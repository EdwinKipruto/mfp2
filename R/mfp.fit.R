# fits an MFP model
#'@import utils
mfp.fit <- function(x, y, weights, offset, cycles,scale,shift,df, center, criterion,
                    xorder,powers, family, method,select, alpha, keep, strata, ftest,
                    control,rownames,nocenter,acdx,verbose){
  #-----------------------------------------------------------------------------
  # STEP 1: Order variable based on xorder
  #-----------------------------------------------------------------------------
  varnames <- names(varorder(x = x, y = y, weights = weights, offset = offset,
                             xorder = xorder, family = family, method = method,
                             strata=strata,control = control,rownames= rownames,
                            nocenter = nocenter)$pvalues)
  xnames <-colnames(x)
  #-----------------------------------------------------------------------------
  #                              STEP 2
  # Make a list of all the initial FP powers that are equal to one. Assign names
  # to scale,alpha, select, and df. Because FP powers are ordered by varnames,
  # the aforementioned vectors should be ordered to match the variables as well.
  #-----------------------------------------------------------------------------
  pp<-setNames(rep(1,ncol(x)), xnames)
  pwrs <- lapply(split(pp, names(pp)), unname)[varnames]
  # Named alpha vector
  alpha <- setNames(alpha, xnames)[varnames]
  # Named select vector
  select <- setNames(select, xnames)[varnames]
  # Force some variables into the model when pvalue criterion is used
  if(!is.null(keep)){
    index <- which(names(select)%in%keep)
    select <- replace(select, list = index, values =rep(1,length(index)))
  }
  # Named df vector
  df <- setNames(df, xnames)[varnames]
  # Named center vector
  center <- setNames(center, xnames)[varnames]
  # adjustment and scaling factors
  shift <- setNames(shift, xnames)[varnames]
  scale <- setNames(scale, xnames)[varnames]
  # logical indicator indicating whether or not a variable has acd
  acdx <- setNames(acdx,xnames)[varnames]
  #-----------------------------------------------------------------------------
  # Check whether acd transformation is required
  #-----------------------------------------------------------------------------
  if(any(acdx)){
  # reset acdx of variables with less than 5 level to false
    acdx <- reset.acd(x, acdx)
  # assign two powers to acd variables (1,NA). The first is for xi, and the
  # second is for acd(xi). NA has been assigned to acd(xi), which will be
  # updated in the MFP cycles.
  acd.variables <- xnames[match(names(which(acdx)), xnames)]
  pwrs.acd <- lapply(1:length(acd.variables), function(x) c(1,NA))
  names(pwrs.acd) <- acd.variables
  # New initial powers with acd having two powers
  pwrs <- modifyList(x=pwrs, val = pwrs.acd)
  # Override df of acd variables by setting them to 4
  df <- replace(df, which(varnames%in%acd.variables), rep(4,length(acd.variables)))
  }

  #-----------------------------------------------------------------------------
  #     STEP 3:  Run the mfp cycles
  #-----------------------------------------------------------------------------
  # Initialize the cycle and a counter that record the stopping point if converged
  j <- 1
  stop.at <- j
  # Set a counter for convergence
  converged <- FALSE
  # Run each cycles and update the initial powers denoted by pwrs
  while(j<=cycles){
    if(verbose) {cat("\nCycle", j)}
    # Each cycle employs pwrs, which are updated after each cycle is completed.
    run.each.cycle <- eachcycle(x = x, y = y,
                                allpowers = pwrs, # acd variables, like FP2, have two powers.
                                df = df,
                                weights = weights,
                                offset = offset,
                                family = family,
                                criterion = criterion,
                                select = select,
                                alpha = alpha,
                                keep = keep,
                                powers = powers,
                                ftest = ftest,
                                control = control,
                                rownames = rownames,
                                strata = strata,
                                nocenter = nocenter,
                                method=method,
                                acdx = acdx, # same length with x
                                verbose = verbose)

    # Estimated powers for the ith cycle
    pwrs.updated <- run.each.cycle
    # Check for convergence
    if(identical(pwrs,pwrs.updated)){
      converged <- TRUE
      pwrs <- pwrs.updated
      stop.at <- j + 1
      cat("\nFractional polynomial fitting algorithm converged after ",j," cycles.","\n")#, "\n\n")
      if (j <= cycles) break
    }else{
      if(j<=cycles){
        # update the powers of the variables at the end of each cycle
        pwrs <- pwrs.updated
        stop.at <- j
        # increment j
        j <- j+1

      }
    }

  }

  # Return warning message if the algorithm failed to converge
  if(!converged){
    stop.at <- stop.at
    warning("No convergence after ", stop.at, " cycles. Results of the last iteration reported", call. = F)
  # Return the mfp model for the last iteration even if not converged
  }
  #=============================================================================
  # Table showing FP Power for each variable---TO MOVE TO mfpa()
  #=============================================================================
  # status: 1 = in the model while 0 = out/removed
  status <- sapply(pwrs, function(x) ifelse(all(is.na(x)),0,1))
  # Final degrees of freedom
  dfx <- sapply(pwrs, df.final)
  # Add "A" to names of acd variables
  xnam <- sapply(names(pwrs), function(x) ifelse(acdx[x], paste0("(A)",x), x))
  # Matrix of FP powers
  mfp.powers <- powermat(pwrs)
  if(criterion=="pvalue"){
    # combine select and alpha vectors into a matrix and cbind with FP powers
    matsel <- do.call(cbind, list(df, select, alpha, status, dfx))
    colnames(matsel) <- c("df.initial","select","alpha", "status", "df.final")
    rownames(matsel) <- xnam
    fp.table <- cbind(matsel, mfp.powers)
  }else{
    if(criterion=="AIC"){
      fp.table <- data.frame(df.initial = df, AIC = rep("AIC", length(pwrs)),
                             status =status, df.final = dfx, mfp.powers)
      rownames(fp.table) <- xnam

    }else{
      fp.table <- data.frame(df.initial = df, AIC = rep("BIC", length(pwrs)),
                             status =status,df.final = dfx, mfp.powers)
      rownames(fp.table) <- xnam

    }
  }

  # Table of scaling and shifting factors
  matssc <- data.frame(shift, scale, center)
  colnames(matssc) <- c("shift","scale","center")
  #=============================================================================
  # Fit the final model with transformed x if nonlinear functions were selected--TO MOVE TO mfpa()
  #=============================================================================
  # Transform x using the final FP powers selected. x has already been shifted and scaled
  X <- xtransform(x = x, power.list = pwrs, center = center, acdx = acdx)

  # Use the transformed x and fit the final model
  # modelfit <- model.fit(x = X, y = y,family = family,weights = weights, offset = offset,
  #                  method = method, strata = strata,control = control,
  #                  rownames = rownames,resid = resid,nocenter = nocenter)$fit
  modelfit <- model.fit(x = X, y = y,family = family,weights = weights, offset = offset,
                        method = method, strata = strata,control = control,
                        rownames = rownames,nocenter = nocenter)
  # untransformed and scaled x
  x <- x[,names(pwrs[!sapply(pwrs, function(x) all(is.na(x)))]), drop = F]
  # variance-covariance matrix
  # if(family!="cox")
  # summary.glm(modelfit)$cov.scaled
  if(family!="cox"){
    family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    fit <- list(coefficients = modelfit$fit$coefficients,
                residuals = modelfit$fit$residuals, fitted.values = modelfit$fit$fitted.values,
                effects =modelfit$fit$effects, R = modelfit$fit$R, rank = modelfit$fit$rank,
                qr =modelfit$fit$qr, family = family, na.action = NULL,
                linear.predictors = modelfit$fit$linear.predictors,
                deviance =  modelfit$fit$deviance, aic =  modelfit$fit$aic,
                null.deviance =  modelfit$fit$null.deviance, weights = modelfit$fit$weights,
                prior.weights =  modelfit$fit$prior.weights, df.residual = modelfit$fit$df.residual,
                  df.null = modelfit$fit$df.null, y = y, fp.table = fp.table,
                shiftscalecenter = matssc,  pwrs = pwrs,acd = acdx, X =X, x = x)
  }else{
    fit <- list(coefficients = modelfit$fit$coefficients,
                residuals = modelfit$fit$residuals, family = "cox",
                linear.predictors = modelfit$fit$linear.predictors,
                var = modelfit$fit$var,loglik = modelfit$fit$loglik,
                score = modelfit$fit$score, means = modelfit$fit$means,
                method <- modelfit$fit$method, class = modelfit$fit$class,
                acd = acdx, y = y, X =X, x = x, n = nrow(y),
                nevent = sum(y[,ncol(y)]),
                #na.action = options()$na.action, # default in coxph. set. not in use anywhere because the user must take care of missing data
                fail = if(is.character(modelfit$fit)){"fail"}, # to work on this later-might not be correct
                fp.table = fp.table,shiftscalecenter = matssc,
                pwrs = pwrs)
  }
    # incase cox fails

  return(fit)
}








