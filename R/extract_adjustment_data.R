#' Function to extract and transform adjustment variables
#' 
#' @param x a matrix of predictors that includes the variable of interest `xi`. 
#' It is assumed that continuous variables have already been shifted and scaled.
#' @param xi name of the continuous predictor for which the FP function will be
#' estimated. There are no binary or two-level variables allowed. All variables
#' except `xi` are referred to as "adjustment variables".
#' @param fp_powers a named list of FP powers of all variables of interest, 
#' including `xi`. Note that these powers are updated during backfitting or MFP 
#' cycles.
#' @param df a numeric vector of degrees of freedom for `xi`.
#' @param powers a set of allowed FP powers.
#' 
#' @details
#' After extracting the adjustment variables this function, using their
#' corresponding FP powers stored in `fp_powers`, transforms them. 
#' This is necessary When evaluating x of interest, as we must account for other
#' variables, which can be transformed or untransformed, depending on the
#' individual powers. It's worth noting that some powers can be NA, indicating
#' that the variable has been left out of the adjustment variables. It also
#' returns the FP data, which is dependent on the degrees of freedom. For example,
#' `df = 2` is equivalent to FP degree one, resulting in the generation of 8 
#' variables. If `acdx` is set to `TRUE`, however, 64 variables are generated.
#' 
#' @return 
#' A list containing the following elements: 
#' 
#' * `powers_fp`: fp powers used for `data_fp`.
#' * `data_fp`: all possible fp transformations for `xi`, see the `data` 
#' component of the output of [generate_transformations_fp()] and 
#' [generate_transformations_acd()]. 
#' * `powers_adj`: fp powers for adjustment variables in `data_adj`.
#' * `data_adj`: adjustment data, i.e. transformed input data for adjustment
#' variables.
transform_data_step <- function(x,
                                    xi,
                                    fp_powers,
                                    df, 
                                    powers,
                                    acdx) {

  # sort x based on the names of the powers
  names_fp_powers <- names(fp_powers)
  x <- x[, names_fp_powers, drop = FALSE]
  
  # extract the matrix of adjustment variables
  vars_adj <- names_fp_powers[!(names_fp_powers %in% xi)]
  x_adj <- x[, vars_adj, drop = FALSE]
  
  # transformations of adjustment variables
  powers_adj <- fp_powers[vars_adj]
  acdx_adj <- unname(acdx[vars_adj])
  
  # generate adjustment data 
  # check whether all adjustment powers = NA
  if (all(is.na(unlist(powers_adj, use.names = FALSE)))) {
    # all adjustment variables were eliminated in MFP backfitting process
    data_adj <- NULL
    powers_adj <- NULL
  } else {
    data_adj <- vector(mode = "list", length = ncol(x_adj))
    
    for (i in 1:ncol(x_adj)) {
      if (acdx_adj[i]) {
        data_adj[[i]] <- transform_vector_acd(
          x = x_adj[, i, drop = TRUE], power = powers_adj[[i]],
          powers = powers
        )
      } else {
        data_adj[[i]] <- transform_vector_fp(
          x = x_adj[, i, drop = TRUE], power = powers_adj[[i]]
        )
      }
    }
    
    # combine into data.frame
    # note some variables may have been extended to more than one column
    # in the loop above due to transformation
    data_adj <- do.call(cbind, data_adj)
    
    # assign arbitrary names to adjustment matrix 
    # the names of adjustment variables at this stage are not relevant
    colnames(data_adj) <- paste0("var_adj", 1:ncol(data_adj))
  }
  
  # generate fp data
  data_xi <- x[, xi, drop = TRUE]
  
  if (length(unique(data_xi)) <= 3 || df == 1) {
    # if a variable has less than 4 levels we do not generate FP data
    # when df = 1 we do not generate FP data because we assume linearity
    data_fp <- data_xi
    powers_fp <- 1
  } else {
    if (acdx[xi]) {
      fpd <- generate_transformations_acd(data_xi, powers = powers)
    } else {
      # note that degree is df / 2
      fpd <- generate_transformations_fp(data_xi, degree = df / 2, powers = powers)
    }
    data_fp <- fpd$data
    powers_fp <- fpd$powers
  }
  
  list(
    powers_fp = powers_fp, 
    data_fp = data_fp, 
    powers_adj = powers_adj,
    data_adj = data_adj
  )
}
