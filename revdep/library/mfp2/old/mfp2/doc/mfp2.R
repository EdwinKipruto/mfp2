## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(
	message = FALSE,
	warning = FALSE,
	tidy = TRUE,
	fig.pos = "h"
)
oldwidth <- options()$width
options(width = 100)
knitr::opts_chunk$set(tidy.opts=list(width.cutoff=75))

library(mfp2)
library(ggplot2)
library(survival)
library(formatR)

## ----include=FALSE--------------------------------------------------------------------------------
# the code in this chunk enables us to truncate the print output for each
# chunk using the `out.lines` option
# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
        
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})

## ----eval=TRUE, echo=FALSE------------------------------------------------------------------------
#==============================================================================
# 8 FP1 FUNCTIONS 
#===============================================================================
x <- seq(0.05, 1.05, length.out = 1000)
funx <- function(x, power){
  ifelse(power == rep(0, length(x)), log(x), x^power)
}
s <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
# transform x using the powers in s
outx <- lapply(s, function(s) funx(x = x, power = s))

# multiply the first 3 with -1 so that the function increase rather than decrease
# due to negative powers. these are betas i,e y = beta*x^p
k <- c(-1, -1, -1, 1, 1, 1, 1, 1)
datx <- matrix(unlist(outx), ncol = length(s))
datax <- datx %*% diag(k)
# standardize the data so that y is in the same range
dataxx <- apply(datax, 2, function(x) (x - min(x)) / (max(x) - min(x)))
colnames(dataxx) <- c(paste0("x", 1:length(s)))
dataxx <- as.data.frame(dataxx)
dataxx$x <- x
width <- 0.5
fig1 <- ggplot(dataxx, aes(x)) + 
  geom_line(aes(y = x1, colour = "x1"), linewidth = width, color = "#006400") + 
  geom_line(aes(y = x2, colour = "x2"), linewidth = width, color ="#ff0000") + 
  geom_line(aes(y = x3, colour = "x3"), linewidth = width, color ="#ffd700") +
  geom_line(aes(y = x4, colour = "x4"), linewidth = width, color ="#ff00ff") +
  geom_line(aes(y = x5, colour = "x5"), linewidth = width, color ="#ffb6c1") +
  geom_line(aes(y = x6, colour = "x6"), linewidth = width, color ="#00ff00") +
  geom_line(aes(y = x7, colour = "x7"), linewidth = width, color ="#0000ff") +
  geom_line(aes(y = x8, colour = "x8"), linewidth = width, color ="#000000") +
  geom_text(aes(x = 0.175, y = 1,size = 5, label = "-2"), color = "#006400") +
  geom_text(aes(x = 0.25, y = 0.9,size = 5, label = "-1"), color = "#ff0000") +
  geom_text(aes(x = 0.33, y = 0.83,size = 5, label = "-0.5"), color = "#ffd700") +
  geom_text(aes(x = 0.4, y = 0.725,size = 5, label = "0"), color = "#ff00ff") +
  geom_text(aes(x = 0.5, y = 0.65,size = 5, label = "0.5"), color = "#ffb6c1") +
  geom_text(aes(x = 0.59, y = 0.575,size = 5, label = "1"), color = "#00ff00") +
  geom_text(aes(x = 0.7, y = 0.475,size = 5, label = "2"), color = "#0000ff") +
  geom_text(aes(x = 0.75, y = 0.4,size = 5, label = "3"), color = "#000000") +
  labs(x="x", y="Fractional polynomial f(x)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title=element_text(size=18),
    axis.text.x = element_blank()) 

# Define functions
f1 <- function(x) 3 - 10*x^2 + 4*x^3
f2 <- function(x) 20 - 15.4*x^2 + 4*x^3
f3 <- function(x) -20 + 6*log(x) + 6*log(x)*log(x)
f4 <- function(x) 20 + 0.3*(x^-2) - 4*(x^-1)
f5 <- function(x) -10 + 5*(x^0.5) + 14*(x^-0.5)
f6 <- function(x) 33 + 19*log(x) - 7*(x^2)
f7 <- function(x) -10 + 10*(x - 1.5) + 10*(x - 1.5)^2

# Create data frames for each function
x <- seq(0.1, 3, by = 0.01)
df1 <- data.frame(x = x, y = f1(x))
df2 <- data.frame(x = x, y = f2(x))
df3 <- data.frame(x = x, y = f3(x))
df4 <- data.frame(x = x, y = f4(x))
df5 <- data.frame(x = x, y = f5(x))
df6 <- data.frame(x = x, y = f6(x))
df7 <- data.frame(x = x, y = f7(x))

# Plot functions
fig2 <- ggplot() +
    geom_line(data = df1, aes(x = x, y = y), color = "#e41a1c") +
    geom_line(data = df2, aes(x = x, y = y), color = "#377eb8") +
    geom_line(data = df3, aes(x = x, y = y), color = "#4daf4a") +
    geom_line(data = df4, aes(x = x, y = y), color = "#984ea3") +
    geom_line(data = df5, aes(x = x, y = y), color = "#ff7f00") +
    geom_line(data = df6, aes(x = x, y = y), color = "#FF00FF") +
    geom_line(data = df7, aes(x = x, y = y), color = "#a65628") +
    scale_x_continuous(limits = c(0, 3)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-25,40)) + theme_classic() +
    geom_text(aes(x = 0.75, y = 2.4,size = 5, label = "(2, 3)"), color = "#e41a1c") +
    geom_text(aes(x = 1.6, y =1.95,size = 5, label = "(2, 3)"), color = "#377eb8")+
    geom_text(aes(x = 0.7, y = -18,size = 5, label = "(0, 0)"), color = "#4daf4a") +
    geom_text(aes(x = 1.25, y = 20,size = 5, label = "(-2, -1)"), color = "#984ea3") +
    geom_text(aes(x = 1.3, y = 12,size = 5, label = "(-0.5, 0.5)"), color = "#ff7f00") +
    geom_text(aes(x = 1.25, y = 29.5,size = 5, label = "(0, 2)"), color = "#FF00FF") +
    geom_text(aes(x = 1, y = -9.5,size = 5, label = "(1, 2)"), color = "#a65628") +
    ylab(" ") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      axis.title = element_text(size=18),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank()
    )

figure <- patchwork::wrap_plots(fig1, fig2, 
                               ncol = 2, nrow = 1, 
                               widths = 8, heights = 2)


## ----fig1, echo=FALSE, fig.cap="Illustration of the flexibility of the FP family. 8 FP1 functions (left panel) and a subset of FP2 functions (right panel). FPs are global functions. See section 1.3.1.4 for a proposed extension of local features", fig.height= 5, fig.width=8----
figure

## ----eval=FALSE-----------------------------------------------------------------------------------
#  # install the package
#  install.packages("mfp2")
#  
#  # load the package
#  library(mfp2)

## ----out.lines = 10-------------------------------------------------------------------------------
# Load prostate data
data("prostate")
head(prostate, 6)

# create covariate matrix x and numeric vector y 
x <- as.matrix(prostate[,2:8])
y <- as.numeric(prostate$lpsa)

## -------------------------------------------------------------------------------------------------
# default interface
fit_default <- mfp2(x, y, 
                    verbose = FALSE)
summary(fit_default)

# formula interface
fit_formula <- mfp2(lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) + fp(cp),
                    data = prostate, 
                    verbose = FALSE)
summary(fit_formula)

## ----eval = TRUE----------------------------------------------------------------------------------
# minimum values for each covariate
apply(x, 2, min)

# shifting values for each covariate
apply(x, 2, find_shift_factor)

## ----eval=TRUE------------------------------------------------------------------------------------
# shift x if nonpositive values exist
shift <- apply(x, 2, find_shift_factor)
xnew <- sweep(x, 2, shift, "+")

# scaling factors
apply(xnew, 2, find_scale_factor)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # Default interface
#  mfp2(x,y,
#       shift = c(0, 0, 1, 0, 0, 0, 0),
#       scale = c(10, 1, 100, 10, 100, 1, 10))
#  
#  # Formula interface
#  mfp2(lpsa ~ fp(age, shift = 0, scale = 10) + svi + fp(pgg45, shift = 1, scale = 100) + fp(cavol, shift = 0, scale = 10) + fp(weight, shift = 0, scale = 100) + fp(bph, shift = 0, scale = 1) + fp(cp, shift = 0, scale = 10),
#       data = prostate)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # Default Interface
#  mfp2(x,y, df = c(4, 1, 4, 4, 4, 4, 4))
#  
#  # Formula Interface
#  mfp2(lpsa ~ fp(age, df = 4) + svi + fp(pgg45, df = 4) + fp(cavol, df = 4) + fp(weight, df = 4) + fp(bph, df = 4) + fp(cp, df = 4),
#       data = prostate)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # Default Interface
#  mfp2(x,y)
#  
#  # Formula Interface, all continuous variables passed to the fp() function
#  mfp2(lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) + fp(cp),
#       data = prostate)
#  
#  # Formula Interface, `cp` not passed to the fp() function
#  mfp2(lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) + cp,
#       data = prostate)
#  
#  # Formula Interface, all covariates not passed to the fp() function
#  mfp2(lpsa ~ age + svi + pgg45 + cavol + weight + bph + cp,
#       data = prostate)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # Default Interface
#  mfp2(x,y,
#       select = rep(0.05, ncol(x)),
#       alpha = rep(0.05, ncol(x)))
#  
#  # Formula Interface
#  mfp2(lpsa ~ fp(age, select = 0.05) + svi + fp(pgg45, select = 0.05) + fp(cavol, select = 0.05) +
#         fp(weight, select = 0.05) + fp(bph, select = 0.05) + fp(cp, select = 0.05),
#       select = 0.05, alpha = 0.05,
#       data = prostate)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # default Interface
#  # Set select to 1 for variables "age" and "svi"
#  mfp2(x,y,
#       select = c(1, 1, 0.05, 0.05, 0.05, 0.05, 0.05),
#       alpha = rep(0.05, ncol(x)))
#  
#  # use keep argument
#  mfp2(x,y,
#       select = rep(0.05, ncol(x)),
#       alpha = rep(0.05, ncol(x)),
#       keep = c("age", "svi"))
#  
#  # formula Interface
#  # use fp() function and set select to 1 for "age" and "svi"
#  mfp2(lpsa ~ fp(age, select = 1) + fp(svi, df = 1, select = 1) + fp(pgg45, select = 0.05) +
#         fp(cavol, select = 0.05) + fp(weight, select = 0.05) + fp(bph, select = 0.05) + fp(cp, select = 0.05),
#       select = 0.05, alpha = 0.05, data = prostate)
#  
#  # use keep argument
#  mfp2(lpsa ~ fp(age, select = 0.05) + svi + fp(pgg45, select = 0.05) + fp(cavol, select = 0.05) +
#         fp(weight, select = 0.05) + fp(bph, select = 0.05) + fp(cp, select = 0.05),
#       select = 0.05, alpha = 0.05,
#       keep = c("age", "svi"),
#       data = prostate)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  # Default Interface
#  mfp2(x, y,
#       criterion = "aic",
#       keep = c("age", "svi"))
#  
#  # Formula Interface
#  mfp2(lpsa ~ fp(age) + fp(svi, df = 1) + fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) +  fp(cp),
#       criterion = "aic",
#       keep = c("age", "svi"),
#       data = prostate)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  # Default Interface
#  mfp2(x, y,
#       criterion = "pvalue",
#       ftest = TRUE)
#  
#  # Formula Interface
#  mfp2(lpsa ~ fp(age) + svi + fp(pgg45) + fp(cavol) + fp(weight) + fp(bph) +  fp(cp),
#       criterion = "pvalue",
#       ftest = TRUE,
#       data = prostate)

## ----eval=FALSE-----------------------------------------------------------------------------------
#  # create a list of power terms for covariates "age" and "cavol"
#  powx <- list(age = c(1, 0.5), cavol = c(1, 0, 1))
#  
#  # Default Interface
#  mfp2(x, y,
#       criterion = "pvalue",
#       powers = powx)
#  
#  # Formula Interface
#  mfp2(lpsa ~ fp(age, powers = c(1 , 0.5) ) + svi + fp(pgg45) + fp(cavol, powers = c(1, 0, 1)) + fp(weight) + fp(bph) +  fp(cp),
#       data = prostate)
#  

## ----out.lines = 103------------------------------------------------------------------------------
fit <- mfp2(x, y, 
            criterion = "pvalue", 
            select = 0.05, alpha = 0.05, 
            ftest = TRUE)

## ----eval= TRUE,  tidy = TRUE, tidy.opts = list(width.cutoff = 80)--------------------------------
fit <- mfp2(x, y, verbose = FALSE)
class(fit)

## ----eval= TRUE-----------------------------------------------------------------------------------
print(fit)

## ----eval= TRUE-----------------------------------------------------------------------------------
# extract only regression coefficients
coef(fit)

# display regression coefficients with other statistics
summary(fit)

## ----out.lines = 10, eval= TRUE-------------------------------------------------------------------
 # extract the first five observations from 'x'
new_observations <- x[1:5, ] 
dim(new_observations)

# make predictions for these new observations.
predict(fit, newdata = new_observations)

## ----eval= TRUE,fig.cap="\\label{fig:fig2}Prostate data. Graphical presentation of results from an FP analysis"----
plots <- fracplot(fit)

# plots is a list
class(plots)

# The plot of cavol is the first element of the list
plots[[1]] + ggplot2::ylab("Partial predictor + residuals")


## ----eval= TRUE, fig.cap="\\label{fig:fig3}Prostate data. Illustration on how to use reference points (ref point = 30)"----
plots <- fracplot(fit, type = "contrasts", ref = list(cavol = 30))
plots[[1]] + ggplot2::ylab("Partial predictor")

## ----eval = TRUE, fig.cap = "Art data. Illustration on how to use terms_seq argument in fracplot. In the original data, the function FP(0,3) was selected, equivalent to $\\beta_0 + \\beta_1log(x5) + \\beta_2x5^3$. This function was plotted using the original values of x5 in the data (right panel) and a sequence of equidistant new data points using the range of the original x5 data (left panel). The two functions are identical up to x5 = 62. However, beyond this point, a linear trend is estimated when the original values are used due to lack of data points between 62.21 and 235. If the new data is used, it results in a curve, clearly depicting the FP2 function."----

# load art data
data("art")

# fit mfp model using art data
xx <- as.matrix(art[,"x5", drop = FALSE])
yy <- as.numeric(art$y)
fit2 <- mfp2(xx, yy, verbose = FALSE)

# generate plot using original data 
plot1 <- fracplot(fit2)
plot1[[1]] <- plot1[[1]] + 
  ggplot2::ylab("Partial predictor + residuals") 

# generate plot using sequence of data
plot2 <- fracplot(fit2, terms_seq = "equidistant")
plot2[[1]] <- plot2[[1]] + 
  ggplot2::ylab("Partial predictor + residuals") 

# combine plots
patchwork::wrap_plots(plot2[[1]], plot1[[1]], ncol = 2, widths = 8, heights = 4)


## ----out.lines = 103------------------------------------------------------------------------------
# load art data
data("art")

# order the levels of the variable such that Levels: 1 < 2 < 3
art$x4 <- factor(art$x4, ordered = TRUE, levels = c(1, 2, 3))

# convert x9 to factor and assign level 1 as reference group
art$x9 <- factor(art$x9, ordered = FALSE, levels = c(1, 2, 3))

# create dummy variables for x4 and x9 based on ordinal an categorical coding and drop 
# the original variable
art <- create_dummy_variables(art, 
                              var_ordinal = c("x4"),
                              var_nominal = c("x9"),
                              drop_variables = TRUE)

# display the first 20 observations
head(art, 10)

## ----out.lines = 103------------------------------------------------------------------------------
# create matrix x and outcome y
x <- as.matrix(art[,-1])
y <- as.numeric(art$y)

# fit mfp using default interface with default parameters
fit <- mfp2(x,y, verbose = FALSE)
fit

## ----out.lines = 103------------------------------------------------------------------------------
# load pima data
data("pima")
head(pima)

# matrix x
x <- as.matrix(pima[,2:9])
# outcome y
y <- as.vector(pima$y)

# fit mfp
fit <- mfp2(x, y, 
            family = "binomial", 
            verbose = FALSE)
fit

## ----out.lines = 103, fig.show='hold', fig.cap="Pima data. Graphical presentation of results from an MFP analysis", echo=FALSE----
plots <- fracplot(fit)

# rename y-label for all plots
plots <- lapply(plots, function(v) {
  v + ggplot2::ylab("Partial pred + residuals")
})

# combine multiple ggplot objects into a single figure
patchwork::wrap_plots(plots, ncol = 2)

## ----out.lines = 103------------------------------------------------------------------------------
# load gbsg data
data("gbsg")

# data preparation 
# create dummy variable for grade using ordinal coding
gbsg <- create_dummy_variables(gbsg, 
                               var_ordinal = "grade",
                               drop_variables = TRUE)

# covariate matrix x 
x <- as.matrix(gbsg[,-c(1, 6, 10, 11)])
head(x, 10)

# use Surv() function to create outcome y
y <- survival::Surv(gbsg$rectime, gbsg$censrec)

# fit mfp and keep hormon in the model
fit1 <- mfp2(x, y, family = "cox",keep = "hormon", 
             control = coxph.control(iter.max = 50),
                          verbose = FALSE)
fit1

## ----out.lines = 103------------------------------------------------------------------------------
# remove nodes and include enodes
x <- as.matrix(gbsg[,-c(1, 5, 10, 11)])

# fit mfp and keep hormon in the model
fit2 <- mfp2(x, y, family = "cox", keep = "hormon",
             powers = list(enodes = c(0.5, 1, 2, 3)), 
             control = coxph.control(iter.max = 50),
             verbose = FALSE)
fit2

## ----out.lines = 103------------------------------------------------------------------------------
# using default interface
fit2 <- mfp2(x[,-7], y, family = "cox", 
             strata = x[,7], 
             verbose = FALSE)
fit2

# using formula interface
fit3 <- mfp2(Surv(rectime,censrec)~ age + meno + size + nodes + pgr + er + grade_1 + 
               grade_2 + strata(hormon),
             family = "cox", 
             data =gbsg, 
             verbose = FALSE)

## ----out.lines = 103------------------------------------------------------------------------------
# Generate artificial data with sigmoid relationship
set.seed(54) 
n <- 500  
x <- matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x") )

# Apply sigmoid transformation to x
sigmoid <- function(x) {
  1 / (1 + exp(-1.7*x))
}

# Generate y with sigmoid relationship to x
y <- as.numeric(20*sigmoid(x) + rnorm(n, mean = 0, sd = 0.5))

## ----out.lines = 103------------------------------------------------------------------------------
# default interface
fit1 <- mfp2(x, y, 
             acdx = "x", 
             verbose = FALSE)

# formula interface
datax <- data.frame(y, x)
fit2 <- mfp2(y ~ fp(x, acdx = TRUE),
             data = datax,
             verbose = FALSE)

# display selected power terms
fit2

# fit usual mfp: bad idea but useful for illustration
fit3 <- mfp2(y ~ fp(x, acdx = FALSE), 
             data = datax, 
             verbose = FALSE)
fit3

## ----echo=FALSE,eval=TRUE-------------------------------------------------------------------------
p1 <- fracplot(fit1, terms = "x")[[1]] + ggplot2::ylab("y")
p2 <- fracplot(fit3, terms = "x")[[1]] + ggplot2::ylab("y")
figurex <- patchwork::wrap_plots(p1, p2, 
                               ncol = 2, nrow = 1, 
                               widths = 8, heights = 2)

## ----echo=FALSE, fig.cap="The selected function using the ACD transformation is FP1(1,1) represented by $\\beta_0 + \\beta_1x^1 +\\beta_2a^1$(left panel), while the usual MFP selected FP2(3,3), equivalent to $\\beta_0 + \\beta_1x^3 +\\beta_2x^3\\log(x)$ (right panel). See text for more details", fig.width=8, fig.height=4----
figurex

## ----eval=TRUE, echo=FALSE----------------------------------------------------
# reset width to the original value
options(width=oldwidth)

