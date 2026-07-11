pkgname <- "ctlr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ctlr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ctl")
### * ctl

flush(stderr()); flush(stdout())

### Name: ctl
### Title: Clinical Tolerance Limits for Assessing Agreement
### Aliases: ctl

### ** Examples

## Don't show: 
if (requireNamespace("withr", quietly = TRUE)) withAutoprint({ # examplesIf
## End(Don't show)
withr::with_tempdir({
   # Tolerance limit plot
   ctl(ctl_dataset1, idvar = "id", ynew = "y1", yref = "y2", intercept = 5, slope = 0.15,
       tlplot = TRUE, plots_to_file = TRUE, results_to_file = TRUE, outputs_path = ".")

   ctl(ctl_dataset2, idvar = "id", ynew = "y1", yref = "y2", intercept = 1, slope = 0.2,
       seed = 11446158, tlplot = TRUE, plots_to_file = TRUE, results_to_file = TRUE,
       outputs_path = ".")

   # Conditional probability of agreement plot
   ctl(ctl_dataset1, idvar = "id", ynew = "y1", yref = "y2", intercept = 5, slope = 0.15,
       nbsimul = 100, cpaplot = TRUE, plots_to_file = TRUE, results_to_file = TRUE,
       outputs_path = ".")

   ctl(ctl_dataset1, idvar = "id", ynew = "y1", yref = "y2", intercept = 0, slope = 0.15,
       nbsimul = 100, cpaplot = TRUE, pointwise = TRUE, plots_to_file = TRUE,
       results_to_file = TRUE, outputs_path = ".")

   ctl(ctl_dataset2, idvar = "id", ynew = "y1", yref = "y2", intercept = 1, slope = 0.2,
       seed = 11446158, nbsimul = 100, cpaplot = TRUE, simultaneous = TRUE)
})
## Don't show: 
}) # examplesIf
## End(Don't show)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
