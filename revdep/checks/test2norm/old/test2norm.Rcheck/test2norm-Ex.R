pkgname <- "test2norm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('test2norm')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("PsychTestData")
### * PsychTestData

flush(stderr()); flush(stdout())

### Name: PsychTestData
### Title: Neuropsychological test data
### Aliases: PsychTestData
### Keywords: datasets

### ** Examples

data(PsychTestData)
test2norm(data = PsychTestData, test = "rawscore",
          test.min = 0, test.max = 36, test.better = "High",
          group.id = "group", control.id = "control",
          demographics = c("age", "male"))



cleanEx()
nameEx("raw2scaled")
### * raw2scaled

flush(stderr()); flush(stdout())

### Name: raw2scaled
### Title: Convert raw neuropsychological test scores to scaled scores.
### Aliases: raw2scaled

### ** Examples

data(PsychTestData)
raw2scaled(data = PsychTestData, test = "rawscore",
           test.min = 0, test.max = 36, test.better = "High",
           group.id = "group", control.id = "control")



cleanEx()
nameEx("score2adjust")
### * score2adjust

flush(stderr()); flush(stdout())

### Name: score2adjust
### Title: Convert neuropsychological test scores to demographically
###   adjusted norms.
### Aliases: score2adjust

### ** Examples

data(PsychTestData)
PsychTestData$scaledscore <- raw2scaled(data=PsychTestData, test="rawscore",
                                        test.min=0, test.max=36,
                                        test.better="High", group.id="group",
                                        control.id="control")[[2]]
score2adjust(data = PsychTestData, test.score = "scaledscore",
             group.id = "group", control.id = "control",
             demographics = c("age", "male"))



cleanEx()
nameEx("test2norm")
### * test2norm

flush(stderr()); flush(stdout())

### Name: test2norm
### Title: Convert raw neuropsychological test scores to demographically
###   adjusted norms.
### Aliases: test2norm

### ** Examples

data(PsychTestData)
test2norm(data = PsychTestData, test = "rawscore",
          test.min = 0, test.max = 36, test.better = "High",
          group.id = "group", control.id = "control",
          demographics = c("age", "male"))



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
