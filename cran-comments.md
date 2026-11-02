## Submission

This is an update of mfp2 from 1.0.1 to 1.1.0. The main changes are new
model families (negative-binomial, Gamma, inverse-Gaussian, multinomial,
ordinal, parametric survival, Fine-Gray and GEE), interaction analysis
with `mfpi()`, spike-at-zero modelling, and `plot.mfp2()`. See NEWS.md.

## Maintainer email change

The maintainer is unchanged (Edwin Kipruto), but the email address has
changed from edwin.kipruto@uniklinik-freiburg.de to
kiprutoedwin2626@gmail.com. I no longer have access to the old address
because I have changed jobs, so I cannot confirm the change from it.
Co-author Willi Sauerbrei (wilhelm.sauerbrei@uniklinik-freiburg.de) can
confirm the change on my behalf if needed.

## Test environments

* Local Windows 10 x64, R 4.6.1: `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 0 notes

CRAN's incoming checks will report a NOTE about the new maintainer email
address; this is explained in the section above.

## Reverse dependencies

Checked 3 reverse dependencies: ctlr, midoc, and test2norm. No new problems
were found.
