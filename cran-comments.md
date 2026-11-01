## Submission

This is an update of mfp2 from version 1.0.1 to 1.1.0. It adds new model
families, interaction analysis with `mfpi()`, spike-at-zero modelling, and a
`plot()` method. See NEWS.md for details.

## Maintainer email change

The maintainer's email address has changed from
edwin.kipruto@uniklinik-freiburg.de to edwin.kipruto@email.uni-freiburg.de.
I no longer have access to the previous address <REASON, e.g. the account was
closed when my position changed>, so I cannot send confirmation from it.
Co-author Willi Sauerbrei (wilhelm.sauerbrei@uniklinik-freiburg.de) can
confirm the change if needed.

## Removed exports

Eleven low-level helper functions are no longer exported (listed in NEWS.md).
They are still available internally. mfp2 has no reverse dependencies on CRAN,
so no other packages are affected.

## Test environments

* local: R 4.x.x, <OS>
* win-builder: R-devel and R-release
* <add others, e.g. R-hub or GitHub Actions: macOS, Ubuntu>

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New maintainer: Edwin Kipruto <edwin.kipruto@email.uni-freiburg.de>
  Old maintainer(s): Edwin Kipruto <edwin.kipruto@uniklinik-freiburg.de>

## Reverse dependencies

There are no reverse dependencies on CRAN.
