# calculates index of the best model when acd transformation is desired
# pvalue-pvalues of: Null vs M1, lin(xi) vs M1, fp1(xi) vs M1, FP1(acd(xi)) vs M1
# and lin(acd(xi)) vs FP1(acd(xi)) in that order. Note M1 = FP1(xi, acd(xi))
best.model.index.acd <- function(pvalue, select, alpha) {
  if (pvalue[1] > select) {
    index.bestmodel <- 1 # Null vs M1: if not sig then NULL model is chosen
  } else {
    # Linear(xi) vs M1: if not sig then choose linear xi is chosen
    if (pvalue[2] > alpha) {
      index.bestmodel <- 2
    } else {
      # FP1(xi) vs M1: if not sig then choose FP1(xi)
      if (pvalue[3] > alpha) {
        index.bestmodel <- 3
      } else {
        # FP1(acd(xi)) vs M1: if not sig then choose FP1(acd(xi))
        if (pvalue[4] > alpha) {
          # since the FP1(axi) vs M1 is not sig, we can check whether linear(axi)
          # is better than FP1(axi)
          if (pvalue[5] > alpha) {
            index.bestmodel <- 6
          } else {
            index.bestmodel <- 4
          }
        } else {
          index.bestmodel <- 5
        }
      }
    }
  }
  index.bestmodel
}
