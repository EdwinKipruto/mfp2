# -----------------------------------------------------------------------------
# Package-wide verbose-output helpers -----------------------------------------
# -----------------------------------------------------------------------------

#' Emit Package Verbose Output
#'
#' Emits fitting-time verbose output using \code{cat()} for speed.
#'
#' This helper is intended for progress, diagnostics, and verbose fitting output
#' inside \code{if (verbose)} blocks. Do not use it in S3 print methods, where
#' direct \code{cat()} and \code{print()} calls are also appropriate.
#'
#' @param ... Objects passed to \code{cat()}.
#' @param appendLF Logical scalar. Whether to append a newline.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
mfp2_message <- function(..., appendLF = TRUE) {
  cat(..., if (appendLF) "\n" else "", sep = "")
  invisible(NULL)
}


#' Emit Multiple Package Verbose Lines
#'
#' Emits a character vector of already formatted verbose-output lines.
#'
#' @param lines Character vector. Each element is one output line.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
mfp2_message_lines <- function(lines) {
  if (is.null(lines) || length(lines) == 0L) {
    return(invisible(NULL))
  }
  
  cat(as.character(lines), sep = "\n")
  cat("\n")
  
  invisible(NULL)
}


#' Capture Printed Output and Emit It as Verbose Output
#'
#' Evaluates an expression, captures its printed output, and emits the captured
#' text using \code{cat()}.
#'
#' @param expr Expression whose printed output should be captured.
#'
#' @return Invisibly returns the evaluated result of \code{expr}.
#'
#' @keywords internal
#' @noRd
mfp2_message_capture <- function(expr) {
  result <- NULL
  
  lines <- utils::capture.output({
    result <- eval.parent(substitute(expr))
  })
  
  mfp2_message_lines(lines)
  
  invisible(result)
}