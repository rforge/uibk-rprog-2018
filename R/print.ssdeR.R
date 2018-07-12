################
# print.ssdeR
################
print.ssdeR <- function (x, digits = max(3, getOption("digits") - 3), ...)
{

  k <- length(x$coefficients)
  kt <- length(x$firststage$estimate)
  ks <- length(x$firststage$estimate)
  mmo <- model.matrix(x, "o")
  mmt <- model.matrix(x, "t")
  mms <- model.matrix(x, "s")

  cfo <- x$coefficients
  cft <- x$firststage$estimate[1:NCOL(mmt)]
  cfs <- x$firststage$estimate[(NCOL(mmt)+1):(NCOL(mmt)+(NCOL(mms) + 2))]

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
                                                        0.85)), "", sep = "\n")

  if (k) {
    cat(paste("Outcome Coefficients:\n", sep = ""))
    print.default(format(cfo[], digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n")
    if (ks) {
      cat(paste("Selection Coefficients:\n", sep = ""))
      print.default(format(cfs[], digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\n")
    }
    if (kt) {
      cat(paste("Treatment Coefficients:\n", sep = ""))
      print.default(format(cft[], digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\n")
    }
  }
  cat("\n")
  if (is.numeric(x$df.residual)) {
    cat(paste("df.residual: ", format(x$df.residual, digits = digits),
              "\n", sep = ""))
  }
  cat("\n")
  invisible(x)
}
