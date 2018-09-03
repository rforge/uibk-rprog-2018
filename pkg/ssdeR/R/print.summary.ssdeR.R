#################
# print.summary.ssdeR
#################
print.summary.ssdeR <- function (x, digits = max(3, getOption("digits") ), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") *
                                                        0.85)), "", sep = "\n")
  cat("\n")
  if (!x$firststage$code) {
    cat("model did not converge\n")
  }else{cat(paste0(x$firststage$message))}
  cat("\n")


  cat(paste("Standardized residuals:\n", sep = ""))
  print(structure(round(as.vector(quantile(x$residuals)),
                        digits = digits), .Names = c("Min", "1Q", "Median",
                                                     "3Q", "Max")))
  if (NROW(x$coefficients$treatment)) {
    cat(paste("\nCoefficients (treatment model):\n", sep = ""))
    printCoefmat(x$coefficients$treatment, digits = digits,
                 signif.legend = FALSE)
  }
  else cat("\nNo coefficients (in treatment model)\n")
  if (NROW(x$coefficients$selection)) {
    cat(paste("\nCoefficients (selection model):\n", sep = ""))
    printCoefmat(x$coefficients$selection, digits = digits,
                 signif.legend = FALSE)
  }
  else cat("\nNo coefficients ( in selection model)\n")
  if (NROW(x$coefficients$outcome)) {
    cat(paste("\nCoefficients (outcome model):\n", sep = ""))
    printCoefmat(x$coefficients$outcome, digits = digits,
                 signif.legend = FALSE)
  }
  else cat("\nNo coefficients ( in outcome model)\n")

  if (NROW(x$coefficients$Aux.Param)) {
    cat(paste("\nAuxiliary Parameters:\n", sep = ""))
    printCoefmat(x$coefficients$Aux.Param, digits = digits,
                 signif.legend = FALSE)
  }
  else cat("\nNo Auxiliary Parameters\n")

  if (getOption("show.signif.stars") & any(do.call("rbind",
                                                   x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
        "\n")
  cat("\n")

  if (is.numeric(x$df.residual)) {
    cat(paste("df.residual: ", format(x$df.residual, digits = digits),
              "\n", sep = ""))
  }

  cat("Log-likelihood:", formatC(x$loglik, digits = digits),
      "on", x$DF1, "Df\n")
  cat("\n")

  cat("AIC:", formatC(x$AIC, digits = digits))
  cat("\n")
  cat("BIC:", formatC(x$BIC, digits = digits))
  cat("\n")

  cat("---",
      "\n")


  cat("Log-likelihood (Bivariate Probit):", formatC(x$firststage$maximum, digits = digits),
      "on", x$DF2, "Df\n")
  cat(paste("Number of iterations in First Stage", x$firststage$type, ":",
            x$firststage$iterations, "\n"))

  invisible(x)
}
