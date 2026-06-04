#' Save capture data and trap positions to CSV files
#'
#' @param captures Integer matrix `[T x P]` of daily captures. Defaults to
#'   the built-in `mosquito_captures` dataset.
#' @param trap_positions Numeric matrix `[P x 2]` of trap coordinates.
#'   Defaults to `mosquito_trap_positions`.
#' @param dir Directory where CSV files will be written. Default `"."`.
#'
#' @return Invisibly returns the paths of the two files written.
#' @export
save_data_csv <- function(
    captures       = NULL,
    trap_positions = NULL,
    dir            = "."
) {
  if (is.null(captures))       captures       <- get("mosquito_captures",       envir = asNamespace("mosquitoSIT"))
  if (is.null(trap_positions)) trap_positions <- get("mosquito_trap_positions", envir = asNamespace("mosquitoSIT"))

  T_days <- nrow(captures)
  P      <- ncol(captures)

  # --- Captures CSV ---
  captures_df <- as.data.frame(captures)
  colnames(captures_df) <- paste0("trap_", seq_len(P))
  captures_df <- cbind(day = seq_len(T_days), captures_df)

  captures_path <- file.path(dir, "mosquito_captures.csv")
  write.csv(captures_df, captures_path, row.names = FALSE)
  message("Captures saved to: ", captures_path)

  # --- Trap positions CSV ---
  traps_df <- as.data.frame(trap_positions)
  colnames(traps_df) <- c("x_coord", "y_coord")
  traps_df <- cbind(trap = seq_len(P), traps_df)

  traps_path <- file.path(dir, "mosquito_trap_positions.csv")
  write.csv(traps_df, traps_path, row.names = FALSE)
  message("Trap positions saved to: ", traps_path)

  invisible(c(captures_path, traps_path))
}


#' Plot prior and posterior distributions for D, LAMBDA and GAMMA
#'
#' Creates a PDF with one panel per parameter, showing the posterior density
#' from MCMC samples overlaid on the prior distribution.
#'
#' @param fit A `CmdStanMCMC` object returned by [run_mosquito_inference()].
#' @param path Path of the output PDF file. Default `"prior_posterior.pdf"`.
#' @param r_piege The fixed trap radius used during inference (used only for
#'   the plot title). Default 3.5.
#'
#' @return Invisibly returns `path`.
#' @export
plot_prior_posterior <- function(fit, path = "prior_posterior.pdf", r_piege = 3.5) {

  draws <- fit$draws(c("D", "LAMBDA", "GAMMA"), format = "matrix")

  params <- list(
    D      = list(samples   = draws[, "D"],
                  prior_fun = function(x) dlnorm(x, meanlog = log(300), sdlog = 0.5),
                  xlab      = "D  (m²/day)",
                  xlim      = c(0, 600)),
    LAMBDA = list(samples   = draws[, "LAMBDA"],
                  prior_fun = function(x) dlnorm(x, meanlog = log(0.2), sdlog = 0.5),
                  xlab      = "LAMBDA  (day⁻¹)",
                  xlim      = c(0, 0.6)),
    GAMMA  = list(samples   = draws[, "GAMMA"],
                  prior_fun = function(x) dlnorm(x, meanlog = log(0.8), sdlog = 0.5),
                  xlab      = "GAMMA",
                  xlim      = c(0, 4))
  )

  pdf(path, width = 10, height = 4)
  on.exit(dev.off())

  # Reserve an outer top margin (oma) for the overall title so it sits
  # ABOVE the individual panel titles instead of overlapping them.
  old_par <- par(mfrow = c(1, 3),
                 mar = c(4, 4, 2.5, 1),
                 oma = c(0, 0, 3, 0))
  on.exit(par(old_par), add = TRUE)

  for (nm in names(params)) {
    p      <- params[[nm]]
    dens   <- density(p$samples)
    x_grid <- seq(p$xlim[1], p$xlim[2], length.out = 500)
    prior  <- p$prior_fun(x_grid)

    ylim <- c(0, max(dens$y, prior) * 1.1)

    plot(dens,
         main = nm,
         xlab = p$xlab,
         ylab = "Density",
         xlim = p$xlim,
         ylim = ylim,
         col  = "steelblue",
         lwd  = 2)

    lines(x_grid, prior, col = "tomato", lwd = 2, lty = 2)

    legend("topright",
           legend = c("Posterior", "Prior"),
           col    = c("steelblue", "tomato"),
           lwd    = 2,
           lty    = c(1, 2),
           bty    = "n")
  }

  mtext(paste0("R_PIEGE fixed at ", r_piege, " m"),
        side = 3, line = 0.5, outer = TRUE, cex = 1, font = 2)

  message("Plot saved to: ", path)
  invisible(path)
}
