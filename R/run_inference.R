#' Prepare Stan data list for the mosquito diffusion-capture model
#'
#' @param captures Integer matrix of shape `[T_days, P]` with daily capture counts.
#' @param trap_positions Numeric matrix of shape `[P, 2]` with (x, y) trap coordinates in metres.
#' @param N_grid Number of grid points per axis (grid is `N_grid x N_grid`). Default 61.
#' @param dt Time step in days. Default 0.08.
#' @param X_0,Y_0 Initial release position in metres. Default `(0, 0)`.
#' @param R0 Radius (m) of the initial Gaussian mosquito density. Default 5.
#' @param n_initial Total initial mosquito count. Default 50 000.
#'
#' @return A named list ready to pass to `cmdstanr`.
#' @export
prepare_stan_data <- function(
    captures,
    trap_positions,
    N_grid    = 61L,
    dt        = 0.08,
    X_0       = 0.0,
    Y_0       = 0.0,
    R0        = 5.0,
    n_initial = 50000
) {
  stopifnot(is.matrix(captures), is.integer(captures))
  stopifnot(is.matrix(trap_positions), ncol(trap_positions) == 2)

  T_days <- nrow(captures)
  P      <- ncol(captures)
  steps  <- as.integer(T_days / dt)

  if (nrow(trap_positions) != P) {
    stop("`trap_positions` must have one row per trap (", P, " traps detected in `captures`).")
  }

  list(
    T          = T_days,
    P          = P,
    CAPTURES   = captures,
    POS_PIEGES = trap_positions,
    N          = as.integer(N_grid),
    steps      = steps,
    dt         = dt,
    X_0        = X_0,
    Y_0        = Y_0,
    R0         = R0,
    n_initial  = n_initial
  )
}


#' Run Bayesian inference for the mosquito diffusion-capture model
#'
#' Compiles and samples from the Stan PDE model that infers the diffusion
#' coefficient (`D`), mortality rate (`LAMBDA`), and trap capture rate
#' (`GAMMA`) from daily mosquito trap data.
#'
#' The Stan model solves a 2D diffusion-reaction PDE on a square grid and
#' computes expected captures per trap per day. Observed counts are modelled
#' as Poisson with these expected values as rates.
#'
#' @param captures Integer matrix `[T_days x P]` of daily capture counts.
#'   Defaults to the built-in `mosquito_captures` dataset (Bouyer et al.).
#' @param trap_positions Numeric matrix `[P x 2]` of (x, y) trap coordinates
#'   in metres. Defaults to `mosquito_trap_positions`.
#' @param N_grid Grid resolution (number of nodes per axis). Default 61.
#' @param dt Time step in days for the PDE solver. Default 0.08.
#' @param X_0,Y_0 Release point coordinates (m). Default `(0, 0)`.
#' @param R0 Initial Gaussian release radius (m). Default 5.
#' @param n_initial Total initial mosquito count. Default 50 000.
#' @param chains Number of MCMC chains. Default 1.
#' @param parallel_chains Chains to run in parallel. Default 1.
#' @param iter_warmup Warmup iterations per chain. Default 150.
#' @param iter_sampling Sampling iterations per chain. Default 150.
#' @param adapt_delta Target acceptance rate. Default 0.90.
#' @param max_treedepth Maximum NUTS tree depth. Default 10.
#' @param output_dir Directory for Stan CSV output files. Default `"."`.
#' @param force_recompile Re-compile the Stan model even if cached. Default `FALSE`.
#'
#' @return A `CmdStanMCMC` object. Key parameters to inspect:
#'   `"D"` (m²/day), `"LAMBDA"` (day⁻¹), `"GAMMA"`, `"R_PIEGE"` (m).
#'
#' @examples
#' \dontrun{
#' fit <- run_mosquito_inference()
#' fit$summary(c("D", "LAMBDA", "GAMMA", "R_PIEGE"))
#' fit$cmdstan_diagnose()
#' }
#'
#' @export
run_mosquito_inference <- function(
    captures        = NULL,
    trap_positions  = NULL,
    N_grid          = 61L,
    dt              = 0.08,
    X_0             = 0.0,
    Y_0             = 0.0,
    R0              = 5.0,
    n_initial       = 50000,
    chains          = 1L,
    parallel_chains = 1L,
    iter_warmup     = 150L,
    iter_sampling   = 150L,
    adapt_delta     = 0.90,
    max_treedepth   = 10L,
    output_dir      = ".",
    force_recompile = FALSE
) {
  if (is.null(captures))       captures       <- get("mosquito_captures",        envir = asNamespace("mosquitoSIT"))
  if (is.null(trap_positions)) trap_positions <- get("mosquito_trap_positions",  envir = asNamespace("mosquitoSIT"))

  stan_file <- system.file("stan", "mosquito_diffusion.stan", package = "mosquitoSIT")
  if (!nzchar(stan_file)) stop("Stan model file not found inside the mosquitoSIT package.")

  mod <- cmdstanr::cmdstan_model(stan_file, force_recompile = force_recompile)

  T_days <- nrow(captures)
  steps  <- as.integer(T_days / dt)
  dx     <- 800 / (N_grid - 1)

  message(sprintf("Grid: N=%d, dx=%.1f m, dt=%.3f days, steps=%d", N_grid, dx, dt, steps))
  message(sprintf("CFL check: dt * 4 * D_max / dx^2 = %.3f  (must be < 1 for stability)",
                  dt * 4 * 500 / dx^2))

  stan_data <- prepare_stan_data(
    captures       = captures,
    trap_positions = trap_positions,
    N_grid         = N_grid,
    dt             = dt,
    X_0            = X_0,
    Y_0            = Y_0,
    R0             = R0,
    n_initial      = n_initial
  )

  mod$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    output_dir      = output_dir
  )
}
