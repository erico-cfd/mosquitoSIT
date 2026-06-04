# =============================================================================
# Irregular cumulative observations
# =============================================================================
# Real collectors do not empty every trap every day. What gets recorded is the
# CUMULATIVE number of captured mosquitoes on the days a trap was actually
# visited, and each trap can have its own (irregular) visit schedule.
#
# Example for one trap over 11 days (X = not collected that day):
#
#   day        1  2  3  4  5  6  7  8  9 10 11
#   cumulative 0  0  1  X  X  X 10 15  X  X 20
#
# This is encoded as a [T x P] matrix of cumulative counts with NA for the
# uncollected days. The functions below turn that into the interval records
# the Stan model expects (one Poisson likelihood term per collection event).
# =============================================================================


#' Build irregular observation records from a cumulative matrix
#'
#' Converts a `[T x P]` matrix of cumulative capture counts (with `NA` on the
#' days a trap was not collected) into the flat interval representation used by
#' the Stan model.
#'
#' For each trap, every observed (non-`NA`) day `d` produces one record: the
#' interval spans from the day after the previous observation up to `d`, and
#' the count is the increase in the cumulative total over that interval.
#'
#' @param cumulative Numeric/integer matrix `[T x P]` of cumulative counts,
#'   with `NA` where the trap was not collected. Cumulative values must be
#'   non-decreasing within each trap.
#'
#' @return A list with `N_obs`, `obs_trap`, `obs_start`, `obs_end`, `obs_count`.
#' @export
build_irregular_obs <- function(cumulative) {
  stopifnot(is.matrix(cumulative))
  T_days <- nrow(cumulative)
  P      <- ncol(cumulative)

  obs_trap  <- integer(0)
  obs_start <- integer(0)
  obs_end   <- integer(0)
  obs_count <- integer(0)

  for (p in seq_len(P)) {
    col      <- cumulative[, p]
    obs_days <- which(!is.na(col))
    if (length(obs_days) == 0) next

    prev_day <- 0L   # day before the experiment started
    prev_cum <- 0    # cumulative captures before the start are 0

    for (d in obs_days) {
      inc <- col[d] - prev_cum
      if (inc < 0) {
        stop("Cumulative counts must be non-decreasing (trap ", p,
             ", day ", d, ": ", col[d], " < previous ", prev_cum, ").")
      }
      obs_trap  <- c(obs_trap,  p)
      obs_start <- c(obs_start, prev_day + 1L)
      obs_end   <- c(obs_end,   d)
      obs_count <- c(obs_count, inc)

      prev_day <- d
      prev_cum <- col[d]
    }
  }

  list(
    N_obs     = length(obs_count),
    obs_trap  = as.integer(obs_trap),
    obs_start = as.integer(obs_start),
    obs_end   = as.integer(obs_end),
    obs_count = as.integer(round(obs_count))
  )
}


#' Prepare Stan data from irregular cumulative observations
#'
#' @param cumulative Numeric matrix `[T x P]` of cumulative captures with `NA`
#'   on uncollected days (see [build_irregular_obs()]).
#' @param trap_positions Numeric matrix `[P x 2]` of (x, y) trap coordinates (m).
#' @param N_grid Grid points per axis. Default 61.
#' @param dt Time step (days). Default 0.08.
#' @param X_0,Y_0 Release position (m). Default `(0, 0)`.
#' @param R0 Initial Gaussian release radius (m). Default 5.
#' @param n_initial Total initial mosquito count. Default 50000.
#' @param r_piege Fixed trap radius (m). Default 3.5.
#'
#' @return A named list ready to pass to `cmdstanr`.
#' @export
prepare_stan_data_irregular <- function(
    cumulative,
    trap_positions,
    N_grid    = 61L,
    dt        = 0.08,
    X_0       = 0.0,
    Y_0       = 0.0,
    R0        = 5.0,
    n_initial = 50000,
    r_piege   = 3.5
) {
  stopifnot(is.matrix(cumulative))
  stopifnot(is.matrix(trap_positions), ncol(trap_positions) == 2)

  T_days <- nrow(cumulative)
  P      <- ncol(cumulative)

  if (nrow(trap_positions) != P) {
    stop("`trap_positions` must have one row per trap (", P,
         " traps detected in `cumulative`).")
  }

  obs <- build_irregular_obs(cumulative)

  .assemble_stan_data(
    obs            = obs,
    T_days         = T_days,
    trap_positions = trap_positions,
    N_grid         = N_grid,
    dt             = dt,
    X_0            = X_0,
    Y_0            = Y_0,
    R0             = R0,
    n_initial      = n_initial,
    r_piege        = r_piege
  )
}


#' Run inference with irregular cumulative observations
#'
#' Same diffusion-capture model as [run_mosquito_inference()], but the
#' likelihood is built from cumulative counts collected on arbitrary,
#' per-trap schedules instead of dense daily counts.
#'
#' @param cumulative Numeric matrix `[T x P]` of cumulative captures with `NA`
#'   on uncollected days.
#' @param trap_positions Numeric matrix `[P x 2]` of trap coordinates (m).
#'   Defaults to the built-in `mosquito_trap_positions`.
#' @param N_grid,dt,X_0,Y_0,R0,n_initial,r_piege Model/grid settings, passed to
#'   [prepare_stan_data_irregular()].
#' @param chains,parallel_chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth,output_dir
#'   Sampler settings, passed to `cmdstanr`.
#' @param force_recompile Re-compile the Stan model. Default `FALSE`.
#'
#' @return A `CmdStanMCMC` object.
#' @export
run_mosquito_inference_irregular <- function(
    cumulative,
    trap_positions  = NULL,
    N_grid          = 61L,
    dt              = 0.08,
    X_0             = 0.0,
    Y_0             = 0.0,
    R0              = 5.0,
    n_initial       = 50000,
    r_piege         = 3.5,
    chains          = 1L,
    parallel_chains = 1L,
    iter_warmup     = 150L,
    iter_sampling   = 150L,
    adapt_delta     = 0.90,
    max_treedepth   = 10L,
    output_dir      = ".",
    force_recompile = FALSE
) {
  if (is.null(trap_positions))
    trap_positions <- get("mosquito_trap_positions", envir = asNamespace("mosquitoSIT"))

  stan_file <- system.file("stan", "mosquito_diffusion.stan", package = "mosquitoSIT")
  if (!nzchar(stan_file)) stop("Stan model file not found inside the mosquitoSIT package.")

  mod <- cmdstanr::cmdstan_model(stan_file, force_recompile = force_recompile)

  stan_data <- prepare_stan_data_irregular(
    cumulative     = cumulative,
    trap_positions = trap_positions,
    N_grid         = N_grid,
    dt             = dt,
    X_0            = X_0,
    Y_0            = Y_0,
    R0             = R0,
    n_initial      = n_initial,
    r_piege        = r_piege
  )

  dx <- 800 / (N_grid - 1)
  message(sprintf("Grid: N=%d, dx=%.1f m, dt=%.3f days, steps=%d",
                  N_grid, dx, dt, stan_data$steps))
  message(sprintf("Irregular observations: %d collection events across %d traps over %d days",
                  stan_data$N_obs, stan_data$P, stan_data$T))

  fit <- mod$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    output_dir      = output_dir
  )

  print_inference_summary(fit, r_piege = r_piege)
  fit
}


#' Build a demo irregular dataset by punching holes in daily captures
#'
#' Takes a dense daily-capture matrix, turns it into a cumulative matrix, then
#' randomly drops a fraction of the collection days per trap (replacing them
#' with `NA`). The last day is always kept so the final total is observed.
#' Useful for testing whether the irregular-interval inference recovers similar
#' parameters to the fully-observed case.
#'
#' @param daily Integer matrix `[T x P]` of daily captures. Defaults to the
#'   built-in `mosquito_captures`.
#' @param keep_prob Probability that any given intermediate day is collected.
#'   Default 0.5.
#' @param seed Optional integer for reproducibility.
#'
#' @return A `[T x P]` cumulative matrix with `NA` on the dropped days.
#' @export
make_irregular_demo <- function(daily = NULL, keep_prob = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(daily))
    daily <- get("mosquito_captures", envir = asNamespace("mosquitoSIT"))

  T_days <- nrow(daily)
  P      <- ncol(daily)

  cumulative <- apply(daily, 2, cumsum)
  storage.mode(cumulative) <- "double"

  for (p in seq_len(P)) {
    for (t in seq_len(T_days - 1L)) {        # never drop the final day
      if (stats::runif(1) > keep_prob) cumulative[t, p] <- NA
    }
  }

  cumulative
}
