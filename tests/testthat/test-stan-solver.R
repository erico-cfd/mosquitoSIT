# This test calls the Stan solver directly (no sampling) using fixed_param=TRUE.
# It verifies the no-trap analytical case: when GAMMA=0, captures must be zero.

test_that("Stan solver: GAMMA=0 produces zero captures", {
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
          message = "CmdStan not installed")

  stan_file <- system.file("stan", "mosquito_diffusion.stan", package = "mosquitoSIT")
  mod <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)

  d <- prepare_stan_data(
    captures       = mosquito_captures,
    trap_positions = mosquito_trap_positions,
    r_piege        = 3.5
  )

  fit <- mod$sample(
    data          = d,
    fixed_param   = TRUE,
    iter_sampling = 1,
    chains        = 1,
    inits         = list(list(
      log_D      = log(200),
      log_LAMBDA = log(0.2),
      log_GAMMA  = log(0.8)
    )),
    refresh = 0
  )

  captures_no_trap <- as.vector(fit$draws("captures_no_trap"))
  expect_true(all(abs(captures_no_trap) < 1e-10),
              label = "all captures_no_trap values should be zero when GAMMA=0")
})
