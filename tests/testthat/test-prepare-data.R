test_that("prepare_stan_data passes r_piege through correctly", {
  d <- prepare_stan_data(
    captures       = mosquito_captures,
    trap_positions = mosquito_trap_positions,
    r_piege        = 3.5
  )
  expect_equal(d$R_PIEGE, 3.5)
})

test_that("prepare_stan_data errors when captures is not integer", {
  bad_captures <- matrix(1.5, nrow = 20, ncol = 21)
  expect_error(
    prepare_stan_data(bad_captures, mosquito_trap_positions),
    regexp = "integer"
  )
})

test_that("prepare_stan_data errors when trap count mismatches captures columns", {
  wrong_traps <- matrix(c(0, 0), nrow = 1)
  expect_error(
    prepare_stan_data(mosquito_captures, wrong_traps),
    regexp = "21"
  )
})
