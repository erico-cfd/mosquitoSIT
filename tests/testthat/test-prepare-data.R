test_that("prepare_stan_data returns correct dimensions with built-in data", {
  d <- prepare_stan_data(
    captures       = mosquito_captures,
    trap_positions = mosquito_trap_positions
  )
  expect_equal(d$T, 20L)
  expect_equal(d$P, 21L)
  expect_equal(d$N, 61L)
  expect_equal(d$steps, as.integer(20 / 0.08))
})

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
  wrong_traps <- matrix(c(0, 0), nrow = 1)   # 1 trap, but captures has 21
  expect_error(
    prepare_stan_data(mosquito_captures, wrong_traps),
    regexp = "21"
  )
})

test_that("built-in data has expected dimensions", {
  expect_equal(dim(mosquito_captures),        c(20L, 21L))
  expect_equal(dim(mosquito_trap_positions),  c(21L, 2L))
  expect_true(is.integer(mosquito_captures))
  expect_true(all(mosquito_captures >= 0L))
})
