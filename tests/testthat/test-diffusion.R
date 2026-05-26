# These tests use .run_diffusion_r(), an internal R mirror of the Stan PDE
# solver. By testing against known analytical solutions we verify that the
# finite-difference scheme is implemented correctly.
#
# Analytical solution for 2D diffusion of an initial Gaussian:
#
#   h0(x,y) = A * exp(-(x²+y²) / R0²)
#
# At time t with pure diffusion (LAMBDA=0, no traps):
#
#   h(x,y,t) = A * R0² / (R0² + 4Dt) * exp(-(x²+y²) / (R0² + 4Dt))
#
# Two key properties we test:
#   1. Total mass is conserved  (no mosquitoes created or destroyed)
#   2. Spatial shape matches the formula above

# Place the single required trap far outside the domain so it captures nothing.
.far_trap <- matrix(c(0, 5000), nrow = 1)

test_that("pure diffusion conserves total mass", {
  # With LAMBDA=0 and an unreachable trap, the total number of mosquitoes
  # must stay equal to n_initial throughout the simulation.
  n_initial <- 50000

  result <- mosquitoSIT:::.run_diffusion_r(
    D             = 100,
    LAMBDA        = 0,
    GAMMA         = 0,
    R_PIEGE       = 1,
    n_initial     = n_initial,
    R0            = 30,
    N             = 61,
    T_days        = 5,
    dt            = 0.08,
    trap_positions = .far_trap
  )

  final_mass <- sum(result$h_final) * result$dx^2

  # Accept up to 1% error from boundary effects and floating-point accumulation
  expect_equal(final_mass, n_initial, tolerance = 0.01 * n_initial)
})


test_that("diffusion spreads as the analytical Gaussian solution", {
  # Parameters chosen so the Gaussian stays well within the grid boundaries.
  D         <- 50
  T_days    <- 5
  R0        <- 30      # initial radius large enough to be resolved on the grid
  n_initial <- 50000
  N         <- 61

  result <- mosquitoSIT:::.run_diffusion_r(
    D             = D,
    LAMBDA        = 0,
    GAMMA         = 0,
    R_PIEGE       = 1,
    n_initial     = n_initial,
    R0            = R0,
    N             = N,
    T_days        = T_days,
    dt            = 0.08,
    trap_positions = .far_trap
  )

  # Analytical solution: effective radius grows as R_eff² = R0² + 4*D*T
  R_eff_sq <- R0^2 + 4 * D * T_days
  A        <- n_initial / (pi * R0^2)           # amplitude of initial Gaussian
  A_t      <- A * R0^2 / R_eff_sq              # amplitude at time T

  xs <- result$xs
  h_analytical <- outer(xs, xs, function(a, b) A_t * exp(-(a^2 + b^2) / R_eff_sq))

  # The numerical and analytical fields should be highly correlated (shape match)
  cor_val <- cor(as.vector(result$h_final), as.vector(h_analytical))
  expect_gt(cor_val, 0.999)

  # The peak value at the centre (0,0) should match within 5%
  centre <- which.min(abs(xs))   # index closest to x=0
  expect_equal(
    result$h_final[centre, centre],
    h_analytical[centre, centre],
    tolerance = 0.05 * h_analytical[centre, centre]
  )
})


test_that("mortality reduces total mass over time", {
  # With LAMBDA > 0 and no traps, total mass must strictly decrease.
  result_alive <- mosquitoSIT:::.run_diffusion_r(
    D = 100, LAMBDA = 0,   GAMMA = 0, R_PIEGE = 1,
    n_initial = 50000, R0 = 30, T_days = 5, dt = 0.08,
    trap_positions = .far_trap
  )
  result_dying <- mosquitoSIT:::.run_diffusion_r(
    D = 100, LAMBDA = 0.3, GAMMA = 0, R_PIEGE = 1,
    n_initial = 50000, R0 = 30, T_days = 5, dt = 0.08,
    trap_positions = .far_trap
  )

  mass_alive <- sum(result_alive$h_final) * result_alive$dx^2
  mass_dying <- sum(result_dying$h_final) * result_dying$dx^2

  expect_gt(mass_alive, mass_dying)
})
