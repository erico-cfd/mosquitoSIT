# Internal R mirror of the Stan solve_edp() function.
# Used only for unit testing — not exported to users.
# Mirrors the finite-difference scheme in inst/stan/mosquito_diffusion.stan.
.run_diffusion_r <- function(
    D, LAMBDA, GAMMA, R_PIEGE,
    n_initial = 50000,
    X_0 = 0, Y_0 = 0, R0 = 5,
    N = 61, T_days = 20, dt = 0.08,
    trap_positions = matrix(c(0, 2000), nrow = 1)
) {
  L_MIN <- -400; L_MAX <- 400
  dx  <- (L_MAX - L_MIN) / (N - 1)
  dx2 <- dx^2
  xs  <- seq(L_MIN, L_MAX, length.out = N)
  P   <- nrow(trap_positions)
  steps          <- as.integer(T_days / dt)
  steps_per_day  <- steps %/% T_days

  # --- Initial Gaussian density, scaled to n_initial mosquitoes ---
  h <- outer(xs - X_0, xs - Y_0, function(a, b) exp(-(a^2 + b^2) / R0^2))
  h <- h * (n_initial / (sum(h) * dx2))

  # --- Pre-compute capture kernels (same rescaling as Stan) ---
  analytical_integral <- GAMMA * pi * R_PIEGE^2
  F_TRAP <- array(0, dim = c(P, N, N))
  for (p in 1:P) {
    for (i in 1:N) {
      for (j in 1:N) {
        d2 <- (xs[i] - trap_positions[p, 1])^2 + (xs[j] - trap_positions[p, 2])^2
        if (d2 <= 10000) F_TRAP[p, i, j] <- GAMMA * exp(-d2 / R_PIEGE^2)
      }
    }
    num_integral <- sum(F_TRAP[p, , ]) * dx2
    if (num_integral > 1e-10)
      F_TRAP[p, , ] <- F_TRAP[p, , ] * (analytical_integral / num_integral)
  }
  V_TRAP <- apply(F_TRAP, c(2, 3), sum)

  # --- Time loop ---
  captures <- matrix(0, P, T_days)
  for (n in seq_len(steps)) {
    jour <- (n - 1L) %/% steps_per_day + 1L
    if (jour > T_days) jour <- T_days

    for (p in 1:P)
      captures[p, jour] <- captures[p, jour] + sum(h * F_TRAP[p, , ]) * dt * dx2

    # Vectorised 2D Laplacian (interior nodes only)
    lap <- h[1:(N-2), 2:(N-1)] + h[3:N,     2:(N-1)] +
           h[2:(N-1), 1:(N-2)] + h[2:(N-1), 3:N    ] -
           4 * h[2:(N-1), 2:(N-1)]

    h_next <- h
    h_next[2:(N-1), 2:(N-1)] <- pmax(0,
      h[2:(N-1), 2:(N-1)] +
      dt * (D * lap / dx2 - (LAMBDA + V_TRAP[2:(N-1), 2:(N-1)]) * h[2:(N-1), 2:(N-1)])
    )
    h <- h_next
  }

  list(h_final = h, captures = captures, dx = dx, xs = xs)
}
