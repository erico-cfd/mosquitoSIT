// =============================================================================
// MOSQUITO DIFFUSION-CAPTURE MODEL v5 — RESOLUTION-INDEPENDENT
// =============================================================================
// The capture Gaussian has radius R ~ 3.5m. On a coarse grid (dx ~ 27m),
// the numerical integral of the Gaussian severely underestimates or 
// overestimates the true integral (depending on trap-to-node alignment).
//
// Fix: after computing F_TRAP[p] = GAMMA * exp(-d²/R²), we RESCALE it so
// that its numerical integral matches the analytical one (GAMMA * pi * R²).
//
// - Coarse grid (dx >> R): correction compensates for under-sampling → correct total
// - Fine grid (dx << R): correction ≈ 1 → natural Gaussian behavior
// =============================================================================

functions {  // This is the edp solver
    matrix solve_edp(real D, real GAMMA, real LAMBDA, real R_PIEGE, 
                     int steps, int T, int N, int P, real dt, real dx, 
                     matrix h_initial, array[,,] real DIST_CARRE,
                     array[,,] int TRAP_MASK) {
        
        matrix[P, T] CAPTURES_SIM = rep_matrix(0, P, T);
        matrix[N, N] h = h_initial; 
        
        real dx2 = square(dx);
        real inv_R2 = 1.0 / square(R_PIEGE);
        int steps_per_day = steps %/% T;

        // ==========================================================
        // PRE-COMPUTE capture kernels with integral correction
        // ==========================================================
        array[P] matrix[N, N] F_TRAP;
        matrix[N, N] V_TRAP = rep_matrix(0, N, N);
        
        // Analytical integral of the Gaussian capture function:
        // integral of GAMMA * exp(-r²/R²) over R² = GAMMA * pi * R²
        real analytical_integral = GAMMA * pi() * square(R_PIEGE);
        
        for (p in 1:P) {
            F_TRAP[p] = rep_matrix(0, N, N);
            real num_sum = 0;
            
            // Step 1: compute raw Gaussian values
            for (i in 1:N) {
                for (j in 1:N) {
                    if (TRAP_MASK[p, i, j] == 1) {
                        real fp = GAMMA * exp(-DIST_CARRE[p, i, j] * inv_R2);
                        F_TRAP[p][i, j] = fp;
                        num_sum += fp;
                    }
                }
            }
            
            // Step 2: compute numerical integral and correction factor
            real num_integral = num_sum * dx2;
            
            // Correction: scale so that sum(F_TRAP[p]) * dx² = analytical_integral
            if (num_integral > 1e-10) {
                real correction = analytical_integral / num_integral;
                F_TRAP[p] = F_TRAP[p] * correction;
            }
            
            // Accumulate total capture rate
            V_TRAP += F_TRAP[p];
        }

        // ==========================================================
        // TIME LOOP
        // ==========================================================
        for (n in 1:steps) { 
            int JOUR = (n - 1) %/% steps_per_day + 1;
            if (JOUR > T) JOUR = T;
            
            // Capture accumulation (vectorized per trap)
            for (p in 1:P) {
                CAPTURES_SIM[p, JOUR] += sum(h .* F_TRAP[p]) * dt * dx2;
            }
            
            // PDE update
            matrix[N, N] h_next = h;
            for (i in 2:(N - 1)) {
                for (j in 2:(N - 1)) {
                    real LAP = h[i-1, j] + h[i+1, j] 
                             + h[i, j-1] + h[i, j+1] 
                             - 4 * h[i, j];
                    real reaction = (LAMBDA + V_TRAP[i, j]) * h[i, j];
                    h_next[i, j] = fmax(0.0, h[i, j] + dt * (D * LAP / dx2 - reaction));
                }
            }
            h = h_next;
        }
        return CAPTURES_SIM;
    }
}

data {
    int<lower=1> T;                          // total number of simulated days
    int<lower=1> P;                          // number of traps
    array[P] vector[2] POS_PIEGES;

    // --- Irregular cumulative observations ---------------------------------
    // Each record k is one collection event: over the days
    // obs_start[k] .. obs_end[k] (inclusive), trap obs_trap[k] accumulated
    // obs_count[k] new mosquitoes. This supports arbitrary, per-trap
    // observation schedules ("holes" in the data). Daily observation is just
    // the special case where every interval has length 1.
    int<lower=1> N_obs;
    array[N_obs] int<lower=1> obs_trap;
    array[N_obs] int<lower=1> obs_start;
    array[N_obs] int<lower=1> obs_end;
    array[N_obs] int<lower=0> obs_count;

    int<lower=2> N;
    int<lower=1> steps;
    real<lower=0> dt;

    real X_0;
    real Y_0;
    real<lower=0> R0;
    real<lower=1> n_initial;

    // R_PIEGE is fixed by the user instead of being derived from GAMMA.
    // This removes the perfect collinearity between GAMMA and R_PIEGE
    // (they shared only one degree of freedom). A typical value is 1.0 m.
    real<lower=0> R_PIEGE;
}

transformed data { 
    real L_MIN = -400;
    real L_MAX = 400;
    real dx = (L_MAX - L_MIN) / (N - 1);

    real SEUIL_DIST_CARRE = 10000.0;

    array[P, N, N] real DIST_CARRE;
    array[P, N, N] int TRAP_MASK;
    matrix[N, N] h_initial = rep_matrix(0, N, N);

    for (p in 1:P) {
        for (i in 1:N) {
            for (j in 1:N) {
                real x = L_MIN + (i - 1) * dx;
                real y = L_MIN + (j - 1) * dx;
                real d2 = square(x - POS_PIEGES[p, 1]) + square(y - POS_PIEGES[p, 2]);
                DIST_CARRE[p, i, j] = d2;
                TRAP_MASK[p, i, j] = (d2 <= SEUIL_DIST_CARRE) ? 1 : 0;
            }
        }
    }

    {
        real total_mass = 0;
        for (i in 1:N) {
            for (j in 1:N) {
                real x = L_MIN + (i - 1) * dx;
                real y = L_MIN + (j - 1) * dx;
                h_initial[i, j] = exp(-(square(x - X_0) + square(y - Y_0)) / square(R0));
                total_mass += h_initial[i, j];
            }
        }
        real scale = n_initial / (total_mass * square(dx));
        for (i in 1:N) {
            for (j in 1:N) {
                h_initial[i, j] *= scale;
            }
        }
    }
}

parameters {
    // Sampling on the log scale: unconstrained reals, positivity guaranteed
    // by exponentiation. No hard bounds needed.
    real log_D;
    real log_LAMBDA;
    real log_GAMMA;
}

transformed parameters {
    real D      = exp(log_D);
    real LAMBDA = exp(log_LAMBDA);
    real GAMMA  = exp(log_GAMMA);

    matrix[P, T] CAPTURES_SIM = solve_edp(D, GAMMA, LAMBDA, R_PIEGE,
                                           steps, T, N, P, dt, dx,
                                           h_initial, DIST_CARRE,
                                           TRAP_MASK);
}

model {
    // Priors on the log scale (equivalent to log-normal priors on D, LAMBDA, GAMMA).
    // Centered at log of our prior beliefs: D~300, LAMBDA~0.2, GAMMA~0.8
    log_D      ~ normal(log(300), 0.5);
    log_LAMBDA ~ normal(log(0.2), 0.5);
    log_GAMMA  ~ normal(log(0.8), 0.5);

    // Likelihood over irregular cumulative intervals.
    // The expected count over an interval is the sum of the daily expected
    // captures inside it. A sum of independent Poisson counts is Poisson with
    // the summed rate, so the cumulative increment is Poisson(sum of rates).
    for (k in 1:N_obs) {
        real rate = 0;
        for (t in obs_start[k]:obs_end[k]) {
            rate += CAPTURES_SIM[obs_trap[k], t];
        }
        obs_count[k] ~ poisson(fmax(1e-6, rate));
    }
}
