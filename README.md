# mosquitoSIT

Bayesian inference for mosquito dispersal models in the context of the **Sterile Insect Technique (SIT)**. Given daily trap capture counts, the package estimates three biological parameters using a 2D diffusion-reaction PDE solved inside a Stan model and sampled with Hamiltonian Monte Carlo (HMC/NUTS).

| Parameter | Meaning | Prior |
|-----------|---------|-------|
| `D` | Diffusion coefficient (m²/day) — how fast mosquitoes spread | Log-normal, median 300 |
| `LAMBDA` | Mortality rate (day⁻¹) — daily death probability | Log-normal, median 0.2 |
| `GAMMA` | Trap capture rate — peak attraction strength of a trap | Log-normal, median 0.8 |
| `R_PIEGE` | Trap effective radius (m) — **fixed by the user** | Not estimated |

All three parameters are sampled on the log scale (`log_D`, `log_LAMBDA`, `log_GAMMA`), which removes the need for hard boundaries and improves sampler efficiency. `R_PIEGE` is fixed as a known input rather than estimated, removing a collinearity issue with `GAMMA`.

---

## Installation

You need R (≥ 4.1) and a working Stan installation.

**Step 1 — Install Stan (once)**
```r
install.packages("cmdstanr", repos = "https://mc-stan.org/r-packages/")
cmdstanr::install_cmdstan()
```

**Step 2 — Install mosquitoSIT**
```r
devtools::install_github("erico-cfd/mosquitoSIT")
```

---

## Usage

### Run with built-in data (Bouyer et al.)

The package includes the trap positions and daily captures from the reference study (21 traps, 20 days).

```r
library(mosquitoSIT)

fit <- run_mosquito_inference(r_piege = 3.5)
```

Expected console output during sampling:
```
Grid: N=61, dx=13.3 m, dt=0.080 days, steps=250
CFL check: dt * 4 * D_max / dx^2 = 0.900  (must be < 1 for stability)
Compiling Stan program...
Running MCMC with 1 chain...
Chain 1 Iteration:   1 / 300 [  0%]  (Warmup)
...
Chain 1 Iteration: 300 / 300 [100%]  (Sampling)
Chain 1 finished in ...s.
```

### Inspect results

```r
fit$summary(c("D", "LAMBDA", "GAMMA"))
```

Expected output (values will vary):
```
  variable  mean  median    sd   mad     q5    q95  rhat ess_bulk ess_tail
  D        163.   163.   5.62  5.93  154.   172.   0.997     95.0     88.1
  LAMBDA     0.220  0.221  0.00735 0.00718  0.208  0.231  1.03      81.7     76.2
  GAMMA      0.783  0.781  0.0348  0.0328   0.731  0.845  1.01      65.6     71.3
```

- **`mean` / `median`** — posterior estimate of the parameter
- **`sd`** — uncertainty (standard deviation of the posterior)
- **`rhat`** — convergence diagnostic: should be close to 1.00 (< 1.01 is good)
- **`ess_bulk`** — effective sample size: higher is better (aim for > 100)

### Run diagnostics

```r
fit$cmdstan_diagnose()
```

This checks for divergences, low E-BFMI, and tree depth warnings. A healthy run looks like:

```
No divergent transitions found.
No saturated tree depths found.
E-BFMI satisfactory for all chains.
```

---

## Run with your own data

### Input format

| Argument | Type | Shape | Description |
|----------|------|-------|-------------|
| `captures` | integer matrix | `[T × P]` | Daily captures: rows = days, columns = traps |
| `trap_positions` | numeric matrix | `[P × 2]` | Trap (x, y) coordinates in metres from release point |
| `r_piege` | numeric | scalar | Fixed trap effective radius in metres |

```r
# Example: 10 days, 5 traps
my_captures <- matrix(as.integer(c(
  0, 0, 0, 2, 0,
  0, 1, 0, 5, 0,
  0, 3, 1, 8, 1,
  1, 4, 2, 6, 2,
  2, 5, 3, 4, 3,
  3, 4, 4, 3, 2,
  2, 3, 3, 2, 1,
  1, 2, 2, 1, 0,
  0, 1, 1, 0, 0,
  0, 0, 0, 0, 0
)), nrow = 10, ncol = 5, byrow = TRUE)

my_traps <- matrix(c(
    0,  100,
  100,    0,
    0, -100,
 -100,    0,
    0,    0
), nrow = 5, ncol = 2, byrow = TRUE)

fit <- run_mosquito_inference(
  captures       = my_captures,
  trap_positions = my_traps,
  r_piege        = 3.5
)
```

### Tuning the sampler

```r
fit <- run_mosquito_inference(
  r_piege         = 3.5,
  chains          = 4,     # run 4 independent chains (better convergence check)
  parallel_chains = 4,     # run them in parallel
  iter_warmup     = 500,   # more warmup for difficult posteriors
  iter_sampling   = 500,
  adapt_delta     = 0.95,  # increase if you get divergent transitions
  output_dir      = "results/"
)
```

### Prepare data without running

If you want to inspect or modify the Stan data list before sampling:

```r
stan_data <- prepare_stan_data(
  captures       = my_captures,
  trap_positions = my_traps,
  r_piege        = 3.5,
  N_grid         = 61,
  dt             = 0.08
)
str(stan_data)
```

---

## Grid and solver parameters

| Argument | Default | Effect |
|----------|---------|--------|
| `r_piege` | 1.0 | Fixed trap effective radius (m) — recommended: 3.5 |
| `N_grid` | 61 | Grid resolution (61 × 61 nodes over 800 m × 800 m) |
| `dt` | 0.08 | PDE time step in days — smaller = more accurate but slower |
| `X_0`, `Y_0` | 0, 0 | Release point coordinates (m) |
| `R0` | 5.0 | Initial Gaussian release radius (m) |
| `n_initial` | 50000 | Total initial mosquito count |

**CFL stability condition:** the solver requires `dt × 4 × D / dx² < 1`. With the defaults (`N_grid=61`, `dt=0.08`) this holds for `D` up to ~500 m²/day. If you coarsen the grid or increase `dt`, verify the CFL check printed at runtime.

---

## Reference

Model based on the methodology described in:

> Bouyer et al. — *Dispersal of Sterile Aedes albopictus males*, see [article](https://umr-astre.pages-forge.inrae.fr/sit-methods/article.pdf).
