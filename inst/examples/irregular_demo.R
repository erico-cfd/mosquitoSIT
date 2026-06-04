# =============================================================================
# Demo: inference with irregular cumulative observations ("holes" in the data)
# =============================================================================
# Run from R with:  source(system.file("examples/irregular_demo.R",
#                                       package = "mosquitoSIT"))
# =============================================================================

library(mosquitoSIT)

# -----------------------------------------------------------------------------
# Example A — a tiny hand-built dataset (the format you describe)
# -----------------------------------------------------------------------------
# Cumulative captures per day, NA where the trap was NOT collected.
# Each trap can have its own irregular schedule.
#
#   trap 1: 0 0 1 X X X 10 15 X X 20
#   trap 2: 0 X X X 5  5  7 X X X 12
NA_ <- NA
cumulative <- cbind(
  trap1 = c(0, 0, 1, NA_, NA_, NA_, 10, 15, NA_, NA_, 20),
  trap2 = c(0, NA_, NA_, NA_, 5, 5, 7, NA_, NA_, NA_, 12)
)

# Inspect how this is turned into collection intervals:
print(build_irregular_obs(cumulative))
# obs_start..obs_end is the day range each increment accumulated over.

# Two traps at arbitrary positions (metres)
traps2 <- matrix(c(0, 50,
                   50, 0), nrow = 2, byrow = TRUE)

# (Uncomment to actually run the sampler on this toy example)
# fit_toy <- run_mosquito_inference_irregular(cumulative, traps2, r_piege = 3.5)


# -----------------------------------------------------------------------------
# Example B — punch holes in the real Bouyer dataset and re-fit
# -----------------------------------------------------------------------------
# This is the key TEST: does the inference still recover similar parameters
# when only ~50% of collection days are available?

cumulative_holed <- make_irregular_demo(keep_prob = 0.5, seed = 1)

fit_irreg <- run_mosquito_inference_irregular(
  cumulative = cumulative_holed,
  r_piege    = 3.5
)

# Compare the posterior to the fully-observed run from main:
#   D ~ 164, LAMBDA ~ 0.219, GAMMA ~ 0.775
fit_irreg$summary(c("D", "LAMBDA", "GAMMA"))
