# Mizer

Mizer is an R package for dynamic multi-species size-spectrum modelling of
fish communities. It tracks the full size distribution of each species and
the plankton resource, computing growth, predation, and mortality from
individual-level physiology.

## Core workflow

```r
library(mizer)

# 1. Create model parameters from a species data frame
params <- newMultispeciesParams(species_params, interaction)

# 2. Find the steady state (sets initial values)
params <- steady(params)

# 3. Calibrate to observed biomasses / yields
params <- calibrateBiomass(params)  # adjusts kappa
params <- matchBiomasses(params)    # adjusts R_max per species
params <- matchGrowth(params)       # adjusts h per species

# 4. Tune density-dependent reproduction
params <- setBevertonHolt(params, reproduction_level = 0.25)

# 5. Project forward in time
sim <- project(params, t_max = 20, effort = 1)

# 6. Analyse results
plot(sim)
getBiomass(sim)
getYield(sim)
plotSpectra(sim)
```

## Key objects

**`MizerParams`** — holds all model parameters. Never modify slots directly.
All setter functions return a new copy: `params <- setFishing(params, ...)`.
Change species parameters with `given_species_params(params) <- value`, which
triggers recalculation of dependent quantities.

**`MizerSim`** — simulation output from `project()`. Arrays: `N(sim)` (time ×
species × size), `NResource(sim)`.

## Species parameters

The `species_params` data frame must have `species` (name) and `w_max`
(maximum weight in grams). Everything else has defaults.

| Column | Meaning |
|--------|---------|
| `w_mat` | Maturity weight (g) |
| `beta` | Preferred predator/prey mass ratio (default ~100) |
| `sigma` | S.d. of lognormal predation kernel (default ~1.3) |
| `h` | Max intake rate coefficient |
| `alpha` | Assimilation efficiency (default 0.6) |
| `erepro` | Reproductive efficiency |
| `R_max` | Beverton-Holt max recruitment |
| `biomass_observed` | Observed biomass for `calibrateBiomass()` |
| `yield_observed` | Observed yield for `matchYields()` |

## Units

Weights in grams, lengths in cm, time in years.

## Plotting

Mizer provides many built-in plotting functions. Always prefer these over
writing custom plotting code.

```r
plot(sim)              # overview of simulation
plotSpectra(sim)       # size spectra
plotBiomass(sim)       # biomass over time
plotYield(sim)         # yield over time
plotGrowthCurves(sim)  # growth curves
plotFMort(sim)         # fishing mortality
```

The return values of most `get...()` functions also have `plot()` methods,
so you can visualise any quantity directly without writing custom plotting code:

```r
plot(getSSB(sim))           # ArrayTimeBySpecies  → time series per species
plot(getTrophicLevel(params)) # ArraySpeciesBySize → curve per species
```

Grep for "plot" in `llms-full.txt` to discover the full list of available
plots before writing any custom code. Grep for a specific function name to
look up its documentation — do not read the whole file.

## Extending mizer

To replace a rate function: `params <- setRateFunction(params, "Encounter", myFun)`.
To add a new ecosystem component: `params <- setComponent(params, "detritus", ...)`.
See https://sizespectrum.org/mizer/articles/extending-mizer.html

## Single-species model and resource dynamics

### `newSingleSpeciesParams()` uses a fixed resource by default

```r
params <- newSingleSpeciesParams(w_max = 1e4, kappa = 0.005, f0 = 0.3)
resource_dynamics(params)  # "resource_constant" -- food is truly fixed!
```

The fixed background is the *only* food source. To enable density-dependent
coupling (resource responds to fish predation), switch explicitly:

```r
params <- setResource(params, resource_rate = 0.1,
                      resource_dynamics = "resource_semichemostat")
```

### `setResource()` with `balance = TRUE` (default) does NOT change feeding level

`balance = TRUE` adjusts `c_R` so that the current `N_R` is exactly at its
semi-chemostat equilibrium. Changing `r_R` with balance left on preserves the
equilibrium feeding level and growth rates — the jam signal is therefore the
same for `r_R = 0.1` and `r_R = 10` when `f0` is fixed. Use different `f0`
values (not different `r_R` values) to get clearly distinct perturbation-decay
curves.

The resource level at equilibrium is `N_R/c_R = r_R/(r_R + mu_R)`. Low `r_R`
causes strong depletion at sizes where predation is heavy (small prey sizes
eaten by abundant small fish). You can inspect this with:

```r
resource_level(params)          # vector over w_full(params)
range(resource_level(params))   # [0.33, 1] with r_R=0.1, e.g.
```

### `resource_level(p) <- rl` adjusts `c_R`, NOT `N_R`

This is a critical pitfall: `getEGrowth()` uses the actual resource abundance
`N_R`, not the capacity `c_R`. Setting `resource_level(p) <- 0.5` only
re-labels what "100%" means; it leaves `N_R` and hence growth rates unchanged.
To actually deplete the food supply and observe its effect on `g(w)`:

```r
# WRONG -- growth rates unchanged:
resource_level(p) <- 0.8

# CORRECT -- actually reduces N_R:
initialNResource(p) <- initialNResource(p) * 0.8
```

### Sensitivity of g(w) near the critical feeding level

With `f0` close to `fc` (default `fc = 0.25`), growth collapses rapidly:

```
g(w) = alpha * h * w^n * (f(w) - fc)
```

After scaling `N_R` by `rl`: `f_new = rl*f0 / (rl*f0 + 1 - f0)`.
Growth goes to zero at `rl_min = (1-f0)*fc / (f0*(1-fc))`.
For `f0 = 0.3, fc = 0.25`: `rl_min ≈ 0.778` — a 22% resource drop halts all
growth. Use depletion values **above** this threshold to show the collapse:

```r
# For f0=0.3, safe depletion range:
depletion_levels <- c(1.0, 0.95, 0.90, 0.85, 0.80)
```

### Resource indexing uses `w_full()`, not `w()`

`NResource(sim)` has dimensions `[n_times, length(w_full(params))]`.
Always use `w_full(params)` to map resource sizes to indices:

```r
wf       <- w_full(params)
prey_idx <- which.min(abs(wf - 0.5))   # index of w = 0.5 g on full grid
NResource(sim)[time_idx, prey_idx]
```

### Resource coupling is weak in the single-species model

With the default broad predation kernel (`sigma = 1.3`) and cannibalism off,
a localised bump in fish at size `w` causes negligible depletion at its
preferred prey size `w/beta`. The predation mortality on resource at `w/beta`
is tiny (`mu_R ~ 1e-7 yr^{-1}`) compared to `r_R`, so the resource there
stays near full capacity. The strong coupling is only at very small sizes
(w ~ 1e-5 g) where `mu_R >> r_R` and the many small fish eat heavily.

Consequence: the "traffic jam" signal in a single-species model is primarily
**kinematic** — perturbations decay slowly because `g(w) ~ f_0 - f_c` is tiny
near `f_c`, not because of self-reinforcing density-dependent resource
depletion. For a true **dynamic instability** (self-amplifying jam), use a
community model (`newCommunityParams()`) where fish eat each other.

### Template for perturbation experiments

```r
# add a Gaussian density bump at one size class
add_bump <- function(params, w_bump = 50, amplitude = 4, log10_sigma = 0.4) {
  N0  <- initialN(params); wv <- w(params)
  bump <- amplitude * exp(-((log10(wv) - log10(w_bump))^2) / (2*log10_sigma^2))
  N0[1, ] <- N0[1, ] * (1 + bump)
  initialN(params) <- N0
  params
}

# run baseline and perturbed pair
run_pair <- function(params, t_max = 30) {
  list(
    base = project(params,             t_max = t_max, effort = 0, progress_bar = FALSE),
    pert = project(add_bump(params),   t_max = t_max, effort = 0, progress_bar = FALSE)
  )
}

# extract jam signal: ratio N_perturbed / N_baseline over time
ratio_at_t <- function(sim_base, sim_pert, t_idx) {
  N(sim_pert)[t_idx, 1, ] / N(sim_base)[t_idx, 1, ]
}
```

### `getFeedingLevel()` return type

Returns `ArrayTimeBySpeciesBySize` with dimensions `[n_times, n_species, n_sizes]`.
Index as `fl[time_idx, species_idx, size_idx]`, not with named dimensions.

## API documentation (local copies)

Concise overview of the mizer API (start here):
/home/gustav/R/x86_64-pc-linux-gnu-library/4.6/mizerAgents/llms.txt

Full API documentation with complete details on every function:
/home/gustav/R/x86_64-pc-linux-gnu-library/4.6/mizerAgents/llms-full.txt
Grep or search this file for specific function names - do not read the whole file.

