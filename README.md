# clonecensorweighting

`clonecensorweighting` is an R package for building reproducible
clone-censor-weighting workflows for target trial emulation.

This repository is meant to be easy to use from a fresh GitHub clone and easy
to maintain as a shared collaboration project. The recommended setup uses:

- `rig` to install and switch to the project R version
- `renv` to restore the same package environment for every collaborator
- GitHub Actions to check reproducibility and package health

## Quick start

### 1. Clone the repository

```sh
git clone https://github.com/CausalInferenceLab/clonecensorweighting.git
cd clonecensorweighting
```

### 2. Use `rig` to install R 4.4.2

This project is pinned to R `4.4.2` for reproducibility.

```sh
rig add 4.4.2
rig default 4.4.2
```

If you already have R `4.4.2`, you can skip `rig add 4.4.2`.

### 3. Restore the project package library with `renv`

Open R in the project directory and run:

```r
install.packages("renv")
renv::restore()
```

`renv::restore()` installs the package versions recorded in `renv.lock`, so
everyone works with the same dependency set.

Note: this repository includes a `.Rprofile` file that activates `renv`
automatically when you open the project in R.

### 4. Install the package locally

From the project root:

```sh
R CMD INSTALL .
```

Then in R:

```r
library(clonecensorweighting)
```

## First example

```r
library(clonecensorweighting)

trial_data <- tibble::tibble(
  id = c(1, 2),
  follow_up = c(10, 12),
  event = c(1, 0),
  treatment = c("A", "B")
)

clones <- clone_censor_weighting(
  data = trial_data,
  id = "id",
  follow_up = "follow_up",
  event = "event",
  treatment = "treatment",
  regimes = c("A", "B")
)

clones

make_surv_response(
  data = trial_data,
  follow_up = "follow_up",
  event = "event"
)
```

## What the package currently provides

The package is still an early, lightweight foundation for future work. Right
now it includes:

- `read_trial_data()` to read trial-style CSV data into a tibble
- `clone_censor_weighting()` to create a starter cloned dataset across regimes
- `make_surv_response()` to build a `survival::Surv()` response object

## Working together on this repository

For collaborative work, the safest pattern is:

1. Clone the repository.
2. Switch to R `4.4.2` with `rig`.
3. Run `renv::restore()`.
4. Make your changes in a branch.
5. Run checks before opening a pull request.

Useful commands:

```sh
R CMD INSTALL .
R CMD check --no-manual .
```

In R:

```r
testthat::test_local()
```

If you add, remove, or upgrade dependencies, update the lockfile from R:

```r
renv::settings$snapshot.type("explicit")
renv::status(dev = TRUE)
renv::snapshot(dev = TRUE)
```

Please commit both code changes and the updated `renv.lock` when dependency
changes are intentional.

## Continuous integration

This repository includes two GitHub Actions workflows:

- `check-reproducible.yaml` runs `R CMD check` with pinned R `4.4.2` and `renv`
- `check-latest.yaml` runs broader checks across operating systems and R versions

Together, these workflows help keep the project reproducible for collaborators
and stable for future users.
