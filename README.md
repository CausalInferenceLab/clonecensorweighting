# clonecensorweighting

Infrastructure for building a CRAN-ready R package around clone-censor-weighting
workflows for target trial emulation.

## Local development

Pin local development to R 4.4.2 with `rig`:

```sh
rig add 4.4.2
rig default 4.4.2
```

Then install `renv` once and restore the project library:

```r
install.packages("renv")
renv::restore()
```

If you need to refresh the lockfile after changing package dependencies:

```r
renv::settings$snapshot.type("explicit")
renv::status(dev = TRUE)
renv::snapshot(dev = TRUE)
```

## CI

The repository includes two workflows:

- `check-reproducible.yaml` runs `R CMD check` against pinned `R 4.4.2` with `renv`.
- `check-latest.yaml` runs broader cross-platform checks against release, devel, and oldrel.
