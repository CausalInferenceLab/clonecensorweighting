# ------------------------------------------------------------
# Clone-Censor-Weight simulated example
# Post-MI antiplatelet initiation window: start within 30 days vs no start within 30 days
# ------------------------------------------------------------
# Required package: survival (base-recommended)
# Optional package used in formulas: splines (comes with R)

expit <- function(x) 1 / (1 + exp(-x))

# Natural treatment initiation model.
# Before day 30, treatment is symptom-driven and more confounded.
# On day 30, initiation can optionally become more protocolized because patients
# attend a scheduled cardiology review visit. For backward compatibility,
# the older \`day30_review_bump\` interface is also supported.

.resolve_common_length <- function(...) {
  lens <- vapply(list(...), length, integer(1))
  n <- max(lens)
  if (!all(lens %in% c(1L, n))) {
    stop("All inputs must have length 1 or a common length.", call. = FALSE)
  }
  n
}

.recycle_to_length <- function(x, n) {
  if (length(x) == n) return(x)
  rep(x, n)
}

.resolve_day30_mode <- function(day30_mode = c("auto", "protocolized", "legacy_bump", "none"),
                                protocolized_day30 = TRUE,
                                day30_review_bump = NULL) {
  day30_mode <- match.arg(day30_mode)
  if (day30_mode == "auto") {
    if (!is.null(day30_review_bump)) {
      return("legacy_bump")
    }
    if (isTRUE(protocolized_day30)) {
      return("protocolized")
    }
    return("none")
  }
  day30_mode
}

compute_lp_start_post_mi <- function(day,
                                     L_now,
                                     age10,
                                     bleed,
                                     frailty,
                                     large_mi,
                                     grace = 30,
                                     day30_mode = c("auto", "protocolized", "legacy_bump", "none"),
                                     protocolized_day30 = TRUE,
                                     day30_review_bump = NULL,
                                     day30_lp_intercept = 0.20,
                                     day30_lp_L = 0.25,
                                     day30_lp_bleed = -0.35,
                                     day30_lp_age10 = -0.05,
                                     day30_lp_frailty = -0.05,
                                     day30_lp_large_mi = 0.15) {
  n <- .resolve_common_length(day, L_now, age10, bleed, frailty, large_mi)
  day <- .recycle_to_length(day, n)
  L_now <- .recycle_to_length(L_now, n)
  age10 <- .recycle_to_length(age10, n)
  bleed <- .recycle_to_length(bleed, n)
  frailty <- .recycle_to_length(frailty, n)
  large_mi <- .recycle_to_length(large_mi, n)

  resolved_mode <- .resolve_day30_mode(
    day30_mode = day30_mode,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump
  )

  natural_lp <- -4.4 +
    0.75 * L_now -
    1.35 * bleed -
    0.15 * age10 -
    0.10 * frailty +
    0.45 * large_mi

  natural_lp[day > grace] <- natural_lp[day > grace] - 0.35
  out <- natural_lp

  is_day30 <- day == grace
  if (any(is_day30)) {
    if (resolved_mode == "protocolized") {
      out[is_day30] <- day30_lp_intercept +
        day30_lp_L * L_now[is_day30] +
        day30_lp_bleed * bleed[is_day30] +
        day30_lp_age10 * age10[is_day30] +
        day30_lp_frailty * frailty[is_day30] +
        day30_lp_large_mi * large_mi[is_day30]
    } else if (resolved_mode == "legacy_bump") {
      bump <- if (is.null(day30_review_bump)) 0 else day30_review_bump
      out[is_day30] <- natural_lp[is_day30] + bump
    }
  }

  out
}

run_tests_compute_lp_start_post_mi <- function(tol = 1e-10, verbose = TRUE) {
  assert_close <- function(observed, expected, label) {
    ok <- isTRUE(all.equal(observed, expected, tolerance = tol, check.attributes = FALSE))
    data.frame(test = label, passed = ok, stringsAsFactors = FALSE)
  }

  rows <- list()
  k <- 0L

  # 1) Pre-grace natural model is unchanged by day-30 settings.
  obs1 <- compute_lp_start_post_mi(
    day = 29,
    L_now = c(-0.5, 0.25),
    age10 = c(0, 1),
    bleed = c(0, 1),
    frailty = c(-0.3, 0.7),
    large_mi = c(0, 1),
    day30_mode = "protocolized",
    day30_lp_intercept = 99
  )
  exp1 <- -4.4 +
    0.75 * c(-0.5, 0.25) -
    1.35 * c(0, 1) -
    0.15 * c(0, 1) -
    0.10 * c(-0.3, 0.7) +
    0.45 * c(0, 1)
  k <- k + 1L; rows[[k]] <- assert_close(obs1, exp1, "pre_grace_uses_natural_model")

  # 2) Post-grace model subtracts 0.35 from the natural LP.
  obs2 <- compute_lp_start_post_mi(
    day = 31,
    L_now = 0.4,
    age10 = -0.5,
    bleed = 1,
    frailty = 0.2,
    large_mi = 1,
    day30_mode = "none"
  )
  exp2 <- (-4.4 + 0.75 * 0.4 - 1.35 * 1 - 0.15 * (-0.5) - 0.10 * 0.2 + 0.45 * 1) - 0.35
  k <- k + 1L; rows[[k]] <- assert_close(obs2, exp2, "post_grace_subtracts_035")

  # 3) Explicit protocolized day-30 model is used on day 30.
  obs3 <- compute_lp_start_post_mi(
    day = 30,
    L_now = c(0, 1),
    age10 = c(0.2, -0.4),
    bleed = c(1, 0),
    frailty = c(0.5, -0.5),
    large_mi = c(0, 1),
    day30_mode = "protocolized",
    day30_lp_intercept = 0.2,
    day30_lp_L = 0.25,
    day30_lp_bleed = -0.35,
    day30_lp_age10 = -0.05,
    day30_lp_frailty = -0.05,
    day30_lp_large_mi = 0.15
  )
  exp3 <- 0.2 +
    0.25 * c(0, 1) -
    0.35 * c(1, 0) -
    0.05 * c(0.2, -0.4) -
    0.05 * c(0.5, -0.5) +
    0.15 * c(0, 1)
  k <- k + 1L; rows[[k]] <- assert_close(obs3, exp3, "day30_protocolized_formula")

  # 4) Backward-compatible legacy bump interface works on day 30.
  obs4 <- compute_lp_start_post_mi(
    day = 30,
    L_now = 0.6,
    age10 = 0.3,
    bleed = 1,
    frailty = -0.1,
    large_mi = 0,
    day30_review_bump = 2.5
  )
  nat4 <- -4.4 + 0.75 * 0.6 - 1.35 * 1 - 0.15 * 0.3 - 0.10 * (-0.1) + 0.45 * 0
  exp4 <- nat4 + 2.5
  k <- k + 1L; rows[[k]] <- assert_close(obs4, exp4, "day30_legacy_bump_formula")

  # 5) Auto mode defaults to protocolized when no legacy bump is supplied.
  obs5 <- compute_lp_start_post_mi(
    day = 30,
    L_now = 0.2,
    age10 = 0.1,
    bleed = 0,
    frailty = 0.4,
    large_mi = 1,
    day30_mode = "auto",
    protocolized_day30 = TRUE,
    day30_lp_intercept = 0.2,
    day30_lp_L = 0.25,
    day30_lp_bleed = -0.35,
    day30_lp_age10 = -0.05,
    day30_lp_frailty = -0.05,
    day30_lp_large_mi = 0.15
  )
  exp5 <- 0.2 + 0.25 * 0.2 - 0.35 * 0 - 0.05 * 0.1 - 0.05 * 0.4 + 0.15 * 1
  k <- k + 1L; rows[[k]] <- assert_close(obs5, exp5, "auto_mode_protocolized_default")

  # 6) Higher L increases LP in both natural and protocolized modes.
  low_L <- compute_lp_start_post_mi(30, L_now = -1, age10 = 0, bleed = 0, frailty = 0, large_mi = 0,
                                    day30_mode = "protocolized")
  high_L <- compute_lp_start_post_mi(30, L_now = 1, age10 = 0, bleed = 0, frailty = 0, large_mi = 0,
                                     day30_mode = "protocolized")
  k <- k + 1L; rows[[k]] <- data.frame(test = "higher_L_increases_lp_on_day30", passed = (high_L > low_L), stringsAsFactors = FALSE)

  res <- do.call(rbind, rows)
  if (isTRUE(verbose)) {
    print(res)
    if (all(res$passed)) {
      message("All compute_lp_start_post_mi() tests passed.")
    } else {
      message("At least one compute_lp_start_post_mi() test failed.")
    }
  }
  res
}

show_compute_lp_start_post_mi_examples <- function() {
  ex <- expand.grid(
    day = c(29, 30, 31),
    L_now = c(-1, 0, 1),
    bleed = c(0, 1),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  ex$age10 <- 0
  ex$frailty <- 0
  ex$large_mi <- 0
  ex$lp_protocolized <- compute_lp_start_post_mi(
    day = ex$day,
    L_now = ex$L_now,
    age10 = ex$age10,
    bleed = ex$bleed,
    frailty = ex$frailty,
    large_mi = ex$large_mi,
    day30_mode = "protocolized"
  )
  ex$p_protocolized <- expit(ex$lp_protocolized)
  ex$lp_legacy_bump <- compute_lp_start_post_mi(
    day = ex$day,
    L_now = ex$L_now,
    age10 = ex$age10,
    bleed = ex$bleed,
    frailty = ex$frailty,
    large_mi = ex$large_mi,
    day30_review_bump = 2.5
  )
  ex$p_legacy_bump <- expit(ex$lp_legacy_bump)
  ex
}

# -----------------------------------------------------------------------------
# 1) Data-generating model
# -----------------------------------------------------------------------------
# Interpretation:
# - Day 0 = hospital discharge after myocardial infarction (MI; heart attack)
# - start_day = day the antiplatelet is started at the END of that day
#               (so the medication can affect outcomes starting the next day)
# - event_day = day of a serious cardiovascular event during follow-up
# - outcome = death or urgent cardiovascular readmission (composite example)
# -----------------------------------------------------------------------------

simulate_post_mi <- function(n = 2000,
                             max_follow = 180,
                             grace = 30,
                             seed = 20260406,
                             protocolized_day30 = TRUE,
                             day30_review_bump = NULL,
                             day30_lp_intercept = 0.20,
                             day30_lp_L = 0.25,
                             day30_lp_bleed = -0.35,
                             day30_lp_age10 = -0.05,
                             day30_lp_frailty = -0.05,
                             day30_lp_large_mi = 0.15) {
  set.seed(seed)

  id <- seq_len(n)

  # Baseline covariates
  age <- pmin(pmax(rnorm(n, mean = 68, sd = 10), 40), 95)
  age10 <- (age - 68) / 10
  bleed <- rbinom(n, 1, 0.20)      # high bleeding risk / contraindication
  frailty <- rnorm(n, 0, 1)        # general frailty / multimorbidity burden
  large_mi <- rbinom(n, 1, 0.35)   # index MI was relatively large/severe

  baseline_list <- vector("list", n)
  daily_list <- vector("list", n)

  for (i in seq_len(n)) {
    # Initial daily instability (higher = more symptoms / physician concern / ischemic instability)
    L_prev <- 0.6 * age10[i] +
      0.9 * bleed[i] +
      0.6 * frailty[i] +
      0.8 * large_mi[i] +
      rnorm(1, 0, 1)

    started <- FALSE
    start_day <- Inf
    event <- 0L
    event_day <- max_follow

    person_days <- vector("list", max_follow)
    keep_days <- 0L

    for (day in seq_len(max_follow)) {
      # Drug is active only if started at the END of a previous day
      active_startday <- as.integer(started && (start_day < day))

      if (day > 1) {
        # Daily instability evolves over time and improves somewhat once treatment is active
        L_prev <- 0.75 * L_prev +
          0.10 * age10[i] +
          0.15 * frailty[i] +
          0.25 * large_mi[i] -
          0.35 * active_startday +
          rnorm(1, mean = 0, sd = 0.6)
      }
      L_now <- L_prev

      # Daily event hazard during the day
      lp_event <- -5.9 +
        0.65 * L_now +
        0.35 * age10[i] +
        0.65 * bleed[i] +
        0.45 * frailty[i] +
        0.45 * large_mi[i] -
        0.75 * active_startday

      p_event <- expit(lp_event)
      y <- rbinom(1, 1, p_event)

      keep_days <- keep_days + 1L
      person_days[[keep_days]] <- data.frame(
        id = id[i],
        day = day,
        L = L_now,
        active_startday = active_startday,
        stringsAsFactors = FALSE
      )

      # Event happens DURING the day
      if (y == 1) {
        event <- 1L
        event_day <- day
        break
      }

      # Treatment can start at the END of the day
      if (!started) {
        lp_start <- compute_lp_start_post_mi(
          day = day,
          L_now = L_now,
          age10 = age10[i],
          bleed = bleed[i],
          frailty = frailty[i],
          large_mi = large_mi[i],
          grace = grace,
          day30_mode = "auto",
          protocolized_day30 = protocolized_day30,
          day30_review_bump = day30_review_bump,
          day30_lp_intercept = day30_lp_intercept,
          day30_lp_L = day30_lp_L,
          day30_lp_bleed = day30_lp_bleed,
          day30_lp_age10 = day30_lp_age10,
          day30_lp_frailty = day30_lp_frailty,
          day30_lp_large_mi = day30_lp_large_mi
        )

        p_start <- expit(lp_start)
        a <- rbinom(1, 1, p_start)
        if (a == 1) {
          started <- TRUE
          start_day <- day
        }
      }
    }

    if (!is.finite(start_day)) {
      start_day <- NA_integer_
    }

    baseline_list[[i]] <- data.frame(
      id = id[i],
      age = age[i],
      age10 = age10[i],
      bleed = bleed[i],
      frailty = frailty[i],
      large_mi = large_mi[i],
      start_day = start_day,
      event = event,
      event_day = event_day,
      stringsAsFactors = FALSE
    )

    daily_list[[i]] <- do.call(rbind, person_days[seq_len(keep_days)])
  }

  baseline <- do.call(rbind, baseline_list)
  daily <- do.call(rbind, daily_list)

  list(
    baseline = baseline,
    daily = daily,
    meta = list(max_follow = max_follow, grace = grace, seed = seed)
  )
}

# -----------------------------------------------------------------------------
# 2) Clone each patient into two strategies
#    arm 1: start_by_30     = initiate within 30 days
#    arm 2: no_start_by_30  = do not initiate within 30 days
# -----------------------------------------------------------------------------
# Timing convention:
# - If actual start_day = t, treatment starts at END of day t
# - Therefore a no_start_by_30 clone is artificially censored at END of day t
# - If a start_by_30 clone has still not started by END of day 30, it is censored at END of day 30
# - If an event occurs on day t and treatment also starts on day t, the event is counted first
#   because the event occurs during the day and the drug starts only at the day's end.
# -----------------------------------------------------------------------------

make_ccw_clones <- function(baseline,
                            daily,
                            max_follow = 180,
                            grace = 30) {

  daily_full <- merge(
    daily,
    baseline[, c("id", "age", "age10", "bleed", "frailty", "large_mi",
                 "start_day", "event_day", "event")],
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )

  out <- vector("list", nrow(baseline) * 2)
  idx <- 0L

  arms <- c("start_by_30", "no_start_by_30")

  for (arm in arms) {
    for (i in seq_len(nrow(baseline))) {
      s <- baseline$start_day[i]
      # event_day is meaningful only if an event actually occurred.
      # Otherwise administrative follow-up at day 180 would be misread as an event.
      e <- if (baseline$event[i] == 1L) baseline$event_day[i] else Inf
      pid <- baseline$id[i]

      # Arm-specific artificial censoring day
      if (arm == "start_by_30") {
        cday <- if (is.na(s) || s > grace) grace else Inf
      } else {
        cday <- if (!is.na(s) && s <= grace) s else Inf
      }

      # Event has priority if it happens on or before the censoring day
      if (e <= cday && e <= max_follow) {
        final_day <- e
        final_event <- 1L
        final_censor <- 0L
      } else {
        final_day <- min(cday, max_follow)
        final_event <- 0L
        final_censor <- as.integer(is.finite(cday) && cday <= max_follow)
      }

      sub <- daily_full[daily_full$id == pid & daily_full$day <= final_day, ]
      if (nrow(sub) == 0L) next

      sub <- sub[order(sub$day), ]
      sub$arm <- arm
      sub$clone <- paste0(pid, "_", arm)
      sub$tstart <- sub$day - 1
      sub$tstop <- sub$day
      sub$event_ccw <- 0L
      sub$censor_ccw <- 0L
      sub$event_ccw[nrow(sub)] <- final_event
      sub$censor_ccw[nrow(sub)] <- final_censor

      if (arm == "start_by_30") {
        # At risk of artificial censoring only up to day 30 and only until the observed start day
        at_risk <- (sub$day <= grace) & (is.na(s) | (sub$day <= s))
      } else {
        # No-start arm can be censored on any day 1..30 if treatment starts that day
        at_risk <- (sub$day <= grace)
      }
      sub$at_risk_censor <- as.integer(at_risk)

      idx <- idx + 1L
      out[[idx]] <- sub
    }
  }

  clone_long <- do.call(rbind, out[seq_len(idx)])
  rownames(clone_long) <- NULL

  clone_long$arm <- factor(clone_long$arm,
                           levels = c("no_start_by_30", "start_by_30"))
  clone_long
}

# -----------------------------------------------------------------------------
# 3) Fit censoring models and create weights
# -----------------------------------------------------------------------------
# Standard CCW uses IPCW.
# Here we compute two choices:
#   (A) w_ipcw = stabilized inverse probability of censoring weights (recommended default)
#   (B) w_ow   = overlap-style weights based on cumulative uncensoring probability
#                (sensitivity analysis, not the standard CCW default estimand)
#
# The asymmetry is important:
# - no_start_by_30 arm can deviate on any day 1..30 if treatment starts
# - start_by_30 arm can deviate only at the end of day 30, and only among clones still untreated then
# -----------------------------------------------------------------------------

add_ccw_weights <- function(clone_long,
                            grace = 30,
                            eps = 1e-6,
                            truncate = c(0.01, 0.99)) {

  dat <- clone_long[order(clone_long$id, clone_long$arm, clone_long$day), ]
  dat$w_ipcw <- 1
  dat$w_ow <- 1

  # ---------- No-start-by-30 arm: day-by-day censoring hazard ----------
  risk_no30 <- subset(dat,
                      arm == "no_start_by_30" &
                        at_risk_censor == 1 &
                        event_ccw == 0)

  if (nrow(risk_no30) > 0L) {
    fit_no30_den <- glm(
      censor_ccw ~ day + I(day^2) + age10 + bleed + frailty + large_mi + L,
      family = binomial(),
      data = risk_no30
    )

    fit_no30_num <- glm(
      censor_ccw ~ day + I(day^2),
      family = binomial(),
      data = risk_no30
    )

    risk_no30$p_den <- pmin(pmax(predict(fit_no30_den, type = "response"), eps), 1 - eps)
    risk_no30$p_num <- pmin(pmax(predict(fit_no30_num, type = "response"), eps), 1 - eps)

    # Probability of staying uncensored through the END of this day
    risk_no30$u_den <- 1 - risk_no30$p_den
    risk_no30$u_num <- 1 - risk_no30$p_num

    # Merge predictions back
    key <- paste(risk_no30$clone, risk_no30$day, sep = "__")
    dat$key_tmp <- paste(dat$clone, dat$day, sep = "__")
    dat$u_den_no30 <- 1
    dat$u_num_no30 <- 1
    m <- match(dat$key_tmp, key)
    has <- !is.na(m)
    dat$u_den_no30[has] <- risk_no30$u_den[m[has]]
    dat$u_num_no30[has] <- risk_no30$u_num[m[has]]
    dat$key_tmp <- NULL
  } else {
    dat$u_den_no30 <- 1
    dat$u_num_no30 <- 1
  }

  # ---------- Start-by-30 arm: censoring only at the END of day 30 ----------
  risk_start30 <- subset(dat,
                         arm == "start_by_30" &
                           day == grace &
                           at_risk_censor == 1 &
                           event_ccw == 0)

  if (nrow(risk_start30) > 0L) {
    fit_start30_den <- glm(
      censor_ccw ~ age10 + bleed + frailty + large_mi + L,
      family = binomial(),
      data = risk_start30
    )

    fit_start30_num <- glm(
      censor_ccw ~ 1,
      family = binomial(),
      data = risk_start30
    )

    risk_start30$p_den <- pmin(pmax(predict(fit_start30_den, type = "response"), eps), 1 - eps)
    risk_start30$p_num <- pmin(pmax(predict(fit_start30_num, type = "response"), eps), 1 - eps)

    risk_start30$u_den <- 1 - risk_start30$p_den
    risk_start30$u_num <- 1 - risk_start30$p_num

    key <- paste(risk_start30$clone, risk_start30$day, sep = "__")
    dat$key_tmp <- paste(dat$clone, dat$day, sep = "__")
    dat$u_den_start30 <- 1
    dat$u_num_start30 <- 1
    m <- match(dat$key_tmp, key)
    has <- !is.na(m)
    dat$u_den_start30[has] <- risk_start30$u_den[m[has]]
    dat$u_num_start30[has] <- risk_start30$u_num[m[has]]
    dat$key_tmp <- NULL
  } else {
    dat$u_den_start30 <- 1
    dat$u_num_start30 <- 1
  }

  # ---------- Build cumulative weights clone by clone ----------
  split_dat <- split(dat, dat$clone)
  weighted_list <- vector("list", length(split_dat))
  j <- 0L

  for (nm in names(split_dat)) {
    d <- split_dat[[nm]]
    d <- d[order(d$day), ]

    G_den <- rep(1, nrow(d))  # cumulative denominator prob of staying uncensored through START of interval
    G_num <- rep(1, nrow(d))  # cumulative numerator prob through START of interval

    # "End of interval" version is used only for the overlap-style sensitivity analysis
    G_den_end <- rep(1, nrow(d))
    G_num_end <- rep(1, nrow(d))

    for (k in seq_len(nrow(d))) {
      if (k == 1L) {
        G_den[k] <- 1
        G_num[k] <- 1
      } else {
        G_den[k] <- G_den_end[k - 1L]
        G_num[k] <- G_num_end[k - 1L]
      }

      if (as.character(d$arm[1]) == "no_start_by_30") {
        u_den_k <- d$u_den_no30[k]
        u_num_k <- d$u_num_no30[k]
      } else {
        u_den_k <- d$u_den_start30[k]
        u_num_k <- d$u_num_start30[k]
      }

      G_den_end[k] <- G_den[k] * u_den_k
      G_num_end[k] <- G_num[k] * u_num_k
    }

    # Standard stabilized IPCW
    d$w_ipcw <- G_num / pmax(G_den, eps)

    # Overlap-style sensitivity analysis on the cumulative uncensoring score.
    # Uses 1 - G_den_end(t): emphasizes clones/intervals with intermediate adherence probability.
    # This is NOT the standard CCW estimand, but a pragmatic sensitivity analysis when inverse
    # weights become unstable.
    d$w_ow <- pmax(1 - G_den_end, eps)

    j <- j + 1L
    weighted_list[[j]] <- d
  }

  dat_w <- do.call(rbind, weighted_list)
  rownames(dat_w) <- NULL

  # Optional truncation of IPCW
  if (!is.null(truncate)) {
    q <- quantile(dat_w$w_ipcw, probs = truncate, na.rm = TRUE)
    dat_w$w_ipcw_trunc <- pmin(pmax(dat_w$w_ipcw, q[1]), q[2])
  } else {
    dat_w$w_ipcw_trunc <- dat_w$w_ipcw
  }

  # Rescale overlap-style weights to mean 1 (for numerical convenience)
  dat_w$w_ow <- dat_w$w_ow / mean(dat_w$w_ow, na.rm = TRUE)

  dat_w
}

# -----------------------------------------------------------------------------
# 4) Weighted Cox regression
# -----------------------------------------------------------------------------
# We cluster on the original patient id because the two clones from the same person
# are dependent.
# -----------------------------------------------------------------------------

fit_ccw_cox <- function(weighted_clone_long,
                        weight_var = c("w_ipcw_trunc", "w_ipcw", "w_ow")) {
  weight_var <- match.arg(weight_var)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  survival::coxph(
    survival::Surv(tstart, tstop, event_ccw) ~ arm + cluster(id),
    data = weighted_clone_long,
    weights = weighted_clone_long[[weight_var]],
    robust = TRUE,
    ties = "efron"
  )
}

# -----------------------------------------------------------------------------
# 5) Naive immortal-time-biased analyses (for contrast only)
# -----------------------------------------------------------------------------

fit_naive_cox <- function(baseline) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  baseline$group_naive <- ifelse(!is.na(baseline$start_day) & baseline$start_day <= 30,
                                 "start_by_30", "no_start_by_30")
  baseline$group_naive <- factor(baseline$group_naive,
                                 levels = c("no_start_by_30", "start_by_30"))

  survival::coxph(
    survival::Surv(event_day, event) ~ group_naive,
    data = baseline,
    ties = "efron"
  )
}

fit_naive_logistic180 <- function(baseline,
                                  grace = 30) {
  dat <- baseline
  dat$arm_naive <- factor(
    ifelse(!is.na(dat$start_day) & dat$start_day <= grace,
           "start_by_30", "no_start_by_30"),
    levels = c("no_start_by_30", "start_by_30")
  )

  fit <- glm(
    event ~ arm_naive,
    family = binomial(),
    data = dat
  )

  structure(
    list(
      fit = fit,
      arm_var = "arm_naive",
      arm_levels = levels(dat$arm_naive)
    ),
    class = "post_mi_logistic180_fit"
  )
}

predict_logistic180_risk <- function(fit_object,
                                     arm_levels = NULL) {
  if (inherits(fit_object, "glm")) {
    stop("Please pass the object returned by fit_naive_logistic180().")
  }

  fit <- fit_object$fit
  arm_var <- fit_object$arm_var

  if (is.null(arm_levels)) {
    arm_levels <- fit_object$arm_levels
  }

  newdat <- data.frame(dummy = seq_along(arm_levels))
  newdat[[arm_var]] <- factor(arm_levels, levels = fit_object$arm_levels)
  newdat$dummy <- NULL

  p <- as.numeric(predict(fit, newdata = newdat, type = "response"))
  out <- data.frame(
    arm = arm_levels,
    risk = p,
    stringsAsFactors = FALSE
  )

  if (all(c("no_start_by_30", "start_by_30") %in% out$arm)) {
    p0 <- out$risk[out$arm == "no_start_by_30"]
    p1 <- out$risk[out$arm == "start_by_30"]
    contrast <- data.frame(
      contrast = c("risk_difference", "risk_ratio"),
      estimate = c(p1 - p0, p1 / pmax(p0, 1e-8)),
      stringsAsFactors = FALSE
    )
  } else {
    contrast <- NULL
  }

  list(risk_180 = out, contrasts = contrast)
}

# -----------------------------------------------------------------------------
# 6) Monte Carlo truth under explicit intervention regimes
# -----------------------------------------------------------------------------
# We directly simulate the two strategies of interest:
#   g_start_by_30     = allow natural start through day 30; if still untreated at
#                       the END of day 30, force treatment start at the END of day 30
#   g_no_start_by_30  = prohibit treatment start through the END of day 30;
#                       from day 31 onward, allow the natural treatment process again
#
# The Monte Carlo truth is therefore on the 180-day RISK scale, not the hazard-ratio scale.
# -----------------------------------------------------------------------------

simulate_regime_outcomes_chunk <- function(age10,
                                           bleed,
                                           frailty,
                                           large_mi,
                                           z0,
                                           zL,
                                           u_event,
                                           u_start,
                                           regime = c("start_by_30", "no_start_by_30"),
                                           max_follow = 180,
                                           grace = 30,
                                           protocolized_day30 = TRUE,
                                           day30_review_bump = NULL,
                                           day30_lp_intercept = 0.20,
                                           day30_lp_L = 0.25,
                                           day30_lp_bleed = -0.35,
                                           day30_lp_age10 = -0.05,
                                           day30_lp_frailty = -0.05,
                                           day30_lp_large_mi = 0.15) {
  regime <- match.arg(regime)

  n <- length(age10)

  L_prev <- 0.6 * age10 +
    0.9 * bleed +
    0.6 * frailty +
    0.8 * large_mi +
    z0

  started <- rep(FALSE, n)
  start_day <- rep(NA_integer_, n)
  event <- integer(n)
  event_day <- rep(max_follow, n)
  alive <- rep(TRUE, n)

  for (day in seq_len(max_follow)) {
    active_startday <- as.integer(started & !is.na(start_day) & (start_day < day))

    if (day > 1L) {
      L_prev <- 0.75 * L_prev +
        0.10 * age10 +
        0.15 * frailty +
        0.25 * large_mi -
        0.35 * active_startday +
        zL[, day - 1L]
    }
    L_now <- L_prev

    # Event during the day
    idx_alive <- which(alive)
    if (length(idx_alive) == 0L) break

    lp_event <- -5.9 +
      0.65 * L_now[idx_alive] +
      0.35 * age10[idx_alive] +
      0.65 * bleed[idx_alive] +
      0.45 * frailty[idx_alive] +
      0.45 * large_mi[idx_alive] -
      0.75 * active_startday[idx_alive]

    p_event <- expit(lp_event)
    y_today <- u_event[idx_alive, day] < p_event

    if (any(y_today)) {
      ids <- idx_alive[y_today]
      event[ids] <- 1L
      event_day[ids] <- day
      alive[ids] <- FALSE
    }

    # Treatment can start at the END of the day, only among those still alive and untreated
    idx_treat <- which(alive & !started)
    if (length(idx_treat) == 0L) next

    lp_start_nat <- compute_lp_start_post_mi(
      day = day,
      L_now = L_now[idx_treat],
      age10 = age10[idx_treat],
      bleed = bleed[idx_treat],
      frailty = frailty[idx_treat],
      large_mi = large_mi[idx_treat],
      grace = grace,
      day30_mode = "auto",
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi
    )

    p_start_nat <- expit(lp_start_nat)
    natural_start <- u_start[idx_treat, day] < p_start_nat
    start_today <- rep(FALSE, length(idx_treat))

    if (regime == "start_by_30") {
      if (day < grace) {
        start_today <- natural_start
      } else if (day == grace) {
        # Force everyone who survives untreated to be started by the END of day 30
        start_today <- rep(TRUE, length(idx_treat))
      } else {
        # After day 30 everyone alive should already be started under this regime
        start_today <- rep(FALSE, length(idx_treat))
      }
    } else {
      # no_start_by_30
      if (day <= grace) {
        start_today <- rep(FALSE, length(idx_treat))
      } else {
        start_today <- natural_start
      }
    }

    if (any(start_today)) {
      ids <- idx_treat[start_today]
      started[ids] <- TRUE
      start_day[ids] <- day
    }
  }

  alive_through_grace <- event_day > grace
  started_by_grace <- !is.na(start_day) & (start_day <= grace)

  list(
    event = event,
    event_day = event_day,
    start_day = start_day,
    alive_through_grace = alive_through_grace,
    started_by_grace = started_by_grace
  )
}

monte_carlo_truth_post_mi <- function(n_mc = 200000,
                                      max_follow = 180,
                                      grace = 30,
                                      chunk_size = 20000,
                                      seed = 20260406,
                                      conf.level = 0.95,
                                      protocolized_day30 = TRUE,
                                      day30_review_bump = NULL,
                                      day30_lp_intercept = 0.20,
                                      day30_lp_L = 0.25,
                                      day30_lp_bleed = -0.35,
                                      day30_lp_age10 = -0.05,
                                      day30_lp_frailty = -0.05,
                                      day30_lp_large_mi = 0.15) {
  set.seed(seed)

  z_alpha <- qnorm(1 - (1 - conf.level) / 2)

  y_start <- integer(n_mc)
  y_no30 <- integer(n_mc)

  alive_start_grace <- 0L
  alive_no30_grace <- 0L
  started_start_grace <- 0L
  started_no30_grace <- 0L

  offset <- 0L
  n_chunks <- ceiling(n_mc / chunk_size)

  for (b in seq_len(n_chunks)) {
    m <- min(chunk_size, n_mc - offset)
    idx <- seq.int(offset + 1L, offset + m)

    age <- pmin(pmax(rnorm(m, mean = 68, sd = 10), 40), 95)
    age10 <- (age - 68) / 10
    bleed <- rbinom(m, 1, 0.20)
    frailty <- rnorm(m, 0, 1)
    large_mi <- rbinom(m, 1, 0.35)

    z0 <- rnorm(m, mean = 0, sd = 1)
    zL <- matrix(rnorm(m * (max_follow - 1L), mean = 0, sd = 0.6),
                 nrow = m, ncol = max_follow - 1L)
    u_event <- matrix(runif(m * max_follow), nrow = m, ncol = max_follow)
    u_start <- matrix(runif(m * max_follow), nrow = m, ncol = max_follow)

    sim_start <- simulate_regime_outcomes_chunk(
      age10 = age10,
      bleed = bleed,
      frailty = frailty,
      large_mi = large_mi,
      z0 = z0,
      zL = zL,
      u_event = u_event,
      u_start = u_start,
      regime = "start_by_30",
      max_follow = max_follow,
      grace = grace,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi
    )

    sim_no30 <- simulate_regime_outcomes_chunk(
      age10 = age10,
      bleed = bleed,
      frailty = frailty,
      large_mi = large_mi,
      z0 = z0,
      zL = zL,
      u_event = u_event,
      u_start = u_start,
      regime = "no_start_by_30",
      max_follow = max_follow,
      grace = grace,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi
    )

    y_start[idx] <- sim_start$event
    y_no30[idx] <- sim_no30$event

    alive_start_grace <- alive_start_grace + sum(sim_start$alive_through_grace)
    alive_no30_grace <- alive_no30_grace + sum(sim_no30$alive_through_grace)
    started_start_grace <- started_start_grace +
      sum(sim_start$started_by_grace & sim_start$alive_through_grace)
    started_no30_grace <- started_no30_grace +
      sum(sim_no30$started_by_grace & sim_no30$alive_through_grace)

    offset <- offset + m
  }

  p_start <- mean(y_start)
  p_no30 <- mean(y_no30)

  se_start <- sqrt(p_start * (1 - p_start) / n_mc)
  se_no30 <- sqrt(p_no30 * (1 - p_no30) / n_mc)

  rd <- p_start - p_no30
  rr <- p_start / pmax(p_no30, 1e-8)

  se_rd <- sd(y_start - y_no30) / sqrt(n_mc)

  p_start_safe <- pmax(p_start, 1e-8)
  p_no30_safe <- pmax(p_no30, 1e-8)
  if_log_rr <- (y_start - p_start) / p_start_safe - (y_no30 - p_no30) / p_no30_safe
  se_log_rr <- sd(if_log_rr) / sqrt(n_mc)

  risk_180 <- data.frame(
    regime = c("start_by_30", "no_start_by_30"),
    risk = c(p_start, p_no30),
    mc_se = c(se_start, se_no30),
    lower = c(p_start - z_alpha * se_start, p_no30 - z_alpha * se_no30),
    upper = c(p_start + z_alpha * se_start, p_no30 + z_alpha * se_no30),
    stringsAsFactors = FALSE
  )

  contrasts <- data.frame(
    contrast = c("risk_difference", "risk_ratio"),
    estimate = c(rd, rr),
    mc_se = c(se_rd, se_log_rr),
    lower = c(rd - z_alpha * se_rd,
              exp(log(rr) - z_alpha * se_log_rr)),
    upper = c(rd + z_alpha * se_rd,
              exp(log(rr) + z_alpha * se_log_rr)),
    stringsAsFactors = FALSE
  )

  diagnostics <- data.frame(
    regime = c("start_by_30", "no_start_by_30"),
    alive_through_grace = c(alive_start_grace, alive_no30_grace),
    started_by_grace_among_alive = c(started_start_grace, started_no30_grace),
    proportion_started_by_grace_among_alive = c(
      started_start_grace / pmax(alive_start_grace, 1),
      started_no30_grace / pmax(alive_no30_grace, 1)
    ),
    stringsAsFactors = FALSE
  )

  list(
    risk_180 = risk_180,
    contrasts = contrasts,
    diagnostics = diagnostics,
    meta = list(
      n_mc = n_mc,
      max_follow = max_follow,
      grace = grace,
      chunk_size = chunk_size,
      seed = seed,
      conf.level = conf.level
    )
  )
}


run_day30_interface_smoke_tests <- function(verbose = TRUE) {
  test_rows <- list()
  k <- 0L

  lp_tests <- run_tests_compute_lp_start_post_mi(verbose = FALSE)
  k <- k + 1L
  test_rows[[k]] <- data.frame(
    test = "compute_lp_start_post_mi_unit_tests",
    passed = all(lp_tests$passed),
    stringsAsFactors = FALSE
  )

  sim <- simulate_post_mi(
    n = 50,
    max_follow = 60,
    grace = 30,
    seed = 123,
    day30_review_bump = 2.5
  )
  k <- k + 1L
  test_rows[[k]] <- data.frame(
    test = "simulate_post_mi_accepts_legacy_day30_review_bump",
    passed = is.list(sim) && all(c("baseline", "daily", "meta") %in% names(sim)) && nrow(sim$baseline) == 50,
    stringsAsFactors = FALSE
  )

  truth <- monte_carlo_truth_post_mi(
    n_mc = 2000,
    max_follow = 60,
    grace = 30,
    chunk_size = 500,
    seed = 321,
    day30_review_bump = 2.5
  )
  k <- k + 1L
  test_rows[[k]] <- data.frame(
    test = "monte_carlo_truth_accepts_legacy_day30_review_bump",
    passed = is.list(truth) && all(c("risk_180", "contrasts", "diagnostics", "meta") %in% names(truth)),
    stringsAsFactors = FALSE
  )

  diag_ok <- with(truth$diagnostics,
                  all(proportion_started_by_grace_among_alive[regime == "start_by_30"] == 1) &&
                    all(proportion_started_by_grace_among_alive[regime == "no_start_by_30"] == 0))
  k <- k + 1L
  test_rows[[k]] <- data.frame(
    test = "monte_carlo_regime_diagnostics_are_exact",
    passed = isTRUE(diag_ok),
    stringsAsFactors = FALSE
  )

  out <- do.call(rbind, test_rows)
  if (isTRUE(verbose)) {
    print(out)
    if (all(out$passed)) {
      message("All day-30 interface smoke tests passed.")
    } else {
      message("At least one day-30 interface smoke test failed.")
    }
  }
  out
}

# -----------------------------------------------------------------------------
# 7) Pooled logistic risk estimation
# -----------------------------------------------------------------------------
# This is the preferred way to estimate 180-day risk so that the estimand matches
# the Monte Carlo truth on the RISK scale.
#
# There are three useful versions:
#   (A) naive observed-group pooled logistic on the original cohort
#       -> immortal time bias
#   (B) clone-censor pooled logistic WITHOUT weights
#       -> informative censoring / selection bias
#   (C) clone-censor pooled logistic WITH censoring weights
#       -> target CCW analysis
#
# For validation against Monte Carlo truth, a saturated time model (time_model = "factor")
# is recommended. It avoids smoothing bias from overly restrictive splines and makes the
# pooled-logistic risk match the day-specific weighted hazard product-limit estimator.
# -----------------------------------------------------------------------------

make_naive_pooled_long <- function(baseline,
                                   daily,
                                   max_follow = 180,
                                   grace = 30) {
  dat <- merge(
    daily,
    baseline[, c("id", "age", "age10", "bleed", "frailty", "large_mi",
                 "start_day", "event_day", "event")],
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )

  dat <- dat[dat$day <= max_follow, ]
  dat <- dat[order(dat$id, dat$day), ]
  dat$arm <- factor(
    ifelse(!is.na(dat$start_day) & dat$start_day <= grace,
           "start_by_30", "no_start_by_30"),
    levels = c("no_start_by_30", "start_by_30")
  )
  dat$event_plr <- as.integer(dat$event == 1L & dat$day == dat$event_day)
  dat$tstart <- dat$day - 1L
  dat$tstop <- dat$day
  rownames(dat) <- NULL
  dat
}

aggregate_pooled_person_period <- function(long_data,
                                         outcome_var,
                                         arm_var = "arm",
                                         weight_var = NULL) {
  dat <- long_data
  dat$.arm <- factor(dat[[arm_var]],
                     levels = c("no_start_by_30", "start_by_30"))

  if (is.null(weight_var)) {
    dat$.w <- 1
  } else {
    dat$.w <- dat[[weight_var]]
  }

  dat$.event_w <- dat[[outcome_var]] * dat$.w

  agg <- stats::aggregate(
    cbind(.w, .event_w) ~ day + .arm,
    data = dat,
    FUN = sum
  )
  names(agg)[names(agg) == ".arm"] <- arm_var
  names(agg)[names(agg) == ".w"] <- "n_risk"
  names(agg)[names(agg) == ".event_w"] <- "events"
  agg$prop_event <- agg$events / pmax(agg$n_risk, 1e-8)
  agg
}

fit_pooled_logistic_model <- function(long_data,
                                      outcome_var,
                                      arm_var = "arm",
                                      weight_var = NULL,
                                      time_df = 5,
                                      time_model = c("factor", "ns"),
                                      eps = 1e-8) {
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("Package 'splines' is required.")
  }

  time_model <- match.arg(time_model)
  agg <- aggregate_pooled_person_period(
    long_data = long_data,
    outcome_var = outcome_var,
    arm_var = arm_var,
    weight_var = weight_var
  )
  agg[[arm_var]] <- factor(agg[[arm_var]],
                           levels = c("no_start_by_30", "start_by_30"))

  if (time_model == "factor") {
    form_txt <- paste0(
      "prop_event ~ ", arm_var,
      " * factor(day)"
    )
  } else {
    form_txt <- paste0(
      "prop_event ~ ", arm_var,
      " * splines::ns(day, df = ", time_df, ")"
    )
  }
  fml <- stats::as.formula(form_txt)

  fit <- glm(
    fml,
    family = quasibinomial(),
    data = agg,
    weights = agg$n_risk
  )

  structure(
    list(
      fit = fit,
      data = agg,
      outcome_var = outcome_var,
      arm_var = arm_var,
      weight_var = weight_var,
      time_df = time_df,
      time_model = time_model,
      eps = eps,
      arm_levels = levels(agg[[arm_var]])
    ),
    class = "post_mi_pooled_fit"
  )
}

predict_pooled_logistic_risk <- function(fit_object,
                                         horizon = 180,
                                         arm_levels = NULL) {
  if (inherits(fit_object, "glm")) {
    stop("Please pass the object returned by fit_pooled_logistic_model().")
  }

  fit <- fit_object$fit
  arm_var <- fit_object$arm_var
  eps <- fit_object$eps

  if (is.null(arm_levels)) {
    arm_levels <- fit_object$arm_levels
  }

  out <- vector("list", length(arm_levels))

  for (j in seq_along(arm_levels)) {
    newdat <- data.frame(day = seq_len(horizon))
    newdat[[arm_var]] <- factor(rep(arm_levels[j], horizon),
                                levels = fit_object$arm_levels)

    p_day <- as.numeric(predict(fit, newdata = newdat, type = "response"))
    p_day <- pmin(pmax(p_day, eps), 1 - eps)

    surv <- cumprod(1 - p_day)
    risk <- 1 - surv

    out[[j]] <- data.frame(
      arm = arm_levels[j],
      day = seq_len(horizon),
      hazard = p_day,
      survival = surv,
      risk = risk,
      stringsAsFactors = FALSE
    )
  }

  daily_risk <- do.call(rbind, out)
  rownames(daily_risk) <- NULL

  risk_180 <- daily_risk[daily_risk$day == horizon, c("arm", "risk")]
  rownames(risk_180) <- NULL

  if (all(c("no_start_by_30", "start_by_30") %in% risk_180$arm)) {
    p0 <- risk_180$risk[risk_180$arm == "no_start_by_30"]
    p1 <- risk_180$risk[risk_180$arm == "start_by_30"]

    contrasts <- data.frame(
      contrast = c("risk_difference", "risk_ratio"),
      estimate = c(p1 - p0, p1 / pmax(p0, eps)),
      stringsAsFactors = FALSE
    )
  } else {
    contrasts <- NULL
  }

  list(
    daily_risk = daily_risk,
    risk_180 = risk_180,
    contrasts = contrasts
  )
}

estimate_naive_pooled_risk <- function(baseline,
                                       daily,
                                       max_follow = 180,
                                       grace = 30,
                                       time_df = 5,
                                       time_model = c("factor", "ns")) {
  time_model <- match.arg(time_model)
  long_naive <- make_naive_pooled_long(
    baseline = baseline,
    daily = daily,
    max_follow = max_follow,
    grace = grace
  )

  fit_obj <- fit_pooled_logistic_model(
    long_data = long_naive,
    outcome_var = "event_plr",
    arm_var = "arm",
    weight_var = NULL,
    time_df = time_df,
    time_model = time_model
  )

  pred <- predict_pooled_logistic_risk(
    fit_object = fit_obj,
    horizon = max_follow
  )

  list(
    long_data = long_naive,
    fit = fit_obj,
    risk_180 = pred$risk_180,
    contrasts = pred$contrasts,
    daily_risk = pred$daily_risk
  )
}

estimate_ccw_pooled_risk <- function(clone_long,
                                     weight_var = NULL,
                                     max_follow = 180,
                                     time_df = 5,
                                     time_model = c("factor", "ns")) {
  time_model <- match.arg(time_model)
  fit_obj <- fit_pooled_logistic_model(
    long_data = clone_long,
    outcome_var = "event_ccw",
    arm_var = "arm",
    weight_var = weight_var,
    time_df = time_df,
    time_model = time_model
  )

  pred <- predict_pooled_logistic_risk(
    fit_object = fit_obj,
    horizon = max_follow
  )

  list(
    fit = fit_obj,
    risk_180 = pred$risk_180,
    contrasts = pred$contrasts,
    daily_risk = pred$daily_risk,
    weight_var = weight_var
  )
}

# Direct weighted day-specific hazard product-limit estimator.
# This is a useful sanity check for the pooled logistic fit.
# With time_model = "factor", the pooled-logistic risk should be essentially identical.
estimate_discrete_time_risk_from_hazards <- function(long_data,
                                                     outcome_var,
                                                     arm_var = "arm",
                                                     weight_var = NULL,
                                                     horizon = 180,
                                                     eps = 1e-8) {
  agg <- aggregate_pooled_person_period(
    long_data = long_data,
    outcome_var = outcome_var,
    arm_var = arm_var,
    weight_var = weight_var
  )
  agg[[arm_var]] <- factor(agg[[arm_var]],
                           levels = c("no_start_by_30", "start_by_30"))

  out <- vector("list", length(levels(agg[[arm_var]])))

  for (j in seq_along(levels(agg[[arm_var]]))) {
    arm_j <- levels(agg[[arm_var]])[j]
    sub <- agg[agg[[arm_var]] == arm_j, c("day", "n_risk", "events")]
    template <- data.frame(day = seq_len(horizon))
    sub <- merge(template, sub, by = "day", all.x = TRUE, sort = TRUE)
    sub$n_risk[is.na(sub$n_risk)] <- 0
    sub$events[is.na(sub$events)] <- 0
    sub$hazard <- ifelse(sub$n_risk > 0,
                         pmin(pmax(sub$events / sub$n_risk, eps), 1 - eps),
                         0)
    sub$survival <- cumprod(1 - sub$hazard)
    sub$risk <- 1 - sub$survival
    sub$arm <- arm_j
    out[[j]] <- sub[, c("arm", "day", "hazard", "survival", "risk")]
  }

  daily_risk <- do.call(rbind, out)
  rownames(daily_risk) <- NULL

  risk_180 <- daily_risk[daily_risk$day == horizon, c("arm", "risk")]
  rownames(risk_180) <- NULL

  p0 <- risk_180$risk[risk_180$arm == "no_start_by_30"]
  p1 <- risk_180$risk[risk_180$arm == "start_by_30"]
  contrasts <- data.frame(
    contrast = c("risk_difference", "risk_ratio"),
    estimate = c(p1 - p0, p1 / pmax(p0, eps)),
    stringsAsFactors = FALSE
  )

  list(
    daily_risk = daily_risk,
    risk_180 = risk_180,
    contrasts = contrasts
  )
}

compare_risk_estimates_to_truth <- function(estimate_object,
                                            truth_object) {
  est_risk <- estimate_object$risk_180
  names(est_risk)[names(est_risk) == "risk"] <- "estimate"

  true_risk <- truth_object$risk_180[, c("regime", "risk")]
  names(true_risk) <- c("arm", "truth")

  risk_comparison <- merge(est_risk, true_risk, by = "arm", all = TRUE)
  risk_comparison$bias <- risk_comparison$estimate - risk_comparison$truth

  if (!is.null(estimate_object$contrasts) && !is.null(truth_object$contrasts)) {
    est_con <- estimate_object$contrasts
    names(est_con)[names(est_con) == "estimate"] <- "estimate"

    true_con <- truth_object$contrasts[, c("contrast", "estimate")]
    names(true_con) <- c("contrast", "truth")

    contrast_comparison <- merge(est_con, true_con, by = "contrast", all = TRUE)
    contrast_comparison$bias <- contrast_comparison$estimate - contrast_comparison$truth
  } else {
    contrast_comparison <- NULL
  }

  list(
    risk = risk_comparison,
    contrasts = contrast_comparison
  )
}

# -----------------------------------------------------------------------------
# 8A) Support diagnostics for the start_by_30 arm
# -----------------------------------------------------------------------------
# Why this matters:
# Under the standard initiation-window CCW estimand used here, clones assigned to
# start_by_30 who are still untreated at the END of day 30 are represented by people
# who naturally start exactly on day 30. If those natural day-30 starters are too rare,
# the start arm has weak empirical support and finite-sample bias can remain even when
# the censoring model is correctly specified.
# -----------------------------------------------------------------------------

diagnose_start_by_30_support <- function(baseline,
                                         weighted_clone_long = NULL,
                                         grace = 30,
                                         weight_var = "w_ipcw") {
  riskset_nat <- subset(
    baseline,
    event_day > grace & (is.na(start_day) | start_day >= grace)
  )

  n_riskset <- nrow(riskset_nat)
  n_start_day30 <- sum(riskset_nat$start_day == grace, na.rm = TRUE)
  prop_start_day30 <- if (n_riskset > 0L) n_start_day30 / n_riskset else NA_real_

  out <- list(
    natural_day30_support = data.frame(
      n_event_free_untreated_through_day30 = n_riskset,
      n_natural_day30_starters = n_start_day30,
      proportion_natural_day30_starters = prop_start_day30,
      stringsAsFactors = FALSE
    )
  )

  if (!is.null(weighted_clone_long)) {
    risk_start30 <- subset(
      weighted_clone_long,
      arm == "start_by_30" &
        day == grace &
        at_risk_censor == 1 &
        event_ccw == 0
    )

    after_day30 <- subset(
      weighted_clone_long,
      arm == "start_by_30" &
        day == (grace + 1)
    )

    ess <- if (nrow(after_day30) > 0L) {
      w <- after_day30[[weight_var]]
      (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
    } else {
      NA_real_
    }

    out$weighted_day30_diagnostics <- data.frame(
      n_start_arm_day30_riskset = nrow(risk_start30),
      n_uncensored_natural_day30_starters = sum(risk_start30$censor_ccw == 0L, na.rm = TRUE),
      median_predicted_uncensor_prob = if (nrow(risk_start30) > 0L) median(risk_start30$u_den_start30, na.rm = TRUE) else NA_real_,
      p10_predicted_uncensor_prob = if (nrow(risk_start30) > 0L) quantile(risk_start30$u_den_start30, 0.10, na.rm = TRUE) else NA_real_,
      p90_predicted_uncensor_prob = if (nrow(risk_start30) > 0L) quantile(risk_start30$u_den_start30, 0.90, na.rm = TRUE) else NA_real_,
      effective_sample_size_day31 = ess,
      stringsAsFactors = FALSE
    )
  }

  out
}

# -----------------------------------------------------------------------------
# 8) Optional helper: translate a fitted Cox model into model-implied 180-day risks
# -----------------------------------------------------------------------------
# WARNING:
# - This does NOT make the Cox HR directly comparable to the Monte Carlo truth.
# - It only gives the risks IMPLIED BY the fitted proportional-hazards model.
# - If proportional hazards is violated, or if the Cox model is used only as a summary HR,
#   these model-implied risks may differ materially from the true intervention-specific risks.
# -----------------------------------------------------------------------------

predict_risk180_from_cox <- function(cox_fit,
                                     horizon = 180) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  bh <- survival::basehaz(cox_fit, centered = FALSE)
  bh <- bh[order(bh$time), ]

  if (all(bh$time < horizon)) {
    H0 <- max(bh$hazard)
  } else {
    H0 <- max(bh$hazard[bh$time <= horizon])
  }

  beta <- unname(stats::coef(cox_fit)[1])

  risk_no30 <- 1 - exp(-H0)
  risk_start <- 1 - exp(-H0 * exp(beta))

  risk_180 <- data.frame(
    arm = c("no_start_by_30", "start_by_30"),
    risk = c(risk_no30, risk_start),
    stringsAsFactors = FALSE
  )

  contrasts <- data.frame(
    contrast = c("risk_difference", "risk_ratio", "hazard_ratio"),
    estimate = c(risk_start - risk_no30,
                 risk_start / pmax(risk_no30, 1e-8),
                 exp(beta)),
    stringsAsFactors = FALSE
  )

  list(risk_180 = risk_180, contrasts = contrasts)
}

# -----------------------------------------------------------------------------
# 9) Example workflow
# -----------------------------------------------------------------------------
# source("ccw_post_mi_example.R")
#
# # Step 1: simulate one observed dataset
# sim <- simulate_post_mi(n = 5000, seed = 20260406)
#
# # Step 2: CCW cloning and censoring weights
# ccw <- make_ccw_clones(sim$baseline, sim$daily, grace = 30, max_follow = 180)
# ccw_w <- add_ccw_weights(ccw, grace = 30)
#
# # Step 3: Monte Carlo truth under explicit interventions
# truth_mc <- monte_carlo_truth_post_mi(
#   n_mc = 200000,
#   max_follow = 180,
#   grace = 30,
#   chunk_size = 20000,
#   seed = 1
# )
# truth_mc$risk_180
# truth_mc$contrasts
#
# # Step 4A: simple naive patient-level logistic regression (biased)
# fit_naive_logit180 <- fit_naive_logistic180(sim$baseline, grace = 30)
# naive_logit180_risk <- predict_logistic180_risk(fit_naive_logit180)
# naive_logit180_risk
#
# # Step 4B: naive pooled logistic on observed groups (immortal time bias)
# naive_plr <- estimate_naive_pooled_risk(
#   baseline = sim$baseline,
#   daily = sim$daily,
#   max_follow = 180,
#   grace = 30,
#   time_df = 5,
#   time_model = "factor"
# )
# naive_plr$risk_180
# naive_plr$contrasts
# compare_risk_estimates_to_truth(naive_plr, truth_mc)
#
# # Step 4C: clone-censor pooled logistic WITHOUT weights (selection bias)
# ccw_unweighted <- estimate_ccw_pooled_risk(
#   clone_long = ccw,
#   weight_var = NULL,
#   max_follow = 180,
#   time_df = 5,
#   time_model = "factor"
# )
# ccw_unweighted$risk_180
# ccw_unweighted$contrasts
# compare_risk_estimates_to_truth(ccw_unweighted, truth_mc)
#
# # Step 4D: clone-censor pooled logistic WITH stabilized IPCW (target analysis)
# ccw_weighted <- estimate_ccw_pooled_risk(
#   clone_long = ccw_w,
#   weight_var = "w_ipcw_trunc",
#   max_follow = 180,
#   time_df = 5,
#   time_model = "factor"
# )
# ccw_weighted$risk_180
# ccw_weighted$contrasts
# compare_risk_estimates_to_truth(ccw_weighted, truth_mc)
#
# # Sanity check: direct weighted day-specific hazard product-limit estimator
# ccw_weighted_hazard <- estimate_discrete_time_risk_from_hazards(
#   long_data = ccw_w,
#   outcome_var = "event_ccw",
#   arm_var = "arm",
#   weight_var = "w_ipcw_trunc",
#   horizon = 180
# )
# ccw_weighted_hazard$risk_180
# compare_risk_estimates_to_truth(ccw_weighted_hazard, truth_mc)
#
# # Optional overlap-style sensitivity analysis
# ccw_overlap <- estimate_ccw_pooled_risk(
#   clone_long = ccw_w,
#   weight_var = "w_ow",
#   max_follow = 180,
#   time_df = 5,
#   time_model = "factor"
# )
# ccw_overlap$risk_180
# ccw_overlap$contrasts
# compare_risk_estimates_to_truth(ccw_overlap, truth_mc)
#
# # Step 5: Cox regression as a secondary summary measure
# fit_ipcw_cox <- fit_ccw_cox(ccw_w, weight_var = "w_ipcw_trunc")
# summary(fit_ipcw_cox)
# exp(coef(fit_ipcw_cox))
# exp(confint(fit_ipcw_cox))
#
# # Model-implied 180-day risks from the Cox PH model
# cox_risk180 <- predict_risk180_from_cox(fit_ipcw_cox, horizon = 180)
# cox_risk180
# compare_risk_estimates_to_truth(cox_risk180, truth_mc)

# -----------------------------------------------------------------------------
# 9A) Support-focused workflow (recommended for package validation)
# -----------------------------------------------------------------------------
# source("ccw_post_mi_example_supported.R")
#
# sim <- simulate_post_mi(
#   n = 5000,
#   max_follow = 180,
#   grace = 30,
#   seed = 20260406,
#   day30_review_bump = 2.5  # backward-compatible legacy day-30 bump
# )
#
# ccw <- make_ccw_clones(sim$baseline, sim$daily, grace = 30, max_follow = 180)
# ccw_w <- add_ccw_weights(ccw, grace = 30, truncate = c(0.01, 0.99))
#
# diagnose_start_by_30_support(sim$baseline, ccw_w, grace = 30, weight_var = "w_ipcw")
#
# truth_mc <- monte_carlo_truth_post_mi(
#   n_mc = 200000,
#   max_follow = 180,
#   grace = 30,
#   chunk_size = 20000,
#   seed = 1,
#   day30_review_bump = 2.5  # backward-compatible legacy day-30 bump
# )
#
# ccw_weighted <- estimate_ccw_pooled_risk(
#   clone_long = ccw_w,
#   weight_var = "w_ipcw_trunc",
#   max_follow = 180,
#   time_model = "factor"
# )
#
# compare_risk_estimates_to_truth(ccw_weighted, truth_mc)


# -----------------------------------------------------------------------------
# 9) Repeated-simulation validation helper
# -----------------------------------------------------------------------------
# Important: compare_risk_estimates_to_truth() compares one simulated dataset to
# the Monte Carlo truth. Its "bias" column is replicate-level estimation error,
# not the estimator's expected bias across repeated samples. Use the helper below
# for actual validation studies.
run_ccw_validation_replicates <- function(n_reps = 200,
                                          n = 5000,
                                          max_follow = 180,
                                          grace = 30,
                                          base_seed = 20260406,
                                          truth_object = NULL,
                                          protocolized_day30 = TRUE,
                                          day30_review_bump = NULL,
                                          day30_lp_intercept = 0.20,
                                          day30_lp_L = 0.25,
                                          day30_lp_bleed = -0.35,
                                          day30_lp_age10 = -0.05,
                                          day30_lp_frailty = -0.05,
                                          day30_lp_large_mi = 0.15,
                                          truncate = c(0.01, 0.99)) {
  if (is.null(truth_object)) {
    truth_object <- monte_carlo_truth_post_mi(
      n_mc = 200000,
      max_follow = max_follow,
      grace = grace,
      chunk_size = 20000,
      seed = base_seed + 9999,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi
    )
  }

  out <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    sim <- simulate_post_mi(
      n = n,
      max_follow = max_follow,
      grace = grace,
      seed = base_seed + r,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi
    )

    ccw <- make_ccw_clones(sim$baseline, sim$daily, max_follow = max_follow, grace = grace)
    ccw_w <- add_ccw_weights(ccw, grace = grace, truncate = truncate)
    est <- estimate_ccw_pooled_risk(
      clone_long = ccw_w,
      weight_var = "w_ipcw_trunc",
      max_follow = max_follow,
      time_model = "factor"
    )

    cmp <- compare_risk_estimates_to_truth(est, truth_object)
    risk_cmp <- cmp$risk
    con_cmp <- cmp$contrasts
    out[[r]] <- data.frame(
      rep = r,
      est_no30 = risk_cmp$estimate[risk_cmp$arm == "no_start_by_30"],
      est_start = risk_cmp$estimate[risk_cmp$arm == "start_by_30"],
      bias_no30 = risk_cmp$bias[risk_cmp$arm == "no_start_by_30"],
      bias_start = risk_cmp$bias[risk_cmp$arm == "start_by_30"],
      rd_est = con_cmp$estimate[con_cmp$contrast == "risk_difference"],
      rd_bias = con_cmp$bias[con_cmp$contrast == "risk_difference"],
      rr_est = con_cmp$estimate[con_cmp$contrast == "risk_ratio"],
      rr_bias = con_cmp$bias[con_cmp$contrast == "risk_ratio"],
      stringsAsFactors = FALSE
    )
  }

  reps <- do.call(rbind, out)
  summary <- data.frame(
    metric = c("risk_no_start", "risk_start", "risk_difference", "risk_ratio"),
    mean_estimate = c(mean(reps$est_no30), mean(reps$est_start), mean(reps$rd_est), mean(reps$rr_est)),
    empirical_sd = c(sd(reps$est_no30), sd(reps$est_start), sd(reps$rd_est), sd(reps$rr_est)),
    mean_error = c(mean(reps$bias_no30), mean(reps$bias_start), mean(reps$rd_bias), mean(reps$rr_bias)),
    rmse = c(
      sqrt(mean(reps$bias_no30^2)),
      sqrt(mean(reps$bias_start^2)),
      sqrt(mean(reps$rd_bias^2)),
      sqrt(mean(reps$rr_bias^2))
    ),
    stringsAsFactors = FALSE
  )

  list(
    truth = truth_object,
    replicate_results = reps,
    summary = summary
  )
}

# -----------------------------------------------------------------------------
# 10) Convenience wrappers for one-shot execution
# -----------------------------------------------------------------------------
# These wrappers make it easy to:
# - run a single simulated dataset and compare estimates with Monte Carlo truth
# - inspect risk, risk difference, risk ratio, and estimate-minus-truth in one place
# - run repeated simulations to estimate actual finite-sample bias

run_single_ccw_validation <- function(n = 5000,
                                      max_follow = 180,
                                      grace = 30,
                                      sim_seed = 20260406,
                                      truth_seed = 1,
                                      truth_n_mc = 200000,
                                      truth_chunk_size = 20000,
                                      protocolized_day30 = TRUE,
                                      day30_review_bump = NULL,
                                      day30_lp_intercept = 0.20,
                                      day30_lp_L = 0.25,
                                      day30_lp_bleed = -0.35,
                                      day30_lp_age10 = -0.05,
                                      day30_lp_frailty = -0.05,
                                      day30_lp_large_mi = 0.15,
                                      truncate = c(0.01, 0.99),
                                      time_model = "factor",
                                      include_naive = TRUE,
                                      include_unweighted = TRUE,
                                      include_overlap = FALSE) {
  sim <- simulate_post_mi(
    n = n,
    max_follow = max_follow,
    grace = grace,
    seed = sim_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  truth_mc <- monte_carlo_truth_post_mi(
    n_mc = truth_n_mc,
    max_follow = max_follow,
    grace = grace,
    chunk_size = truth_chunk_size,
    seed = truth_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  ccw <- make_ccw_clones(sim$baseline, sim$daily, max_follow = max_follow, grace = grace)
  ccw_w <- add_ccw_weights(ccw, grace = grace, truncate = truncate)
  support <- diagnose_start_by_30_support(sim$baseline, ccw_w, grace = grace, weight_var = "w_ipcw")

  out <- list(
    sim = sim,
    truth = truth_mc,
    support = support
  )

  if (isTRUE(include_naive)) {
    naive_plr <- estimate_naive_pooled_risk(
      baseline = sim$baseline,
      daily = sim$daily,
      max_follow = max_follow,
      grace = grace,
      time_model = time_model
    )
    out$naive_pooled <- naive_plr
    out$naive_pooled_compare <- compare_risk_estimates_to_truth(naive_plr, truth_mc)
  }

  if (isTRUE(include_unweighted)) {
    ccw_unweighted <- estimate_ccw_pooled_risk(
      clone_long = ccw,
      weight_var = NULL,
      max_follow = max_follow,
      time_model = time_model
    )
    out$ccw_unweighted <- ccw_unweighted
    out$ccw_unweighted_compare <- compare_risk_estimates_to_truth(ccw_unweighted, truth_mc)
  }

  ccw_weighted <- estimate_ccw_pooled_risk(
    clone_long = ccw_w,
    weight_var = "w_ipcw_trunc",
    max_follow = max_follow,
    time_model = time_model
  )
  out$ccw_weighted <- ccw_weighted
  out$ccw_weighted_compare <- compare_risk_estimates_to_truth(ccw_weighted, truth_mc)

  if (isTRUE(include_overlap)) {
    ccw_overlap <- estimate_ccw_pooled_risk(
      clone_long = ccw_w,
      weight_var = "w_ow",
      max_follow = max_follow,
      time_model = time_model
    )
    out$ccw_overlap <- ccw_overlap
    out$ccw_overlap_compare <- compare_risk_estimates_to_truth(ccw_overlap, truth_mc)
  }

  make_method_rows <- function(method_name, cmp) {
    risk_rows <- data.frame(
      method = method_name,
      type = "risk",
      label = cmp$risk$arm,
      estimate = cmp$risk$estimate,
      truth = cmp$risk$truth,
      estimate_minus_truth = cmp$risk$bias,
      stringsAsFactors = FALSE
    )
    con_rows <- data.frame(
      method = method_name,
      type = "contrast",
      label = cmp$contrasts$contrast,
      estimate = cmp$contrasts$estimate,
      truth = cmp$contrasts$truth,
      estimate_minus_truth = cmp$contrasts$bias,
      stringsAsFactors = FALSE
    )
    rbind(risk_rows, con_rows)
  }

  summary_rows <- list()
  k <- 0L
  if (!is.null(out$naive_pooled_compare)) {
    k <- k + 1L; summary_rows[[k]] <- make_method_rows("naive_pooled", out$naive_pooled_compare)
  }
  if (!is.null(out$ccw_unweighted_compare)) {
    k <- k + 1L; summary_rows[[k]] <- make_method_rows("ccw_unweighted", out$ccw_unweighted_compare)
  }
  k <- k + 1L; summary_rows[[k]] <- make_method_rows("ccw_weighted", out$ccw_weighted_compare)
  if (!is.null(out$ccw_overlap_compare)) {
    k <- k + 1L; summary_rows[[k]] <- make_method_rows("ccw_overlap", out$ccw_overlap_compare)
  }

  out$summary_table <- do.call(rbind, summary_rows)
  out
}

print_single_ccw_validation <- function(x, digits = 4) {
  if (!is.list(x) || is.null(x$truth) || is.null(x$ccw_weighted_compare)) {
    stop("Please pass the object returned by run_single_ccw_validation().", call. = FALSE)
  }

  cat("\n=== Monte Carlo truth (180-day risk) ===\n")
  print(round(x$truth$risk_180, digits))

  cat("\n=== Support diagnostics for start_by_30 arm ===\n")
  print(round(x$support$natural_day30_support, digits))
  print(round(x$support$weighted_day30_diagnostics, digits))

  if (!is.null(x$summary_table)) {
    cat("\n=== Estimates, truth, and estimate-minus-truth ===\n")
    tbl <- x$summary_table
    num_cols <- vapply(tbl, is.numeric, logical(1))
    tbl[num_cols] <- lapply(tbl[num_cols], round, digits = digits)
    print(tbl)
  }

  invisible(x)
}

# -----------------------------------------------------------------------------
# 10) One-shot helper: run one full CCW analysis and compare with truth
# -----------------------------------------------------------------------------
run_one_ccw_analysis <- function(n = 5000,
                                 max_follow = 180,
                                 grace = 30,
                                 data_seed = 20260406,
                                 truth_seed = 1,
                                 n_mc = 200000,
                                 chunk_size = 20000,
                                 protocolized_day30 = TRUE,
                                 day30_review_bump = NULL,
                                 day30_lp_intercept = 0.20,
                                 day30_lp_L = 0.25,
                                 day30_lp_bleed = -0.35,
                                 day30_lp_age10 = -0.05,
                                 day30_lp_frailty = -0.05,
                                 day30_lp_large_mi = 0.15,
                                 truncate = c(0.01, 0.99),
                                 weight_var = "w_ipcw_trunc",
                                 time_model = "factor",
                                 include_naive = TRUE,
                                 include_unweighted = TRUE,
                                 include_cox = TRUE) {
  sim <- simulate_post_mi(
    n = n,
    max_follow = max_follow,
    grace = grace,
    seed = data_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  truth <- monte_carlo_truth_post_mi(
    n_mc = n_mc,
    max_follow = max_follow,
    grace = grace,
    chunk_size = chunk_size,
    seed = truth_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  ccw <- make_ccw_clones(sim$baseline, sim$daily, max_follow = max_follow, grace = grace)
  ccw_w <- add_ccw_weights(ccw, grace = grace, truncate = truncate)
  support <- diagnose_start_by_30_support(sim$baseline, ccw_w, grace = grace, weight_var = "w_ipcw")

  weighted_est <- estimate_ccw_pooled_risk(
    clone_long = ccw_w,
    weight_var = weight_var,
    max_follow = max_follow,
    time_model = time_model
  )
  weighted_cmp <- compare_risk_estimates_to_truth(weighted_est, truth)

  out <- list(
    args = list(
      n = n,
      max_follow = max_follow,
      grace = grace,
      data_seed = data_seed,
      truth_seed = truth_seed,
      n_mc = n_mc,
      chunk_size = chunk_size,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi,
      truncate = truncate,
      weight_var = weight_var,
      time_model = time_model
    ),
    support = support,
    truth = truth,
    weighted = list(
      estimate = weighted_est,
      compare = weighted_cmp
    ),
    sim = sim,
    ccw = ccw,
    ccw_w = ccw_w
  )

  if (isTRUE(include_naive)) {
    naive_est <- estimate_naive_pooled_risk(
      baseline = sim$baseline,
      daily = sim$daily,
      max_follow = max_follow,
      grace = grace,
      time_model = time_model
    )
    out$naive <- list(
      estimate = naive_est,
      compare = compare_risk_estimates_to_truth(naive_est, truth)
    )
  }

  if (isTRUE(include_unweighted)) {
    unweighted_est <- estimate_ccw_pooled_risk(
      clone_long = ccw,
      weight_var = NULL,
      max_follow = max_follow,
      time_model = time_model
    )
    out$unweighted <- list(
      estimate = unweighted_est,
      compare = compare_risk_estimates_to_truth(unweighted_est, truth)
    )
  }

  if (isTRUE(include_cox)) {
    cox_fit <- fit_ccw_cox(ccw_w, weight_var = weight_var)
    cox_risk180 <- predict_risk180_from_cox(cox_fit, horizon = max_follow)
    out$cox <- list(
      fit = cox_fit,
      risk_180 = cox_risk180,
      compare = compare_risk_estimates_to_truth(cox_risk180, truth)
    )
  }

  out
}

# -----------------------------------------------------------------------------
# 11) Compact printer for one-shot output
# -----------------------------------------------------------------------------
print_one_ccw_analysis <- function(x, digits = 4) {
  fmt <- function(obj) {
    if (is.data.frame(obj)) return(round_df(obj, digits = digits))
    obj
  }

  cat("\n=== Support diagnostics ===\n")
  print(fmt(x$support$natural_day30_support))
  print(fmt(x$support$weighted_day30_diagnostics))

  cat("\n=== Monte Carlo truth: 180-day risks ===\n")
  print(fmt(x$truth$risk_180))
  cat("\n=== Monte Carlo truth: contrasts ===\n")
  print(fmt(x$truth$contrasts))

  cat("\n=== Weighted CCW estimate ===\n")
  print(fmt(x$weighted$estimate$risk_180))
  print(fmt(x$weighted$estimate$contrasts))
  cat("\n=== Weighted CCW vs truth (one-run error; NOT repeated-sampling bias) ===\n")
  print(fmt(x$weighted$compare$risk))
  print(fmt(x$weighted$compare$contrasts))

  if (!is.null(x$naive)) {
    cat("\n=== Naive pooled logistic vs truth ===\n")
    print(fmt(x$naive$estimate$risk_180))
    print(fmt(x$naive$estimate$contrasts))
    print(fmt(x$naive$compare$risk))
    print(fmt(x$naive$compare$contrasts))
  }

  if (!is.null(x$unweighted)) {
    cat("\n=== Unweighted CCW vs truth ===\n")
    print(fmt(x$unweighted$estimate$risk_180))
    print(fmt(x$unweighted$estimate$contrasts))
    print(fmt(x$unweighted$compare$risk))
    print(fmt(x$unweighted$compare$contrasts))
  }

  if (!is.null(x$cox)) {
    cat("\n=== Cox-implied 180-day risks vs truth ===\n")
    print(fmt(x$cox$risk_180$risk_180))
    print(fmt(x$cox$risk_180$contrasts))
    print(fmt(x$cox$compare$risk))
    print(fmt(x$cox$compare$contrasts))
  }

  invisible(x)
}

# Utility: round numeric columns in a data.frame for compact printing
round_df <- function(df, digits = 4) {
  out <- df
  is_num <- vapply(out, is.numeric, logical(1))
  out[is_num] <- lapply(out[is_num], round, digits = digits)
  out
}

# -----------------------------------------------------------------------------
# 10) One-shot helper: run one full validation analysis and collect bias/RD/RR
# -----------------------------------------------------------------------------
# This is the easiest entry point when you want, in one call, to get:
# - Monte Carlo truth
# - support diagnostics
# - naive and CCW estimates
# - risk / risk-difference / risk-ratio comparisons to truth
#
# IMPORTANT:
# - The "bias" columns in the returned tables are replicate-level estimation
#   errors for THIS simulated dataset, not repeated-sampling expected bias.
# - For actual estimator bias, use run_ccw_validation_replicates().
# -----------------------------------------------------------------------------

run_post_mi_validation_once <- function(n = 5000,
                                        max_follow = 180,
                                        grace = 30,
                                        seed = 20260406,
                                        n_mc = 200000,
                                        truth_seed = 1,
                                        chunk_size = 20000,
                                        protocolized_day30 = TRUE,
                                        day30_review_bump = NULL,
                                        day30_lp_intercept = 0.20,
                                        day30_lp_L = 0.25,
                                        day30_lp_bleed = -0.35,
                                        day30_lp_age10 = -0.05,
                                        day30_lp_frailty = -0.05,
                                        day30_lp_large_mi = 0.15,
                                        truncate = c(0.01, 0.99),
                                        include_overlap = TRUE,
                                        include_cox = TRUE) {
  sim <- simulate_post_mi(
    n = n,
    max_follow = max_follow,
    grace = grace,
    seed = seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  truth_mc <- monte_carlo_truth_post_mi(
    n_mc = n_mc,
    max_follow = max_follow,
    grace = grace,
    chunk_size = chunk_size,
    seed = truth_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  ccw <- make_ccw_clones(
    baseline = sim$baseline,
    daily = sim$daily,
    grace = grace,
    max_follow = max_follow
  )

  ccw_w <- add_ccw_weights(
    clone_long = ccw,
    grace = grace,
    truncate = truncate
  )

  support <- diagnose_start_by_30_support(
    baseline = sim$baseline,
    weighted_clone_long = ccw_w,
    grace = grace,
    weight_var = "w_ipcw"
  )

  est_naive_logit <- predict_logistic180_risk(
    fit_object = fit_naive_logistic180(sim$baseline, grace = grace)
  )

  est_naive_pooled <- estimate_naive_pooled_risk(
    baseline = sim$baseline,
    daily = sim$daily,
    max_follow = max_follow,
    grace = grace,
    time_model = "factor"
  )

  est_ccw_unweighted <- estimate_ccw_pooled_risk(
    clone_long = ccw,
    weight_var = NULL,
    max_follow = max_follow,
    time_model = "factor"
  )

  est_ccw_weighted <- estimate_ccw_pooled_risk(
    clone_long = ccw_w,
    weight_var = "w_ipcw_trunc",
    max_follow = max_follow,
    time_model = "factor"
  )

  est_ccw_overlap <- NULL
  if (isTRUE(include_overlap)) {
    est_ccw_overlap <- estimate_ccw_pooled_risk(
      clone_long = ccw_w,
      weight_var = "w_ow",
      max_follow = max_follow,
      time_model = "factor"
    )
  }

  fit_cox <- NULL
  est_cox_risk180 <- NULL
  if (isTRUE(include_cox)) {
    fit_cox <- fit_ccw_cox(ccw_w, weight_var = "w_ipcw_trunc")
    est_cox_risk180 <- predict_risk180_from_cox(fit_cox, horizon = max_follow)
  }

  estimates <- list(
    naive_logit180 = est_naive_logit,
    naive_pooled = est_naive_pooled,
    ccw_unweighted = est_ccw_unweighted,
    ccw_weighted = est_ccw_weighted,
    ccw_overlap = est_ccw_overlap,
    cox_risk180 = est_cox_risk180
  )
  estimates <- estimates[!vapply(estimates, is.null, logical(1))]

  comparisons <- lapply(estimates, compare_risk_estimates_to_truth, truth_object = truth_mc)

  risk_summary <- do.call(
    rbind,
    lapply(names(comparisons), function(method_name) {
      x <- comparisons[[method_name]]$risk
      if (is.null(x) || nrow(x) == 0L) return(NULL)
      data.frame(
        method = method_name,
        arm = x$arm,
        estimate = x$estimate,
        truth = x$truth,
        bias = x$bias,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(risk_summary) <- NULL

  contrast_summary <- do.call(
    rbind,
    lapply(names(comparisons), function(method_name) {
      x <- comparisons[[method_name]]$contrasts
      if (is.null(x) || nrow(x) == 0L) return(NULL)
      data.frame(
        method = method_name,
        contrast = x$contrast,
        estimate = x$estimate,
        truth = x$truth,
        bias = x$bias,
        stringsAsFactors = FALSE
      )
    })
  )
  if (!is.null(contrast_summary)) {
    rownames(contrast_summary) <- NULL
  }

  list(
    meta = list(
      n = n,
      max_follow = max_follow,
      grace = grace,
      seed = seed,
      n_mc = n_mc,
      truth_seed = truth_seed,
      chunk_size = chunk_size,
      protocolized_day30 = protocolized_day30,
      day30_review_bump = day30_review_bump,
      day30_lp_intercept = day30_lp_intercept,
      day30_lp_L = day30_lp_L,
      day30_lp_bleed = day30_lp_bleed,
      day30_lp_age10 = day30_lp_age10,
      day30_lp_frailty = day30_lp_frailty,
      day30_lp_large_mi = day30_lp_large_mi,
      truncate = truncate
    ),
    support = support,
    truth = truth_mc,
    estimates = estimates,
    comparisons = comparisons,
    risk_summary = risk_summary,
    contrast_summary = contrast_summary,
    fit_cox = fit_cox,
    sim = sim,
    ccw = ccw,
    ccw_w = ccw_w
  )
}

print_post_mi_validation_once <- function(result_object,
                                          digits = 4) {
  stopifnot(is.list(result_object))

  cat("\n=== Support diagnostics ===\n")
  print(round_df_numeric(result_object$support$natural_day30_support, digits = digits), row.names = FALSE)
  if (!is.null(result_object$support$weighted_day30_diagnostics)) {
    print(round_df_numeric(result_object$support$weighted_day30_diagnostics, digits = digits), row.names = FALSE)
  }

  cat("\n=== Monte Carlo truth: 180-day risks ===\n")
  print(round_df_numeric(result_object$truth$risk_180, digits = digits), row.names = FALSE)

  cat("\n=== Monte Carlo truth: contrasts ===\n")
  print(round_df_numeric(result_object$truth$contrasts, digits = digits), row.names = FALSE)

  cat("\n=== Method-by-arm comparison to truth ===\n")
  print(round_df_numeric(result_object$risk_summary, digits = digits), row.names = FALSE)

  cat("\n=== Method-by-contrast comparison to truth ===\n")
  print(round_df_numeric(result_object$contrast_summary, digits = digits), row.names = FALSE)

  invisible(result_object)
}


round_df_numeric <- function(x, digits = 4) {
  if (is.null(x)) return(x)
  out <- x
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  out
}


# -----------------------------------------------------------------------------
# 10) Convenience wrappers: one-shot analysis and validation summary
# -----------------------------------------------------------------------------

run_ccw_one_shot <- function(n = 5000,
                             n_mc = 200000,
                             max_follow = 180,
                             grace = 30,
                             sim_seed = 20260406,
                             truth_seed = 1,
                             protocolized_day30 = TRUE,
                             day30_review_bump = NULL,
                             day30_lp_intercept = 0.20,
                             day30_lp_L = 0.25,
                             day30_lp_bleed = -0.35,
                             day30_lp_age10 = -0.05,
                             day30_lp_frailty = -0.05,
                             day30_lp_large_mi = 0.15,
                             truncate = c(0.01, 0.99),
                             include_naive = TRUE,
                             include_unweighted = TRUE,
                             include_overlap = FALSE,
                             include_cox = TRUE) {
  sim <- simulate_post_mi(
    n = n,
    max_follow = max_follow,
    grace = grace,
    seed = sim_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  ccw <- make_ccw_clones(
    baseline = sim$baseline,
    daily = sim$daily,
    grace = grace,
    max_follow = max_follow
  )

  ccw_w <- add_ccw_weights(
    clone_long = ccw,
    grace = grace,
    truncate = truncate
  )

  support <- diagnose_start_by_30_support(
    baseline = sim$baseline,
    weighted_clone_long = ccw_w,
    grace = grace,
    weight_var = "w_ipcw"
  )

  truth_mc <- monte_carlo_truth_post_mi(
    n_mc = n_mc,
    max_follow = max_follow,
    grace = grace,
    chunk_size = 20000,
    seed = truth_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  out <- list(
    support = support,
    truth = truth_mc,
    objects = list(sim = sim, ccw = ccw, ccw_w = ccw_w),
    comparisons = list()
  )

  if (isTRUE(include_naive)) {
    naive_plr <- estimate_naive_pooled_risk(
      baseline = sim$baseline,
      daily = sim$daily,
      max_follow = max_follow,
      grace = grace,
      time_model = "factor"
    )
    out$objects$naive_plr <- naive_plr
    out$comparisons$naive <- compare_risk_estimates_to_truth(naive_plr, truth_mc)
  }

  if (isTRUE(include_unweighted)) {
    ccw_unweighted <- estimate_ccw_pooled_risk(
      clone_long = ccw,
      weight_var = NULL,
      max_follow = max_follow,
      time_model = "factor"
    )
    out$objects$ccw_unweighted <- ccw_unweighted
    out$comparisons$ccw_unweighted <- compare_risk_estimates_to_truth(ccw_unweighted, truth_mc)
  }

  ccw_weighted <- estimate_ccw_pooled_risk(
    clone_long = ccw_w,
    weight_var = "w_ipcw_trunc",
    max_follow = max_follow,
    time_model = "factor"
  )
  out$objects$ccw_weighted <- ccw_weighted
  out$comparisons$ccw_weighted <- compare_risk_estimates_to_truth(ccw_weighted, truth_mc)

  if (isTRUE(include_overlap)) {
    ccw_overlap <- estimate_ccw_pooled_risk(
      clone_long = ccw_w,
      weight_var = "w_ow",
      max_follow = max_follow,
      time_model = "factor"
    )
    out$objects$ccw_overlap <- ccw_overlap
    out$comparisons$ccw_overlap <- compare_risk_estimates_to_truth(ccw_overlap, truth_mc)
  }

  if (isTRUE(include_cox)) {
    fit_ipcw_cox <- fit_ccw_cox(ccw_w, weight_var = "w_ipcw_trunc")
    cox_risk180 <- predict_risk180_from_cox(fit_ipcw_cox, horizon = max_follow)
    out$objects$fit_ipcw_cox <- fit_ipcw_cox
    out$objects$cox_risk180 <- cox_risk180
    out$comparisons$cox <- compare_risk_estimates_to_truth(cox_risk180, truth_mc)
  }

  out
}

run_ccw_validation_summary <- function(n_reps = 100,
                                       n = 5000,
                                       n_mc = 200000,
                                       max_follow = 180,
                                       grace = 30,
                                       sim_seed = 20260406,
                                       truth_seed = 1,
                                       protocolized_day30 = TRUE,
                                       day30_review_bump = NULL,
                                       day30_lp_intercept = 0.20,
                                       day30_lp_L = 0.25,
                                       day30_lp_bleed = -0.35,
                                       day30_lp_age10 = -0.05,
                                       day30_lp_frailty = -0.05,
                                       day30_lp_large_mi = 0.15,
                                       truncate = c(0.01, 0.99)) {
  truth_mc <- monte_carlo_truth_post_mi(
    n_mc = n_mc,
    max_follow = max_follow,
    grace = grace,
    chunk_size = 20000,
    seed = truth_seed,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi
  )

  val <- run_ccw_validation_replicates(
    n_reps = n_reps,
    n = n,
    max_follow = max_follow,
    grace = grace,
    base_seed = sim_seed,
    truth_object = truth_mc,
    protocolized_day30 = protocolized_day30,
    day30_review_bump = day30_review_bump,
    day30_lp_intercept = day30_lp_intercept,
    day30_lp_L = day30_lp_L,
    day30_lp_bleed = day30_lp_bleed,
    day30_lp_age10 = day30_lp_age10,
    day30_lp_frailty = day30_lp_frailty,
    day30_lp_large_mi = day30_lp_large_mi,
    truncate = truncate
  )

  val
}



# 한 번의 simulated dataset + Monte Carlo truth + naive/CCW 결과를 한 번에 생성
res_once <- run_post_mi_validation_once(
  n = 5000,
  max_follow = 180,
  grace = 30,
  seed = 20260406,
  n_mc = 200000,
  truth_seed = 1,
  chunk_size = 20000,
  protocolized_day30 = TRUE,
  day30_review_bump = NULL,
  day30_lp_intercept = 0.20,
  day30_lp_L = 0.25,
  day30_lp_bleed = -0.35,
  day30_lp_age10 = -0.05,
  day30_lp_frailty = -0.05,
  day30_lp_large_mi = 0.15,
  truncate = c(0.01, 0.99),
  include_overlap = FALSE,
  include_cox = TRUE
)

# 보기 좋게 출력
print_post_mi_validation_once(res_once)

# 핵심 표만 따로 보기
res_once$truth$risk_180
res_once$truth$contrasts

res_once$risk_summary
res_once$contrast_summary

# weighted CCW만 따로 보기
subset(res_once$risk_summary, method == "ccw_weighted")
subset(
  res_once$contrast_summary,
  method == "ccw_weighted" &
    contrast %in% c("risk_difference", "risk_ratio")
)

# support 진단도 같이 확인
res_once$support$natural_day30_support
res_once$support$weighted_day30_diagnostics


### 아래 코드는 반복
# source("ccw_post_mi_example_protocolized.R")

# val <- run_ccw_validation_replicates(
#   n_reps = 100,
#   n = 5000,
#   max_follow = 180,
#   grace = 30,
#   base_seed = 20260406,
#   protocolized_day30 = TRUE,
#   day30_review_bump = NULL,
#   day30_lp_intercept = 0.20,
#   day30_lp_L = 0.25,
#   day30_lp_bleed = -0.35,
#   day30_lp_age10 = -0.05,
#   day30_lp_frailty = -0.05,
#   day30_lp_large_mi = 0.15,
#   truncate = c(0.01, 0.99)
# )

# val$truth$risk_180
# val$truth$contrasts
# val$summary
# head(val$replicate_results)

# # Weighted CCW 간단하게

# source("ccw_post_mi_example_protocolized.R")

# res_once <- run_post_mi_validation_once(
#   n = 5000,
#   seed = 20260406,
#   n_mc = 200000,
#   truth_seed = 1,
#   protocolized_day30 = TRUE,
#   day30_review_bump = NULL,
#   day30_lp_intercept = 0.20,
#   day30_lp_L = 0.25,
#   day30_lp_bleed = -0.35,
#   day30_lp_age10 = -0.05,
#   day30_lp_frailty = -0.05,
#   day30_lp_large_mi = 0.15,
#   include_overlap = FALSE,
#   include_cox = FALSE
# )

# subset(res_once$risk_summary, method == "ccw_weighted")
# subset(
#   res_once$contrast_summary,
#   method == "ccw_weighted" &
#     contrast %in% c("risk_difference", "risk_ratio")
# )

# # compute_lp_start_post_mi 가 제대로 바뀌어는지 테스트
# source("ccw_post_mi_example_protocolized.R")

# run_tests_compute_lp_start_post_mi()
# run_day30_interface_smoke_tests()

# compute_lp_start_post_mi(
#   day = 30,
#   L_now = c(-1, 0, 1),
#   age10 = 0,
#   bleed = 0,
#   frailty = 0,
#   large_mi = 0,
#   day30_mode = "protocolized"
# )

# compute_lp_start_post_mi(
#   day = 30,
#   L_now = c(-1, 0, 1),
#   age10 = 0,
#   bleed = 0,
#   frailty = 0,
#   large_mi = 0,
#   day30_review_bump = 2.5
# )

# show_compute_lp_start_post_mi_examples()

# # 가장 중요한 출력
# res_once$support$natural_day30_support
# res_once$support$weighted_day30_diagnostics
# subset(res_once$risk_summary, method == "ccw_weighted")
# subset(res_once$contrast_summary, method == "ccw_weighted")
