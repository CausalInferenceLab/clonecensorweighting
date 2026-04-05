test_that("clone_censor_weighting expands rows across regimes", {
  trial_data <- tibble::tibble(
    id = c(1, 2),
    follow_up = c(10, 12),
    event = c(1, 0),
    treatment = c("A", "B")
  )

  result <- clone_censor_weighting(
    data = trial_data,
    id = "id",
    follow_up = "follow_up",
    event = "event",
    treatment = "treatment",
    regimes = c("A", "B")
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 4)
  expect_equal(result$.censored, c(0L, 1L, 1L, 0L))
  expect_equal(result$.weight, rep(1, 4))
})

test_that("make_surv_response returns a Surv object", {
  trial_data <- tibble::tibble(
    follow_up = c(5, 8),
    event = c(1, 0)
  )

  surv_response <- make_surv_response(
    data = trial_data,
    follow_up = "follow_up",
    event = "event"
  )

  expect_s3_class(surv_response, "Surv")
  expect_equal(length(surv_response), 2)
})
