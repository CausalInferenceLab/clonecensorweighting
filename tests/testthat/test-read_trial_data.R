test_that("read_trial_data reads csv input into a tibble", {
  csv_file <- tempfile(fileext = ".csv")
  writeLines(
    c(
      "id,follow_up,event,treatment",
      "1,10,1,A",
      "2,12,0,B"
    ),
    con = csv_file
  )

  result <- read_trial_data(csv_file)

  expect_s3_class(result, "tbl_df")
  expect_equal(names(result), c("id", "follow_up", "event", "treatment"))
  expect_equal(nrow(result), 2)
})
