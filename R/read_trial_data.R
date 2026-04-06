#' Read trial-style data from CSV
#'
#' @param file Path to a CSV file.
#' @param show_col_types Passed to [readr::read_csv()].
#'
#' @return A tibble.
#' @export
read_trial_data <- function(file, show_col_types = FALSE) {
  tibble::as_tibble(
    readr::read_csv(file, show_col_types = show_col_types)
  )
}
