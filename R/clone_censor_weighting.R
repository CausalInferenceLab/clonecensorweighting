#' Create a starter clone-censor-weighting dataset
#'
#' Expands each observation into one row per treatment strategy and flags
#' immediate deviations from the observed baseline treatment as censored clones.
#' This is a light-weight scaffold for package development rather than a full
#' causal inference implementation.
#'
#' @param data A data frame containing one row per participant.
#' @param id The name of the participant identifier column.
#' @param follow_up The name of the follow-up time column.
#' @param event The name of the event indicator column.
#' @param treatment The name of the observed treatment column.
#' @param regimes Optional vector of treatment strategies to clone. When `NULL`,
#'   unique non-missing observed treatments are used.
#'
#' @return A tibble with one row per participant-strategy combination and the
#'   additional columns `.clone_id`, `.regime`, `.censored`, and `.weight`.
#' @export
clone_censor_weighting <- function(
    data,
    id,
    follow_up,
    event,
    treatment,
    regimes = NULL
) {
  .assert_data_frame(data)

  required_columns <- c(id, follow_up, event, treatment)
  .assert_required_columns(data, required_columns)

  tbl <- tibble::as_tibble(data)

  if (!is.numeric(tbl[[follow_up]])) {
    stop("`follow_up` must refer to a numeric column.", call. = FALSE)
  }

  if (!all(stats::na.omit(tbl[[event]]) %in% c(0, 1))) {
    stop("`event` must contain only 0/1 values.", call. = FALSE)
  }

  if (is.null(regimes)) {
    regimes <- sort(unique(stats::na.omit(tbl[[treatment]])))
  }

  if (length(regimes) == 0) {
    stop("`regimes` must contain at least one treatment strategy.", call. = FALSE)
  }

  id_values <- tbl[[id]]
  observed_treatment <- tidyr::replace_na(as.character(tbl[[treatment]]), ".missing")
  regime_values <- as.character(regimes)

  clones <- purrr::map(
    regime_values,
    function(regime) {
      clone <- dplyr::mutate(
        tbl,
        .clone_id = paste(id_values, regime, sep = "::"),
        .regime = regime,
        .censored = as.integer(observed_treatment != regime),
        .weight = 1
      )

      clone[, c(
        ".clone_id",
        ".regime",
        ".censored",
        ".weight",
        setdiff(names(clone), c(".clone_id", ".regime", ".censored", ".weight"))
      )]
    }
  )

  cloned_tbl <- dplyr::bind_rows(clones)
  row_order <- order(cloned_tbl[[id]], cloned_tbl[[follow_up]], cloned_tbl[[".regime"]])

  tibble::as_tibble(cloned_tbl[row_order, , drop = FALSE])
}

#' Construct a survival response
#'
#' @param data A data frame with follow-up and event columns.
#' @param follow_up The name of the follow-up time column.
#' @param event The name of the event indicator column.
#'
#' @return An object of class `"Surv"`.
#' @export
make_surv_response <- function(data, follow_up, event) {
  .assert_data_frame(data)
  .assert_required_columns(data, c(follow_up, event))

  survival::Surv(time = data[[follow_up]], event = data[[event]])
}
