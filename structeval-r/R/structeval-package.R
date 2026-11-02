#' structeval: structure-aware evaluation of extracted fields
#'
#' Scores extracted fields against reference labels ([read_oncoshot()],
#' [score()]), reports metrics with intervals matched to each field's
#' structure ([metrics()]), tests a threshold, a paired difference and the
#' balance of omissions and commissions ([test_threshold()],
#' [test_paired()], [test_asymmetry()]), and plans sample sizes
#' ([n_threshold()], [n_paired()], [min_certifiable()], the `power_*()`
#' functions and [simulate_counts()]).
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(".data")
