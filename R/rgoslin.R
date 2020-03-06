#' @useDynLib rgoslin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Check the provided lipid name against the built-in grammars.
#' Will return FALSE if none of the parsers was able to parse the provided name successfully.
#' @param lipidName the lipid name to check.
#' @export
isValidLipidName <- function(lipidName) {
  return(rcpp_is_valid_lipid_name(lipidName))
}

#' Parse the provided lipid name and return structural information.
#' Will return an empty vector if none of the parsers was able to parse the provided name successfully.
#' @param lipidName the lipid name to check.
#' @export
parseLipidName <- function(lipidName) {
  return(rcpp_parse_lipid_name(lipidName))
}