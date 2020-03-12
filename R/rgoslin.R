#' @useDynLib rgoslin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Check lipid name.
#'
#' \code{isValidLipidName} checks the provided lipid name against the built-in grammars.
#' Will return FALSE if none of the parsers was able to parse the provided name successfully or if any error was raised.
#' @param lipidName The lipid name to check.
#' @examples 
#' isValidLipidName("PC 32:1")
#' isValidLipidName("PC(32:1)")
#' isValidLipidName("PCX(32:1)")
#' @return TRUE if the lipidName could be parsed, FALSE otherwise.
#' @export
isValidLipidName <- function(lipidName) {
  tryCatch(
    {
      return(rcpp_is_valid_lipid_name(lipidName))
    }, error = function(err) {
      message(paste("Could not parse the provided lipid name", lipidName," with any of the available parsers!"))
      return(FALSE)
    }
  )
}

#' Parse lipid name.
#'
#' \code{parseLipidName} reads the provided lipid name and returns structural information as a named vector.
#' Will return a vector with the "Grammar" element set to "NOT_PARSEABLE" if none of the parsers was able to parse the provided name successfully.
#' If any error was raised, returns an empty vector.
#' @param lipidName The lipid name to parse.
#' @examples 
#' parseLipidName("PC 32:1")
#' parseLipidName("LPC 18:1_16:0") 
#' parseLipidName("TG(18:1_18:0_16:1)")
#' @return Named vector with details of the lipid, empty vector otherwise.
#' @export
parseLipidName <- function(lipidName) {
  tryCatch(
    {
      return(rcpp_parse_lipid_name(lipidName))
    }, error = function(err) {
      message(paste("Could not parse the provided lipid name", lipidName," with any of the available parsers!"))
      return(c())
    }
  )
}