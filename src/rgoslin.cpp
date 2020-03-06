#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <cppgoslin/cppgoslin/cppgoslin.h>

using namespace Rcpp;

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
bool rcpp_is_valid_lipid_name(std::string lipid_name) {
    /* create instance of lipid parser containing several grammars */
    LipidParser lipid_parser;

    /* parsing lipid name into a lipid container data structure */
    LipidAdduct* lipidAdduct = lipid_parser.parse(lipid_name);

    /* checking if parsing was successful, otherwise lipid reference remains NULL */
    bool isValidLipidName = false;
    if (lipidAdduct != NULL){
        isValidLipidName = true;
    }

    return isValidLipidName;
}

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
SEXP rcpp_parse_lipid_name(std::string lipid_name) {
    /* create instance of lipid parser containing several grammars */
    LipidParser lipid_parser;
    
    /* parsing lipid name into a lipid container data structure */
    LipidAdduct* lipidAdduct = lipid_parser.parse(lipid_name);
    std::map<std::string, std::string> lipidDetails;
    if (lipidAdduct != NULL){
        LipidSpeciesInfo info = lipidAdduct->lipid->info;
        std::string nativeLevelName = lipidAdduct->get_lipid_string(info.level);
        std::string category = lipidAdduct->get_lipid_string(CATEGORY);
        std::string species = lipidAdduct->get_lipid_string(SPECIES);
        lipidDetails["originalName"] = lipid_name;
        lipidDetails["nativeLevelName"] = nativeLevelName;
        lipidDetails["category"] = category;
        lipidDetails["species"] = species;
    }
    return Rcpp::wrap(lipidDetails);
}