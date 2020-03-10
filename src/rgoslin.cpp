#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <cppgoslin/cppgoslin/cppgoslin.h>
#include <typeinfo>

using namespace Rcpp;

/**
 * Join elements of a string vector with the provided delimiter
 */
std::string join(std::vector<std::string>& vector, const char* delimiter) {
    std::ostringstream out;
    if (!vector.empty())
    {
        std::copy(vector.begin(), vector.end() - 1, std::ostream_iterator<string>(out, delimiter));
        out << vector.back();
    }
    return out.str();
}

/**
 *  Return a textual representation for the provided lipid level
 */
std::string get_lipid_level_str(LipidLevel& level) {
    switch (level){
        case UNDEFINED_LEVEL:
            return "UNDEFINED";
        case CATEGORY:
            return "CATEGORY";
        case CLASS:
            return "CLASS";
        case SPECIES:
            return "SPECIES";
        case MOLECULAR_SUBSPECIES:
            return "MOLECULAR_SUBSPECIES";
        case STRUCTURAL_SUBSPECIES:
            return "STRUCTURAL_SUBSPECIES";
        case ISOMERIC_SUBSPECIES:
            return "ISOMERIC_SUBSPECIES";
        default:
            return "UNDEFINED";
    }
}

/** 
 * Uses the cppGoslin library's LipidParser->parse(lipid_name) method to check,
 * whether the provided name is valid in any of the supported parsers/grammars.
 * Returns true if at least one parser supports the lipid name, false if no parser 
 * supports the lipid name.
 */
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
bool rcpp_is_valid_lipid_name(std::string lipid_name) {
    /* create instance of lipid parser containing several grammars */
    LipidParser lipid_parser;
    LipidAdduct* lipidAdduct;
    /* parsing lipid name into a lipid container data structure */
    try {
        lipidAdduct = lipid_parser.parse(lipid_name);
        /* checking if parsing was successful, otherwise lipid reference remains NULL */
        bool isValidLipidName = lipidAdduct != NULL;
        if(lipidAdduct) {
            delete lipidAdduct;
        }
        return isValidLipidName;
    } catch (LipidException &e){
        warning("Parsing of lipid name '" +lipid_name+"' caused an exception: "+ e.what());
        return false;
    }
}

/** 
 * Reimplements the cppGoslin library's LipidParser->parse(lipid_name) method.
 * It adds information about the parser/grammar that was used to parse the lipid name,
 * if parsing was successful.
 */
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
SEXP rcpp_parse_lipid_name(std::string lipid_name) {
    CharacterVector lipidDetails;
    /* parsing lipid name into a lipid container data structure */
    LipidAdduct* lipidAdduct;
    try {
        /* create instance of lipid parser containing several grammars */
        LipidParser lipid_parser;
        std::string grammar = "NA";
        for (auto parser : lipid_parser.parser_list) {
            try {
                LipidAdduct *lipid = parser->parse(lipid_name);
                if (lipid){
                    lipidAdduct = lipid;
                    grammar = parser->grammar_name;
                    break;
                }
            } catch (LipidException &e){
                forward_exception_to_r(e);
            }
        }

        lipidDetails["Normalized Name"] = "NA";
        lipidDetails["Original Name"] = "NA";
        lipidDetails["Grammar"] = "NA";
        lipidDetails["Lipid Maps Category"] = "NA";
        lipidDetails["Lipid Maps Main Class"] = "NA";
        lipidDetails["Functional Class Abbr"] = "NA";
        lipidDetails["Functional Class Synonyms"] = "NA";
        lipidDetails["Level"] = "NA";
        lipidDetails["Total C"] = "NA";
        lipidDetails["Total OH"] = "NA";
        lipidDetails["Total DB"] = "NA";
        lipidDetails["FA1 Position"] = "NA";
        lipidDetails["FA1 C"] = "NA";
        lipidDetails["FA1 OH"] = "NA";
        lipidDetails["FA1 DB"] = "NA";
        lipidDetails["FA1 Bond Type"] = "NA";
        lipidDetails["FA2 Position"] = "NA";
        lipidDetails["FA2 C"] = "NA";
        lipidDetails["FA2 OH"] = "NA";
        lipidDetails["FA2 DB"] = "NA";
        lipidDetails["FA2 Bond Type"] = "NA";
        lipidDetails["LCB Position"] = "NA";
        lipidDetails["LCB C"] = "NA";
        lipidDetails["LCB OH"] = "NA";
        lipidDetails["LCB DB"] = "NA";
        lipidDetails["LCB Bond Type"] = "NA";
        lipidDetails["FA3 Position"] = "NA";
        lipidDetails["FA3 C"] = "NA";
        lipidDetails["FA3 OH"] = "NA";
        lipidDetails["FA3 DB"] = "NA";
        lipidDetails["FA3 Bond Type"] = "NA";
        lipidDetails["FA4 Position"] = "NA";
        lipidDetails["FA4 C"] = "NA";
        lipidDetails["FA4 OH"] = "NA";
        lipidDetails["FA4 DB"] = "NA";
        lipidDetails["FA4 Bond Type"] = "NA";
    
        if (lipidAdduct){
            std::string grammar = "NA";
            std::string originalName = "NA";
            std::string nativeLevelName = "NA";
            std::string lipidMapsCategory = "NA";
            std::string lipidMapsMainClass = "NA";
            std::string species = "NA";
            std::string headGroup = "NA";
            std::string headGroupSynonyms = "NA";
            std::string level = "NA";
            std::string totalC = "NA";
            std::string totalOH = "NA";
            std::string totalDB = "NA";
            
            lipidMapsCategory = lipidAdduct->get_lipid_string(CATEGORY);
            lipidMapsMainClass = lipidAdduct->get_lipid_string(CLASS);
            species = lipidAdduct->get_lipid_string(SPECIES);
            if(lipidAdduct->lipid) {
                LipidSpecies* lipid = lipidAdduct->lipid;
                Rcout << "Lipid object is defined" << "\n";
                LipidSpeciesInfo info = (*lipid).info;
                Rcout << "Lipid species info is defined" << "\n";
                // nativeLevelName = lipidAdduct->get_lipid_string(info.level);
                // Rcout << "Lipid string is defined" << "\n";
                // lipidMapsMainClass = lipidAdduct->lipid->get_class_name();
                // Rcout << "Lipid maps main class is defined" << "\n";
                // headGroup = lipid->head_group;
                // Rcout << "Lipid head group is defined" << "\n";
                // LipidClassMeta lcMeta = lipid_classes.at(lipid->get_class(headGroup));
                // Rcout << "Lipid class meta is defined" << "\n";
                // headGroupSynonyms = "[" + join(lcMeta.synonyms, ", ") + "]";
                // level = get_lipid_level_str(info.level);
                // Rcout << "Lipid level is defined" << "\n";
                // totalC = info.num_carbon;
                // Rcout << "Lipid total carbon is defined" << "\n";
                // totalOH = info.num_hydroxyl;
                // Rcout << "Lipid total hydroxyl is defined" << "\n";
                // totalDB = info.num_double_bonds;
                // Rcout << "Lipid total double bonds is defined" << "\n";
            }
    //         lipidDetails["Species Name"] = species;
    //         lipidDetails["Functional Class Abbr"] = "[" + headGroup + "]";
    //         LipidClassMeta lcMeta = lipid_classes.at(lipidAdduct->lipid->get_class(headGroup));
    //         std:string synonyms = "[" + join(lcMeta.synonyms, ", ") + "]";
    //         lipidDetails["Functional Class Synonyms"] = synonyms;
    //         lipidDetails["Level"] = get_lipid_level_str(info.level);
    //         lipidDetails["Total C"] = info.num_carbon;
    //         lipidDetails["Total OH"] = info.num_hydroxyl;
    //         lipidDetails["Total DB"] = info.num_double_bonds;
    //             //species;
    // // Normalized Name	Original Name	Grammar	Lipid Maps Category	Lipid Maps Main Class	Functional Class Abbr	Functional Class Synonyms	Level	Total #C	Total #OH	Total #DB	FA1 Position	FA1 #C	FA1 #OH	FA1 #DB	FA1 Bond Type	FA2 Position	FA2 #C	FA2 #OH	FA2 #DB	FA2 Bond Type	LCB Position	LCB #C	LCB #OH	LCB #DB	LCB Bond Type	FA3 Position	FA3 #C	FA3 #OH	FA3 #DB	FA3 Bond Type	FA4 Position	FA4 #C	FA4 #OH	FA4 #DB	FA4 Bond Type        
            lipidDetails["Normalized Name"] = nativeLevelName;
            lipidDetails["Original Name"] = lipid_name;
            lipidDetails["Grammar"] = grammar;
            lipidDetails["Lipid Maps Category"] = lipidMapsCategory;
            lipidDetails["Lipid Maps Main Class"] = lipidMapsMainClass;
            lipidDetails["Functional Class Abbr"] = "[" + headGroup + "]";
            lipidDetails["Functional Class Synonyms"] = headGroupSynonyms;
            lipidDetails["Level"] = headGroupSynonyms;
            lipidDetails["Total C"] = totalC;
            lipidDetails["Total OH"] = totalOH;
            lipidDetails["Total DB"] = totalDB;
            
            
        } else {
            lipidDetails["Normalized Name"] = "NA";
            lipidDetails["Original Name"] = lipid_name;
            lipidDetails["Grammar"] = "NA";
        }
        if(lipidAdduct) {
            delete lipidAdduct;
        }
        return lipidDetails;
    } catch(LipidException &e) {
        warning("Parsing of lipid name '" +lipid_name+"' caused an exception: "+ e.what());
        return lipidDetails;
    }
}