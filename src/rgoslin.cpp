#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <cppgoslin.h>
#include <typeinfo>

using namespace Rcpp;
using namespace goslin;

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
 * Returns the lipid name for the requested level or NA_STRING if that is not possible, or lipidAdduct is FALSE.
 */
std::string get_lipid_name_for_level_with_warnings(LipidAdduct* lipidAdduct, LipidLevel level, bool withWarning) {
    String chr_na = NA_STRING;
    if(lipidAdduct) {
        LipidSpecies* lipid = lipidAdduct->lipid;
        LipidSpeciesInfo info = (*lipid).info;
        try {
            return lipidAdduct->get_lipid_string(level);
        } catch(LipidException &e) {
            if(withWarning) {
                warning("Lipid '"+lipidAdduct->get_lipid_string(info.level)+"' with native level '"+get_lipid_level_str(info.level)+"' can not generate name for more specific level '" +get_lipid_level_str(level)+"'!");
            }
            return chr_na;
        }
    }
    return chr_na;
}

/**
 * Returns the lipid name for the requested level or NA_STRING if that is not possible, or lipidAdduct is FALSE.
 */
std::string get_lipid_name_for_level(LipidAdduct* lipidAdduct, LipidLevel level) {
    return get_lipid_name_for_level_with_warnings(lipidAdduct, level, false);
}

/**
 * Creates the CharacterVector for the provided adduct. Will be filled with default values,
 * if the lipidAdduct is false. The grammar argument must be the name of a valid grammar, or 
 * the NA_STRING.
 */
SEXP handle_lipid(LipidAdduct* lipidAdduct, std::string lipid_name, std::string grammar) {
    CharacterVector lipidDetails;
    String chr_na = NA_STRING;
    lipidDetails["Normalized Name"] = chr_na;
    lipidDetails["Original Name"] = chr_na;
    lipidDetails["Grammar"] = chr_na;
    lipidDetails["Lipid Maps Category"] = chr_na;
    lipidDetails["Lipid Maps Main Class"] = chr_na;
    lipidDetails["Species Name"] = chr_na;
    lipidDetails["Molecular Subspecies Name"] = chr_na;
    lipidDetails["Structural Subspecies Name"] = chr_na;
    lipidDetails["Isomeric Subspecies Name"] = chr_na;
    lipidDetails["Functional Class Abbr"] = chr_na;
    lipidDetails["Functional Class Synonyms"] = chr_na;
    lipidDetails["Level"] = chr_na;
    lipidDetails["Total C"] = chr_na;
    lipidDetails["Total OH"] = chr_na;
    lipidDetails["Total DB"] = chr_na;
    lipidDetails["Mass"] = chr_na;
    lipidDetails["Sum Formula"] = chr_na;
    lipidDetails["FA1 Position"] = chr_na;
    lipidDetails["FA1 C"] = chr_na;
    lipidDetails["FA1 OH"] = chr_na;
    lipidDetails["FA1 DB"] = chr_na;
    lipidDetails["FA1 Bond Type"] = chr_na;
    lipidDetails["FA2 Position"] = chr_na;
    lipidDetails["FA2 C"] = chr_na;
    lipidDetails["FA2 OH"] = chr_na;
    lipidDetails["FA2 DB"] = chr_na;
    lipidDetails["FA2 Bond Type"] = chr_na;
    lipidDetails["LCB Position"] = chr_na;
    lipidDetails["LCB C"] = chr_na;
    lipidDetails["LCB OH"] = chr_na;
    lipidDetails["LCB DB"] = chr_na;
    lipidDetails["LCB Bond Type"] = chr_na;
    lipidDetails["FA3 Position"] = chr_na;
    lipidDetails["FA3 C"] = chr_na;
    lipidDetails["FA3 OH"] = chr_na;
    lipidDetails["FA3 DB"] = chr_na;
    lipidDetails["FA3 Bond Type"] = chr_na;
    lipidDetails["FA4 Position"] = chr_na;
    lipidDetails["FA4 C"] = chr_na;
    lipidDetails["FA4 OH"] = chr_na;
    lipidDetails["FA4 DB"] = chr_na;
    lipidDetails["FA4 Bond Type"] = chr_na;
    
    if (lipidAdduct){
        // grammar = parser->grammar_name;
        std::string originalName = chr_na;
        std::string nativeLevelName = chr_na;
        std::string lipidMapsCategory = chr_na;
        std::string lipidMapsMainClass = chr_na;
        std::string species = chr_na;
        std::string headGroup = chr_na;
        std::string headGroupSynonyms = chr_na;
        std::string level = chr_na;
        std::string totalC = chr_na;
        std::string totalOH = chr_na;
        std::string totalDB = chr_na;
        std::string mass = chr_na;
        std::string formula = chr_na;
        
        lipidMapsCategory = lipidAdduct->get_lipid_string(CATEGORY);
        lipidMapsMainClass = lipidAdduct->get_lipid_string(CLASS);
        species = lipidAdduct->get_lipid_string(SPECIES);
        LipidSpecies* lipid = lipidAdduct->lipid;
        if(lipid) {
            LipidSpeciesInfo info = (*lipid).info;
            nativeLevelName = lipidAdduct->get_lipid_string(info.level);
            lipidMapsMainClass = lipid->get_class_name();
            headGroup = lipid->head_group;
            LipidClassMeta lcMeta = LipidClasses::get_instance().lipid_classes.at(lipid->get_class(headGroup));
            headGroupSynonyms = "[" + join(lcMeta.synonyms, ", ") + "]";
            level = get_lipid_level_str(info.level);
            std::ostringstream cbuffer; 
            cbuffer << info.num_carbon;
            totalC = cbuffer.str();
            std::ostringstream ohbuffer;
            ohbuffer << info.num_hydroxyl;
            totalOH = ohbuffer.str();
            std::ostringstream dbbuffer;
            dbbuffer << info.num_double_bonds;
            totalDB = dbbuffer.str();
            std::ostringstream massbuffer;
            massbuffer << lipidAdduct->get_mass();
            mass = massbuffer.str();
            std::ostringstream formulabuffer;
            formulabuffer << lipidAdduct->get_sum_formula();
            formula = formulabuffer.str();
            
            // Normalized Name	Original Name	Grammar	Lipid Maps Category	Lipid Maps Main Class	Functional Class Abbr	Functional Class Synonyms	Level	Total #C	Total #OH	Total #DB	FA1 Position	FA1 #C	FA1 #OH	FA1 #DB	FA1 Bond Type	FA2 Position	FA2 #C	FA2 #OH	FA2 #DB	FA2 Bond Type	LCB Position	LCB #C	LCB #OH	LCB #DB	LCB Bond Type	FA3 Position	FA3 #C	FA3 #OH	FA3 #DB	FA3 Bond Type	FA4 Position	FA4 #C	FA4 #OH	FA4 #DB	FA4 Bond Type        
            lipidDetails["Normalized Name"] = nativeLevelName;
            lipidDetails["Original Name"] = lipid_name;
            lipidDetails["Grammar"] = grammar;
            lipidDetails["Lipid Maps Category"] = lipidMapsCategory;
            lipidDetails["Lipid Maps Main Class"] = lipidMapsMainClass;
            lipidDetails["Species Name"] = species;
            lipidDetails["Molecular Subspecies Name"] = get_lipid_name_for_level(lipidAdduct, MOLECULAR_SUBSPECIES);
            lipidDetails["Structural Subspecies Name"] = get_lipid_name_for_level(lipidAdduct, STRUCTURAL_SUBSPECIES);
            lipidDetails["Isomeric Subspecies Name"] = get_lipid_name_for_level(lipidAdduct, ISOMERIC_SUBSPECIES);
            lipidDetails["Functional Class Abbr"] = "[" + headGroup + "]";
            lipidDetails["Functional Class Synonyms"] = headGroupSynonyms;
            lipidDetails["Level"] = level;
            lipidDetails["Total C"] = totalC;
            lipidDetails["Total OH"] = totalOH;
            lipidDetails["Total DB"] = totalDB;
            lipidDetails["Mass"] = mass;
            lipidDetails["Sum Formula"] = formula;
            
            int faCnt = 1;
            for(FattyAcid* fap:lipid->get_fa_list()) {
                FattyAcid fa = (*fap);
                std::ostringstream prefs;
                if(fa.lcb) {
                    prefs << "LCB ";
                } else {
                    prefs << "FA" << faCnt << " ";
                }
                string prefix = prefs.str();
                std::ostringstream pos;
                pos << fa.position;
                lipidDetails[prefix + "Position"] = pos.str();
                std::ostringstream nc;
                nc << fa.num_carbon;
                lipidDetails[prefix + "C"] = nc.str();
                std::ostringstream noh;
                noh << fa.num_hydroxyl;
                lipidDetails[prefix + "OH"] = noh.str();
                std::ostringstream ndb;
                ndb << fa.num_double_bonds;
                lipidDetails[prefix + "DB"] = ndb.str();
                string fa_bond_type = "ESTER";
                switch(fa.lipid_FA_bond_type){
                case (UNDEFINED_FA): fa_bond_type = "UNDEFINED"; break;
                case (ESTER): fa_bond_type = "ESTER"; break;
                case (ETHER_PLASMANYL): fa_bond_type = "ETHER_PLASMANYL"; break;
                case (ETHER_PLASMENYL): fa_bond_type = "ETHER_PLASMENYL"; break;
                case (NO_FA): fa_bond_type = "NO_FA"; break;
                default: warning("Unknown bond type in FA " + prefix + " of " + lipid_name); fa_bond_type = "UNKNOWN";
                }
                lipidDetails[prefix + "Bond Type"] = fa_bond_type;
                ++faCnt;
            }
        }
        delete lipidAdduct;
    } else {
        lipidDetails["Normalized Name"] = chr_na;
        lipidDetails["Original Name"] = lipid_name;
        lipidDetails["Grammar"] = "NOT_PARSEABLE";
        warning("Could not find matching parser for grammar '" +grammar+"'! Valid ones are 'Goslin', 'GoslinFragments', 'SwissLipids', 'LipidMaps', 'HMDB'");
    }
    return lipidDetails;
}

/** 
 * Allows to define the specific grammar to use to parse a lipid.
 * It adds information about the parser/grammar that was used to parse the lipid name,
 * if parsing was successful.
 * 
 * The grammar argument can be one of Goslin, GoslinFragments, SwissLipids, LipidMaps, HMDB.
 */
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
SEXP rcpp_parse_lipid_name_with_grammar(std::string lipid_name, std::string target_grammar) {
    CharacterVector lipidDetails;
    LipidParser lipid_parser;
    LipidAdduct* lipidAdduct = NULL;
    for (auto parser : lipid_parser.parser_list) {
        try {
            if(target_grammar.compare(parser->grammar_name)==0) {
                lipidAdduct = parser->parse(lipid_name, false);
                if (lipidAdduct){
                    lipidDetails = handle_lipid(lipidAdduct, lipid_name, parser->grammar_name);
                }
            }
        } catch (LipidException &e){
            forward_exception_to_r(e);
        }
    }
    return lipidDetails;
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
    String chr_na = NA_STRING;
    /* parsing lipid name into a lipid container data structure */
    LipidAdduct* lipidAdduct = NULL;
    try {
        /* create instance of lipid parser containing several grammars */
        LipidParser lipid_parser;
        std::string grammar = chr_na;
        for (auto parser : lipid_parser.parser_list) {
            try {
                LipidAdduct *lipid = parser->parse(lipid_name, false);
                if (lipid){
                    lipidAdduct = lipid;
                    if (parser) {
                        grammar = parser->grammar_name;
                    }
                    break;
                }
            } catch (LipidException &e){
                forward_exception_to_r(e);
            }
        }
        lipidDetails = handle_lipid(lipidAdduct, lipid_name, grammar);
        return lipidDetails;
    } catch(LipidException &e) {
        warning("Parsing of lipid name '" +lipid_name+"' caused an exception: "+ e.what());
        return lipidDetails;
    }
}