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
        LipidSpeciesInfo *info = lipidAdduct->lipid->info;
        try {
            return lipidAdduct->get_lipid_string(level);
        } catch(LipidException &e) {
            if(withWarning) {
                warning("Lipid '"+lipidAdduct->get_lipid_string(info->level)+"' with native level '"+get_lipid_level_str(info->level) + "' can not generate name for more specific level '" +get_lipid_level_str(level)+"'!");
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
    String chr_na = NA_STRING;
    
    DataFrame lipidDetails = DataFrame::create();
    lipidDetails.push_back(chr_na, "Normalized.Name");
    lipidDetails.push_back(chr_na, "Original.Name");
    lipidDetails.push_back(chr_na, "Grammar");
    lipidDetails.push_back(chr_na, "Adduct");
    lipidDetails.push_back(chr_na, "Adduct.Charge");
    lipidDetails.push_back(chr_na, "Lipid.Maps.Category");
    lipidDetails.push_back(chr_na, "Lipid.Maps.Main.Class");
    lipidDetails.push_back(chr_na, "Species.Name");
    lipidDetails.push_back(chr_na, "Molecular.Subspecies.Name");
    lipidDetails.push_back(chr_na, "Structural.Subspecies.Name");
    lipidDetails.push_back(chr_na, "Isomeric.Subspecies.Name");
    lipidDetails.push_back(chr_na, "Functional.Class.Abbr");
    lipidDetails.push_back(chr_na, "Functional.Class.Synonyms");
    lipidDetails.push_back(chr_na, "Level");
    lipidDetails.push_back(chr_na, "Total.C");
    lipidDetails.push_back(chr_na, "Total.OH");
    lipidDetails.push_back(chr_na, "Total.DB");
    lipidDetails.push_back(chr_na, "Mass");
    lipidDetails.push_back(chr_na, "Sum.Formula");
    lipidDetails.push_back(chr_na, "FA1.Position");
    lipidDetails.push_back(chr_na, "FA1.C");
    lipidDetails.push_back(chr_na, "FA1.OH");
    lipidDetails.push_back(chr_na, "FA1.DB");
    lipidDetails.push_back(chr_na, "FA1.Bond.Type");
    lipidDetails.push_back(chr_na, "FA1.DB.Positions");
    lipidDetails.push_back(chr_na, "FA2.Position");
    lipidDetails.push_back(chr_na, "FA2.C");
    lipidDetails.push_back(chr_na, "FA2.OH");
    lipidDetails.push_back(chr_na, "FA2.DB");
    lipidDetails.push_back(chr_na, "FA2.Bond.Type");
    lipidDetails.push_back(chr_na, "FA2.DB.Positions");
    lipidDetails.push_back(chr_na, "LCB.Position");
    lipidDetails.push_back(chr_na, "LCB.C");
    lipidDetails.push_back(chr_na, "LCB.OH");
    lipidDetails.push_back(chr_na, "LCB.DB");
    lipidDetails.push_back(chr_na, "LCB.Bond.Type");
    lipidDetails.push_back(chr_na, "FA3.Position");
    lipidDetails.push_back(chr_na, "FA3.C");
    lipidDetails.push_back(chr_na, "FA3.OH");
    lipidDetails.push_back(chr_na, "FA3.DB");
    lipidDetails.push_back(chr_na, "FA3.Bond.Type");
    lipidDetails.push_back(chr_na, "FA3.DB.Positions");
    lipidDetails.push_back(chr_na, "FA4.Position");
    lipidDetails.push_back(chr_na, "FA4.C");
    lipidDetails.push_back(chr_na, "FA4.OH");
    lipidDetails.push_back(chr_na, "FA4.DB");
    lipidDetails.push_back(chr_na, "FA4.Bond.Type");
    lipidDetails.push_back(chr_na, "FA4.DB.Positions");
    
    if (lipidAdduct){
        // grammar = parser->grammar_name;
        std::string originalName = chr_na;
        std::string nativeLevelName = chr_na;
        std::string adductString = chr_na;
        int adductCharge = 0;
        std::string lipidMapsCategory = chr_na;
        std::string lipidMapsMainClass = chr_na;
        std::string species = chr_na;
        std::string headGroup = chr_na;
        std::string headGroupSynonyms = chr_na;
        std::string level = chr_na;
        int totalC = 0;
        int totalOH = 0;
        int totalDB = 0;
        double mass = 0;
        std::string formula = chr_na;
        
        Adduct* adduct = lipidAdduct->adduct;
        lipidMapsCategory = lipidAdduct->lipid->get_lipid_string(CATEGORY);
        lipidMapsMainClass = lipidAdduct->lipid->get_lipid_string(CLASS);
        species = lipidAdduct->lipid->get_lipid_string(SPECIES);
        LipidSpecies* lipid = lipidAdduct->lipid;
        if(lipid) {
            LipidSpeciesInfo *info = lipid->info;
            nativeLevelName = lipidAdduct->lipid->get_lipid_string(info->level);
            if(adduct) {
                adductString = lipidAdduct->adduct->get_lipid_string();
                adductCharge = lipidAdduct->adduct->get_charge();
            }
            lipidMapsMainClass = lipid->headgroup->get_class_name();
            headGroup = lipid->headgroup->headgroup;
            LipidClassMeta lcMeta = LipidClasses::get_instance().lipid_classes.at(lipid->headgroup->get_class(headGroup));
            headGroupSynonyms = "[" + join(lcMeta.synonyms, ", ") + "]";
            level = get_lipid_level_str(info->level);
            totalC = info->num_carbon;
            totalOH = ((info->functional_groups->find("OH") != info->functional_groups->end()) ? info->functional_groups->at("OH").size() : 0);
            totalDB = info->double_bonds->get_num();
            mass = lipidAdduct->get_mass();
            formula = lipidAdduct->get_sum_formula();
            

            // Normalized Name	Original Name	Grammar	Lipid Maps Category	Lipid Maps Main Class	Functional Class Abbr	Functional Class Synonyms	Level	Total #C	Total #OH	Total #DB	FA1 Position	FA1 #C	FA1 #OH	FA1 #DB	FA1 Bond Type	FA2 Position	FA2 #C	FA2 #OH	FA2 #DB	FA2 Bond Type	LCB Position	LCB #C	LCB #OH	LCB #DB	LCB Bond Type	FA3 Position	FA3 #C	FA3 #OH	FA3 #DB	FA3 Bond Type	FA4 Position	FA4 #C	FA4 #OH	FA4 #DB	FA4 Bond Type        
            lipidDetails["Normalized.Name"] = nativeLevelName;
            lipidDetails["Original.Name"] = lipid_name;
            lipidDetails["Adduct"] = adductString;
            lipidDetails["Adduct.Charge"] = adductCharge;
            lipidDetails["Grammar"] = grammar;
            lipidDetails["Lipid.Maps.Category"] = lipidMapsCategory;
            lipidDetails["Lipid.Maps.Main.Class"] = lipidMapsMainClass;
            lipidDetails["Species.Name"] = species;
            lipidDetails["Molecular.Subspecies.Name"] = get_lipid_name_for_level(lipidAdduct, MOLECULAR_SUBSPECIES);
            lipidDetails["Structural.Subspecies.Name"] = get_lipid_name_for_level(lipidAdduct, STRUCTURAL_SUBSPECIES);
            lipidDetails["Isomeric.Subspecies.Name"] = get_lipid_name_for_level(lipidAdduct, ISOMERIC_SUBSPECIES);
            lipidDetails["Functional.Class.Abbr"] = "[" + headGroup + "]";
            lipidDetails["Functional.Class.Synonyms"] = headGroupSynonyms;
            lipidDetails["Level"] = level;
            lipidDetails["Total.C"] = totalC;
            lipidDetails["Total.OH"] = totalOH;
            lipidDetails["Total.DB"] = totalDB;
            lipidDetails["Mass"] = mass;
            lipidDetails["Sum.Formula"] = formula;
            int faCnt = 1;
            for(FattyAcid* fap : lipid->get_fa_list()) {
                string prefix = (fap->lcb ? "LCB" : "FA") + std::to_string(faCnt) + ".";
                lipidDetails[prefix + "Position"] = fap->position;
                lipidDetails[prefix + "C"] = fap->num_carbon;
                lipidDetails[prefix + "OH"] = ((fap->functional_groups->find("OH") != fap->functional_groups->end()) ? fap->functional_groups->at("OH").size() : 0);
                lipidDetails[prefix + "DB"] = lipidAdduct->get_mass();
                string fa_bond_type = "ESTER";
                switch(fap->lipid_FA_bond_type){
                case (UNDEFINED_FA): fa_bond_type = "UNDEFINED"; break;
                case (ESTER): fa_bond_type = "ESTER"; break;
                case (ETHER_PLASMANYL): fa_bond_type = "ETHER_PLASMANYL"; break;
                case (ETHER_PLASMENYL): fa_bond_type = "ETHER_PLASMENYL"; break;
                case (NO_FA): fa_bond_type = "NO_FA"; break;
                default: warning("Unknown bond type in FA " + prefix + " of " + lipid_name); fa_bond_type = "UNKNOWN";
                }
                lipidDetails[prefix + "Bond.Type"] = fa_bond_type;
                std::ostringstream dbPos;
                dbPos << "[";
                std::vector<std::string> dbPosPairs;
                for (auto kv : fap->double_bonds->double_bond_positions){
                    dbPosPairs.push_back(std::to_string(kv.first) + kv.second);
                }
                dbPos << join(dbPosPairs, ", ");
                dbPos << "]";
                lipidDetails[prefix + "DB.Positions"] = dbPos.str();
                ++faCnt;
            }
        }
        delete lipidAdduct;
    } else {
        lipidDetails["Normalized.Name"] = chr_na;
        lipidDetails["Original.Name"] = lipid_name;
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
                    return handle_lipid(lipidAdduct, lipid_name, parser->grammar_name);
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
    String chr_na = NA_STRING;
    /* parsing lipid name into a lipid container data structure */
    LipidAdduct* lipidAdduct = NULL;
    try {
        /* create instance of lipid parser containing several grammars */
        LipidParser lipid_parser;
        lipidAdduct = lipid_parser.parse(lipid_name);
        return handle_lipid(lipidAdduct, lipid_name, (lipidAdduct != 0) ? String(lipid_parser.lastSuccessfulParser->grammar_name) : chr_na);
    } catch(LipidException &e) {
        warning("Parsing of lipid name '" +lipid_name+"' caused an exception: "+ e.what());
        CharacterVector lipidDetails;
        return lipidDetails;
    }
}

