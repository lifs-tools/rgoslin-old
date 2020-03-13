#include "FattyAcid.h"

FattyAcid::FattyAcid(string _name, int _num_carbon, int _num_double_bonds, int _num_hydroxyl, LipidFaBondType _lipid_FA_bond_type, bool _lcb, int _position){
    name = _name;
    position = _position;
    num_carbon = _num_carbon;
    num_double_bonds = _num_double_bonds;
    num_hydroxyl = _num_hydroxyl;
    lipid_FA_bond_type = _lipid_FA_bond_type;
    lcb = _lcb;
    
    if (num_carbon < 2){
        throw ConstraintViolationException("FattyAcid must have at least 2 carbons! Got " + to_string(num_carbon));
    }
    
    if (position < -1){
        throw ConstraintViolationException("FattyAcid position must be greater or equal to -1 (undefined) or greater or equal to 0 (0 = first position)!");
    }
    
    if (num_hydroxyl < 0){
        throw ConstraintViolationException("FattyAcid must have at least 0 hydroxy groups! Got " + to_string(num_hydroxyl));
    }
}

FattyAcid::FattyAcid(FattyAcid* fa){
    name = fa->name;
    position = fa->position;
    num_carbon = fa->num_carbon;
    num_hydroxyl = fa->num_hydroxyl;
    lipid_FA_bond_type = fa->lipid_FA_bond_type;
    lcb = fa->lcb;
}

FattyAcid::~FattyAcid(){

}

string FattyAcid::suffix(LipidFaBondType lipid_FA_bond_type){
    switch(lipid_FA_bond_type){
        case (ETHER_PLASMANYL): return "a";
        case (ETHER_PLASMENYL): return "p";
        default: return "";
    }
}


string FattyAcid::to_string(bool special_case, LipidLevel level){
    stringstream s;
    
    string lipid_suffix = suffix(lipid_FA_bond_type);
    if (special_case && lipid_suffix.length() > 0) s << "O-";
    s << num_carbon << ":" << num_double_bonds;
    
    if (num_hydroxyl) s << ";" << num_hydroxyl;
    s << lipid_suffix;
    return s.str();
}
