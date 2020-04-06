/*
MIT License

Copyright (c) 2020 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de
                   Nils Hoffmann  -  nils.hoffmann {at} isas.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include "LipidSpecies.h"

#include <iostream>
using namespace std;

LipidSpecies::LipidSpecies(string _head_group, LipidCategory _lipid_category, LipidClass _lipid_class, LipidSpeciesInfo *lipid_species_info){
    head_group = _head_group;
    lipid_category = (_lipid_category != NO_CATEGORY) ? _lipid_category : get_category(head_group);
    lipid_class = (_lipid_class != NO_CLASS) ? _lipid_class : get_class(head_group);
    if (lipid_species_info != NULL){
        info.level = lipid_species_info->level;
        info.num_carbon = lipid_species_info->num_carbon;
        info.num_hydroxyl = lipid_species_info->num_hydroxyl;
        info.num_double_bonds = lipid_species_info->num_double_bonds;
        info.lipid_FA_bond_type = lipid_species_info->lipid_FA_bond_type;
    }
    else {
        info.level = UNDEFINED_LEVEL;
    }
    use_head_group = false;
}


LipidSpecies::~LipidSpecies(){
    for (auto fa : fa_list){
        delete fa;
    }
}


string LipidSpecies::get_lipid_string(LipidLevel level){
    switch (level){
            
        default:
            throw RuntimeException("LipidSpecies does not know how to create a lipid string for level " + to_string(level));
            
        case UNDEFINED_LEVEL:
            throw RuntimeException("LipidSpecies does not know how to create a lipid string for level " + to_string(level));
            
        case CATEGORY:
            return get_category_string(lipid_category);
            
        case CLASS:
            return (!use_head_group ? get_class_string(lipid_class) : head_group);
            
        case NO_LEVEL:
        case SPECIES:
            if (!validate()){
                throw ConstraintViolationException("No fatty acly chain information present for lipid " + get_class_string(lipid_class));
            }
            stringstream st;
            st << (!use_head_group ? get_class_string(lipid_class) : head_group);
            
            if (info.num_carbon > 0){
                bool special_case = (lipid_class == PC) | (lipid_class == LPC) | (lipid_class == PE) | (lipid_class == LPE);
                st << ((lipid_category !=  ST) ? " " : "/");
                st << ((FattyAcid)info).to_string(special_case);
            }
            
            return st.str();
    }
    return "";
}

LipidCategory LipidSpecies::get_category(string _head_group){
    if (!StringCategory.size()){
        for (const auto& kvp : lipid_classes){
            LipidCategory category = kvp.second.lipid_category;
            for (auto hg : kvp.second.synonyms){
                StringCategory.insert(pair<string, LipidCategory>(hg, category));
            }
        }
    }
    
    
    auto cat = StringCategory.find(_head_group);
    return (cat != StringCategory.end()) ? StringCategory.at(_head_group) : UNDEFINED;
}


LipidLevel LipidSpecies::get_lipid_level(){
    return SPECIES;
}


LipidClass LipidSpecies::get_class(string _head_group){
    if (!StringClass.size()){
        for (auto kvp : lipid_classes){
            LipidClass l_class = kvp.first;
            for (auto hg : kvp.second.synonyms){
                StringClass.insert({hg, l_class});
            }
        }
    }
    
    auto cl = StringClass.find(_head_group);
    return (cl != StringClass.end()) ? cl->second : UNDEFINED_CLASS;
}


string LipidSpecies::get_class_string(LipidClass _lipid_class){
    if (!ClassString.size()){
        for (auto kvp : lipid_classes){
            ClassString.insert({kvp.first, kvp.second.synonyms.at(0)});
        }
    }
    auto cl = ClassString.find(_lipid_class);
    return (cl != ClassString.end()) ? ClassString.at(_lipid_class) : "UNDEFINED";
}


string LipidSpecies::get_class_name(){
    return lipid_classes.at(lipid_class).class_name;
}


string LipidSpecies::get_category_string(LipidCategory _lipid_category){
    return CategoryString.at(_lipid_category);
}

vector<FattyAcid*> LipidSpecies::get_fa_list(){
    return fa_list;
}


bool LipidSpecies::validate(){
    if (use_head_group) return true;
    if (lipid_classes.find(lipid_class) == lipid_classes.end()) return false;
    return lipid_classes.at(lipid_class).max_num_fa == 0 || (lipid_classes.at(lipid_class).max_num_fa > 0 && info.num_carbon >= 2);
}

