#ifndef LIPID_ISOMERIC_SUBSPECIES_H
#define LIPID_ISOMERIC_SUBSPECIES_H

#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidStructuralSubspecies.h"
#include "cppgoslin/domain/IsomericFattyAcid.h"
#include <sstream>
#include <vector>

using namespace std;

class LipidIsomericSubspecies : public LipidStructuralSubspecies {
public:
    LipidIsomericSubspecies(string head_group, vector<FattyAcid*>* fa = NULL);
    ~LipidIsomericSubspecies();
    string build_lipid_isomeric_substructure_name();
    string get_lipid_string(LipidLevel level = NO_LEVEL);
    LipidLevel get_lipid_level();
            
};

#endif /* LIPID_ISOMERIC_SUBSPECIES_H */
