#ifndef FATTY_ACID_H
#define FATTY_ACID_H

#include <string>
#include "cppgoslin/domain/LipidExceptions.h"
#include "cppgoslin/domain/LipidEnums.h"
#include <sstream>

using namespace std;

class FattyAcid {
public:
    string name;
    int position;
    int num_carbon;
    int num_hydroxyl;
    int num_double_bonds;
    LipidFaBondType lipid_FA_bond_type;
    bool lcb;

    FattyAcid(string _name, int _num_carbon, int _num_double_bonds, int _num_hydroxyl, LipidFaBondType _lipid_FA_bond_type, bool _lcb, int _position);
    FattyAcid(FattyAcid* fa);
    virtual ~FattyAcid();
    virtual string to_string(bool special_case, LipidLevel level = NO_LEVEL);
    static string suffix(LipidFaBondType _lipid_FA_bond_type);
};
#endif /* FATTY_ACID_H */
