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


#include "LipidIsomericSubspecies.h"

LipidIsomericSubspecies::LipidIsomericSubspecies(Headgroup* _headgroup, vector<FattyAcid*> *_fa) : LipidStructuralSubspecies(_headgroup, _fa) {            
    info.level = ISOMERIC_SUBSPECIES;
}

LipidIsomericSubspecies::~LipidIsomericSubspecies(){
    
}



LipidLevel LipidIsomericSubspecies::get_lipid_level(){
    return ISOMERIC_SUBSPECIES;
}




string LipidIsomericSubspecies::get_lipid_string(LipidLevel level){
    switch(level){
        case NO_LEVEL:
        case ISOMERIC_SUBSPECIES:
            return LipidMolecularSubspecies::build_lipid_subspecies_name(ISOMERIC_SUBSPECIES);
            
        case SPECIES:
        case STRUCTURAL_SUBSPECIES:
        case MOLECULAR_SUBSPECIES:
        case CATEGORY:
        case CLASS:
            return LipidStructuralSubspecies::get_lipid_string(level);
    
        default:
            throw IllegalArgumentException("LipidIsomericSubspecies does not know how to create a lipid string for level " + std::to_string(level));
    }
}
