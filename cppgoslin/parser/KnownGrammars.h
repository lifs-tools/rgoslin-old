#ifndef KNOWN_GRAMMARS_H
#define KNOWN_GRAMMARS_H

#include <string>
#include "cppgoslin/parser/Parser.h"

using namespace std;

static const string goslin_grammar = " \
/* \
 * MIT License \
 *  \
 * Copyright (c) 2017 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de \
 *                    Bing Peng   -   bing.peng {at} isas.de \
 *                    Nils Hoffmann  -  nils.hoffmann {at} isas.de \
 * \
 * Permission is hereby granted, free of charge, to any person obtaining a copy \
 * of this software and associated documentation files (the 'Software'), to deal \
 * in the Software without restriction, including without limitation the rights \
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell \
 * copies of the Software, and to permit persons to whom the Software is \
 * furnished to do so, subject to the following conditions:; \
 *  \
 * The above copyright notice and this permission notice shall be included in all \
 * copies or substantial portions of the Software. \
 *  \
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR \
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, \
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE \
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER \
 * LIABILITY, WHether IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, \
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE \
 * SOFTWARE. \
*/ \
 \
/* This is a BNF / ANTLR4 grammar for lipid subspecies identifiers following \
 * J.K. Pauling et al. 2017, PLoS One, 12(11):e0188394. \
 */ \
 \
grammar Goslin; \
 \
/* first rule is always start rule */ \
lipid : lipid_eof EOF; \
lipid_eof : lipid_pure | lipid_pure adduct_info; \
lipid_pure : gl | pl | sl | cholesterol | mediatorc; \
 \
 \
/* adduct information */ \
adduct_info : '[M' adduct ']' charge charge_sign | adduct_separator '[M' adduct ']' charge charge_sign; \
adduct : '+H' | '+2H' | '+NH4' | '-H' | '-2H' | '+HCOO' | '+CH3COO' | charge_sign arbitrary_adduct; \
arbitrary_adduct : adduct4 | adduct4 adduct4; \
adduct4 : adduct2 | adduct2 adduct2; \
adduct2 : character | character character; \
 \
fa2 : fa2_unsorted | fa2_sorted; \
fa2_unsorted: fa DASH fa | fa UNDERSCORE fa; \
fa2_sorted: fa SLASH fa | fa BACKSLASH fa; \
fa3 : fa3_unsorted | fa3_sorted; \
fa3_unsorted: fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa; \
fa3_sorted: fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa; \
fa4 : fa4_unsorted | fa4_sorted; \
fa4_unsorted: fa DASH fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa UNDERSCORE fa; \
fa4_sorted: fa SLASH fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa BACKSLASH fa; \
 \
/* glycerolipid rules */ \
gl : mgl | dgl | sgl | tgl; \
mgl : hg_mgl_full headgroup_separator fa; \
dgl : hg_dgl_full headgroup_separator gl_species | hg_dgl_full headgroup_separator dgl_subspecies; \
sgl : hg_sgl_full headgroup_separator gl_species | hg_sgl_full headgroup_separator dgl_subspecies; \
tgl : hg_tgl_full headgroup_separator gl_species | hg_tgl_full headgroup_separator tgl_subspecies; \
gl_species : fa; \
dgl_subspecies : fa2; \
tgl_subspecies : fa3; \
 \
hg_mgl_full : hg_mgl | hg_mgl heavy_hg; \
hg_dgl_full : hg_dgl | hg_dgl heavy_hg; \
hg_sgl_full : hg_sgl | hg_sgl heavy_hg; \
hg_tgl_full : hg_tgl | hg_tgl heavy_hg; \
 \
hg_mgl : 'MAG'; \
hg_dgl : 'DAG'; \
hg_sgl : 'MGDG' | 'DGDG' | 'SQDG'; \
hg_tgl : 'TAG'; \
 \
 \
 \
/* phospholipid rules */ \
pl : lpl | dpl | pl_o | cl | mlcl; \
pl_o : lpl_o | dpl_o; \
lpl : hg_lplc headgroup_separator fa; \
lpl_o : hg_lpl_oc plasmalogen_separator fa; \
dpl : hg_plc headgroup_separator pl_species | hg_plc headgroup_separator pl_subspecies; \
dpl_o : hg_pl_oc plasmalogen_separator pl_species | hg_pl_oc plasmalogen_separator pl_subspecies; \
mlcl : hg_mlclc headgroup_separator pl_species | hg_mlclc headgroup_separator mlcl_subspecies; \
cl : hg_clc headgroup_separator pl_species | hg_clc headgroup_separator cl_subspecies; \
 \
pl_species : fa; \
pl_subspecies : fa2; \
mlcl_subspecies : fa3; \
cl_subspecies : fa4; \
 \
heavy_hg : heavy; \
 \
hg_clc : hg_cl | hg_cl heavy_hg; \
hg_cl : 'CL'; \
hg_mlclc : hg_mlcl | hg_mlcl heavy_hg; \
hg_mlcl : 'MLCL'; \
hg_plc : hg_pl | hg_pl heavy_hg; \
hg_pl : 'BMP' | 'CDPDAG' | 'DMPE' | 'MMPE' | 'PA' | 'PC' | 'PE' | 'PEt' | 'PG' | 'PI' | 'PIP' | 'PIP2' | 'PIP3' | 'PS'; \
hg_lplc : hg_lpl | hg_lpl heavy_hg; \
hg_lpl : 'LPA' | 'LPC' | 'LPE' | 'LPG' | 'LPI' | 'LPS'; \
hg_lpl_oc : hg_lpl_o ' O' | hg_lpl_o heavy_hg ' O'; \
hg_lpl_o : 'LPC' | 'LPE'; \
hg_pl_oc : hg_pl_o ' O' | hg_pl_o heavy_hg ' O'; \
hg_pl_o : 'PC' | 'PE'; \
 \
 \
 \
/* sphingolipid rules */ \
sl : lsl | dsl; \
lsl : hg_lslc headgroup_separator lcb; \
dsl : hg_dslc headgroup_separator sl_species | hg_dslc headgroup_separator sl_subspecies; \
 \
sl_species : lcb; \
sl_subspecies : lcb sorted_fa_separator fa; \
 \
hg_lslc : hg_lsl | hg_lsl heavy_hg; \
hg_lsl : 'LCB' | 'LCBP' | 'LHexCer' | 'LSM'; \
hg_dslc : hg_dsl | hg_dsl heavy_hg; \
hg_dsl : 'Cer' | 'CerP' | 'EPC' | 'GB3' | 'GB4' | 'GD3' | 'GM3' | 'GM4' | 'Hex3Cer' | 'Hex2Cer' | 'HexCer' | 'IPC' | 'M(IP)2C' | 'MIPC' | 'SHexCer' | 'SM'; \
 \
 \
 \
/* cholesterol lipids */ \
cholesterol : chc | che; \
chc : ch | ch heavy_hg; \
ch : 'Ch' | 'Cholesterol'; \
che : hg_chec headgroup_separator fa; \
hg_chec : hg_che | hg_che heavy_hg; \
hg_che : 'ChE' | 'CE'; \
 \
 \
/* mediator lipids */ \
mediatorc : mediator | mediator heavy_hg; \
mediator : '10-HDoHE' | '11-HDoHE' | '11-HETE' | '11,12-DHET' | '11(12)-EET'| '12-HEPE' | '12-HETE' | '12-HHTrE' | '12-OxoETE' | '12(13)-EpOME' | '13-HODE' | '13-HOTrE' | '14,15-DHET' | '14(15)-EET' | '14(15)-EpETE' | '15-HEPE' | '15-HETE' | '15d-PGJ2' | '16-HDoHE' | '16-HETE' | '18-HEPE' | '5-HEPE' | '5-HETE' | '5-HpETE' | '5-OxoETE' | '5,12-DiHETE' | '5,6-DiHETE' | '5,6,15-LXA4' | '5(6)-EET' | '8-HDoHE' | '8-HETE' | '8,9-DHET' | '8(9)-EET' | '9-HEPE' | '9-HETE' | '9-HODE' | '9-HOTrE' | '9(10)-EpOME' | 'AA' | 'alpha-LA' | 'DHA' | 'EPA' | 'Linoleic acid' | 'LTB4' | 'LTC4' | 'LTD4' | 'Maresin 1' | 'Palmitic acid' | 'PGB2' | 'PGD2' | 'PGE2' | 'PGF2alpha' | 'PGI2' | 'Resolvin D1' | 'Resolvin D2' | 'Resolvin D3' | 'Resolvin D5' | 'tetranor-12-HETE' | 'TXB1' | 'TXB2' | 'TXB3'; \
 \
 \
 \
/* generic rules */ \
fa : fa_pure | fa_pure heavy_fa | fa_pure ether | fa_pure ether heavy_fa; \
heavy_fa : heavy; \
fa_pure : carbon carbon_db_separator db | carbon carbon_db_separator db db_hydroxyl_separator hydroxyl; \
ether : 'a' | 'p'; \
lcb : lcb_pure | lcb_pure heavy_lcb; \
heavy_lcb : heavy; \
lcb_pure : carbon carbon_db_separator db db_hydroxyl_separator hydroxyl | old_hydroxyl carbon carbon_db_separator db; \
carbon : number; \
db : db_count | db_count round_open_bracket db_position round_close_bracket; \
db_count : number; \
db_position : number cistrans | number cistrans | db_position db_position_separator db_position; \
cistrans : 'E' | 'Z'; \
hydroxyl : number; \
old_hydroxyl : 'd' | 't'; \
number :  digit; \
digit : '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | digit digit; \
 \
 \
 \
 \
 \
heavy : '(+' isotopes  ')'; \
isotopes : isotopes isotopes | isotope; \
isotope : '[' isotope_number ']' isotope_element isotope_count | '[' isotope_number ']' isotope_element; \
isotope_number : number; \
isotope_element : element; \
isotope_count : number; \
element : 'C' | 'H' | 'O' | 'N' | 'P' | 'S'; \
 \
/* separators */ \
SPACE : ' '; \
COLON : ':'; \
SEMICOLON : ';'; \
DASH : '-'; \
UNDERSCORE : '_'; \
SLASH : '/'; \
BACKSLASH : '\\\\'; \
COMMA: ','; \
ROB: '('; \
RCB: ')'; \
FRAGMENT_SEPARATOR : ' - '; \
 \
sorted_fa_separator : SLASH | BACKSLASH; \
adduct_separator : SPACE; \
unsorted_fa_separator : DASH | UNDERSCORE; \
plasmalogen_separator : headgroup_separator | DASH; \
headgroup_separator : SPACE; \
carbon_db_separator : COLON; \
db_hydroxyl_separator : SEMICOLON; \
db_position_separator : COMMA; \
round_open_bracket : ROB; \
round_close_bracket : RCB; \
 \
character : 'a' | 'b' | 'c' | 'd' | 'e' | 'f' | 'g' | 'h' | 'i' | 'j' | 'k' | 'l' | 'm' | 'n' | 'o' | 'p' | 'q' | 'r' | 's' | 't' | 'u' | 'v' | 'w' | 'x' | 'y' | 'z' | 'A' | 'B' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'J' | 'K' | 'L' | 'M' | 'N' | 'O' | 'P' | 'Q' | 'R' | 'S' | 'T' | 'U' | 'V' | 'W' | 'X' | 'Y' | 'Z' |'0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'; \
charge : '1' | '2' | '3' | '4'; \
charge_sign : '-' | '+'; \
\
";

static const string goslin_fragment_grammar = "/* \
 * MIT License \
 *  \
 * Copyright (c) 2017 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de \
 *                    Bing Peng   -   bing.peng {at} isas.de \
 *                    Nils Hoffmann  -  nils.hoffmann {at} isas.de \
 * \
 * Permission is hereby granted, free of charge, to any person obtaining a copy \
 * of this software and associated documentation files (the 'Software'), to deal \
 * in the Software without restriction, including without limitation the rights \
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell \
 * copies of the Software, and to permit persons to whom the Software is \
 * furnished to do so, subject to the following conditions:; \
 *  \
 * The above copyright notice and this permission notice shall be included in all \
 * copies or substantial portions of the Software. \
 *  \
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR \
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, \
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE \
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER \
 * LIABILITY, WHether IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, \
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE \
 * SOFTWARE. \
*/ \
 \
/* This is a BNF / ANTLR4 grammar for lipid subspecies identifiers following \
 * J.K. Pauling et al. 2017, PLoS One, 12(11):e0188394. \
 */ \
 \
grammar GoslinFragments; \
 \
/* first rule is always start rule */ \
lipid : lipid_eof EOF; \
lipid_eof : just_lipid | just_lipid FRAGMENT_SEPARATOR fragment_name; \
just_lipid : lipid_pure | lipid_pure adduct_info; \
lipid_pure : gl | pl | sl | cholesterol | mediatorc; \
 \
/* rules for fragments */ \
fragment_name : frag_char; \
frag_char : frag_char frag_char | character | sign; \
 \
/* adduct information */ \
adduct_info : '[M' adduct ']' charge charge_sign | adduct_separator '[M' adduct ']' charge charge_sign; \
adduct : '+H' | '+2H' | '+NH4' | '-H' | '-2H' | '+HCOO' | '+CH3COO' | charge_sign arbitrary_adduct; \
arbitrary_adduct : adduct4 | adduct4 adduct4; \
adduct4 : adduct2 | adduct2 adduct2; \
adduct2 : character | character character; \
 \
fa2 : fa2_unsorted | fa2_sorted; \
fa2_unsorted: fa DASH fa | fa UNDERSCORE fa; \
fa2_sorted: fa SLASH fa | fa BACKSLASH fa; \
fa3 : fa3_unsorted | fa3_sorted; \
fa3_unsorted: fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa; \
fa3_sorted: fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa; \
fa4 : fa4_unsorted | fa4_sorted; \
fa4_unsorted: fa DASH fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa UNDERSCORE fa; \
fa4_sorted: fa SLASH fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa BACKSLASH fa; \
 \
/* glycerolipid rules */ \
gl : mgl | dgl | sgl | tgl; \
mgl : hg_mgl_full headgroup_separator fa; \
dgl : hg_dgl_full headgroup_separator gl_species | hg_dgl_full headgroup_separator dgl_subspecies; \
sgl : hg_sgl_full headgroup_separator gl_species | hg_sgl_full headgroup_separator dgl_subspecies; \
tgl : hg_tgl_full headgroup_separator gl_species | hg_tgl_full headgroup_separator tgl_subspecies; \
gl_species : fa; \
dgl_subspecies : fa2; \
tgl_subspecies : fa3; \
 \
hg_mgl_full : hg_mgl | hg_mgl heavy_hg; \
hg_dgl_full : hg_dgl | hg_dgl heavy_hg; \
hg_sgl_full : hg_sgl | hg_sgl heavy_hg; \
hg_tgl_full : hg_tgl | hg_tgl heavy_hg; \
 \
hg_mgl : 'MAG'; \
hg_dgl : 'DAG'; \
hg_sgl : 'MGDG' | 'DGDG' | 'SQDG'; \
hg_tgl : 'TAG'; \
 \
 \
 \
/* phospholipid rules */ \
pl : lpl | dpl | pl_o | cl | mlcl; \
pl_o : lpl_o | dpl_o; \
lpl : hg_lplc headgroup_separator fa; \
lpl_o : hg_lpl_oc plasmalogen_separator fa; \
dpl : hg_plc headgroup_separator pl_species | hg_plc headgroup_separator pl_subspecies; \
dpl_o : hg_pl_oc plasmalogen_separator pl_species | hg_pl_oc plasmalogen_separator pl_subspecies; \
mlcl : hg_mlclc headgroup_separator pl_species | hg_mlclc headgroup_separator mlcl_subspecies; \
cl : hg_clc headgroup_separator pl_species | hg_clc headgroup_separator cl_subspecies; \
 \
pl_species : fa; \
pl_subspecies : fa2; \
mlcl_subspecies : fa3; \
cl_subspecies : fa4; \
 \
heavy_hg : heavy; \
 \
hg_clc : hg_cl | hg_cl heavy_hg; \
hg_cl : 'CL'; \
hg_mlclc : hg_mlcl | hg_mlcl heavy_hg; \
hg_mlcl : 'MLCL'; \
hg_plc : hg_pl | hg_pl heavy_hg; \
hg_pl : 'BMP' | 'CDPDAG' | 'DMPE' | 'MMPE' | 'PA' | 'PC' | 'PE' | 'PEt' | 'PG' | 'PI' | 'PIP' | 'PIP2' | 'PIP3' | 'PS'; \
hg_lplc : hg_lpl | hg_lpl heavy_hg; \
hg_lpl : 'LPA' | 'LPC' | 'LPE' | 'LPG' | 'LPI' | 'LPS'; \
hg_lpl_oc : hg_lpl_o ' O' | hg_lpl_o heavy_hg ' O'; \
hg_lpl_o : 'LPC' | 'LPE'; \
hg_pl_oc : hg_pl_o ' O' | hg_pl_o heavy_hg ' O'; \
hg_pl_o : 'PC' | 'PE'; \
 \
 \
 \
/* sphingolipid rules */ \
sl : lsl | dsl; \
lsl : hg_lslc headgroup_separator lcb; \
dsl : hg_dslc headgroup_separator sl_species | hg_dslc headgroup_separator sl_subspecies; \
 \
sl_species : lcb; \
sl_subspecies : lcb sorted_fa_separator fa; \
 \
hg_lslc : hg_lsl | hg_lsl heavy_hg; \
hg_lsl : 'LCB' | 'LCBP' | 'LHexCer' | 'LSM'; \
hg_dslc : hg_dsl | hg_dsl heavy_hg; \
hg_dsl : 'Cer' | 'CerP' | 'EPC' | 'GB3' | 'GB4' | 'GD3' | 'GM3' | 'GM4' | 'Hex3Cer' | 'Hex2Cer' | 'HexCer' | 'IPC' | 'M(IP)2C' | 'MIPC' | 'SHexCer' | 'SM'; \
 \
 \
 \
/* cholesterol lipids */ \
cholesterol : chc | che; \
chc : ch | ch heavy_hg; \
ch : 'Ch' | 'Cholesterol'; \
che : hg_chec headgroup_separator fa; \
hg_chec : hg_che | hg_che heavy_hg; \
hg_che : 'ChE' | 'CE'; \
 \
 \
/* mediator lipids */ \
mediatorc : mediator | mediator heavy_hg; \
mediator : '10-HDoHE' | '11-HDoHE' | '11-HETE' | '11,12-DHET' | '11(12)-EET'| '12-HEPE' | '12-HETE' | '12-HHTrE' | '12-OxoETE' | '12(13)-EpOME' | '13-HODE' | '13-HOTrE' | '14,15-DHET' | '14(15)-EET' | '14(15)-EpETE' | '15-HEPE' | '15-HETE' | '15d-PGJ2' | '16-HDoHE' | '16-HETE' | '18-HEPE' | '5-HEPE' | '5-HETE' | '5-HpETE' | '5-OxoETE' | '5,12-DiHETE' | '5,6-DiHETE' | '5,6,15-LXA4' | '5(6)-EET' | '8-HDoHE' | '8-HETE' | '8,9-DHET' | '8(9)-EET' | '9-HEPE' | '9-HETE' | '9-HODE' | '9-HOTrE' | '9(10)-EpOME' | 'AA' | 'alpha-LA' | 'DHA' | 'EPA' | 'Linoleic acid' | 'LTB4' | 'LTC4' | 'LTD4' | 'Maresin 1' | 'Palmitic acid' | 'PGB2' | 'PGD2' | 'PGE2' | 'PGF2alpha' | 'PGI2' | 'Resolvin D1' | 'Resolvin D2' | 'Resolvin D3' | 'Resolvin D5' | 'tetranor-12-HETE' | 'TXB1' | 'TXB2' | 'TXB3'; \
 \
 \
 \
/* generic rules */ \
fa : fa_pure | fa_pure heavy_fa | fa_pure ether | fa_pure ether heavy_fa; \
heavy_fa : heavy; \
fa_pure : carbon carbon_db_separator db | carbon carbon_db_separator db db_hydroxyl_separator hydroxyl; \
ether : 'a' | 'p'; \
lcb : lcb_pure | lcb_pure heavy_lcb; \
heavy_lcb : heavy; \
lcb_pure : carbon carbon_db_separator db db_hydroxyl_separator hydroxyl | old_hydroxyl carbon carbon_db_separator db; \
carbon : number; \
db : db_count | db_count round_open_bracket db_position round_close_bracket; \
db_count : number; \
db_position : number cistrans | number cistrans | db_position db_position_separator db_position; \
cistrans : 'E' | 'Z'; \
hydroxyl : number; \
old_hydroxyl : 'd' | 't'; \
number :  digit; \
digit : '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | digit digit; \
 \
 \
 \
 \
 \
heavy : '(+' isotopes  ')'; \
isotopes : isotopes isotopes | isotope; \
isotope : '[' isotope_number ']' isotope_element isotope_count | '[' isotope_number ']' isotope_element; \
isotope_number : number; \
isotope_element : element; \
isotope_count : number; \
element : 'C' | 'H' | 'O' | 'N' | 'P' | 'S'; \
 \
/* separators */ \
SPACE : ' '; \
COLON : ':'; \
SEMICOLON : ';'; \
DASH : '-'; \
UNDERSCORE : '_'; \
SLASH : '/'; \
BACKSLASH : '\\\\'; \
COMMA: ','; \
ROB: '('; \
RCB: ')'; \
FRAGMENT_SEPARATOR : ' - '; \
 \
sorted_fa_separator : SLASH | BACKSLASH; \
adduct_separator : SPACE; \
unsorted_fa_separator : DASH | UNDERSCORE; \
plasmalogen_separator : headgroup_separator | DASH; \
headgroup_separator : SPACE; \
carbon_db_separator : COLON; \
db_hydroxyl_separator : SEMICOLON; \
db_position_separator : COMMA; \
round_open_bracket : ROB; \
round_close_bracket : RCB; \
 \
character : 'a' | 'b' | 'c' | 'd' | 'e' | 'f' | 'g' | 'h' | 'i' | 'j' | 'k' | 'l' | 'm' | 'n' | 'o' | 'p' | 'q' | 'r' | 's' | 't' | 'u' | 'v' | 'w' | 'x' | 'y' | 'z' | 'A' | 'B' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'J' | 'K' | 'L' | 'M' | 'N' | 'O' | 'P' | 'Q' | 'R' | 'S' | 'T' | 'U' | 'V' | 'W' | 'X' | 'Y' | 'Z' |'0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'; \
charge : '1' | '2' | '3' | '4'; \
charge_sign : '-' | '+'; \
sign : '-' | '+' | ' ' | ',' | '(' | ')' | '[' | ']' | ':' | '*' | ';'; \
 \
";

static const string lipid_maps_grammar = "//////////////////////////////////////////////////////////////////////////////// \
// MIT License \
//  \
// Copyright (c) 2017 Dominik Kopczynski   -   dominik.kopczynski {at} isas.de \
//                    Bing Peng   -   bing.peng {at} isas.de \
//                    Nils Hoffmann  -  nils.hoffmann {at} isas.de \
// \
// Permission is hereby granted, free of charge, to any person obtaining a copy \
// of this software and associated documentation files (the 'Software'), to deal \
// in the Software without restriction, including without limitation the rights \
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell \
// copies of the Software, and to permit persons to whom the Software is \
// furnished to do so, subject to the following conditions: \
//  \
// The above copyright notice and this permission notice shall be included in all \
// copies or substantial portions of the Software. \
//  \
// THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR \
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, \
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE \
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER \
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, \
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE \
// SOFTWARE. \
//////////////////////////////////////////////////////////////////////////////// \
 \
//// This is a BNF grammer for lipid subspecies identifiers followed by \
//// J.K. Pauling et al. 2017, PLoS One, 12(11):e0188394.  \
 \
grammar LipidMaps; \
 \
/* first rule is always start rule */ \
lipid: lipid_rule EOF; \
lipid_rule: lipid_mono | lipid_mono isotope; \
lipid_mono: lipid_pure | lipid_pure isoform; \
lipid_pure: pure_fa | gl | pl | sl | cholesterol | mediator; \
isoform: square_open_bracket 'rac' square_close_bracket | square_open_bracket 'iso' number square_close_bracket; \
isotope: SPACE round_open_bracket element number round_close_bracket | DASH round_open_bracket element number round_close_bracket | DASH element number; \
element: 'd'; \
 \
 \
/* pure fatty acid */ \
pure_fa: fa; \
 \
fa2 : fa2_unsorted | fa2_sorted; \
fa2_unsorted: fa DASH fa | fa UNDERSCORE fa; \
fa2_sorted: fa SLASH fa | fa BACKSLASH fa; \
fa3 : fa3_unsorted | fa3_sorted; \
fa3_unsorted: fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa; \
fa3_sorted: fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa; \
fa4 : fa4_unsorted | fa4_sorted; \
fa4_unsorted: fa DASH fa DASH fa DASH fa | fa UNDERSCORE fa UNDERSCORE fa UNDERSCORE fa; \
fa4_sorted: fa SLASH fa SLASH fa SLASH fa | fa BACKSLASH fa BACKSLASH fa BACKSLASH fa; \
lcb_fa_sorted: lcb SLASH fa | lcb BACKSLASH fa; \
 \
/* glycerolipid rules */ \
gl: sgl | tgl; \
sgl: hg_sglc sgl_species | hg_sglc sgl_subspecies; \
sgl_species: round_open_bracket fa round_close_bracket | fa; \
sgl_subspecies: round_open_bracket fa2 round_close_bracket | fa2; \
tgl: hg_glc tgl_species | hg_glc tgl_subspecies; \
tgl_species: round_open_bracket fa round_close_bracket | fa; \
tgl_subspecies: round_open_bracket fa3 round_close_bracket | fa3; \
 \
hg_sglc: hg_sgl | hg_sgl headgroup_separator; \
hg_sgl: 'MGDG' | 'DGDG' | 'SQDG' | 'SQMG' | 'DG'; \
hg_glc: hg_gl | hg_gl headgroup_separator; \
hg_gl: 'MG' | 'DG' | 'TG'; \
 \
 \
 \
 \
/* phospholipid rules */ \
pl: lpl | dpl | cl | fourpl; \
lpl: hg_lplc round_open_bracket fa_lpl round_close_bracket | hg_lplc fa_lpl; \
fa_lpl: fa | fa2; \
dpl: hg_ddpl dpl_species | hg_ddpl dpl_subspecies; \
dpl_species: round_open_bracket fa round_close_bracket | fa; \
dpl_subspecies: round_open_bracket fa2 round_close_bracket | fa2; \
cl: hg_clc cl_species | hg_clc cl_subspecies; \
cl_species: round_open_bracket fa round_close_bracket | fa; \
cl_subspecies: round_open_bracket '1\'-' square_open_bracket fa2 square_close_bracket ',3\'-' square_open_bracket fa2 square_close_bracket round_close_bracket | hg_clc '1\'-' square_open_bracket fa2 square_close_bracket ',3\'-' square_open_bracket fa2 square_close_bracket; \
fourpl: hg_fourplc round_open_bracket fa4 round_close_bracket | hg_fourplc fa4; \
 \
hg_ddpl: hg_dplc pip_position | hg_dplc; \
 \
hg_clc: hg_cl | hg_cl headgroup_separator; \
hg_cl: 'CL'; \
hg_dplc: hg_dpl | hg_dpl headgroup_separator; \
hg_dpl: 'LBPA' | 'CDP-DG' | 'DMPE' | 'MMPE' | 'PA' | 'PC' | 'PE' | 'PEt' | 'PG' | 'PI' | 'PIP' | 'PIP2' | 'PIP3' | 'PS' | 'PIM1' | 'PIM2' | 'PIM3' | 'PIM4' | 'PIM5' | 'PIM6' | 'Glc-DG' | 'PGP' | 'PE-NMe2' | 'AC2SGL' | 'DAT' | 'PE-NMe' | 'PT' | 'Glc-GP' | 'NAPE'; \
hg_lplc: hg_lpl | hg_lpl headgroup_separator; \
hg_lpl: 'LysoPC' | 'LPC' | 'LysoPE' | 'LPE' | 'LPIM1' | 'LPIM2' | 'LPIM3' | 'LPIM4' | 'LPIM5' | 'LPIM6' | 'CPA' | 'LPA'; \
hg_fourplc: hg_fourpl | hg_fourpl headgroup_separator; \
hg_fourpl: 'PAT16' | 'PAT18'; \
pip_position: square_open_bracket pip_pos square_close_bracket; \
pip_pos: pip_pos COMMA pip_pos | number '\''; \
 \
 \
 \
/* sphingolipid rules */ \
sl: lsl | dsl | sphingoxine; \
lsl: hg_lslc round_open_bracket lcb round_close_bracket | hg_lslc lcb; \
dsl: hg_dslc dsl_species | hg_dslc dsl_subspecies; \
dsl_species: round_open_bracket lcb round_close_bracket | lcb; \
dsl_subspecies: round_open_bracket lcb_fa_sorted round_close_bracket | lcb_fa_sorted; \
 \
sphingoxine: sphingoxine_pure | sphingoxine_var; \
sphingoxine_pure: sphingosine_name | sphinganine_name; \
sphingoxine_var: ctype headgroup_separator sphingosine_name | ctype headgroup_separator sphinganine_name; \
sphingosine_name: 'Sphingosine' | 'So' | 'Sphingosine-1-phosphate'; \
sphinganine_name: 'Sphinganine' | 'Sa' | 'Sphinganine-1-phosphate'; \
ctype: 'C' number; \
 \
hg_dslc: hg_dsl | hg_dsl headgroup_separator; \
hg_dsl: 'Cer' | 'CerP' | 'EPC' | 'GB3' | 'GB4' | 'GD3' | 'GM3' | 'GM4' | 'Hex3Cer' | 'Hex2Cer' | 'HexCer' | 'IPC' | 'M(IP)2C' | 'MIPC' | 'SHexCer' | 'SM' | 'PE-Cer' | 'PI-Cer' | 'GlcCer' | 'FMC-5' | 'FMC-6' | 'LacCer' | 'GalCer' | '(3\'-sulfo)Galbeta-Cer'; \
 \
hg_lslc: hg_lsl | hg_lsl headgroup_separator; \
hg_lsl: 'SPH' | 'S1P' | 'HexSph' | 'SPC' | 'SPH-P'; \
 \
 \
 \
/* cholesterol lipids */ \
cholesterol: chc | chec; \
chc: ch | ch headgroup_separator; \
ch: 'Cholesterol'; \
chec: che | che headgroup_separator | che_fa; \
che: fa headgroup_separator hg_che; \
che_fa: hg_che round_open_bracket fa round_close_bracket; \
hg_che: 'Cholesteryl ester' | 'Cholesterol ester' | 'CE'; \
 \
 \
/* mediator lipids */ \
mediator: mediator_var | mediator_const; \
mediator_var: mediator_prefix mediator_name_separator mediator_var_names | '(+/-)-' mediator_prefix mediator_name_separator mediator_var_names; \
mediator_prefix: mediator_numbers | mediator_prefix mediator_prefix; \
mediator_numbers: mediator_numbers mediator_separator mediator_numbers | mediator_number_pure | mediator_number_pure mediator_pos; \
mediator_number_pure: number | round_open_bracket number round_close_bracket; \
mediator_pos: 'R' | 'S'; \
mediator_var_names: mediator_var_name | mediator_oxo '-' mediator_var_name; \
mediator_var_name: 'HHTrE' | 'EpOME' | 'HODE' | 'HOTrE' | 'DHET' | 'EET' | 'EpETE' | 'HEPE' | 'HETE' | 'PGJ2' | 'HDoHE' | 'HpETE' | 'ETE' | 'DiHETE' | 'LXA4'; \
mediator_const: 'Arachidonic acid' | 'Arachidonic Acid' | 'alpha-LA' | 'DHA' | 'EPA' | 'Linoleic acid' | 'LTB4' | 'LTC4' | 'LTD4' | 'Maresin 1' | 'Palmitic acid' | 'PGB2' | 'PGD2' | 'PGE2' | 'PGF2alpha' | 'PGI2' | 'Resolvin D1' | 'Resolvin D2' | 'Resolvin D3' | 'Resolvin D5' | 'tetranor-12-HETE' | 'TXB1' | 'TXB2' | 'TXB3'; \
mediator_oxo: 'Oxo' | 'oxo'; \
 \
 \
 \
/* generic rules */ \
fa: fa_unmod | fa_unmod fa_mod; \
fa_unmod: ether fa_pure | fa_pure; \
fa_mod: round_open_bracket modification round_close_bracket; \
modification: modification ',' modification | number mod_text | mod_text; \
mod_text: mod_text mod_text | square_open_bracket mod_text square_close_bracket | 'OH' | 'Ke' | 'S' | 'OOH' | 'My' | 'Me' | 'R' | 'Br' | 'CHO' | 'COOH' | 'Cp' | '(R)'; \
ether: 'P-' | 'O-'; \
fa_pure: carbon carbon_db_separator db | carbon carbon_db_separator db db_hydroxyl_separator hydroxyl; \
lcb_fa: carbon carbon_db_separator db; \
lcb: hydroxyl_lcb lcb_fa | lcb_fa; \
carbon: number; \
db: db_count | db_count round_open_bracket db_position round_close_bracket; \
db_count: number; \
db_position: number | number cistrans | db_position db_position_separator db_position; \
cistrans: 'E' | 'Z'; \
hydroxyl: number; \
hydroxyl_lcb: 'm' | 'd' | 't'; \
number:  '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9' | number number; \
 \
/* separators */ \
SPACE : ' '; \
COLON : ':'; \
SEMICOLON : ';'; \
DASH : '-'; \
UNDERSCORE : '_'; \
SLASH : '/'; \
BACKSLASH : '\\\\'; \
COMMA: ','; \
ROB: '('; \
RCB: ')'; \
SOB: '['; \
SCB: ']'; \
 \
fa_separator: UNDERSCORE | SLASH | BACKSLASH | DASH; \
headgroup_separator: SPACE; \
carbon_db_separator: COLON; \
db_hydroxyl_separator: SEMICOLON; \
db_position_separator: COMMA; \
mediator_separator: COMMA; \
mediator_name_separator: DASH; \
round_open_bracket: ROB; \
round_close_bracket: RCB; \
square_open_bracket: SOB; \
square_close_bracket: SCB; \
 \
";


#endif /* KNOWN_GRAMMARS_H */
