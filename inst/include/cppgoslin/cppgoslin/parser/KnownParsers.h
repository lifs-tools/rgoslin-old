#ifndef KNOWN_PARSERS_H
#define KNOWN_PARSERS_H


#include "cppgoslin/parser/GoslinFragmentParserEventHandler.h"
#include "cppgoslin/parser/GoslinParserEventHandler.h"
#include "cppgoslin/parser/LipidMapsParserEventHandler.h"
#include "cppgoslin/parser/SwissLipidsParserEventHandler.h"
#include "cppgoslin/parser/BaseParserEventHandler.h"

class GoslinParser : public Parser<LipidAdduct*> {
public:
    GoslinParser();
    ~GoslinParser();
};


class GoslinFragmentParser : public Parser<LipidAdduct*> {
public:
    GoslinFragmentParser();
    ~GoslinFragmentParser();
};


class LipidMapsParser : public Parser<LipidAdduct*> {
public:
    LipidMapsParser();
    ~LipidMapsParser();
};


class SwissLipidsParser : public Parser<LipidAdduct*> {
public:
    SwissLipidsParser();
    ~SwissLipidsParser();
};


class LipidParser {
public:
    vector<Parser<LipidAdduct*>*> parser_list;
    
    LipidParser();
    ~LipidParser();
    LipidAdduct* parse(string lipid_name);
};      

#endif /* KNOWN_PARSERS_H */
