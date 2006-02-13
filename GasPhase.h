#if !defined(GASPHASE_H_INCLUDED)
#define GASPHASE_H_INCLUDED

#include "NumKeyword.h"
#include "NameDouble.h"
#define EXTERNAL extern
#include "global.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector

#include "char_star.h"

class cxxGasPhase : public cxxNumKeyword
{

public:
        cxxGasPhase();
        cxxGasPhase(struct gas_phase *);
        ~cxxGasPhase();

        enum GP_TYPE {
                GP_PRESSURE     = 0,
                GP_VOLUME       = 1
        };

        struct gas_phase *cxxGasPhase2gas_phase();

        struct gas_comp *cxxGasPhaseComp2gas_comp();

        //void dump_xml(std::ostream& os, unsigned int indent = 0)const;

        void dump_raw(std::ostream& s_oss, unsigned int indent)const;

        void read_raw(CParser& parser);


protected:
        cxxNameDouble gasPhaseComps;
	GP_TYPE type;
        double total_p;
        double volume;
public:
        //static std::map<int, cxxGasPhase>& map;

};

#endif // !defined(GASPHASE_H_INCLUDED)
