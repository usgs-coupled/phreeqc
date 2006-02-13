#if !defined(SSASSEMBLAGESS_H_INCLUDED)
#define SSASSEMBLAGESS_H_INCLUDED

#include "NameDouble.h"
#define EXTERNAL extern
#include "global.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector

#include "char_star.h"

class cxxSSassemblageSS 
{

public:
        cxxSSassemblageSS();
        cxxSSassemblageSS(struct s_s *);
        ~cxxSSassemblageSS();

        enum SS_PARAMETER_TYPE {
                SS_PARM_NONE                          = -1,
                SS_PARM_A0_A1                         = 0,
                SS_PARM_GAMMAS                        = 1,
                SS_PARM_DIST_COEF                     = 2,
                SS_PARM_MISCIBILITY                   = 3,
                SS_PARM_SPINODAL                      = 4,
                SS_PARM_CRITICAL                      = 5,
                SS_PARM_ALYOTROPIC                    = 6,
                SS_PARM_DIM_GUGG                      = 7,
                SS_PARM_WALDBAUM                      = 8,
                SS_PARM_MARGULES                      = 9
        };

        static struct s_s *cxxSSassemblageSS2s_s(std::list<cxxSSassemblageSS>& el);

        void dump_xml(std::ostream& os, unsigned int indent = 0)const;

        void dump_raw(std::ostream& s_oss, unsigned int indent)const;

        void read_raw(CParser& parser);

protected:
        char *name;
        //std::list<cxxSSassemblageSS> ppAssemblageSS;
	cxxNameDouble comps;
        //double total_moles;
        //double dn;
        double a0, a1;
        double ag0, ag1;
        //bool s_s_in;
        bool miscibility;
        //bool spinodal;
        //double tk;
	double xb1, xb2;
        //SS_PARAMETER_TYPE type;
        //double p[4];

public:

};

#endif // !defined(SSASSEMBLAGESS_H_INCLUDED)
