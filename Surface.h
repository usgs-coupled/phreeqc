#if !defined(SURFACE_H_INCLUDED)
#define SURFACE_H_INCLUDED

#include "NumKeyword.h"
#define EXTERNAL extern
#include "global.h"
#include <cassert> // assert
#include <map>     // std::map
#include <string>  // std::string
#include <list>    // std::list
#include <vector>  // std::vector

#include "char_star.h"
#include "SurfComp.h"
#include "SurfCharge.h"

class cxxSurface : public cxxNumKeyword
{

public:
        cxxSurface();
        cxxSurface(struct surface *);
        ~cxxSurface();

        struct surface *cxxSurface2surface();

        struct surf_comp *cxxSurfComp2surf_comp();

        void dump_xml(std::ostream& os, unsigned int indent = 0)const;

        void dump_raw(std::ostream& s_oss, unsigned int indent)const;

        void read_raw(CParser& parser);

        bool get_related_phases(void);

        bool get_related_rate(void);


protected:
        std::list<cxxSurfComp> surfComps;
        std::list<cxxSurfCharge> surfCharges;
        bool diffuse_layer;
        bool edl;
        bool only_counter_ions;
        bool donnan;
        double thickness;
        //double debye_units;
        //int transport;

public:
        //static std::map<int, cxxSurface>& map;

};

#endif // !defined(SURFACE_H_INCLUDED)
