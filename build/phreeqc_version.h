#define PROGRAM PHREEQC
#define PHREEQC_VER_MAJOR      2
#define PHREEQC_VER_MINOR     11
#define PHREEQC_VER_PATCH      0  
#define PHREEQC_VER_REVISION   0

#define APR_STRINGIFY(n) APR_STRINGIFY_HELPER(n)
#define APR_STRINGIFY_HELPER(n) #n

/** Version number */
#define PHREEQC_VER_NUM        APR_STRINGIFY(PHREEQC_VER_MAJOR) \
                           "." APR_STRINGIFY(PHREEQC_VER_MINOR) \
                           "." APR_STRINGIFY(PHREEQC_VER_PATCH) \
                           "." APR_STRINGIFY(PHREEQC_VER_REVISION)



#define PRODUCT_NAME       APR_STRINGIFY(PROGRAM) \
                       "-" APR_STRINGIFY(PHREEQC_VER_MAJOR) \
                       "." APR_STRINGIFY(PHREEQC_VER_MINOR)
