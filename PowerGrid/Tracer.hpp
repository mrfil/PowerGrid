
/*****************************************************************************

    File Name   [Tracer.hpp]

    Synopsis    [Tracer class for inserting nvtx range markers using RAII idiom. 
                Look at NVIDIA dev blog for more details.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/10/2018]

*****************************************************************************/
#ifndef PowerGrid_Tracer_hpp
#define PowerGrid_Tracer_hpp

#ifdef USE_NVTX
class Tracer {
public:
    Tracer(const char* name)
    {
        nvtxRangePushA(name);
    }
    ~Tracer()
    {
        nvtxRangePop();
    }
};
#define RANGE(name) Tracer uniq_name_using_macros(name);
#else
#define RANGE(name)
#endif

#endif
