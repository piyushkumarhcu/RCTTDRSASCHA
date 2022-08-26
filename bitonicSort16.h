
#ifndef _BITONICSORT16_H_
#define _BITONICSORT16_H_

#include <iostream>
#include "ap_int.h"
#include "algo_top.h"
#include "algo_top_parameters.h"
#define Nbclusters 16

using namespace std;

class GreaterSmaller{
public:
    Cluster greater, smaller;
};

void bitonicSort16(Cluster in[Nbclusters], Cluster out[Nbclusters]);

#endif

