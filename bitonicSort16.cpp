
#include "bitonicSort16.h"

//Main CAE block (compare and exchange)

GreaterSmaller AscendDescend(Cluster &x, Cluster &y){
//#pragma HLS PIPELINE II=9
#pragma HLS INLINE
    GreaterSmaller s;

    s.greater = (x.clusterEnergy() > y.clusterEnergy()) ? x : y;
    s.smaller = (x.clusterEnergy() > y.clusterEnergy()) ? y : x;

    return s;
}

void FourinSmallFir(Cluster &x0, Cluster &x1, Cluster &x2, Cluster &x3,
					Cluster &y0, Cluster &y1, Cluster &y2, Cluster &y3){
//#pragma HLS PIPELINE II=9
#pragma HLS INLINE
    GreaterSmaller res, res1;
        res = AscendDescend(x0, x2);
        y0 = res.smaller; y2 = res.greater;

        res1 = AscendDescend(x1, x3);
        y1 = res1.smaller; y3 = res1.greater;
}

void FourinGreatFir(Cluster &x0, Cluster &x1, Cluster &x2, Cluster &x3,
					Cluster &y0, Cluster &y1, Cluster &y2, Cluster &y3){
//#pragma HLS PIPELINE II=9
#pragma HLS INLINE
    GreaterSmaller res, res1;
        res = AscendDescend(x0, x2);
        y0 = res.greater; y2 = res.smaller;

        res1 = AscendDescend(x1, x3);
        y1 = res1.greater; y3 = res1.smaller;
}

void EightinSmallFir(Cluster &x0, Cluster &x1,Cluster &x2, Cluster &x3, Cluster &x4, Cluster &x5,
                Cluster &x6, Cluster &x7, Cluster &y0, Cluster &y1, Cluster &y2, Cluster &y3, Cluster &y4, Cluster &y5, Cluster &y6, Cluster &y7){
//#pragma HLS PIPELINE II=9
#pragma HLS INLINE
        GreaterSmaller res, res1, res2, res3;
        res = AscendDescend(x0, x4);
        y0 = res.smaller; y4 = res.greater;

        res1 = AscendDescend(x1, x5);
        y1 = res1.smaller; y5 = res1.greater;

        res2 = AscendDescend(x2, x6);
        y2 = res2.smaller; y6 = res2.greater;

        res3 = AscendDescend(x3, x7);
        y3 = res3.smaller; y7 = res3.greater;

}

void EightinGreatFir(Cluster &x0, Cluster &x1, Cluster &x2, Cluster &x3, Cluster &x4, Cluster &x5,
                Cluster &x6, Cluster &x7, Cluster &y0, Cluster &y1, Cluster &y2, Cluster &y3, Cluster &y4, Cluster &y5, Cluster &y6,Cluster &y7){
//#pragma HLS PIPELINE II=9
#pragma HLS INLINE
        GreaterSmaller res, res1, res2, res3;
        res = AscendDescend(x0, x4);
        y0 = res.greater; y4 = res.smaller;

        res1 = AscendDescend(x1, x5);
        y1 = res1.greater; y5 = res1.smaller;

        res2 = AscendDescend(x2, x6);
        y2 = res2.greater; y6 = res2.smaller;

        res3 = AscendDescend(x3, x7);
        y3 = res3.greater; y7= res3.smaller;

}

void bitonicSort16(Cluster in[Nbclusters], Cluster out[Nbclusters]){
//#pragma HLS PIPELINE II=9
#pragma HLS ARRAY_PARTITION variable=in
#pragma HLS ARRAY_PARTITION variable=out

GreaterSmaller result, result1;

Cluster a[Nbclusters], b[Nbclusters], c[Nbclusters], d[Nbclusters], e[Nbclusters], f[Nbclusters], g[Nbclusters], h[Nbclusters], l[Nbclusters];
#pragma HLS ARRAY_PARTITION variable=a
#pragma HLS ARRAY_PARTITION variable=b
#pragma HLS ARRAY_PARTITION variable=c
#pragma HLS ARRAY_PARTITION variable=d
#pragma HLS ARRAY_PARTITION variable=e
#pragma HLS ARRAY_PARTITION variable=f
#pragma HLS ARRAY_PARTITION variable=g
#pragma HLS ARRAY_PARTITION variable=h
#pragma HLS ARRAY_PARTITION variable=l

//Starting of first stage

for(loop i=0; i<Nbclusters/4; i++){
    #pragma HLS unroll
    result = AscendDescend(in[4*i], in[4*i+1]);
    a[4*i] = result.smaller; a[4*i+1] = result.greater; 
}

for(loop i=0; i<Nbclusters/4; i++){
    #pragma HLS unroll
    result1 = AscendDescend(in[4*i+2], in[4*i+3]);
    a[4*i+2] = result1.greater; a[4*i+3] = result1.smaller; 
}

//Starting of second stage

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    FourinSmallFir(a[8*i], a[8*i+1], a[8*i+2], a[8*i+3], b[8*i], b[8*i+1], b[8*i+2], b[8*i+3]);
    FourinGreatFir(a[8*i+4], a[8*i+5], a[8*i+6], a[8*i+7], b[8*i+4], b[8*i+5], b[8*i+6], b[8*i+7]); 
}

//Starting of third stage

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    result = AscendDescend(b[8*i], b[8*i+1]);
    c[8*i] = result.smaller; c[8*i+1] = result.greater; 
 
    result1 = AscendDescend(b[8*i+2], b[8*i+3]);
    c[8*i+2] = result1.smaller; c[8*i+3] = result1.greater;
}

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    result = AscendDescend(b[8*i+4], b[8*i+5]);
    c[8*i+4] = result.greater; c[8*i+5] = result.smaller; 

    result1 = AscendDescend(b[8*i+6], b[8*i+7]);
    c[8*i+6] = result1.greater; c[8*i+7] = result1.smaller;
}

//Starting of fourth stage

EightinSmallFir(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
EightinGreatFir(c[8], c[9], c[10], c[11], c[12], c[13], c[14], c[15], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);

//Starting fifth stage

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    FourinSmallFir(d[4*i], d[4*i+1], d[4*i+2], d[4*i+3], e[4*i], e[4*i+1], e[4*i+2], e[4*i+3]);
}

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    FourinGreatFir(d[4*i+8], d[4*i+9], d[4*i+10], d[4*i+11], e[4*i+8], e[4*i+9], e[4*i+10], e[4*i+11]);
}

//Starting sixth stage

for(loop i=0; i<Nbclusters/4; i++)
{
    #pragma HLS unroll
    result = AscendDescend(e[2*i], e[2*i+1]);
    f[2*i] = result.smaller; f[2*i+1] = result.greater;
}

for(loop i=0; i<Nbclusters/4; i++){
    #pragma HLS unroll
    result1 = AscendDescend(e[2*i+8], e[2*i+9]);
    f[2*i+8] = result1.greater; f[2*i+9] = result1.smaller;
}

//Starting seventh stage

for(loop i=0; i<Nbclusters/2; i++){
    #pragma HLS unroll
    result = AscendDescend(f[i], f[i+8]);
    g[i] = result.smaller; g[i+8] = result.greater;
}

//starting eighth stage

for(loop i=0; i<Nbclusters/8; i++){
    #pragma HLS unroll
    EightinSmallFir(g[8*i], g[8*i+1], g[8*i+2], g[8*i+3], g[8*i+4], g[8*i+5], g[8*i+6], g[8*i+7],
                h[8*i], h[8*i+1], h[8*i+2], h[8*i+3], h[8*i+4], h[8*i+5], h[8*i+6], h[8*i+7]);
}

//starting ninth stage

for(loop i=0; i<Nbclusters/4; i++){
    #pragma HLS unroll
    FourinSmallFir(h[4*i], h[4*i+1], h[4*i+2], h[4*i+3], l[4*i], l[4*i+1], l[4*i+2], l[4*i+3]);
}

//starting tenth stage

for(loop i=0; i<Nbclusters/2; i++){
    #pragma HLS unroll
    result = AscendDescend(l[2*i], l[2*i+1]);
    out[2*i] = result.smaller; out[2*i+1] = result.greater;
}

//out[0] = in[0] ;

}


