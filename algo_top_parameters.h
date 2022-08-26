/**
 * algo_top_parameters.h
 *
 * Defines the algorithm related configuration parameters.
 */

#ifndef __ALGO_TOP_PARAMETERS_H__
#define __ALGO_TOP_PARAMETERS_H__

/** Common algorithm definitions, do not remove **/
// Un-comment to have algo as a passthrough for testing purposes
//#define ALGO_PASSTHROUGH

// Number of data words per processing cycle/frame
const int N_WORDS_PER_FRAME	= 9;
// Word has 64 bits
const int N_SHORTS_PER_WORD = 4;

/** Algorithm specific parameters **/
#define TOWERS_IN_ETA 3
#define TOWERS_IN_PHI 4

/** More common algorithm definitions, do not remove **/
// N_INPUT_LINKS and N_OUTPUT_LINKS need to be numerically specified so that wrapper generators work.
// -- N_INPUT_LINKS = TOWERS_IN_ETA*TOWERS_IN_PHI
//#define N_INPUT_LINKS	26
// 12x3 ECAL + 1 HCAL
#define N_INPUT_LINKS	32
// -- N_OUTPUT_LINKS = TOWERS_IN_PHI*2 //For now we have each eta slice in two output links
#define N_OUTPUT_LINKS  2	

// 16G inputs
#define N_INPUT_WORDS_PER_FRAME 9
// 25G outputs
#define N_OUTPUT_WORDS_PER_FRAME 9


#endif
 /* !__ALGO_TOP_PARAMETERS_H__ */


