#undef EXTERN
#undef SET_IT
#undef SET_IT_2
#undef SET_IT_3

#ifdef DEFINE_REPL_VARIABLES
#define EXTERN
#define SET_IT = {0}
#define SET_IT_2 = 0
#define SET_IT_3 = 2
#else
#define EXTERN extern
#define SET_IT
#define SET_IT_2
#define SET_IT_3
#endif

#include "sim_types.h"
#include "gen_def.h"

EXTERN promo_form promo_adds[20];
EXTERN int promoter[MAX_PROMO][PROMO_INFO] SET_IT;				// Information about each template:
										// First index is promoter ID: assuming
#ifndef PROMO_INDEXING								// slow growth, there will be never be
#define PROMO_INDEXING								// more than two copies of the tx unit
										// Second index:
#define ON_OFF 1                                                                // --[1] Current on/off (1/0) status
#define LOC_ON_OFF 2                                                            // --[2] Location in action list of currently-
#define PARENT_PROMO 3                                                          //       available on or off action
#define NEWEST_RNAP 4                                                           // --[3] Location in action list of possible
#define OLDEST_RNAP 5                                                           //       loading event
#define LAST_BP 6                                                               // --[4] Identity of youngest RNAP on template

#endif

EXTERN float mu_promo_on, mu_promo_off;                                         // Const parameters assoc. with promo on/off

EXTERN int num_promoters;                                                       // Tracking replication
EXTERN int num_promo_adds;
EXTERN int birth_promo;
EXTERN int occupied_promo_slot[MAX_PROMO];

EXTERN float repl_time;

EXTERN int pconverter[33];                                                      // Used for standardizing promoter
EXTERN int branch[3];                                                           // indexing

EXTERN int rnap_bound SET_IT_2;
EXTERN int ribo_bound SET_IT_2;	
EXTERN int next_promo SET_IT_3;        	                                       


#ifndef REPL_FNS
#define REPL_FNS

void Start_from_zero(float fract);

#endif




