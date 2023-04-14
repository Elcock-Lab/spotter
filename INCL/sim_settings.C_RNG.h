#undef EXTERN
#undef SET_IT
#undef SET_IT_2

#ifdef DEFINE_SIM_SET_VARIABLES
#define EXTERN
#define SET_IT = {0}
#define SET_IT_2 = 0
#else
#define EXTERN extern
#define SET_IT
#define SET_IT_2
#endif

#include "sim_types.C_RNG.h"
#include "gen_def.C_RNG.h"

EXTERN int QUIET_MODE;

EXTERN int num_to_recycle;
EXTERN int discard[1000][2];

EXTERN int tx_size;
EXTERN int tsl_start[MAX_GENES], tsl_stop[MAX_GENES];                           // Starts and stops for each gene
EXTERN int pol_width, ribo_width;                                               // For processing clashes/collisions of
EXTERN int pol_ribo_width;                                                      // RNAP/RNAP, ribo/ribo, and RNAP/ribo
EXTERN int degrade_rnap_width, degrade_ribo_width;                              // For determining if degradation can start
EXTERN int gene_num, tx_end;                                                    // Renamed for free sharing within file

EXTERN int plasmid;
EXTERN int plasmid_linearized;
EXTERN int chrom;

EXTERN int dna_length;
EXTERN int pol_topo_width;
EXTERN int pol_gyrase_width;

EXTERN int num_free_rnap;
EXTERN int num_free_ribo;
EXTERN int MAX_PROTEINS;

EXTERN int start_pos, stop_pos, min_dist;                                      
EXTERN int stable;
EXTERN int translated;
EXTERN int num_transcripts SET_IT_2;
EXTERN float mean_proteins_per_cell[MAX_GENES];
EXTERN int num_proteins_made[MAX_GENES] SET_IT;                                 // For logging completed proteins
EXTERN int init_proteins[MAX_GENES] SET_IT;
EXTERN int product_added[MAX_GENES] SET_IT;                                     // Zeroed out for each update interval
EXTERN int split;                                                               // Genomic position dividing stable
										// and mRNA (2 cases)
EXTERN float  temp_off;

EXTERN int failed_loads;
EXTERN int max_RNAP_load_gap;

EXTERN int no_super_pause[20000];
EXTERN float no_super_dur[20000];
EXTERN int say_it;

