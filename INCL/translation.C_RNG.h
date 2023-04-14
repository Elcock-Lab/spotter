#undef EXTERN
#undef SET_IT

#ifdef DEFINE_TSL_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.C_RNG.h"
#include "gen_def.C_RNG.h"

EXTERN int **rna;                                                               // Information about extant RNAs
EXTERN int MAX_RNA;                                                             // Size is MAX_RNA*RNA_INFO (see below)
										// First index is an arbitary identifier
#ifndef RNA_INDEXING								
#define RNA_INDEXING								// Second index:

#define SOURCE_PROMO 1                                                          // --[1] Identity of source promoter (tu)
#define SOURCE_RNAP 2                                                           // --[2] Identity of RNAP generating RNA
#define FIVE_END 3                                                              // --[3] Location of the 5' end of the RNA
#define THREE_END 4                                                             // --[4] Location of the 3' end of the RNA
#define LOC_START_DEG 143                                                       // --[5]...[27] On/off (1/0) status of tsl's
#define LOC_CONT_DEG 144                                                        // --[28]...[50] Location of loading action
#define BIRTH_TIME 145                                                          //               for respective tsl
#define MATURE 146                                                              // --[51]...[73] Youngest ribo on resp. tsl
#define TOTAL_RIBO 147                                                          // --[74]...[96] Oldest ribo on resp. tsl
										// --[97]...[119] Blocked/unblocked (1/0)
										// 		  status of respective tsl
#define MASTER_INDEX 149							// --[120]...[142] Number of ribosomes loaded
										// --[143] Location of initiate degrade action
										// --[144] Location of continue degrade action
										// --[145] Birth time of the RNA (approx; sec.)
										// --[146] Mature? (y/n->1/0)
#endif	
										// Information about ribosomes in pool
EXTERN int **ribosome;                                                          // Size is MAX_RIBO*RIBO_INFO (see below)
EXTERN int MAX_RIBO_PER_RNA;                                                    // First index is a ribosome identifier
EXTERN int MAX_RIBO;                                                            

#ifndef RIBO_INDEXING								// Second index:
#define RIBO_INDEXING
										// --[1] Identity of the associated RNA
#define SOURCE_RNA 1                                                            // --[2] Position on RNA template
#define ARRIVAL 5                                                               // --[3] Identity of ribosome -1 upstream
#define ABORT_POSSIBLE 7                                                        // --[4] Identity of ribosome +1 downstream
#define IN_CONTACT 8                                                            // --[5] Arrival time (as step)--collisions
#define SOURCE_GENE 9                                                           // --[6] Location on action list of currently-
										//       available move
										// --[7] CSAT abort available?
										// --[8] Location on action list of abort
										// --[9] Identity of associated gene
										// --[10] Status (place to indicate "paused")
										// --[11] Location on action list of (un)paus
#endif										// Note: if "status" incl., expand RIBO_INFO


EXTERN int TSL_ON_OFF[MAX_GENES];                                               // Indices for the RNA object: these will be
EXTERN int LOC_RIBO_LOAD[MAX_GENES];                                            // determined below and treated as constants
EXTERN int NEWEST_RIBO[MAX_GENES];
EXTERN int OLDEST_RIBO[MAX_GENES];
EXTERN int LOC_RIBO_MOVE[MAX_GENES];
EXTERN int NUM_RIBO[MAX_GENES];
EXTERN int TSL_BLOCK[MAX_GENES];
EXTERN int **rna_strip;                                                         // For storing current positions of ribosomes
										// on each copy of the RNA template

EXTERN float  mu_ribo_load[MAX_GENES];                                          // Ribosome loading (constant, by gene)
EXTERN float  basal_ribo_load[MAX_GENES];                                       // For concentration-dependent adjustment
EXTERN float  *mu_ribo_move;                                                    // Position-specific rates for ribosome
EXTERN float  mu_start_degrade, mu_continue_degrade;                            // Rates for starting/continuing degradation
EXTERN float  mu_ribo_csat;                                                     // For collision-stimulated termination

EXTERN ribo_rate_form *ribo_dwell;                                              // Gives access passim

EXTERN int first_free_rna;                                                      // Location of available objects; global
EXTERN int first_free_ribo;

EXTERN int *rna_block_grid;                                                     // For testing tsl start site blockage

EXTERN int num_final_rna;
EXTERN int *rna_to_log;
EXTERN rna_log_form *final_rna_log;

#ifndef TSL_FNS
#define TSL_FNS

int Determine_RNAP_riboload_effect(double time, int rnap_index,
                                   int orig_three_end, int just_made);

int Check_for_ribosome_stacking(int ribo_index);

int Check_for_downstream_ribosome(int ribo_index);

int Add_collision_stimulated_termination_event(double time,int ribo_index);

int Load_ribosome(double time, int gene_index, int rna_index);

int Advance_ribosome(double time, int ribo_index, int gene_index,
		     int rna_index, int *num_proteins_made,
		     int MAX_PROTEINS, float **proteins_made,
		     int *num_codons, float *mean_elongation_rate);

int Add_amino_acid(double time, int ribo_index);

int Remove_stalled_ribosome(double time, int ribo_index);

int Reschedule_ribosome_translocation(double time, int ribo_index);

int Attempt_ribosome_move(double time, int ribo_index,
			  float p_push_ribo_forward,
                          float  p_ribo_rnap_push, float p_rnap_rnap_push,
                          int *num_proteins_made, int MAX_PROTEINS,
                          float **proteins_made, int *num_codons,
			  float *mean_elongation_rate);

int Start_RNA_degradation(double time, int rna_index);

int Degrade_first_RNA_nt(double time, int rna_index);

int Continue_RNA_degradation(double time, int rna_index, int *rna_log_ctr,
                             int **rna_log);

#endif

