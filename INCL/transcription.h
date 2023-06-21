#undef EXTERN
#undef SET_IT

#ifdef DEFINE_TX_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.h"
#include "gen_def.h"

EXTERN int rnap[MAX_RNAP][RNAP_INFO] SET_IT;					// Information about RNA polymerases
										// First index is an identifier
#ifndef RNAP_INDICES_DEFINED
#define RNAP_INDICES_DEFINED
										// Second index:
#define POSITION 2								// --[1] Identity of the associated promoter 
#define MINUS_ONE 3								// --[2] Position on the DNA template
#define PLUS_ONE 4								// --[3] Identity of RNAP -1 upstream
#define STATE 6									// --[4] Identity of RNAP +1 downstream
#define LOC_MOVE 6								// --[5] Status (Onpath pre/post,p1,p2,p3)
#define CLOSEST_RIBO 7
#define LOC_PAUSE 7								// --[6] Location on action list of currently-
#define LOC_ABORT 8								//       available move
#define RIBO_TETHER 8
#define RNAP_RNA 9								// --[7] Location on action list of (un)pause
										// --[8] Location on action list of abort
										// --[9] Associated RNA object
#define BARRIER_AHEAD 10
#define BARRIER_BEHIND 11
#define TOPO_UPSTREAM 12
#define TOPO_DOWNSTREAM 13

#define SUPER_INDEX 14

#endif

EXTERN int *rnap_search[MAX_RNAP];
EXTERN int **dna_strip;								// For storing current positions of RNAP
										// on each copy of the DNA template
EXTERN float  mu_promo_load;                                                    // Const parameter assoc. with RNAP loading
EXTERN float  basal_promo_load;                                                 // For concentration-dep adj of rates
EXTERN float  *mu_rnap_move;                                                    // Position-specific rates for RNAP
EXTERN float  mu_closed_to_open[MAX_PROMO], mu_unbind_rnap;                     // For multistep transcription initiation
EXTERN float orig_mu_closed_to_open;
EXTERN float  min_rnap_dwell;                                                   // A max rate used for an RNAP in contact
										// with an upstream RNAP and so always forward-
										// disposed

EXTERN rnap_rate_form *rnap_dwell;                                              // Gives access passim

EXTERN double mu_pause[5];                                                       // Collected inv rates for P1,P2,P3 entry
EXTERN double mu_exit[5];                                                        // Collected inv rates for P1,P2,P3 exit
EXTERN double mu_rnap_add;                                                       // If distinct ribo pause states added,
										// change to mu_rnap_pause,mu_rnap_exit
EXTERN float orig_p_closed_to_open;
EXTERN float p_closed_to_open[MAX_PROMO];                                       // Pre-calculated rate-based p's
EXTERN float p_p1_to_p2;                                                        // to enable sampling from uniform
EXTERN float p_p2_to_p3;                                                        // rather then exp distribution at
EXTERN float mu_net_promo;                                                      // rxn forks; others (below) pos-specif
EXTERN float mu_net_p1;                                                         // dwelltimes are net for firing of
EXTERN float mu_net_p2;                                                         // either of a pair of chosen rxns
										// at a particular stage

EXTERN int BUBBLE_ADJ;
EXTERN int bubble[2];

EXTERN int DNA_SCRUNCH;
EXTERN int hybrid_length;
EXTERN int TSS_adj;

EXTERN int NO_PAUSE_MODEL;                                                      // If on, only RNAP velocity is affected
										// by supercoiling--no pauses occur;
										// NB: not compatible with backtracking

EXTERN int ANTICASCADE_ON;                                                      // Determines whether blocked RNAPs
										// can enter state P1

EXTERN int RNAP_ANTIPAUSE;                                                      // Determines whether in-contact
EXTERN int RIBO_ANTIPAUSE;                                                      // RNAP(ribosome) affects access to
										// pause state

EXTERN int BACKTRACK_ON;
EXTERN int MAX_BACKTRACK;                                                       // Max distance from start in P2/3
EXTERN int EFFECTOR_STATE;                                                      // Used to determine which abutting
EXTERN float bt_adj;                                                            // RNAPs can cause pause exit in
										// upstream RNAPs

EXTERN int RIBO_RNAP_RANGE;                                                     // Max number of nt's allowed between
										// ribosome and RNAP for interaction

EXTERN int TETHER_LENGTH;                                                       // Maximum allowable RNAP-ribosome
										// distance after formation of
										// expressome

EXTERN int TX_INIT_SC_FX;
EXTERN float k_unbind;
EXTERN float k_open;

EXTERN int first_free_rnap;

EXTERN float mu_term;

EXTERN last_rate_form last_rates[MAX_RNAP];

EXTERN float mu_rnap_surveillance;
EXTERN float mu_rnap_eval;

#ifndef TX_FNS
#define TX_FNS

void Get_position_without_shutoff();

int Find_first_free(int **look);

int Find_first_free_RNAP();

int Check_for_RNAP_RNAP_contacts(int rnap_index, int updown);

int Find_RNAP_nearest_ribosome_distance(int rnap_index);

int Check_for_RNAP_ribosome_contacts(int rnap_index, int range);

int Check_ribosome_tether(int rnap_index, int move);

int Determine_dominant_trailer(int rnap_index, float rnap_cf,
                               float ribo_cf);

int Advance_pause_state(double time, int rnap_index);

int Process_pretranslocated_state(double time, int rnap_index);

int Exit_off_pathway_state(double time, int rnap_index);

int Update_RNAP_rates(double time, int rnap_index);

void Freeze_RNAP(double time, int rnap_index);

void Unfreeze_RNAP(double time, int rnap_index);

int Determine_off_pathway_move(double time,int rnap_index);

int Load_RNAP(double time, int promoter_index);

int Activate_RNAP(double time, int rnap_index);

void Generate_RNAP_reaction_set(int rnap_index,int *rxn_set);

int Process_terminating_RNAP(double time, int rnap_index);

int Terminate_transcription(double time, int rnap_index);

int Schedule_RNAP_checkup(double time, int rnap_index);

int Evaluate_transcriptional_progress(double time,int rnap_index);

int Remove_RNAP(double time, int rnap_index);

int Rectify_downstream_RNAP(double time, int rnap_index);

int Add_RNA_nt(double time, int rnap_index, int MAX_PROTEINS,
               float **proteins_made);

int Advance_RNAP(double time, int rnap_index, int MAX_PROTEINS,
                 float **proteins_made, int displace);

int Process_blocked_RNAP(double time, int rnap_index);

int Attempt_RNAP_move(double time, int rnap_index,
                      float rnap_push_prob, int MAX_PROTEINS,
                      float **proteins_made, int direction);

#endif
