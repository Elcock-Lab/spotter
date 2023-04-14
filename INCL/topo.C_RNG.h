#undef EXTERN
#undef SET_IT

#ifdef DEFINE_TOPO_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.C_RNG.h"
#include "gen_def.C_RNG.h"


EXTERN int topo[MAX_TOPO][TOPO_INFO] SET_IT;                                    // Information about topoisomerases;
										// will share SOURCE_PROMO,MINUS_ONE,
#ifndef TOPO_OBJ_INDICES							// PLUS_ONE, and POSITION with RNAPs
#define TOPO_OBJ_INDICES

#define CYCLES_LEFT 5                                                           // Number of cycles remaining in burst
#define TOPO_TYPE 6                                                             // As advertised
#define PROMO_BLOCKER 7                                                         // Indicates topo prevents RNAP loading
#define POS_BLOCK 11                                                            // UPSTREAM/MIDDLE/DOWNSTREAM wrt tx unit
#define CLOSEST_UP 12                                                           // Indicates this topo is closest to an RNAP
#define CLOSEST_DOWN 13                                                         // As above, but downstream
#define BURST_STATE 14                                                          // Topoisomerase in burst mode? (ON/OFF)

#endif

EXTERN int *topo_search[MAX_TOPO];
EXTERN int first_free_topo;
EXTERN float rate_topoI;
EXTERN float rate_gyrase;
EXTERN float gyrase_min_sigma;

EXTERN int num_topo[MAX_PROMO];
EXTERN int *topo_finder[MAX_PROMO];
EXTERN int curr_topo_list[MAX_PROMO];
EXTERN int topo_master_list[MAX_PROMO][2][MAX_TOPO];
EXTERN int starter_up[MAX_PROMO],starter_down[MAX_PROMO];

EXTERN int num_upstream_update;
EXTERN int num_downstream_update;
EXTERN int upstream_list[MAX_TOPO];
EXTERN int downstream_list[MAX_TOPO];

EXTERN int promo_topo_blocked[MAX_PROMO];
EXTERN int topoIA_width,gyrase_width;
EXTERN int topo_rnap_width[3];
EXTERN int topo_repressor_width[3];
EXTERN int topo_dist_check[5];
EXTERN float mu_topo[5];
EXTERN float mu_topo_binding[3];
EXTERN float mu_topo_unbind[3];

EXTERN float p_rnap_topo_eject[2];

EXTERN int num_burst_rate_distr;
EXTERN float burst_rate_distr[25];
EXTERN float burst_set[25];
EXTERN float mu_topo_intraburst[MAX_TOPO];
EXTERN float mu_topo_interburst[5];
EXTERN float mu_topoIA_intraburst;
EXTERN float mean_lagtime_topoIA;
EXTERN float mu_topo_lagtime[MAX_TOPO];
EXTERN float mean_burst_size;
EXTERN float mean_topo_run_length[5];
EXTERN double g_alpha, g_beta;

EXTERN int inter_rnap_topo;
EXTERN int guided_topo_binding;
EXTERN int guided_gyrase_binding;
EXTERN float relative_gyrase_ds_affinity;
EXTERN int sep_topoIA_interrnap_rate;
EXTERN float relative_topoIA_ds_binding;
EXTERN float rel_extra_topo;

#ifndef TOPO_FNS
#define TOPO_FNS

void Unbind_topoisomerase(int topo_index);

int Find_closest_topoisomerase(int rnap_index, int side);

int Check_upstream(int left_check, int right_check, int check_site);

int Check_downstream(int left_check, int right_check, int check_site);

int Check_middle(int left_check, int right_check, int check_site);

void Update_topoisomerase_RNAP_connections();

void Unbind_topoisomerase_on_schedule(double time, int topo_index);

void Update_topoisomerase_lk(float new_lk, float new_sigma,
			     int new_len, int left_check,
			     int right_check, int loc_check,
			     int promo_index, int rnap_done);

void Load_burst_distribution();

void Determine_mean_lagtime(int topo_index);

void Switch_burst_state(double time, int topo_index);

int Define_plasmid_wraparound_region(int checksite, int promo_index,
                                     int dir_ext, int *wrap_ends);

void Perform_topoisomerase_cycle(double time, int topo_index);

int Attempt_topoisomerase_binding(double time, int topo_type,
                                  int topo_index, int min_spacing);

#endif
