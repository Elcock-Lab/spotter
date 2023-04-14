#undef EXTERN
#undef SET_IT

#ifdef DEFINE_SUPER_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.h"
#include "gen_def.h"

EXTERN int start_barrier;
EXTERN int stop_barrier;
EXTERN float start_sigma;
EXTERN float resist_exp;
EXTERN float min_resist;
EXTERN float max_resist;
EXTERN float implicit_f;
EXTERN float *rot_resist;
EXTERN int rnap_sep[MAX_RNAP][2];
EXTERN float lk[MAX_RNAP][2];
EXTERN float sigma[MAX_RNAP][2];
EXTERN float twist[MAX_RNAP][2];
EXTERN float loading_sigma[MAX_RNAP];
EXTERN float helical_repeat[MAX_RNAP][2];
EXTERN float starting_helical_repeat;
EXTERN float relaxed_lk[15000000];
EXTERN float start_lk_up;
EXTERN float start_lk_down;
EXTERN float promo_sigma[MAX_PROMO];
EXTERN float shear_modulus;
EXTERN float torque_factor;
EXTERN float torque_zero_pause_p;
EXTERN float torque_zero_pause_dur;
EXTERN float rnap_p1_effect;
EXTERN float ribo_p1_effect;
EXTERN float rnap_e1_effect;
EXTERN float ribo_e1_effect;
EXTERN float P23_release_torque;
EXTERN int trailer_list[2];
EXTERN float min_upstr_sigma;

EXTERN int rnap_relax;
EXTERN int ribo_resist_model;
EXTERN float relax_factor;
EXTERN float mean_interribo_dist;
EXTERN float flat_resist_param;
EXTERN double complex_resist[MAX_RNAP];
EXTERN int curr_rnap_rot[MAX_RNAP];

EXTERN int writhe_partition;

EXTERN int supercoiling_off;

EXTERN float open_lk[MAX_PROMO];
EXTERN float init_sigma[MAX_PROMO];
EXTERN float curr_init_lk_up[MAX_PROMO];
EXTERN float curr_init_lk_down[MAX_PROMO];
EXTERN float relaxed_open_lk;

#ifndef SUPER_FNS
#define SUPER_FNS

void Assign_rotation_resistance();

float Calculate_net_torque(int rnap_index,int relax_call);

double Torque_dependent_mu_forward(int rnap_index);

double Torque_dependent_mu_enter_P1(int rnap_index);

double Torque_dependent_mu_exit(int rnap_index);

float Apply_180bp_partition_plot(int rnap_index, int side);

float Apply_435bp_partition_plot(int rnap_index, int side);

float Apply_885bp_partition_plot(int rnap_index, int side);

void Distribute_twist_and_writhe(int rnap_index, int side);

void Process_supercoiling_changes(int rnap_index, int move,
                                  int side_adj);

void Update_RNAP_complex_resistance(int rnap_index);

int Update_relaxation_rates(double time, int rnap_index);

int Relax_RNAP_complex(double time, int rnap_index);

#endif
