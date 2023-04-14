#undef EXTERN
#undef SET_IT
#undef SET_IT_2
#ifdef DEFINE_REG_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.h"
#include "gen_def.h"

EXTERN float promo_kill;                                                        // For emulating promoter shutoff
										// experiments
EXTERN float warm_up_time;
EXTERN int pseudo_shutoff;
EXTERN float t_pseudo_shut;

EXTERN float lk_barrier[MAX_PROMO][4];
EXTERN float relxd_lk[25][5];
EXTERN float barrier_relaxed_open_lk;
EXTERN float no_barr_relaxed_open_lk;
EXTERN float frac[25][5];
EXTERN int far_to_near_dist;
EXTERN int promo_loc[25], promo_barrier_left[20],promo_barrier_right[20];
EXTERN int end_loc[25];
EXTERN float frac_occupied[20];
EXTERN int shutoff_done;

#ifndef REG_FNS
#define REG_FNS

int Activate_promoter(double time, int *promo_cycle, int promo_index,
                      float promoter_on_log[][MAX_PROMO_CYCLES][2]);

int Deactivate_promoter(double time, int *promo_cycle, int promo_index,
                      	float promoter_on_log[][MAX_PROMO_CYCLES][2],
                        int shutoff);

void Shut_down_all_promoters(double time, int *promo_cycle,
                             float promoter_on_log[][MAX_PROMO_CYCLES][2],
                             int promo_index);

#endif
