#undef EXTERN
#undef SET_IT

#ifdef DEFINE_RXN_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "reaction_macros.C_RNG.h"
#include "sim_types.C_RNG.h"
#include "gen_def.C_RNG.h"

EXTERN reaction_form reactions[MAX_REACTIONS] SET_IT;
EXTERN tree_form reaction_queue[MAX_REACTIONS];

EXTERN int reaction_index[MAX_REACTIONS];                                       // Index for tracking reactions in queue
EXTERN int queue_size;                                                          // Current number of reactions in queue

#ifndef RXN_FNS
#define RXN_FNS

void Swap_queue_positions(int i, int j);

void Update_reaction_position(int position);

void Update_reaction_queue(int reaction, double time);

void Add_reaction_to_queue(int reaction, double time);

void Order_queue(int i);

void Build_initial_reaction_queue(int length);

void Remove_reaction(int reaction);

void Initialize_reaction_times();

void Assign_reaction_information();

int Uniform_int(int max_num);

double Sample_gamma_distribution(double g_alpha, double g_beta);

double Sample_gaussian_distribution (double sigma);

double Sample_uniform_distribution();

double Find_firing_time(double mu);

void Add_kymograph_reactions(double time);

void Add_pseudoreactions_to_queue(
			int num_fish_times, float *fish_times,
			int num_seq_times, float *seq_times,
			float snap_time, float  repl_time);


#endif
