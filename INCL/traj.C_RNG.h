#ifndef TRAJ_H
#define TRAJ_H

#include "sim_types.C_RNG.h"

int Generate_trajectory(int num_traj, float sim_stop,
			float fract, float k_on, float k_off,
			float k_loading, rnap_rate_form *rnap_dwell,
			ribo_rate_form *ribo_dwell, float snap_time,
			int tx_length, float *seq_window,
			float *sample_window, gene_form *gene_set,
			float RNA_lifetime, char *traj_name,
			int start_seed, int num_genes, int traj_index,
			float time_step, float k_continue_degrade,
                        float p_rnap_rnap_push, float p_ribo_rnap_push,
			float p_ribo_ribo_push, int pol_size,
			int ribo_size, int degrade_width,
			float mean_ribo_rate, float min_dwell,
			float generation_time, float b_time, float c_time,
                        float d_time, float k_unloading, float k_to_open,
			float k_CSAT, float mean_RNAP_conc,
			float mean_ribo_conc, int *coord, int stable_RNA,
                        int read_state, char state_file[][200],
			float birth_vol, int tsl_yn, int *gene_id,
			double mu_rnap_p2, double mu_rnap_p3,
			double mu_rnap_e1, double mu_rnap_e2,
			double mu_rnap_e3, double mu_rnap_N, double mu_rnap_c,
                        double mu_rnap_d, int min_effector,
			int rnap_antipause, int ribo_antipause,
                        int ribo_rnap_range, int max_backtrack,
			float t_promo_shutoff, char checklist[][200]);

#endif

