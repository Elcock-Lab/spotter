#ifndef IO_H
#define IO_H

#include "sim_types.h"

int Read_basic_parameters(char *name,char name_set[][200],float *time_set,
                          float *param_set, gene_form *gene_set,
                          int *coord, int *stable, int *state_file,
                          int *tsl_yn);

int Read_master_input_file(char *name, char name_set[][200], float  *time_set, float  *param_set,
                           gene_form *gene_set, int *coord, int *stable, int *state_file, int *tsl_yn);


int Get_tx_length(char *name);

int Read_RNAP_rates(char *name, rnap_rate_form *rate);

int Read_ribo_rates(char *name, ribo_rate_form *rate);

void Read_invariant_rates(double *mu_rnap_p2, double *mu_rnap_p3,
                          double *mu_rnap_e1, double *mu_rnap_e2,
                          double *mu_rnap_e3, double *mu_rnap_N,
                          double *mu_rnap_c, double *mu_rnap_d,
                          int *min_effector, int *rnap_antipause,
                          int *ribo_antipause, int *ribo_rnap_range,
                          int *max_backtrack, char *name);

void Read_gene_identities(int unit, int *gene_id);

void Read_additional_parameters(char *name);

void Print_assigned_values(float  sim_start, float  sim_stop, float  *seq_window, float  *sample_window,
                           float  k_loading, float  k_on, float  k_off, float  RNA_lifetime,
                           float  prot_lifetime, int num_genes, gene_form *gene_set, int tx_length,
                           char name_set[][200], float  snap_time, float  fract, int num_traj,
                           float generation_time, float b_time, float c_time, float d_time,
                           float k_unloading, float k_to_open, float k_continue_degrade,
                           float p_rnap_rnap_push,float p_ribo_rnap_push,float p_ribo_ribo_push,
                           float k_CSAT, int pol_width, int ribo_width, int degrade_width,
                           float mean_ribo_rate, float min_dwell, int *coord, float birth_vol,
                           int stable_RNA, int read_state, float mean_RNAP_conc, float mean_ribo_conc,
                           char state_file[][200], int tsl_yn);

void Write_rna_ribosome_snapshot(int num_traj, int max_ribo_per_rna,
				 int **RNA_info,int num_mRNA_snapshot,
				 int traj_index);

void Write_rnap_dna_snapshot(int num_promoters, int tx_length,
                             int **rnap_snapshot, int traj_index);

void Write_pseudoseq_files(int tx_length, int *rnap_seq, int *ribo_seq, int traj_index);

void Write_summary_file(char *traj_name, int start_seed, float  sim_stop, int num_traj, float  k_on,
                        float  k_off, float  k_loading, float  *seq_window, float  *sample_window,
                        float  snap_time, int tx_length, float  fract, float  RNA_lifetime,
                        gene_form *gene_set, float  *protein_mean, float  *protein_sd, float *snapshot_ribo_mean,
                        float  snapshot_ribo_sd, float  *mean_ribos_through_rna, float  *sd_ribos_through_rna,
                        float  fish_final[][5], float  fish_final_sd[][5], float  *fish_pool,
                        float  *fish_pool_sd, int num_genes, int traj_index, float  time_step,
                        float  k_continue_degrade, float  p_rnap_rnap_push, float  p_ribo_rnap_push,
                        float  p_ribo_ribo_push, int pol_width, int ribo_width, int degrade_width,
                        float  mean_ribo_rate, float  min_dwell, float  *mean_elongation_rate,
                        int *coord, int stable, int translated, int read_state, float generation_time,
                        float b_time, float c_time, float d_time, float k_unloading, float k_to_open,
                        float k_CSAT, float birth_vol, float mean_RNAP_conc, float mean_ribo_conc,
                        int num_transcripts, int *num_proteins_made, float *mean_proteins_per_cell,
			float *mean_ribo_gene);

void Write_frame(double time, int promoter_index, FILE *fp);

void Prepare_for_no_regulation_simulation_without_file(float *param_set);

void Prepare_for_no_supercoiling_simulation_without_file();

#endif

