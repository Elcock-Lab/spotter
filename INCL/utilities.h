#undef EXTERN
#undef SET_IT

#ifdef DEFINE_UTIL_VARIABLES
#define EXTERN
#define SET_IT = {0}
#else
#define EXTERN extern
#define SET_IT
#endif

#include "sim_types.h"
#include "gen_def.h"


EXTERN int num_to_get;
EXTERN int num_started;
EXTERN int num_tx;
EXTERN int num_tx_finished;                                                     // Used for calculation of tx
EXTERN int num_tx_unfinished;							// elongation rates
EXTERN float *tx_rates;                                                   
EXTERN float *preshut_rates;
EXTERN float *postshut_rates;
EXTERN float *postshut_rates_gene;
EXTERN float *postshut_rates_beyond;
EXTERN float *tx_times;
EXTERN float *final_start_time;
EXTERN float *final_stop_time;
EXTERN float *gene_time;
EXTERN float *final_gene_time;
EXTERN float mean_sd_tx[2];
EXTERN float finished_mean_sd[2];
EXTERN float unfinished_mean_sd[2];
EXTERN float *last_rnap_time;
EXTERN float **three_times;
EXTERN int *num_three_times;
EXTERN int *shutoff_position;
EXTERN int *final_shut_pos;
EXTERN int *final_pos;

EXTERN int radio_incr;								// Global variables for analysis
EXTERN int radio_index;								// of RNA abundances
EXTERN int radio_loc[2];
EXTERN int radio_probe[50000][2];
EXTERN int fish_five_set[20000];
EXTERN int fish_three_set[20000];

EXTERN int fish_five_start;
EXTERN int fish_five_stop;
EXTERN int fish_three_start;
EXTERN int fish_three_stop;

EXTERN int record_on;								// Storage for kymograph data
EXTERN float sim_endtime;							// and rotation visualization
EXTERN int num_kymo_times;
EXTERN float kymo_times[10001];
EXTERN int num_kymo_loc;
EXTERN int kymo_loc[1001];
EXTERN kymo_form kymo_data[10001][1001];
EXTERN int first_rnap_engaged[10001][3];

EXTERN int REPORT_DWT;								// For dwell time densities
EXTERN char dwt_name[100];
EXTERN int dwt_window;
EXTERN int dwt_start,dwt_stop;
EXTERN int num_dwt;
EXTERN float dwt[1000000];
EXTERN int curr_dwt_window[MAX_RNAP];
EXTERN float last_window_entry[MAX_RNAP];

EXTERN int ONE_RNAP_RUN;
EXTERN int RNAP_run_started;
EXTERN int second_load_OK;

EXTERN int SINGMOL_TRACE;							// For making single-molecule
EXTERN char trace_name[100];							// traces
EXTERN int num_trace;
EXTERN int trace_RNAP,trace_pos[MAX_FISH_TIMES][2];
EXTERN int trace_3end[MAX_FISH_TIMES][2];
EXTERN float trace_time[MAX_FISH_TIMES];
EXTERN float acq_freq;

EXTERN int MAKE_KYMOGRAPH;							// Specifications for kymographs
EXTERN int kymo_min_pos;
EXTERN int kymo_max_pos;
EXTERN int kymo_window;
EXTERN float kymo_start;
EXTERN float kymo_duration;
EXTERN float kymo_increment;

EXTERN int frame;								// Storage for data necessary for
EXTERN int sim_total_rnap;							// standard movies
EXTERN int sim_total_ribosomes;
EXTERN int sim_total_rna;

EXTERN int partial_proteins[5000][20];

EXTERN int *active_riboseq;                                                     // Constantly-updated pseudoseq snapshots;
EXTERN int *active_rnapseq;                                                     // Accessed at seq intervals

#ifndef UTIL_FNS
#define UTIL_FNS

void Add_kymograph_info(double time, int current_reaction, int kymo_index);

void Add_FISH_data(double time, int check_stage, int *fish_probe,
                   int fish_total[][MAX_FISH_PROBES],
		   int current_reaction);

void Add_pseudoseq_data(double time, int check_stage, int *rnap_seq,
			int *ribo_seq, int promo_index, int current_reaction,
			FILE *fp);

void Log_nascent_proteins(double time, int current_reaction);

void Determine_translation_efficiencies(int rna_log_ctr,
					int **rna_log,
					int num_genes,
                                        float *mean_ribos_through_rna,
					float *sample_window);

void Determine_protein_production(int *num_proteins_made, int MAX_PROTEINS,
                                  float **proteins_made,
				  float *sample_window,
                                  int num_ten_min_windows,
				  int num_genes, float *protein_mean);


void Take_snapshot(double time, int promoter_index, int traj,
		   int *num_mRNA_snapshot, int MAX_RIBO_PER_RNA,
		   int **RNA_info, float *snapshot_ribo_mean,
		   int **rnap_snapshot, int current_reaction);

#endif
