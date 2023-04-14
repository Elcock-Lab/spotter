#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define DEFINE_SIM_SET_VARIABLES
#include "INCL/sim_settings.C_RNG.h"

#include "INCL/sim_types.C_RNG.h"
#include "INCL/gen_def.C_RNG.h"
#include "INCL/io.C_RNG.h"
#include "INCL/reaction_manager.C_RNG.h"
#include "INCL/replication.C_RNG.h"
#include "INCL/regulation.C_RNG.h"
#include "INCL/transcription.C_RNG.h"
#include "INCL/translation.C_RNG.h"
#include "INCL/supercoiling.C_RNG.h"
#include "INCL/topo.C_RNG.h"
#include "INCL/utilities.C_RNG.h"
#include "INCL/traj.C_RNG.h"


int Generate_trajectory(int num_traj, float sim_stop,
                        float fract, float k_on, float k_off,
                        float k_loading, rnap_rate_form *rnap_sent,
                        ribo_rate_form *ribo_sent, float snap_time,
                        int tx_length, float *seq_window,
                        float *sample_window, gene_form *gene_set,
                        float RNA_lifetime, char *traj_name,
                        int start_seed, int num_genes, int traj_index,
                        float time_step, float k_continue_degrade,
                        float p_rnap_rnap_push, float p_ribo_rnap_push,
                        float p_ribo_ribo_push, int pol_size,
                        int ribo_size, int degrade_width,
                        float mean_ribo_all, float min_dwell,
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
                        float t_promo_shutoff, char checklist[][200])

{


	int i,j,k,traj,ctr;
	int gene_index;
	int object_index;
	int current_reaction=0;
	int reaction_id;
	int result;
	int edge;
	int kymo_step;
	float t;
	int promo_cycle[MAX_PROMO] = {0};
	float promoter_on_log[MAX_PROMO][MAX_PROMO_CYCLES][2];
	float snapshot_ribo_mean[MAX_GENES];
	float snapshot_ribo_sd = 0.0;
	float fish_times[MAX_FISH_TIMES];
	int num_fish_times=0;
	float seq_times[MAX_FISH_TIMES];
	int num_seq_times=0;
	float seq_check = 0.0;
	int fish_probe[MAX_FISH_PROBES];
	int num_ten_min_windows;
	int num_mRNA_snapshot = 0;
	int **rnap_snapshot;
	int *rnap_seq;
	int *ribo_seq;
	float protein_mean[MAX_GENES];
	float protein_sd[MAX_GENES];
	float mean_ribos_through_rna[MAX_GENES];
	float sd_ribos_through_rna[MAX_GENES];

	int num_codons[MAX_GENES];
	float mean_elongation_rate[MAX_GENES];

	int milleresque[10000] = {0};
	float time;

	double curr_time;							// Current simulation time (from last rxn)
        float p_ribo_forward_push;						// Probability of productive collision

	float dwell_adjust;
	float dwell_sum;
	float dwell_loss;
	int dwell_ct;
	float gene_tsl_info[MAX_GENES];

	int *overlap_zone;
	float *adj_factor;

	int num_tsl_info=0;
	int tot_fish_check = 0;
	int fish_total[MAX_FISH_TIMES][MAX_FISH_PROBES]={0};
	float fish_pool[MAX_FISH_PROBES];
	float fish_pool_sd[MAX_FISH_PROBES];
	float fish_check = 0.0;							// Integers (in sec) to deal with possible
	float snap_check = 0.0;							// imprecision in evaluating float  equality
	int rna_log_ctr = 0;							// Logs giving birth times of RNAs and the
										// number of ribosomes total load during life
	int MAX_TOT_RNA;							// Used in defining VLA below
	float *ribo_temp;
	float rna_max;
	float gene_max_ribo;
	float ribo_max=0.0;
	float mean_ribo_rate[MAX_GENES];

	float k_closed_to_open=1.0;						// TEMP: will come from input

	float fish_final[5][5], fish_final_sd[5][5];				// Vestigial: makes compatible with oldstyle

	char gen_output[200];

	translated = tsl_yn;
	rna_max = fract * sim_stop;
	rna_max *= (k_loading/(k_loading+k_unloading));				


	if (k_closed_to_open < 5.0) {						// If stable RNA, this rate is treated as
		rna_max *= k_closed_to_open;					// infinite (cf Liang,...,Dennis, JMB '99)
	}									// and a k_unbind will be chosen to give
	if (stable_RNA == 0 || translated == HYBRID) {				// an appropriate net tx initiation rate
		MAX_TOT_RNA = rna_max * 20;

		MAX_TOT_RNA *= 4.0;

		MAX_PROTEINS = MAX_TOT_RNA * 200.0;				// Pretty conservative...
		if (translated != HYBRID) {
			rna_max *= (RNA_lifetime/sim_stop);					
		}
	}
	else {									// Not necessary for stable transcripts
		MAX_TOT_RNA = 10;
	}
	printf("MAXIMUM TOTAL RNA IS: %d\n",MAX_TOT_RNA);

	rna_max *= 50.0;			

	MAX_RNA = (int) rna_max;						
	if (MAX_RNA < 25) {							
		MAX_RNA = 25;
	}
	if (MAX_RNA > 10000) {							
		MAX_RNA = 10000;
	}
	if (stable_RNA == 1 && translated != HYBRID) {
		MAX_PROTEINS = MAX_RNA;
	}
	
	printf("WILL USE MAX RNA VALUE OF %d\n",MAX_RNA);
	printf("WILL USE MAX PROTEIN VALUE OF %d\n",MAX_PROTEINS);
	rna = malloc(MAX_RNA*sizeof(*rna));
	for (i = 0; i < MAX_RNA; i++) {
		rna[i] = calloc(RNA_INFO,sizeof(*rna[i]));
	}

	for (i = 1; i <= num_genes; i++) {					// Define start, stop, and loading rate
		tsl_start[i] = gene_set[i].start;				// for each gene on the transcription unit
		tsl_stop[i] = gene_set[i].stop;
		if (gene_set[i].k_init > 0) {
			mu_ribo_load[i] = 1/gene_set[i].k_init;
		}
		else {
			mu_ribo_load[i] = 9999999999.9;
		}
		mean_ribo_rate[i] = gene_set[i].mean_rate;
		TSL_ON_OFF[i] = 4 + i;
		LOC_RIBO_LOAD[i] = 27 + i;
		NEWEST_RIBO[i] = 50 + i;
		OLDEST_RIBO[i] = 73 + i;
		TSL_BLOCK[i] = 96 + i;
		NUM_RIBO[i] = 119 + i;
		num_codons[i] = ((tsl_stop[i] - tsl_start[i]) + 1)/3;
		printf("TSL START %d IS %d\n", i, tsl_start[i]);
		gene_max_ribo = mu_ribo_load[i] * (mean_ribo_rate[i]);
		printf("GENE MAX RIBO IS: %f\n",gene_max_ribo);
		gene_max_ribo = ((tsl_stop[i]+1)-tsl_start[i])/
                                gene_max_ribo;
		printf("GENE MAX RIBO IS: %f\n",gene_max_ribo);
		ribo_max += gene_max_ribo;
	}
	MAX_RIBO = (int) ribo_max;
	printf("MAX RIBO IS: %d\n", MAX_RIBO);
	MAX_RIBO_PER_RNA = ((int) ((double) tx_length/30.0)) + 5;
	printf("WILL USE MAX RIBO PER RNA VALUE OF %d\n",
               MAX_RIBO_PER_RNA);
	if (MAX_RIBO < 20) {
		printf("ADJUSTING: MAX RIBO IS: %d\n", MAX_RIBO);
		MAX_RIBO = 20;
	}
	MAX_RIBO *= MAX_RNA;
	if (MAX_RIBO > 20000) {
		MAX_RIBO = 20000;
	}
	if (stable_RNA == 1 && tsl_yn != HYBRID) {
		MAX_RIBO = 5;
	}

	MAX_RIBO = 20000;							// BIGMEM

	printf("WILL USE MAX RIBO VALUE OF %d\n",MAX_RIBO);

	ribosome = malloc(MAX_RIBO*sizeof(*ribosome));
	for (i = 0; i < MAX_RIBO; i++) {
		ribosome[i] = calloc(RIBO_INFO,sizeof(*ribosome[i]));
	}

	for (i = 0; i < MAX_RNAP; i++) {
		rnap_search[i] = &rnap[i][0];
	}
	for (i = 0; i < MAX_TOPO; i++) {
		topo_search[i] = &topo[i][0];
	}

	memset(num_topo,0,sizeof(num_topo));
	memset(topo_master_list,0,sizeof(topo_master_list));
	memset(curr_topo_list,0,sizeof(curr_topo_list));
	memset(starter_up,0,sizeof(starter_up));
	memset(starter_down,0,sizeof(starter_down));

	for (i = 1; i < MAX_PROMO; i++) {
		topo_finder[i] = &topo_master_list[i][0][0];
	}
	Load_burst_distribution();

	num_upstream_update = 0;
	num_downstream_update = 0;
	memset(upstream_list,0,sizeof(upstream_list));
	memset(downstream_list,0,sizeof(downstream_list));

        if (MAX_RNA > 10000) {
                MAX_RNA = 10000;
                printf("SETTING MAXIMUM NUMBER OF RNA EXTANT IN"
                        "SIMULATION TO CEILING OF 10000\n");
        }
        if (MAX_TOT_RNA > 100000) {
                MAX_TOT_RNA = 100000;
                printf("SETTING MAXIMUM NUMBER OF RNA MADE THROUGHOUT"
                        "SIMULATION TO CEILING OF 100000\n");
        }
        if (MAX_PROTEINS > 2000000) {
                MAX_PROTEINS = 2000000;
                printf("SETTING MAXIMUM NUMBER OF PROTEINS CREATED"
                        "THROUHGOUT SIMULATION TO CEILING OF 2000000\n");
        }

        int **RNA_info;
        int **rna_log;
        float **proteins_made;

        RNA_info = malloc(MAX_RNA*sizeof(*RNA_info));
        for (i = 0; i < MAX_RNA; i++) {
                RNA_info[i] = calloc(MAX_RIBO_PER_RNA,sizeof(*RNA_info[i]));
        }
        rna_log = malloc(MAX_TOT_RNA*sizeof(*rna_log));
        for (i = 0; i < MAX_TOT_RNA; i++) {
                rna_log[i] = calloc(MAX_GENES,sizeof(*rna_log[i]));
        }
        proteins_made = malloc((MAX_GENES+1)*sizeof(*proteins_made));
        for (i = 0; i < MAX_GENES+1; i++) {
                proteins_made[i] =
                   calloc(MAX_PROTEINS,sizeof(*proteins_made[i]));
        }

	tx_size = tx_length+50;
	dna_strip = malloc(MAX_PROMO*sizeof(*dna_strip));				
	rnap_snapshot = malloc(MAX_PROMO*sizeof(*rnap_snapshot));				
	for (i = 0; i < MAX_PROMO; i++) {						
		dna_strip[i] = calloc((tx_size+1),sizeof(*dna_strip[i]));
		rnap_snapshot[i] = calloc((tx_size+1),
                                          sizeof(*rnap_snapshot[i]));
	}
	rna_strip = malloc(MAX_RNA*sizeof(*rna_strip));
	for (i = 0; i < MAX_RNA; i++) {
		rna_strip[i] = calloc((tx_size+1),sizeof(*rna_strip[i]));
		for (j = 1; j <= tx_size; j++) {
			rna_strip[i][j] = -1;
		}
	}
	mu_rnap_move = malloc((tx_size+1)*sizeof(*mu_rnap_move));
	mu_ribo_move = malloc((tx_size+1)*sizeof(*mu_ribo_move));
	rna_block_grid = calloc((tx_size+1),sizeof(*rna_block_grid));
	active_riboseq = calloc((tx_size+1),sizeof(*active_riboseq));
	active_rnapseq = calloc((tx_size+1),sizeof(*active_rnapseq));
	rnap_seq = calloc((tx_size+1),sizeof(*rnap_seq));
	ribo_seq = calloc((tx_size+1),sizeof(*ribo_seq));
	ribo_temp = malloc((tx_size+1)*sizeof(*ribo_temp));

	num_to_get = 200;
	num_tx = 0;
	num_started = 0;
	tx_times = calloc(MAX_RNAP,sizeof(*tx_times));
	gene_time = calloc(MAX_RNAP,sizeof(*gene_time));
	shutoff_position = calloc(MAX_RNAP,sizeof(*shutoff_position));
	tx_rates = calloc(MAX_TOT_RNA,sizeof(*tx_rates));
	final_shut_pos = calloc(MAX_TOT_RNA,sizeof(*final_shut_pos));
	final_pos = calloc(MAX_TOT_RNA,sizeof(*final_shut_pos));
	preshut_rates = calloc(MAX_TOT_RNA,sizeof(*preshut_rates));
	postshut_rates = calloc(MAX_TOT_RNA,sizeof(*postshut_rates));
	postshut_rates_gene =
	   calloc(MAX_TOT_RNA,sizeof(*postshut_rates_gene));
	postshut_rates_beyond =
	   calloc(MAX_TOT_RNA,sizeof(*postshut_rates_beyond));
	final_start_time =
	   calloc(MAX_TOT_RNA,sizeof(*final_start_time));
	final_stop_time =
	   calloc(MAX_TOT_RNA,sizeof(*final_stop_time));
	final_gene_time =
	   calloc(MAX_TOT_RNA,sizeof(*final_gene_time));
	memset(mean_sd_tx,0,sizeof(mean_sd_tx));
	memset(finished_mean_sd,0,sizeof(finished_mean_sd));
	memset(unfinished_mean_sd,0,sizeof(unfinished_mean_sd));

	num_three_times = calloc((tx_size+1),sizeof(*num_three_times));
	three_times = malloc((tx_size+1)*sizeof(*three_times));
	for (i = 0; i < tx_size + 1; i++) {
		three_times[i] = calloc(20001,sizeof(*three_times[i]));
	}
	last_rnap_time = calloc(MAX_RNAP,sizeof(*last_rnap_time));

	radio_incr = 10;
	radio_loc[0] = 433;
	radio_loc[1] = 3107;
	memset(radio_probe,0,sizeof(radio_probe));
	memset(fish_five_set,0,sizeof(fish_five_set));
	memset(fish_three_set,0,sizeof(fish_three_set));

        num_final_rna = 0;
        rna_to_log = calloc(MAX_RNA,sizeof(*rna_to_log));
        final_rna_log = calloc(MAX_TOT_RNA,sizeof(rna_log_form));
		 
	overlap_zone = calloc((tx_size+1),sizeof(*overlap_zone));
	adj_factor = malloc((tx_size+1)*sizeof(*adj_factor));

	gene_num = num_genes;							// Trivial renamings to allow
	tx_end = tx_length;							// global use

	start_pos = coord[0];
	stop_pos = coord[1];
	stable = stable_RNA;
	if (translated == HYBRID) {
		split = (tsl_stop[gene_num-1]+tsl_start[gene_num])/2;
	}	

	RNAP_run_started = 0;
	ONE_RNAP_RUN = 0;
	second_load_OK = 0;

	SINGMOL_TRACE = 0;
	trace_RNAP = 0;
	num_trace = 0;
	acq_freq = 5.0;

	MAKE_KYMOGRAPH = 0;
	kymo_min_pos = 100;
	kymo_max_pos = 450;
	kymo_window = 500;
	kymo_start = 50.0;
	kymo_duration = 50.0;
	kymo_increment = 0.01;

	mu_promo_on = 1/k_on;
	mu_promo_off = 1/k_off;
	mu_promo_load = 1/k_loading;						// Currently set under assumptions that
	mu_unbind_rnap = 1/k_unloading;						// closed->open rate is 1 s^-1, fractional
	orig_mu_closed_to_open = 1/k_to_open;					// promoter occupancy is 0.45, so that
	basal_promo_load = mu_promo_load;					// the effective rate of k_tx is 0.45 s^-1
	memcpy(basal_ribo_load,mu_ribo_load,sizeof(mu_ribo_load));		// k_on for RNAP is 45 s^-1 (input file)
										// See notes and timestep version for
										// justification
	mu_start_degrade = RNA_lifetime;
	mu_continue_degrade = 1/k_continue_degrade;
	if (k_CSAT < 0.000001) {						// Note that at present csat comes from
		mu_ribo_csat = INFINITY;					// p_ribo_ribo; this could be hijacked
	}									// if necessary
	else {
		mu_ribo_csat = 1/k_CSAT;
	}

	p_ribo_forward_push = p_ribo_ribo_push;

	ribo_dwell = ribo_sent;
	rnap_dwell = rnap_sent;
	

	for (i = 1; i <= num_genes; i++) {
		mean_ribo_rate[i] = 1.0/mean_ribo_rate[i];
		printf("MEAN RIBO RATE: %f\n",mean_ribo_rate[i]);
	}
	memset(gene_tsl_info,0,sizeof(gene_tsl_info));
	for (i = 1; i <= num_genes; i++) {
		dwell_ct = 0;
		for (j = tsl_start[i]; j <= tsl_stop[i]; j+= 3) {
			ribo_temp[j] = 0.0;
			for (k = 0; k < 3; k++) {
				ribo_temp[j] += ribo_dwell[j+k].exp;
			}
			gene_tsl_info[i] += ribo_temp[j];
			dwell_ct++;
		}
		gene_tsl_info[i] /= ((double) dwell_ct);
	}

	for (i = 1; i <= num_genes; i++) {
		dwell_sum = 0.0;
		dwell_loss = 0.0;
		dwell_adjust = gene_tsl_info[i];
		dwell_adjust /= (mean_ribo_rate[i]*3.0);
		for (j = tsl_start[i]; j <= tsl_stop[i]; j+= 3) {
			ribo_temp[j] /= dwell_adjust;
			if (ribo_temp[j] < min_dwell) {
				dwell_loss += (min_dwell-ribo_temp[j]);
			}
			else {
				dwell_sum += ribo_temp[j];
			}
		}
		dwell_sum = (dwell_sum - dwell_loss)/dwell_sum;
		for (j = tsl_start[i]; j <= tsl_stop[i]; j+= 3) {
			if (ribo_temp[j] < min_dwell) {
				ribo_temp[j] = min_dwell;
			}
			else {
				ribo_temp[j] *= dwell_sum;
			}
		}
		for (j = tsl_start[i]; j <= tsl_stop[i]; j+= 3) {
			if (ribo_temp[j] > ribo_dwell[j].codon - min_dwell) {
				ribo_dwell[j].mu_add = ribo_dwell[j].codon -
						       min_dwell;
				ribo_dwell[j].mu_f = ribo_temp[j] - 
						     ribo_dwell[j].mu_add;
			}
			else {
				ribo_dwell[j].mu_f = min_dwell;
				ribo_dwell[j].mu_add = ribo_temp[j] - 
						       ribo_dwell[j].mu_f;
				if (ribo_dwell[j].mu_add < 0.0) {
					ribo_dwell[j].mu_add = 0.00001;
				}
			}
		}
						     
	}
	for (i = 1; i <= num_genes; i++) {
		gene_tsl_info[i] /= mean_ribo_rate[i];
	}

	for (i = 1; i <= tx_length; i++) {
		rnap_dwell[i].p_on_f  = rnap_dwell[i].mu_f /			// Enabling sampling
				        (rnap_dwell[i].mu_f +			// from uniform, not
					 rnap_dwell[i].mu_p);			// exp distrib;
		rnap_dwell[i].p_off_f = rnap_dwell[i].mu_f_off /		// NB: passed as
					(rnap_dwell[i].mu_f_off +		// *rates*
					 rnap_dwell[i].mu_b_off);
		rnap_dwell[i].mu_net_pre = rnap_dwell[i].mu_f +
					   rnap_dwell[i].mu_p;
		rnap_dwell[i].mu_net_off = rnap_dwell[i].mu_f_off +
				  	   rnap_dwell[i].mu_b_off; 
		rnap_dwell[i].mu_f = 1/rnap_dwell[i].mu_f;			// Convert rates
		rnap_dwell[i].mu_b = 1/rnap_dwell[i].mu_b;			// to dt's
		rnap_dwell[i].mu_p = 1/rnap_dwell[i].mu_p;
		rnap_dwell[i].mu_f_off = 1/rnap_dwell[i].mu_f_off;
		rnap_dwell[i].mu_b_off = 1/rnap_dwell[i].mu_b_off;
		rnap_dwell[i].mu_net_pre = 1/rnap_dwell[i].mu_net_pre;
		rnap_dwell[i].mu_net_off = 1/rnap_dwell[i].mu_net_off;
	}

	mu_pause[OFF_P1] = mu_rnap_p2;						// Inv rate P1->P2
	mu_pause[OFF_P2] = mu_rnap_p3;						// Inv rate P2->P3
	mu_exit[OFF_P1] = mu_rnap_e1;						// Exit from P1
	mu_exit[OFF_P2] = mu_rnap_e2;						// Exit from P2
	mu_exit[OFF_P3] = mu_rnap_e3;						// Exit from P3
	mu_rnap_add = mu_rnap_c;						// Nt addition

	orig_p_closed_to_open = k_to_open/(k_to_open+k_unloading);
	p_p1_to_p2 = (1/mu_rnap_p2)/((1/mu_rnap_p2)+(1/mu_rnap_e1)); 
	p_p2_to_p3 = (1/mu_rnap_p3)/((1/mu_rnap_p3)+(1/mu_rnap_e2));
	mu_net_promo = 1/(k_to_open + k_unloading);
	mu_net_p1 = 1/((1/mu_rnap_p2)+(1/mu_rnap_e1));
	mu_net_p2 = 1/((1/mu_rnap_p3)+(1/mu_rnap_e2));

	EFFECTOR_STATE = min_effector;
	RNAP_ANTIPAUSE = rnap_antipause;
	RIBO_ANTIPAUSE = ribo_antipause;
	RIBO_RNAP_RANGE = ribo_rnap_range;
	MAX_BACKTRACK = max_backtrack;

	BACKTRACK_ON = 1;		

	TETHER_LENGTH = 5;
	
	printf("P P1 -> P2: %f\n",p_p1_to_p2);
	printf("P P2 -> P3: %f\n",p_p2_to_p3);
	printf("NET MU STATE P1: %f\n",mu_net_p1);
	printf("NET MU STATE P2: %f\n",mu_net_p2);
	printf("MU EXIT STATE P3: %f\n",mu_exit[OFF_P3]);

	promo_kill = t_promo_shutoff;
	warm_up_time = 0.0;
	repl_time = 99999.9;

	memset(no_super_pause,0,sizeof(no_super_pause));

	BUBBLE_ADJ = OFF;
	bubble[UPSTREAM]   = 0;
	bubble[DOWNSTREAM] = 0;

	DNA_SCRUNCH = OFF;
	TSS_adj = 0;

	rnap_relax = OFF;
	writhe_partition = OFF;
	supercoiling_off = OFF;
	ribo_resist_model = IMPLICIT_RNAP_POS;
	mean_interribo_dist = 100.0;
	ANTICASCADE_ON = OFF;
	NO_PAUSE_MODEL = OFF;
	BACKTRACK_ON = OFF;
	QUIET_MODE = OFF;
	max_RNAP_load_gap = 9999999;
	resist_exp = 1.5;

	memset(lk,0,sizeof(lk));
	rot_resist = malloc((tx_size+1)*sizeof(*rot_resist));

	min_resist = 1.0;
	max_resist = 1.0;

	k = 5 - ((strcmp(checklist[6],"NOT SUPPLIED") != 0) * 4);
	for (i = 1; i <= k; i++) {
	   if (strcmp(checklist[i],"NOT SUPPLIED") != 0) {
		Read_additional_parameters(checklist[i]);
	   }
	   else {
		if (i == 3) {
		   Prepare_for_no_supercoiling_simulation_without_file();
		}
	   }
	}
	Assign_rotation_resistance();

	if (NO_PAUSE_MODEL) {
		for (i = 1; i <= tx_length; i++) {
			rnap_dwell[i].mu_f = rnap_dwell[i].net;                 // FOR single-step model
		}								// without pauses
	}			

	printf("SOME STUFF YOU MIGHT WANT TO KNOW:\n");
	printf("WRITHE PARITITONING STATUS (1/0 = ON/OFF): %d\n",
		writhe_partition);
	printf("NO SUPERCOILING STATUS (1/0 = ON/OFF): %d\n",
		supercoiling_off);
	printf("RNAP COMPLEX RELAXATION STATUS (1/0 = ON/OFF): %d\n",
		rnap_relax);
	printf("RELAXATION FACTOR (MU FACTOR): %f\n",relax_factor);
	printf("RIBOSOME RESISTANCE MODEL: %d\n",ribo_resist_model);
	printf("MEAN INTERRIBO DIST (IMPLICIT MODEL): %f\n",
		mean_interribo_dist);
	printf("PAUSE-FREE MODEL (1/0 = ON/OFF): %d\n",NO_PAUSE_MODEL);
	printf("ANTI-PAUSE CASCADE (1/0 = ON/OFF): %d\n",ANTICASCADE_ON);
	printf("BACKTRACKING MODEL (1/0 = ON/OFF): %d\n",BACKTRACK_ON);
	printf("FLAT RESISTANCE PARAMETER: %f\n",flat_resist_param);
	printf("BT ADJ: %f\n",bt_adj);
	
	sim_stop += warm_up_time;
	promo_kill += warm_up_time;

	sim_endtime = sim_stop;

	for (i = fish_five_start; i <= fish_five_stop; i += 40) {
		fish_five_set[i] = 1;
	}
	for (i = fish_three_stop; i >= fish_three_start; i -= 40) {
		fish_three_set[i] = 1;
	}


	pseudo_shutoff = 0;
	if (promo_kill > sim_stop) {
		pseudo_shutoff++;
		promo_kill = t_pseudo_shut;
	}

	if (plasmid_linearized) {						// Override if necessary:
		start_sigma = 0.0;						// linearized always
	}									// starts relaxed

	stop_barrier += tx_length;
	starting_helical_repeat = (1.0 - (0.3*start_sigma))*10.5;
	for (i = 1; i < 15000000; i++) {
		relaxed_lk[i] = ((float) i)/10.5;
	}
	shear_modulus = 300;							//pN/nm^2, per Heberling et al.
	torque_factor = (shear_modulus*PI*PI)/10.5;
	start_lk_up = (1.0 + start_sigma)*
	              (((float) (1 - start_barrier))/10.5);
	start_lk_down = (1.0 + start_sigma)*
	              (((float) (stop_barrier - 1))/10.5);
	printf("STARTING SIGMA IS %f\n",start_sigma);
	printf("****ASSIGNING START LK'S OF %f {UP} AND %f {DOWN}\n",
		start_lk_up,start_lk_down);

	relaxed_open_lk = (((float) (1 - start_barrier))/10.5);
	relaxed_open_lk += (((float) (stop_barrier - 1))/10.5);
	printf("RELAXED OPEN LK: %f\n",relaxed_open_lk);

	for (i = 1; i < MAX_PROMO; i++) {
		init_sigma[i] = start_sigma;
		curr_init_lk_up[i] = start_lk_up;
		curr_init_lk_down[i] = start_lk_down;
		open_lk[i] = start_lk_up + start_lk_down;
		printf("PROMO %d: INITIAL SIGMA %f; CURR LK UP %f; "
		       "CURR LK DOWN: %f; STARTING OPEN LK: %f\n",
			i,init_sigma[i],curr_init_lk_up[i],
			curr_init_lk_down[i],open_lk[i]);
	}

	promo_barrier_left[OFF]  = start_barrier;
	promo_barrier_right[OFF] = stop_barrier;

	for (i = 1; i < MAX_PROMO; i++) {
		p_closed_to_open[i] = orig_p_closed_to_open;
		mu_closed_to_open[i] = orig_mu_closed_to_open;
	}

	torque_zero_pause_p = 0.0128;						// Per Heberling et al.
	torque_zero_pause_dur = 0.55812;
	rnap_p1_effect = 0.01;
	ribo_p1_effect = 0.01;
	rnap_e1_effect = 0.01;
	ribo_e1_effect = 0.01;
	

	memset(fish_times, 0, sizeof(fish_times));
	memset(seq_times, 0, sizeof(seq_times));
	memset(discard,0,sizeof(discard));
	memset(mean_elongation_rate,0,sizeof(mean_elongation_rate));
	memset(mean_proteins_per_cell,0,sizeof(mean_proteins_per_cell));
	memset(snapshot_ribo_mean,0,sizeof(snapshot_ribo_mean));
	
	pol_width = pol_size;							// Set clash-checking widths
	ribo_width = ribo_size;
	pol_ribo_width = (int) (((float ) pol_width +
		(float ) ribo_width)/2.0 + 0.5);
	degrade_rnap_width = (pol_width/2) + degrade_width;
	degrade_ribo_width = (ribo_width/2) + degrade_width;			// Room for degradation apparatus--
										// can change

	for (i = 1; i <= num_genes; i++) {
		edge = tsl_start[i] - ribo_width;
		if (edge < 1) {
			edge = 1;
		}
		for (j = edge; j <= tsl_start[i] + ribo_width; j++) {
			rna_block_grid[j] = i;
		}
	}
		
	fish_probe[1] = 1;
	fish_probe[2] = tx_length/2;
	fish_probe[3] = tx_length;

        record_on = 0;
        ctr = 0;
        t = 0.0;
	kymo_step = (int) (kymo_duration/kymo_increment);
        for (i = 1; i <= kymo_step; i++) {
                t += kymo_increment;
                ctr++;
                kymo_times[ctr] = t;
        }
        num_kymo_times = ctr;

        ctr = 0;
        for (i = 1; i <= kymo_window; i++) {
                ctr++;
                kymo_loc[ctr] = i;
        }
        num_kymo_loc = ctr;
        for (i = 1; i <= num_kymo_times; i++) {
                for (j = 1; j <= num_kymo_loc; j++) {
                        kymo_data[i][j].RNAP_bound = 0;
                        kymo_data[i][j].RNAP_state = 0;
                        kymo_data[i][j].RNAP_rot = 0;
                }
        }
        memset(first_rnap_engaged,0,sizeof(first_rnap_engaged));
	memset(curr_rnap_rot,0,sizeof(curr_rnap_rot));

	ctr = 0;
	t = sample_window[0];
	for (i = 1; i < 100; i++) {
		if (t > sample_window[1]) {
			break;
		}
		ctr++;
		fish_times[ctr] = t;
		t += 60.0;
	}
	num_fish_times = ctr;
	num_ten_min_windows = ((int) (sample_window[1] + 0.001));
	num_ten_min_windows -= ((int) (sample_window[0] + 0.0001));
	num_ten_min_windows /= 600;

	ctr = 0;
	t = seq_window[0];
	for (i = 1; i < 50000; i++) {
		if (t > seq_window[1]) {
			break;
		}
		ctr++;
		seq_times[ctr] = t;
		t += (1.0/acq_freq);
	}
	num_seq_times = ctr;

	traj = 1;
	num_traj = 1;

	temp_off = sim_stop+1.0;
	Initialize_reaction_times();
	Assign_reaction_information();
        printf("         %d:%f  QUEUE SIZE %d\n",
		reaction_queue[1].reaction,
		reaction_queue[1].reaction_time,queue_size);
        printf("  %d:%f        %d:%f\n",
		reaction_queue[2].reaction,
		reaction_queue[2].reaction_time,
		reaction_queue[3].reaction,
		reaction_queue[3].reaction_time);
        printf("%d:%f %d:%f %d:%f %d:%f\n",reaction_queue[4].reaction,
		reaction_queue[4].reaction_time,
		reaction_queue[5].reaction,
		reaction_queue[5].reaction_time,
		reaction_queue[6].reaction,
		reaction_queue[6].reaction_time,
		reaction_queue[7].reaction,
		reaction_queue[7].reaction_time);

	birth_promo = 1;
	num_promoters = birth_promo;
	Start_from_zero(fract);


	printf("REACTION ZERO IS %d:%f\n",
		reaction_queue[0].reaction,
		reaction_queue[0].reaction_time);
        printf("         %d:%f  QUEUE SIZE %d\n",
		reaction_queue[1].reaction,
		reaction_queue[1].reaction_time,queue_size);
        printf("  %d:%f        %d:%f\n",
		reaction_queue[2].reaction,
		reaction_queue[2].reaction_time,
		reaction_queue[3].reaction,
		reaction_queue[3].reaction_time);
        printf("%d:%f %d:%f %d:%f %d:%f\n",reaction_queue[4].reaction,
		reaction_queue[4].reaction_time,
		reaction_queue[5].reaction,
		reaction_queue[5].reaction_time,
		reaction_queue[6].reaction,
		reaction_queue[6].reaction_time,
		reaction_queue[7].reaction,
		reaction_queue[7].reaction_time);

	Add_pseudoreactions_to_queue(num_fish_times, fish_times,
					num_seq_times, seq_times,
					snap_time, repl_time);

	curr_time = 0.000;
	ctr = 1;

	for (i = 0; i < MAX_RNAP; i++) {
		if (ribo_resist_model == FLAT_RESISTANCE) {
			complex_resist[i] = flat_resist_param;
		}
		else {
			complex_resist[i] = 1.0/(0.12 + 0.05);			// RNAP from Stokes; DNA from
		}								// Tripathi (0.05 chi)
	}


	dna_length = stop_barrier - start_barrier;
	mu_topo_binding[TOPO_IA] = 1.0/rate_topoI;
	mu_topo_binding[GYRASE] = 1.0/rate_gyrase;


	topo_rnap_width[TOPO_IA] = ((pol_width + topoIA_width)/2)+1;
	topo_rnap_width[GYRASE] = ((pol_width + gyrase_width)/2)+1;
	topo_dist_check[TOPO_IA + TOPO_IA] = topoIA_width;
	topo_dist_check[TOPO_IA + GYRASE] = (topoIA_width+gyrase_width)/2;
	topo_dist_check[GYRASE + GYRASE] = gyrase_width;
	memset(promo_topo_blocked,0,sizeof(promo_topo_blocked));


	Add_reaction_to_queue(TOPO_I_BINDING_RXN(1),
	   Find_firing_time(mu_topo_binding[TOPO_IA]));
	Add_reaction_to_queue(GYRASE_BINDING_RXN(1),
	   Find_firing_time(mu_topo_binding[GYRASE]));

	memset(partial_proteins,0,sizeof(partial_proteins));
	shutoff_done= 0;

	frame = 0;
	sim_total_rna = 0;
	sim_total_rnap = 0;
	sim_total_ribosomes = 0;
	FILE *fp_for_xtc;
	fp_for_xtc = fopen("tx_tsl_traj.for_xtc.txt","w");

	printf("FISH CHECK AT %f; SEQ CHECK AT %f; SNAPSHOT AT %f\n",
               fish_check,seq_check,snap_check);

	failed_loads = 0;

	say_it = 0;

	printf("BEGINNING SIMULATION\n");

	while (curr_time < sim_stop) {
		curr_time = reaction_queue[1].reaction_time;
		current_reaction = reaction_queue[1].reaction;
		reaction_id = reactions[current_reaction].reaction_id;
		object_index = reactions[current_reaction].object_index;
		if (reaction_id == PROMO_ON) {
			result = Activate_promoter(curr_time,promo_cycle,
						   object_index,
                                                   promoter_on_log);
		}
		else if (reaction_id == PROMO_OFF) {
			result = Deactivate_promoter(curr_time,
						     promo_cycle,
						     object_index,
                                                     promoter_on_log, 0);
		}
		else if (reaction_id == LOAD_RNAP) {
			result = Load_RNAP(curr_time,object_index);
		}
		else if (reaction_id == LOAD_RIBO) {
			gene_index =
			 reactions[current_reaction].gene_index;
			result =
			 Load_ribosome(curr_time,gene_index,object_index);
		}
		else if (reaction_id == RNAP_F_TRANSL) {
			result = Attempt_RNAP_move(curr_time, 
						   object_index,
                                                   p_rnap_rnap_push,
						   MAX_PROTEINS,
                                                   proteins_made,
						   FORWARD);
		}
		else if (reaction_id == RNAP_B_TRANSL) {
			result = Attempt_RNAP_move(curr_time,
						   object_index,
                                                   p_rnap_rnap_push,
						   MAX_PROTEINS,
                                                   proteins_made,
						   REVERSE);
		}
		else if (reaction_id == RNAP_NT_ADD) {
			result = Add_RNA_nt(curr_time, object_index,
					    MAX_PROTEINS,proteins_made);
		}
		else if (reaction_id == RNAP_OFF_FOR) {
			result = Attempt_RNAP_move(curr_time,
						   object_index,
                                                   p_rnap_rnap_push,
						   MAX_PROTEINS,
                                                   proteins_made,
						   FORWARD);
		}
		else if (reaction_id == RNAP_OFF_REV) {
			result = Attempt_RNAP_move(curr_time,
						   object_index,
                                                   p_rnap_rnap_push,
						   MAX_PROTEINS,
                                                   proteins_made,
						   REVERSE);
		}
		else if (reaction_id == RNAP_ENTER_PAUSE) {
			result = Advance_pause_state(
					curr_time,object_index);
		}
		else if (reaction_id == RNAP_EXIT_PAUSE) {
			result = Exit_off_pathway_state(
					curr_time,object_index);
		}


		else if (reaction_id == RNAP_TERMINATION) {
			result = Terminate_transcription(
					curr_time,object_index);
		}

		else if (reaction_id == RNAP_SURVEILLANCE) {
			result = Schedule_RNAP_checkup(
					curr_time,object_index);
		}
		else if (reaction_id == RNAP_EVALUATION) {
			result = Evaluate_transcriptional_progress(
					curr_time,object_index);
		}
		else if (reaction_id == RNAP_RELAX) {
			result = Relax_RNAP_complex(
					curr_time,object_index);
		}
			
		else if (reaction_id == RIBO_DECODE) {
			result = Add_amino_acid(curr_time, object_index);
		}
		else if (reaction_id == RIBO_TRANSLOC) {
			result = Attempt_ribosome_move(curr_time,
						       object_index,
                                                       p_ribo_forward_push,
                                                       p_ribo_rnap_push,
                                                       p_rnap_rnap_push,
                                                       num_proteins_made,
						       MAX_PROTEINS,
                                                       proteins_made,
						       num_codons,
                                                       mean_elongation_rate);
		}
		else if (reaction_id == COLLISION_ABORT) {
			result = Remove_stalled_ribosome(
					curr_time,object_index);
		}
		else if (reaction_id == START_DEGRADE) {
			result = Start_RNA_degradation(
					curr_time,object_index);
		}
		else if (reaction_id == CONT_DEGRADE) {
			result = Continue_RNA_degradation(
					curr_time,object_index,
                                        &rna_log_ctr,rna_log);
		}
		else if (reaction_id == ACTIVATE_RNAP) {
			result = Activate_RNAP(curr_time,object_index);
		}
		else if (reaction_id == REMOVE_RNAP) {
			result = Remove_RNAP(curr_time,object_index);
		}
		else if (reaction_id == ADD_FISH_DATA) {
			tot_fish_check++;
			Add_FISH_data(curr_time,tot_fish_check,fish_probe,
				      fish_total,current_reaction);		
		}
		else if (reaction_id == ADD_SEQ_DATA) {
			Add_pseudoseq_data(curr_time,ctr,rnap_seq,ribo_seq,
                                           1,current_reaction,fp_for_xtc);
		}
		else if (reaction_id == TAKE_SNAPSHOT) {
			Take_snapshot(curr_time,1,traj,&num_mRNA_snapshot,
                                      MAX_RIBO_PER_RNA, RNA_info,
				      snapshot_ribo_mean, rnap_snapshot,
				      current_reaction);
			Write_rna_ribosome_snapshot(num_traj,
						    MAX_RIBO_PER_RNA,
						    RNA_info,
                                    		    num_mRNA_snapshot,
						    traj_index);
			Write_rnap_dna_snapshot(num_promoters,tx_length,
						rnap_snapshot,traj_index);
			num_mRNA_snapshot = 0;
			memset(snapshot_ribo_mean, 0, 
				sizeof(snapshot_ribo_mean));
		}
		else if (reaction_id == KILL_PROMOTERS) {
			if (pseudo_shutoff) {
				Get_position_without_shutoff();
			}
		}
		else if (reaction_id == TOPO_I_BINDING) {
			Attempt_topoisomerase_binding(curr_time,TOPO_IA,1,
						topo_rnap_width[TOPO_IA]);
		}
		else if (reaction_id == GYRASE_BINDING) {
			Attempt_topoisomerase_binding(curr_time,GYRASE,1,
						topo_rnap_width[GYRASE]);
		}
		else if (reaction_id == TOPOISOMERASE_UNBINDING) {
			Unbind_topoisomerase_on_schedule(curr_time,
							 object_index);
		}
		else if (reaction_id == TOPOISOMERASE_CLEAVAGE) {
			Perform_topoisomerase_cycle(curr_time,
						    object_index);
		}
		else if (reaction_id == TOPOISOMERASE_STATE_CHANGE) {
			Switch_burst_state(curr_time,object_index);
		}
                else if (reaction_id == KYMOGRAPH) {
                        Add_kymograph_info(curr_time,current_reaction,
                                           object_index);
                }
		else if (reaction_id == LOG_NASCENT_PROTEINS) {
			Log_nascent_proteins(curr_time,current_reaction);
		}
				
		sim_stop = sim_endtime;
		result++;
	}		

// POST-SIM PROESSING

	printf("FINISHED WITH SIMULATIION\n");
	memset(mean_ribos_through_rna,0,sizeof(mean_ribos_through_rna));
	memset(sd_ribos_through_rna,0,sizeof(sd_ribos_through_rna));
	memset(protein_mean,0,sizeof(protein_mean));
	memset(protein_sd,0,sizeof(protein_sd));
	memset(fish_pool,0,sizeof(fish_pool));
	memset(fish_pool_sd,0,sizeof(fish_pool_sd));

	printf("LOGGED %d COMPLETE RNAS\n",rna_log_ctr);
	for (i = 1; i <= rna_log_ctr; i++) {
		printf("%10d", rna_log[i][0]);
		for (j = 1; j <= num_genes; j++) {
			printf("%10d",rna_log[i][j]);
		}
		printf("\n");
	}
	printf("RELATIVE RIBOSEQ-DERIVED GENE TRANSLATION:\n");
	for (i = 1; i <= num_tsl_info; i++) {
		printf("GENE%5d:%10.5f\n",i,gene_tsl_info[i]);
	}

	printf("CHECKING SOME STUFF...\n");
	printf("TOTAL FISH CHECKS: %d\n",tot_fish_check);
	printf("FISH TOTALS FIRST: %d %d %d\n",
		fish_total[1][1],fish_total[1][2],fish_total[1][3]);
	printf("FISH TOTALS PENULT: %d %d %d\n",
		fish_total[tot_fish_check - 1][1],
		fish_total[tot_fish_check - 1][2],
		fish_total[tot_fish_check - 1][3]);
	printf("FISH TOTALS LAST: %d %d %d\n",
		fish_total[tot_fish_check][1],
		fish_total[tot_fish_check][2],
		fish_total[tot_fish_check][3]);

	Determine_translation_efficiencies(rna_log_ctr, rna_log,
					   num_genes,mean_ribos_through_rna,
					   sample_window);

	Determine_protein_production(num_proteins_made, MAX_PROTEINS,
				     proteins_made, sample_window,
				     num_ten_min_windows,num_genes,
                                     protein_mean);
	
	for (i = 1; i <= num_genes; i++) {
		mean_elongation_rate[i] /=
			(num_proteins_made[i] - init_proteins[i]);
	}
	for (i = 1; i <= tot_fish_check; i++) {
		for (j = 1; j <= 3; j++) {
			fish_pool[j] += ((float ) fish_total[i][j]);
		}
	}
	for (i = 1; i <= 3; i++) {
		fish_pool[i] /= ((float ) tot_fish_check);
	}
	for (i = 1; i <= tot_fish_check; i++) {
		for (j = 1; j <= 3; j++) {
			fish_pool_sd[j] += ((fish_total[i][j] -
					     fish_pool[j]) *
                                            (fish_total[i][j] -
					     fish_pool[j]));
		}
	}
	for (i = 1; i <= 3; i++) {
		fish_pool_sd[i] /= ((float ) tot_fish_check);
		fish_pool_sd[i] = sqrt(fish_pool_sd[i]);
	}

	num_tx_finished = num_tx;

	for (i = 1; i < MAX_RNAP; i++) {
		if (rnap[i][POSITION] && rnap[i][STATE] != CLOSED) {
			num_tx++;
			final_pos[num_tx] = rnap[i][POSITION];
	                tx_rates[num_tx] = ((double) rnap[i][POSITION]) /
	                                   (sim_stop -
					    tx_times[i]);
			final_shut_pos[num_tx] = shutoff_position[i];
			final_start_time[num_tx] = tx_times[i];
			final_gene_time[num_tx] = gene_time[i];
			final_stop_time[num_tx] = sim_stop;
		}
	}

	for (i = 1; i <= num_tx; i++) {
		mean_sd_tx[0] += tx_rates[i];
		if (i <= num_tx_finished) {
			finished_mean_sd[0] += tx_rates[i];
		}
		else {
			unfinished_mean_sd[0] += tx_rates[i];
		}
		if (final_shut_pos[i] && final_shut_pos[i] > tx_end) {
			preshut_rates[i] = tx_rates[i];
			postshut_rates[i] = -999.99;
			postshut_rates_gene[i] = -999.99;
			postshut_rates_beyond[i] = -999.99;
		}
		else if (final_shut_pos[i]) {
			preshut_rates[i] =
			   ((double) final_shut_pos[i]) /
			   (promo_kill - final_start_time[i]);
			postshut_rates[i] =
			   ((double) (final_pos[i] - final_shut_pos[i])) /
			   (final_stop_time[i] - promo_kill);
			if (final_gene_time[i] > promo_kill) {
				postshut_rates_gene[i] = 
				   ((double) (tsl_stop[1] -
				    final_shut_pos[i]))/
				   (final_gene_time[i] - promo_kill);
				postshut_rates_beyond[i] =
				   ((double) (final_pos[i] -
				    tsl_stop[1])) /
				   (final_stop_time[i]-final_gene_time[i]);
			}
			else if (final_gene_time[i] > 0.01) {
				postshut_rates_gene[i] = -999.99;
				postshut_rates_beyond[i]=postshut_rates[i];
			}
			else {
				postshut_rates_gene[i] =postshut_rates[i];
				postshut_rates_beyond[i] = -999.99;
			}
		}
		else {
			preshut_rates[i] = -999.99;
			postshut_rates[i] = tx_rates[i];
			if (final_gene_time[i] > 0.01) {
				postshut_rates_gene[i] =
				   ((double) tsl_stop[1]) /
				   (final_gene_time[i] -
				    final_start_time[i]);
				postshut_rates_beyond[i] =
				   ((double) (final_pos[i] -
				    tsl_stop[1]))/
				   (final_stop_time[i]-final_gene_time[i]);
			}
			else {
				postshut_rates_gene[i] = tx_rates[i];
				postshut_rates_beyond[i] = -999.99;
			}
		}
	}
			
        FILE *fp_tx;
        sprintf(gen_output,"tx_rate_info_%d.txt",traj_index);
        fp_tx = fopen(gen_output,"w");

	for (i = 1; i <= num_tx; i++) {
		fprintf(fp_tx,"%15.5f%15d%15.5f%15.5f%15.5f"
			"%15.5f%15d%15.5f%15.5f%15.5f\n",
			tx_rates[i],final_shut_pos[i],preshut_rates[i],
			postshut_rates[i],postshut_rates_gene[i],
			postshut_rates_beyond[i],final_pos[i],
			final_start_time[i] - warm_up_time,
			final_gene_time[i] - warm_up_time,
			final_stop_time[i] - warm_up_time);
	}

	mean_sd_tx[0] /= ((double) num_tx);
	for (i = 1; i <= num_tx; i++) {
		mean_sd_tx[1] += ((tx_rates[i] - mean_sd_tx[0]) *
				  (tx_rates[i] - mean_sd_tx[0]));
	}
	mean_sd_tx[1] /= ((double) num_tx);
	mean_sd_tx[1] = sqrt(mean_sd_tx[1]);
	fprintf(fp_tx,"MEAN (SD) TX ELONGATION:%10f %10f\n",mean_sd_tx[0],
		mean_sd_tx[1]);
	
	finished_mean_sd[0] /= ((double) num_tx_finished);
	for (i = 1; i <= num_tx_finished; i++) {
		finished_mean_sd[1] += ((tx_rates[i]-finished_mean_sd[0])*
					(tx_rates[i]-finished_mean_sd[0]));
	}
	finished_mean_sd[1] /= ((double) num_tx_finished);
	finished_mean_sd[1] = sqrt(finished_mean_sd[1]);
	fprintf(fp_tx,"MEAN (SD) TX ELONGATION COMPLETE TX:%10f %10f\n",
		finished_mean_sd[0],finished_mean_sd[1]);

	num_tx_unfinished = num_tx - num_tx_finished;
	unfinished_mean_sd[0] /= ((double) num_tx_unfinished);
	for (i = num_tx_finished + 1; i <= num_tx; i++) {
	    unfinished_mean_sd[1] += ((tx_rates[i]-unfinished_mean_sd[0])*
				      (tx_rates[i]-unfinished_mean_sd[0]));
	}
	unfinished_mean_sd[1] /= ((double) num_tx_unfinished);
	unfinished_mean_sd[1] = sqrt(unfinished_mean_sd[1]);
	fprintf(fp_tx,"MEAN (SD) TX ELONGATION INCOMPLETE TX:%10f %10f\n",
		unfinished_mean_sd[0],unfinished_mean_sd[1]);

	mean_sd_tx[0] = 0;
	mean_sd_tx[1] = 0;
	for (i = 1; i <= num_tx; i++) {
		mean_sd_tx[0] += (1.0/(tx_rates[i]/tx_end));
	}
	mean_sd_tx[0] /= ((double) num_tx);
	for (i = 1; i <= num_tx; i++) {
		mean_sd_tx[1] += (((1.0/(tx_rates[i]/tx_end)) -
				    mean_sd_tx[0]) *
				  ((1.0/(tx_rates[i]/tx_end)) -
				    mean_sd_tx[0]));
	}
	mean_sd_tx[1] /= ((double) num_tx);
	mean_sd_tx[1] = sqrt(mean_sd_tx[1]);
	fprintf(fp_tx,"MEAN (SD) TX ELONGATION:%10f %10f\n",tx_end/mean_sd_tx[0],
		tx_end/mean_sd_tx[1]);
	fclose(fp_tx);

	printf("RADIO PROBES!\n");
        sprintf(gen_output,"probe_based_RNA_abundance_%d.txt",traj_index);
        FILE *fp_probe;
        fp_probe = fopen(gen_output,"w");
	for (i = radio_incr;i <= (int) (sim_stop + 0.5);i += radio_incr) {
		if (i <= (int) (warm_up_time + 0.01)) {
			continue;
		}
		radio_probe[i][0] += radio_probe[i - radio_incr][0];
		radio_probe[i][1] += radio_probe[i - radio_incr][1];
		fprintf(fp_probe,"%12d%12d%12d\n",
			i - ((int) (warm_up_time + 0.01)),
			radio_probe[i][0],radio_probe[i][1]);
	}
	fclose(fp_probe);


	for (i = 1; i <= gene_num; i++) {
		mean_proteins_per_cell[i] /= ((float) num_fish_times);
		printf("MEAN PROTEINS PER CELL GENE %d: %f\n",
			i,mean_proteins_per_cell[i]);
	}

	for (i = 1; i <= num_proteins_made[1]; i++) {
		time = proteins_made[1][i];
		radio_index =
		   ((((int) time)/radio_incr)+1)*radio_incr;
		milleresque[radio_index]++;
	}
	printf("MILLERESQUE!\n");
        sprintf(gen_output,"protein_log_%d.txt",traj_index);
        FILE *fp_mill;
        fp_mill = fopen(gen_output,"w");
	for (i = radio_incr; i <= (int) (sim_stop + 0.5); i += radio_incr) {
		if (i <= (int) (warm_up_time + 0.01)) {
			continue;
		}
		milleresque[i] += milleresque[i - radio_incr];
		fprintf(fp_mill,"%12d%12d",
			i - ((int) (warm_up_time+0.01)), milleresque[i]);
		for (j = 1; j <= 9; j++) {
			fprintf(fp_mill,"%12d",partial_proteins[i][j]);
		}
		fprintf(fp_mill,"\n");
	}
	fclose(fp_mill);

        sprintf(gen_output,"rna_info_log_%d.txt",traj_index);
        FILE *fp_log;
        fp_log = fopen(gen_output,"w");
        for (i = 1; i <= sim_total_rna; i++) {
                fprintf(fp_log,"%10d%10.3f%10.3f%10.3f%10.3f",
                    i,final_rna_log[i].start_tx, final_rna_log[i].end_tx,
                    final_rna_log[i].start_deg,final_rna_log[i].end_deg);
                for (j = 1; j <= num_genes; j++) {
                        fprintf(fp_log,"%10d%10d%10d%10d",
                                final_rna_log[i].ribo_start_nascent[j],
                                final_rna_log[i].ribo_end_nascent[j],
                                final_rna_log[i].ribo_start_mature[j],
                                final_rna_log[i].ribo_end_mature[j]);
                }
                fprintf(fp_log,"\n");
        }
        fclose(fp_log);

	current_reaction = reaction_queue[1].reaction;
	printf("I THINK THE CURRENT RXN IS: %d\n",current_reaction);
	Take_snapshot(curr_time,1,traj,&num_mRNA_snapshot,
                      MAX_RIBO_PER_RNA,RNA_info,snapshot_ribo_mean,
                      rnap_snapshot,current_reaction);
	printf("FINAL SNAPSHOT HAS BEEN TAKEN...\n");
	Write_rna_ribosome_snapshot(num_traj,MAX_RIBO_PER_RNA,
				    RNA_info,num_mRNA_snapshot,
				    traj_index);
	Write_rnap_dna_snapshot(num_promoters,tx_length,
				rnap_snapshot,traj_index);
	Write_pseudoseq_files(tx_length,rnap_seq,ribo_seq,traj_index);
	Write_summary_file(traj_name,start_seed,sim_stop,num_traj,k_on,k_off,k_loading,
                           seq_window,sample_window,snap_time,tx_length,fract,RNA_lifetime,
                           gene_set,protein_mean,protein_sd,snapshot_ribo_mean,
                           snapshot_ribo_sd,mean_ribos_through_rna,sd_ribos_through_rna,
                           fish_final,fish_final_sd,fish_pool,fish_pool_sd,num_genes,
                           traj_index,time_step,k_continue_degrade,p_rnap_rnap_push,
                           p_ribo_rnap_push,p_ribo_ribo_push,pol_width,ribo_width,
                           degrade_width,mean_ribo_all,min_dwell,mean_elongation_rate,
                           coord,stable,translated,read_state,generation_time,b_time,
                           c_time,d_time,k_unloading,k_to_open,k_CSAT,birth_vol,
                           mean_RNAP_conc,mean_ribo_conc,num_transcripts,num_proteins_made,
			   mean_proteins_per_cell,mean_ribo_rate);



	if (REPORT_DWT) {
	 FILE *fp_dwt;
	 sprintf(dwt_name,"dwt_dump_%d.txt",traj_index);
	 fp_dwt = fopen(dwt_name,"w");
	 for (i = 1; i <= num_dwt; i++) {
		fprintf(fp_dwt,"%10.4f\n",dwt[i]);
	 }
	 fclose(fp_dwt);
	}

	if (SINGMOL_TRACE) {
	 FILE *fp_trace;
	 sprintf(trace_name,"single_mol_trace_%d.txt",traj_index);
	 fp_trace = fopen(trace_name,"w");
         for (i = 1; i <= num_trace; i++) {
                fprintf(fp_trace,"%10.3f%10d%10d%10d%10d\n",
                        trace_time[i],trace_pos[i][0],
                        trace_3end[i][0],trace_pos[i][1],
                        trace_3end[i][1]);
         }
	 fclose(fp_trace);
	}
	  
	fclose(fp_for_xtc);

	printf("END RIBOSOME POSITIONS:\n");
	for (i = 1; i <= 4; i++) {
		printf("RIBOSOME %d: %d\n",i,ribosome[i][POSITION]);
	} 
	for (i = 1; i <= gene_num; i++) {
		printf("GENE %d ID: %d\n",i,gene_id[i]);
	} 

        if (record_on) {
                FILE *fp_kymo;
                fp_kymo = fopen("kymo_dump_rotate","w");
                for (i = 1; i <= num_kymo_times; i++) {
                        for (j = 1; j <= num_kymo_loc; j++) {
                                fprintf(fp_kymo,"%10.3f%10d%10.3f%10d%10d%10d",
                                        kymo_times[i],
                                        kymo_loc[j],kymo_data[i][j].sigma,
                                        kymo_data[i][j].RNAP_bound,
                                        kymo_data[i][j].RNAP_state,
                                        kymo_data[i][j].RNAP_rot);
                                if (j == 1) {
                                        fprintf(fp_kymo,"%10d%10d%10d",
                                                first_rnap_engaged[i][0],
                                                first_rnap_engaged[i][1],
                                                first_rnap_engaged[i][2]);
                                }
                                fprintf(fp_kymo,"\n");
                        }
                }
                fclose(fp_kymo);
        }

	for (i = 0; i < MAX_PROMO; i++) {
		free(dna_strip[i]);
		free(rnap_snapshot[i]);
	}
	free(dna_strip);
	free(rnap_snapshot);
	for (i = 0; i < MAX_RNA; i++) {
		free(rna_strip[i]);
		free(rna[i]);
	}
	free(rna);
	free(rna_strip);
	for (i = 0; i < MAX_RIBO; i++) {
		free(ribosome[i]);
	}
	free(ribosome);
	free(mu_rnap_move);
	free(mu_ribo_move);
	free(rna_block_grid);
	free(active_riboseq);
	free(active_rnapseq);
	free(rnap_seq);
	free(ribo_seq);
	free(ribo_temp);
	free(overlap_zone);
	free(adj_factor);

	return(43);
}
	











