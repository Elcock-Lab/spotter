#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define DEFINE_UTIL_VARIABLES
#include "INCL/utilities.C_RNG.h"

#include "INCL/sim_types.C_RNG.h"
#include "INCL/gen_def.C_RNG.h"
#include "INCL/sim_settings.C_RNG.h"
#include "INCL/io.C_RNG.h"
#include "INCL/reaction_manager.C_RNG.h"
#include "INCL/replication.C_RNG.h"
#include "INCL/regulation.C_RNG.h"
#include "INCL/transcription.C_RNG.h"
#include "INCL/translation.C_RNG.h"
#include "INCL/supercoiling.C_RNG.h"
#include "INCL/topo.C_RNG.h"


void Add_kymograph_info(double time,int current_reaction, int kymo_index)
{
        int i,check,which_promo,site;
        int barrier_check,up_check,down_check;
        int temp_plus=0,temp_minus=0;
        int spacing = 0,test_len;
        int first_rnap,last_rnap,first_tx_rnap=0;
        int min_spacing = 0;
        float test_lk;

	which_promo = 1;
	first_tx_rnap = promoter[which_promo][NEWEST_RNAP];
	if (rnap[first_tx_rnap][STATE] == CLOSED) {
        	first_tx_rnap = rnap[first_tx_rnap][PLUS_ONE];
	}
	if (first_tx_rnap) {
        	first_rnap_engaged[kymo_index][0] = first_tx_rnap;
        	first_rnap_engaged[kymo_index][1] =
					rnap[first_tx_rnap][POSITION];
        	first_rnap_engaged[kymo_index][2] =
					curr_rnap_rot[first_tx_rnap];
	}
	first_rnap = first_tx_rnap;
	last_rnap = promoter[which_promo][OLDEST_RNAP] *
	    (rnap[promoter[which_promo][NEWEST_RNAP]][STATE] != CLOSED);

	for (i = 1; i <= num_kymo_loc; i++) {
            site = kymo_loc[i];

            if (site < 1) {
                if (promoter[which_promo][NEWEST_RNAP] &&
                    rnap[promoter[which_promo][NEWEST_RNAP]][STATE] !=
                    CLOSED) {
                        temp_plus =
                           promoter[which_promo][NEWEST_RNAP];
                }
                else {                                                          // This should give the proper
                        temp_plus =                                             // right boundary whether or
                           rnap[promoter[which_promo][NEWEST_RNAP]]             // not an RNAP is present
                               [PLUS_ONE];
                }
            }
            else if (site > tx_end) {
                if (promoter[which_promo][OLDEST_RNAP] &&
                    rnap[promoter[which_promo][OLDEST_RNAP]][STATE] !=
                    CLOSED) {
                        temp_minus =
                           promoter[which_promo][OLDEST_RNAP];
                }
                else {
                        temp_minus =
                           rnap[promoter[which_promo][OLDEST_RNAP]]
                               [MINUS_ONE];
                }
            }
            else {
                check = promoter[which_promo][NEWEST_RNAP];
                while (check && rnap[check][POSITION] <= site) {
                 if (rnap[check][STATE] != CLOSED) {
                  temp_minus = check;
                  kymo_data[kymo_index][rnap[check][POSITION]].RNAP_bound=      // At 1 bp resoultion, site and
                        check;                                                  // real position are the same;
                  kymo_data[kymo_index][rnap[check][POSITION]].RNAP_state=      // will be redundant but OK
                        rnap[check][STATE];
                  kymo_data[kymo_index][rnap[check][POSITION]].RNAP_rot=
                        curr_rnap_rot[check];
                 }
                        check = rnap[check][PLUS_ONE];
                }
                if (temp_minus) {
                        temp_plus =
                           rnap[temp_minus][PLUS_ONE];
                }
                else {
                        if (!(promoter[which_promo][NEWEST_RNAP]) ||
                            rnap[promoter[which_promo][NEWEST_RNAP]][STATE]
                            == CLOSED) {
                                temp_plus =
                                  rnap[promoter[which_promo][NEWEST_RNAP]]
                                      [PLUS_ONE];
                        }
                        else {
                                temp_plus =
                                   promoter[which_promo][NEWEST_RNAP];
                        }
                }
            }

            if (temp_minus && temp_plus) {
                barrier_check =
                   rnap[temp_minus][BARRIER_AHEAD];
                up_check =
                   rnap[temp_minus][POSITION];
                down_check =
                   rnap[temp_plus][POSITION];
                test_lk = lk[temp_minus][DOWNSTREAM];
                test_len = rnap_sep[temp_minus][DOWNSTREAM];
		if (barrier_check <= 0) {
                   if (inter_rnap_topo) {
                        spacing = ((site - up_check >= min_spacing) &&          // Use if interRNAP OK
                                   (down_check - site >= min_spacing));
                   }
                   else {
                        spacing = 0;	                                        // If activity between RNAPs
                   }								// not allowed, topo cannot
		}                                                       	// be placed
             }
             else if (temp_minus) {

                barrier_check =                                                             
                   rnap[temp_minus][BARRIER_AHEAD];                               
                up_check =
                   rnap[temp_minus][POSITION];
                test_lk = lk[temp_minus][DOWNSTREAM];
                test_len = rnap_sep[temp_minus][DOWNSTREAM];
		if (barrier_check <= 0) {
                   if (!(plasmid)) {
                                spacing = ((site -
                                    rnap[temp_minus][POSITION]
                                    >= min_spacing) &&
                                   (stop_barrier - site >=
                                    min_spacing));
                   }
                   else {
                                spacing = (site-rnap[temp_minus][POSITION]      // No problem with stop
                                           >= min_spacing);                     // barrier in plasmid
                   }
		}
            }
            else if (temp_plus) {

                barrier_check =
                   rnap[temp_plus][BARRIER_BEHIND];
                down_check =
                   rnap[temp_plus][POSITION];
                test_lk = lk[temp_plus][UPSTREAM];
                test_len = rnap_sep[temp_plus][UPSTREAM];
		if (barrier_check <= 0) {
                   if (!(plasmid)) {
                                spacing =
                                   ((rnap[temp_plus][POSITION]
                                   - site >= min_spacing) &&
                                   (site - start_barrier >= min_spacing));
                   }
                   else {
                                spacing = (rnap[temp_plus][POSITION] -site      // No problem with start
                                           >= min_spacing);                     // barrier in plasmid
                  }
	       }
            }
            else {
               if (!(plasmid)) {
                                spacing =
                                  ((stop_barrier - site >= min_spacing) &&
                                  (site - start_barrier >= min_spacing));
               }
               else {                                                  		// Good to go in plasmid:
                                spacing = 1;                                    // no need to check start
               }                                                       		// or stop barriers
               test_lk = open_lk[which_promo];
               test_len = stop_barrier - start_barrier;
            }

            kymo_data[kymo_index][i].sigma =
               (test_lk - relaxed_lk[test_len])/relaxed_lk[test_len];

	}

        Update_reaction_queue(current_reaction,INFINITY);
        Remove_reaction(current_reaction);

        if (kymo_index == num_kymo_times) {
                sim_endtime = time + 0.0001;
        }
	spacing++;
	first_rnap--;
	last_rnap--;
        return;
}

void Add_FISH_data(double time,int check_stage,int *fish_probe,
                   int fish_total[][MAX_FISH_PROBES],int current_reaction)
{
	int i;
	printf("WORKING ON FISH AT TIME %f...\n",time);
	for (i = 1; i < MAX_RNA; i++) {
		if (rna[i][SOURCE_PROMO] == 0) {
			continue;
		}
		if (rna[i][FIVE_END] == fish_probe[1]) {
			fish_total[check_stage][1]++;
		}
		if (rna[i][FIVE_END] <= fish_probe[2] && rna[i][THREE_END] >=
                    fish_probe[2]) {
			fish_total[check_stage][2]++;
		}
		if (rna[i][THREE_END] == tx_end) {
			fish_total[check_stage][3]++;
		}			
	}
	for (i = 1; i <= gene_num; i++) {
		mean_proteins_per_cell[i] += ((float) num_proteins_made[i]);
	}
	Update_reaction_queue(current_reaction,INFINITY);
	Remove_reaction(current_reaction);
	return;
}
	

void Add_pseudoseq_data(double time,int check_stage,int *rnap_seq,int *ribo_seq,
                        int promo_index,int current_reaction, FILE *fp)
{
	int i,j;
	for (i = 1; i < MAX_RIBO; i++) {
		if (ribosome[i][1] != 0) {
			j = Check_for_ribosome_stacking(i);
			if (j == 0) {
				active_riboseq[ribosome[i][POSITION]]++;
			}
			else {
//				printf("NO DICE: RIBOSOME WAS STACKED!\n");
			}
		}
	}
	for (i = 1; i <= tx_end; i++) {
		ribo_seq[i] += active_riboseq[i];
//		for (j = 1; j <= num_promoters; j++) {
//			rnap_seq[i] += (dna_strip[j][i] > 0);
//		}
	}
        for (i = 1; i < MAX_RNAP; i++) {
                if (rnap[i][POSITION]) {
                        rnap_seq[rna[rnap[i][RNAP_RNA]][THREE_END]]++;
                }
        }
        if (SINGMOL_TRACE) {
                num_trace++;
                trace_time[num_trace] = time;
                trace_pos[num_trace][0] = rnap[trace_RNAP][POSITION];
                trace_3end[num_trace][0] =
                   rna[rnap[trace_RNAP][RNAP_RNA]][THREE_END];
                trace_pos[num_trace][1] = rnap[trace_RNAP+1][POSITION];
                trace_3end[num_trace][1] =
                   rna[rnap[trace_RNAP+1][RNAP_RNA]][THREE_END];
        }
	memset(active_riboseq,0,(tx_size+1)*sizeof(*active_riboseq));
	Update_reaction_queue(current_reaction,INFINITY);
	Remove_reaction(current_reaction);

	frame++;
	Write_frame(time,promo_index,fp);

	return;
}
	
void Log_nascent_proteins(double time, int current_reaction)
{
	int i,tenth,time_index;
	int frac_done;
	time_index = (int) (time + 0.01);
	tenth = (tsl_stop[1] - tsl_start[1]) + 1;
	tenth /= 10;

	printf("TAKING CENSUS OF PARTIAL PROTEINS AT TIME %f; "
		"TENTH IS %d\n",
		time,tenth);

	for (i = 1; i < MAX_RIBO; i++) {
		if (ribosome[i][1] != 0) {
			frac_done = (ribosome[i][POSITION]-tsl_start[1]) /
				     tenth;
			partial_proteins[time_index][frac_done]++;
		}
	}							
	Update_reaction_queue(current_reaction,INFINITY);
	Remove_reaction(current_reaction);
	return;
}

void Determine_translation_efficiencies(int rna_log_ctr,
					int **rna_log,
					int num_genes,
                                        float *mean_ribos_through_rna,
					float *sample_window)
{
	int i,j,ctr=0,a_quo,ad_quem;
	a_quo = (int) sample_window[0];
	ad_quem = (int) sample_window[1];
	for (i = 1; i <= rna_log_ctr; i++) {
		if (rna_log[i][0] >= a_quo && rna_log[i][0] <= ad_quem) {
			ctr++;
			for (j = 1; j <= num_genes; j++) {
				mean_ribos_through_rna[j] +=
					((float ) rna_log[i][j]);
			}
		}
	}
	for (i = 1; i <= num_genes; i++) {
		mean_ribos_through_rna[i] /= ((float ) ctr);
	}
	return;
}
		
void Determine_protein_production(int *num_proteins_made,
				  int MAX_PROTEINS,
                                  float **proteins_made,
				  float *sample_window,
                                  int num_ten_min_windows,
				  int num_genes,
				  float *protein_mean)
{
	int i,j,k;
	float  ii;
	for (i = 0; i < num_ten_min_windows; i++) {
		ii = (float ) i;
		for (j = 1; j <= num_genes; j++) {
			for (k = 1; k <= num_proteins_made[j]; k++) {
			    if (proteins_made[j][k] >= 
				   sample_window[0] + (ii * 600) &&
                                proteins_made[j][k] <=
				   sample_window[0] + ((ii + 1) * 600)) {
					protein_mean[j] += 1.0;
				}
			}
		}
	}
	for (i = 1; i <= num_genes; i++) {
		protein_mean[i] /= ((float ) num_ten_min_windows);
		protein_mean[i] *= 6.0;
	}
	return;
}

void Take_snapshot(double time, int promoter_index, int traj,
                   int *num_mRNA_snapshot, int MAX_RIBO_PER_RNA,
                   int **RNA_info, float *snapshot_ribo_mean,
                   int **rnap_snapshot, int current_reaction)
{
	int i,j,k,ctr=0;
	for (i = 1; i <= tx_end; i++) {
		for (j = 1; j <= num_promoters; j++) {
			rnap_snapshot[j][i] = dna_strip[j][i];
		}
	}
	for (i = 1; i < MAX_RNA; i++) {
		if (rna[i][1] == 0) {
			continue;
		}
		*num_mRNA_snapshot += 1;
		RNA_info[*num_mRNA_snapshot][0] = i;			
		RNA_info[*num_mRNA_snapshot][1] = rna[i][FIVE_END];
		RNA_info[*num_mRNA_snapshot][2] = rna[i][THREE_END];
		RNA_info[*num_mRNA_snapshot][3] = (rna[i][MATURE]  == 0);
		ctr = 0;
		for (j = 1; j <= tx_end; j++) {
			if (rna_strip[i][j] > 0) {
				ctr++;
				RNA_info[*num_mRNA_snapshot][4 + ctr] = j;
			}
		}
		for (j = 1; j <= gene_num; j++) {
			for (k = tsl_start[j]; k <= tsl_stop[j]; k++) {
				if (rna_strip[i][k] > 0) {
					snapshot_ribo_mean[j] += 1.0;
				}
			}
		}
		RNA_info[*num_mRNA_snapshot][4] = ctr;
	}
	printf("TAKING SNAPSHOT AT %f WITH %d RNAs\n",
		time,*num_mRNA_snapshot);
	for (i = 1; i <= gene_num; i++) {
		snapshot_ribo_mean[i] /= ((float ) *num_mRNA_snapshot);
	}
	Update_reaction_queue(current_reaction,INFINITY);
	Remove_reaction(current_reaction);
	return;
}
		
