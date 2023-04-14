#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define DEFINE_TSL_VARIABLES
#include "INCL/translation.C_RNG.h"

#include "INCL/sim_types.C_RNG.h"
#include "INCL/gen_def.C_RNG.h"
#include "INCL/sim_settings.C_RNG.h"
#include "INCL/reaction_manager.C_RNG.h"
#include "INCL/replication.C_RNG.h"
#include "INCL/transcription.C_RNG.h"
#include "INCL/supercoiling.C_RNG.h"
#include "INCL/utilities.C_RNG.h"


int Determine_RNAP_riboload_effect(double time, int rnap_index,
				   int orig_three_end, int just_made)
{
	int i,start,rna_index;
	double next_fire;
	rna_index = rnap[rnap_index][RNAP_RNA];
	if (translated == ON || translated == HYBRID) {
		if (translated == HYBRID) {
			start = gene_num;
		}
		else {
			start = 1;
		}
		for (i = start; i <= gene_num; i++) {
			if (tsl_start[i] > orig_three_end -
			    pol_ribo_width) {					// Don't waste time on
				break;						// untranscribed binding
			}							// sites...
			if (rna[rna_index][TSL_ON_OFF[i]] == KILLED) {		
				continue;					// ...or degraded sites
			}
			if ((just_made || rnap[rnap_index][STATE] ==		// This covers transcripts
			    PRE_TRANSLOC) && rna[rna_index][THREE_END] ==	// newly extended to open
			    tsl_start[i] + pol_ribo_width) {			// a ribosome binding site;
				rna[rna_index][TSL_ON_OFF[i]] = ON;		// fn is only accessed with
				next_fire = time + Find_firing_time(		// this state on nt add
				   mu_ribo_load[i]);
				Add_reaction_to_queue(RIBO_LOAD_RXN(
				   rna_index,i),next_fire);
				break;
			}
			if (rnap[rnap_index][STATE] == OFF_P2 ||		// This covers the effects
			    rnap[rnap_index][STATE] == OFF_P3) {		// of off-pathway movement:
				if (rnap[rnap_index][POSITION] >=
				    tsl_start[i] + pol_ribo_width) {
				     if (rna[rna_index][TSL_ON_OFF[i]] ==	// If movement without effect,
				         ON) {					// leave loading on...
					  continue;
				     }
				     else {					// If it frees up the ribo
					  rna[rna_index][TSL_ON_OFF[i]] =	// site, add back loading 
					    ON;					// rxns. Note that this will
					  next_fire = time +			// never generate a new site
						      Find_firing_time(		// b/c of ban on forward
							mu_ribo_load[i]);	// tracking
					  Add_reaction_to_queue(
					     RIBO_LOAD_RXN(rna_index,i),
					        next_fire);
				     }
				}
				else {
				     if (rna[rna_index][TSL_ON_OFF[i]] ==	// If backtracking does not
				         OFF) {					// affect site, no change...
					  continue;
				     }
				     else {					// If it results in the
					  rna[rna_index][TSL_ON_OFF[i]] =	// occlusion of the site,
					    OFF;				// remove loading rxns
					  Update_reaction_queue(
					     RIBO_LOAD_RXN(rna_index,i),
					        INFINITY);
					  Remove_reaction(
					     RIBO_LOAD_RXN(rna_index,i));
				     }
				}
			}
		}
	}
	return(1);
}

int Check_for_ribosome_stacking(int ribo_index)
{
	int i,blocker=0;
	int rna_index,start,stop;
	rna_index = ribosome[ribo_index][SOURCE_RNA];
	for (i = (ribosome[ribo_index][POSITION] + ribo_width) - 1; i <=
             ribosome[ribo_index][POSITION] + ribo_width + 9; i++) {
		if (rna_strip[rna_index][i] > 0) {
			blocker = rna_strip[rna_index][i];			
			break;							
		}
	}
	if (blocker > 0) {
		return(blocker);
	}
	start = MAX_VAL(0,(ribosome[ribo_index][POSITION] -
			   ribo_width - 9));
	stop = MAX_VAL(0,((ribosome[ribo_index][POSITION] - 
			   ribo_width) + 1));
	for (i = start; i <= stop; i++) {
		if (rna_strip[rna_index][i] > 0) {
			blocker = rna_strip[rna_index][i];			
			break;							
		}
	}
	return(blocker);
}

int Check_for_downstream_ribosome(int ribo_index)
{
	int i,blocker=0;
	int rna_index;
	rna_index = ribosome[ribo_index][SOURCE_RNA];
	for (i = (ribosome[ribo_index][POSITION] + ribo_width) - 1; i <=
             ribosome[ribo_index][POSITION] + ribo_width + 3; i++) {
		if (rna_strip[rna_index][i] > 0) {
			blocker = rna_strip[rna_index][i];			// Should only need to check these because
			break;							// of loading restrictions
		}
	}
	return(blocker);
}

int Add_collision_stimulated_termination_event(double time, int ribo_index)
{
	double next_fire;
	ribosome[ribo_index][ABORT_POSSIBLE] = ON;
        if (!(reaction_index[CSAT_TERM_RXN(ribo_index)] < 0)) {
                return(1);
        }
	next_fire = time + Find_firing_time(mu_ribo_csat);
	Add_reaction_to_queue(CSAT_TERM_RXN(ribo_index),next_fire);
	return(1);
}

int Load_ribosome(double time, int gene_index, int rna_index)
{
	int rnap_index,result,blocker=0;
	double next_fire;
	next_fire = time + Find_firing_time(mu_ribo_load[gene_index]);
	Update_reaction_queue(RIBO_LOAD_RXN(rna_index,gene_index),
	   next_fire);
	if (rna[rna_index][TSL_ON_OFF[gene_index]] < ON) {			// To cover "KILLED" status from degr
		return(0);
	}
	if (rna[rna_index][TSL_BLOCK[gene_index]] > 0) {			// To cover occlusion of RBS
		return(0);
	}
	ribosome[first_free_ribo][SOURCE_RNA] = rna_index;
	ribosome[first_free_ribo][SOURCE_GENE] = gene_index;
	ribosome[first_free_ribo][POSITION] = tsl_start[gene_index];
	ribosome[first_free_ribo][MINUS_ONE] = 0;

	sim_total_ribosomes++;
	ribosome[first_free_ribo][SUPER_INDEX] = sim_total_ribosomes;

        if (!(rna[rna_index][MATURE])) {
          final_rna_log[
           rna[rna_index][MASTER_INDEX]].ribo_start_nascent[gene_index]++;
        }
        else {
          final_rna_log[
           rna[rna_index][MASTER_INDEX]].ribo_start_mature[gene_index]++;
        }

	ribosome[first_free_ribo][ARRIVAL] = (int) (time/0.001);		// For calculating elongation rate:
										// keeps as int; ms precision should
	blocker = Check_for_downstream_ribosome(first_free_ribo);		// be adequate
	rnap_index = rna[rna_index][SOURCE_RNAP];
	if (rna[rna_index][NEWEST_RIBO[gene_index]] != 0) {
	      ribosome[first_free_ribo][PLUS_ONE] =
		   rna[rna_index][NEWEST_RIBO[gene_index]];
	      ribosome[rna[rna_index][NEWEST_RIBO[gene_index]]][MINUS_ONE]
		    = first_free_ribo;
	}
	else {
		ribosome[first_free_ribo][PLUS_ONE] = 0;
		rna[rna_index][OLDEST_RIBO[gene_index]] = first_free_ribo;

		if (rna[rna_index][OLDEST_RIBO[gene_index - 1]] == 0 ||		// This is to deal with the case
		    ribosome[rna[rna_index][OLDEST_RIBO[gene_index - 1]]]	// where genes are so overlapped
		      [POSITION] < ribosome[first_free_ribo][POSITION]) {	// that a ribosome can initiate
			rnap[rnap_index][CLOSEST_RIBO] = first_free_ribo;	// on an downstream gene *behind*

			if (!(QUIET_MODE)) {
			  printf("RIBOSOME %d LOADED; NEAREST RNAP IS %d "
				 "AT %d\n",
			  first_free_ribo,rnap_index,
			  rnap[rnap_index][POSITION]);
			}

			if (Check_ribosome_tether(rnap_index,0)) {		// ribosomes on the upstream gene
				rnap[rnap_index][RIBO_TETHER] =
				   first_free_ribo;

			    if (!(QUIET_MODE)) {
				printf("LOADED RIBOSOME %d FORMS COMPLEX "
					"WITH RNAP %d\n", 
					first_free_ribo,rnap_index);
			    }

			}
			else {

			    if (!(QUIET_MODE)) {
				printf("LOADED RIBOSOME %d DID NOT FORM "
					"COMPLEX WITH RNAP %d\n",
					first_free_ribo,rnap_index);
			    }

			}

		}
	}
	rna[rna_index][TSL_BLOCK[gene_index]]++;
	rna[rna_index][NUM_RIBO[gene_index]]++;
	rna[rna_index][TOTAL_RIBO]++;

	if (rnap_relax && ribo_resist_model == EXPLICIT_RIBOSOMES &&
	    !(rna[rna_index][MATURE])) {
		Update_RNAP_complex_resistance(rnap_index);
		result = Update_relaxation_rates(time,rnap_index);		// Adjust rotation rates for
										// new ribosome (not nec. for
	}									// RNA nt) b/c covered on 
										// forward translocation

	rna[rna_index][NEWEST_RIBO[gene_index]] = first_free_ribo;
	ribosome[first_free_ribo][STATE] = POST_TRANSLOC;
	next_fire = time + Find_firing_time(
		  ribo_dwell[ribosome[first_free_ribo][POSITION]].mu_add);
	Add_reaction_to_queue(RIBO_DECODE_RXN(first_free_ribo),next_fire);
	rna_strip[rna_index][tsl_start[gene_index]] = first_free_ribo;
	ribo_bound++;


	if (rna[rna_index][MATURE] == 0 && rnap[rnap_index][STATE]
	    == OFF_P1 && RIBO_ANTIPAUSE) {
		result = Update_RNAP_rates(time,rnap_index);			// Now accounting for torque
	}
			
	if (rna[rna_index][MATURE] == 0 &&
	    ((rnap[rnap_index][STATE] == OFF_P2 ||				// Trying ribosome
	      rnap[rnap_index][STATE] == OFF_P3) &&				// accelerant...
	     rnap[rnap_index][POSITION] ==
	     rna[rna_index][THREE_END]) && RIBO_ANTIPAUSE) {
		result = Update_RNAP_rates(time,rnap_index);
	}
	first_free_ribo = Find_first_free(ribosome);
	blocker++;
	result++;
	return(1);
}

int Advance_ribosome(double time, int ribo_index, int gene_index,
		     int rna_index, int *num_proteins_made,
		     int MAX_PROTEINS, float **proteins_made,
		     int *num_codons, float *mean_elongation_rate)
{
	int block_value;
	int blocker=0,result;
	int rnap_index;
	int same_gene_ribo,next_gene_ribo;
	float  elong_rate;
	double next_fire;
	block_value = rna_block_grid[ribosome[ribo_index][POSITION]];		// Determine whether original
	rna_strip[rna_index][ribosome[ribo_index][POSITION]] = 0;		// and new positions occlude
	rna[rna_index][TSL_BLOCK[block_value]] -= (block_value > 0);  		// a translation start using
	ribosome[ribo_index][POSITION] += 3;					// master block grid
	block_value = rna_block_grid[ribosome[ribo_index][POSITION]];
	rnap_index = rna[rna_index][SOURCE_RNAP];

        if (!(reaction_index[CSAT_TERM_RXN(ribo_index)] < 0)) {
                Update_reaction_queue(CSAT_TERM_RXN(ribo_index),INFINITY);
                Remove_reaction(CSAT_TERM_RXN(ribo_index));
        }

        if (ribosome[ribo_index][POSITION] > tsl_stop[gene_index] ||
            (rna[rna_index][MATURE] &&                                          // Covering possibility of
             ribosome[ribo_index][POSITION] > rna[rna_index][THREE_END])) {     // ribosome moving past 3'
                ribosome[ribosome[ribo_index][MINUS_ONE]][PLUS_ONE] = 0;        // end of terminated RNA:
                rna[rna_index][OLDEST_RIBO[gene_index]] =                       // NOTE currently moving
                   ribosome[ribo_index][MINUS_ONE];                             // smoothly off RNA end
		if (rna[rna_index][NEWEST_RIBO[gene_index]] ==
		    ribo_index) {
			rna[rna_index][NEWEST_RIBO[gene_index]] = 0;				
		}
		if (!(rna[rna_index][MATURE]) &&
		    ribo_index == rnap[rnap_index][CLOSEST_RIBO]) {		// Reassign closest ribo if
			same_gene_ribo = ribosome[ribo_index][MINUS_ONE];	// necessary; deal with
			next_gene_ribo = 					// possibility of overlap;
			   rna[rna_index][OLDEST_RIBO[gene_index + 1]];		// if neither of the checked
			if (rnap[rnap_index][POSITION] -			// ribosomes exists, CLOSEST
			    ribosome[same_gene_ribo][POSITION] <		// will automatically be 0,
			    rnap[rnap_index][POSITION] -			// obviating the tethering
			    ribosome[next_gene_ribo][POSITION]) {		// check below
				rnap[rnap_index][CLOSEST_RIBO] =
				   same_gene_ribo;
			}
			else {
				rnap[rnap_index][CLOSEST_RIBO] =
				   next_gene_ribo;
			}
			if (rnap[rnap_index][CLOSEST_RIBO] &&			
			    rnap[rnap_index][RIBO_TETHER]) {			// If the RNAP isn't already
				rnap[rnap_index][RIBO_TETHER] = 0;		// tethered, it can't be
				if (Check_ribosome_tether(rnap_index,0)) {	// yoked to a more distant
					rnap[rnap_index][RIBO_TETHER] =		// ribosome; covers
					   rnap[rnap_index][CLOSEST_RIBO];	// *extreme* interxn dist's
				}
			}

			if (!(QUIET_MODE)) {
				printf("RIBOSOME %d FINISHES GENE; "
					"CLOSEST RIBOSOME IS NOW %d AND "
					"TETHER IS NOW %d\n", 
					ribo_index, 
					rnap[rnap_index][CLOSEST_RIBO],
					rnap[rnap_index][RIBO_TETHER]);
			}

		}
			

		Update_reaction_queue(RIBO_MOVE_RXN(ribo_index),INFINITY);
		Remove_reaction(RIBO_MOVE_RXN(ribo_index));

                if (ribosome[ribo_index][POSITION] > tsl_stop[gene_index]) {
                  elong_rate = (float ) num_codons[gene_index];
                  elong_rate /= (time - (0.001 *
                                ((float ) ribosome[ribo_index][ARRIVAL])));
                  mean_elongation_rate[gene_index] += elong_rate;
                  num_proteins_made[gene_index]++;
                  product_added[gene_index]++;
                  proteins_made[gene_index][num_proteins_made[gene_index]]
                     = time;
		  rna[rna_index][TOTAL_RIBO]--;
		  if (!(rna[rna_index][MATURE]) && rnap_relax &&
		      ribo_resist_model == EXPLICIT_RIBOSOMES) {
		     Update_RNAP_complex_resistance(rnap_index);
		     Update_relaxation_rates(time,rnap_index);
		  }
                  if (!(rna[rna_index][MATURE])) {
                     final_rna_log[
                      rna[rna_index][
                        MASTER_INDEX]].ribo_end_nascent[gene_index]++;
                  }
                  else {
                     final_rna_log[
                      rna[rna_index][
                        MASTER_INDEX]].ribo_end_mature[gene_index]++;
                  }
		     
                }

		memset(ribosome[ribo_index],0,
		   RIBO_INFO*sizeof(*ribosome[ribo_index]));
		ribo_bound--;
		first_free_ribo = ribo_index;
	}
	else {
		rna[rna_index][TSL_BLOCK[block_value]] += (block_value
		   > 0);  
		rna_strip[rna_index][ribosome[ribo_index][POSITION]]
		   = ribo_index;
		ribosome[ribo_index][ABORT_POSSIBLE] = OFF;

		ribosome[ribo_index][STATE] = POST_TRANSLOC;
		Update_reaction_queue(RIBO_MOVE_RXN(ribo_index),INFINITY);
		Remove_reaction(RIBO_MOVE_RXN(ribo_index));
		next_fire = time + Find_firing_time(
		   ribo_dwell[ribosome[ribo_index][POSITION]].mu_add);
		Add_reaction_to_queue(RIBO_DECODE_RXN(ribo_index),
		   next_fire);


		blocker = Check_for_downstream_ribosome(ribo_index);

		if (!(rna[rna_index][MATURE]) && ribo_index ==
		    rnap[rnap_index][CLOSEST_RIBO]) {
			if (!(rnap[rnap_index][RIBO_TETHER])) {
				if (Check_ribosome_tether(rnap_index,0)) {
					rnap[rnap_index][RIBO_TETHER] =
					   ribo_index;

				   if (!(QUIET_MODE)) {
					printf("RIBOSOME %d MOVED TO %d; "
						"COMPLEX WITH RNAP %d AT "
						"%d FORMED\n", 
						ribo_index,
						ribosome[
						  ribo_index][POSITION],
						rnap_index,rnap[
						  rnap_index][POSITION]);
				   }

				}
			}
			if (RIBO_ANTIPAUSE && rnap[rnap_index][STATE] ==
			    OFF_P1) {
				result = Update_RNAP_rates(time,		// Includes torque check
					    rnap_index);
			}
			if (RIBO_ANTIPAUSE && ((rnap[rnap_index][STATE] ==
			    OFF_P2 || rnap[rnap_index][STATE] == OFF_P3)
			    && rnap[rnap_index][POSITION] ==
			    rna[rna_index][THREE_END])) {
				result = Update_RNAP_rates(time,
					    rnap_index);
			}
		}
	}
	blocker++;
	result++;
	return(1);
}

int Add_amino_acid(double time, int ribo_index)
{
	double next_fire;
	ribosome[ribo_index][STATE] = PRE_TRANSLOC;
	Update_reaction_queue(RIBO_DECODE_RXN(ribo_index),INFINITY);
	Remove_reaction(RIBO_DECODE_RXN(ribo_index));
	next_fire = time +
		    ribo_dwell[ribosome[ribo_index][POSITION]].mu_f;
	Add_reaction_to_queue(RIBO_MOVE_RXN(ribo_index),next_fire);
	return(1);
}

int Remove_stalled_ribosome(double time, int ribo_index)
{
	int gene_index, rna_index,block_value;
	gene_index = ribosome[ribo_index][SOURCE_GENE];				// Account for loss of ribosome
	rna_index = ribosome[ribo_index][SOURCE_RNA];				// in RNA information files and
	block_value = rna_block_grid[ribosome[ribo_index][POSITION]];		// eliminate gaps in +1/-1
	rna[rna_index][TSL_BLOCK[block_value]] -= (block_value > 0);  		// ribosome slots, then recycle
	rna_strip[rna_index][ribosome[ribo_index][POSITION]] = 0;
	if (ribo_index == rna[rna_index][OLDEST_RIBO[gene_index]]) {
		rna[rna_index][OLDEST_RIBO[gene_index]] =
		   ribosome[ribo_index][MINUS_ONE];	
	} 											
	if (ribo_index == rna[rna_index][NEWEST_RIBO[gene_index]]) {				
		rna[rna_index][NEWEST_RIBO[gene_index]] =
		   ribosome[ribo_index][PLUS_ONE];	
	}
	ribosome[ribosome[ribo_index][PLUS_ONE]][MINUS_ONE] =
	   ribosome[ribo_index][MINUS_ONE]; 
	ribosome[ribosome[ribo_index][MINUS_ONE]][PLUS_ONE] =
	   ribosome[ribo_index][PLUS_ONE];
	if (!(reaction_index[RIBO_MOVE_RXN(ribo_index)] < 0)) {
		Update_reaction_queue(RIBO_MOVE_RXN(ribo_index),INFINITY);
		Remove_reaction(RIBO_MOVE_RXN(ribo_index));
	}
	if (!(reaction_index[RIBO_DECODE_RXN(ribo_index)] < 0)) {
		Update_reaction_queue(RIBO_DECODE_RXN(ribo_index),
		   INFINITY);
		Remove_reaction(RIBO_DECODE_RXN(ribo_index));
	}
	if (!(reaction_index[CSAT_TERM_RXN(ribo_index)] < 0)) {
		Update_reaction_queue(CSAT_TERM_RXN(ribo_index),INFINITY);
		Remove_reaction(CSAT_TERM_RXN(ribo_index));
	}
	memset(ribosome[ribo_index],0,
	       RIBO_INFO*sizeof(*ribosome[ribo_index]));
	ribo_bound--;
	first_free_ribo = ribo_index;
	return(1);
}

int Reschedule_ribosome_translocation(double time, int ribo_index)
{
	double next_fire;
	next_fire = time + Find_firing_time(
	   ribo_dwell[ribosome[ribo_index][POSITION]].mu_f);
	Update_reaction_queue(RIBO_MOVE_RXN(ribo_index),next_fire);
	return(1);
}

int Attempt_ribosome_move(double time, int ribo_index,
			  float p_push_ribo_forward,
			  float p_ribo_rnap_push, float p_rnap_rnap_push,
                          int *num_proteins_made, int MAX_PROTEINS,
                          float **proteins_made, int *num_codons,
			  float *mean_elongation_rate)
{
	int gene_index, rna_index,rnap_index;
	int ok_dist,result;
	int blocker = 0, move_open=0;
	float  push;
	gene_index = ribosome[ribo_index][SOURCE_GENE];
	rna_index = ribosome[ribo_index][SOURCE_RNA];
	blocker = Check_for_downstream_ribosome(ribo_index);
	if (blocker != 0) {
                if (mu_ribo_csat > 0) {
                 Add_collision_stimulated_termination_event(time,blocker);
                }
		push = Sample_uniform_distribution();
		if (push < p_push_ribo_forward) {
			move_open = 0;
			result = Attempt_ribosome_move(
				   time,blocker,p_push_ribo_forward,
                                   p_ribo_rnap_push,p_rnap_rnap_push,
                                   num_proteins_made,MAX_PROTEINS,
                                   proteins_made,num_codons,
                                   mean_elongation_rate);
			move_open = Check_for_downstream_ribosome(
				       ribo_index);
			if (move_open == 0) {
				result = Advance_ribosome(
					 time,ribo_index,gene_index,
					 rna_index,num_proteins_made,
					 MAX_PROTEINS,proteins_made,
					 num_codons,mean_elongation_rate);
			}
			else {
				Reschedule_ribosome_translocation(
				   time,ribo_index);
			}
		}
		else {
			Reschedule_ribosome_translocation(time,
			   ribo_index);
		}
	}	
	else {
		rnap_index = rna[rna_index][SOURCE_RNAP];
		ok_dist = 3 - RIBO_RNAP_RANGE;					// Adjustment to distance
		if (rna[rna_index][MATURE] == 1 ||				// allows testing FX
                    (abs((ribosome[ribo_index][POSITION] + ok_dist) -		// of ribosome
		    rnap[rnap_index][POSITION])) >= pol_ribo_width) {		// "springiness"
			result = Advance_ribosome(time,ribo_index,
				   gene_index,rna_index,num_proteins_made,
				   MAX_PROTEINS,proteins_made,num_codons,
                                   mean_elongation_rate);
		}
		else {
			push = Sample_uniform_distribution();
			if (push <= p_ribo_rnap_push) {
				result = Attempt_RNAP_move(time,
				         rnap_index,p_rnap_rnap_push,
                                         MAX_PROTEINS,proteins_made,
					 FORWARD);
				if ((abs((ribosome[ribo_index][POSITION] +
				     ok_dist)-rnap[rnap_index][POSITION]))
				     >= pol_ribo_width) {
					result = Advance_ribosome(time,
					  ribo_index,gene_index,rna_index,
					  num_proteins_made,MAX_PROTEINS,
					  proteins_made,num_codons,
                                          mean_elongation_rate);
				}
				else {
					Reschedule_ribosome_translocation(
					   time,ribo_index);
				}
			}
			else {
				Reschedule_ribosome_translocation(
				   time,ribo_index);
			}
		}
	}
	result++;
	return(1);			
}				

int Start_RNA_degradation(double time, int rna_index)
{
	double next_fire;
	Update_reaction_queue(INIT_DEGR_RXN(rna_index),INFINITY);
	Remove_reaction(INIT_DEGR_RXN(rna_index));
	next_fire = time + Find_firing_time(mu_continue_degrade);
	Add_reaction_to_queue(CONT_DEGR_RXN(rna_index),next_fire);
	final_rna_log[rna[rna_index][MASTER_INDEX]].start_deg = time;
	return(1);
}

int Degrade_first_RNA_nt(double time, int rna_index)
{
	int offset=0;
	double next_fire;
	if (translated == HYBRID) {
		offset += (split - 1);
	}
	if (rna[rna_index][MATURE] == 1 ||
            rnap[rna[rna_index][SOURCE_RNAP]][POSITION] >
	    degrade_rnap_width + offset + 2) {
		if (ribosome[rna[rna_index][NEWEST_RIBO[1]]][POSITION]== 0
		   || ribosome[rna[rna_index][NEWEST_RIBO[1]]][POSITION] >
		      degrade_ribo_width + offset + 2) {
			rna[rna_index][FIVE_END] = 2 + offset;
			if (rna[rna_index][FIVE_END] >= tsl_start[1] -
			    degrade_ribo_width) {
                                if (!(reaction_index[RIBO_LOAD_RXN(
				      rna_index,1)] < 0)) {			// ADJ SO THAT ALL STATES
                                        Update_reaction_queue(			// (FROM GRE ETC S/ BE OK)...
					   RIBO_LOAD_RXN(rna_index,1),
					     INFINITY);
                                        Remove_reaction(
					   RIBO_LOAD_RXN(rna_index,1));
                                }
                                rna[rna_index][TSL_ON_OFF[1]] = KILLED;
			}
			next_fire = time + Find_firing_time(
						mu_continue_degrade);
			Update_reaction_queue(CONT_DEGR_RXN(rna_index),
			   next_fire);
			rna_strip[rna_index][1+offset] = -1;
			return(1);
		}
	}
	next_fire = time + Find_firing_time(mu_continue_degrade);
	Update_reaction_queue(CONT_DEGR_RXN(rna_index),next_fire);
	return(0);
}

int Continue_RNA_degradation(double time, int rna_index, int *rna_log_ctr,
			     int **rna_log)
{
	int i, success;
	double next_fire;
	if (rna[rna_index][FIVE_END] == 1 || (translated == HYBRID &&
	    rna[rna_index][FIVE_END] == split)) {
		success = Degrade_first_RNA_nt(time,rna_index);
		return(success);
	}
	if (rna[rna_index][MATURE] == 1 ||
            (abs(rnap[rna[rna_index][SOURCE_RNAP]][POSITION] -
		 (rna[rna_index][FIVE_END] + 1)) > degrade_rnap_width)) {
		if (rna_strip[rna_index][rna[rna_index][FIVE_END] + 1 +
		    degrade_ribo_width] == 0 ||		
                    rna_strip[rna_index][rna[rna_index][FIVE_END] + 1 +
		    degrade_ribo_width] == -1) {
			rna_strip[rna_index][rna[rna_index][FIVE_END]]=-1;
			rna[rna_index][FIVE_END]++;

			if (rnap_relax && !(rna[rna_index][MATURE])) {
			   Update_RNAP_complex_resistance(
			      rna[rna_index][SOURCE_RNAP]);
			   success = Update_relaxation_rates(
			      time,rna[rna_index][SOURCE_RNAP]);
			}

			if (fish_five_set[rna[rna_index][FIVE_END]]) {
				radio_index =
				 ((((int) time)/radio_incr)+1)*radio_incr;
				radio_probe[radio_index][0]--;
			}
			if (fish_three_set[rna[rna_index][FIVE_END]]) {
				radio_index =
				 ((((int) time)/radio_incr)+1)*radio_incr;
				radio_probe[radio_index][1]--;
			}

			for (i = 1; i <= gene_num; i++) {
				if (rna[rna_index][TSL_ON_OFF[i]] ==
				    KILLED) {
					continue;
				} 
				if (rna[rna_index][FIVE_END] >=
				    tsl_start[i] - degrade_ribo_width) {	// ADJ SO THAT ALL STATES
                                     if (!(reaction_index[RIBO_LOAD_RXN(	// (FROM GRE ETC S/ BE OK)...
					    rna_index,i)] < 0)) {
					     Update_reaction_queue(
					      RIBO_LOAD_RXN(rna_index,i),
					      INFINITY);
					     Remove_reaction(
					      RIBO_LOAD_RXN(rna_index,i));
					}
					rna[rna_index][TSL_ON_OFF[i]] =
						KILLED;
				}
			}
                        if (rna[rna_index][FIVE_END] > tx_end ||
                            (rna[rna_index][MATURE] &&                                                        
                             rna[rna_index][FIVE_END] >
			     rna[rna_index][THREE_END])) {                     

				Update_reaction_queue(
				     CONT_DEGR_RXN(rna_index),INFINITY);
				Remove_reaction(
				     CONT_DEGR_RXN(rna_index));
				*rna_log_ctr += 1;
				rna_log[*rna_log_ctr][0] =
				     rna[rna_index][BIRTH_TIME];
				for (i = 1; i <= gene_num; i++) {
					rna_log[*rna_log_ctr][i] =
					    rna[rna_index][NUM_RIBO[i]];
				}
                                final_rna_log[
                                   rna[rna_index][MASTER_INDEX]].end_deg =
                                        time;
				memset(rna[rna_index], 0,
					RNA_INFO*sizeof(*rna[rna_index]));
				first_free_rna = rna_index;
			}
			else {
				next_fire = time +
					    Find_firing_time(
						mu_continue_degrade);
				Update_reaction_queue(
				    CONT_DEGR_RXN(rna_index),next_fire);
			}
			return(1);
		}
		else {
			next_fire = time +
				    Find_firing_time(mu_continue_degrade);
			Update_reaction_queue(
			    CONT_DEGR_RXN(rna_index),next_fire);
			return(0);
		}
	}
	else {
		next_fire = time + Find_firing_time(mu_continue_degrade);
		Update_reaction_queue(CONT_DEGR_RXN(rna_index),next_fire);
		return(0);
	}
}

