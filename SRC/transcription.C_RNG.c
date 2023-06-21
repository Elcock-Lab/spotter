#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define DEFINE_TX_VARIABLES
#include "INCL/transcription.C_RNG.h"

#include "INCL/sim_types.C_RNG.h"
#include "INCL/gen_def.C_RNG.h"
#include "INCL/sim_settings.C_RNG.h"
#include "INCL/reaction_manager.C_RNG.h"
#include "INCL/replication.C_RNG.h"
#include "INCL/regulation.C_RNG.h"
#include "INCL/translation.C_RNG.h"
#include "INCL/supercoiling.C_RNG.h"
#include "INCL/topo.C_RNG.h"
#include "INCL/utilities.C_RNG.h"
	

void Get_position_without_shutoff()
{
	int i;
	for (i = 1; i < MAX_RNAP; i++) {
		if (rnap[i][POSITION]) {
			shutoff_position[i] = rnap[i][POSITION];
		}
	}

	Update_reaction_queue(PROMO_KILL_RXN,INFINITY);
	Remove_reaction(PROMO_KILL_RXN);
	return;
}
	
int Find_first_free(int **look)							
{							
	int i;									// Compatible "topo_search" 
	for (i = 1; i <= 1000000; i++) {					// and "rnap_search" 
		if (look[i][1] == 0) {						// allow general use
			break;
		}
	}
	return(i);
}	

int Find_first_free_RNAP()
{
	int i;
	for (i = 1; i <= 1000000; i++) {
		if (rnap[i][1] == 0) {
			break;
		}
	}
	return(i);
}
	

int Check_for_RNAP_RNAP_contacts(int rnap_index, int updown)
{
	if (rnap[rnap_index][updown] != 0 &&
	    abs(rnap[rnap_index][POSITION] -
	    rnap[rnap[rnap_index][updown]][POSITION]) == pol_width) {
		return(1);
	}
	return(0);
}

int Find_RNAP_nearest_ribosome_distance(int rnap_index)
{
	if (!(rnap[rnap_index][CLOSEST_RIBO])) {
		return (-1000);
	}
	return ((rnap[rnap_index][POSITION] -
		 ribosome[rnap[rnap_index][CLOSEST_RIBO]][POSITION]) -
		pol_ribo_width);
}

int Check_for_RNAP_ribosome_contacts(int rnap_index,
						   int range)
{
	int gap;
	gap = Find_RNAP_nearest_ribosome_distance(rnap_index);
	if (gap < -999 || gap > range) {
		return(0);
	}
	return(1);
}

int Check_ribosome_tether(int rnap_index, int move)
{
	return(0);
/*
	if (Find_RNAP_nearest_ribosome_distance(rnap_index) + move >
	    TETHER_LENGTH) {
		return(0);
	}
	return(1);
*/
}

int Determine_dominant_trailer(int rnap_index, float rnap_cf,
			       float ribo_cf)
{
	int rnap_behind,rnap_trailer=0,ribo_trailer=0;
	rnap_behind = rnap[rnap_index][MINUS_ONE];
	if (rnap_behind &&
	    RNAP_ANTIPAUSE &&
	    rnap[rnap_behind][STATE] <= EFFECTOR_STATE &&
	    Check_for_RNAP_RNAP_contacts(rnap_index,MINUS_ONE)) {
		rnap_trailer++;
	}
	if (RIBO_ANTIPAUSE &&
	    Check_for_RNAP_ribosome_contacts(rnap_index,0)) {
		ribo_trailer++;
	}
	trailer_list[0] = rnap_trailer;
	trailer_list[1] = ribo_trailer;
	if (rnap_trailer) {
		if (ribo_trailer && ribo_cf > rnap_cf) {
			return(TRAILING_RIBO);
		}
		else {
			return(TRAILING_RNAP);
		}
	}
	else {
		if (ribo_trailer) {
			return(TRAILING_RIBO);
		}
		else {
			return(INDEPENDENT);
		}
	}
} 

int Advance_pause_state(double time, int rnap_index)
{
	int result;
	double next_fire,exit_fire,addl_fire,exit_mu;
	double term_fire,torque_adj,addl_mu;

	rnap[rnap_index][STATE]++;
	exit_mu = Torque_dependent_mu_exit(rnap_index);			
	if (exit_mu < 0) {
		result = Exit_off_pathway_state(time,rnap_index);
		return(1);
	}
	exit_fire = time + Find_firing_time(exit_mu);
	if (rnap[rnap_index][STATE] < OFF_P3) {					// Now keeping BOTH
		if (rnap[rnap_index][STATE] == OFF_P1) {			// pause advance AND
			next_fire = time +					// exit rxns for P1 and
				    Find_firing_time(mu_pause[OFF_P1]);		// P2;
			Add_reaction_to_queue(
			   RNAP_EXIT_RXN(rnap_index),exit_fire);
			if (exit_fire < INFINITY) {
				last_rates[rnap_index].last_e1_mu =
				   exit_mu;
			}
			else {
				last_rates[rnap_index].last_e1_mu =
				   -999.99;
			}
			if (!(reaction_index[RNAP_ON_F_RXN(rnap_index)]
			    < 0)) {
				Update_reaction_queue(
				   RNAP_ON_F_RXN(rnap_index),INFINITY);
			}
			term_fire = time + Find_firing_time(mu_term);		// Add'l stall-dependent
			Add_reaction_to_queue(					// termination rxn
			   RNAP_TERM_RXN(rnap_index),term_fire);
		}
		else {								// Decision about whether
			next_fire = time +					// or not ke2 and ke3 are
				    Find_firing_time(mu_pause[OFF_P2]);		// senstive to supercoiling
			Update_reaction_queue(					// in torque exit fn above
			   RNAP_EXIT_RXN(rnap_index),exit_fire);
			addl_fire = INFINITY;
			if (BACKTRACK_ON &&
			    rnap[rnap_index][POSITION] -
			    rna[rnap[rnap_index][RNAP_RNA]][FIVE_END]		// Check for adequate room to
			    > 15) {						// backtrack
			     addl_mu = rnap_dwell[
			      rnap[rnap_index][POSITION]].mu_b_off;
			     if (!(supercoiling_off)) {
			      torque_adj =
				Torque_dependent_mu_forward(rnap_index);
			      if (torque_adj < 0.001) {
				addl_mu /= 2.0;
			      }
			      else if (torque_adj > 1.999) {
				addl_mu = INFINITY;
			      }
			      else {
				torque_adj = 2.0 - torque_adj;
				addl_mu /= torque_adj;
			      }
			     }
			     addl_fire = time + Find_firing_time(
				addl_mu);
			}
			Add_reaction_to_queue(
	   		   RNAP_OFF_B_RXN(rnap_index),addl_fire);
			Add_reaction_to_queue(
	   		   RNAP_OFF_F_RXN(rnap_index),INFINITY);
		}
		if (!(reaction_index[RNAP_PAUSE_RXN(rnap_index)] < 0)) {
			Update_reaction_queue(RNAP_PAUSE_RXN(rnap_index),
			   next_fire);
		}
		else {
			Add_reaction_to_queue(RNAP_PAUSE_RXN(rnap_index),
			   next_fire);
		}	
	}
	else {
		Update_reaction_queue(
		   RNAP_PAUSE_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_PAUSE_RXN(rnap_index));
		Update_reaction_queue(
		   RNAP_EXIT_RXN(rnap_index),exit_fire);
	}
	result++;
	return(1);
}

int Process_pretranslocated_state(double time,int rnap_index)
{
	int result;
	double pause_fire,f_fire,p1_mu,f_mu;

	if (!(NO_PAUSE_MODEL)) {
	  if (!(ANTICASCADE_ON) ||
	    !(Check_for_RNAP_RNAP_contacts(rnap_index,PLUS_ONE))) {
		f_fire = time + Find_firing_time(
			 rnap_dwell[rnap[rnap_index][POSITION]].mu_f);
		p1_mu = Torque_dependent_mu_enter_P1(rnap_index);
	  }
	  else {
		f_fire = INFINITY;
		p1_mu = INFINITY;
	  } 

	  if (p1_mu > 0) {
		pause_fire = time + Find_firing_time(p1_mu);
	  }
	  else {
		result = Advance_pause_state(time,rnap_index);
		return(1);
	  }
	}
	else {
	  f_mu = Torque_dependent_mu_forward(rnap_index);
	  if (f_mu < 0.001) {
		f_fire = INFINITY;
		last_rates[rnap_index].last_f_mu = -999.99;
	  }
	  else {
		f_mu = rnap_dwell[rnap[rnap_index][POSITION]].mu_f/f_mu;
		f_fire = time + Find_firing_time(f_mu);
		last_rates[rnap_index].last_f_mu = f_mu;
	  }
	  p1_mu = INFINITY;
	  pause_fire = INFINITY;
	}

	if (p1_mu < INFINITY) {
		last_rates[rnap_index].last_p1_mu = p1_mu;
	}
	else {									// Negative value indicates
		last_rates[rnap_index].last_p1_mu = -999.99;			// that RNAP starts with
	}									// kp1 = 0
	if (!(reaction_index[RNAP_PAUSE_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_PAUSE_RXN(rnap_index),
		   pause_fire);
	}
	else {
		Add_reaction_to_queue(RNAP_PAUSE_RXN(rnap_index),
		   pause_fire);
	}	
	if (!(reaction_index[RNAP_ON_F_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_ON_F_RXN(rnap_index),
		   f_fire);
	}
	else {
		Add_reaction_to_queue(RNAP_ON_F_RXN(rnap_index),
		   f_fire);
	}
	result++;
	return(1);
}

int Exit_off_pathway_state(double time, int rnap_index)
{
	int rna_index,orig_three_end,result;
	rna_index = rnap[rnap_index][RNAP_RNA];					// This covers escapes
	orig_three_end = rna[rna_index][THREE_END];				// with and without
	rna[rna_index][THREE_END] = rnap[rnap_index][POSITION];			// greAB cutting

	if (orig_three_end != rna[rna_index][THREE_END]) {
		num_three_times[orig_three_end]++;
		three_times[orig_three_end][num_three_times[orig_three_end]] =
		   time - last_rnap_time[rnap_index];
		last_rnap_time[rnap_index] = time;
	}
		
	if (!(reaction_index[RNAP_EXIT_RXN(rnap_index)] < 0)) {		
		Update_reaction_queue(RNAP_EXIT_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_EXIT_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_OFF_B_RXN(rnap_index)] < 0)) {		
		Update_reaction_queue(RNAP_OFF_B_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_OFF_B_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_OFF_F_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_OFF_F_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_OFF_F_RXN(rnap_index));
	}
	rnap[rnap_index][STATE] = OFF_P1;					// Prevents inadvertent
	result = Determine_RNAP_riboload_effect(time,rnap_index,		// gene reactivation
		     orig_three_end,0);
	rnap[rnap_index][STATE] = PRE_TRANSLOC;
	if (!(reaction_index[RNAP_TERM_RXN(rnap_index)] < 0)) {			// Remove stall-dependent
		Update_reaction_queue(RNAP_TERM_RXN(rnap_index),INFINITY);	// termination rxn
		Remove_reaction(RNAP_TERM_RXN(rnap_index));
	}
	result = Process_pretranslocated_state(time,rnap_index);		// Should handle irrepsective
	result++;								// of cause: fn removes pause
	return(1);								// rxn
}

int Update_RNAP_rates(double time, int rnap_index)
{
	int reaction,result;
	double next_fire,new_mu, *old_mu, *time_adj;
	double torque_adj;

	result = Update_relaxation_rates(time,rnap_index);

	new_mu = Torque_dependent_mu_enter_P1(rnap_index);
	old_mu = &last_rates[rnap_index].last_p1_mu;
	time_adj = &last_rates[rnap_index].p1_time_adj;
	reaction = RNAP_PAUSE_RXN(rnap_index);
	if (!(rnap_index) ||							// Check for existence
	    rnap[rnap_index][STATE] == CLOSED ||				// of RNAP and status
	    rnap[rnap_index][STATE] == POST_TRANSLOC) {				// affected by torque
		return(0);
	}

	if (rnap[rnap_index][STATE] == PRE_TRANSLOC) {
	   if (!(NO_PAUSE_MODEL)) {

		if (new_mu < -1.0) {

			result = Advance_pause_state(time,rnap_index);

			last_rates[rnap_index].last_p1_mu = -999.9;
			return(1);
		}
		
		if (ANTICASCADE_ON &&
		    !(reaction_queue[reaction_index[				// Add forward time for
		      RNAP_ON_F_RXN(rnap_index)]].reaction_time			// RNAP newly liberated
		    < INFINITY) &&						// from downstream
		    !(Check_for_RNAP_RNAP_contacts(				// blockade
		       rnap_index,PLUS_ONE))) {
			Update_reaction_queue(RNAP_ON_F_RXN(rnap_index),
			 time +
			 Find_firing_time(
			   rnap_dwell[rnap[rnap_index][POSITION]].mu_f));
		}
	   }
	   else {
		reaction = RNAP_ON_F_RXN(rnap_index);
		torque_adj = Torque_dependent_mu_forward(rnap_index);
		if (torque_adj < 0.001) {
			new_mu = INFINITY;
		}
		else {
			new_mu =
			  rnap_dwell[rnap[rnap_index][POSITION]].mu_f / 
			  torque_adj;
		}
		old_mu = &last_rates[rnap_index].last_f_mu;
		time_adj = &last_rates[rnap_index].f_time_adj;
	   }
		
	}
	if (rnap[rnap_index][STATE] == OFF_P1) {
		reaction = RNAP_EXIT_RXN(rnap_index);
		new_mu = Torque_dependent_mu_exit(rnap_index);
		if (new_mu < -1.0) {
			result = Exit_off_pathway_state(time,rnap_index);
			last_rates[rnap_index].last_e1_mu = -999.9;
			return(1);
		}
		old_mu = &last_rates[rnap_index].last_e1_mu;
		time_adj = &last_rates[rnap_index].e1_time_adj;
	}

	if ((rnap[rnap_index][STATE] == OFF_P2 ||				// Here, assume that
	    rnap[rnap_index][STATE] == OFF_P3) &&				// torque does not affect
	    rnap[rnap_index][POSITION] == 					// ability of trailing
	    rna[rnap[rnap_index][RNAP_RNA]][THREE_END]) {			// RNAP or ribosome to 
		new_mu = Torque_dependent_mu_exit(rnap_index);			// promote transition to
		if (new_mu < -1.0) {						// on-pathway state; note
			result = Exit_off_pathway_state(time,rnap_index);	// that if torque prevents
			return(result);						// forward movement, though,
		}								// stalls are likely to
		else {								// recur.
			return(0);
		}
	}

	if (rnap[rnap_index][STATE] == OFF_P2 || 				// No adjustment to rates
	    rnap[rnap_index][STATE] == OFF_P3) {				// unless at start pos'n
		return(0);
	}


	next_fire =
	   reaction_queue[reaction_index[reaction]].reaction_time;
	if (new_mu > 9999.99) {							// If first zero rate:
		if (next_fire < INFINITY) {					// Old mu retained to keep
			*time_adj = next_fire - time;				// last non-zero rate
		}
		next_fire = INFINITY;
	}
	else {
		if (next_fire < INFINITY) {
			next_fire -= time;
			next_fire *= (new_mu/(*old_mu));			// Per G&B: new/old
			next_fire += time;					// mu b/c old/new
			*old_mu = new_mu;					// propensities
		}
		else if (*old_mu > -1.0) {					// If a pre-zero-rate
			next_fire = *time_adj;					// time exists, use it
			next_fire *= (new_mu/(*old_mu));			// to find tau
			next_fire += time;
			*old_mu = new_mu;
		}
		else {								// If first non-zero
			next_fire = time + Find_firing_time(new_mu);		// rate, draw new time
			*old_mu = new_mu;
		}
	}
	Update_reaction_queue(reaction,next_fire);
	return(1);
}

void Freeze_RNAP(double time, int rnap_index)
{
	int result;
	if (rnap[rnap_index][STATE] == CLOSED) {				// If "unfreezing" is 
		result = Remove_RNAP(time,rnap_index);				// implemented, will need
		return;
	}
	if (!(reaction_index[CLOSED_TO_OPEN_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(CLOSED_TO_OPEN_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(CLOSED_TO_OPEN_RXN(rnap_index));
	}
	if (!(reaction_index[UNBIND_RNAP_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(UNBIND_RNAP_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(UNBIND_RNAP_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_ON_F_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_ON_F_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_ON_F_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_ON_B_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_ON_B_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_ON_B_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_ADD_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_ADD_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_ADD_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_PAUSE_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_PAUSE_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_PAUSE_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_EXIT_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_EXIT_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_EXIT_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_OFF_F_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_OFF_F_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_OFF_F_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_OFF_B_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_OFF_B_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_OFF_B_RXN(rnap_index));
	}
	result++;
	rnap[rnap_index][STATE] = OFF_P3;
	rna[rnap[rnap_index][RNAP_RNA]][THREE_END] = 0;				// Assigned to prevent
	return;									// reactivation; FIX
}

void Unfreeze_RNAP(double time, int rnap_index)
{
	rnap[rnap_index][STATE] = PRE_TRANSLOC;
	rna[rnap[rnap_index][RNAP_RNA]][THREE_END] =
	   rnap[rnap_index][POSITION];
	Process_pretranslocated_state(time,rnap_index);
	return;
}

int Determine_off_pathway_move(double time,int rnap_index)
{
	double next_fire,torque_adjust;
	int rna_index, f_ok = 0, b_ok = 0;
	rna_index = rnap[rnap_index][RNAP_RNA];

	if (rnap[rnap_index][POSITION]+1 <= rna[rna_index][THREE_END]) {	// Offpathway transloc
		f_ok = 1;							// beyond 3' not allowed
	}
	if (rnap[rnap_index][POSITION] -					// Check for complete
	    rna[rnap[rnap_index][RNAP_RNA]][FIVE_END] > 15) {			// bubble post backtrack
		b_ok = 1;							// Note: no clash check
	}									// here (below)
	if (!(BACKTRACK_ON)){
		f_ok = 0;							// Switch backtracking
		b_ok = 0;							// on/off for testing
		return(0);
	}
	if (f_ok || b_ok) {
		if (supercoiling_off) {
			torque_adjust = 1.0;
		}
		else {
			torque_adjust =
			  Torque_dependent_mu_forward(rnap_index);
		}
	}
	next_fire = INFINITY;
	if (f_ok) {
	   if (torque_adjust < 0.001) {
		next_fire = INFINITY;
	   }
	   else {
		next_fire = time + Find_firing_time(
		    rnap_dwell[rnap[rnap_index][POSITION]].mu_f_off /
		    torque_adjust) * bt_adj;
	   }
	}
	Update_reaction_queue(RNAP_OFF_F_RXN(rnap_index),next_fire);
	torque_adjust = 2.0 - torque_adjust;
	next_fire = INFINITY;
	if (b_ok) {
	   if (torque_adjust < 0.001) {
		next_fire = INFINITY;
	   }
	   else {
		next_fire = time + Find_firing_time(
		    rnap_dwell[rnap[rnap_index][POSITION]].mu_b_off /
		    torque_adjust) * bt_adj;
	   }
	}
	Update_reaction_queue(RNAP_OFF_B_RXN(rnap_index),next_fire);
	return(1);
}

int Load_RNAP(double time, int promoter_index)
{
	double next_fire;
	float try_it;
        float promo_sc,sc_init_adj=1.0,p_init_adj;
        float k_sc_open;
	next_fire = time + Find_firing_time(mu_promo_load);

	if (!(ONE_RNAP_RUN && RNAP_run_started)) {
	   Update_reaction_queue(RNAP_LOAD_RXN(promoter_index),next_fire);
	}
	else {
	   Update_reaction_queue(RNAP_LOAD_RXN(promoter_index),INFINITY);
	   return(0);
	}

	if (promo_topo_blocked[promoter_index]) {				// Check for topoisomerase
										// preventing loading
	if (!(QUIET_MODE)) {
		printf("PROMOTER BLOCKED BY TOPOISOMERASE AT TIME %f\n",
			time);
	}

		return(0);
	}


	if (promoter[promoter_index][NEWEST_RNAP] == 0 ||
            abs(rnap[promoter[promoter_index][NEWEST_RNAP]][POSITION]-1)>=
	    pol_width + TSS_adj) {
		rnap[first_free_rnap][SOURCE_PROMO] = promoter_index;				
		rnap[first_free_rnap][POSITION] = 1;			
		rnap[first_free_rnap][MINUS_ONE] = 0;
		rnap[first_free_rnap][STATE] = CLOSED;
		if (promoter[promoter_index][NEWEST_RNAP] != 0) {
		   rnap[first_free_rnap][PLUS_ONE] =
			   promoter[promoter_index][NEWEST_RNAP];
		   rnap[promoter[promoter_index][NEWEST_RNAP]][MINUS_ONE]=
			 first_free_rnap;
		}
		else {
		   rnap[first_free_rnap][PLUS_ONE] = 0;
		   promoter[promoter_index][OLDEST_RNAP] = first_free_rnap;
		}
		promoter[promoter_index][NEWEST_RNAP] = first_free_rnap;
	
		p_init_adj = p_closed_to_open[promoter_index];
                if (TX_INIT_SC_FX) {
                  if (rnap[first_free_rnap][PLUS_ONE]) {
                    promo_sc = sigma[rnap[first_free_rnap][PLUS_ONE]]
                                    [UPSTREAM];
                  }
                  else {
                    promo_sc = (open_lk[promoter_index]/
                                relaxed_lk[stop_barrier - start_barrier])
                                - 1.0;
                  }
                  sc_init_adj = exp(-122.0*promo_sc);
                  if (sc_init_adj > 3.0) {
                    sc_init_adj = 3.0;
                  }
                  if (sc_init_adj > 1.0/p_closed_to_open[promoter_index]) {
                    sc_init_adj = (1.0/p_closed_to_open[promoter_index]) -
                                  0.0001;
                  }
                  k_sc_open = (p_closed_to_open[promoter_index] *
                               k_unbind) /
                              ((1.0 / sc_init_adj) -
                               p_closed_to_open[promoter_index]);
                  next_fire = time + Find_firing_time(
                                        1.0/(k_unbind + k_sc_open));
		  p_init_adj *= sc_init_adj;
                }
                else {
                  next_fire = time + Find_firing_time(mu_net_promo);
                }

		try_it = Sample_uniform_distribution();
		if (try_it < p_init_adj) {
			Add_reaction_to_queue(
			  CLOSED_TO_OPEN_RXN(first_free_rnap),next_fire);
		}
		else if (failed_loads > max_RNAP_load_gap) {
			Add_reaction_to_queue(
			   CLOSED_TO_OPEN_RXN(first_free_rnap),next_fire);
		}
		else {
			Add_reaction_to_queue(
			   UNBIND_RNAP_RXN(first_free_rnap),next_fire);
			if (p_closed_to_open[promoter_index] > 0.001) {
				failed_loads++;
			}
		}
		dna_strip[promoter_index][1] = first_free_rnap;
		rnap_bound++;
		first_free_rnap = Find_first_free(rnap_search);
		return(1);
	}
	return(0);
}



int Activate_RNAP(double time, int rnap_index)
{
	int i,result,promo_index,rnap_ahead;
	int last_rnap;
	int check_topo;
	int whole_sep,left_sep,right_sep,wrapped=0;
	double next_fire;
	float old_lk,lk_whole,upstr_share,downstr_share;

	Update_reaction_queue(CLOSED_TO_OPEN_RXN(rnap_index),INFINITY);
	Remove_reaction(CLOSED_TO_OPEN_RXN(rnap_index));
	rnap[rnap_index][STATE] = PRE_TRANSLOC;
	rnap[rnap_index][RNAP_RNA] = first_free_rna;

	sim_total_rnap++;
	rnap[rnap_index][SUPER_INDEX] = sim_total_rnap;
	sim_total_rna++;
	rna[first_free_rna][MASTER_INDEX] = sim_total_rna;

	final_rna_log[sim_total_rna].start_tx = time;

	rna[first_free_rna][SOURCE_PROMO]= rnap[rnap_index][SOURCE_PROMO];
	if (stable == OFF) {
		next_fire = time + Find_firing_time(mu_start_degrade);
		Add_reaction_to_queue(INIT_DEGR_RXN(first_free_rna),
		   next_fire);
	}

	next_fire = time + Find_firing_time(mu_rnap_surveillance);
	Add_reaction_to_queue(RNAP_SURVEILLANCE_RXN(rnap_index),next_fire);

	Add_reaction_to_queue(RNAP_RELAX_RXN(rnap_index),INFINITY);		// Net torque always = 0
	last_rates[rnap_index].last_assist_mu = -9999.9;			// so put in queue but not
	last_rates[rnap_index].last_resist_mu = -9999.9;			// scheduled at loading
	complex_resist[rnap_index] = 1.0/(0.12+0.05);

	rna[first_free_rna][SOURCE_RNAP] = rnap_index;
	rna[first_free_rna][FIVE_END] = 1;
	rna[first_free_rna][THREE_END] = 1;
	rna[first_free_rna][BIRTH_TIME] = (int) (time + 0.5);
	rna_strip[first_free_rna][1] = 0;
	first_free_rna = Find_first_free(rna);

	curr_dwt_window[rnap_index] = -1;
	RNAP_run_started++;
	trace_RNAP = rnap_index;

	promo_index = rnap[rnap_index][SOURCE_PROMO];
	rnap_ahead = rnap[rnap_index][PLUS_ONE];
	last_rnap = promoter[promo_index][OLDEST_RNAP];

	curr_rnap_rot[rnap_index] = curr_rnap_rot[rnap_ahead];


	if (rnap_ahead) {

		if (!(QUIET_MODE)) {
			printf("ACTIVATING RNAP %d; "
				"THERE IS ANOTHER RNAP IN FRONT: %d\n",
				rnap_index,rnap[rnap_index][PLUS_ONE]);
		}

		rnap_sep[rnap_index][DOWNSTREAM] =
		   rnap[rnap[rnap_index][PLUS_ONE]][POSITION] - 
		   1;
		rnap_sep[rnap_index][UPSTREAM] =
		   rnap_sep[rnap_ahead][UPSTREAM] -
		   rnap_sep[rnap_index][DOWNSTREAM];

		if (!(QUIET_MODE)) {
			printf("RNAP SEP UP IS %d, DOWN IS %d\n",
				rnap_sep[rnap_index][UPSTREAM],
				rnap_sep[rnap_index][DOWNSTREAM]);
		}

		loading_sigma[rnap_index] =
		   sigma[rnap[rnap_index][PLUS_ONE]][UPSTREAM];

		if (!(QUIET_MODE)) {
			printf("LOADING SIGMA IS %f\n",
				loading_sigma[rnap_index]);
			printf("OLDER RNAP LK UP IS %f; DOWN IS %f\n",
			lk[rnap[rnap_index][PLUS_ONE]][UPSTREAM],
			lk[rnap[rnap_index][PLUS_ONE]][DOWNSTREAM]);
		}

		if (BUBBLE_ADJ) {
		 rnap_sep[rnap_index][DOWNSTREAM] -=
		    (bubble[UPSTREAM] + bubble[DOWNSTREAM]);
		  old_lk = lk[rnap_ahead][UPSTREAM];
		  upstr_share = ((float) rnap_sep[rnap_index][UPSTREAM]) /
				((float) (rnap_sep[rnap_index][UPSTREAM] +
				       rnap_sep[rnap_index][DOWNSTREAM]));
		  upstr_share *= old_lk;
		  downstr_share = old_lk - upstr_share;
		  lk[rnap_index][UPSTREAM] = upstr_share;
		  lk[rnap_index][DOWNSTREAM] = downstr_share;
		  rnap_sep[rnap_ahead][UPSTREAM] =
			rnap_sep[rnap_index][DOWNSTREAM];
		  lk[rnap_ahead][UPSTREAM] = lk[rnap_index][DOWNSTREAM];
		}
		
		else {
		 lk[rnap_index][UPSTREAM] =
		   lk[rnap[rnap_index][PLUS_ONE]][UPSTREAM] *			// Note that adjustments to
		   (((float) rnap_sep[rnap_index][UPSTREAM])/			// the oldest-RNAP Lk are not
		   ((float) rnap_sep[rnap[rnap_index][PLUS_ONE]][UPSTREAM]));	// made here for cases where
		 lk[rnap[rnap_index][PLUS_ONE]][UPSTREAM] -=			// it is linked to newly-
		   lk[rnap_index][UPSTREAM];					// activated RNAP: these are
		 lk[rnap_index][DOWNSTREAM] =					// made in the call to
		   lk[rnap[rnap_index][PLUS_ONE]][UPSTREAM];			// "Process" below
		 rnap_sep[rnap[rnap_index][PLUS_ONE]][UPSTREAM] =
		   rnap_sep[rnap_index][DOWNSTREAM];
		}

		if (!(QUIET_MODE)) {
			printf("LK UP IS %f; LK DOWN IS %f\n",
			lk[rnap_index][UPSTREAM],
			lk[rnap_index][DOWNSTREAM]);
		}

	}
	else {
		whole_sep = stop_barrier -	
			    start_barrier;
		if (plasmid) {
			left_sep  = whole_sep;
			right_sep = whole_sep;
		}
		else {
			left_sep  = 1 - start_barrier;
			right_sep = stop_barrier - 1;
		}	
		lk_whole  = open_lk[promo_index];

		if (!(wrapped)) {
			rnap_sep[rnap_index][UPSTREAM]   = left_sep;
			rnap_sep[rnap_index][DOWNSTREAM] = right_sep;
			if (BUBBLE_ADJ) {
				rnap_sep[rnap_index][UPSTREAM]   -=
					bubble[UPSTREAM];
				rnap_sep[rnap_index][DOWNSTREAM] -=
					bubble[DOWNSTREAM];
                                if (plasmid) {			                // If plasmid without RNAP
                                     rnap_sep[rnap_index][UPSTREAM]   -=        // ahead, both directions
                                        bubble[DOWNSTREAM];                     // adjusted by full bubble
                                     rnap_sep[rnap_index][DOWNSTREAM] -=        // length
                                        bubble[UPSTREAM];
                                }
				whole_sep -=
				 (bubble[UPSTREAM] + bubble[DOWNSTREAM]);
			}
			lk[rnap_index][UPSTREAM]   = lk_whole *
			  (((float) left_sep)/((float) whole_sep));
			lk[rnap_index][DOWNSTREAM]   = lk_whole *
			  (((float) right_sep)/((float) whole_sep));
		}
		init_sigma[promo_index] = (lk[rnap_index][UPSTREAM]/
		   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]) - 1.0;
		curr_init_lk_up[promo_index]=lk[rnap_index][UPSTREAM];
		curr_init_lk_down[promo_index]=lk[rnap_index][DOWNSTREAM];
		loading_sigma[rnap_index] = init_sigma[promo_index];
	}

	if (!(QUIET_MODE)) {
		printf("ACTIVATING RNAP: ASSIGNING START LK'S OF %f {UP} "
			"AND %f {DOWN}\n",
			lk[rnap_index][UPSTREAM],
			lk[rnap_index][DOWNSTREAM]);
	}

	Process_supercoiling_changes(rnap_index,0,BOTH);

	result = Process_pretranslocated_state(time,rnap_index);
	Update_RNAP_rates(time,rnap[rnap_index][PLUS_ONE]);			// Only downstream RNAP
										// needs to be dealt
										// with (b/c at promoter)

	rnap[rnap_index][TOPO_UPSTREAM] =					// Current auto-assignments
	  starter_up[promo_index];						// determined on topo binding
	topo[starter_up[promo_index]][CLOSEST_UP] = rnap_index;
	if (rnap[rnap[rnap_index][PLUS_ONE]][TOPO_UPSTREAM] ==
	    starter_up[promo_index]) {
		rnap[rnap[rnap_index][PLUS_ONE]][TOPO_UPSTREAM] = 0;
	}
	if (!(rnap[rnap_index][PLUS_ONE]) ||
	    rnap[rnap[rnap_index][PLUS_ONE]][POSITION] >
	    topo[starter_down[promo_index]][POSITION]) {
		rnap[rnap_index][TOPO_DOWNSTREAM] =					
		  starter_down[promo_index];
		topo[starter_down[promo_index]][CLOSEST_DOWN] =
		   rnap_index;
	}

	for (i = 1; i <= num_topo[promo_index]; i++) {
		check_topo = topo_finder[promo_index][i];
		if (topo[check_topo][POS_BLOCK] == UPSTREAM) {
			topo[check_topo][PLUS_ONE] = rnap_index;
		}
		else {
			if (!(topo[check_topo][MINUS_ONE])) {
				topo[check_topo][MINUS_ONE] =
				   rnap_index;
			}
		}
	}

	if (!(QUIET_MODE)) {
		printf("ACTIVATED RNAP %d ASSIGNED TOPO DOWN %d AND TOPO "
			"UP %d\n",
			rnap_index,rnap[rnap_index][TOPO_DOWNSTREAM],
			rnap[rnap_index][TOPO_UPSTREAM]);
	}
		

	num_started++;
	tx_times[rnap_index] = time;
	last_rnap_time[rnap_index] = time;
	last_rnap--;
	printf("RNAP %d BEGINS TRANSCRIPTION ON PROMOTER %d AT TIME %f\n",
		rnap_index,rnap[rnap_index][SOURCE_PROMO],time);
	failed_loads = 0;
	printf("RESETTING ACTIVATION ATTEMPTS TO ZERO\n");
	result++;
	return(1);
}
	

void Generate_RNAP_reaction_set(int rnap_index,int *rxn_set)
{
        rxn_set[1] = UNBIND_RNAP_RXN(rnap_index);
        rxn_set[2] = RNAP_ON_F_RXN(rnap_index);
        rxn_set[3] = RNAP_ON_B_RXN(rnap_index);
        rxn_set[4] = RNAP_ADD_RXN(rnap_index);
        rxn_set[5] = RNAP_PAUSE_RXN(rnap_index);
        rxn_set[6] = RNAP_EXIT_RXN(rnap_index);
        rxn_set[7] = RNAP_OFF_F_RXN(rnap_index);
        rxn_set[8] = RNAP_OFF_B_RXN(rnap_index);
        return;
}

int Process_terminating_RNAP(double time, int rnap_index)
{
        int i,rnap_behind,rnap_ahead,loc;
        int ref_rnap = 0,newest_rnap,oldest_rnap,plasmid_linked_rnap=0;
        int oldest_leaving=0,newest_leaving=0;
        int promo_index;
        int rxn_set[10];
        float temp_sigma=0.0;

	if (!(QUIET_MODE)) {
		printf("PROCESSING RNAP TERMINATION OF %d AT TIME %f\n",
			rnap_index,time);
	}

        dna_strip[
         rnap[rnap_index][SOURCE_PROMO]][rnap[rnap_index][POSITION]] = 0;
        rnap_behind  = rnap[rnap_index][MINUS_ONE];
        rnap_ahead   = rnap[rnap_index][PLUS_ONE];
        loc = rnap[rnap_index][POSITION];
        promo_index = rnap[rnap_index][SOURCE_PROMO];
        newest_rnap =
           promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP];
        if (rnap_index == newest_rnap) {
                promoter[promo_index][NEWEST_RNAP] = rnap_ahead;
        }
        if (rnap[newest_rnap][STATE] == CLOSED) {
                newest_rnap = rnap[newest_rnap][PLUS_ONE];
        }
        if (rnap_index == newest_rnap) {
                newest_leaving++;
        }
        oldest_rnap = promoter[promo_index][OLDEST_RNAP];
        if (rnap_index == oldest_rnap) {
                promoter[promo_index][OLDEST_RNAP] = rnap_behind;
                oldest_leaving++;
        }
        if (rnap_behind) {
                rnap[rnap_behind][PLUS_ONE] = rnap[rnap_index][PLUS_ONE];
        }
        if (rnap_ahead) {
                rnap[rnap_ahead][MINUS_ONE] = rnap[rnap_index][MINUS_ONE];
        }
        rnap_behind *= (rnap[rnap_behind][STATE] != CLOSED);
        if (rnap_behind && !(rnap[rnap_behind][TOPO_DOWNSTREAM])) {
                rnap[rnap_behind][TOPO_DOWNSTREAM] =                            // If RNAP behind and it doesn't
                   rnap[rnap_index][TOPO_DOWNSTREAM];                           // have a blocking topoisomerase,
                topo[rnap[rnap_behind][TOPO_DOWNSTREAM]][CLOSEST_DOWN] =        // it inherits that of the RNAP
                   rnap_behind;                                                 // departing...
        }
        if (rnap_ahead && !(rnap[rnap_ahead][TOPO_UPSTREAM])) {                 // ...as the RNAP ahead inherits
                rnap[rnap_ahead][TOPO_UPSTREAM] =                               // an upstream topoisomerase if
                   rnap[rnap_index][TOPO_UPSTREAM];                             // it doesn't have one already
                topo[rnap[rnap_ahead][TOPO_UPSTREAM]][CLOSEST_UP] =
                   rnap_ahead;
        }
        rnap_behind *= (rnap[rnap_index][BARRIER_BEHIND] <= 0);                 // Checks for inter-RNAP linkage
        rnap_ahead  *= (rnap[rnap_index][BARRIER_AHEAD] <= 0);
        if (rnap_behind) {
                rnap_sep[rnap_behind][DOWNSTREAM] +=
                 rnap_sep[rnap_index][DOWNSTREAM];
                lk[rnap_behind][DOWNSTREAM] +=
                 lk[rnap_index][DOWNSTREAM];
		if (BUBBLE_ADJ) {
			rnap_sep[rnap_behind][DOWNSTREAM] +=
			  (bubble[UPSTREAM] + bubble[DOWNSTREAM]);
		}
                if (rnap[rnap_index][BARRIER_AHEAD] > 0) {
                        rnap[rnap_behind][BARRIER_AHEAD] =
                        rnap[rnap_index][BARRIER_AHEAD];
                }
                ref_rnap = rnap_behind;
                Process_supercoiling_changes(rnap_behind,0,DOWNSTREAM);
        }
        if (rnap_ahead) {
                if (!(ref_rnap)) {                                              // Changes will have been
                     rnap_sep[rnap_ahead][UPSTREAM] +=                          // made already if no
                       rnap_sep[rnap_index][UPSTREAM];                          // barrier between RNAPs
                     lk[rnap_ahead][UPSTREAM] +=                                // ahead and behind
                       lk[rnap_index][UPSTREAM];
		     if (BUBBLE_ADJ) {
			rnap_sep[rnap_ahead][UPSTREAM] +=
			  (bubble[UPSTREAM] + bubble[DOWNSTREAM]);
		     }
                     if (rnap[rnap_index][BARRIER_BEHIND] > 0) {
                        rnap[rnap_ahead][BARRIER_BEHIND] =
                        rnap[rnap_index][BARRIER_BEHIND];
                     }
                     ref_rnap = rnap_ahead;
                     Process_supercoiling_changes(rnap_ahead,0,UPSTREAM);
                }
        }

        if (plasmid) {
                if (oldest_leaving && !(newest_leaving) &&
                    rnap[rnap_index][BARRIER_AHEAD] <= 0 &&
                    rnap[newest_rnap][BARRIER_BEHIND] <= 0) {
                        plasmid_linked_rnap = newest_rnap;
                        if (!(rnap_behind)) {              		        // Not necessary if linked
                                rnap_sep[newest_rnap][UPSTREAM] +=              // RNAP is RNAP behind; if so,
                                   rnap_sep[rnap_index][UPSTREAM];              // it will have been handled
                                lk[newest_rnap][UPSTREAM] +=                    // above
                                   lk[rnap_index][UPSTREAM];
				if (BUBBLE_ADJ) {
				  rnap_sep[newest_rnap][UPSTREAM] +=
				   (bubble[UPSTREAM]+bubble[DOWNSTREAM]);
				}
                                Process_supercoiling_changes(
                                   newest_rnap,0,UPSTREAM);
                                if (!(ref_rnap)) {
                                        ref_rnap = newest_rnap;
                                }
                        }
                }
                if (newest_leaving && !(oldest_leaving) &&
                    rnap[rnap_index][BARRIER_BEHIND] <= 0 &&
                    rnap[oldest_rnap][BARRIER_AHEAD] <= 0) {
                        plasmid_linked_rnap = oldest_rnap;
                        if (!(rnap_ahead)) {		                        // Similarly, if linked
                                rnap_sep[oldest_rnap][DOWNSTREAM] +=            // RNAP is RNAP ahead,
                                   rnap_sep[rnap_index][DOWNSTREAM];            // adjustment not needed
                                lk[oldest_rnap][DOWNSTREAM] +=
                                   lk[rnap_index][DOWNSTREAM];
				if (BUBBLE_ADJ) {
				  rnap_sep[oldest_rnap][DOWNSTREAM] +=
				   (bubble[UPSTREAM]+bubble[DOWNSTREAM]);
				}
                                Process_supercoiling_changes(
                                   oldest_rnap,0,DOWNSTREAM);
                                if (!(ref_rnap)) {
                                        ref_rnap = oldest_rnap;
                                }
                        }
                }
        }

        if (rnap_behind) {
        	Update_RNAP_rates(time,rnap_behind);
        }
        if (rnap_ahead) {
        	Update_RNAP_rates(time,rnap_ahead);
        }
        if (plasmid_linked_rnap &&
            !(plasmid_linked_rnap == rnap_behind ||
            plasmid_linked_rnap == rnap_ahead)) {
        	Update_RNAP_rates(time,plasmid_linked_rnap);
        }

        if (rnap_behind) {
                Update_topoisomerase_lk(lk[rnap_behind][DOWNSTREAM],
                                        sigma[rnap_behind][DOWNSTREAM],
                                        rnap_sep[rnap_behind][DOWNSTREAM],
                                        rnap[rnap_behind][POSITION],
                                        rnap[rnap_behind][POSITION] +
                                        rnap_sep[rnap_behind][DOWNSTREAM],
                                        MIDDLE,promo_index,rnap_index);
        }
        else if (rnap_ahead) {
                Update_topoisomerase_lk(lk[rnap_ahead][UPSTREAM],
                                        sigma[rnap_ahead][UPSTREAM],
                                        rnap_sep[rnap_ahead][UPSTREAM],
                                        rnap[rnap_ahead][POSITION] -
                                        rnap_sep[rnap_ahead][UPSTREAM],
                                        rnap[rnap_ahead][POSITION],
                                        MIDDLE,promo_index,rnap_index);
        }
        else {
            if (!(ref_rnap)) {
                temp_sigma =							// Account for
                   (open_lk[promo_index]/relaxed_open_lk) - 1.0;		// equilibration of upstream
                Update_topoisomerase_lk(open_lk[promo_index],			// and downstream sigma even
                                        temp_sigma,				// if no RNAPs left
                                        stop_barrier - start_barrier,
                                        start_barrier,
                                        stop_barrier,
                                        MIDDLE, promo_index,
                                        rnap_index);
            }
        }

        if (plasmid_linked_rnap && plasmid_linked_rnap == newest_rnap) {
                Update_topoisomerase_lk(lk[newest_rnap][UPSTREAM],              // Adjust topoisiomerase binding
                                        temp_sigma,                             // from start to newest RNAP
                                        rnap_sep[newest_rnap][UPSTREAM],
                                        start_barrier,
                                        rnap[newest_rnap][POSITION],
                                        MIDDLE,promo_index,rnap_index);
                if (!(rnap_behind)) {
                        Update_topoisomerase_lk(lk[newest_rnap][UPSTREAM],      // Check region from right
                                        temp_sigma,                             // edge to stop barrier as well;
                                        rnap_sep[newest_rnap][UPSTREAM],        // if RNAP behind exists, this
                                        rnap[rnap_index][POSITION] -            // will have been done above
                                        rnap_sep[rnap_index][UPSTREAM],
                                        stop_barrier,
                                        MIDDLE,promo_index,rnap_index);
                }
        }
        if (plasmid_linked_rnap && plasmid_linked_rnap == oldest_rnap) {        // Adjustment from oldest RNAP
                Update_topoisomerase_lk(lk[oldest_rnap][DOWNSTREAM],            // to transcript end
                                        temp_sigma,
                                        rnap_sep[oldest_rnap][DOWNSTREAM],
                                        rnap[oldest_rnap][POSITION],
                                        stop_barrier,
                                        MIDDLE,promo_index,rnap_index);
                if (!(rnap_ahead)) {                                            // Check region from start
                        Update_topoisomerase_lk(lk[oldest_rnap][DOWNSTREAM],    // barrier to RNAP downstream
                                        temp_sigma,                             // from departing RNAP or
                                        rnap_sep[oldest_rnap][DOWNSTREAM],      // left edge (always the
                                        start_barrier,                          // latter, else covered above)
                                        rnap[rnap_index][POSITION] +
                                        rnap_sep[rnap_index][DOWNSTREAM],
                                        MIDDLE,promo_index,rnap_index);
                }
        }

        rnap_sep[rnap_index][UPSTREAM]   = 0;
        rnap_sep[rnap_index][DOWNSTREAM] = 0;
        lk[rnap_index][UPSTREAM]   = 0.0;
        lk[rnap_index][DOWNSTREAM] = 0.0;

        Generate_RNAP_reaction_set(rnap_index,rxn_set);                         // Purge of all rxns for
        for (i = 1; i <= 8; i++) {                                              // the terminating RNAP

		printf("TERM RXN: %d; INDEX: %d; SCHEDULED TIME: %f\n",
			rxn_set[i],reaction_index[i],
			reaction_queue[reaction_index[i]].reaction_time);

                if (!(reaction_index[rxn_set[i]] < 0)) {
                        Update_reaction_queue(rxn_set[i],INFINITY);
                        Remove_reaction(rxn_set[i]);
                }
        }

	loc++;
        final_rna_log[rna[rnap[
           rnap_index][RNAP_RNA]][MASTER_INDEX]].end_tx = time * -1.0;
        memset(rnap[rnap_index],0,RNAP_INFO*sizeof(*rnap[rnap_index]));
        rnap_bound--;
        first_free_rnap = rnap_index;
        return(1);
}

int Terminate_transcription(double time, int rnap_index)
{

	if (!(QUIET_MODE)) {
		printf("TERMINATING TRANSCRIPTION OF RNAP %d AT TIME %f\n",
		rnap_index,time);
	}

        int rna_index;
        rna_index = rnap[rnap_index][RNAP_RNA];
        rna[rna_index][MATURE] = 1;                                             // Note "MATURE" also used
        Process_terminating_RNAP(time,rnap_index);                              // for terminated RNAs
	if (!(reaction_index[RNAP_TERM_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_TERM_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_TERM_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_EVAL_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_EVAL_RXN(rnap_index),INFINITY);
		Remove_reaction(RNAP_EVAL_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_SURVEILLANCE_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_SURVEILLANCE_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_SURVEILLANCE_RXN(rnap_index));
	}
	if (!(reaction_index[RNAP_RELAX_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(RNAP_RELAX_RXN(rnap_index),
				      INFINITY);
		Remove_reaction(RNAP_RELAX_RXN(rnap_index));
	}

        return(1);
}

			
int Schedule_RNAP_checkup(double time, int rnap_index)
{
	last_rates[rnap_index].checkstart = time;
	last_rates[rnap_index].checkpos = rnap[rnap_index][POSITION];
	Add_reaction_to_queue(RNAP_EVAL_RXN(rnap_index),time+mu_rnap_eval);
	Update_reaction_queue(RNAP_SURVEILLANCE_RXN(rnap_index),INFINITY);
	Remove_reaction(RNAP_SURVEILLANCE_RXN(rnap_index));
	return(1);
}

int Evaluate_transcriptional_progress(double time,int rnap_index)
{
	double mfd_cf_rate, next_fire;
	Update_reaction_queue(RNAP_EVAL_RXN(rnap_index),INFINITY);
	Remove_reaction(RNAP_EVAL_RXN(rnap_index));
	mfd_cf_rate = ((double) (rnap[rnap_index][POSITION] -
				last_rates[rnap_index].checkpos)) /
		      (time - last_rates[rnap_index].checkstart);
	if (mfd_cf_rate > 7.0) {
		next_fire = time + Find_firing_time(mu_rnap_surveillance);
		Add_reaction_to_queue(RNAP_SURVEILLANCE_RXN(rnap_index),
				      next_fire);
		return(0);
	}
	Terminate_transcription(time, rnap_index);
	return(1);
}

int Remove_RNAP(double time, int rnap_index)
{
	dna_strip[rnap[rnap_index][SOURCE_PROMO]][
		rnap[rnap_index][POSITION]] = 0;
	rnap[rnap[rnap_index][PLUS_ONE]][MINUS_ONE] = 0;
	if (promoter[rnap[rnap_index][SOURCE_PROMO]][OLDEST_RNAP] == 
	    rnap_index) {
		promoter[rnap[rnap_index][SOURCE_PROMO]][OLDEST_RNAP] = 0;
	}
	promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP] =
	   rnap[rnap_index][PLUS_ONE];
	if (!(reaction_index[UNBIND_RNAP_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(UNBIND_RNAP_RXN(rnap_index),
				      INFINITY);				// Conditional to deal with
		Remove_reaction(UNBIND_RNAP_RXN(rnap_index));			// forced removal
	}									

	if (!(reaction_index[CLOSED_TO_OPEN_RXN(rnap_index)] < 0)) {
		Update_reaction_queue(CLOSED_TO_OPEN_RXN(rnap_index),
		   INFINITY);
		Remove_reaction(CLOSED_TO_OPEN_RXN(rnap_index));
	}

	memset(rnap[rnap_index],0,RNAP_INFO*sizeof(*rnap[rnap_index]));
	rnap_bound--;
	first_free_rnap = rnap_index;
	return(1);
}

int Rectify_downstream_RNAP(double time, int rnap_index)			// Vestigial: will need to be
{										// changed to be used...
	int reaction;
	float orig_dwell,next_fire;
	reaction = RNAP_MOVE_RXN(rnap_index);
	orig_dwell = mu_rnap_move[rnap[rnap_index][POSITION]];
	next_fire = reaction_queue[reaction_index[reaction]].reaction_time;
	next_fire -= time;
	next_fire *= (min_rnap_dwell/orig_dwell);
	next_fire += time;
	Update_reaction_queue(RNAP_MOVE_RXN(rnap_index),next_fire);
	return(1);
}

int Add_RNA_nt(double time, int rnap_index, int MAX_PROTEINS,
               float **proteins_made)
{
	int i,start,result,loc,rna_index,promo_index;
	int new_old,newest_rnap;
	int checkster,wrap_ends[2];
	int windex;
	double next_fire;
	float temp_sigma;
	if (!(reaction_index[RNAP_ADD_RXN(rnap_index)] < 0)) {			// Check allows direct call
		Update_reaction_queue(RNAP_ADD_RXN(rnap_index),INFINITY);	// at transcript end without
		Remove_reaction(RNAP_ADD_RXN(rnap_index));			// rxn assignment
	}
	if (stable == ON) {
		start = gene_num;
		if (translated == HYBRID) {
			start--;
		}
		for (i = start; i >= 1; i--) {
			if (rnap[rnap_index][POSITION] == tsl_stop[i]+1) {	// Using proteins category for
				num_proteins_made[i]++;				// stable RNA genes
				proteins_made[i][num_proteins_made[i]] =
				   time;
				product_added[i]++;
				break;
			}
		}
	}
	if (rnap[rnap_index][POSITION] == split) {
		rna[rnap[rnap_index][RNAP_RNA]][THREE_END] = split-1;
		rna[rnap[rnap_index][RNAP_RNA]][MATURE] = 1;
		rnap[rnap_index][RNAP_RNA] = first_free_rna;
		rna[first_free_rna][SOURCE_PROMO] =
		   rnap[rnap_index][SOURCE_PROMO];
		rna[first_free_rna][SOURCE_RNAP] = rnap_index;
		rna[first_free_rna][FIVE_END] = split;
		rna[first_free_rna][THREE_END] = split;
		rna[first_free_rna][BIRTH_TIME] = (int) (time + 0.5);
		next_fire = time + Find_firing_time(mu_start_degrade);
		Add_reaction_to_queue(INIT_DEGR_RXN(first_free_rna),
		   next_fire);
		first_free_rna = Find_first_free(rna);
	}	

	if (rnap[rnap_index][POSITION] == tsl_stop[1]) {
		gene_time[rnap_index] = time;
	}

	loc = rnap[rnap_index][POSITION];
	num_three_times[loc - 1]++;
	three_times[loc - 1][num_three_times[loc - 1]] =
	   time - last_rnap_time[rnap_index];
	last_rnap_time[rnap_index] = time;
		
	rnap[rnap_index][STATE] = PRE_TRANSLOC;
	loc = rnap[rnap_index][POSITION];
	rna_index = rnap[rnap_index][RNAP_RNA];
	if (rnap[rnap_index][POSITION] <= tx_end) {
		rna[rnap[rnap_index][RNAP_RNA]][THREE_END] =
		   rnap[rnap_index][POSITION];

		if (!(QUIET_MODE)) {
		  printf("JUST ADDED NT FOR RNAP AT %d AT TIME %f\n",
			rnap_index,time);
		}

		result = Process_pretranslocated_state(time,rnap_index);	// Above no longer needed:
										// contacts will checked
										// on processing

		if (!(QUIET_MODE)) {
		  printf("DONE WITH PROCESSING OF "
			 "PRETRANSLOCATED STATE\n");
		}
			

		Determine_RNAP_riboload_effect(time,rnap_index,
		   rna[rna_index][THREE_END],1);

		if (fish_five_set[rnap[rnap_index][POSITION]]) {
			radio_index =
			   ((((int) time)/radio_incr)+1)*radio_incr;
			radio_probe[radio_index][0]++;
		}
		if (fish_three_set[rnap[rnap_index][POSITION]]) {
			radio_index =
			   ((((int) time)/radio_incr)+1)*radio_incr;
			radio_probe[radio_index][1]++;
		}

		if (rnap_relax) {
			Update_RNAP_complex_resistance(rnap_index);
		}

		if (REPORT_DWT) {
			windex = (rnap[rnap_index][POSITION] - dwt_start)/
				 dwt_window;
			if (windex > -1 &&
			    rnap[rnap_index][POSITION] < dwt_stop) {
				if (windex !=
				    curr_dwt_window[rnap_index]) {
				    if (windex > 0) {
				       num_dwt++;
				       dwt[num_dwt] = time -
					last_window_entry[rnap_index];
				    }
				    curr_dwt_window[rnap_index] = windex;
				    last_window_entry[rnap_index] = time;
				}
			}
		}

		return(1);
	}
	else {										
		rnap[rnap[rnap_index][MINUS_ONE]][PLUS_ONE] = 0;	
		promoter[rnap[rnap_index][SOURCE_PROMO]][OLDEST_RNAP] =
		   rnap[rnap_index][MINUS_ONE];
		if (promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP]
		    == rnap_index) {
		     promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP]
		        = 0;
		}

		if (!(reaction_index[RNAP_PAUSE_RXN(rnap_index)] < 0)) {
			Update_reaction_queue(RNAP_PAUSE_RXN(rnap_index),
					      INFINITY);
			Remove_reaction(RNAP_PAUSE_RXN(rnap_index));
		}

		newest_rnap =
		   promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP];
		promo_index = rnap[rnap_index][SOURCE_PROMO];
		if (rnap[newest_rnap][STATE] == CLOSED) {
			newest_rnap = rnap[newest_rnap][PLUS_ONE];
		}
		if (rnap[rnap_index][MINUS_ONE] &&
		    rnap[rnap[rnap_index][MINUS_ONE]][STATE] != CLOSED &&
		    rnap[rnap[rnap_index][MINUS_ONE]][BARRIER_AHEAD] <=0){

		  rnap_sep[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] +=
		     rnap_sep[rnap_index][DOWNSTREAM];
		  if (BUBBLE_ADJ) {
		     rnap_sep[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] += 
			(bubble[UPSTREAM] + bubble[DOWNSTREAM]);
		  }
		  lk[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] +=
		     lk[rnap_index][DOWNSTREAM];
		  Process_supercoiling_changes(rnap[rnap_index][MINUS_ONE],	// Will also adjust newest
		     0,DOWNSTREAM);						// RNAP if upstream and
		  Update_RNAP_rates(time,rnap[rnap_index][MINUS_ONE]);		// departing RNAPs in contact
                  if (rnap[rnap_index][MINUS_ONE] != newest_rnap) {
                        Update_RNAP_rates(time,newest_rnap);
                  }
		}

		lk[rnap_index][UPSTREAM] = 0.0;					// To avoid counting on
		lk[rnap_index][DOWNSTREAM] = 0.0;				// on transient loadings
		rnap_sep[rnap_index][UPSTREAM] = 0;				// without activation
		rnap_sep[rnap_index][DOWNSTREAM] = 0;

		new_old = rnap[rnap_index][MINUS_ONE];
		if (new_old && rnap[new_old][STATE] != CLOSED &&
		    rnap[new_old][BARRIER_AHEAD] <= 0){
			Update_topoisomerase_lk(
					    lk[new_old][DOWNSTREAM],
					    sigma[new_old][DOWNSTREAM],
					    rnap_sep[new_old][DOWNSTREAM],	// Should be able to use sep
					    rnap[new_old][POSITION],		// here: already accounts for
					    stop_barrier,			// add'l length in wraparound
					    MIDDLE, promo_index,
					    rnap_index);
			if (plasmid) {
			   checkster = Define_plasmid_wraparound_region(
					    rnap[new_old][POSITION],
					    promo_index,FORWARD,
					    wrap_ends);
			   Update_topoisomerase_lk(
					    lk[new_old][DOWNSTREAM],
					    sigma[new_old][DOWNSTREAM],
					    rnap_sep[new_old][DOWNSTREAM],
					    wrap_ends[0],wrap_ends[1],
					    MIDDLE, promo_index,0);
			}
		}
		else {
			temp_sigma =
			  (open_lk[promo_index]/relaxed_open_lk) - 1.0;
			Update_topoisomerase_lk(
					    open_lk[promo_index],
					    temp_sigma,
					    stop_barrier - start_barrier,
					    start_barrier,
					    stop_barrier,
					    MIDDLE, promo_index,
					    rnap_index);
		}

		num_tx++;
		tx_rates[num_tx] = ((double) tx_end) / 
				   (time - tx_times[rnap_index]);
		final_pos[num_tx] = tx_end;
		final_start_time[num_tx] = tx_times[rnap_index];
		final_gene_time[num_tx] = gene_time[rnap_index];
		final_stop_time[num_tx] = time;
		final_shut_pos[num_tx] = shutoff_position[rnap_index];
		if (time < promo_kill) {
			final_shut_pos[num_tx] = 99999;
		}
		shutoff_position[rnap_index] = 0;
		
		if (!(reaction_index[RNAP_TERM_RXN(rnap_index)] < 0)) {
			Update_reaction_queue(RNAP_TERM_RXN(rnap_index),
					      INFINITY);
			Remove_reaction(RNAP_TERM_RXN(rnap_index));
		}
		if (!(reaction_index[RNAP_EVAL_RXN(rnap_index)] < 0)) {
			Update_reaction_queue(RNAP_EVAL_RXN(rnap_index),
					      INFINITY);
			Remove_reaction(RNAP_EVAL_RXN(rnap_index));
		}
		if (!(reaction_index[RNAP_SURVEILLANCE_RXN(rnap_index)]
		    < 0)) {
			Update_reaction_queue(
			  RNAP_SURVEILLANCE_RXN(rnap_index),INFINITY);
			Remove_reaction(
			  RNAP_SURVEILLANCE_RXN(rnap_index));
		}
		if (!(reaction_index[RNAP_RELAX_RXN(rnap_index)]
		    < 0)) {
			Update_reaction_queue(
			  RNAP_RELAX_RXN(rnap_index),INFINITY);
			Remove_reaction(
			  RNAP_RELAX_RXN(rnap_index));
		}


		printf("RNAP %d FINISHED TRANSCRIPTION AT TIME %f; "
			"RATE IS: %f\n",
			rnap_index,time,tx_rates[num_tx]);
		printf("%d RNAP RUNS COMPLETE\n",num_tx);

                final_rna_log[rna[
                 rnap[rnap_index][RNAP_RNA]][MASTER_INDEX]].end_tx = time;

		num_transcripts++;						// Deal with completion of transcript:
//		Update_reaction_queue(RNAP_MOVE_RXN(rnap_index),INFINITY);	// RNAP in i-1 position no longer has
//		Remove_reaction(RNAP_MOVE_RXN(rnap_index));			// an i+1 neighbor; the oldest RNAP on
										// the promoter is now RNAP i-1; and the
		rna[rnap[rnap_index][RNAP_RNA]][MATURE] = 1;			// moving of this RNAP needs to be
		memset(rnap[rnap_index],0,					// removed from the reaction pool
		       RNAP_INFO*sizeof(*rnap[rnap_index]));
		rnap_bound--;
		checkster++;
		result++;
		first_free_rnap = rnap_index;					// Eliminate RNAP as object of concern
		return(0);							// on RNA and recycle now-empty RNAP	
	}
}

int Advance_RNAP(double time, int rnap_index, int MAX_PROTEINS,
                 float **proteins_made, int displace)
{
	int result;
	int checkster,wrap_ends[2];
	int first_rnap,last_rnap;
	double next_fire;
	dna_strip[rnap[rnap_index][SOURCE_PROMO]][
		rnap[rnap_index][POSITION]] = 0;
	rnap[rnap_index][POSITION] += displace;
	rna_strip[rnap[rnap_index][RNAP_RNA]][
		rnap[rnap_index][POSITION]] = 0;

	first_rnap= promoter[rnap[rnap_index][SOURCE_PROMO]][NEWEST_RNAP];
	if (rnap[first_rnap][STATE] == CLOSED) {
		first_rnap = rnap[first_rnap][PLUS_ONE];
	}
	last_rnap = promoter[rnap[rnap_index][SOURCE_PROMO]][OLDEST_RNAP];


	if (rnap[rnap_index][POSITION] <= tx_end) {
		dna_strip[rnap[rnap_index][SOURCE_PROMO]][
			rnap[rnap_index][POSITION]] = rnap_index;
		if (rnap[rnap_index][STATE] == PRE_TRANSLOC) {
			rnap[rnap_index][STATE] = POST_TRANSLOC;
			Update_reaction_queue(RNAP_ON_F_RXN(rnap_index),
			   INFINITY);
		   if (!(NO_PAUSE_MODEL)) {
			next_fire = time + Find_firing_time(mu_rnap_add);
			if (no_super_pause[rnap[rnap_index][POSITION]]) {
				next_fire += Find_firing_time(
				 no_super_dur[rnap[rnap_index][POSITION]]);
			}
			Add_reaction_to_queue(RNAP_ADD_RXN(rnap_index),
			   next_fire);
		   }
		   else {							// Immediate add in
			result = Add_RNA_nt(time,rnap_index,MAX_PROTEINS,	// no-pause model
					proteins_made);
		   }
		}
		if (rnap[rnap_index][STATE] == OFF_P2 ||			// P1 immobile and not included
		    rnap[rnap_index][STATE] == OFF_P3) {
			Determine_RNAP_riboload_effect(time,rnap_index,
		   	   rna[rnap[rnap_index][RNAP_RNA]][THREE_END],0);
			Determine_off_pathway_move(time,rnap_index);
			if (displace < 0 &&					// Check for backtracking into
			    !(rnap[rnap_index][RIBO_TETHER]) &&			// ribosome-RNAP tether (omit
			    rnap[rnap_index][CLOSEST_RIBO]) {			// check if tether already
				if (Check_ribosome_tether(rnap_index,0)) {	// in place or no ribosome)
					rnap[rnap_index][RIBO_TETHER] =
					  rnap[rnap_index][CLOSEST_RIBO];

				}
			}
		}	
	

		if (!(QUIET_MODE)) {
			printf("JUST ADVANCED RNAP %d AT TIME %f\n",
				rnap_index,time);
		}

		Process_supercoiling_changes(rnap_index,displace,BOTH);
                if (!(rnap[rnap_index][STATE] == OFF_P2 ||
                      rnap[rnap_index][STATE] == OFF_P3)) {
                        Update_reaction_queue(
                           RNAP_PAUSE_RXN(rnap_index),INFINITY);
                }
		Update_RNAP_rates(time,rnap_index);
		Update_RNAP_rates(time,rnap[rnap_index][MINUS_ONE]);		// Calls now cover BOTH effects
		Update_RNAP_rates(time,rnap[rnap_index][PLUS_ONE]);		// handled in section commented
										// out AND torque-dependent
										// rate changes in +/-1 pos'ns
		if (!(QUIET_MODE)) {
			printf("DONE WITH +/- 1 UPDATES\n");
		}

		Update_topoisomerase_lk(lk[rnap_index][DOWNSTREAM],		// Remove topoisomerases as
					sigma[rnap_index][DOWNSTREAM],		// necessary in intervals within
					rnap_sep[rnap_index][DOWNSTREAM],	// which sigma has changed; all
					rnap[rnap_index][POSITION],		// parameters can be borrowed
					rnap[rnap_index][POSITION] +		// from the RNAP advanced
					rnap_sep[rnap_index][DOWNSTREAM],
					MIDDLE,
					rnap[rnap_index][SOURCE_PROMO],
					0);
		if (plasmid && rnap_index == last_rnap &&			// Overshoot above should not be
		    rnap[rnap_index][BARRIER_AHEAD] <= 0) {			// problematic; topoisomerase
			checkster = Define_plasmid_wraparound_region(		// adjustments in wraparound are
					rnap[rnap_index][POSITION],		// picked up here
					rnap[rnap_index][SOURCE_PROMO],
					FORWARD,wrap_ends);
			Update_topoisomerase_lk(
					lk[rnap_index][DOWNSTREAM],
					sigma[rnap_index][DOWNSTREAM],
					rnap_sep[rnap_index][DOWNSTREAM],
					wrap_ends[0],wrap_ends[1],
					MIDDLE,
					rnap[rnap_index][SOURCE_PROMO],
					0);
			if (rnap[first_rnap][BARRIER_BEHIND] <= 0 &&		// Update newest RNAP if linked
			    rnap[rnap_index][MINUS_ONE] != first_rnap &&	// and not already done
			    rnap_index != first_rnap) {
				Update_RNAP_rates(time,first_rnap);
			}
		}
		wrap_ends[0] = 0;
		wrap_ends[1] = 0;
		Update_topoisomerase_lk(lk[rnap_index][UPSTREAM],
					sigma[rnap_index][UPSTREAM],
					rnap_sep[rnap_index][UPSTREAM],
					rnap[rnap_index][POSITION] -
					rnap_sep[rnap_index][UPSTREAM],
					rnap[rnap_index][POSITION],
					MIDDLE,
					rnap[rnap_index][SOURCE_PROMO],
					0);
		if (plasmid && rnap_index == first_rnap &&
		    rnap[rnap_index][BARRIER_BEHIND] <= 0) {
			checkster = Define_plasmid_wraparound_region(
					rnap[rnap_index][POSITION],
					rnap[rnap_index][SOURCE_PROMO],
					REVERSE,wrap_ends);
			Update_topoisomerase_lk(
					lk[rnap_index][UPSTREAM],
					sigma[rnap_index][UPSTREAM],
					rnap_sep[rnap_index][UPSTREAM],
					wrap_ends[0],wrap_ends[1],
					MIDDLE,
					rnap[rnap_index][SOURCE_PROMO],
					0);
			if (rnap[last_rnap][BARRIER_AHEAD] <= 0 &&		// Similarly, update oldest
			    rnap[rnap_index][PLUS_ONE] != last_rnap &&		// RNAP if necessary
			    rnap_index != last_rnap) {
				Update_RNAP_rates(time,last_rnap);
			}
		}

		if (MAKE_KYMOGRAPH &&
		    !(record_on) &&
                    time > kymo_start &&
                    sim_endtime - time > kymo_duration &&
                    rnap[rnap_index][POSITION] >= kymo_min_pos &&
                    rnap[rnap_index][POSITION] <= kymo_max_pos) {
                        record_on++;
                        Add_kymograph_reactions(time);
                }

		return(1);
	}
	else {
		Update_reaction_queue(RNAP_ON_F_RXN(rnap_index),INFINITY);	// Should always be on pathway:
		Remove_reaction(RNAP_ON_F_RXN(rnap_index));			// Immediately perform pseudo-add
		result = Add_RNA_nt(time,rnap_index,MAX_PROTEINS,		// at transcript end
		   proteins_made);
		checkster++;
		result++;
		return(0);
	}
}

int Process_blocked_RNAP(double time, int rnap_index)
{
	int result = 0;								// Below no longer needed b/c
	result = Process_pretranslocated_state(time,rnap_index);		// of autocheck on processing
/*
	if (rnap[rnap_index][MINUS_ONE] != 0 &&
	    rnap[rnap[rnap_index][MINUS_ONE]][STATE] <= EFFECTOR_STATE &&
	    Check_for_RNAP_RNAP_contacts(rnap_index,MINUS_ONE) == 1) {
		result += (Process_pretranslocated_state(time,rnap_index,
		   TRAILING_RNAP));
	}
	if (Check_for_RNAP_ribosome_contacts(rnap_index,TETHER_LENGTH)
	    == 1) {								// Note that this MUST be
		result += (Process_pretranslocated_state(time,rnap_index,	// called (no ribo check if
		   TRAILING_RIBO));						// forward above); note that if
	}									// BOTH ribosomes and RNAPs are
	if (result == 0) {							// effectors, this check is
		result = Process_pretranslocated_state(time,rnap_index,		// potentially redundant; it's
		   INDEPENDENT);						// preserved for indep. testing
	}									// NB: The commented out option
*/
/*										// below MAY be required to
	if (result == 0) {							// prevent cascading pauses
		result = Process_pretranslocated_state(time,rnap_index,		// behind the lead RNAP....
		   AUTO_FORWARD);
	}
*/
	return(result);
}

int Attempt_RNAP_move(double time, int rnap_index,
                      float rnap_push_prob, int MAX_PROTEINS,
                      float **proteins_made, int direction)
{
	int result=0,new_index,loc,rna_index;
	int check,rnap_check,ribo_check=0,rna_check=0,displace;
	int topo_check=0;
	int tether_check=1;
	int topo_id;
	float  push;
	loc = rnap[rnap_index][POSITION];
	rna_index = rnap[rnap_index][RNAP_RNA];
	if (direction == FORWARD) {
		check = PLUS_ONE;
		displace = 1;
	}
	else {
		check = MINUS_ONE;
		displace = -1;
	}
	rnap_check = (rnap[rnap_index][check] == 0 ||
	   Check_for_RNAP_RNAP_contacts(rnap_index,check) == 0);
	if (rnap_check) {
		if (rnap[rnap_index][RIBO_TETHER]) {				// Assuming for now that
//		if (rnap[rnap_index][RIBO_TETHER] &&				// off-pathway RNAPs may
//		    rnap[rnap_index][STATE]  == PRE_TRANSLOC) {			// also be tethered
			tether_check=Check_ribosome_tether(rnap_index,1);

		}
		ribo_check = ((direction == FORWARD && tether_check) ||
		   (direction == REVERSE && 
		    Check_for_RNAP_ribosome_contacts(rnap_index,1)		// Check for adequate room to
		     == 0));							// b.t. 1 nt relative to ribo
		if (ribo_check) {
			if (direction == FORWARD) {				// Currently preventing forward
				if (rnap[rnap_index][STATE] ==			// tracking beyond RNA 3' end...
				    PRE_TRANSLOC ||
				    rna[rna_index][THREE_END] >=
				    loc + displace) {
					rna_check = 1;
				}
			}
			else {
				if ((loc + displace) -				// ...and requiring a complete
				    rna[rna_index][FIVE_END] > 14 &&		// bubble while enforcing a
				    rna[rna_index][THREE_END] -			// maximum backtrack depth in
				    (loc + displace) <= MAX_BACKTRACK) {	// reverse movement
					rna_check = 1;
				}
			}
		}
		if (ribo_check && rna_check) {
			topo_id = rnap[rnap_index][TOPO_DOWNSTREAM];
			if (displace < 0) {
				topo_id = rnap[rnap_index][TOPO_UPSTREAM];
			}
			topo_check = (!(topo_id) || abs((loc + displace) -
			   topo[topo_id][POSITION]) >=
			   topo_rnap_width[topo[topo_id][TOPO_TYPE]]);
			if (!(topo_check)) {
				if (Sample_uniform_distribution() <
				    p_rnap_topo_eject[
				      topo[topo_id][TOPO_TYPE]]) {

					printf("RNAP %d ON PROMOTER %d "
						"ATTEMPTING TO MOVE TO %d "
						"REMOVES TOPOISOMERASE %d "
						"AT %d\n", rnap_index,
						rnap[
						  rnap_index][SOURCE_PROMO],
						rnap[rnap_index][POSITION]
						 + displace, topo_id,
						topo[topo_id][POSITION]);

					Unbind_topoisomerase_on_schedule(
					   time,topo_id);
					topo_check = 1;
				}
				else {

					printf("RNAP %d ON PROMOTER %d "
						"ATTEMPTING TO MOVE TO %d "
						"BLOCKED BY TOPOISOMERASE %d AT %d\n",
						rnap_index,rnap[rnap_index][SOURCE_PROMO],
						rnap[rnap_index][POSITION] + displace,
						topo_id,topo[topo_id][POSITION]);

				}
			}

		}


	}
	if (rnap_check && ribo_check && rna_check && topo_check) {
		result = Advance_RNAP(time,rnap_index,MAX_PROTEINS,
			   proteins_made,displace);
		return(1);
	}

	else if (rnap[rnap_index][STATE] == PRE_TRANSLOC) {			// Recursive push attempts on downstream
		push = Sample_uniform_distribution();				// RNAPs(low prob if many stacked)
		if (push <= rnap_push_prob) {					// This is vestigial stuff changed only
			new_index = rnap[rnap_index][PLUS_ONE];			// to be compatible with new rates/fns;
			result = Attempt_RNAP_move(time,new_index,		// it will be overhauled to be used....
				    rnap_push_prob,MAX_PROTEINS,
				    proteins_made,direction);
			if (rnap[rnap_index][PLUS_ONE] == 0 ||
			    abs((rnap[rnap_index][POSITION]+1) -
			    rnap[rnap[rnap_index][PLUS_ONE]][POSITION])
			    >= pol_width) {
				result = Advance_RNAP(time,rnap_index,
					    MAX_PROTEINS,proteins_made,
					    displace);
			}
			else {
				result = Process_blocked_RNAP(
					   time,rnap_index);
			}
		}
		else {								// Assign new rxn time if blocked
			result = Process_blocked_RNAP(time,rnap_index);
		}
	}
	else {
		Determine_off_pathway_move(time,rnap_index);
	}
	return(result);
}

