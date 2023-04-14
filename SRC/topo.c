#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define DEFINE_TOPO_VARIABLES
#include "INCL/topo.h"

#include "INCL/sim_types.h"
#include "INCL/gen_def.h"
#include "INCL/sim_settings.h"
#include "INCL/reaction_manager.h"
#include "INCL/replication.h"
#include "INCL/regulation.h"
#include "INCL/transcription.h"
#include "INCL/translation.h"
#include "INCL/supercoiling.h"
#include "INCL/utilities.h"


void Unbind_topoisomerase(int topo_index)
{										// Fn covers unbinding events
										// due to either scheduled
	Update_reaction_queue(TOPO_UNBINDING_RXN(topo_index),			// departure or changes in
			      INFINITY);					// local lk; the former does
	Remove_reaction(TOPO_UNBINDING_RXN(topo_index));			// not require separate
	if (!(reaction_index[TOPO_CLEAVAGE_RXN(topo_index)] < 0)) {		// invocation of the update fn
		Update_reaction_queue(TOPO_CLEAVAGE_RXN(topo_index),
				      INFINITY);
		Remove_reaction(TOPO_CLEAVAGE_RXN(topo_index));
	}
	if (!(reaction_index[TOPO_SWITCH_STATE_RXN(topo_index)] < 0)) {
		Update_reaction_queue(TOPO_SWITCH_STATE_RXN(topo_index),
				      INFINITY);
		Remove_reaction(TOPO_SWITCH_STATE_RXN(topo_index));
	}

	if (!(QUIET_MODE)) {
		printf("TOPOISOMERASE %d (TYPE %d) LEAVING PROMOTER %d "
			"AT POSITION %d: IT WAS CLOSEST DOWNSTREAM TOPO "
			"FOR RNAP %d AND CLOSEST UPSTREAM TOPO FOR %d; "
			"PLUS/MINUS ONE RNAP ARE %d/%d\n",
			topo_index, topo[topo_index][TOPO_TYPE],
			topo[topo_index][SOURCE_PROMO], 
			topo[topo_index][POSITION], 
			topo[topo_index][CLOSEST_DOWN],
			topo[topo_index][CLOSEST_UP],
			topo[topo_index][PLUS_ONE],
			topo[topo_index][MINUS_ONE]);
	}

	if (topo[topo_index][CLOSEST_UP]) {
		num_upstream_update++;
		upstream_list[num_upstream_update] =
		  topo[topo_index][PLUS_ONE];
		rnap[topo[topo_index][PLUS_ONE]][TOPO_UPSTREAM] = 0;
	}
	if (topo[topo_index][CLOSEST_DOWN]) {
		num_downstream_update++;
		downstream_list[num_downstream_update] =
		  topo[topo_index][MINUS_ONE];
		rnap[topo[topo_index][MINUS_ONE]][TOPO_DOWNSTREAM] = 0;
	}
	if (topo[topo_index][PROMO_BLOCKER]) {					// Update availability of 
		promo_topo_blocked[topo[topo_index][SOURCE_PROMO]]--;		// promoter using
	}									// association stored in topo

	memset(topo[topo_index],0,
	       TOPO_INFO*sizeof(*topo[topo_index]));
	first_free_topo = Find_first_free(topo_search);
	return;
}


	
int Find_closest_topoisomerase(int rnap_index, int side)
{
	int i,promo_index,side_factor,adj_rnap;
	int check, which = 0, dist, curr_min =99999;
	promo_index = rnap[rnap_index][SOURCE_PROMO];
	side_factor = (side == DOWNSTREAM ? 1 : -1);
	adj_rnap = rnap[rnap_index][MINUS_ONE + (side == DOWNSTREAM)];
	if (adj_rnap) {								// Topoisomerase only relevant for
		curr_min = abs(rnap[adj_rnap][POSITION] -			// movement if closer than closest
			       rnap[rnap_index][POSITION]);			// (up)downstream RNAP
	}
	for (i = 1; i <= num_topo[promo_index]; i++) {
		check = topo_master_list[promo_index]
					[curr_topo_list[promo_index]][i];
		dist = (topo[check][POSITION] -
			rnap[rnap_index][POSITION]) * side_factor;
		if (dist > 0 && dist < curr_min) {
			curr_min = dist;
			which = check;
		}
	}
	return(which);
}
			
		

int Check_upstream(int left_check, int right_check,
				 int check_site)
{
	return(check_site < left_check);
}

int Check_downstream(int left_check, int right_check,
				 int check_site)
{
	return(check_site > right_check);
}

int Check_middle(int left_check, int right_check,
				 int check_site)
{
	return(check_site > left_check && check_site < right_check);

}

void Update_topoisomerase_RNAP_connections()
{
	int i,which;
	for (i = 1; i <= num_upstream_update; i++) {
		which = Find_closest_topoisomerase(upstream_list[i],
						   UPSTREAM);
		if (which) {
			topo[which][CLOSEST_UP] = upstream_list[i];
			rnap[upstream_list[i]][TOPO_UPSTREAM] = which;
		}
	}
	for (i = 1; i <= num_downstream_update; i++) {
		which = Find_closest_topoisomerase(downstream_list[i],
						   DOWNSTREAM);
		if (which) {
			topo[which][CLOSEST_DOWN] = downstream_list[i];
			rnap[downstream_list[i]][TOPO_DOWNSTREAM] = which;
		}
	}
	return;
}
		
void Unbind_topoisomerase_on_schedule(double time, int topo_index)
{
	int i, promo_index, num_topo_temp = 0, next_slot;
	int lost_starter_up=0,lost_starter_down=0;
	int test_dist,testmin=99999;

	promo_index = topo[topo_index][SOURCE_PROMO];
	next_slot = (curr_topo_list[promo_index] + 1) % 2;
	num_upstream_update = 0;
	num_downstream_update = 0;
	if (topo_index == starter_up[promo_index]) {
		lost_starter_up++;
		starter_up[promo_index] = 0;
	}
	else if (topo_index == starter_down[promo_index]) {
		lost_starter_down++;
		starter_down[promo_index] = 0;
	}
	for (i = 1; i <= num_topo[promo_index]; i++) {
		if (topo_finder[promo_index][i] == topo_index) {
			continue;
		}
		num_topo_temp++;
		topo_master_list[promo_index][next_slot][num_topo_temp] =
		   topo_finder[promo_index][i];
		if (lost_starter_down) {
			test_dist = topo[topo_finder[promo_index][i]]
					   [POSITION] - 1;
			if (test_dist >= 0 && test_dist < testmin) {
				testmin = test_dist;
				starter_down[promo_index] =
				   topo_finder[promo_index][i];
			}
		}
		else if (lost_starter_up) {
			test_dist = 1 - 
				    topo[topo_finder[promo_index][i]]
					   [POSITION];
			if (test_dist >= 1 && test_dist < testmin) {
				testmin = test_dist;
				starter_up[promo_index] =
				   topo_finder[promo_index][i];
			}
		}
			
	}
	Unbind_topoisomerase(topo_index);
	curr_topo_list[promo_index] = next_slot;
	num_topo[promo_index] = num_topo_temp;
	topo_finder[promo_index] =
	  &topo_master_list[promo_index][next_slot][0];
	Update_topoisomerase_RNAP_connections();
	return;
}
	

void Update_topoisomerase_lk(float new_lk, float new_sigma, int new_len,
			     int left_check, int right_check,
			     int loc_check, int promo_index,
			     int rnap_done)
{
	int i,which;
	int site,in_range,next_slot;
	int num_topo_temp = 0, min_dist_up=99999,min_dist_down=99999;
	int temp_up,temp_down;
	float max_topo_lk;
	int (*check_fn[])(int, int, int) = {Check_upstream,Check_downstream,
					    NULL, Check_middle};


	if (!(QUIET_MODE)) {
		printf("IN UPDATE TOPOISOMERASE LK FN! CHECKING %d "
			"INTERVAL FROM %d TO %d ON PROMOTER %d\n",
			loc_check,left_check,right_check,promo_index);
	}

	max_topo_lk = relaxed_lk[new_len] - 1.0;
	next_slot = (curr_topo_list[promo_index] + 1) % 2;
	num_upstream_update = 0;
	num_downstream_update = 0;
	starter_up[promo_index] = 0;
	starter_down[promo_index] = 0;
	for (i = 1; i <= num_topo[promo_index]; i++) {
		which = topo_finder[promo_index][i];
		site = topo[which][POSITION];
		temp_down = site - 1;
		temp_up = 1 - site;
		if (rnap_done) {						// Processing adjustments on end
			if (topo[which][MINUS_ONE] == rnap_done) {
				if (rnap[rnap[rnap_done][MINUS_ONE]][STATE]
				    != CLOSED) {
					topo[which][MINUS_ONE] =
					  rnap[rnap_done][MINUS_ONE];
				}
				else {
					topo[which][MINUS_ONE] = 0;
				}
			}
			if (topo[which][PLUS_ONE] == rnap_done) {
				topo[which][PLUS_ONE] =
				  rnap[rnap_done][PLUS_ONE];
			} 
		}
		in_range =
		  (*check_fn[loc_check])(left_check,right_check,site);
		if (in_range && topo[which][TOPO_TYPE] == TOPO_IA &&
		    new_lk > max_topo_lk) {

			if (!(QUIET_MODE)) {
				printf("TOPO IA NUMBER %d AT LOCATION "
					"%d IS BEING REMOVED!\n",
					which,topo[which][POSITION]);
			}	

			Unbind_topoisomerase(which);
			continue;
		}
		num_topo_temp++;
		topo_master_list[promo_index][next_slot][num_topo_temp] =
		  which;
		if (temp_down >= 0 && temp_down < min_dist_down) {
			min_dist_down = temp_down;
			starter_down[promo_index] = which;
		}
		if (temp_up >= 1 && temp_up < min_dist_up) {
			min_dist_up = temp_up;
			starter_up[promo_index] = which;
		}
	}
	num_topo[promo_index] = num_topo_temp;
	curr_topo_list[promo_index] = next_slot;
	topo_finder[promo_index] =
	  &topo_master_list[promo_index][next_slot][0];

	Update_topoisomerase_RNAP_connections();
}
		
void Load_burst_distribution()
{
	int i;
	float interburst_time;
	float mean_burst_time,mean_intraburst_rate;
	num_burst_rate_distr = 20;
	burst_rate_distr[0] = -1.00;						// Fixed: from distrib
	burst_rate_distr[1] = 0.084;						// in Ashley_Osherhoff, NAR 2017
	burst_rate_distr[2] = 0.139;						// Fig. S2
	burst_rate_distr[3] = 0.178;
	burst_rate_distr[4] = 0.207;
	burst_rate_distr[5] = 0.217;
	burst_rate_distr[6] = 0.243;
	burst_rate_distr[7] = 0.260;
	burst_rate_distr[8] = 0.284;
	burst_rate_distr[9] = 0.321;
	burst_rate_distr[10] = 0.354;
	burst_rate_distr[11] = 0.378;
	burst_rate_distr[12] = 0.405;
	burst_rate_distr[13] = 0.421;
	burst_rate_distr[14] = 0.450;
	burst_rate_distr[15] = 0.476;
	burst_rate_distr[16] = 0.515;
	burst_rate_distr[17] = 0.618;
	burst_rate_distr[18] = 0.691;
	burst_rate_distr[19] = 0.777;
	burst_rate_distr[20] = 1.001;
	mean_burst_size = 6.2;
	mean_intraburst_rate = 107.0;						// From Ashley et al.
	g_alpha = 2.1;
	g_beta = 114.0;
	mean_burst_time = mean_burst_size/mean_intraburst_rate;
	mu_topo_interburst[TOPO_IA] = 9999.9;					// Rate effectively zero
	mu_topo_interburst[GYRASE]  = 2.0;
	mean_topo_run_length[TOPO_IA] = 20.0;
	mean_topo_run_length[GYRASE] = 6.2;
	mu_topoIA_intraburst = 1.0/3.3;			
//	mean_lagtime_topoIA = 5.0;						// From Terekhova,Marko,Mondragon
	for (i = 1; i <= num_burst_rate_distr; i++) {
		burst_set[i] = 50.0/((float) i);
		interburst_time = burst_set[i];
		burst_set[i] = (50.0 - burst_set[i])/106.0;
		interburst_time -= burst_set[i];
		burst_set[i] /= mean_burst_time;
		if (burst_set[i] > 0.01) {
			burst_set[i] = interburst_time/burst_set[i];
		}
		else {
			burst_set[i] = 50.0;
		}

		if (!(QUIET_MODE)) {
			printf("BURST SET[%d]: %f\n",i,burst_set[i]);
		}

	}
	mu_topo_unbind[GYRASE] = 2.5;						// From Stracy,Sherratt,Zawadki,NAR 2019
	mu_topo_unbind[TOPO_IA] = 1/0.043;					// From Tiwari,...,Tse-Dinh,Darici,
										// Biochem Biophys Res Comm 2014
	return;
}

void Determine_mean_lagtime(int topo_index)
{
	int i;
	float draw;
	draw = gsl_rng_uniform(r);
	i = num_burst_rate_distr;
	while (draw < burst_rate_distr[i]) {
		i--;
	}
	mu_topo_lagtime[topo_index] = burst_set[i + 1];
	return;
}
			
void Switch_burst_state(double time, int topo_index)
{

	if (!(QUIET_MODE)) {
		printf("IN SWITCH FUNCTION AT TIME %f\n",time);
	}

	int topo_type;
	float run_length,next_switch=-999.9;
	double next_fire;
	topo[topo_index][BURST_STATE] =
	   (topo[topo_index][BURST_STATE] + 1) % 2;
	topo_type = topo[topo_index][TOPO_TYPE];
	if (topo[topo_index][BURST_STATE] == ON) {
		if (!(reaction_index[TOPO_SWITCH_STATE_RXN(topo_index)]		// Post-burst switch
		    < 0)) {							// time determined
			Update_reaction_queue(					// by run length
			   TOPO_SWITCH_STATE_RXN(topo_index),INFINITY);
			Remove_reaction(
			   TOPO_SWITCH_STATE_RXN(topo_index));
		}
		run_length =
		   Find_firing_time(mean_topo_run_length[topo_type]);
		if (topo_type == TOPO_IA) {
			mu_topo_intraburst[topo_index] = 
			   mu_topoIA_intraburst;
		}
		else {
			run_length /= 2.0;					// B/c 2 coils/cycle;
			mu_topo_intraburst[topo_index] =			// Mu = 1/rate, where
			   2.0 * gsl_ran_gamma(r,g_alpha,1.0/g_beta);		// rates follow inverse
		}								// gamma distribution
		next_fire = time +
		   Find_firing_time(mu_topo_intraburst[topo_index]);
		topo[topo_index][CYCLES_LEFT] = (int) (run_length+0.5);
	}
	else {
		topo[topo_index][CYCLES_LEFT] = 0;
		next_fire = time +
		   Find_firing_time(mu_topo_interburst[topo_type]);
		if (topo_type == GYRASE) {
			Determine_mean_lagtime(topo_index);
		}
		else {
			mu_topo_lagtime[topo_index] = mean_lagtime_topoIA;
		}

		if (!(QUIET_MODE)) {
			printf("MU LAGTIME IS %f\n",
				mu_topo_lagtime[topo_index]);
		}

		next_switch = time +
		   Find_firing_time(mu_topo_lagtime[topo_index]);		// Determined at loading
		Add_reaction_to_queue(TOPO_SWITCH_STATE_RXN(topo_index),	// for topo; above for
		   next_switch);						// gyrase with draw
	}
	if (!(reaction_index[TOPO_CLEAVAGE_RXN(topo_index)] < 0)) {
		Update_reaction_queue(TOPO_CLEAVAGE_RXN(topo_index),
				      next_fire);
	}
	else {
		Add_reaction_to_queue(TOPO_CLEAVAGE_RXN(topo_index),
				      next_fire);
	}

	if (!(QUIET_MODE)) {
		printf("TOPOISOMERASE %d: TYPE: %d; BURST STATE: %d; "
			"NEXT CLEAVAGE RXN: %f; STATE CHANGE TIME %f\n",
			topo_index,topo[topo_index][TOPO_TYPE],
			topo[topo_index][BURST_STATE],next_fire,
			next_switch);
	}

	return;
}
	
int Define_plasmid_wraparound_region(int checksite, int promo_index,
				     int dir_ext, int *wrap_ends)
{
	int check_rnap;								// Fn creates a second interval
	if (dir_ext == FORWARD) {						// to be checked for topo binding
		wrap_ends[0] = start_barrier;
		check_rnap = promoter[promo_index][NEWEST_RNAP];
		if (rnap[check_rnap][STATE] == CLOSED) {			// Will only be accessed if
			check_rnap = rnap[check_rnap][PLUS_ONE];		// an active RNAP is in place
		}
		if (!(check_rnap) ||					
		    rnap[check_rnap][BARRIER_BEHIND] > 0) {
			wrap_ends[1] = start_barrier;
			return(0);						// Reports identity (absence of) 
		}								// connected RNAP for calls where
		wrap_ends[1] = rnap[check_rnap][POSITION];			// this is needed
		return(check_rnap);
	}
	if (dir_ext == REVERSE) {
		wrap_ends[1] = stop_barrier;
		check_rnap = promoter[promo_index][OLDEST_RNAP];
		if (!(check_rnap) ||
		    rnap[check_rnap][STATE] == CLOSED ||
		    rnap[check_rnap][BARRIER_AHEAD] > 0) {
			wrap_ends[0] = stop_barrier;
			return(0);
		}
		wrap_ends[0] = rnap[check_rnap][POSITION];
		return(check_rnap);
	}
	return(0);
}
			

void Perform_topoisomerase_cycle(double time, int topo_index)
{

	if (!(QUIET_MODE)) {
		printf("IN CLEAVAGE RXN FOR TOPO %d AT TIME %f\n"
			,topo_index, time);
	}

	int site,promo_index;
	int left_bracket,right_bracket,left_dist,right_dist;
	int rnap_up,rnap_down,new_len=0;
	int plasmid_linked_rnap = 0;
	int wrap_ends[2] = {0};
	float old_lk,new_lk,new_sigma,check_sigma=0.0;
	float *adj_loc;
	double next_fire;
	site = topo[topo_index][POSITION];
	promo_index = topo[topo_index][SOURCE_PROMO];

	left_bracket = start_barrier;
	right_bracket = stop_barrier;
	adj_loc = &open_lk[promo_index];
	left_dist  =  site - left_bracket;
	right_dist =  right_bracket - site;
	rnap_up    =  topo[topo_index][MINUS_ONE];
	rnap_up    *= (site - rnap[rnap_up][POSITION] < left_dist);
	rnap_down  =  topo[topo_index][PLUS_ONE];
	rnap_down  *= (rnap[rnap_down][POSITION] - site < right_dist);
	if (rnap_up) {
		left_bracket = rnap[rnap_up][POSITION];
		adj_loc = &lk[rnap_up][DOWNSTREAM];
		new_len = rnap_sep[rnap_up][DOWNSTREAM];
		check_sigma = sigma[rnap_up][DOWNSTREAM];
	}
	if (rnap_down) {
		right_bracket = rnap[rnap_down][POSITION];
		adj_loc = &lk[rnap_down][UPSTREAM];
		new_len = rnap_sep[rnap_down][UPSTREAM];
		check_sigma = sigma[rnap_down][UPSTREAM];
	}

	if (plasmid &&								// Account for wraparound
	    ((left_bracket == start_barrier &&					// in plasmid if necessary
	      right_bracket != stop_barrier) ||
	     (right_bracket == stop_barrier &&
	      left_bracket != start_barrier))) {
		if (right_bracket == stop_barrier) {				// Extend search for linked
			plasmid_linked_rnap = 					// RNAP beyond stop...
			   Define_plasmid_wraparound_region(site,
			      promo_index,FORWARD,wrap_ends);
		}
		if (left_bracket == start_barrier) {				// ...and start sites of 
			plasmid_linked_rnap =					// chromosomal equivalent
			   Define_plasmid_wraparound_region(site,
			      promo_index,REVERSE,wrap_ends);
		}
	}
			
	new_lk = *adj_loc;

	if (!(rnap_up || rnap_down)) {
		if (plasmid_linked_rnap) {					// If connected, get 
			if (right_bracket == stop_barrier) {			// relevant linking number
			       new_len =					// from wraparound RNAP
				rnap_sep[plasmid_linked_rnap][UPSTREAM];
			       adj_loc =
				&lk[plasmid_linked_rnap][UPSTREAM];
			       check_sigma =
				sigma[plasmid_linked_rnap][UPSTREAM];
			}
			if (left_bracket == start_barrier) {
			       new_len =
				rnap_sep[plasmid_linked_rnap][DOWNSTREAM];
			       adj_loc =
				&lk[plasmid_linked_rnap][DOWNSTREAM];
			       check_sigma =
				sigma[plasmid_linked_rnap][DOWNSTREAM];
			}
		}
		else {
			new_len = right_bracket - left_bracket;
			check_sigma = (new_lk/relaxed_lk[new_len]) - 1.0;
		}
	}

	new_lk = *adj_loc;

	if (!(QUIET_MODE)) {
		printf("TOPO TYPE IS %d; LOCATION IS %d; CHECK SIGMA IS: "
			"%f; LK IS %f\n",
			topo[topo_index][TOPO_TYPE],
			topo[topo_index][POSITION], check_sigma,new_lk);
	}

	old_lk = new_lk;
	if (plasmid_linearized && (left_bracket == start_barrier ||		// No change to Lk in free
	    right_bracket == stop_barrier)) {					// ends of linearized plasmid
	}
	else if (topo[topo_index][TOPO_TYPE] == TOPO_IA) {			// Check for introducing
		new_lk += 1.0;							// (+) supercoiling not req'd:
	}									// will have been removed if
	else {									// complete cycle not possible
		if (check_sigma >= gyrase_min_sigma &&
		    new_lk > 2.0) {
			if (check_sigma > 0) {
				new_lk -= 2.0;
			}
			else if (gsl_rng_uniform(r) < 0.1) {

			if (!(QUIET_MODE)) {
				printf("UNLIKELIER INTRODUCTION OF (-) "
					"SUPERCOILS! CONGRATULATIONS!\n");
			}

				new_lk -= 2.0;
			}
		}
	}
	if (new_lk != old_lk) {

		if (!(QUIET_MODE)) {
			printf("TOPO CYCLE ALLOWED: CHANGING LK FROM %f "
				"TO %f\n",
				old_lk,new_lk);
		}

		new_sigma = 0.0;
		if (!(wrap_ends[0] == 0 && wrap_ends[1] == 0)) {		// In this case, wraparound is
			if (!(rnap_up || rnap_down) &&				// relevant but length cannot
			    !(plasmid_linked_rnap)) {				// be obtained from a reference
				new_len += (wrap_ends[1] - wrap_ends[0]);	// RNAP: supplement using
			}							// interval created above
			Update_topoisomerase_lk(new_lk,new_sigma,new_len,
					wrap_ends[0],wrap_ends[1],MIDDLE,
					promo_index,0);
		}
		new_sigma = (new_lk/relaxed_lk[new_len]) - 1.0;			// Can be dropped
		Update_topoisomerase_lk(new_lk,new_sigma,new_len,
					left_bracket,right_bracket,MIDDLE,
					promo_index,0);

		if (!(plasmid_linearized)) {					// If ring or chromosomal DNA,
			open_lk[promo_index] += (new_lk - old_lk);		// open_lk always updated b/c it may
		}								// be used in partitioning Lk
										// in the absence of RNAPs

		*adj_loc = new_lk;						// Note that change to the pointer is
										// redundant on an open operon

		if (rnap_down) {						// Must check first b/c of order above
		     Process_supercoiling_changes(rnap_down,0,UPSTREAM);	
		}								// *SHOULD* be OK because lk will be
		else if (rnap_up) {						// changed in the fn called for the
		     Process_supercoiling_changes(rnap_up,0,DOWNSTREAM);	// first/last RNAP if connected
		}
		else if (plasmid_linked_rnap) {					// Make separate call if not changed
		     Process_supercoiling_changes(plasmid_linked_rnap,		// with up/down calls
						  0,BOTH);
		}
		Update_RNAP_rates(time,rnap_up);
		Update_RNAP_rates(time,rnap_down);
		if (plasmid_linked_rnap) {					// Updates rates if cleavage between
		     Update_RNAP_rates(time,plasmid_linked_rnap);		// topologically-linked first and last
		}								// RNAPs {will not be incl. in up/down)

		if (plasmid) {							// Account for any changes so that
			lk_barrier[promo_index][UPSTREAM] =			// in plasmid up- and downstream may
			   lk_barrier[promo_index][DOWNSTREAM];			// be used interchangeably
		}
		if (!(topo[topo_index][SOURCE_PROMO])) {			// Covers possibility
			return;							// that topoisomerase is
		}								// dislodged post-cleavage
	}
	else {

		if (!(QUIET_MODE)) {
			printf("TOPO CYCLE VERBOTEN!\n");
		}

	}

	if (topo[topo_index][BURST_STATE] == ON) {
		topo[topo_index][CYCLES_LEFT]--;

		if (!(QUIET_MODE)) {
			printf("TOPO IN BURST STATE: CYCLES NOW "
				"REMAINING: %d\n",
				topo[topo_index][CYCLES_LEFT]);
		}

		if (topo[topo_index][CYCLES_LEFT] > 0) {
			next_fire = time +
			 Find_firing_time(mu_topo_intraburst[topo_index]);
		}
		else {
			Switch_burst_state(time,topo_index);
			return;
		}
	}
	else {
		next_fire = time +
			    Find_firing_time(mu_topo_interburst[
				topo[topo_index][TOPO_TYPE]]);
	}	
	Update_reaction_queue(TOPO_CLEAVAGE_RXN(topo_index), next_fire);	// Reaction already in queue
	return;
}
	

int Attempt_topoisomerase_binding(double time, int topo_type,
				  int topo_index, int min_spacing)
{
	int i,check,which_promo,site;
	int barrier_check,up_check,down_check;
	int temp_plus=0,temp_minus=0;
	int spacing = 0,test_len;
	int try_site,check_site,check_type;
	int rnap_position,adj_topo;
	int first_rnap,last_rnap;
	double next_fire;
	float test_lk,mu_topo_exit;
	next_fire = time + Find_firing_time(mu_topo_binding[topo_type]);
	if (topo_type == TOPO_IA) {
		Update_reaction_queue(TOPO_I_BINDING_RXN(1),next_fire);

		if (!(QUIET_MODE)) {
			printf("ATTEMPTING TOPO I RXN AT TIME %f\n",time);
		}

	}
	if (topo_type == GYRASE) {
		Update_reaction_queue(GYRASE_BINDING_RXN(1),next_fire);

		if (!(QUIET_MODE)) {
			printf("ATTEMPTING GYRASE RXN AT TIME %f\n",time);
		}

	}
	site = gsl_rng_uniform_int(r,num_promoters*dna_length);
	which_promo = (site / dna_length);
	which_promo++;

	if (!(QUIET_MODE)) {
		printf("DRAW IS %d; PROMO IS %d\n",site,which_promo);
	}

	site = site % dna_length;
	site++;

	if (!(QUIET_MODE)) {
		printf("SITE ADJUSTED TO %d\n",site);
	}

	site += start_barrier;

	if (!(QUIET_MODE)) {
		printf("SITE ADJUSTED FURTHER TO %d\n",site);
	}

	first_rnap = promoter[which_promo][NEWEST_RNAP];
	if (rnap[first_rnap][STATE] == CLOSED) {
		first_rnap = rnap[first_rnap][PLUS_ONE];
	}
	last_rnap = promoter[which_promo][OLDEST_RNAP];

	if (topo_type == TOPO_IA &&
	    abs(1 - site) < topo_rnap_width[TOPO_IA] + 2) { 
		return(0);
	}
	
	if (guided_topo_binding && topo_type == TOPO_IA && site > -150) {
		return(0);
	}
	if (guided_gyrase_binding && topo_type == GYRASE &&			// REVISE THESE FOR PLASMID
	    site < tx_end + 150) {						// WILL NEED TO BE CAREFUL
		return(0);							// TO MATCH SIZE OF "UPSTREAM"
	}									// REGION IN PLASMID WITH
	if (sep_topoIA_interrnap_rate && topo_type == TOPO_IA &&		// EQUIVALENT IN CHROMOSOME
	    site > -150) {							// (EVEN THOUGH GYRASE/TOPO
										// EFFECTS SHOULD CANCEL...)
		if (site >= rnap[first_rnap][POSITION] &&
		    site <= rnap[last_rnap][POSITION]) {
			if (gsl_rng_uniform(r) >
			    relative_topoIA_ds_binding) {
				return(0);
			}
		}
		else {								// Reduced binding rates for
			if (gsl_rng_uniform(r) >				// sites neither in privileged
			    relative_topoIA_ds_binding *			// upstream region nor 
			    rel_extra_topo) {					// associated with (between)
				return(0);					// transcribing RNAPs
			}
		}
	}
		


	if (abs(site - rnap[promoter[which_promo][NEWEST_RNAP]][POSITION])	// Quick check of attempt
	    < topo_rnap_width[topo_type]) {					// against position of any
		return(0);							// closed-conformation RNAP
	}									// (only open used below)

	if (site < 1) {
		if (promoter[which_promo][NEWEST_RNAP] &&
		    rnap[promoter[which_promo][NEWEST_RNAP]][STATE] !=
		    CLOSED) {
			temp_plus =
			   promoter[which_promo][NEWEST_RNAP];
		}
		else {								// This should give the proper
			temp_plus =						// right boundary whether or
			   rnap[promoter[which_promo][NEWEST_RNAP]]		// not an RNAP is present
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

	if (!(QUIET_MODE)) {
		printf("INITIAL TOPO BEHIND %d AT %d; AHEAD: %d AT %d\n",
			temp_minus,rnap[temp_minus][POSITION],temp_plus,
			rnap[temp_plus][POSITION]);
		printf("BARRIER CHECK RNAP BEHIND (%d): BARRRIER BEHIND: "
			"%d; BARRIER AHEAD: %d\n",
			temp_minus,rnap[temp_minus][BARRIER_BEHIND],
			rnap[temp_minus][BARRIER_AHEAD]);
		printf("BARRIER CHECK RNAP AHEAD (%d): BARRRIER BEHIND: "
			"%d; BARRIER AHEAD: %d\n",
			temp_plus,rnap[temp_plus][BARRIER_BEHIND],
			rnap[temp_plus][BARRIER_AHEAD]);
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
			spacing = ((site - up_check >= min_spacing) &&		// Use if interRNAP OK
				   (down_check - site >= min_spacing));
		   }
		   else {
			spacing = 0;						// If activity between RNAPs
		   }								// not allowed, topo cannot
		}								// be placed
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
			spacing = (site-rnap[temp_minus][POSITION]		// No problem with stop
				   >= min_spacing);				// barrier in plasmid
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
			spacing = (rnap[temp_plus][POSITION] -site		// No problem with start
			   	   >= min_spacing);				// barrier in plasmid
		   }
		}
	}
	else {
		if (!(plasmid)) {
		  	spacing =
			  ((stop_barrier - site >= min_spacing) &&
			  (site - start_barrier >= min_spacing));
		}
		else {								// Good to go in plasmid:
			spacing = 1;						// no need to check start
		}								// or stop barriers
		test_lk = open_lk[which_promo];
		test_len = stop_barrier - start_barrier;
	}

	if (!(QUIET_MODE)) {
		printf("FINAL TOPO CHECKS: MINUS ONE %d; PLUS ONE %d; "
			"SPACING %d\n",
			temp_minus,temp_plus,spacing);
	}

	if (spacing) {

		if (!(QUIET_MODE)) {
			printf("NOW CHECKING AGAINST %d TOPOISOMERASES "
				"ON PROMOTER %d\n",
				num_topo[which_promo],which_promo);
		}

		if (plasmid && site < 1) {
			try_site = stop_barrier + site +			// Account for plasmid
				   abs(start_barrier);				// wrap-around
		}
		else {
			try_site = site;
		}

		for (i = 1; i <= num_topo[which_promo]; i++) {			// Check for occlusion by
			check_site =
			   topo[topo_finder[which_promo][i]][POSITION];		// other topoisomerases on
			check_type =
			   topo[topo_finder[which_promo][i]][TOPO_TYPE];	// tx unit (should be cheaper
			if (plasmid && check_site < 1) {
				check_site = stop_barrier + check_site +
					     abs(start_barrier);
			}
			if (abs(try_site - check_site) <			// than checking DNA grid)
			    topo_dist_check[topo_type + check_type]) {
				spacing = 0;
				break;
			}
		}
	}
	if (spacing) {

		if (!(QUIET_MODE)) {
			printf("SPACING CHECKS PASSED!\n");
			printf("FIRST FREE TOPO IS %d\n",first_free_topo);
		}

		if (topo_type == TOPO_IA) {					// Topoisomerase IA binding
			if (test_lk + 1.0 > relaxed_lk[test_len]) {		// prevented if insufficient
				return(0);					// negative supercoiling; currently
			}							// cutoff is ability to perform
		}								// at least one cycle

		if (!(QUIET_MODE)) {
			printf("PREASSIGN\n");
		}

		topo[first_free_topo][TOPO_TYPE] = topo_type;
		topo[first_free_topo][POSITION] = site;
		topo[first_free_topo][SOURCE_PROMO] = which_promo;
		topo[first_free_topo][PLUS_ONE] = temp_plus;
		topo[first_free_topo][MINUS_ONE] = temp_minus;
		topo[first_free_topo][BURST_STATE] = -1;

		if (!(QUIET_MODE)) {
			printf("POSTASSIGN\n");
		}

		num_topo[which_promo]++;
		topo_finder[which_promo][num_topo[which_promo]] =
		   first_free_topo;

		if (!(QUIET_MODE)) {
			printf("POSTASSIGN2\n");
		}

		if (site < 1) {
			topo[first_free_topo][POS_BLOCK] = UPSTREAM;
		}
		else if (site < tx_end) {
			topo[first_free_topo][POS_BLOCK] = MIDDLE;
		}
		else {
			topo[first_free_topo][POS_BLOCK] = DOWNSTREAM;
		}

		if (!(QUIET_MODE)) {
			printf("PRESWITCH\n");
		}	

		Switch_burst_state(time,first_free_topo);

		if (!(QUIET_MODE)) {
			printf("POSTSWITCH\n");
		}

		mu_topo_exit = mu_topo_unbind[topo_type];
		if (topo_type == GYRASE && site < tx_end + 1) {			// Preference for downstream
			mu_topo_exit /= relative_gyrase_ds_affinity;		// binding sites; cf. Sutormin
		}								// et al., NAR 2019
		next_fire = time + Find_firing_time(mu_topo_exit);
		Add_reaction_to_queue(TOPO_UNBINDING_RXN(first_free_topo),
				      next_fire);

		if (!(QUIET_MODE)) {
			printf("UNBINDING SCHEDULED FOR %f\n",next_fire);
		}

		rnap_position = rnap[temp_minus][POSITION];
		adj_topo = rnap[temp_minus][TOPO_DOWNSTREAM];
		if (temp_minus && (!(adj_topo) ||
		    site - rnap_position <
		    topo[adj_topo][POSITION] - rnap_position)) {
			rnap[temp_minus][TOPO_DOWNSTREAM] =
			   first_free_topo;
			topo[first_free_topo][CLOSEST_DOWN] = temp_minus;	// CHECK CLOSEST CATEGORIES!!!
			topo[adj_topo][CLOSEST_DOWN] = 0;
		}
		rnap_position = rnap[temp_plus][POSITION];
		adj_topo = rnap[temp_plus][TOPO_UPSTREAM];
		if (temp_plus && (!(adj_topo) ||
		    rnap_position - site <
		    rnap_position - topo[adj_topo][POSITION])) {
			rnap[temp_plus][TOPO_UPSTREAM] = first_free_topo;
			topo[first_free_topo][CLOSEST_UP] = temp_plus;
			topo[adj_topo][CLOSEST_UP] = 0;
		}

		if (!(QUIET_MODE)) {
			printf("CURRENT DOWNSTREAM TOPO AUTO ASSIGN %d: "
				"TYPE %d AT POSITION %d\n",
				starter_down[which_promo],
				topo[starter_down[which_promo]][TOPO_TYPE],
				topo[starter_down[which_promo]][POSITION]);
		}

		if (site >= 1 && (!(starter_down[which_promo]) ||
		    site - 1 <
		    topo[starter_down[which_promo]][POSITION] - 1)) {
			starter_down[which_promo] = first_free_topo;

			if (!(QUIET_MODE)) {
				printf("DOWNSTREAM TOPO REPLACED BY %d\n",
					first_free_topo);
			}

		}

		if (!(QUIET_MODE)) {
			printf("CURRENT UPSTREAM TOPO AUTO ASSIGN %d: "
				"TYPE %d AT POSITION %d\n",
				starter_up[which_promo],
				topo[starter_up[which_promo]][TOPO_TYPE],
				topo[starter_up[which_promo]][POSITION]);
		}

		if (site <= 0 && (!(starter_up[which_promo]) ||
		    1 - site <
		    1 - topo[starter_up[which_promo]][POSITION])) {
			starter_up[which_promo] = first_free_topo;

			if (!(QUIET_MODE)) {
				printf("UPSTREAM TOPO REPLACED BY %d\n",
					first_free_topo);
			}

		}
		if (abs(site - 1) < topo_rnap_width[topo_type]) {		// Should be able to use 0,1,2 to
			promo_topo_blocked[which_promo]++;			// indicate number of topo's
			topo[first_free_topo][PROMO_BLOCKER]++;			// blocking promoter--subtract
		}								// on unbinding of topoisomerase
		first_free_topo = Find_first_free(topo_search);
	}

	if (!(QUIET_MODE)) {
		printf("ABOUT TO EXIT TOPO BINDING FN\n");
	}

	return(spacing);
}

