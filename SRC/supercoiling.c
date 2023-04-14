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

#define DEFINE_SUPER_VARIABLES
#include "INCL/supercoiling.h"

#include "INCL/sim_types.h"
#include "INCL/gen_def.h"
#include "INCL/sim_settings.h"
#include "INCL/reaction_manager.h"
#include "INCL/replication.h"
#include "INCL/transcription.h"
#include "INCL/translation.h"


void Assign_rotation_resistance()
{
	int i;
	float adj;

	for (i = 1; i <= tx_end; i++) {
		rot_resist[i] = pow((float) i,	resist_exp);
	}
	adj = max_resist - min_resist;
	for (i = 1; i <= tx_end; i++) {
		rot_resist[i] = min_resist +
				adj*(rot_resist[i]/rot_resist[tx_end]);
	}

	return;
}

float Calculate_net_torque(int rnap_index,int relax_call)
{
	float upstream_contrib,downstream_contrib;
	if (supercoiling_off) {
		return(0.0);
	}
	upstream_contrib   = -1.0 * torque_factor * 
			     twist[rnap_index][UPSTREAM];
	downstream_contrib = torque_factor *
			     twist[rnap_index][DOWNSTREAM];
	if (!(relax_call) &&
	    twist[rnap_index][UPSTREAM] <= min_upstr_sigma) {
		return(11.1);
	}
	return(rot_resist[rnap[rnap_index][POSITION]] *
	       (upstream_contrib + downstream_contrib));
}
	
double Torque_dependent_mu_forward(int rnap_index)
{
        double torque,adj_mu;
        torque = (double) Calculate_net_torque(rnap_index,0);

        adj_mu =  -0.0002 * (pow(torque,5.0)) +
                   0.0008 * (pow(torque,4.0)) +
                   0.0041 * (pow(torque,3.0)) -
                    0.035 * (pow(torque,2.0)) -
                   0.2166 * torque +
                  22.3574;
        adj_mu /= 22.3574;
        if (adj_mu > 2.0) {
                adj_mu = 2.0;
        }
        if (adj_mu < 0.0001) {
                return(0.0);
        }

        return(adj_mu);

}

	
double Torque_dependent_mu_enter_P1(int rnap_index)
{
	int cause;
	double torque,rel_pause, p_pause, new_kp;
	double old_kp, old_kf, adj = 1.0;
	cause = Determine_dominant_trailer(rnap_index,1/rnap_p1_effect,
					   1/ribo_p1_effect);
	torque = (double) Calculate_net_torque(rnap_index,0);
	if (torque < -10.5) {
		return (INFINITY);
	}
	if (ANTICASCADE_ON &&
	    Check_for_RNAP_RNAP_contacts(rnap_index,PLUS_ONE)) {		// Assumes that P1
		return (INFINITY);						// requires partial
	}									// forward translocation
										// and is unavailable to
										// a blocked RNAP
	if (torque > 11) {
		return(-9999.99);
	}
	cause++;
	if (torque <= 7.5) {
		rel_pause = (0.0001*torque*torque) + (0.0022*torque) +
			   0.0128;
	}
	else {
		rel_pause = (0.0453*torque*torque) - (0.5621*torque) +
			   1.7032;
	}	
	rel_pause /= torque_zero_pause_p;
	p_pause = 1.0 - rnap_dwell[rnap[rnap_index][POSITION]].p_on_f;
	rel_pause *= adj;
	if (rel_pause * p_pause >= 1.0) {
		return(-9999.99);
	}
	else {
		old_kp = 1.0/rnap_dwell[rnap[rnap_index][POSITION]].mu_p;	// Recalculate kp:
		old_kf = 1.0/rnap_dwell[rnap[rnap_index][POSITION]].mu_f;	// ((relative*kp*kf)/
		new_kp = (rel_pause * old_kf * old_kp *				// (kp+kf)) / 
		   rnap_dwell[rnap[rnap_index][POSITION]].mu_net_pre) /		// (1 - ((relative*kp)/
		   (1.0 - (rel_pause * old_kp *					// (kp+kf)))
		   rnap_dwell[rnap[rnap_index][POSITION]].mu_net_pre));		// (Note: net_pre = 
		return(1.0/new_kp);						//  1/(kf+kp))
	}
}

double Torque_dependent_mu_exit(int rnap_index)
{
	int trailer;
	double torque,relative_dur,adj=1.0;
	trailer = Determine_dominant_trailer(rnap_index,rnap_e1_effect,
					     ribo_e1_effect);
	if (rnap[rnap_index][STATE] == OFF_P1) {				// Assuming for now that
		torque = (double) Calculate_net_torque(rnap_index,0);		// exit rates from P2
		if (torque <= -10.0) {						// and P3 are independent
			return(-9999.99);					// of torque; changes in
		}								// pause duration under
		if (torque >= 10.0) {						// torque come only from
			return(INFINITY);					// entry into elemental
		}								// pause...
		relative_dur = 0.1914*(tan(((PI/23.0)*torque)+0.1469))
			       + 0.5298;
		relative_dur /= torque_zero_pause_dur;

		if (trailer == INDEPENDENT) {
			adj = 1.0;
		}
		relative_dur *= adj;
	}
	else {									// Will not get here unless
		if ((trailer_list[0] != 0 && RNAP_ANTIPAUSE) ||			// RNAP has already passed
		    (trailer_list[1] != 0 && RIBO_ANTIPAUSE)) {			// back-in-register test;
			return(-9999.99);					// effect is currently 
		}								// all-or-nothing...
		torque = (double) Calculate_net_torque(rnap_index,0);
		if (torque < P23_release_torque) {
			return(-9999.99);
		}
		relative_dur = 1.0;
	}
	return(mu_exit[rnap[rnap_index][STATE]]*relative_dur);
}		

float Apply_180bp_partition_plot(int rnap_index, int side)
{
	float twist_contrib;
	if (sigma[rnap_index][side] >= -0.05 &&
	    sigma[rnap_index][side] <=  0.05) {
		twist_contrib = sigma[rnap_index][side];
	}
	else if (sigma[rnap_index][side] < 0.0 &&
		 sigma[rnap_index][side] > -0.09) {
		twist_contrib = (1.0 -
		  (0.5 * ((sigma[rnap_index][side] -
		  -0.05)/-0.04)));
		twist_contrib *= sigma[rnap_index][side];
	}
	else if (sigma[rnap_index][side] > 0.0 &&
		 sigma[rnap_index][side] < 0.07) {
		twist_contrib = (1.0 -
		  (0.5 * ((sigma[rnap_index][side] -
		  0.05)/0.02)));
		twist_contrib *= sigma[rnap_index][side];
	}
	else {
		twist_contrib = 0.5 * sigma[rnap_index][side];
	}
	return(twist_contrib);
}

float Apply_435bp_partition_plot(int rnap_index, int side)
{
	float twist_contrib;
	if (sigma[rnap_index][side] >= -0.03 &&
	    sigma[rnap_index][side] <=  0.02) {
		twist_contrib = sigma[rnap_index][side];
	}
	else if (sigma[rnap_index][side] < 0.0 &&
		 sigma[rnap_index][side] > -0.09) {
		twist_contrib = (1.0 -
		  (0.67 * ((sigma[rnap_index][side] -
		  -0.03)/-0.06)));
		twist_contrib *= sigma[rnap_index][side];
	}
	else if (sigma[rnap_index][side] > 0.0 &&
		 sigma[rnap_index][side] < 0.09) {
		twist_contrib = (1.0 -
		  (0.67 * ((sigma[rnap_index][side] -
		  0.02)/0.07)));
		twist_contrib *= sigma[rnap_index][side];
	}
	else {
		twist_contrib = sigma[rnap_index][side] * 0.33;
	}
	return(twist_contrib);
}

float Apply_885bp_partition_plot(int rnap_index, int side)
{
	float twist_contrib;
	if (sigma[rnap_index][side] >= -0.02 &&
	    sigma[rnap_index][side] <=  0.01) {
		twist_contrib = sigma[rnap_index][side] * 0.5;
	}
	else if (sigma[rnap_index][side] > -0.11 &&
		 sigma[rnap_index][side] < -0.02) {
		twist_contrib = 0.5 -
		  (0.2 * ((sigma[rnap_index][side] - -0.02) / -0.09));
		twist_contrib *= sigma[rnap_index][side];
	}
	else if (sigma[rnap_index][side] < 0.11 &&
		 sigma[rnap_index][side] > 0.01) {
		twist_contrib = 0.5 -
		  (0.2 * ((sigma[rnap_index][side] - 0.01) / 0.1));
		twist_contrib *= sigma[rnap_index][side];
	}
	else {
		twist_contrib = sigma[rnap_index][side] * 0.3;
	}
	return(twist_contrib);
}

void Distribute_twist_and_writhe(int rnap_index, int side)
{
	float frac_shorter,frac_longer;
	if (!(writhe_partition)) {
		twist[rnap_index][side] = sigma[rnap_index][side];
		helical_repeat[rnap_index][side] = 10.5 *
		   (1.0 - twist[rnap_index][side]);
		return;
	}
	
	if (rnap_sep[rnap_index][side] <= 150.0) {				// Twist here is additional
		twist[rnap_index][side] = sigma[rnap_index][side];		// twist relative to relaxed
	}									// B-DNA
	else if (rnap_sep[rnap_index][side] <= 180.0) {
		frac_longer = ((float) (rnap_sep[rnap_index][side]-150))/
			       30.0;
		frac_shorter = 1 - frac_longer;
		twist[rnap_index][side] =
		  (frac_shorter * sigma[rnap_index][side]) +
		  (frac_longer *
		   Apply_180bp_partition_plot(rnap_index,side));
	}
	else if (rnap_sep[rnap_index][side] <= 435.0) {
		frac_longer = ((float) (rnap_sep[rnap_index][side]-180))/
			       255.0;
		frac_shorter = 1 - frac_longer;
		twist[rnap_index][side] =
		  (frac_shorter *
		   Apply_180bp_partition_plot(rnap_index,side)) +
		  (frac_longer *
		   Apply_435bp_partition_plot(rnap_index,side));
	}
	else if (rnap_sep[rnap_index][side] <= 885.0) {
		frac_longer = ((float) (rnap_sep[rnap_index][side]-435))/
			       450.0;
		frac_shorter = 1 - frac_longer;
		twist[rnap_index][side] =
		  (frac_shorter *
		   Apply_435bp_partition_plot(rnap_index,side)) +
		  (frac_longer *
		   Apply_885bp_partition_plot(rnap_index,side));
	}
	else {
		twist[rnap_index][side] =
		  Apply_885bp_partition_plot(rnap_index,side);
	}

	helical_repeat[rnap_index][side] = 10.5 *
	   (1.0 - twist[rnap_index][side]);
	return;
}	

void Process_supercoiling_changes(int rnap_index, int move,
				  int side_adj)
{
	int newest_rnap,oldest_rnap,check_rnap;
	int single_rnap=0,connected=0;
	int skip_upstream_lk=0,skip_downstream_lk=0;
	int promo_index;
	int up_adj=0,down_adj=0;
	float twist_add,tot_link=0.0;
	float torque;
	float relax_add;

                promo_index = rnap[rnap_index][SOURCE_PROMO];
                oldest_rnap = promoter[promo_index][OLDEST_RNAP];
                newest_rnap = promoter[promo_index][NEWEST_RNAP];
                if (rnap[newest_rnap][STATE] == CLOSED) {
                        newest_rnap = rnap[newest_rnap][PLUS_ONE];
                }
		if (plasmid &&
		    rnap[newest_rnap][BARRIER_BEHIND] <= 0 &&
		    rnap[oldest_rnap][BARRIER_AHEAD] <= 0) {
			connected++;
		}
		if (rnap_index == newest_rnap &&
		    rnap_index == oldest_rnap) {
			single_rnap++;
		}
		twist_add = (1.0/10.5);
		if (move < 0.5) {
			twist_add *= (1.0 - 
			   rot_resist[rnap[rnap_index][POSITION]]);
			twist_add *=
			   (1.0 + twist[rnap_index][UPSTREAM]);
		}
		else {
			twist_add *= (1.0 - 
			   rot_resist[rnap[rnap_index][POSITION]]);
			twist_add *=
			   (1.0 + twist[rnap_index][DOWNSTREAM]);
		}

		if (single_rnap && connected) {					// Movement of a single,
			skip_upstream_lk++;					// connected RNAP on a 
			skip_downstream_lk++;					// plasmid template leaves
		}								// spacing and Lk unchanged
		else {
			if (!(DNA_SCRUNCH) ||
			    rnap[rnap_index][POSITION] > hybrid_length) {
				rnap_sep[rnap_index][UPSTREAM] += move;
			}
			if (DNA_SCRUNCH &&					// Extended bubble reduced
			    rnap[rnap_index][POSITION] ==			// on dissociation of open
			    hybrid_length + 1) {				// complex from promoter;
				rnap_sep[rnap_index][UPSTREAM] +=		// transcription from this
				  ((hybrid_length - 1)*(move==1));		// point forward is with
			}							// standard (constant) bubble
			rnap_sep[rnap_index][DOWNSTREAM] -= move;
				
			
			if (plasmid_linearized &&				// Spacing and Lk changed
			    rnap_index == newest_rnap &&			// for newest RNAP, but in
			    rnap[rnap_index][BARRIER_BEHIND] <= 0) {		// the absence of a barrier
				skip_upstream_lk++;				// sigma is constant
				lk[rnap_index][UPSTREAM] = relaxed_lk[		// (relaxed) upstream
				   rnap_sep[rnap_index][UPSTREAM]];
				sigma[rnap_index][UPSTREAM] = 0.0;
			}
			if (plasmid_linearized &&				// So too for the oldest 
			    rnap_index == oldest_rnap &&			// RNAP downstream
			    rnap[rnap_index][BARRIER_AHEAD] <= 0) {
				skip_downstream_lk++;
				lk[rnap_index][DOWNSTREAM] = relaxed_lk[
				   rnap_sep[rnap_index][DOWNSTREAM]];
				sigma[rnap_index][DOWNSTREAM] = 0.0;
			}
		}
	
		relax_add = 1.0/105.0;						// Section preventing
										// time-wasting oscillations
		if ((side_adj == ASSIST_TORQUE_RELAX ||				// across zero-torque states
		     side_adj == RESIST_TORQUE_RELAX) &&			// at low viscosity/drag
		    !(skip_upstream_lk && skip_downstream_lk)) {
		    if (side_adj == ASSIST_TORQUE_RELAX) {
		     curr_rnap_rot[rnap_index]--;
		     if (skip_downstream_lk) {					// For last RNAP on open
		      if (/*lk[rnap_index][UPSTREAM] > 				// template: if relaxation
			  relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] &&	*/	// changes sign of sigma from
			  lk[rnap_index][UPSTREAM] - relax_add <		// + to - upstream, reduce dLk
			  relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]) {		// addition to difference
			   relax_add = fabs(lk[rnap_index][UPSTREAM] -		// between LkO and Lk (sigma
			    relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]);	// is always zero downstream)
		       }
		      }	
		      else if (skip_upstream_lk) {				// For first RNAP on open
		       if (/*lk[rnap_index][DOWNSTREAM] <			// template: if relaxation
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]] && */	// changes sign downstream
			   lk[rnap_index][DOWNSTREAM] + relax_add >		// from - to +, reduce dLk
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]) {	// to difference from LkO
			    relax_add = fabs(					// (sigma zero upstream)
			     relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]] -
			     lk[rnap_index][DOWNSTREAM]);
		       }
		      }
		      else {							// For all RNAPs within closed
		       if (/*lk[rnap_index][UPSTREAM]/				// topological domains and
			   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] >		// RNAPs bracketed by two RNAPs
			   lk[rnap_index][DOWNSTREAM]/				// on open templates: if Lk add 
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]		// changes sign of torque (i.e.,
			   &&	*/						// assistive -> resistive), 
			   (lk[rnap_index][UPSTREAM] - relax_add)/		// reduce dLk by assuming Lk
			   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] <		// equilibrates withing domain,
			   (lk[rnap_index][DOWNSTREAM] + relax_add)/		// resulting in zero net torque.
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]) {	// Note that |dLk| can be used
			    relax_add = (lk[rnap_index][UPSTREAM] +		// here b/c it's always assumed
					 lk[rnap_index][DOWNSTREAM])/		// to be positive below
					relaxed_lk[
					 rnap_sep[rnap_index][UPSTREAM] +	// B/c relaxation conserves Lk
					 rnap_sep[rnap_index][DOWNSTREAM]];	// w/in domain centered on RNAP,
			    relax_add = fabs(					// the same dLk will balance
			     relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]*	// sigma on either side of RNAP
			     relax_add - lk[rnap_index][DOWNSTREAM]);
		       }
		      }
		     }
		    else {							// Under resistive torque:
		     curr_rnap_rot[rnap_index]++;
		     if (skip_downstream_lk) {					// For last RNAP on open
		      if (/*lk[rnap_index][UPSTREAM] < 				// template: if relaxation
			  relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] &&	*/	// changes sign of sigma from
			  lk[rnap_index][UPSTREAM] + relax_add >		// - to + upstream, proceed
			  relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]) {		// as above
			   relax_add = fabs(lk[rnap_index][UPSTREAM] -
			    relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]);
		       }
		      }	
		      else if (skip_upstream_lk) {				// For first RNAP on open
		       if (/*lk[rnap_index][DOWNSTREAM] >			// template: if relaxation
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]] && */	// changes sign downstream
			   lk[rnap_index][DOWNSTREAM] - relax_add <		// from + to -, etc.
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]) {
			    relax_add = fabs(
			     relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]] -
			     lk[rnap_index][DOWNSTREAM]);
		       }
		      }
		      else {							// For all RNAPs within closed
		       if (/*lk[rnap_index][UPSTREAM]/				// topological domains and
			   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] >		// RNAPs bracketed by two RNAPs
			   lk[rnap_index][DOWNSTREAM]/				// on open templates: if Lk add 
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]		// changes sign of torque (i.e.,
			   &&	*/						// resistive -> assistive), etc. 
			   (lk[rnap_index][UPSTREAM] + relax_add)/
			   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]] >
			   (lk[rnap_index][DOWNSTREAM] - relax_add)/
			   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]) {
			    relax_add = (lk[rnap_index][UPSTREAM] +
					 lk[rnap_index][DOWNSTREAM])/
					relaxed_lk[
					 rnap_sep[rnap_index][UPSTREAM] +
					 rnap_sep[rnap_index][DOWNSTREAM]];
			    relax_add = fabs(
			     relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]*
			     relax_add - lk[rnap_index][DOWNSTREAM]);
			}
		       }
		      }
		}
			   

		if (!(skip_upstream_lk)) {
			if (side_adj == ASSIST_TORQUE_RELAX) {
				lk[rnap_index][UPSTREAM] -=
				   relax_add;
			}
			else if (side_adj == RESIST_TORQUE_RELAX) {
				lk[rnap_index][UPSTREAM] +=
				   relax_add;
			}
			else {						
				lk[rnap_index][UPSTREAM] +=
				   (twist_add * ((float) move));
			}
		}
		if (!(skip_downstream_lk)) {
			if (side_adj == ASSIST_TORQUE_RELAX) {
				lk[rnap_index][DOWNSTREAM] +=
				   relax_add;
			}
			else if (side_adj == RESIST_TORQUE_RELAX) {
				lk[rnap_index][DOWNSTREAM] -=
				   relax_add;
			}
			else {
				lk[rnap_index][DOWNSTREAM] -=
				   (twist_add * ((float) move));
			}
		}


		sigma[rnap_index][UPSTREAM] = -1.0 +
		   (lk[rnap_index][UPSTREAM]/
		   relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]);
		sigma[rnap_index][DOWNSTREAM] = -1.0 +
		   (lk[rnap_index][DOWNSTREAM]/
		   relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]);

	Distribute_twist_and_writhe(rnap_index,UPSTREAM);
	Distribute_twist_and_writhe(rnap_index,DOWNSTREAM);

	if (move || side_adj == BOTH || side_adj == UPSTREAM ||
	    side_adj == ASSIST_TORQUE_RELAX ||
	    side_adj == RESIST_TORQUE_RELAX) {
		up_adj++;
	}
	if (move || side_adj == BOTH || side_adj == DOWNSTREAM ||
	    side_adj == ASSIST_TORQUE_RELAX ||
	    side_adj == RESIST_TORQUE_RELAX) {
		down_adj++;
	}

      if (up_adj) {
	if (rnap[rnap_index][MINUS_ONE] &&
	    rnap[rnap_index][BARRIER_BEHIND] <= 0) {
		rnap_sep[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] =
		   rnap_sep[rnap_index][UPSTREAM];
		lk[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] =
		   lk[rnap_index][UPSTREAM];
		sigma[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] =
		   sigma[rnap_index][UPSTREAM];
		helical_repeat[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] =
		   helical_repeat[rnap_index][UPSTREAM];
		twist[rnap[rnap_index][MINUS_ONE]][DOWNSTREAM] = 
		   twist[rnap_index][UPSTREAM];
	}
	if (rnap_index == newest_rnap && connected) {
                rnap_sep[oldest_rnap][DOWNSTREAM] =
                   rnap_sep[rnap_index][UPSTREAM];
                lk[oldest_rnap][DOWNSTREAM] =
                   lk[rnap_index][UPSTREAM];
                helical_repeat[oldest_rnap][DOWNSTREAM] =
                   helical_repeat[rnap_index][UPSTREAM];
                twist[oldest_rnap][DOWNSTREAM] =
                   twist[rnap_index][UPSTREAM];
        }
      }
      if (down_adj) {
	if (rnap[rnap_index][PLUS_ONE] &&
	    rnap[rnap_index][BARRIER_AHEAD] <= 0) {
		rnap_sep[rnap[rnap_index][PLUS_ONE]][UPSTREAM] =
		   rnap_sep[rnap_index][DOWNSTREAM];
		lk[rnap[rnap_index][PLUS_ONE]][UPSTREAM] =
		   lk[rnap_index][DOWNSTREAM];
		sigma[rnap[rnap_index][PLUS_ONE]][UPSTREAM] =
		   sigma[rnap_index][DOWNSTREAM];
		helical_repeat[rnap[rnap_index][PLUS_ONE]][UPSTREAM] =
		   helical_repeat[rnap_index][DOWNSTREAM];
		twist[rnap[rnap_index][PLUS_ONE]][UPSTREAM] = 
		   twist[rnap_index][DOWNSTREAM];
	}
	if (rnap_index == oldest_rnap && connected) {
                rnap_sep[newest_rnap][UPSTREAM] =
                   rnap_sep[rnap_index][DOWNSTREAM];
                lk[newest_rnap][UPSTREAM] =
                   lk[rnap_index][DOWNSTREAM];
                helical_repeat[newest_rnap][UPSTREAM] =
                   helical_repeat[rnap_index][DOWNSTREAM];
                twist[newest_rnap][UPSTREAM] =
                   twist[rnap_index][DOWNSTREAM];
        }
      }

	if (!(QUIET_MODE)) {
		if (rnap[rnap_index][SOURCE_PROMO] == 1 &&
		    (side_adj != RESIST_TORQUE_RELAX && side_adj !=
		     ASSIST_TORQUE_RELAX)) {
			printf("RNAP %d AT %d; HELICAL RPT: %f; ROTATION "
				"FRACTION: %f; TWIST ADD: %f; DIST "
				"UPSTREAM: %d; SIGMA UPSTREAM %f; TWIST "
				"UPSTREAM: %f; DIST DOWNSTREAM: %d; "
				"SIGMA DOWNSTREAM: %f; TWIST "
				"DOWNSTREAM: %f\n",
				rnap_index, rnap[rnap_index][POSITION],
				helical_repeat[rnap_index][DOWNSTREAM],
				rot_resist[rnap[rnap_index][POSITION]],
				twist_add,rnap_sep[rnap_index][UPSTREAM],
				sigma[rnap_index][UPSTREAM],
				twist[rnap_index][UPSTREAM],
				rnap_sep[rnap_index][DOWNSTREAM],
				sigma[rnap_index][DOWNSTREAM],
				twist[rnap_index][DOWNSTREAM]);
			printf("LK UPSTR: %f\n",
			    lk[rnap_index][UPSTREAM]);
			printf("LK DOWNSTR: %f\n",
			    lk[rnap_index][DOWNSTREAM]);
			printf("CF UPSTR: %f\n",
			    relaxed_lk[rnap_sep[rnap_index][UPSTREAM]]);
			printf("CF DOWNSTR: %f\n",
			    relaxed_lk[rnap_sep[rnap_index][DOWNSTREAM]]);
			check_rnap = oldest_rnap;
			while (check_rnap) {
				tot_link += lk[check_rnap][UPSTREAM];
				check_rnap = rnap[check_rnap][MINUS_ONE];
			}
			printf("TOTAL LINKING NUMBER: %f\n",tot_link);
			torque = Calculate_net_torque(rnap_index,0);
			printf("NET TORQUE RNAP %d: %f; RNAP BEHIND (%d): "
				"%f; RNAP AHEAD (%d): %f\n",
				rnap_index,torque,
				rnap[rnap_index][MINUS_ONE],
				Calculate_net_torque(
				  rnap[rnap_index][MINUS_ONE],0),
				rnap[rnap_index][PLUS_ONE],
				Calculate_net_torque(
				  rnap[rnap_index][PLUS_ONE],0));
			if (NO_PAUSE_MODEL) {
				printf("MU F FOR RNAP %d IS: %f; "
					"ADJUSTED FOR TORQUE: %f\n",
					rnap_index,
					rnap_dwell[
					  rnap[rnap_index][POSITION]].mu_f,
					Torque_dependent_mu_forward(
					  rnap_index));
			}
			else {
				printf("MU ENTER P1 UNADJ FOR RNAP %d IS: "
					"%f; ADJUSTED FOR TORQUE: %f\n",
					rnap_index,
					rnap_dwell[
					  rnap[rnap_index][POSITION]].mu_p,
					Torque_dependent_mu_enter_P1(
					  rnap_index));
				printf("UNADJ PROB P1 ENTRY: %f; mu_f %f; "
					"mu_p %f\n",
					1.0 - 
					  rnap_dwell[
					   rnap[rnap_index][POSITION]].p_on_f,
					rnap_dwell[
					  rnap[rnap_index][POSITION]].mu_f,
					rnap_dwell[
					  rnap[rnap_index][POSITION]].mu_p);
				printf("RELATIVE OFF PATHWAY EXIT "
					"ADJ: %f\n",
					Torque_dependent_mu_exit(
					  rnap_index));
				printf("\n");
			}
		}
	}
}

void Update_RNAP_complex_resistance(int rnap_index)
{
	int rna_index;
	double rna_contrib,ribo_contrib;
	if (ribo_resist_model == FLAT_RESISTANCE) {
		complex_resist[rnap_index] = flat_resist_param;
		return;
	}
	rna_index = rnap[rnap_index][RNAP_RNA];
	ribo_contrib = (double) rna[rna_index][TOTAL_RIBO];
	rna_contrib = (double) (rna[rna_index][THREE_END] -
			       rna[rna_index][FIVE_END]);
	if (ribo_resist_model == IMPLICIT_RNAP_POS) {
		ribo_contrib = ((double) rnap[rnap_index][POSITION])/
			        mean_interribo_dist;
	}
	else if (ribo_resist_model == IMPLICIT_EXTANT_RNA) {
		ribo_contrib = rna_contrib/mean_interribo_dist;
	}
	else if (ribo_resist_model == EXPLICIT_RIBOSOMES) {
		ribo_contrib = (double) rna[rna_index][TOTAL_RIBO];
	}
	else {									// RNA-only sphere; compare
		ribo_contrib = 0.0;						// only with RNAP "sphere"
	}	
	rna_contrib = pow(rna_contrib*0.34, resist_exp);
	rna_contrib *= 0.0005;
	complex_resist[rnap_index] = 0.12 +ribo_contrib*0.96;
	if (rna_contrib > complex_resist[rnap_index]) {
		complex_resist[rnap_index] = rna_contrib;
	} 
	complex_resist[rnap_index] =
		1.0/(complex_resist[rnap_index]+0.05);
	return;
}

int Update_relaxation_rates(double time, int rnap_index)
{
	int reaction;
	double torque,next_fire;
	double new_mu;
	if (!(rnap_relax)) {
		return(0);
	}
	if (!(rnap_index) || rnap[rnap_index][STATE] == CLOSED) {
		return(0);
	}
	reaction = RNAP_RELAX_RXN(rnap_index);
	torque = (double) Calculate_net_torque(rnap_index,1);
	if (fabs(torque) < 0.01) {
		Update_reaction_queue(reaction,INFINITY);
		last_rates[rnap_index].last_resist_mu = -9999.9;
		last_rates[rnap_index].last_assist_mu = -9999.9;
		return(0);
	}
	new_mu = fabs(torque*complex_resist[rnap_index]);
	new_mu /= (2*PI);
	new_mu = ((1.0/105.0)/new_mu);
	new_mu *= relax_factor;
	if (torque > 0 && last_rates[rnap_index].last_resist_mu > 0) {
	       next_fire =
		  reaction_queue[reaction_index[reaction]].reaction_time;
	       next_fire -= time;
	       next_fire *= (new_mu/last_rates[rnap_index].last_resist_mu);
	       next_fire += time;
	}
	else if (torque < 0 && last_rates[rnap_index].last_assist_mu > 0) {
	       next_fire =
		  reaction_queue[reaction_index[reaction]].reaction_time;
	       next_fire -= time;
	       next_fire *= (new_mu/last_rates[rnap_index].last_assist_mu);
	       next_fire += time;
	}
	else {
	       next_fire = time + Find_firing_time(new_mu);
	}
	if (torque > 0) {
		last_rates[rnap_index].last_resist_mu = new_mu;
		last_rates[rnap_index].last_assist_mu = -9999.9;
	}
	else {
		last_rates[rnap_index].last_resist_mu = -9999.9;
		last_rates[rnap_index].last_assist_mu = new_mu;
	}
	Update_reaction_queue(reaction,next_fire);
	return(1);
}	

int Relax_RNAP_complex(double time, int rnap_index)
{
	double torque;
	last_rates[rnap_index].last_resist_mu = -9999.9;
	last_rates[rnap_index].last_assist_mu = -9999.9;
	torque = (double) Calculate_net_torque(rnap_index,1);

	if (torque > 0) {
		Process_supercoiling_changes(
		   rnap_index,0,RESIST_TORQUE_RELAX);
	}
	else {
		Process_supercoiling_changes(
		   rnap_index,0,ASSIST_TORQUE_RELAX);
	}
	Update_RNAP_rates(time,rnap_index);
	Update_RNAP_rates(time,rnap[rnap_index][MINUS_ONE]);
	Update_RNAP_rates(time,rnap[rnap_index][PLUS_ONE]);
	return(1);
}

