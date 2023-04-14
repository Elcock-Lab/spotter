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

#define DEFINE_REG_VARIABLES
#include "INCL/regulation.h"

#include "INCL/sim_types.h"
#include "INCL/gen_def.h"
#include "INCL/sim_settings.h"
#include "INCL/reaction_manager.h"
#include "INCL/replication.h"
#include "INCL/transcription.h"
#include "INCL/translation.h"
#include "INCL/supercoiling.h"
#include "INCL/topo.h"
#include "INCL/utilities.h"

	
int Activate_promoter(double time, int *promo_cycle, int promo_index,
                      float promoter_on_log[][MAX_PROMO_CYCLES][2])
{
	double next_fire;
	promoter[promo_index][ON_OFF] = ON;
	promo_cycle[promo_index]++;
	promoter_on_log[promo_index][promo_cycle[promo_index]][0] = time;
	Update_reaction_queue(promo_index,INFINITY);				// Here we assume that
	next_fire = time + Find_firing_time(mu_promo_off);			// promoters > 1 will have
	Update_reaction_queue(PROMO_OFF_RXN(promo_index), next_fire);		// these rxns added as they
	next_fire = time + Find_firing_time(mu_promo_load);			// come into being so they
	Update_reaction_queue(RNAP_LOAD_RXN(promo_index), next_fire);		// need only be updated
	return(1);
}

int Deactivate_promoter(double time, int *promo_cycle, int promo_index,
                        float promoter_on_log[][MAX_PROMO_CYCLES][2],
			int shutoff)
{
	double next_fire;
	promoter[promo_index][ON_OFF] = OFF;
	promoter_on_log[promo_index][promo_cycle[promo_index]][1] = time;
	Update_reaction_queue(PROMO_OFF_RXN(promo_index),INFINITY);
	Update_reaction_queue(RNAP_LOAD_RXN(promo_index),INFINITY);
	if (!(shutoff)) {
		next_fire = time + Find_firing_time(mu_promo_on);
	}
	else {
		printf("SHUTTING DOWN PROMOTER %d AT TIME %f\n",
			promo_index,time);
		next_fire = INFINITY;
	}
	Update_reaction_queue(promo_index,next_fire);
	return(1);
}

void Shut_down_all_promoters(double time, int *promo_cycle, 
                             float promoter_on_log[][MAX_PROMO_CYCLES][2],
			     int promo_index)
{
	int i,result;
	for (i = 1; i <= num_promoters; i++) {
		result = Deactivate_promoter(time,promo_cycle,i,
			    promoter_on_log,1);
	}
	for (i = 1; i < MAX_RNAP; i++) {
		if (rnap[i][POSITION]) {
			shutoff_position[i] = rnap[i][POSITION];
		}
	}
	Update_reaction_queue(PROMO_KILL_RXN,INFINITY);
	Remove_reaction(PROMO_KILL_RXN);
	result++;
	return;
}
		
