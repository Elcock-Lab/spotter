#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define DEFINE_REPL_VARIABLES
#include "INCL/replication.C_RNG.h"

#include "INCL/sim_types.C_RNG.h"
#include "INCL/gen_def.C_RNG.h"
#include "INCL/sim_settings.C_RNG.h"
#include "INCL/reaction_manager.C_RNG.h"
#include "INCL/regulation.C_RNG.h"
#include "INCL/transcription.C_RNG.h"
#include "INCL/translation.C_RNG.h"
#include "INCL/topo.C_RNG.h"



void Start_from_zero(float fract)
{
	int i,add=0;
	float x;
	for (i = 1; i <= birth_promo; i++) {
		promoter[i][ON_OFF] = OFF;					// First promoter is present at start:
		promoter[i][LOC_ON_OFF] = 1;					// Add as first element of action list and
		promoter[i][LAST_BP] = tx_end;					// draw to determine if on or off at
										// start of trajectory
		x = Sample_uniform_distribution();
		printf("DRAW WAS %f; TARGET %f\n",x,fract);
		reaction_queue[add+1].reaction = i;
		reaction_queue[add+2].reaction = PROMO_OFF_RXN(i);
		reaction_queue[add+3].reaction = RNAP_LOAD_RXN(i);
		reaction_index[i] = add + 1;
		reaction_index[PROMO_OFF_RXN(i)] = add + 2;
		reaction_index[RNAP_LOAD_RXN(i)] = add + 3;
		if (x <= fract) {
			reaction_queue[add + 2].reaction_time =
				Find_firing_time(mu_promo_off);
			reaction_queue[add + 3].reaction_time =
				Find_firing_time(mu_promo_load);

			promoter[i][ON_OFF] = ON;
		}
		else {
			reaction_queue[add + 1].reaction_time =
				Find_firing_time(mu_promo_on);
		}
		add += 3;
	}

	Build_initial_reaction_queue(3*birth_promo);


	first_free_rnap = 1;
	first_free_ribo = 1;
	first_free_rna = 1;
	first_free_topo = 1;
	return;
}
