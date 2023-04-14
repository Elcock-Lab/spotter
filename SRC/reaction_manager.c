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

#define DEFINE_RXN_VARIABLES
#include "INCL/reaction_manager.h"

#include "INCL/sim_types.h"
#include "INCL/gen_def.h"
#include "INCL/sim_settings.h"
#include "INCL/replication.h"
#include "INCL/regulation.h"
#include "INCL/translation.h"
#include "INCL/utilities.h"


void Swap_queue_positions(int i, int j)
{
        int temp_reaction;
        double temp_time;
        reaction_index[reaction_queue[i].reaction] = j;
        reaction_index[reaction_queue[j].reaction] = i;
        temp_reaction = reaction_queue[i].reaction;
        temp_time = reaction_queue[i].reaction_time;
        reaction_queue[i].reaction = reaction_queue[j].reaction;
        reaction_queue[i].reaction_time = reaction_queue[j].reaction_time;
        reaction_queue[j].reaction = temp_reaction;
        reaction_queue[j].reaction_time = temp_time;
        return;
}

void Update_reaction_position(int position)
{
        int min_child;
        if (reaction_queue[position].reaction_time <
            reaction_queue[PARENT(position)].reaction_time) {
                Swap_queue_positions(position,PARENT(position));
                Update_reaction_position(PARENT(position));
        }
        else if (MIN_TIME(
	  reaction_queue[FIRST_CHILD(position)].reaction_time,
          reaction_queue[SECOND_CHILD(position)].reaction_time) <
          reaction_queue[position].reaction_time) {
            min_child =
		FIRST_CHILD(position) +
                MIN_CHILD(
		  reaction_queue[FIRST_CHILD(position)].reaction_time,
                  reaction_queue[SECOND_CHILD(position)].reaction_time);
            Swap_queue_positions(position,min_child);
            Update_reaction_position(min_child);
        }
        else {
            return;
        }
}

void Update_reaction_queue(int reaction, double time)
{
        reaction_queue[reaction_index[reaction]].reaction_time = time;
        Update_reaction_position(reaction_index[reaction]);
        return;
}

void Add_reaction_to_queue(int reaction, double time)
{
        queue_size++;
        reaction_queue[queue_size].reaction = reaction;
        reaction_queue[queue_size].reaction_time = time; 
        reaction_index[reaction] = queue_size;
        Update_reaction_position(queue_size);
        return;
}

void Order_queue(int i)
{
        int min;
        int first_child,second_child;
        min = i;
        first_child = FIRST_CHILD(i);
        second_child = SECOND_CHILD(i);
        if (first_child <= queue_size) {
                if (reaction_queue[first_child].reaction_time <
                    reaction_queue[i].reaction_time) {
                        min = first_child;
                }
        }
        if (second_child <= queue_size) {
                if (reaction_queue[second_child].reaction_time <
                    reaction_queue[min].reaction_time) {
                        min = second_child;
                }
        }
        if (min != i) {
                Swap_queue_positions(i,min);
                Order_queue(min);
        }
        return;
}

void Build_initial_reaction_queue(int length)
{
        int i;
        queue_size = length;
        for (i = queue_size/2; i >= 1; i--) {
                Order_queue(i);
        }
        return;
}

void Remove_reaction(int reaction)
{
	int i;
	i = reaction_index[reaction];
	Swap_queue_positions(i,queue_size);
	queue_size--;
	Update_reaction_position(i);
	reaction_index[reaction] = -1;
	return;
}

void Initialize_reaction_times()
{
	int i;
	reaction_queue[0].reaction = 0;
	reaction_queue[0].reaction_time = -INFINITY;
	for (i = 1; i < MAX_REACTIONS; i++) {
		reaction_queue[i].reaction = 0;
		reaction_queue[i].reaction_time = INFINITY;
		reaction_index[i] = -1;
	}
	return;
}

void Assign_reaction_information()
{
	int i,j,k;
	for (i = 1; i <= MAX_PROMO; i++) {
		j=i;
		reactions[j].object_type = PROMOTER;
		reactions[j].object_index = i;
		reactions[j].reaction_id = PROMO_ON;
		j=PROMO_OFF_RXN(i);
		reactions[j].object_type = PROMOTER;
		reactions[j].object_index = i;
		reactions[j].reaction_id = PROMO_OFF;
		j=RNAP_LOAD_RXN(i);
		reactions[j].object_type = PROMOTER;
		reactions[j].object_index = i;
		reactions[j].reaction_id = LOAD_RNAP;
	}
	for (i = 1; i <= MAX_RNA; i++) {
		j=INIT_DEGR_RXN(i);
		reactions[j].object_type = RNA;
		reactions[j].object_index = i;
		reactions[j].reaction_id = START_DEGRADE;
		j=CONT_DEGR_RXN(i);
		reactions[j].object_type = RNA;
		reactions[j].object_index = i;
		reactions[j].reaction_id = CONT_DEGRADE;
		for (j = 1; j <= gene_num; j++) {
			k = RIBO_LOAD_RXN(i,j);
			reactions[k].object_type = RNA;
			reactions[k].object_index = i;
			reactions[k].gene_index = j;
			reactions[k].reaction_id = LOAD_RIBO;
		}
	}		
	for (i = 1; i <= MAX_RNAP; i++) {
		j=CLOSED_TO_OPEN_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = ACTIVATE_RNAP;
		j=UNBIND_RNAP_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = REMOVE_RNAP;
		j=RNAP_ON_F_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_F_TRANSL;
		j=RNAP_ON_B_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_B_TRANSL;
		j=RNAP_ADD_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_NT_ADD;
		j=RNAP_PAUSE_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_ENTER_PAUSE;
		j=RNAP_EXIT_RXN(i);		
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_EXIT_PAUSE;
		j=RNAP_OFF_F_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_OFF_FOR;
		j=RNAP_OFF_B_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_OFF_REV;
		j=RNAP_TERM_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_TERMINATION;
		j=RNAP_SURVEILLANCE_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_SURVEILLANCE;
		j=RNAP_EVAL_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_EVALUATION;
		j=RNAP_RELAX_RXN(i);
		reactions[j].object_type = RNAP;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RNAP_RELAX;
	}
	for (i = 1; i <= MAX_RIBO; i++) {
		j=RIBO_MOVE_RXN(i);
		reactions[j].object_type = RIBOSOME;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RIBO_TRANSLOC;
		j=RIBO_DECODE_RXN(i);
		reactions[j].object_type = RIBOSOME;
		reactions[j].object_index = i;
		reactions[j].reaction_id = RIBO_DECODE;
		j=CSAT_TERM_RXN(i);
		reactions[j].object_type = RIBOSOME;
		reactions[j].object_index = i;
		reactions[j].reaction_id = COLLISION_ABORT;
	}
	for (i = 1; i <= MAX_FISH_TIMES; i++) {
		j=ADD_SEQ_RXN(i);
		reactions[j].object_type = 0;
		reactions[j].object_index = i;
		reactions[j].reaction_id = ADD_SEQ_DATA;
		j=ADD_FISH_RXN(i);
		reactions[j].object_type = 0;
		reactions[j].object_index = i;
		reactions[j].reaction_id = ADD_FISH_DATA;
	}
	for (i = 1; i <= 100; i++) {
		j=NASCENT_PROTEIN_COUNT_RXN(i);
		reactions[j].object_type = 0;
		reactions[j].object_index = i;
		reactions[j].reaction_id = LOG_NASCENT_PROTEINS;
	}
	j=SNAPSHOT_RXN(1);
	reactions[j].object_type = 0;
	reactions[j].object_index = 1;
	reactions[j].reaction_id = TAKE_SNAPSHOT;
	for (i = 1; i <= 3600; i++) {
		j=UPDATE_RXN(i);
		reactions[j].object_type = 0;
		reactions[j].object_index = 1;
		reactions[j].reaction_id = UPDATE_CONCENTRATIONS;
	}
	for (i = 1; i <= MAX_PROMO; i++) {
		j=ADD_PROMO_RXN(i);
		reactions[j].object_type = PROMOTER;
		reactions[j].object_index = i;
		reactions[j].reaction_id = ADD_PROMOTER;
	}
	
	
	j=PROMO_KILL_RXN;
	reactions[j].object_type = 0;
	reactions[j].object_index = 1;
	reactions[j].reaction_id = KILL_PROMOTERS;

	j=TOPO_I_BINDING_RXN(1);
	reactions[j].object_type = 0;
	reactions[j].object_index = 1;
	reactions[j].reaction_id = TOPO_I_BINDING;
	
	j=GYRASE_BINDING_RXN(1);
	reactions[j].object_type = 0;
	reactions[j].object_index = 1;
	reactions[j].reaction_id = GYRASE_BINDING;

	for (i = 1; i <= MAX_TOPO; i++) {
		j=TOPO_UNBINDING_RXN(i);
		reactions[j].object_type = TOPOISOMERASE;
		reactions[j].object_index = i;
		reactions[j].reaction_id = TOPOISOMERASE_UNBINDING;
		j=TOPO_CLEAVAGE_RXN(i);
		reactions[j].object_type = TOPOISOMERASE;
		reactions[j].object_index = i;
		reactions[j].reaction_id = TOPOISOMERASE_CLEAVAGE;
		j=TOPO_SWITCH_STATE_RXN(i);
		reactions[j].object_type = TOPOISOMERASE;
		reactions[j].object_index = i;
		reactions[j].reaction_id = TOPOISOMERASE_STATE_CHANGE;
	}
	
        for (i = 1; i <= 10000; i++) {
                j = KYMO_RXN(i);
                reactions[j].object_type = 0;
                reactions[j].object_index = i;
                reactions[j].reaction_id = KYMOGRAPH;
        }

	
	return;
}

void Add_pseudoreactions_to_queue(int num_fish_times, float *fish_times,
                                  int num_seq_times, float *seq_times,
                                  float snap_time, float repl_time)
{
	int i;
	for (i = 1; i <= num_seq_times; i++) {
		Add_reaction_to_queue(ADD_SEQ_RXN(i),seq_times[i]);
	}
	for (i = 1; i <= num_fish_times; i++) {
		Add_reaction_to_queue(ADD_FISH_RXN(i),fish_times[i]);
	}
	Add_reaction_to_queue(SNAPSHOT_RXN(1),snap_time);

	for (i = 1; i < 24; i++) {
		Add_reaction_to_queue(NASCENT_PROTEIN_COUNT_RXN(i),
		   ((float) (i * 10)));
	}
	Add_reaction_to_queue(NASCENT_PROTEIN_COUNT_RXN(24),239.999);

	Add_reaction_to_queue(PROMO_KILL_RXN,promo_kill);

	return;
}
	
void Add_kymograph_reactions(double time)
{
        int i;
        for (i = 1; i <= num_kymo_times; i++) {
                Add_reaction_to_queue(KYMO_RXN(i), time + kymo_times[i]);
        }
}

double Find_firing_time(double mu)
{
	return (gsl_ran_exponential(r,mu));
}
