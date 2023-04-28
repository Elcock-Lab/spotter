#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "INCL/RNAP_rate_generator.h"
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"


int main (int argc, char *argv[])
{
        int i,num_bp;

	int GENERATE_SEQUENCE=1;
	int USE_NET_SEQ=1;
	int FLAT_SEQ=0;
	int RANDOM_SEQ=0;
	int REVERSE_COMPLEMENT=0;
	int NORMALIZE_RATES=1;
	int tx_start=1,tx_stop=100,rand_seed=-1;
	int new_start,new_stop;
	int bubble_adj=0,hybrid_adj=0;
	double mean_elong=30.0;
	char flat_repeat[100],param_file[500],seed_label[500];
	char scrunch[200],window[200];
	char seq_file[500],dwell_file[500],rev_file[500];

	sprintf(flat_repeat,"N/A");
	sprintf(seed_label,"N/A");
	sprintf(rev_file,"N/A");
	sprintf(dwell_file,"NOT GIVEN; WILL USE NET-SEQ");
	sprintf(seq_file,"NOT GIVEN; WILL GENERATE");
	sprintf(param_file,
	  "dekker_derived_rnap_paramaters.kc_500.no_advanced_pauses");
	sprintf(scrunch,"no_scrunch");
	sprintf(window,"no_slide");
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i],"-sequence") == 0) {
			i++;
			if (strcmp(argv[i],"random") == 0) {
			   RANDOM_SEQ++;
			   if (rand_seed < 0) {
			      sprintf(seed_label,
				"NOT GIVEN; WILL USE SYSTEM TIME\n");
			      strcat(seed_label,
				"***SAVE aligned_sequence.txt AND USE\n");
			      strcat(seed_label,
				"***START AMD STOP FOR THE SEQUENCE\n");
			      strcat(seed_label,
				"***GENERATED (BELOW) TO RERUN");
			   }
			   strcat(seq_file," A RANDOM SEQUENCE");
			}
			else if (strcmp(argv[i],"flat") == 0) {
				FLAT_SEQ++;
				sprintf(flat_repeat,"A (DEFAULT)");
				strcat(seq_file," A FLAT SEQUENCE");
			}
			else {
				GENERATE_SEQUENCE--;
				strcpy(seq_file,argv[i]);
			}
		}
		if (strcmp(argv[i],"-dwell") == 0) {
			USE_NET_SEQ--;
			NORMALIZE_RATES--;
			i++;
			strcpy(dwell_file,argv[i]);
		}
		if (strcmp(argv[i],"-rate") == 0) {
			NORMALIZE_RATES++;
			i++;
			mean_elong = strtod(argv[i],NULL);
		}
		if (strcmp(argv[i],"-reverse") == 0) {
			REVERSE_COMPLEMENT++;
			if (GENERATE_SEQUENCE) {
				GENERATE_SEQUENCE--;
			}
			i++;
			strcpy(rev_file,argv[i]);
		}
		if (strcmp(argv[i],"-tx_start") == 0) {
			i++;
			tx_start = strtod(argv[i],NULL);
		}
		if (strcmp(argv[i],"-tx_stop") == 0) {
			i++;
			tx_stop = strtod(argv[i],NULL);
		}
		if (strcmp(argv[i],"-seed") == 0) {
			i++;
			rand_seed = strtod(argv[i],NULL);
			sprintf(seed_label,"%d",rand_seed);
		}
		if (strcmp(argv[i],"-repeat") == 0) {
			i++;
			strcpy(flat_repeat,argv[i]);
		}
		if (strcmp(argv[i],"-param") == 0) {
			i++;
			strcpy(param_file,argv[i]);
		}
		if (strcmp(argv[i],"-bubble") == 0) {
			i++;
			bubble_adj = strtod(argv[i],NULL) - 12;
		}
		if (strcmp(argv[i],"-hybrid") == 0) {
			i++;
			hybrid_adj = strtod(argv[i],NULL) - 9;
		}
		if (strcmp(argv[i],"-scrunch") == 0) {
			sprintf(scrunch,"scrunch");
		}
		if (strcmp(argv[i],"-slide") == 0) {
			sprintf(window,"slide");
		}
	}

	char status_check[4][20];
	sprintf(status_check[0],"OFF");
	sprintf(status_check[1],"ON");
	sprintf(status_check[2],"ON");
	sprintf(status_check[3],"ON");
	printf("GENERATE SEQUENCE IS %s\n",status_check[GENERATE_SEQUENCE]);
	printf("RANDOM SEQUENCE IS %s\n",status_check[RANDOM_SEQ]);
	printf("FLAT SEQUENCE IS %s\n",status_check[FLAT_SEQ]);
	printf("SEQUENCE FILE: %s\n",seq_file);
	printf("USE NET-SEQ IS %s\n",status_check[USE_NET_SEQ]);
	printf("DWELL FILE: %s\n",dwell_file);
	printf("REVERSE COMPLEMENT IS %s\n",status_check[REVERSE_COMPLEMENT]);
	printf("REVERSE FILE: %s\n",rev_file);
	printf("NORMALIZE RATES IS %s\n",status_check[NORMALIZE_RATES]);
	printf("MEAN ELONG RATE IS %f NT/S\n",mean_elong);
	printf("TX START IS %d\n",tx_start);
	printf("TX STOP IS %d\n",tx_stop);
	printf("RANDOM SEED IS %s\n",seed_label);
	printf("FLAT REPEAT CHARACTER IS %s\n",flat_repeat);
	printf("PARAMETER FILE: %s\n",param_file);

	if (GENERATE_SEQUENCE) {
		Generate_sequence_file(FLAT_SEQ,RANDOM_SEQ,flat_repeat,
					rand_seed,tx_start,tx_stop);
		sprintf(seq_file,"seq.request");
		tx_start += 19;
		tx_stop += 19;
		printf("WITHIN GENERATED SEQUENCE:\n");
		printf("TX START IS %d\n",tx_start);
		printf("TX STOP IS %d\n",tx_stop);
	}

	if (REVERSE_COMPLEMENT) {
		num_bp = Print_reverse_complement(rev_file);
		sprintf(seq_file,"rev_compl.seq");
		new_start = (tx_stop > tx_start ? tx_stop : tx_start);
		new_stop = (tx_start < tx_stop ? tx_start : tx_stop);
		tx_start = (num_bp - new_start) + 1;
		tx_stop = (num_bp - new_stop) + 1;
		printf("AFTER SEQUENCE REVERSAL:\n");
		printf("TX START IS %d\n",tx_start);
		printf("TX STOP IS %d\n",tx_stop);
	}

	if (!(bubble_adj)) {
		printf("USING DEFAULT BUBBLE SIZE OF 12 NT\n");
	}
	else {
		printf("USING USER-SPECIFIED BUBBLE SIZE OF: %d\n",
			12 + bubble_adj);
	}
	if (!(hybrid_adj)) {
		printf("USING DEFAULT HYBRID LENGTH OF 9 NT\n");
	}
	else {
		printf("USING USER-SPECIFIED HYBRID LENGTH OF: %d\n",
			9 + hybrid_adj);
	}

	Assign_hybridization_energies(seq_file, tx_start, tx_stop,
				      bubble_adj, hybrid_adj,
				      scrunch, window);

	if (USE_NET_SEQ) {
		Assign_NETseq_dwelltimes("aligned_sequence.txt");
		sprintf(dwell_file,"NET_seq_derived_dwelltimes.txt");
	}

	if (NORMALIZE_RATES) {
		Normalize_rates(dwell_file, mean_elong);
		sprintf(dwell_file,"normalized_dwelltime.seq");
	}

	Calculate_final_rates(
		param_file,
		"bubble_plus_hybrid_energies.plus_rates.forward.seq.txt",
		dwell_file);

	return(0);
}

