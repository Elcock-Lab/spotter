#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "INCL/ribosome_rate_generator.h"
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

int main (int argc, char *argv[])
{
        int i,syst_ret;
        char command[1000];

	int STANDARD_SEQUENCE=1;
	int GENEFILE_PROVIDED=0;
	int DWELL_FROM_DECODE=1;
	int FLAT_DWELL=0;
	int FLAT_DECODE=0;
	double flat_decode_time=100.0;
	char seq_file[500],dwell_file[500],gene_file[500],decode_file[200];

	sprintf(seq_file,"aligned_sequence.txt");
	sprintf(decode_file,"fluitt_based_transloc_decoding_split.txt");
	sprintf(dwell_file,"make_flat_dwell");

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i],"-sequence") == 0) {
			i++;
			strcpy(seq_file,argv[i]);
			if (strcmp(argv[i],"aligned_sequence.txt") != 0) {
				STANDARD_SEQUENCE--;
			}
		}
		if (strcmp(argv[i],"-genefile") == 0) {
			i++;
			GENEFILE_PROVIDED++;
			strcpy(gene_file,argv[i]);
		}
		if (strcmp(argv[i],"-dwell") == 0) {
			i++;
			if (strcmp(argv[i],"flat") == 0) {
				FLAT_DWELL++;
				DWELL_FROM_DECODE--;
			}
			else if (strcmp(argv[i],"use_decode") == 0) {
				DWELL_FROM_DECODE++;
			}
			else {
				DWELL_FROM_DECODE--;
				strcpy(dwell_file,argv[i]);
			}
		}
		if (strcmp(argv[i],"-decode") == 0) {
			i++;
			if (strcmp(argv[i],"flat") == 0) {
				FLAT_DECODE++;
				strcpy(decode_file,"make_flat_decode");
			}
			else {
				strcpy(decode_file,argv[i]);
			}
		}
		if (strcmp(argv[i],"-dt") == 0) {
			i++;
			flat_decode_time = strtod(argv[i],NULL);
		}
	}

	if (!(GENEFILE_PROVIDED)) {
		printf("CAN'T PROCEED WITHOUT GENE BOUNDARY FILE\n");
		printf("QUITTING...\n");
		return(0);
	}
	printf("GENE FILE: %s\n",gene_file);

	if (STANDARD_SEQUENCE) {
		printf("USING ALIGNED OUTPUT FILE FROM RNAP FILEMAKING\n");
		Trim_sequence("aligned_sequence.txt",12,5);
		strcpy(seq_file,"trimmed.seq");
	}
	else {
		printf("USING SEQUENCE FILE %s; ",seq_file);
		printf("CONFIRM THAT SEQUENCE BEGINS AT TSS "
			"AND STOPS AT TX END\n");
	}

	printf("FILE OR ACTION FOR DWELL TIMES: %s\n",dwell_file);
	printf("FILE OR ACTION FOR DECODING TIMES: %s\n",decode_file);

	Add_codon_based_decoding_info(seq_file, gene_file, dwell_file,
				      decode_file,0,0, flat_decode_time);

	if (DWELL_FROM_DECODE) {
	   printf("ASSIGNING DWELL TIMES BASED ON DECODING RATES...\n");
	   Assign_decode_only(gene_file,"ribo_rates.seq.full_info");
	   sprintf(command,"cp decoding_only.ribo_rates ribo_rates.seq.full_info");
	   syst_ret = system(command);
	   if (syst_ret == -1) {
		printf("FAILED TO EXECUTE %s\n",command);
	   }
	}

	return(0);
}

