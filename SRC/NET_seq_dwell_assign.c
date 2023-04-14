#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/RNAP_rate_generator.h"
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

#define MAX_BP 5000000
#define GENOME_SIZE 4639675


void Read_NET_seq_energy_function(double energies[][4])
{
	energies[0][0]  = -0.07501;
	energies[0][1]  = -0.07864;						// Values hard-coded
	energies[0][2]  = -0.10185;						// to eliminate file
	energies[0][3]  =  0.25550;						// reading

	energies[1][0]  = -0.13312;
	energies[1][1]  = -0.16128;
	energies[1][2]  = -0.04832;
	energies[1][3]  =  0.34272;

	energies[2][0]  = -0.04656;
	energies[2][1]  = -0.13685;
	energies[2][2]  =  0.11128;
	energies[2][3]  =  0.07213;

	energies[3][0]  =  0.14590;
	energies[3][1]  = -0.16868;
	energies[3][2]  =  0.03490;
	energies[3][3]  = -0.01212;

	energies[4][0]  = -0.00981;
	energies[4][1]  =  0.06710;
	energies[4][2]  = -0.04808;
	energies[4][3]  = -0.00921;

	energies[5][0]  =  0.12656;
	energies[5][1]  =  0.00473;
	energies[5][2]  = -0.11642;
	energies[5][3]  = -0.01487;

	energies[6][0]  =  0.16421;
	energies[6][1]  = -0.14502;
	energies[6][2]  = -0.14727;
	energies[6][3]  =  0.12808;

	energies[7][0]  = -0.07833;
	energies[7][1]  =  0.12342;
	energies[7][2]  =  0.04609;
	energies[7][3]  = -0.09118;

	energies[8][0]  = -0.11130;
	energies[8][1]  =  0.28267;
	energies[8][2]  = -0.06189;
	energies[8][3]  = -0.10947;

	energies[9][0]  = -0.10574;
	energies[9][1]  = -0.09055;
	energies[9][2]  = -0.11474;
	energies[9][3]  =  0.31103;

	energies[10][0] = -0.25472;
	energies[10][1] =  0.26130;
	energies[10][2] =  0.39936;
	energies[10][3] = -0.40594;

	energies[11][0] = -0.33028;
	energies[11][1] = -0.14111;
	energies[11][2] = -0.00640;
	energies[11][3] =  0.47778;

	energies[12][0] = -0.23531;
	energies[12][1] = -0.02632;
	energies[12][2] =  0.14515;
	energies[12][3] =  0.11647;

	energies[13][0] = -0.15616;
	energies[13][1] = -0.09431;
	energies[13][2] =  0.24282;
	energies[13][3] =  0.00765;

	energies[14][0] = -0.14802;
	energies[14][1] =  0.10002;
	energies[14][2] =  0.14771;
	energies[14][3] = -0.09971;

	energies[15][0] =  0.05242;
	energies[15][1] = -0.07661;
	energies[15][2] = -0.08774;
	energies[15][3] =  0.11192;

	return;
}

void Assign_relative_dwelltimes(char *dna_seq, int *seq_range,
				double energies[][4], int aligned_start,
				int aligned_stop)
{
	int i,j,look,ctr;
	double energy = 0.0;

	FILE *fp;
	fp = fopen("NET_seq_derived_dwelltimes.txt","w");
       	for (i = aligned_start; i <= aligned_stop; i++) {
                energy = 0.0;
                ctr = -1;
                for (j = seq_range[0]; j <= seq_range[1]; j++) {
                        ctr++;
                        look = i + j;
                        if (dna_seq[look] == 'A') {
                                energy += energies[ctr][0];
                        }
                        else if (dna_seq[look] == 'T') {
                                energy += energies[ctr][1];
                        }
                        else if (dna_seq[look] == 'C') {
                                energy += energies[ctr][2];
                        }
                        else {
                                energy += energies[ctr][3];
                        }
                }
		fprintf(fp,"%10d%10.5f\n",i - 12, exp(energy));
	}
	fclose(fp);
	return;
}

void Assign_NETseq_dwelltimes(char *seq_name)
{
        int num_bp;
	int seq_range[2];
	int aligned_start,aligned_stop;
	double energies[20][4];
        char dna_seq[MAX_BP];

        printf("READING SEQUENCE: %s\n", seq_name);
        num_bp = Read_DNA_sequence(seq_name, dna_seq);
        printf("FOUND %d BP IN SEQUENCE PROVIDED\n", num_bp);
	
	aligned_start = 13;
	aligned_stop  = num_bp - 5;

	Read_NET_seq_energy_function(energies);

	seq_range[0] = -10;
	seq_range[1] = 5;
	
	Assign_relative_dwelltimes(dna_seq, seq_range, energies,
				   aligned_start, aligned_stop);

	return;
}

