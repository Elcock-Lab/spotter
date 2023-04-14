#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/RNAP_rate_generator.h"
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

#define MAX_BP 5000000
#define GENOME_SIZE 4639675

typedef struct {
	char seq[5];
	double dg;
} energy_form;

energy_form dna_dna_nearest[17];
energy_form rna_dna_nearest[17]; 


void Align_sequence_with_TSS(char *orig_seq, char *final_seq,
			     int num_bp, int start_pos, int stop_pos,
			     int *aligned_start, int *aligned_stop)
{
	int i,ctr=0;
	char padder[5];
	sprintf(padder,"ATCG");
	FILE *fp;
	fp = fopen("aligned_sequence.txt","w");
	for (i = 1; i <= (stop_pos - start_pos) + 13; i++) {
		if ((start_pos + i) - 13 < 1) {					// Pad if sequence provided
			final_seq[i] = padder[i % 4];				// does not cover region 
		}								// upstream of TSS; otherwise
		else {								// use sequence
			final_seq[i] = orig_seq[(start_pos + i) - 13];
		}
		ctr++;
		fprintf(fp,"%c",final_seq[i]);
		if (ctr % 70 == 0) {
			fprintf(fp,"\n");
		}
	}
	for (i = 1; i <= 5; i++) {
		if (num_bp - stop_pos >= i) {
			final_seq[ctr + i] = orig_seq[stop_pos + i];
		}
		else {
			final_seq[ctr + i] = padder[i % 4];
		}
		fprintf(fp,"%c",final_seq[ctr + i]);
		if ((ctr + i) % 70 == 0) {
			fprintf(fp,"\n");
		}
	}
	final_seq[ctr + 6] = '\0';
	*aligned_start = 13;
	*aligned_stop  = 13 + (stop_pos - start_pos);
	fclose(fp);
	return;
}
			
			
void Read_DNA_DNA_nearest_neighbor_energies(energy_form *dna_dna_nearest,
                                            double *dna_init_GC,
					    double *dna_init_AT)
{
	int ctr=0,post_ctr=10;

	dna_dna_nearest[1].dg = -1.00;						// Avoid read-in for
	sprintf(dna_dna_nearest[1].seq,"AA");					// portability
	dna_dna_nearest[2].dg = -0.88;
	sprintf(dna_dna_nearest[2].seq,"AT");
	dna_dna_nearest[3].dg = -0.58;
	sprintf(dna_dna_nearest[3].seq,"TA");
	dna_dna_nearest[4].dg = -1.45;
	sprintf(dna_dna_nearest[4].seq,"CA");
	dna_dna_nearest[5].dg = -1.44;
	sprintf(dna_dna_nearest[5].seq,"GT");
	dna_dna_nearest[6].dg = -1.28;
	sprintf(dna_dna_nearest[6].seq,"CT");
	dna_dna_nearest[7].dg = -1.30;
	sprintf(dna_dna_nearest[7].seq,"GA");
	dna_dna_nearest[8].dg = -2.17;
	sprintf(dna_dna_nearest[8].seq,"CG");
	dna_dna_nearest[9].dg = -2.24;
	sprintf(dna_dna_nearest[9].seq,"GC");
	dna_dna_nearest[10].dg = -1.84;
	sprintf(dna_dna_nearest[10].seq,"GG");
	*dna_init_GC = 0.98;
	*dna_init_AT = 1.03;

	for (ctr = 1; ctr <= 10; ctr++) {					// Redundant entries
		if (strcmp(dna_dna_nearest[ctr].seq,"AA") == 0) {		// explicitly included
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"TT");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
		if (strcmp(dna_dna_nearest[ctr].seq,"CA") == 0) {
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"TG");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
		if (strcmp(dna_dna_nearest[ctr].seq,"GT") == 0) {
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"AC");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
		if (strcmp(dna_dna_nearest[ctr].seq,"CT") == 0) {
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"AG");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
		if (strcmp(dna_dna_nearest[ctr].seq,"GA") == 0) {
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"TC");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
		if (strcmp(dna_dna_nearest[ctr].seq,"GG") == 0) {
			post_ctr++;
			sprintf(dna_dna_nearest[post_ctr].seq,"CC");
			dna_dna_nearest[post_ctr].dg =
				dna_dna_nearest[ctr].dg;
		}
	}
	return;
}

void Read_RNA_DNA_nearest_neighbor_energies(energy_form *rna_dna_nearest,
                                            double *rna_dna_init)
{
	rna_dna_nearest[1].dg   = -1.0;
	sprintf(rna_dna_nearest[1].seq, "AA");
	rna_dna_nearest[2].dg   = -2.1;
	sprintf(rna_dna_nearest[2].seq, "AC");					// Given explicitly to
	rna_dna_nearest[3].dg   = -1.8;						// avoid file-reading;
	sprintf(rna_dna_nearest[3].seq, "AG");					// also converted to DNA
	rna_dna_nearest[4].dg   = -0.9;						// for ease of comparison
	sprintf(rna_dna_nearest[4].seq, "AT");
	rna_dna_nearest[5].dg   = -0.9;
	sprintf(rna_dna_nearest[5].seq, "CA");
	rna_dna_nearest[6].dg   = -2.1;
	sprintf(rna_dna_nearest[6].seq, "CC");
	rna_dna_nearest[7].dg   = -1.7;
	sprintf(rna_dna_nearest[7].seq, "CG");
	rna_dna_nearest[8].dg   = -0.9;
	sprintf(rna_dna_nearest[8].seq, "CT");
	rna_dna_nearest[9].dg   = -1.3;
	sprintf(rna_dna_nearest[9].seq, "GA");
	rna_dna_nearest[10].dg  = -2.7;
	sprintf(rna_dna_nearest[10].seq,"GC");
	rna_dna_nearest[11].dg  = -2.9;
	sprintf(rna_dna_nearest[11].seq,"GG");
	rna_dna_nearest[12].dg  = -1.1;
	sprintf(rna_dna_nearest[12].seq,"GT");
	rna_dna_nearest[13].dg  = -0.6;
	sprintf(rna_dna_nearest[13].seq,"TA");
	rna_dna_nearest[14].dg  = -1.5;
	sprintf(rna_dna_nearest[14].seq,"TC");
	rna_dna_nearest[15].dg  = -1.6;
	sprintf(rna_dna_nearest[15].seq,"TG");
	rna_dna_nearest[16].dg  = -0.2;
	sprintf(rna_dna_nearest[16].seq,"TT");

	*rna_dna_init = 3.1;

	return;
}

double Find_energy_parameter(char *seq, energy_form *nearest)
{
	int i;
	for (i = 1; i <= 16; i++) {
		if (seq[0] == nearest[i].seq[0] &&
		    seq[1] == nearest[i].seq[1]) {
			return(nearest[i].dg);
		}
	}
	return(0.0);
}

double Calculate_energies(int site, char *dna_seq, int dir,
			  double dna_init_GC, double dna_init_AT,
			  double rna_dna_init, double *transloc_energy,
			  double *rna_hybrid, int scrunch,
			  int sliding_window)
{
	int i,ctr = -1,upstream,start_hybr;
	double bubble_energy = 0.0, hybrid_energy = 0.0;
	double post_bubble =0.0;
	char sequence[25];
	*transloc_energy = 0.0;
	if (dir == 0) {
		if (scrunch) {
			upstream = site - 10 - scrunch;
		}
		else {
			upstream = site - 11;
		}
		for (i = upstream - 1; i <= site + 2; i++) {
			ctr++;
			sequence[ctr] = dna_seq[Wraparound(i)];
		}
		ctr--;
	}
	else {
		for (i = site + 11; i >= site - 1; i--) {
			ctr++;
			if (dna_seq[Wraparound(i)] == 'A') {
				sequence[ctr] = 'T';
			}
			else if (dna_seq[Wraparound(i)] == 'T') {
				sequence[ctr] = 'A';
			}
			else if (dna_seq[Wraparound(i)] == 'G') {
				sequence[ctr] = 'C';
			}
			else {
				sequence[ctr] = 'G';
			}
		}
	}
	if (sequence[0] == 'C' || sequence[0] == 'G') {
		bubble_energy -= dna_init_GC;
		if (scrunch || !(sliding_window)) {
			post_bubble -= dna_init_GC;
		}
	}
	else {
		bubble_energy -= dna_init_AT;
		if (scrunch || !(sliding_window)) {
			post_bubble -= dna_init_AT;
		}
	}
	if (!(scrunch) && sliding_window) {
		if (sequence[1] == 'C' || sequence[1] == 'G') {
			post_bubble -= dna_init_GC;
		}
	
		else {
			post_bubble -= dna_init_AT;
		}
	}
	if (sequence[ctr] == 'C' || sequence[ctr] == 'G') {
		bubble_energy -= dna_init_GC;
	}
	else {
		bubble_energy -= dna_init_AT;
	}
	if (sequence[ctr+1] == 'C' || sequence[ctr+1] == 'G') {
		post_bubble -= dna_init_GC;
	}
	else {
		post_bubble -= dna_init_AT;
	}
	hybrid_energy += rna_dna_init;
	*transloc_energy = hybrid_energy;
	char *str = sequence;
	for (i = 1; i <= ctr-1; i++) {
		str++;
		bubble_energy += 
			   Find_energy_parameter(str,dna_dna_nearest);
		if (scrunch || !(sliding_window) || i >= 2) {
			post_bubble +=
			   	Find_energy_parameter(str,dna_dna_nearest);
		}
		if (scrunch) {
			start_hybr = 11;
		}
		else {
			start_hybr = 3;
		}
		if (i >= start_hybr && i < ctr - 1 &&
		    (!(scrunch) || scrunch > 1)) {
			hybrid_energy +=
				Find_energy_parameter(str,rna_dna_nearest);
			if (i > 3 || (!sliding_window)) {
				*transloc_energy +=
					Find_energy_parameter(
						str,rna_dna_nearest);
			}
		}
	}
	str++;
	post_bubble += Find_energy_parameter(str,dna_dna_nearest);
	*rna_hybrid = hybrid_energy;
	bubble_energy *= -1.0;
	post_bubble *= -1.0;
	bubble_energy += (1.0 * ((double) scrunch));
	post_bubble += ((scrunch > 0) * ((double) (scrunch + 1)));
	*transloc_energy = *transloc_energy + post_bubble;
	return(hybrid_energy + bubble_energy);
}	

void Write_energies(int aligned_start, int aligned_stop, 
		    double *energy, double *pre, double *post,
		    int dir)
{
	int i,j,offset;
	double barrier;
	double off_energy[4], on_energy[4];
	double off_step_f,off_step_b,on_step_f,on_step_b;
	char name[200];
	if (dir == 0) {
	   sprintf(name,
	      "bubble_plus_hybrid_energies.plus_rates.forward.seq.txt");
	   offset = 1;
	}
	else {
	   sprintf(name,
	      "bubble_plus_hybrid_energies.plus_rates.reverse.seq.txt");
	   offset = -1;
	}
	FILE *fp;
	fp = fopen(name,"w");
	for (i = aligned_start; i <= aligned_stop; i++) {
		off_energy[1] = energy[Wraparound(i - offset)];
		off_energy[2] = energy[i];
		off_energy[3] = energy[Wraparound(i + offset)];
		if (i <= aligned_start + 0) {
			off_energy[1] = energy[i];
			off_energy[2] = energy[i];
			off_energy[3] = energy[i];
		}	
		on_energy[1]  = pre[i];
		on_energy[2]  = post[i];
		for (j = 1; j <= 3; j++) {
			off_energy[j] *= 1.62333;				// Convert to kbT at 310K
			on_energy[j]  *= 1.62333;
		}
		barrier = (off_energy[2] + off_energy[3])/2.0;
		barrier += 5.0;
		off_step_f = 1000000.0 *
			(exp(-1.0 * (barrier - off_energy[2])));
		barrier = (off_energy[1] + off_energy[2])/2.0;
		barrier += 5.0;
		off_step_b = 1000000.0 *
			(exp(-1.0 * (barrier - off_energy[2])));
		barrier = (on_energy[1] + on_energy[2])/2.0;
		barrier += 5.0;
		on_step_f = 1000000.0 *
			(exp(-1.0 * (barrier - on_energy[1])));
		on_step_b = 1000000.0 *
			(exp(-1.0 * (barrier - on_energy[2])));

		fprintf(fp,"%10d%12.3f%12.3f%12.3f%12.3f"
			   "%12.3f%12.3f%12.3f\n",
			i - 12,energy[i],pre[i],post[i],on_step_f,
			on_step_b,off_step_f,off_step_b);
	}
	fclose(fp);
	return;
}	

void Assign_hybridization_energies(char *seq_name, int start_pos,
				   int stop_pos, char *scruncher,
				   char *slider)
{
        int i,downstream_gap,ctr;
	int num_bp;
	int aligned_start,aligned_stop;
	int SCRUNCH;
	int SLIDING_WINDOW;
	double transition_cf,diff;
	double rna_dna_init, dna_init_GC, dna_init_AT;
	double transloc_energy,rna_hybrid;
	double f_pre[MAX_BP];
	double f_post[MAX_BP];
	double forward[MAX_BP];							// 0(1) -> for(rev)
        char pre_seq[MAX_BP], dna_seq[MAX_BP];					// in bidirectional fnctns

        printf("READING SEQUENCE: %s\n", seq_name);
        num_bp = Read_DNA_sequence(seq_name, pre_seq);
        printf("FOUND %d BP IN SEQUENCE PROVIDED\n", num_bp);
        if (num_bp < stop_pos) {
                printf("ERROR: SEQUENCE SHORTER THAN EXPECTED LENGTH\n");
		return;
        }
	if (start_pos < 13) {
		printf("INCOMPLETE SEQUENCE FOR INITIAL %d BP(S):\n",
			13 - start_pos);
		printf("ADDING RANDOM BP UPSTREAM TO "
			"COVER THESE POSITIONS...\n");
		printf("STOP AND EXTEND SEQUENCE BY THIS AMOUNT "
			"IF THIS IS UNACCEPTABLE\n");
	}
	downstream_gap = 5 - (num_bp - stop_pos);
	if (downstream_gap > 0) {
		printf("INCOMPLETE SEQUENCE FOR LAST %d BP(S):\n",
			downstream_gap);
		printf("ADDING RANDOM BP DOWNSTREAM TO "
			"COVER THESE POSITIONS...\n");
		printf("STOP AND EXTEND SEQUENCE BY THIS AMOUNT "
			"IF THIS IS UNACCEPTABLE\n");
	}

	Align_sequence_with_TSS(pre_seq,dna_seq,
				num_bp,start_pos,stop_pos,
				&aligned_start,&aligned_stop);

	Read_DNA_DNA_nearest_neighbor_energies(dna_dna_nearest,
					       &dna_init_GC,&dna_init_AT);
	Read_RNA_DNA_nearest_neighbor_energies(rna_dna_nearest,
						&rna_dna_init);

	SCRUNCH = (scruncher[0] == 's' || scruncher[0] == 'S');
	if (SCRUNCH) {
		printf("USING DNA SCRUNCHING MODEL\n");
	}

	SLIDING_WINDOW = (slider[0] == 's' || slider[0] == 'S');
	if (SLIDING_WINDOW) {
		printf("USING CONSTANT-SIZE SLIDING BUBBLE\n");
	}
	
	ctr = 0;
	for (i = aligned_start - 1; i <= aligned_start + 7; i++) {
		forward[i] = Calculate_energies(i,dna_seq,0,dna_init_GC,
						dna_init_AT,rna_dna_init,
						&transloc_energy,
						&rna_hybrid,ctr,
						SLIDING_WINDOW);
		f_pre[i] = forward[i];
		f_post[i] = transloc_energy;
		if (SCRUNCH) {
			ctr++;
		}
	}
	for (i = aligned_start + 8; i <= aligned_stop + 1; i++) {
		forward[i] = Calculate_energies(i,dna_seq,0,dna_init_GC,
						dna_init_AT, rna_dna_init,
						&transloc_energy,
						&rna_hybrid,0,
						SLIDING_WINDOW);
		f_pre[i] = forward[i];
		f_post[i] = transloc_energy;
	}
	if (SCRUNCH) {
		transition_cf = Calculate_energies(aligned_start+8,
						   dna_seq,0,dna_init_GC,
						   dna_init_AT,
						   rna_dna_init,
						   &transloc_energy,
						   &rna_hybrid,9,
						   SLIDING_WINDOW);
		printf("SCRUNCHED/UNSCRUNCHED CF: %f %f\n",
			transition_cf,f_post[aligned_start+8]);
		diff = f_post[aligned_start+8] - transition_cf;
		for (i = aligned_start - 1; i <= aligned_start + 7; i++) {
			forward[i] += diff;
			f_pre[i] += diff;
			f_post[i] += diff;
		}
		forward[aligned_start+8] = f_post[aligned_start+8];
		f_pre[aligned_start+8] = f_post[aligned_start+8];
	}

	Write_energies(aligned_start,aligned_stop,forward,f_pre,f_post,0);

	return;
}

