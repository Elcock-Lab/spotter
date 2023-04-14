#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "INCL/sim_types.h"
#include "INCL/io.h"
#include "INCL/traj.h"

int main (int argc, char *argv[])
{
	int i;
	int start_seed;					// Seed for RNG
	float  sim_start, sim_stop;			// Trajectory start and stop times in seconds
	int coord[2];					// Genomic position to tx unit
	float  generation_time;				// Doubling time in simulations
	float  b_time, c_time, d_time;			// Cell-cycle paramters
	int num_traj;					// Number of trajectories to be run
	float  seq_window[2];				// Beginning and end of interval for generating "seq" data
	float  sample_window[2];			// Beginning and end of interval for <mRNA>, <protein_ct>, etc.
	float  snap_time;				// Time "snapshot" of DNA, mRNA, and ribosomes is taken
	float  time_step;				// Trajectory time-step

	rnap_rate_form *rnap_dwell;			// Position-specific RNAP rate info (see header for form)
	ribo_rate_form *ribo_dwell;			// Position-specific ribosome rate info (see header)

	double mu_rnap_p2=INFINITY;			// (Currently) invariant rates for RNAP rxns:				
	double mu_rnap_p3=INFINITY;			//   p2,p3: rates of entry into add'l offpath states
	double mu_rnap_e1=0.565;			//   e1,e2,e2: rates of exit from offpath states
	double mu_rnap_e2=4.149;			//   kN,kd: Rate of nucleotide assoc/dissoc
	double mu_rnap_e3=100.0;			//   kc: Rate of nucleotide add'n
	double mu_rnap_N;				// Rates are those in Janissen-Dekker Biorxiv paper
	double mu_rnap_c=0.002;
	double mu_rnap_d;

	float t_promo_kill;				// Time at which promoters are shut off

	int min_effector;				// Determines which RNAs states upstream effects
	int rnap_antipause=0;				// RNAP capable of freeing upstream P1 from pause (0/1) 
	int ribo_antipause=0;				// Ribo capable of freeing upstream P1 from pause (0/1)
	int ribo_rnap_range;				// Spacing between ribo/RNAP req'd for mvmnt/effect
	int max_backtrack;				// Largest number of nt's RNAP may backtrack

	int tx_length;					// Length of transcription unit (determined from dwelltime files)
	float  k_on;					// k_on for promoter availability
	float  k_off;					// k_off for promoter availability
	float  fract;					// Fraction of time promoter is on
	float  k_loading;				// Rate of RNAP binding at transcription start_site
	float  k_unloading;				// Rate of RNAP unbinding at transcription start site
	float  k_to_open;				// Rate of closed-to-open isomerization for bound RNAP
	float  RNA_lifetime;				// Mean lifetime of mRNA (sec)
	float  k_continue_degrade;			// Rate at which degradation proceeds post-initiaion(bases/sec)
	float  prot_lifetime;				// Mean lifetime of protein (sec)
	float  mean_ribo_rate;				// Used for normalizing riboseq by gene
	float  min_dwell;				// The minimum allowable dwelltime(sec)
	float  p_rnap_rnap_push;			// Probability of upstream RNAP's pushing downstream RNAP 
	float  p_ribo_rnap_push;			// Probability of ribosome pushing downstream RNAP
	float  p_ribo_ribo_push;			// Probability of ribosome of pushing downstream ribosome
	float  k_CSAT;					// Rate of tsl termination for downstream ribosome in contact
	float  birth_vol;				// Size of cell at birth(fL)
	float  mean_RNAP_conc;				// Concentration of RNAPs assumed in k_loading listed
	float  mean_ribo_conc;				// Concentration of ribosomes assumed in k_init's listed
	int stable_RNA;					// Indicates whether mRNA(0) or t/rRNA(1)
	int tsl_yn;
	int read_state;					// Indicates whether simulation begin anew (0) or read state(1)
	int pol_width;					// Width (in basepairs) of RNAP
	int ribo_width;					// Width (in bases) of ribosome
	int degrade_width;				// Width (in bases) assigned to degradation apparatus
	int gene_check,num_genes=0;			// Number of independently-translated genes on transcript
	gene_form gene_set[50];				// Parameters connected with each gene on the tx unit (see io.h)
	float  time_set[20], param_set[25]; 		// Arrays for collecting information from the input file
	char name_set[7][200];				// File names to be used (name[0] is an arbitrary identifier,
	char input_name[200];				// and name[1] and [2] contain names for dwelltime files
	char state_file[2][200];			// Contain DNA/RNAP info([0]) and RNA/ribosome info([1])
	int test_match;
	int traj_index;
	int unit;
	int genecheck=0,gene_id[25]={0};
	int onefile_only=0;

	char singlefile[200];
	char txfile[200],tslfile[200],superfile[200];
	char sysfile[200],regfile[200],checklist[10][200];	
	char genefile[200],namecheck[10][200];

	sprintf(txfile,"NOT SUPPLIED");
	sprintf(tslfile,"NOT SUPPLIED");
	sprintf(superfile,"NOT SUPPLIED");
	sprintf(regfile,"NOT SUPPLIED");
	sprintf(sysfile,"NOT SUPPLIED");
	sprintf(singlefile,"NOT SUPPLIED");
	sprintf(namecheck[1],"TRANSCRIPTION");
	sprintf(namecheck[2],"TRANSLATION");
	sprintf(namecheck[3],"SUPERCOILING");
	sprintf(namecheck[4],"REGULATION");
	sprintf(namecheck[5],"SYSTEM");
	num_traj = 1;
	traj_index = 1;
	start_seed = 1234;	
	t_promo_kill = 99999.9;
	param_set[9]  = -0.01;
	param_set[10] = -0.01;
	param_set[11] = -0.01;

        for (i = 1; i < argc; i++) {
                if (strcmp(argv[i],"-txfile") == 0) {
                        i++;
			strcpy(txfile,argv[i]);
		}
                if (strcmp(argv[i],"-tslfile") == 0) {
                        i++;
			strcpy(tslfile,argv[i]);
		}
                if (strcmp(argv[i],"-superfile") == 0) {
                        i++;
			strcpy(superfile,argv[i]);
		}
                if (strcmp(argv[i],"-sysfile") == 0) {
                        i++;
			strcpy(sysfile,argv[i]);
		}
                if (strcmp(argv[i],"-regfile") == 0) {
                        i++;
			strcpy(regfile,argv[i]);
		}
                if (strcmp(argv[i],"-singlefile") == 0) {
                        i++;
			onefile_only++;
			strcpy(singlefile,argv[i]);
		}
                if (strcmp(argv[i],"-seed") == 0) {
                        i++;
			start_seed = strtol(argv[i],NULL,10);
		}
                if (strcmp(argv[i],"-shut") == 0) {
                        i++;
			t_promo_kill = strtod(argv[i],NULL);
		}
                if (strcmp(argv[i],"-traj") == 0) {
                        i++;
			traj_index = strtol(argv[i],NULL,10);
		}
                if (strcmp(argv[i],"-genefile") == 0) {
                        i++;
			strcpy(genefile,argv[i]);
			genecheck++;
		}
	}
	if (!(onefile_only)) {
	 strcpy(checklist[1],txfile);
	 strcpy(checklist[2],tslfile);
	 strcpy(checklist[3],superfile);
	 strcpy(checklist[4],regfile);
	 strcpy(checklist[5],sysfile);
	}
	else {
	 printf("USING A SINGLE INPUT FILE...\n");
	 printf("...MAKE SURE ALL PARAMETERS FOR ALL MODULES "
		"ARE INCLUDED IN THIS FILE!!!\n");
	 strcpy(checklist[1],singlefile);
	 strcpy(checklist[2],singlefile);
	 strcpy(checklist[3],singlefile);
	 strcpy(checklist[4],singlefile);
	 strcpy(checklist[5],singlefile);
	}
	strcpy(checklist[6],singlefile);

	printf("INPUT FILES:\n");
	for (i = 1; i <= 5; i++) {
	 printf("%s FILE: %s\n",namecheck[i],checklist[i]);
	 if (strcmp(singlefile,"NOT SUPPLIED") == 0 || i < 2) {
	  if (strcmp(checklist[i],"NOT SUPPLIED") != 0) {
	   gene_check = Read_basic_parameters(checklist[i],name_set,
					      time_set,param_set,
                                              gene_set,coord,&stable_RNA,
					      &read_state, &tsl_yn);
	   num_genes = (gene_check > num_genes ? gene_check : num_genes);
	  }
	  else {
	   if (i == 1 || i == 2 || i == 5) {
		printf("...%s FILE REQUIRED; CANNOT PROCEED...QUITTING!\n",
			namecheck[i]);
		return(0);
	   }
	   else {
		printf("...WILL SUPPLY DEFAULT PARAMETERS FOR SIMULATION");
		printf(" WITHOUT %s\n", namecheck[i]);
		if (i == 4) {
			Prepare_for_no_regulation_simulation_without_file(
				param_set);
		}
	   }
	  }
	 }
	}
	printf("\n");

	for (i = 1; i <= 10; i++) {
		time_set[i] *= 60.0;			// Convert all times into seconds
	}
	generation_time = time_set[1];
	b_time = time_set[2];
	c_time = time_set[3];
	d_time = time_set[4];
	sim_start = 0.0;	
	sim_stop = time_set[5];
	seq_window[0] = time_set[6];
	seq_window[1] = time_set[7];
	sample_window[0] = time_set[8];
	sample_window[1] = time_set[9];
	snap_time = time_set[10];
	time_step = time_set[11];

	k_on = param_set[1];
	k_off = param_set[2];
	fract = k_on/(k_on + k_off);
	k_loading = param_set[3];
	k_unloading = param_set[4];
	k_to_open = param_set[5];
	RNA_lifetime = param_set[6];
	k_continue_degrade = param_set[7];
	prot_lifetime = param_set[8];
	p_rnap_rnap_push = param_set[9];
	p_ribo_rnap_push = param_set[10];
	p_ribo_ribo_push = param_set[11];
	k_CSAT = param_set[12];
	pol_width = (int) (param_set[13] + 0.00001);
	ribo_width = (int) (param_set[14] + 0.00001);
	degrade_width = (int) (param_set[15] + 0.00001);
	mean_ribo_rate = param_set[16];
	min_dwell = param_set[17];
	birth_vol = param_set[18];
	mean_RNAP_conc = param_set[19];
	mean_ribo_conc = param_set[20];
	strcpy(state_file[0],name_set[4]);
	strcpy(state_file[1],name_set[5]);

	strcpy(input_name,name_set[2]);
	tx_length = Get_tx_length(input_name);

	rnap_dwell = malloc((tx_length+1)*sizeof(*rnap_dwell));
	ribo_dwell = malloc((tx_length+1)*sizeof(*ribo_dwell));

	tx_length = Read_RNAP_rates(input_name,rnap_dwell);

	strcpy(input_name,name_set[3]);

	test_match = Read_ribo_rates(input_name,ribo_dwell);

	if (test_match != tx_length) {
		printf("ERROR! LENGTH OF RNAP (%d) AND RIBOSOME (%d) DWELLTIME FILES DO NOT MATCH\n",
			tx_length,test_match);
		return(0);
	}

	Read_invariant_rates(&mu_rnap_p2,&mu_rnap_p3,&mu_rnap_e1,&mu_rnap_e2,
                             &mu_rnap_e3,&mu_rnap_N,&mu_rnap_c,&mu_rnap_d,
			     &min_effector,&rnap_antipause,&ribo_antipause,
			     &ribo_rnap_range,&max_backtrack,checklist[1]);

	Print_assigned_values(sim_start,sim_stop,seq_window,sample_window,k_loading,k_on,k_off,RNA_lifetime,
                              prot_lifetime,num_genes,gene_set,tx_length,name_set,snap_time, fract, num_traj,
                              generation_time,b_time,c_time,d_time,k_unloading,k_to_open,k_continue_degrade,
                              p_rnap_rnap_push,p_ribo_rnap_push,p_ribo_ribo_push,k_CSAT,pol_width,ribo_width,
                              degrade_width,mean_ribo_rate,min_dwell,coord,birth_vol,stable_RNA,read_state,
                              mean_RNAP_conc,mean_ribo_conc,state_file,tsl_yn);

        gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r,start_seed);

	sscanf(name_set[1],"iter_%d",&unit);
	if (genecheck) {
	  Read_gene_identities(unit,gene_id);
	}

	strcpy(input_name, name_set[1]);
	i = Generate_trajectory(r, num_traj, sim_stop, fract, k_on, k_off, k_loading, rnap_dwell, ribo_dwell,
                                snap_time, tx_length, seq_window, sample_window, gene_set, RNA_lifetime,
                                input_name, start_seed, num_genes, traj_index, time_step, k_continue_degrade,
                                p_rnap_rnap_push, p_ribo_rnap_push, p_ribo_ribo_push, pol_width, ribo_width,
                                degrade_width,mean_ribo_rate,min_dwell,generation_time,b_time,c_time,d_time,
                                k_unloading,k_to_open,k_CSAT,mean_RNAP_conc,mean_ribo_conc,coord,stable_RNA,
                                read_state,state_file,birth_vol,tsl_yn,gene_id,mu_rnap_p2,mu_rnap_p3,
				mu_rnap_e1,mu_rnap_e2,mu_rnap_e3,mu_rnap_N,mu_rnap_c,mu_rnap_d,min_effector,
				rnap_antipause,ribo_antipause,ribo_rnap_range,max_backtrack,t_promo_kill,
				checklist);

	free(rnap_dwell);
	free(ribo_dwell);

	return(0);
}
