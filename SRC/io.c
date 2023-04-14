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

#include "INCL/sim_types.h"
#include "INCL/gen_def.h"
#include "INCL/sim_settings.h"
#include "INCL/reaction_manager.h"
#include "INCL/replication.h"
#include "INCL/regulation.h"
#include "INCL/transcription.h"
#include "INCL/translation.h"
#include "INCL/supercoiling.h"
#include "INCL/topo.h"
#include "INCL/utilities.h"


int Read_basic_parameters(char *name,char name_set[][200],float *time_set,
			  float *param_set, gene_form *gene_set,
			  int *coord, int *stable, int *state_file,
			  int *tsl_yn)
{
	int ctr = 0,num_genes=0,frac_yn=0,which_gene,val;
	float fraction=0.0,value;
	char line[200], check[40],junk[200],on_off[100];
	
	FILE *fp;
	fp = fopen(name,"r");
        while (fgets(line, 200, fp) != NULL) {
		ctr++;
		sprintf(junk,"%.30s",line);
		sscanf(junk,"%s",check);
		char *str = line;
		str += 39;
		if (strcmp(check,"promoter_name") == 0) {
			sscanf(str,"%s",name_set[1]);
		}
		if (strcmp(check,"rnap_dwell_file") == 0) {
			sscanf(str,"%s",name_set[2]);
		}
		if (strcmp(check,"ribo_dwell_file") == 0) {
			sscanf(str,"%s",name_set[3]);
		}
		if (strcmp(check,"dna_state_file") == 0) {
			sscanf(str,"%s",name_set[4]);
		}
		if (strcmp(check,"rna_state_file") == 0) {
			sscanf(str,"%s",name_set[5]);
		}
		if (strcmp(check,"genomic_start(bp)") == 0) {
			sscanf(str,"%d",&coord[0]);
		}
		if (strcmp(check,"genomic_stop(bp)") == 0) {
			sscanf(str,"%d",&coord[1]);
		}
		if (strcmp(check,"stable_RNA(y/n)") == 0) {
			sscanf(str,"%s",on_off);
			*stable = (on_off[0] == 'y' || on_off[0] == 'Y');
		}	
		if (strcmp(check,"translated(y/n)") == 0) {
			sscanf(str,"%s",on_off);
			*tsl_yn = (on_off[0] == 'y' || on_off[0] == 'Y');
			*tsl_yn += (2 * (check[0] == 'h' ||
					 check[0] == 'H'));
		}	
		if (strcmp(check,"read_state(y/n)") == 0) {
			sscanf(str,"%s",on_off);
			*state_file = (on_off[0] == 'y' ||
				       on_off[0] == 'Y');
		}	
		if (strcmp(check,"generation_time(min)") == 0) {
			sscanf(str,"%f",&time_set[1]);
		}
		if (strcmp(check,"B_time(min)") == 0) {
			sscanf(str,"%f",&time_set[2]);
		}
		if (strcmp(check,"C_time(min)") == 0) {
			sscanf(str,"%f",&time_set[3]);
		}
		if (strcmp(check,"C_time(min)") == 0) {
			sscanf(str,"%f",&time_set[4]);
		}
		if (strcmp(check,"total_sim_time(min)") == 0) {
			sscanf(str,"%f",&time_set[5]);
		}
		if (strcmp(check,"seq_window_start(min)") == 0) {
			sscanf(str,"%f",&time_set[6]);
		}
		if (strcmp(check,"seq_window_stop(min)") == 0) {
			sscanf(str,"%f",&time_set[7]);
		}
		if (strcmp(check,"sample_window_start(min)") == 0) {
			sscanf(str,"%f",&time_set[8]);
		}
		if (strcmp(check,"sample_window_stop(min)") == 0) {
			sscanf(str,"%f",&time_set[9]);
		}
		if (strcmp(check,"snapshot_time(min)") == 0) {
			sscanf(str,"%f",&time_set[10]);
		}
		if (strcmp(check,"timestep(sec)") == 0) {
			sscanf(str,"%f",&time_set[11]);
		}
		if (strcmp(check,"k_on(promoter)") == 0) {
			sscanf(str,"%f",&param_set[1]);
		}
		if (strcmp(check,"use_frac_on(y/n)") == 0) {
			sscanf(str,"%s",on_off);
			frac_yn = (on_off[0] == 'y' || on_off[0] == 'Y');
		}
		if (strcmp(check,"frac_on") == 0) {
			sscanf(str,"%f",&fraction);
		}
		if (strcmp(check,"k_off(promoter)") == 0) {
			if (frac_yn) {						// Get k_off from fraction on if specified
				param_set[2] = ((1 - fraction) *
						param_set[1])/fraction;
			}
			else {
				sscanf(str,"%f",&param_set[2]);
			}
		}
		if (strcmp(check,"k_binding(RNAP)") == 0) {
			sscanf(str,"%f",&param_set[3]);
		}
		if (strcmp(check,"k_unbinding(RNAP)") == 0) {
			sscanf(str,"%f",&param_set[4]);
		}
		if (strcmp(check,"k_closed_to_open(RNAP)") == 0) {
			sscanf(str,"%f",&param_set[5]);
		}
		if (strcmp(check,"mean_mRNA_lifetime(s)") == 0) {
			sscanf(str,"%f",&param_set[6]);
		}
		if (strcmp(check,"k_continue_degr(bases/s)") == 0) {
			sscanf(str,"%f",&param_set[7]);
		}
		if (strcmp(check,"mean_protein_lifetime(s)") == 0) {
			sscanf(str,"%f",&param_set[8]);
		}
		if (strcmp(check,"p_rnap_rnap_push") == 0) {
			sscanf(str,"%f",&param_set[9]);
		}
		if (strcmp(check,"p_ribosome_rnap_push") == 0) {
			sscanf(str,"%f",&param_set[10]);
		}
		if (strcmp(check,"p_ribosome_ribosome_push") == 0) {
			sscanf(str,"%f",&param_set[11]);
		}
		if (strcmp(check,"k_ribo_coll_abort(/sec)") == 0) {
			sscanf(str,"%f",&param_set[12]);
		}
		if (strcmp(check,"RNAP_width(bp)") == 0) {
			sscanf(str,"%f",&param_set[13]);
		}
		if (strcmp(check,"ribosome_width(bases)") == 0) {
			sscanf(str,"%f",&param_set[14]);
		}
		if (strcmp(check,"degradation_width(bases)") == 0) {
			sscanf(str,"%f",&param_set[15]);
		}
		if (strcmp(check,"mean_tsl_rate(bases/sec)") == 0) {
			sscanf(str,"%f",&param_set[16]);
		}
		if (strcmp(check,"min_ribo_dwell(sec/base)") == 0) {
			sscanf(str,"%f",&param_set[17]);
		}
		if (strcmp(check,"birth_volume(fL)") == 0) {
			sscanf(str,"%f",&param_set[18]);
		}
		if (strcmp(check,"mean[RNAP](uM)") == 0) {
			sscanf(str,"%f",&param_set[19]);
		}
		if (strcmp(check,"mean[ribosome](uM)") == 0) {
			sscanf(str,"%f",&param_set[20]);
		}
		if (strcmp(check,"num_RBS") == 0) {
			sscanf(str,"%d",&num_genes);
		}
		if (strcmp(check,"RBS_loc(rel_to_tx_start)") == 0) {
			sscanf(str,"%d %d",&which_gene,&val);
			gene_set[which_gene].start = val;
		}
		if (strcmp(check,"tsl_stop(rel_to_tx_start)") == 0) {
			sscanf(str,"%d %d",&which_gene,&val);
			gene_set[which_gene].stop = val;
		}
		if (strcmp(check,"k_loading_RBS") == 0) {
			sscanf(str,"%d %f",&which_gene,&value);
			gene_set[which_gene].k_init = value;
		}
		if (strcmp(check,"mean_rate_tsl(bases/sec)") == 0) {
			sscanf(str,"%d %f",&which_gene,&value);
			gene_set[which_gene].mean_rate = value;
		}
	}
	fclose(fp);
	return(num_genes);
}		
	

int Read_master_input_file(char *name, char name_set[][200],
			   float *time_set, float *param_set,
                           gene_form *gene_set, int *coord, int *stable,
			   int *state_file, int *tsl_yn)
{
	int ctr = 0,ribo_ctr,num_genes;
	float  fraction=0.0;
	char line[200], check[10],junk[200];
	
	FILE *fp;
	fp = fopen(name,"r");
        while (fgets(line, 200, fp) != NULL) {
		ctr++;
		if (ctr < 4) {
			sscanf(line,"%s %s", junk, name_set[ctr]);
		}
		if (ctr >= 4 && ctr <= 5) {
			sscanf(line,"%s %d", junk, &coord[ctr % 2]);
		}
		if (ctr >= 6 && ctr <= 8) {
			sscanf(line,"%s %s", junk, check);
			if (ctr == 6) {
				*stable = (check[0] == 'y' ||
					   check[0] == 'Y');
			}
			else if (ctr == 7) {
				*tsl_yn = (check[0] == 'y' ||
					   check[0] == 'Y');
				*tsl_yn += (2*(check[0] == 'h' ||
					       check[0] == 'H'));
			}
			else {
				*state_file = (check[0] == 'y' ||
					       check[0] == 'Y');
			}
		}
		if (ctr >= 9 && ctr <= 10) {
			sscanf(line,"%s %s", junk, name_set[ctr-5]);
		}

		if (ctr >= 11 && ctr <= 21) {
			sscanf(line,"%s %f ", junk, &time_set[ctr - 10]);
		}

		if (ctr == 22) {
			sscanf(line,"%s %f ", junk, &param_set[ctr - 21]);
		}
		if (ctr == 23) {
			sscanf(line,"%s %s", junk, check);
		}
		if (ctr == 24) {
			sscanf(line,"%s %f ", junk, &fraction);
		}
		if (ctr == 25) {
			if (strcmp(check,"y") == 0) {					// Get k_off from fraction on if specified
				param_set[2] = ((1 - fraction) *
						param_set[1])/fraction;
			}
			else {
				sscanf(line,"%s %f ", junk, &param_set[2]);
			}
		}
		if (ctr >= 26 && ctr <= 43) {
			sscanf(line,"%s %f ", junk, &param_set[ctr - 23]);
		}
		if (ctr == 44) {
			sscanf(line,"%s %d", junk, &num_genes);
		}
		if (ctr > 44) {
			ribo_ctr = ((ctr - 45) / 4) + 1;
			if ((ctr - 45) % 4 == 0) {
				sscanf(line, "%s %d",
					junk, &gene_set[ribo_ctr].start);
			}
			if ((ctr - 45) % 4 == 1) {
				sscanf(line, "%s %d",
					junk, &gene_set[ribo_ctr].stop);
			}
			if ((ctr - 45) % 4 == 2) {
				sscanf(line, "%s %f ",
					junk, &gene_set[ribo_ctr].k_init);
			}
			if ((ctr - 45) % 4 == 3) {
				sscanf(line, "%s %f ",
					junk, &gene_set[ribo_ctr].mean_rate);
			}
		}
	}
	fclose(fp);
	return(num_genes);		
}

int Get_tx_length(char *name)
{
	int ctr = 0;
	char line[200];

	FILE *fp;
	fp = fopen(name, "r");
	while (fgets(line,200,fp) != NULL) {
		ctr++;
	}
	fclose(fp);
	return(ctr);
}

int Read_RNAP_rates(char *name, rnap_rate_form *rate)
{
	int ctr = 0,bp;
	char line[200];

	FILE *fp;
	fp = fopen(name, "r");
	while (fgets(line,200,fp) != NULL) {
		ctr++;
		sscanf(line,"%d %lf %lf %lf %lf %lf %lf", &bp,
			&rate[ctr].net, &rate[ctr].mu_f,
			&rate[ctr].mu_b, &rate[ctr].mu_p,
			&rate[ctr].mu_f_off, &rate[ctr].mu_b_off);
	}
	fclose(fp);
	return(ctr);
}

int Read_ribo_rates(char *name, ribo_rate_form *rate)
{
	int ctr = 0,bp;
	float fraction;
	char codon[10],line[200];

	FILE *fp;
	fp = fopen(name, "r");
	while (fgets(line,200,fp) != NULL) {
		ctr++;
		sscanf(line,"%d %lf %s %lf %f", &bp,
			&rate[ctr].exp,codon,&rate[ctr].codon,
			&fraction);
		rate[ctr].codon /= 1000.0;
	}
	fclose(fp);
	return(ctr);
}
			
void Read_invariant_rates(double *mu_rnap_p2, double *mu_rnap_p3,
			  double *mu_rnap_e1, double *mu_rnap_e2,
			  double *mu_rnap_e3, double *mu_rnap_N,
			  double *mu_rnap_c, double *mu_rnap_d,
			  int *min_effector, int *rnap_antipause,
			  int *ribo_antipause, int *ribo_rnap_range,
			  int *max_backtrack, char *name)
{
	char check[200],line[200];
	double value;
	FILE *fp;
	fp = fopen(name, "r");
	while (fgets(line,200,fp) != NULL) {
		sscanf(line,"%s %lf",check,&value);
		if (strcmp(check,"kp2") == 0) {
			*mu_rnap_p2 = 1.0/value;
		}
		if (strcmp(check,"kp3") == 0) {
			*mu_rnap_p3 = 1.0/value;
		}
		if (strcmp(check,"ke1") == 0) {
			*mu_rnap_e1 = 1.0/value;
		}
		if (strcmp(check,"ke2") == 0) {
			*mu_rnap_e2 = 1.0/value;
		}
		if (strcmp(check,"ke3") == 0) {
			*mu_rnap_e3 = 1.0/value;
		}
		if (strcmp(check,"kN") == 0) {
			*mu_rnap_N = 1.0/value;
		}
		if (strcmp(check,"kc") == 0) {
			*mu_rnap_c = 1.0/value;
		}
		if (strcmp(check,"kd") == 0) {
			*mu_rnap_d = 1.0/value;
		}
		if (strcmp(check,"min_effector") == 0) {
			*min_effector = (int) (value + 0.01);
		}
		if (strcmp(check,"rnap_antipause") == 0) {
			*rnap_antipause = (int) (value + 0.01);
		}
		if (strcmp(check,"ribo_antipause") == 0) {
			*ribo_antipause = (int) (value + 0.01);
		}
		if (strcmp(check,"ribo_rnap_range") == 0) {
			*ribo_rnap_range = (int) (value + 0.01);
		}
		if (strcmp(check,"max_backtrack") == 0) {
			*max_backtrack = (int) (value + 0.01);
		}
	}
	fclose(fp);
}


void Read_gene_identities(int unit, int *gene_id)
{
	int i,check;
	char line[200];
	FILE *fp;
	fp = fopen("gene_id.all.txt","r");
	while (fgets(line,200,fp) != NULL) {
		sscanf(line,"%d",&check);
		if (check > unit) {
			break;
		}
		if (check == unit) {
			char *str = line;
			str += 5;
			sscanf(str,"%d",&check);
			str += 5;
			for (i = 1; i <= check; i++) {
				sscanf(str,"%d",&gene_id[i]);
				str += 5;
			}
		}
	}
	fclose(fp);
	return;
}

void Read_additional_parameters(char *name)
{
        int pause_loc;
        float pause_dur;
        char label[15],modeltype[15],on_off[20];
        char linearized_check[20];
        char line[200];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
                sprintf(label,"%.12s",line);
                char *str = line;
                str += 30;
                if (strcmp(label,"starting_sig") == 0) {
                        sscanf(str,"%f",&start_sigma);
                }
                if (strcmp(label,"loc_start_ba") == 0) {
                        sscanf(str,"%d",&start_barrier);
                }
                if (strcmp(label,"loc_stop_bar") == 0) {
                        sscanf(str,"%d",&stop_barrier);
                }
                if (strcmp(label,"length_expon") == 0) {
                        sscanf(str,"%f",&resist_exp);
                }
                if (strcmp(label,"min_rot_resi") == 0) {
                        sscanf(str,"%f",&min_resist);
                }
                if (strcmp(label,"max_rot_resi") == 0) {
                        sscanf(str,"%f",&max_resist);
                }
                if (strcmp(label,"implict_forc") == 0) {
                        sscanf(str,"%f",&implicit_f);
                }
                if (strcmp(label,"rate_topoiso") == 0) {
                        sscanf(str,"%f",&rate_topoI);
                }
                if (strcmp(label,"rate_gyrase(") == 0) {
                        sscanf(str,"%f",&rate_gyrase);
                }
                if (strcmp(label,"relative_gyr") == 0) {
                        sscanf(str,"%f",&relative_gyrase_ds_affinity);
                }
                if (strcmp(label,"relative_top") == 0) {
                        sep_topoIA_interrnap_rate = 0;
                        sscanf(str,"%f",&relative_topoIA_ds_binding);
                        if (relative_topoIA_ds_binding < 0.99) {
                                printf("RELATIVE TOPO IA DS BINDING %f\n",
                                        relative_topoIA_ds_binding);
                                sep_topoIA_interrnap_rate = 1;
                                relative_topoIA_ds_binding *=
                                   (((float) (-150 - start_barrier))/
                                    ((float) (stop_barrier + tx_end)));
                                printf("ADJUSTED FOR LENGTH %f\n",
                                        relative_topoIA_ds_binding);
                        }
                }
                if (strcmp(label,"guided_topo_") == 0) {
                        sscanf(str,"%s",on_off);
                        if (on_off[0] == 'y' || on_off[0] == 'Y') {
                                guided_topo_binding = 1;
                                printf("GUIDED TOPOS ON\n");
                        }
                        else {
                                guided_topo_binding = 0;
                                printf("GUIDED TOPOS OFF\n");
                        }
                }
                if (strcmp(label,"guided_gyras") == 0) {
                        sscanf(str,"%s",on_off);
                        if (on_off[0] == 'y' || on_off[0] == 'Y') {
                                guided_gyrase_binding = 1;
                                printf("GUIDED GYRASES ON\n");
                        }
                        else {
                                guided_gyrase_binding = 0;
                                printf("GUIDED GYRASES OFF\n");
                        }
                }
                if (strcmp(label,"inter_rnap_b") == 0) {
                        sscanf(str,"%s",on_off);
                        if (on_off[0] == 'y' || on_off[0] == 'Y') {
                                inter_rnap_topo = 1;
                                printf("INTERRNAP TOPOS OK\n");
                        }
                        else {
                                inter_rnap_topo = 0;
                                printf("NO INTERRNAP TOPOS\n");
                        }
                }
                if (strcmp(label,"lagtime_topo") == 0) {
                        sscanf(str,"%f",&mean_lagtime_topoIA);
                }
                if (strcmp(label,"p_rnap_topoI") == 0) {
                        sscanf(str,"%f",&p_rnap_topo_eject[TOPO_IA]);
                }
                if (strcmp(label,"p_rnap_gyras") == 0) {
                        sscanf(str,"%f",&p_rnap_topo_eject[GYRASE]);
                }
                if (strcmp(label,"max_torque_P") == 0) {
                        sscanf(str,"%f",&P23_release_torque);
                }
                if (strcmp(label,"max_gyrase_s") == 0) {
                        sscanf(str,"%f",&gyrase_min_sigma);
                }
                if (strcmp(label,"topo_IA_widt") == 0) {
                        sscanf(str,"%d",&topoIA_width);
                }
                if (strcmp(label,"gyrase_width") == 0) {
                        sscanf(str,"%d",&gyrase_width);
                }
                if (strcmp(label,"shutoff_equi") == 0) {
                        sscanf(str,"%f",&t_pseudo_shut);
                }
                if (strcmp(label,"tx_on_plasmi") == 0) {
                        sscanf(str,"%s",linearized_check);
                        plasmid =
                           (linearized_check[0] == 'y' ||
                            linearized_check[1] == 'Y');
                }
                if (strcmp(label,"linearized_p") == 0) {                        // Exculsive categories for
                        sscanf(str,"%s",linearized_check);                      // tx on plasmid, linearized
                        plasmid_linearized =                                    // plasmid, and chromosome;
                           (linearized_check[0] == 'y' ||                       // Note that "linearized" is
                            linearized_check[1] == 'Y');                        // ignored if "plasmid" has
                        plasmid_linearized *= plasmid;                          // not been specified
                        plasmid -= plasmid_linearized;                          // explicitly
                        chrom = !(plasmid || plasmid_linearized);
                        if (plasmid) {
                                printf("TX ON PLASMID\n");
                        }
                        else if (plasmid_linearized) {
                                printf("TX ON LINEARIZED PLASMID\n");
                        }
                        else if (chrom) {
                                printf("TX ON CHROMOSOME\n");
                        }
                }
                if (strcmp(label,"5_prime_star") == 0) {
                        sscanf(str,"%d",&fish_five_start);
                }
                if (strcmp(label,"5_prime_stop") == 0) {
                        sscanf(str,"%d",&fish_five_stop);
                }
                if (strcmp(label,"3_prime_star") == 0) {
                        sscanf(str,"%d",&fish_three_start);
                }
                if (strcmp(label,"3_prime_stop") == 0) {
                        sscanf(str,"%d",&fish_three_stop);
                }
                if (strcmp(label,"k_term_stall") == 0) {
                        sscanf(str,"%f",&mu_term);
                        mu_term = 1.0/mu_term;
                }
                if (strcmp(label,"max_RNAP_loa") == 0) {
                         sscanf(str,"%d",&max_RNAP_load_gap);
                }
                if (strcmp(label,"rel_extra_to") == 0) {
                        sscanf(str,"%f",&rel_extra_topo);
                }
                if (strcmp(label,"min_upstream") == 0) {
                        sscanf(str,"%f",&min_upstr_sigma);
                }
                if (strcmp(label,"mfd_surveill") == 0) {
                        sscanf(str,"%f",&mu_rnap_surveillance);
                }
                if (strcmp(label,"mfd_evaluati") == 0) {
                        sscanf(str,"%f",&mu_rnap_eval);
                }
                if (strcmp(label,"sc_insensiti") == 0) {
                        sscanf(str,"%d %f",&pause_loc,&pause_dur);
                        no_super_pause[pause_loc]++;
                        no_super_dur[pause_loc] = pause_dur;
                }
                if (strcmp(label,"RNAP_relax_f") == 0) {
                        sscanf(str,"%f",&relax_factor);
                        if (relax_factor > 0.001) {
                                rnap_relax = ON;
                                relax_factor = 1.0/relax_factor;                // Invert b/c used as mu adj
                        }
                }
                if (strcmp(label,"writhe_parti") == 0) {
                        sscanf(str,"%s",on_off);
                        writhe_partition = (on_off[0] == 'Y' ||
                                            on_off[0] == 'y');
                }
                if (strcmp(label,"no_supercoil") == 0) {
                        sscanf(str,"%s",on_off);
                        supercoiling_off = (on_off[0] == 'Y' ||
                                            on_off[0] == 'y');
                }
                if (strcmp(label,"ribosome_res") == 0) {
                        sscanf(str,"%s",on_off);
                        sprintf(modeltype,"%.13s",on_off);
                        if (strcmp(modeltype,"implicit_RNA_") == 0) {
                                ribo_resist_model = IMPLICIT_EXTANT_RNA;
                        }
                        if (strcmp(modeltype,"implicit_RNAP") == 0) {
                                ribo_resist_model = IMPLICIT_RNAP_POS;
                        }
                        if (strcmp(modeltype,"explicit_ribo") == 0) {
                                ribo_resist_model = EXPLICIT_RIBOSOMES;
                        }
                        if (strcmp(modeltype,"none_RNA_only") == 0) {
                                ribo_resist_model = NO_RIBO_RNA_ONLY;
                        }
                        if (strcmp(modeltype,"flat_resistan") == 0) {
                                ribo_resist_model = FLAT_RESISTANCE;
                        }
                }
                if (strcmp(label,"implicit_int") == 0) {
                        sscanf(str,"%f",&mean_interribo_dist);
                }
                if (strcmp(label,"flat_resist_") == 0) {
                        sscanf(str,"%f",&flat_resist_param);
                }
                if (strcmp(label,"backtracking") == 0) {
                        sscanf(str,"%s",on_off);
                        BACKTRACK_ON = (on_off[0] == 'Y' ||
                                        on_off[0] == 'y');
                }
                if (strcmp(label,"bt_adjustmen") == 0) {
                        sscanf(str,"%f",&bt_adj);
                }
                if (strcmp(label,"no_pause_mod") == 0) {
                        sscanf(str,"%s",on_off);
                        NO_PAUSE_MODEL = (on_off[0] == 'Y' ||
                                          on_off[0] == 'y');
                }
                if (strcmp(label,"anticascade_") == 0) {
                        sscanf(str,"%s",on_off);
                        ANTICASCADE_ON = (on_off[0] == 'Y' ||
                                          on_off[0] == 'y');
                }
                if (strcmp(label,"quiet_mode_o") == 0) {
                        sscanf(str,"%s",on_off);
                        QUIET_MODE = (on_off[0] == 'Y' ||
                                      on_off[0] == 'y');
                }
                if (strcmp(label,"report_dwell") == 0) {
                        sscanf(str,"%s",on_off);
                        REPORT_DWT = (on_off[0] == 'Y' ||
                                      on_off[0] == 'y');
                }
                if (strcmp(label,"DWT_window(n") == 0) {
                        sscanf(str,"%d",&dwt_window);
                }
                if (strcmp(label,"DWT_start_lo") == 0) {
                        sscanf(str,"%d",&dwt_start);
                }
                if (strcmp(label,"DWT_stop_loc") == 0) {
                        sscanf(str,"%d",&dwt_stop);
                }
                if (strcmp(label,"one_RNAP_run") == 0) {
                        sscanf(str,"%s",on_off);
                        ONE_RNAP_RUN = (on_off[0] == 'Y' ||
                                        on_off[0] == 'y');
                }
                if (strcmp(label,"singlemol_tr") == 0) {
                        sscanf(str,"%s",on_off);
                        SINGMOL_TRACE = (on_off[0] == 'Y' ||
                                         on_off[0] == 'y');
                }
                if (strcmp(label,"seq/trace_fr") == 0) {
                        sscanf(str,"%f",&acq_freq);
                }
                if (strcmp(label,"make_kymogra") == 0) {
                        sscanf(str,"%s",on_off);
                        MAKE_KYMOGRAPH = (on_off[0] == 'Y' ||
                                          on_off[0] == 'y');
                }
                if (strcmp(label,"kymo_duratio") == 0) {
                        sscanf(str,"%f",&kymo_duration);
                }
                if (strcmp(label,"kymo_start_t") == 0) {
                        sscanf(str,"%f",&kymo_start);
                }
                if (strcmp(label,"kymo_increme") == 0) {
                        sscanf(str,"%f",&kymo_increment);
                }
                if (strcmp(label,"kymo_min_pos") == 0) {
                        sscanf(str,"%d",&kymo_min_pos);
                }
                if (strcmp(label,"kymo_max_pos") == 0) {
                        sscanf(str,"%d",&kymo_max_pos);
                }
                if (strcmp(label,"kymo_window(") == 0) {
                        sscanf(str,"%d",&kymo_window);
                }
                if (strcmp(label,"bubble_adjus") == 0) {
                        sscanf(str,"%s",on_off);
                        BUBBLE_ADJ = (on_off[0] == 'Y' ||
                                      on_off[0] == 'y');
                }
                if (strcmp(label,"upstream_bub") == 0) {
                        sscanf(str,"%d",&bubble[UPSTREAM]);
                }
                if (strcmp(label,"downstream_b") == 0) {
                        sscanf(str,"%d",&bubble[DOWNSTREAM]);
                }
                if (strcmp(label,"DNA_scrunchi") == 0) {
                        sscanf(str,"%s",on_off);
                        DNA_SCRUNCH = (on_off[0] == 'Y' ||
                                       on_off[0] == 'y');
                }
                if (strcmp(label,"RNA/DNA_hybr") == 0) {
                        sscanf(str,"%d",&hybrid_length);
                }
                if (strcmp(label,"TSS_size_adj") == 0) {
                        sscanf(str,"%d",&TSS_adj);
                }


        }
        fclose(fp);
        return;
}


void Print_assigned_values(float sim_start, float sim_stop,
			   float *seq_window, float *sample_window,
                           float k_loading, float k_on, float k_off,
			   float RNA_lifetime, float prot_lifetime,
			   int num_genes, gene_form *gene_set,
			   int tx_length, char name_set[][200],
			   float snap_time, float fract, int num_traj,
                           float generation_time, float b_time,
			   float c_time, float d_time, float k_unloading,
			   float k_to_open, float k_continue_degrade,
                           float p_rnap_rnap_push,float p_ribo_rnap_push,
			   float p_ribo_ribo_push, float k_CSAT,
			   int pol_size, int ribo_size, int degrade_width,
                           float mean_ribo_rate, float min_dwell,
			   int *coord, float birth_vol, int stable_RNA,
			   int read_state, float mean_RNAP_conc,
			   float mean_ribo_conc, char state_file[][200],
			   int tsl_yn)
{
	int i;
        printf("TESTING:\n");
	printf("COORDINATES OF TX UNIT: %d TO %d\n", coord[0], coord[1]);
	printf("GENERATION TIME: %f\n", generation_time);
	printf("B TIME %f\n",b_time);
	printf("C TIME %f\n",c_time);
	printf("D TIME %f\n",d_time);
	printf("NUMBER OF TRAJECTORIES TO BE RUN: %d\n", num_traj);
        printf("SIM START AND STOP: %f %f\n", sim_start, sim_stop);
        printf("SEQ WINDOW START AND STOP: %f %f\n",
		seq_window[0], seq_window[1]);
        printf("SAMPLE WINDOW START AND STOP: %f %f\n",
		sample_window[0], sample_window[1]);
	printf("TIME OF DNA/mRNA SNAPSHOT: %f\n", snap_time);
        printf("k_on: %f\n", k_on);
        printf("k_off: %f\n", k_off);
	printf("FRACTION OF TIME PROMOTER IS ON: %f\n", fract);
        printf("k_binding: %f\n", k_loading);
        printf("k_unbinding: %f\n", k_unloading);
        printf("k_closed_to_open: %f\n", k_to_open);
	printf("k_cont_degr: %f\n",k_continue_degrade);
	printf("prob rnap/rnap push: %f\n",p_rnap_rnap_push);
	printf("prob ribo/rnap push: %f\n",p_ribo_rnap_push);
	printf("prob ribo/ribo push: %f\n",p_ribo_ribo_push);
	printf("k_CSAT: %f\n",k_CSAT);
        printf("RNA lifetime: %f\n", RNA_lifetime);
        printf("prot lifetime: %f\n", prot_lifetime);
	printf("mean ribo rate: %f\n",mean_ribo_rate);
	printf("min ribo dwell %f\n",min_dwell);
	printf("RNAP width: %d\n",pol_size);
	printf("ribosome width: %d\n", ribo_size);
	printf("degradation_width: %d\n", degrade_width);
	printf("READ STATE FILE?: %d\n", read_state);
	printf("FILES SENT: %s %s\n", state_file[0], state_file[1]);
	printf("STABLE RNA?: %d\n", stable_RNA);
	printf("TRANSLATED?: %d\n", tsl_yn);
	printf("BIRTH VOLUME(fL): %f\n", birth_vol);
	printf("MEAN RNAP CONCENTRATION: %f\n", mean_RNAP_conc);
	printf("MEAN RIBO CONCENTRATION: %f\n", mean_ribo_conc);
        printf("NUM GENES ON TRANSCRIPT: %d\n", num_genes);
        for (i = 1; i <= num_genes; i++) {
                printf("GENE %d: START: %d; STOP %d; k_init %f\n",
			i,gene_set[i].start, gene_set[i].stop,
			gene_set[i].k_init);
        }
        printf("NAME LIST: %s %s %s\n",
		name_set[1], name_set[2], name_set[3]);
        printf("LENGTH OF TX UNIT IS %d\n", tx_length);
	return;
}

void Write_rna_ribosome_snapshot(int num_traj, int max_ribo_per_rna,
				 int **RNA_info, int num_mRNA_snapshot,
				 int traj_index)
{
        int i,j,k;
        char name[200];
        FILE *fp;

        sprintf(name,"rna_ribosome_snapshot_%d.txt", traj_index);
        fp = fopen(name, "w");
        for (i = 1; i <= num_traj; i++) {
               	for (j = 1; j <= num_mRNA_snapshot; j++) {
                       	fprintf(fp,"%7d%7d%7d%7d%7d%7d", i,j,RNA_info[j][1],
			   RNA_info[j][2],RNA_info[j][3],RNA_info[j][4]);
                       	for (k = 1; k <= RNA_info[j][4]; k++) {
                               	fprintf(fp,"%7d", RNA_info[j][k + 4]);
                       	}
                       	fprintf(fp,"\n");
               	}
        }
        fclose(fp);

        return;
}

void Write_rnap_dna_snapshot(int num_promoters, int tx_length,
			     int **rnap_snapshot, int traj_index)
{
        int i,j,rnap_ctr;
        char name[200];
        FILE *fp;

        sprintf(name,"dna_rnap_snapshot_%d.txt", traj_index);
        fp = fopen(name, "w");
        for (i = 1; i <= num_promoters; i++) {
                rnap_ctr = 0;
                for (j = 1; j <= tx_length; j++) {
                        if (rnap_snapshot[i][j] > 0) {
                                rnap_ctr++;
                        }
                }
                fprintf(fp,"%7s%7d%7d"," ",i,rnap_ctr);
                for (j = 1; j <= tx_length; j++) {
                        if (rnap_snapshot[i][j] > 0) {
                                fprintf(fp,"%7d", j);
                        }
                }
                fprintf(fp,"\n");
        }
        fclose(fp);
        return;
}

void Write_pseudoseq_files(int tx_length, int *rnap_seq,
			   int *ribo_seq, int traj_index)
{
        int i;
        char name[200], name2[200];
        FILE *fp;
        FILE *fp2;

        sprintf(name,"pseudo_RNAPseq_%d.txt", traj_index);
        sprintf(name2,"pseudo_riboseq_%d.txt", traj_index);
        fp = fopen(name, "w");
        fp2 = fopen(name2, "w");
        for (i = 1; i <= tx_length; i++) {
                fprintf(fp,"%10d%10d\n",i,rnap_seq[i]);
                fprintf(fp2,"%10d%10d\n",i,ribo_seq[i]);
        }
        fclose(fp);
        fclose(fp2);
        return;
}

void Write_summary_file(char *traj_name, int start_seed, float sim_stop,
			int num_traj, float  k_on, float k_off,
			float k_loading, float *seq_window,
			float *sample_window, float snap_time,
			int tx_length, float fract, float RNA_lifetime,
                        gene_form *gene_set, float  *protein_mean,
			float *protein_sd, float *snapshot_ribo_mean,
                        float snapshot_ribo_sd,
			float *mean_ribos_through_rna,
			float *sd_ribos_through_rna,
                        float fish_final[][5], float fish_final_sd[][5],
			float *fish_pool, float *fish_pool_sd,
			int num_genes, int traj_index, float  time_step,
                        float k_continue_degrade, float p_rnap_rnap_push,
			float p_ribo_rnap_push, float p_ribo_ribo_push,
			int pol_width, int ribo_width, int degrade_width,
                        float mean_ribo_rate, float min_dwell,
			float *mean_elongation_rate, int *coord,
			int stable, int translated, int read_state,
			float generation_time, float b_time, float c_time,
			float d_time, float k_unloading, float k_to_open,
                        float k_CSAT, float birth_vol, float mean_RNAP_conc,
			float mean_ribo_conc, int num_transcripts,
			int *num_proteins_made,
			float *mean_proteins_per_cell,
			float *mean_ribo_gene)
{
        int i;
        char name[200];
	char answer[20];
        FILE *fp;

        sprintf(name,"%s_summary_%d.txt",traj_name, traj_index);
        fp = fopen(name,"w");

        fprintf(fp,"%-40s%20s\n","TRANSCRIPTION UNIT NAME",traj_name);
	fprintf(fp,"%-40s%20d\n","STARTS AT",coord[0]);
	fprintf(fp,"%-40s%20d\n","ENDS AT",coord[1]);
	strcpy(answer,"n");
	if (stable == 1) {
		strcpy(answer,"y");
	}
	fprintf(fp,"%-40s%20s\n","STABLE RNA",answer);
	strcpy(answer,"n");
	if (translated == 1) {
		strcpy(answer,"y");
	}
	if (translated == 2) {
		strcpy(answer,"hybrid");
	}
	fprintf(fp,"%-40s%20s\n","TRANSLATED",answer);
	strcpy(answer,"n");
	if (read_state == 1) {
		strcpy(answer,"y");
	}
	fprintf(fp,"%-40s%20s\n","INITIAL STATE FILE",answer);
        fprintf(fp,"%-40s%20d\n","SIMULATIONS RUN", num_traj);
        fprintf(fp,"%-40s%20d\n", "START SEED", start_seed);
	fprintf(fp,"%-40s%20d\n","GENERATION TIME(MIN)",
		(int) (generation_time/60.0));
	fprintf(fp,"%-40s%20d\n","B TIME(MIN)",(int) (b_time/60.0));
	fprintf(fp,"%-40s%20d\n","C TIME(MIN)",(int) (c_time/60.0));
	fprintf(fp,"%-40s%20d\n","D TIME(MIN)",(int) (d_time/60.0));
        fprintf(fp,"%-40s%20d\n","SIMULATION LENGTH(MIN)",
		((int) (sim_stop + 0.0001))/60 );
        fprintf(fp,"%-40s%20f\n","TIMESTEP(SEC)", time_step);
        fprintf(fp,"%-40s%20d\n","PSEUDOSEQ START(MIN)",
		((int) (seq_window[0] + 0.0001))/60 );
        fprintf(fp,"%-40s%20d\n","PSEUDOSEQ STOP(MIN)",
		((int) (seq_window[1] + 0.0001))/60 );
        fprintf(fp,"%-40s%20d\n","SAMPLE WINDOW START(MIN)",
		((int) (sample_window[0] + 0.0001))/60 );
        fprintf(fp,"%-40s%20d\n","SAMPLE WINDOW STOP(MIN)",
		((int) (sample_window[1] + 0.0001))/60 );
        fprintf(fp,"%-40s%20d\n","RNA SNAPSHOT TIME(MIN)",
		((int) (snap_time + 0.0001))/60 );
        fprintf(fp,"%-40s%20f\n","k_on(PROMOTER)", k_on);
        fprintf(fp,"%-40s%20f\n","k_off(PROMOTER)", k_off);
        fprintf(fp,"%-40s%20f\n","k_on(RNAP)", k_loading);
        fprintf(fp,"%-40s%20f\n","k_off(RNAP)", k_unloading);
        fprintf(fp,"%-40s%20f\n","k_open(RNAP)", k_to_open);
        fprintf(fp,"%-40s%20f\n","FRACTIONAL AVAILABILITY", fract);
        fprintf(fp,"%-40s%20f\n","MEAN TRANSLATION RATE", mean_ribo_rate);
        fprintf(fp,"%-40s%20f\n","MINIMUM DWELLTIME", min_dwell);
        fprintf(fp,"%-40s%20f\n","MEAN RNA LIFETIME(SEC)",RNA_lifetime);
        fprintf(fp,"%-40s%20f\n","RATE DEGR PROGRESS (BASES/SEC)",
		k_continue_degrade);
        fprintf(fp,"%-40s%20f\n","PROB RNAP/RNAP PUSH",p_rnap_rnap_push);
        fprintf(fp,"%-40s%20f\n","PROB RIBOSOME/RNAP PUSH",
		p_ribo_rnap_push);
        fprintf(fp,"%-40s%20f\n","PROB RIBOSOME/RIBOSOME PUSH",
		p_ribo_ribo_push);
	fprintf(fp,"%-40s%20f\n","RATE CSAT",k_CSAT);
	fprintf(fp,"%-40s%20f\n","BIRTH VOLUME(fL)",birth_vol);
	fprintf(fp,"%-40s%20f\n","MEAN [RNAP-free]",mean_RNAP_conc);
	fprintf(fp,"%-40s%20f\n","MEAN [ribo-free]",mean_ribo_conc);
        fprintf(fp,"%-40s%20d\n","RNAP WIDTH (BP)", pol_width);
        fprintf(fp,"%-40s%20d\n","RIBOSOME WIDTH (BASES)", ribo_width);
        fprintf(fp,"%-40s%20d\n","DEGRADE WIDTH (BASES)", degrade_width);
        fprintf(fp,"%-40s%20d\n","TRANSCRIPT LENGTH", tx_length);
        fprintf(fp,"%-40s%20d\n","NUMBER OF GENES ON TRANSCRIPT",
		num_genes);
        for (i = 1; i <= num_genes; i++) {
                fprintf(fp,"START GENE%4d%26s%20d\n",
			i," ",gene_set[i].start);
                fprintf(fp,"STOP GENE%5d%26s%20d\n",
			i," ",gene_set[i].stop);
                fprintf(fp,"k_init GENE%3d%26s%20f\n",
			i," ",gene_set[i].k_init);
        }
        for (i = 1; i <= num_genes; i++) {
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"MEAN PROTEINS PER CELL", "GENE", i,
			mean_proteins_per_cell[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"SD PROTEINS PER CELL", "GENE", i, 0.0);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"MEAN PRODUCTION (PROTEINS/HR)", "GENE",
			i, protein_mean[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"SD PRODUCTION (PROTEINS/HR)", "GENE",
			i, protein_sd[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"MEAN TRANSLATION EFFICIENCY", "GENE",
			i, mean_ribos_through_rna[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"SD TRANSLATION EFFICIENCY", "GENE",
			 i, sd_ribos_through_rna[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"ASSIGNED ELONGATION RATE", "GENE",
			i, (1/mean_ribo_gene[i])/3.0);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"MEAN ELONGATION RATE", "GENE",
			i, mean_elongation_rate[i]);
                fprintf(fp,"%-30s%s%3d   %20f\n",
			"SD ELONGATION RATE", "GENE", i, 0.0);
                fprintf(fp,"%-30s%s%3d   %20d\n",
			"TOTAL PROTEINS PRODUCED", "GENE",
			i, num_proteins_made[i]);
                fprintf(fp,"%-30s%s%3d   %20d\n",
			"SD PROTEINS PRODUCED", "GENE", i, 0);
        	fprintf(fp,"%-30s%s%3d   %20f\n",
			"MEAN RIBOSOMES/RNA (SNAPSHOT)", "GENE",
			i, snapshot_ribo_mean[i]);
        	fprintf(fp,"%-30s%s%3d   %20f\n",
			"SD RIBOSOMES/RNA (SNAPSHOT)", "GENE", i, 0.0);
        }
        fprintf(fp,"%-40s%20f\n",
		"MEAN FISH SIGNAL (5' END)", fish_pool[1]);
        fprintf(fp,"%-40s%20f\n",
		"SD FISH SIGNAL (5' END)", fish_pool_sd[1]);
        fprintf(fp,"%-40s%20f\n",
		"MEAN FISH SIGNAL (MIDGENE)", fish_pool[2]);
        fprintf(fp,"%-40s%20f\n",
		"SD FISH SIGNAL (MIDGENE)", fish_pool_sd[2]);
        fprintf(fp,"%-40s%20f\n",
		"MEAN FISH SIGNAL (3' END)", fish_pool[3]);
        fprintf(fp,"%-40s%20f\n",
		"SD FISH SIGNAL (3' END)", fish_pool_sd[3]);
	fprintf(fp,"%-40s%20d\n",
		"TOTAL TRANSCRIPTS PRODUCED", num_transcripts);
	fprintf(fp,"%-40s%20d\n",
		"SD TRANSCRIPTS PRODUCED", 0);

        fclose(fp);
        return;
}
	
void Write_frame(double time, int promoter_index, FILE *fp)
{
        int i,j,ctr=0,rnap_ctr = 0,rna_ctr=0;
        int ribo_positions[100];
        for (i = 1; i <= tx_end; i++) {
                if (dna_strip[promoter_index][i] > 0) {
                        rnap_ctr++;
                }
        }
        fprintf(fp,"%7d%10.3f%7d%7d",frame,time,0,rnap_ctr);
        if (rnap_ctr > 0) {
                for (i = 1; i <= tx_end; i++) {
                        if (dna_strip[promoter_index][i] > 0) {
                                fprintf(fp,"%7d",rnap[dna_strip[
				  promoter_index][i]][SUPER_INDEX]);
                                fprintf(fp,"%7d", i);
                        }
                }
        }
        fprintf(fp,"\n");
        for (i = 1; i < MAX_RNA; i++) {
                if (rna[i][1] == 0) {
                        continue;
                }
                rna_ctr++;
                fprintf(fp, "%7d", frame);
                fprintf(fp, "%10.3f", time);
                fprintf(fp, "%7d", rna[i][MASTER_INDEX]);
                fprintf(fp, "%7d", rna[i][FIVE_END]);
                fprintf(fp, "%7d", rna[i][THREE_END]);
                fprintf(fp, "%7d", (rna[i][MATURE]  == 0));
                ctr = 0;
                for (j = 1; j <= tx_end; j++) {
                        if (rna_strip[i][j] > 0) {
                                ctr++;
                                ribo_positions[ctr] = j;
                        }
                }
                fprintf(fp, "%7d", ctr);
                for (j = 1; j <= ctr; j++) {
                        fprintf(fp,"%7d", ribosome[rna_strip[i][
					ribo_positions[j]]][SUPER_INDEX]);
                        fprintf(fp,"%7d", ribo_positions[j]);
                }
                fprintf(fp,"\n");
        }
        fprintf(fp,"%7d%10.3f%7d%7d%7d%7d",frame,time, -1,
		num_proteins_made[1],
                promoter[promoter_index][ON_OFF],
		gene_num - 1);
        for (i = 2; i <= gene_num; i++) {
                fprintf(fp,"%7d",num_proteins_made[i]);
        }
        fprintf(fp,"\n");
        return;
}
	
void Prepare_for_no_regulation_simulation_without_file(float *param_set)
{
	param_set[1] = 0.007;
	param_set[2] = 0.0;
	t_pseudo_shut = 90.0;
	return;
}

void Prepare_for_no_supercoiling_simulation_without_file()
{
	supercoiling_off = 1;
	plasmid_linearized = ON;
	start_sigma = 0.0;
	start_barrier = -500;
	stop_barrier = 500;
	resist_exp = 0.0;
	min_resist = 1.0;
	max_resist = 1.0;
	rate_topoI = 0.0;
	rate_gyrase = 0.0;
	P23_release_torque = -10.0;
}
	














