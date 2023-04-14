#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/zoomin_tagalong_tcl.h"

#define MAX_ATOMS 10000
#define MAX_FRAMES 20000
#define MAX_RNAP 20
#define PI 3.141592654

int last_RNAP_state[25];


int Read_RNAP_infofile(char *name,int DNA_length)
{
        int ctr=0,rnap_ctr=0,dom_ctr=0;
	int dom_start=0,dom_end;
	int pos,which_RNAP,state;
	int super_write_ct=0,rnap_write_ct=0;
        double sc,last_sc=999999.9,sc_color;
        char line[150];

	FILE *fp_sc_out;
	FILE *fp_rnap_out;
	fp_sc_out = fopen("supercoiling_changes.tmp","w");
	fp_rnap_out = fopen("rnap_state_changes.tmp","w");

	fprintf(fp_sc_out,"proc load_super_traj_array {} {\n");
	fprintf(fp_sc_out," global lines\n");
	fprintf(fp_sc_out," global numlines\n");

	fprintf(fp_rnap_out,"proc load_rnap_traj_array {} {\n");
	fprintf(fp_rnap_out," global lines2\n");
	fprintf(fp_rnap_out," global numlines2\n");


        FILE *fp;
        fp = fopen(name, "r");
        while (fgets(line, 150, fp) != NULL) {
		char *str = line;
		str += 10;
		sscanf(str,"%d %lf",&pos, &sc);
		if (pos == 1) {
			ctr++;
			rnap_ctr = 0;
			dom_ctr = 0;
			last_sc = 999999.9;
		}
		if (sc != last_sc) {
			sc_color = last_sc/0.1;
			if (sc_color < -1.0) {
				sc_color = -1.0;
			}
			if (sc_color > 1.0) {
				sc_color = 1.0;
			}
			if (pos != 1) {
				dom_end = pos - 1;
				super_write_ct++;
				fprintf(fp_sc_out,
				  " set lines(%d) \"%10d%10d%10d%10d%10.3f\"\n",
				  super_write_ct,ctr,dom_ctr,dom_start,
				  dom_end,sc_color);
			}
			last_sc = sc;
			dom_ctr++;
			dom_start = pos;
		}
		if (pos == DNA_length) {
			sc_color = sc/0.1;
			if (sc_color < -1.0) {
				sc_color = -1.0;
			}
			if (sc_color > 1.0) {
				sc_color = 1.0;
			}
			super_write_ct++;
			fprintf(fp_sc_out,
				" set lines(%d) \"%10d%10d%10d%10d%10.3f\"\n",
				super_write_ct,ctr,dom_ctr,dom_start,DNA_length,sc_color);
		}
				
		str += 20;
		sscanf(str,"%d %d",&which_RNAP,&state);
		if (state == 0) {
			state = 1;
		}
		if (which_RNAP && last_RNAP_state[which_RNAP] != state) {
			rnap_ctr++;
			rnap_write_ct++;
			fprintf(fp_rnap_out,
				" set lines2(%d) \"%10d%10d%10d%10d\"\n",
				rnap_write_ct,ctr, rnap_ctr, which_RNAP, state);
			last_RNAP_state[which_RNAP] = state;
		}
	}

	fprintf(fp_sc_out," set numlines %d\n}\n",super_write_ct);
	fprintf(fp_rnap_out," set numlines2 %d\n}\n",rnap_write_ct);

	fclose(fp);
	fclose(fp_sc_out);
	fclose(fp_rnap_out);
	return(ctr);
}

void Customize_tcl_script(char *label, int DNA_length)
{
	int syst_ret;
	char tcl_template_string[1000000];
	char command[1000],name[500];

	Assign_tcl_template_string(tcl_template_string);
	FILE *fp;
	fp = fopen("temp.temp.tcl","w");
	fprintf(fp,"%s\n",tcl_template_string);
	fclose(fp);

	sprintf(name,"zoomin_viewer.%s.tcl",label);

	sprintf(command,"sed -i \'s/TRAJIC/%s/g\' temp.temp.tcl",
		label);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }
	sprintf(command,"sed -i \'s/1000 /%d /g\' temp.temp.tcl",
		DNA_length * 2);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }
	sprintf(command,"sed -i \'s/1001/%d /g\' temp.temp.tcl",
		DNA_length * 2 + 1);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command, "cat temp.temp.tcl supercoiling_changes.tmp "
			 "rnap_state_changes.tmp > %s",name);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command, "rm temp.temp.tcl; "
			 "rm supercoiling_changes.tmp; "
			 "rm rnap_state_changes.tmp");
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command, "mkdir ZOOMIN_MOVIE_MATERIALS.%s", label);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"cp zoomin_movie.%s* ZOOMIN_MOVIE_MATERIALS.%s;"
			"cp %s ZOOMIN_MOVIE_MATERIALS.%s",
			label,label,name,label);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	return;
}

void Process_rnap_and_supercoiling_changes_for_tcl(char *name,int DNA_length)
{
	int i,num_frames;

	for (i = 1; i < 25; i++) {
		last_RNAP_state[i] = -1;
	}
	num_frames = Read_RNAP_infofile(name,DNA_length);
	printf("FOUND %d FRAMES IN POSITION FILE\n",num_frames);

	
	return;
}
	
