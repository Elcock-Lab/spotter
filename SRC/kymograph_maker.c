#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ATOMS 10000
#define MAX_FRAMES 20000
#define MAX_RNAP 50
#define MAX_DNA_LEN 1000
#define PI 3.141592654

int RNAP_pos[MAX_FRAMES][MAX_RNAP];
int RNAP_pos_coarse[MAX_FRAMES][MAX_RNAP];
int RNAP_pos_2X_coarse[MAX_FRAMES][MAX_RNAP];
int RNAP_rot[MAX_FRAMES][MAX_RNAP];
double sigma[MAX_FRAMES][MAX_DNA_LEN];
double sigma_coarse[MAX_FRAMES][MAX_DNA_LEN];
double min_sigma,max_sigma,sigma_extr;
int first_rnap[MAX_FRAMES][3];
double final_grid[1000][1000];
int RNAP_counted[MAX_DNA_LEN];

int num_rnap=0;
int num_times;
double time,time_cutoff,time_list[1000];
int start_list[50][2];
int pos_list[50][1000];

typedef struct {
        int orig_id;
        int start_time;
        int pos_list[1000];
} RNAP_form;
RNAP_form RNAP_list[50];

int Compare(const void * a, const void * b)
{

  RNAP_form *RNAP_A = (RNAP_form *)a;
  RNAP_form *RNAP_B = (RNAP_form *)b;

  return ( RNAP_B->start_time - RNAP_A->start_time );
}


int Read_RNAP_posfile(char *name)
{
        int ctr=0;
	int pos,which_RNAP,rot;
        double sig;
        char line[150];

        FILE *fp;
        fp = fopen(name, "r");
        while (fgets(line, 150, fp) != NULL) {
		char *str = line;
		str += 10;
		sscanf(str,"%d %lf",&pos, &sig);
		if (pos == 1) {
			ctr++;
		}
		sigma[ctr][pos] = sig;
		str += 20;
		sscanf(str,"%d",&which_RNAP);
		RNAP_pos[ctr][which_RNAP] = pos;
		str += 20;
		sscanf(str,"%d",&rot);
		RNAP_rot[ctr][which_RNAP] = rot;
		if (pos == 1) {
			str += 10;
			sscanf(str,"%d %d %d",&first_rnap[ctr][0],
				&first_rnap[ctr][1],&first_rnap[ctr][2]);
		}
	}
	fclose(fp);
	return(ctr);
}

void Make_custom_R_script_for_viewing(char *label, int num_dna, double time)
{
	int i,syst_ret;
	char script_name[500],colorset[30][200];
	char dir_name[500],command[5000];

	sprintf(colorset[1],"dodgerblue4");
	sprintf(colorset[2],"purple");
	sprintf(colorset[3],"magenta");
	sprintf(colorset[4],"hotpink");
	sprintf(colorset[5],"sienna");
	sprintf(colorset[6],"gray48");
	sprintf(colorset[7],"darkgreen");
	sprintf(colorset[8],"maroon");
	sprintf(colorset[9],"indianred4");
	sprintf(colorset[10],"ivory");
	sprintf(colorset[11],"plum");
	sprintf(colorset[12],"darkred");
	sprintf(colorset[13],"darkolivegreen1");
	sprintf(colorset[14],"darkorchid");
	sprintf(colorset[0],"black");

	sprintf(script_name,"custom_script_for_%s.R",label);
	FILE *fp;
	fp = fopen(script_name,"w");
	fprintf(fp,"options(\"install.lock\"=FALSE)\n");
	fprintf(fp,"install.packages(\"ggplot2\")\n");
	fprintf(fp,"library(ggplot2)\n");
	fprintf(fp,"install.packages(\"reshape\")\n");
	fprintf(fp,"library(reshape)\n");
	fprintf(fp,"install.packages(\"this.path\")\n");
	fprintf(fp,"library(this.path)\n");
	fprintf(fp,"setwd(this.dir())\n");
	fprintf(fp,"getwd()\n\n");
	fprintf(fp,"# Set directory to location of kymograph file if"
		   " this fails to set it properly -"
		   " use setwd(\"FOLDERNAME\")\n;"
		   "# then run remaining commands\n\n");
	fprintf(fp,"kymo_sigma_%s <- "
		   "read.table(\"sigma_info_kymograph.%s.txt\", "
		   "header = TRUE)\n",
		    label,label);
	fprintf(fp,"kymo_rnap_%s <- "
		   "read.fwf(\"rnap_info_for_kymograph.%s.txt\", "
		   "widths = c(8,",
		   label,label);
	for (i = 1; i <= num_rnap; i++) {
		fprintf(fp,"8");
		if (i != num_rnap) {
			fprintf(fp,",");
		}
		else {
			fprintf(fp,"))\n");
		}
	}
	fprintf(fp,"colnames(kymo_rnap_%s) <- c(\"JUNK\",",label);
	for (i = num_rnap; i >= 1; i--) {
		fprintf(fp,"\"RNAP_%d\"",i);
		if (i != 1) {
			fprintf(fp,",");
		}
		else {
			fprintf(fp,")\n");
		}
	}
	fprintf(fp,"kymo_rnap_%s.melty <- "
		   "melt(kymo_rnap_%s,id=\"JUNK\")\n",
		   label,label);
	fprintf(fp,"ggplot(data = kymo_sigma_%s) + "
		   "geom_tile(aes(x=time,y=position,fill=sigma)) + "
		   "scale_fill_gradient2(low = \"#FF0000\",mid = \"#FFFFFF\", "
		   "high = \"#1E90FF\",limits=c(-0.05,0.05)) + "
		   "geom_line(data = kymo_rnap_%s.melty, "
		   "aes(x=JUNK,y=value,color=variable),lwd=1.1) + "
		   "scale_x_continuous(limits = c(0,%d), expand = c(0,0)) + "
		   "scale_y_continuous(limits = c(0,%d), expand = c(0,0)) + "
		   "scale_color_manual(values = c(",
		   label,label,(int) time, num_dna);
	for (i = 1; i <= num_rnap; i++) {
		fprintf(fp,"\"RNAP_%d\" = \"%s\" ",i,colorset[i % 15]);
		if (i != num_rnap) {
			fprintf(fp,",");
		}
		else {
			fprintf(fp,")) + coord_flip() + "
				   "guides(color = guide_legend(reverse=TRUE))\n");
		}
	}
	fclose(fp);

	sprintf(dir_name,"KYMOGRAPH_INFO.%s",label);

	sprintf(command,"mkdir %s",dir_name);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"mv sigma_info_kymograph.%s.txt %s; "
			"mv rnap_info_for_kymograph.%s.txt %s; "
			"mv custom_script_for_%s.R %s", 
			label,dir_name,label,dir_name,label,dir_name);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	return;
}
	

void Build_final_RNAP_list()
{
	int i,j;
	int num_tx_starts,start_loc[10];

	qsort(RNAP_list,50,sizeof(RNAP_form),Compare);

        for (i = 0; i < 50; i++) {
                if (RNAP_list[i].start_time > -10000) {
                        num_rnap++;
                }
        }

        for (i = 1; i <= num_rnap; i++) {
                num_tx_starts = 0;
                memset(start_loc,0,sizeof(start_loc));
                for (j = 1; j <= num_times; j++) {
                        if (RNAP_list[i].pos_list[j] > 0 &&
                            (!(RNAP_list[i].pos_list[j - 1] > 0) ||
                             j == 1)) {
                                num_tx_starts++;
                                start_loc[num_tx_starts] = j;
                        }
                }
                if (num_tx_starts > 1) {
                        num_rnap++;
                        for (j = 1; j <= num_times; j++) {
                                if (j < start_loc[2]) {
                                        RNAP_list[num_rnap].pos_list[j] =0;
                                }
                                else {
                                        RNAP_list[num_rnap].pos_list[j] =
                                         RNAP_list[i].pos_list[j];
                                        RNAP_list[i].pos_list[j] = 0;
                                }
                        }
                        RNAP_list[num_rnap].start_time = start_loc[2];
                }
        }

        qsort(RNAP_list,50,sizeof(RNAP_form),Compare);

        printf("ORDER OF RNAPS: %d %d %d %d %d %d %d %d %d %d\n",
                RNAP_list[0].orig_id, RNAP_list[1].orig_id, 
		RNAP_list[2].orig_id, RNAP_list[3].orig_id, 
		RNAP_list[4].orig_id, RNAP_list[5].orig_id,
		RNAP_list[6].orig_id, RNAP_list[7].orig_id,
		RNAP_list[8].orig_id, RNAP_list[9].orig_id);

        printf("WILL CREATE TABLE FOR PLOTTING %d RNAPs\n",num_rnap);


	return;
}	


int main (int argc, char *argv[])
{
	int i,j;
	int num_dna,num_frames;
	int bp_per_point,frames_per_point;
	double time_per_frame;
	int frame_check;
	char label[200],name[200];
	char sigma_name[500],rnap_name[500];

	strcpy(name,argv[1]);
	num_dna = strtol(argv[2],NULL,10);
	bp_per_point = strtol(argv[3],NULL,10);
	frames_per_point = strtol(argv[4],NULL,10);
	time_per_frame = strtod(argv[5],NULL);
	time_cutoff = strtod(argv[6],NULL);
	strcpy(label,argv[7]);

	memset(sigma,0,sizeof(sigma));
	memset(sigma_coarse,0,sizeof(sigma_coarse));
	memset(final_grid,0,sizeof(final_grid));

	for (i = 1; i < MAX_FRAMES; i++) {
		for (j = 1; j < MAX_RNAP; j++) {
			RNAP_pos[i][j] = -100;
			RNAP_pos_coarse[i][j] = -100;
			RNAP_pos_2X_coarse[i][j] = -100;
			RNAP_rot[i][j] = 0;
		}
	}

	num_frames = Read_RNAP_posfile(name);
	printf("FOuND %d FRAMES IN POSITION FILE\n",num_frames);
	printf("NUM DNA POS %d; BP PER POINT %d; FRAMES PER POINT %d; TIME PER FRAME %f\n",
		num_dna, bp_per_point, frames_per_point, time_per_frame);

	for (i = 1; i <= num_frames; i++) {
		for (j = 1; j <= num_dna; j++) {							// Find mean in bin for sigma
			sigma_coarse[i][((j-1) / bp_per_point) + 1] += sigma[i][j];
		}
		for (j = 1; j <= num_dna / bp_per_point; j++) {
			sigma_coarse[i][j] /= ((double) bp_per_point);
		}
		for (j = 1; j < MAX_RNAP; j++) {							// Just coarsen RNAP position
			RNAP_pos_coarse[i][j] = ((RNAP_pos[i][j]-1) /  bp_per_point) + 1;
		}
	}

	memset(RNAP_counted,0,sizeof(RNAP_counted));
	for (i = 1; i <= num_frames; i++) {
		for (j = 1; j <= num_dna / bp_per_point; j++) {
			final_grid[((i-1) / frames_per_point) + 1][j] += sigma_coarse[i][j];
		}
		for (j = 1; j < MAX_RNAP; j++) {
			if (RNAP_pos_coarse[i][j] > 0) {
				RNAP_pos_2X_coarse[((i-1) / frames_per_point) + 1][j] += 
					RNAP_pos_coarse[i][j];
				RNAP_counted[j]++;
			}
		}
		if (i % frames_per_point == 0) {
			for (j = 1; j <= num_dna / bp_per_point; j++) {
				final_grid[((i-1) / frames_per_point) + 1][j] /= 
					((double) frames_per_point);
			}
			for (j = 1; j < MAX_RNAP; j++) {
				if (RNAP_counted[j]) {
				   RNAP_pos_2X_coarse[((i-1) / frames_per_point) + 1][j] /=
					RNAP_counted[j];
				}
			}
			memset(RNAP_counted,0,sizeof(RNAP_counted));
		}
	}

        for (i = 0; i < 50; i++) {
                RNAP_list[i].start_time = -99900;
                RNAP_list[i].orig_id = i;
        }

	for (i = 1; i <= num_frames / frames_per_point; i++) {
		time = ((double) i) * time_per_frame *
			((double) frames_per_point);
		if (time > time_cutoff) {
			num_times = i - 1;
			break;
		}
		time_list[i] = time;
		frame_check = ((i - 1) * frames_per_point) +
			      (frames_per_point / 2);
		for (j = 1; j < 25; j++) {					
			if (RNAP_pos[frame_check][j] > 0) {
				RNAP_list[j].pos_list[i] =
				   (((RNAP_pos[frame_check][j] - 1) /
				      bp_per_point) + 1) * bp_per_point;
			}
			else {
				RNAP_list[j].pos_list[i] = 0;
			}
			if (RNAP_list[j].start_time < -1000 &&
			    RNAP_list[j].pos_list[i] > 0) {
				if (i == 1) {
				   RNAP_list[j].start_time =
					RNAP_list[j].pos_list[i] / 5;
				   RNAP_list[j].start_time *= -1;
				}
				else {
				   RNAP_list[j].start_time = i;
				}
			}
		}
		num_times = i;
	}

	min_sigma = 9999.9;
	max_sigma = -9999.9;
	FILE *fp;
	sprintf(sigma_name,"sigma_info_kymograph.%s.txt",label);
	fp = fopen(sigma_name,"w");
	fprintf(fp," time  position  sigma\n");			
	for (i = 1; i <= num_frames / frames_per_point; i++) {
		for (j = 1; j <=  num_dna / bp_per_point; j++) {
			fprintf(fp,"%8.3f%8d%8.3f\n",
				((double) i) * time_per_frame *
				((double) frames_per_point),
				j * bp_per_point, final_grid[i][j] );
			if (final_grid[i][j] < min_sigma) {
				min_sigma = final_grid[i][j];
			}
			if (final_grid[i][j] > max_sigma) {
				max_sigma = final_grid[i][j];
			}
		}
	}
	fclose(fp);
	
	Build_final_RNAP_list();

	sigma_extr = (fabs(min_sigma) > fabs(max_sigma) ? 
		      fabs(min_sigma) : fabs(max_sigma));
	sigma_extr = (sigma_extr >= 0.05 ? sigma_extr : 0.05);

	printf("MIN COARSENED SIGMA IS %f; MAX COARSENED SIGMA IS: %f\n"
		"SCALE WILL RANGE FROM %f TO %f\n",
		min_sigma,max_sigma,sigma_extr,-1.0*sigma_extr);	
	
        FILE *fp_rnap;
	sprintf(rnap_name,"rnap_info_for_kymograph.%s.txt",label);
        fp_rnap = fopen(rnap_name,"w");
        for (i = 1; i <= num_times; i++) {
                fprintf(fp_rnap,"%8.3f",time_list[i]);
                for (j = 0; j < num_rnap; j++) {
                        if (RNAP_list[j].pos_list[i] > 0) {
                                fprintf(fp_rnap,"%8d",
					RNAP_list[j].pos_list[i]);
                        }
                        else {
                                if (RNAP_list[j].pos_list[i+1] > 0) {
                                        fprintf(fp_rnap,"%8d",0);
                                }
                                else if (RNAP_list[j].pos_list[i-1] > 0 &&
                                         RNAP_list[j].pos_list[i-1] < 500) {
                                        fprintf(fp_rnap,"%8d",500);
                                }
                                else {
                                        fprintf(fp_rnap,"%8s"," ");
                                }
                        }
                }
                fprintf(fp_rnap,"\n");
        }
        fclose(fp_rnap);
		
	Make_custom_R_script_for_viewing(label,num_dna,time_cutoff);
	
	return(0);
}
	
