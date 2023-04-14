#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int total_ribo;
int ribo_list[5000];

int total_rna;
int rna_list[100];


typedef struct {
	int rna;
	int num_obs;
	int pos[10000];
	double time[10000];
} ribo_form;
ribo_form ribosome[5000];

typedef struct {
	int num_obs;
	int length[100000];
	double time[100000];
} rna_form;
rna_form rna[100];


int Read_tx_tsl_trajectory(char *name)
{
	int i,j;
	int frame,identity,num,num2,num3;
	int max_rna = 0;
	int match,ribo,pos;
	double time;
	char line[1000];

	FILE *fp;
	fp = fopen(name,"r");
	while (fgets(line,1000,fp) != NULL) {
        	sscanf(line,"%d %lf %d %d %d",
			&frame,&time,&identity,&num,&num2);
		if (identity > 0) {
			if (identity > max_rna) {
				max_rna = identity;
			}
			rna[identity].num_obs++;
			rna[identity].time[rna[identity].num_obs] = time;
			rna[identity].length[rna[identity].num_obs] =
				(num2 - num) + 1;

			char *str = line;
			str += 45;
			sscanf(str,"%d",&num3);
			if (num3 > 0) {
				str -= 7;
				for (i = 1; i <= num3; i++) {
				  str += 14;
				  sscanf(str,"%d %d",&ribo,&pos);
				  match = 0;
				  for (j = 1; j <= total_ribo; j++) {
					if (ribo == ribo_list[j]) {
						match = j;
						break;
					}
				  }
				  if (!(match)) {
					total_ribo++;
					ribo_list[total_ribo] = ribo;
					match = total_ribo;
					ribosome[match].rna = identity;
				  }
				  ribosome[match].num_obs++;
				  ribosome[match].time[ribosome[match].num_obs] =
					time;
				  ribosome[match].pos[ribosome[match].num_obs] =
					pos;
				}	
			}
		}
	}
	fclose(fp);


	return(max_rna);
}
			
void Write_R_script(char *r_str, int num_rna)
{
        sprintf(r_str,
                "\n"
                "options(\"install.lock\"=FALSE)\n"
                "install.packages(\"ggplot2\")\n"
                "library(ggplot2)\n"
                "install.packages(\"reshape\")\n"
                "library(reshape)\n"
                "install.packages(\"plotly\")\n"
                "library(plotly)\n"
                "install.packages(\"this.path\")\n"
                "library(this.path)\n"
                "setwd(this.dir())\n"
                "\n"
                "# Change to directory where the ribotrace is located if "
                "script fails here with setwd(\"DIRECTORY NAME\")\n"
                "# and run remaining commands\n"
                "tsl <- read.table(\"rna_and_ribotrace.for_R.txt\","
                "header=TRUE)\n"
                "all_ribos <- tsl[grepl(\"RIBO\",tsl$id),]\n"
                "all_rna <- tsl[!grepl(\"RIBO\",tsl$id),]\n"
                "rna <- list()\n"
                "ribos <- list()\n"
                "for (i in 1:%d) {\n"
                "  rna[[i]] <- subset(all_rna,RNA==i)\n"
                "  ribos[[i]] <- subset(all_ribos, RNA==i)\n"
                "}\n"
                "\n"
                "# The section below makes plots for each of the "
                "individual RNAs;\n"
                "# Display by typing \"rna_plot[[N]]\" where N "
                "is the number of the RNA you want to display \n"
                "rna_plot <- list()\n"
                "for (i in 1:%d) {\n"
                "  rna_plot[[i]] <- plot_ly(rna[[i]], x= ~time, y = ~pos,"
                "type = 'scatter',mode='lines',name=~id)\n"
                "  rna_plot[[i]] <- rna_plot[[i]] %%>%% "
                "add_trace(data=ribos[[i]],x = ~time, y = ~pos, split = ~is,"
                " type = 'scatter', mode ='lines',color=~id)\n"
                "}\n"
                "\n"
                "# The section below will make a 3D plot showing "
                "all of the RNAs;\n"
                "# Display by typing \"multiRNA_3d_plot\"\n"
                "multiRNA_3d_plot <- plot_ly(tsl, x = ~time, y = ~RNA, z = ~pos, "
                "split = ~id, type = 'scatter3d', mode = 'lines', "
                "line = list(width = 2), color = ~id)\n"
                "multiRNA_3d_plot <- multiRNA_3d_plot %%>%% "
                "layout(showlegend=FALSE)\n"
                "multiRNA_3d_plot <- multiRNA_3d_plot %%>%% "
                "add_trace(data=all_rna, x = ~time, "
                "y = ~RNA, z = ~pos, split = ~id, type = 'scatter3d', "
                "mode = 'lines', "
                "line = list(width = 5), color = ~RNA)\n"
                "multiRNA_3d_plot <- multiRNA_3d_plot %%>%% hide_colorbar()\n"
                "\n",num_rna,num_rna);

                return;
}


int main (int argc, char *argv[])
{
	int i,j;
	int num_rna;
	int max_rna;
	char r_str[100000],dir_name[1000];
	char traj_name[1000],label[1000];
	char dir_label[1000],command[1000];

	strcpy(traj_name,argv[1]);

	max_rna = strtol(argv[2],NULL,10);

	strcpy(dir_label,argv[3]);

	total_ribo = 0;
	memset(ribo_list,0,sizeof(ribo_list));
	memset(ribosome,0,sizeof(ribosome));
	memset(rna,0,sizeof(rna));

	num_rna = Read_tx_tsl_trajectory(traj_name);
	printf("FOUND %d RNA(S) IN TRAJECTORY\n",num_rna);
	printf("FOUND %d TOTAL RIBOSOMES IN TRAJECTORY\n",total_ribo);
	printf("WILL PRINT DATA FOR %d RNA(S) AND ASSOCIATED RIBOSOMES\n",
		(max_rna < num_rna ? max_rna : num_rna));

	if (num_rna < max_rna) {
		max_rna = num_rna;
	}

	FILE *fp;
	fp = fopen("rna_and_ribotrace.for_R.txt","w");
	fprintf(fp,"%10s%10s%10s%10s\n","time","pos","RNA","id");
	for (i = 1; i <= max_rna; i++) {
		sprintf(label,"RNA_%d",i);
		for (j = 1; j <= rna[i].num_obs; j++) {
			fprintf(fp,"%10.5f%10d%10d%10s\n",
				rna[i].time[j], rna[i].length[j],
				i, label);
		}
	}
	for (i = 1; i <= total_ribo; i++) {
		if (ribosome[i].rna > max_rna) {
			continue;
		}
		sprintf(label,"RIBO_%d",i);
		for (j = 1; j <= ribosome[i].num_obs; j++) {
			fprintf(fp,"%10.5f%10d%10d%10s\n",
				ribosome[i].time[j], ribosome[i].pos[j],
				ribosome[i].rna, label);
		}
	}
				
	fclose(fp);

	Write_R_script(r_str,max_rna);

	FILE *fp_r;
	fp_r = fopen("rna_and_ribosome_plotter.R","w");
	fprintf(fp_r,"%s",r_str);
	fclose(fp_r);

	sprintf(dir_name,"RNA_TSL_PLOTS.%s",dir_label);
	sprintf(command,"mkdir %s; "
			"mv rna_and_ribotrace.for_R.txt %s; "
			"mv rna_and_ribosome_plotter.R %s",
			dir_name, dir_name, dir_name);
	system(command);
	

	return(0);
}

