#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

char orig_line[20000][200];
char final_line[20000][200];
double val[20000];
int gene[25][2];

int Read_gene_file(char *name)
{
        int ctr=0,start,stop;
        char line[200],gene_id[50];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
                ctr++;
                sscanf(line,"%s %d %d",gene_id,&start,&stop);
		gene[ctr][0] = start;
		gene[ctr][1] = stop;
        }
        fclose(fp);
        return(ctr);
}

int Read_ribo_file(char *name)
{
	int ctr = 0;
        char line[200];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
                ctr++;
		char *str = line;
		str += 40;
		sscanf(str,"%lf",&val[ctr]);
		val[ctr] /= 1000.0;
		strcpy(orig_line[ctr],line);
		strcpy(final_line[ctr],line);
	}
	fclose(fp);
	return(ctr);
}

void Assign_decode_only(char *gene_name, char *rate_name)
{
	int i,j,k;
	int num_reads,num_genes;

	num_genes = Read_gene_file(gene_name);
	printf("FOUND %d GENE(S) IN GENE FILE\n",
		num_genes);

	num_reads = Read_ribo_file(rate_name);
	printf("FOUND %d LINES IN RIBOSOME RATE FILE\n",num_reads);

	for (i = 1; i <= num_genes; i++) {
		for (j = gene[i][0]; j <= gene[i][1]; j += 3) {
			for (k = 0; k <= 2; k++) {
				char *str = orig_line[j + k];
				str += 30;
				sprintf(final_line[j +k],
				  "%10d%20.10f%s", j + k,
				  val[j]/3.0,str);
			}
		}
	}
	
	FILE *fp;
	fp = fopen("decoding_only.ribo_rates","w"); 
	for (i = 1; i <= num_reads; i++) {
		fprintf(fp,"%s",final_line[i]);
	}
	fclose(fp);

	return;
}




