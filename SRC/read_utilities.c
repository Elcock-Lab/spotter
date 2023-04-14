#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int Read_wig_file(char *name, double *reads)
{
        int ctr=0,loc;
	double count;
        char line[200];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
                ctr++;
		sscanf(line,"%d %lf",&loc,&count);
		reads[ctr] = count;
	}
	fclose(fp);
	return(ctr);
}

void Write_normalized_reads(char *name, double *norm, double adj,
			    int num_bp)
{
	int i;
	FILE *fp;
	fp = fopen(name,"w");
	for (i = 1; i <= num_bp; i++) {
		fprintf(fp,"%10d%20.10f\n",i,norm[i]*adj);
	}
	fclose(fp);
	return;
}

void Normalize_rates(char *dwellfile, double nt_per_sec)
{
	int i;
	int num_bp;
	double f_dwell[5000000];
	double adj,check;

	num_bp = Read_wig_file(dwellfile,f_dwell);
	printf("FOUND %d ENTRIES IN DWELL TIME FILE\n",num_bp);

	adj = 0.0;
	for (i = 1; i <= num_bp; i++) {
		adj += f_dwell[i];
	}
	adj = ((double) (num_bp))/adj;
	printf("TARGET GENE TRANSCRIBED AT %f NT PER SEC UNADJ\n",adj);
	adj /= nt_per_sec;

	check = 0.0;
	for (i = 1; i <= num_bp; i++) {
		check += (f_dwell[i] * adj);
	}
	check = ((double) (num_bp))/check;
	printf("TARGET GENE NOW TRANSCRIBED AT %f NT PER SEC\n",check);

	Write_normalized_reads("normalized_dwelltime.seq",
				f_dwell,adj,num_bp);

	return;
}



















