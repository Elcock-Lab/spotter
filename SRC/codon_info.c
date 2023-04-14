#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

#define MAX_BP 5000000

typedef struct {
	double tot_time;
	double frac_decode;
} codon_form;

codon_form codon[5][5][5];

int f_gene_ct[5000000];
char f_gene[5000000][2][10];


int Read_gene_info_file(char *name, int displace)
{
        int i,ctr=0,start,stop;
        int *gene_ct;
        char line[200],gene[50];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
                ctr++;
                sscanf(line,"%s %d %d",gene,&start,&stop);
		start += displace;
		stop  += displace;
                gene_ct = f_gene_ct;
                for (i = start; i <= stop; i += 3) {
                	gene_ct[i]++;
                	strcpy(f_gene[i][gene_ct[i] - 1],gene);
                }
        }
        fclose(fp);
        return(ctr);
}

int Read_codon_info_file(char *name)
{
	int i,ctr=0,finder[4];
	double tot,frac;
	char line[200],check[10];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
		ctr++;
		sscanf(line,"%s %lf %lf",check,&frac,&tot);
		for (i = 0; i <= 2; i++) {
			if (check[i] == 'A') {
				finder[i] = 1;
			}				
			else if (check[i] == 'U') {
				finder[i] = 2;
			}				
			else if (check[i] == 'C') {
				finder[i] = 3;
			}
			else {
				finder[i] = 4;
			}
		}
		codon[finder[0]][finder[1]][finder[2]].frac_decode = frac;
		codon[finder[0]][finder[1]][finder[2]].tot_time = tot;
	}
	fclose(fp);
	return(ctr);
}

int Assign_default_codon_information()
{
	codon[1][1][1].frac_decode = 0.901693492;
	codon[1][1][3].frac_decode = 0.949563039;
	codon[1][1][4].frac_decode = 0.937144967;
	codon[1][1][2].frac_decode = 0.915120149;
	
	codon[1][3][1].frac_decode = 0.941444327;
	codon[1][3][3].frac_decode = 0.956103821;
	codon[1][3][4].frac_decode = 0.937031931;
	codon[1][3][2].frac_decode = 0.829943887;
	
	codon[1][4][1].frac_decode = 0.968241032;
	codon[1][4][3].frac_decode = 0.956967995;
	codon[1][4][4].frac_decode = 0.986524289;
	codon[1][4][2].frac_decode = 0.919508875;
	
	codon[1][2][1].frac_decode = 0.936399477;
	codon[1][2][3].frac_decode = 0.934929839;
	codon[1][2][4].frac_decode = 0.970750293;
	codon[1][2][2].frac_decode = 0.913586325;
	
	codon[3][1][1].frac_decode = 0.955669842;
	codon[3][1][3].frac_decode = 0.959039279;
	codon[3][1][4].frac_decode = 0.95611657;
	codon[3][1][2].frac_decode = 0.971006834;
	 
	codon[3][3][1].frac_decode = 0.969129615;
	codon[3][3][3].frac_decode = 0.904088758;
	codon[3][3][4].frac_decode = 0.962246218;
	codon[3][3][2].frac_decode = 0.951685084;
	
	codon[3][4][1].frac_decode = 0.77728441;
	codon[3][4][3].frac_decode = 0.792123915;
	codon[3][4][4].frac_decode = 0.988995904;
	codon[3][4][2].frac_decode = 0.722175658;
	
	codon[3][2][1].frac_decode = 0.976862589;
	codon[3][2][3].frac_decode = 0.958330729;
	codon[3][2][4].frac_decode = 0.80484463;
	codon[3][2][2].frac_decode = 0.96879954;
	
	codon[4][1][1].frac_decode = 0.873686847;
	codon[4][1][3].frac_decode = 0.936160455;
	codon[4][1][4].frac_decode = 0.790206998;
	codon[4][1][2].frac_decode = 0.898230527;
	
	codon[4][3][1].frac_decode = 0.908696804;
	codon[4][3][3].frac_decode = 0.983283307;
	codon[4][3][4].frac_decode = 0.800810114;
	codon[4][3][2].frac_decode = 0.771657003;
	
	codon[4][4][1].frac_decode = 0.981968927;
	codon[4][4][3].frac_decode = 0.854427579;
	codon[4][4][4].frac_decode = 0.914298041;
	codon[4][4][2].frac_decode = 0.778131451;
	  
	codon[4][2][1].frac_decode = 0.91513395;
	codon[4][2][3].frac_decode = 0.968437013;
	codon[4][2][4].frac_decode = 0.819108889;
	codon[4][2][2].frac_decode = 0.708442451;
	
	codon[2][1][1].frac_decode = 0.762633554;
	codon[2][1][3].frac_decode = 0.90965522;
	codon[2][1][4].frac_decode = 0.95202488;
	codon[2][1][2].frac_decode = 0.851834672;
	 
	codon[2][3][1].frac_decode = 0.91486655;
	codon[2][3][3].frac_decode = 0.967740588;
	codon[2][3][4].frac_decode = 0.931677419;
	codon[2][3][2].frac_decode = 0.837853018;
	
	codon[2][4][1].frac_decode = 0.789170586;
	codon[2][4][3].frac_decode = 0.938858703;
	codon[2][4][4].frac_decode = 0.948018758;
	codon[2][4][2].frac_decode = 0.902246091;
	
	codon[2][2][1].frac_decode = 0.966270822;
	codon[2][2][3].frac_decode = 0.962161261;
	codon[2][2][4].frac_decode = 0.865787087;
	codon[2][2][2].frac_decode = 0.940641681;
	
	codon[1][1][1].tot_time    = 34.58570631;
	codon[1][1][3].tot_time    = 67.41088116;
	codon[1][1][4].tot_time    = 54.09272449;
	codon[1][1][2].tot_time    = 40.05662097;
	
	codon[1][3][1].tot_time    = 58.06440012;
	codon[1][3][3].tot_time    = 77.45548938;
	codon[1][3][4].tot_time    = 53.99562128;
	codon[1][3][2].tot_time    = 19.99340068;
	
	codon[1][4][1].tot_time    = 107.0563745;
	codon[1][4][3].tot_time    = 79.01095949;
	codon[1][4][4].tot_time    = 252.3057879;
	codon[1][4][2].tot_time    = 42.24068169;
	
	codon[1][2][1].tot_time    = 53.4586795;
	codon[1][2][3].tot_time    = 52.25129252;
	codon[1][2][4].tot_time    = 116.2404806;
	codon[1][2][2].tot_time    = 39.34562418;
	
	codon[3][1][1].tot_time    = 76.69722324;
	codon[3][1][3].tot_time    = 83.00635171;
	codon[3][1][4].tot_time    = 77.47799162;
	codon[3][1][2].tot_time    = 117.2690154;
	
	codon[3][3][1].tot_time    = 110.1379206;
	codon[3][3][3].tot_time    = 35.44944198;
	codon[3][3][4].tot_time    = 90.05720332;
	codon[3][3][2].tot_time    = 70.37164273;
	
	codon[3][4][1].tot_time    = 15.26610688;
	codon[3][4][3].tot_time    = 16.35589779;
	codon[3][4][4].tot_time    = 308.9758559;
	codon[3][4][2].tot_time    = 12.23794855;
	
	codon[3][2][1].tot_time    = 146.9481609;
	codon[3][2][3].tot_time    = 81.59490032;
	codon[3][2][4].tot_time    = 17.42201614;
	codon[3][2][2].tot_time    = 108.972753;
	
	codon[4][1][1].tot_time    = 26.91722838;
	codon[4][1][3].tot_time    = 53.25852475;
	codon[4][1][4].tot_time    = 16.20645093;
	codon[4][1][2].tot_time    = 33.40883962;
	
	codon[4][3][1].tot_time    = 37.23856515;
	codon[4][3][3].tot_time    = 203.3895051;
	codon[4][3][4].tot_time    = 17.06913976;
	codon[4][3][2].tot_time    = 14.88988079;
	 
	codon[4][4][1].tot_time    = 188.5633746;
	codon[4][4][3].tot_time    = 23.3560723;
	codon[4][4][4].tot_time    = 39.67237184;
	codon[4][4][2].tot_time    = 15.32438921;
	
	codon[4][2][1].tot_time    = 40.063135;
	codon[4][2][3].tot_time    = 107.72111;
	codon[4][2][4].tot_time    = 18.79583786;
	codon[4][2][2].tot_time    = 11.66150564;
	
	codon[2][1][1].tot_time    = 14.32384422;
	codon[2][1][3].tot_time    = 37.63360765;
	codon[2][1][4].tot_time    = 70.87006836;
	codon[2][1][2].tot_time    = 22.94733901;
	
	codon[2][3][1].tot_time    = 39.93729843;
	codon[2][3][3].tot_time    = 105.3955962;
	codon[2][3][4].tot_time    = 49.76392801;
	codon[2][3][2].tot_time    = 20.9686296;
	
	codon[2][4][1].tot_time    = 16.12678198;
	codon[2][4][3].tot_time    = 55.60889571;
	codon[2][4][4].tot_time    = 65.40821067;
	codon[2][4][2].tot_time    = 34.78121775;
	 
	codon[2][2][1].tot_time    = 100.8029317;
	codon[2][2][3].tot_time    = 89.85500296;
	codon[2][2][4].tot_time    = 25.33288291;
	codon[2][2][2].tot_time    = 57.27925005;

	return(64);
}
	
void Add_codon_info(char *rna_sequence, double *reads, char *name,
		    char gene_label[][2][10], int off_1, int off_2,
		    int num_bp, int FLAT_DWELL, int FLAT_DECODE,
		    double flat_decode_time)
{
	int i,j;
	int finder[4];
	char label[5];

	FILE *fp;
	fp = fopen(name,"w");

	for (i = 1 + off_1; i <= num_bp - off_2; i++) {
		for (j = i; j <= i + 2; j++) {
			if (rna_sequence[j] == 'A') {
				finder[(j-i)+1] = 1;
			}
			else if (rna_sequence[j] == 'U') {
				finder[(j-i)+1] = 2;
			}
			else if (rna_sequence[j] == 'C') {
				finder[(j-i)+1] = 3;
			}
			else {
				finder[(j-i)+1] = 4;
			}
			label[j-i] = rna_sequence[j];
		}
		label[3] = '\0';
		if (FLAT_DWELL) {
			reads[i] = 0.011111;
		}
		if (FLAT_DECODE) {
			codon[finder[1]][finder[2]][finder[3]].tot_time =
				flat_decode_time;
			codon[finder[1]][finder[2]][finder[3]].frac_decode =
			  (codon[finder[1]][finder[2]][finder[3]].tot_time-3.4)/
			  codon[finder[1]][finder[2]][finder[3]].tot_time;
		}
		fprintf(fp,"%10d%20.10f%10s%10.5f%10.5f%10s%10s\n",
			i - off_1, reads[i], label,
			codon[finder[1]][finder[2]][finder[3]].tot_time,
			codon[finder[1]][finder[2]][finder[3]].frac_decode,
			gene_label[i][0], gene_label[i][1]);
	}
	fclose(fp);
	return;
}

void Add_codon_based_decoding_info(char *seq_name, char *gene_name,
				   char *dwell_name, char *codon_name,
				   int start, int stop,
				   double flat_decode_time)
{
	int num_bp;
	int num_reads,num_codon,num_genes;
	double f_reads[5000000];
	char dna_seq[MAX_BP];
	char rna_sequence[100000];

	int FLAT_DWELL = 0;
	int FLAT_DECODE = 0;

	memset(f_reads,0,sizeof(f_reads));

	num_bp = Read_DNA_sequence(seq_name, dna_seq);
	printf("FOUND %d BASE-PAIRS IN SEQUENCE FILE\n",
		num_bp);

	num_genes = Read_gene_info_file(gene_name,start);
	printf("FOUND %d GENE(S) IN GENE FILE\n",
		num_genes);

	if (strcmp(dwell_name,"make_flat_dwell") != 0) {
		num_reads = Read_wig_file(dwell_name,f_reads);
		printf("FOUND %d LINES IN DWELL TIME FILE\n",num_reads);
	}
	else {
		FLAT_DWELL++;
	}

	if (strcmp(codon_name,"make_flat_decode") != 0) {
		if (strcmp("fluitt_based_transloc_decoding_split.txt",
		    codon_name) != 0) {
			num_codon = Read_codon_info_file(codon_name);
		}
		else {
			num_codon = Assign_default_codon_information();
		}
		printf("FOUND INFORMATION FOR %d CODONS\n",num_codon);
	}
	else {
		FLAT_DECODE++;
	}

	Make_RNA_sequence(dna_seq, rna_sequence,1,num_bp,0);
	Add_codon_info(rna_sequence, f_reads, "ribo_rates.seq.full_info",
		       f_gene,start,stop, num_bp,FLAT_DWELL, FLAT_DECODE,
		       flat_decode_time);
	
	return;
}




