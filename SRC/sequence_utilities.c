#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int Read_DNA_sequence (char *dna_name, char *dna_seq)
{
        int c,ctr = 0;
        int info_line=0;

        dna_seq[0] = 'Z';
        FILE *fp;
        fp = fopen(dna_name, "r");
        while ((c = fgetc(fp)) != EOF) {
                if (c == '>') {
                        info_line++;
                }
                if (c == '\n' && info_line) {
                        info_line--;
                }
                if (!(info_line) &&
                    (c == 'A' || c == 'C' || c == 'T' || c == 'G')) {
                        ctr++;
                        dna_seq[ctr] = c;
                }
        }
        dna_seq[ctr + 1] = '\0';
        fclose(fp);

        return ctr;
}

void Generate_sequence_file(int FLAT_SEQ, int RANDOM_SEQ,
                            char *flat_repeat, int rand_seed,
                            int tx_start, int tx_stop)
{
        int i,ctr=0;
        char nt,lett_set[5];
        sprintf(lett_set,"ATCG");
        if (RANDOM_SEQ) {
                if (rand_seed > 0) {
                        srand(rand_seed);
                }
                else {
                        srand(time(NULL));
                }
        }
        FILE *fp;
        fp = fopen("seq.request","w");
        for (i = tx_start - 19; i <= tx_stop + 20; i++) {
                nt = flat_repeat[0];
                if (RANDOM_SEQ) {
                        nt = lett_set[rand() % 4];
                }
                fprintf(fp,"%c",nt);
                ctr++;
                if (ctr % 70 == 0) {
                        fprintf(fp,"\n");
                }
        }
        fclose(fp);
        return;
}

int Print_reverse_complement (char *name)
{
        int i,num_bp,ctr=0;
        char dna_seq[100000];

        printf("READING SEQUENCE: %s\n", name);
        num_bp = Read_DNA_sequence(name, dna_seq);
        printf("FOUND %d BP IN SEQUENCE TO BE REVERSED\n", num_bp);
        FILE *fp;
        fp = fopen("rev_compl.seq","w");
        for (i = num_bp; i >= 1; i--) {
                if (dna_seq[i] == 'A') {
                        fprintf(fp,"T");
                }
                else if (dna_seq[i] == 'T') {
                        fprintf(fp,"A");
                }
                else if (dna_seq[i] == 'C') {
                        fprintf(fp,"G");
                }
                else {
                        fprintf(fp,"C");
                }
                ctr++;
                if (ctr % 70 == 0 && i != 1) {
                        fprintf(fp,"\n");
                }
        }
        fclose(fp);
        return(num_bp);
}

void Make_RNA_sequence(char *dna_seq, char *rna_sequence, int start,
		       int stop, int rev)
{
        int i,ctr=0;

        rna_sequence[0] = 'Z';

        if (rev == 0) {
                for (i = start; i <= stop; i++) {
                        ctr++;
                        if (dna_seq[i] == 'T') {
                                rna_sequence[ctr] = 'U';
                        }
                        else {
                                rna_sequence[ctr] = dna_seq[i];
                        }
                }
        }
        else {
                for (i = stop; i >= start; i--) {
                        ctr++;
                        if (dna_seq[i] == 'T') {
                                rna_sequence[ctr] = 'A';
                        }
                        else if (dna_seq[i] == 'A') {
                                rna_sequence[ctr] = 'U';
                        }
                        else if (dna_seq[i] == 'G') {
                                rna_sequence[ctr] = 'C';
                        }
                        else {
                                rna_sequence[ctr] = 'G';
                        }
                }
        }

        rna_sequence[ctr + 1] = '\0';

        return;
}


void Trim_sequence(char *seq_name, int start_pad, int stop_pad)
{
        int i,num_bp,written=0;
        char dna_seq[100000];

        num_bp = Read_DNA_sequence(seq_name,dna_seq);
        printf("FOUND %d BASE-PAIRS IN SEQUENCE FILE\n",
                num_bp);

        FILE *fp_out;
        fp_out = fopen("trimmed.seq","w");
        for (i = 1; i <= num_bp; i++) {
                if (i > start_pad && i <= num_bp - stop_pad) {
                        fprintf(fp_out,"%c",dna_seq[i]);
                        written++;
                }
                if (written && written % 70 == 0) {
                        fprintf(fp_out,"\n");
                }
        }
        fclose(fp_out);

        return;
}


int Wraparound(int check)
{
        int adj,GENOME_SIZE=4639675;
        adj=check;
        if (adj < 1) {
                adj += GENOME_SIZE;
        }
        if (adj > GENOME_SIZE) {
                adj -= GENOME_SIZE;
        }
        return(adj);
}

