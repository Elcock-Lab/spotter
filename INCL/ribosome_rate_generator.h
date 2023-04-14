#ifndef RIBOFILE_GFN_H
#define RIBOFILE_GEN_H

int Read_gene_file(char *name);

int Read_ribo_file(char *name);

void Assign_decode_only(char *gene_name, char *rate_name);

int Read_gene_info_file(char *name, int displace);

int Read_codon_info_file(char *name);

int Assign_default_codon_information();

void Add_codon_info(char *rna_sequence, double *reads, char *name,
                    char gene_label[][2][10], int off_1, int off_2,
                    int num_bp, int FLAT_DWELL, int FLAT_DECODE);

void Add_codon_based_decoding_info(char *seq_name, char *gene_name,
                                   char *dwell_name, char *codon_name,
                                   int start, int stop,
                                   double flat_decode_time);

#endif
