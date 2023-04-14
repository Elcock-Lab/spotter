#ifndef SEQ_UTIL_H
#define SEQ_UTIL_H

int Read_DNA_sequence (char *dna_name, char *dna_seq);

void Generate_sequence_file(int FLAT_SEQ, int RANDOM_SEQ,
                            char *flat_repeat, int rand_seed,
                            int tx_start, int tx_stop);

int Print_reverse_complement (char *name);

void Make_RNA_sequence(char *dna_seq, char *rna_sequence, int start,
                       int stop, int rev);

void Trim_sequence(char *seq_name, int start_pad, int stop_pad);

int Wraparound(int check);

#endif
