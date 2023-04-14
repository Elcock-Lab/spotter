#ifndef RNAP_FILEGEN_H
#define RNAP_FILEGEM_H


void Align_sequence_with_TSS(char *orig_seq, char *final_seq,
                             int num_bp, int start_pos, int stop_pos,
                             int *aligned_start, int *aligned_stop);

double Calculate_energies(int site, char *dna_seq, int dir,
                          double dna_init_GC, double dna_init_AT,
                          double rna_dna_init, double *transloc_energy,
                          double *rna_hybrid, int scrunch,
                          int sliding_window);

void Write_energies(int aligned_start, int aligned_stop,
                    double *energy, double *pre, double *post,
                    int dir);

void Assign_hybridization_energies(char *seq_name, int start_pos,
                                   int stop_pos, char *scruncher,
                                   char *slider);

void Read_NET_seq_energy_function(double energies[][4]);

void Assign_relative_dwelltimes(char *dna_seq, int *seq_range,
                                double energies[][4], int aligned_start,
                                int aligned_stop);

void Assign_NETseq_dwelltimes(char *seq_name);

void Read_parameters(char *name);

void Assign_default_parameters();

double Inter_or_extrapolate(int num_trials, double target);

double Find_dwell_time(double kp);

int Find_offpathway_rate(double target);

void Calculate_final_rates(char *param_name, char *rate_name,
                           char *dwell_name);

#endif

