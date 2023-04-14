#ifndef READ_UTIL_H
#define READ_UTIL_H

int Read_wig_file(char *name, double *reads);

void Write_normalized_reads(char *name, double *norm, double adj,
                            int num_bp);

void Normalize_rates(char *dwellfile, double nt_per_sec);

#endif
