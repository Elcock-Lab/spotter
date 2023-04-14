#ifndef ZOOMIN_TRAJ_H
#define ZOOMIN_TRAH_H

int Read_RNAP_infofile(char *name,int DNA_length);

void Customize_tcl_script(char *label, int DNA_length);

void Process_rnap_and_supercoiling_changes_for_tcl(char *name, int DNA_length);

int Assign_coordinates_to_DNA_template(int *tx_unit, double coord[][3],
                                       int bonds[][2]);

void Write_initial_pdb(char *label, int num_dna, int num_rnap,int num_RNA,
                       int num_ribo, double coord[][3]);

void Write_bonds(char *label, int num_bonds, int num_atoms,int bonds[][2]);

void Make_initial_structure_and_bonds(char *label, int DNA_length);

#endif
