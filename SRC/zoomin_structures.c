#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/vector_utilities.h"

#define MAX_ATOMS 1000000
#define PI 3.141592654


int Assign_coordinates_to_DNA_template(int *tx_unit, double coord[][3],
				       int bonds[][2])
{
	int i,j,ctr=0;
	double test[3],test2[3],post_rot[3];
	double angle, rot_axis[3];

	test[0] = 0.0;
	test[1] = 2.695;
	test[2] = -2.238;
	test2[0] = 1.087;
	test2[1] = -2.857;
	test2[2] = 2.041;
	rot_axis[0] = 1.0;
	rot_axis[1] = 0.0;
	rot_axis[2] = 0.0;
	angle = 360.0/10.5;

	for (i = tx_unit[0]; i <= tx_unit[1]; i++) {
		ctr++;
		coord[ctr][0] = test[0] + ((double) ctr);
		Rotate_around_arbitrary_axis(
		   rot_axis,((double) (ctr - 1)) * angle,test,post_rot);
		coord[ctr][1] = post_rot[1];
		coord[ctr][2] = post_rot[2];

		j = (tx_unit[1] * 2) + 1;
		j -= ctr;
		coord[j][0] = test2[0] + ((double) ctr);
		Rotate_around_arbitrary_axis(
		   rot_axis,((double) (ctr - 1)) * angle,test2,post_rot);
		coord[j][1] = post_rot[1];
		coord[j][2] = post_rot[2];
	}
	for (i = 2; i <= tx_unit[1]; i++) {
		bonds[i - 1][0] = i - 1;
		bonds[i - 1][1] = i;
	}
	for (i = tx_unit[1] + 2; i <= tx_unit[1] * 2; i++) {
		bonds[i - 2][0] = i - 1;
		bonds[i - 2][1] = i;
	}
	return(ctr);
}

void Write_initial_pdb(char *label, int num_dna, int num_rnap,int num_RNA,
		       int num_ribo, double coord[][3])
{
	int i;
	char name[200];
	
	sprintf(name,"zoomin_movie.%s.pdb",label);
	FILE *fp;
	fp = fopen(name,"w");
	for (i = 1; i <= num_dna * 2; i++) {
		fprintf(fp, "ATOM%7d  P   DNA A%4d%12.3f%8.3f%8.3f  "
			"1.00  0.00\n",
			i,1,coord[i][0],coord[i][1],coord[i][2]);
	}
	for (i = num_dna * 2 + 1; i <= num_dna * 2 + num_rnap + num_rnap; i++) {
		fprintf(fp, "ATOM%7d  C   DNA B%4d%12.3f%8.3f%8.3f  "
			"1.00  0.00\n",
			i,2,coord[i][0],coord[i][1],coord[i][2]);
	}
	fclose(fp);
	return;
}

void Write_bonds(char *label, int num_bonds, int num_atoms,int bonds[][2])
{
	int i;
	char name[200];
	
	sprintf(name,"zoomin_movie.%s.psf",label);
	FILE *fp;
	fp = fopen(name,"w");
        fprintf(fp,"PSF whole-genome\n\n       1 !NTITLE\n REMARKS "
		   "psf for tx/transl snapshot\n\n");
        fprintf(fp,"%8d !NATOM\n", num_atoms);
        for (i = 1; i <= num_atoms; i++) {
//                fprintf(fp,"%8d%2s%5d%3s%4s%5s%11.5f%14.4f%12d\n",
//                i,"U",i,"DNA","N","N",0.0,1000.0,0);
                fprintf(fp,"%8d%2s%5d    DNA  N    NH3   -0.300000       "
			   "14.0070           0\n",i,"U",1);
        }
        fprintf(fp,"\n%8d !NBOND\n", num_bonds);
        for (i = 1; i <= num_bonds; i++) {
                fprintf(fp,"%8d%8d", bonds[i][0], bonds[i][1]);
                if (i % 4 == 0 || i == num_bonds) {
                        fprintf(fp,"\n");
                }
        }
        fclose(fp);
        return;
}

void Make_initial_structure_and_bonds(char *label, int DNA_length)
{
	int i;
	int num_bonds=0,bonds[1000000][2];
	int num_dna=0,num_rnap=0,num_RNA=0,num_ribo=0;
	int tx_unit[2], num_atoms;
	double coord[MAX_ATOMS][3];

	tx_unit[0] = 1;
	tx_unit[1] = DNA_length; 

	num_dna = Assign_coordinates_to_DNA_template(tx_unit,coord,bonds);

	num_rnap = 15;

	for (i = 1; i <= num_rnap; i++) {
		coord[num_dna * 2 + i][0] = -100.0;
		coord[num_dna * 2 + num_rnap + i][0] = -100.0;
		coord[num_dna * 2 + i][1] = 0.0;
		coord[num_dna * 2 + num_rnap + i][1] = 15.0;
		coord[num_dna * 2 + i][2] = 0.0;
		coord[num_dna * 2 + num_rnap + i][2] = 0.0;
	}
	num_bonds = num_dna - 1;
	num_bonds *= 2;
	num_atoms = num_dna + num_dna + num_rnap + num_rnap;
	num_RNA = 0;
	num_ribo = 0;
	Write_initial_pdb(label,num_dna,num_rnap,num_RNA,num_ribo,coord);
	Write_bonds(label,num_bonds,num_atoms,bonds);

	return;
}
	
