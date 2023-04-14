#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/xdrfile.h"
#include "INCL/xdrfile_xtc.h"
#include "INCL/zoomin_trajectory.h"
#include "INCL/vector_utilities.h"

#define MAX_ATOMS 10000
#define MAX_FRAMES 20000
#define MAX_RNAP 50
#define MAX_DNA_LEN 1000
#define PI 3.141592654

int RNAP_pos[MAX_FRAMES][MAX_RNAP];
int RNAP_rot[MAX_FRAMES][MAX_RNAP];
double sigma[MAX_FRAMES][MAX_DNA_LEN];
int first_rnap[MAX_FRAMES][3];

double width_scaling;

int Read_pdb(char *name, double coord[][3])
{
        int ctr=0;
        double x,y,z;
        char line[150], coord_check[10];

        FILE *fp;
        fp = fopen(name, "r");
        while (fgets(line, 150, fp) != NULL) {
		sprintf(coord_check, "%.4s", line);
		if (strcmp(coord_check,"ATOM") != 0) {
			continue;
		}
                ctr++;
                char *str = line;
                str += 30;
                sprintf(coord_check, "%.8s", str);
                sscanf(coord_check, "%lf", &x);
                str += 8;
                sprintf(coord_check, "%.8s", str);
                sscanf(coord_check, "%lf", &y);
                str += 8;
                sprintf(coord_check, "%.8s", str);
                sscanf(coord_check, "%lf", &z);
                coord[ctr][0] = x;
                coord[ctr][1] = y;
                coord[ctr][2] = z;
        }
        fclose(fp);
        return(ctr);
}

int Read_RNAP_posfile(char *name)
{
        int ctr=0;
	int pos,which_RNAP,rot;
        double sig;
        char line[150];

        FILE *fp;
        fp = fopen(name, "r");
        while (fgets(line, 150, fp) != NULL) {
		char *str = line;
		str += 10;
		sscanf(str,"%d %lf",&pos, &sig);
		if (pos == 1) {
			ctr++;
		}
		sigma[ctr][pos] = sig;
		str += 20;
		sscanf(str,"%d",&which_RNAP);
		RNAP_pos[ctr][which_RNAP] = pos;
		str += 20;
		sscanf(str,"%d",&rot);
		RNAP_rot[ctr][which_RNAP] = rot;
		if (pos == 1) {
			str += 10;
			sscanf(str,"%d %d %d",&first_rnap[ctr][0],
				&first_rnap[ctr][1],&first_rnap[ctr][2]);
		}
	}
	fclose(fp);
	return(ctr);
}

void Write_frame(XDRFILE *xfp_out, int frame, float time, int num_atoms, matrix box, float prec,
                 double coord[][3])
{
	int i,j,shit;
	rvec xrvec[MAX_ATOMS];

	for (i = 1; i <= num_atoms; i++) {
		for (j = 0; j < 3; j++) {
			xrvec[i-1][j] = coord[i][j]/10.0;
		}
	}
	shit = write_xtc(xfp_out, num_atoms, frame, time, box, xrvec, prec);
	shit++;
	return;
}

void Make_xtc_file(int num_atoms, double coord[][3],double rot_angle,
		   double *transloc, int num_frames, int start, int stop,
		   int num_DNA, int num_rnap, char *label)
{
	int i,j;
	int frame=0;
	int frame_ctr=0;
	float time=0.0;
	float prec = 1000.0;
	matrix box;
	char vizname[1000];
	double diameter_adj;
	double test1[3], test2[3];
	double rot_axis[3],rnap_angle;
	double DNA_start1[3],DNA_start2[3],last_DNA1[3],last_DNA2[3],DNA_angle;

	rot_axis[0] = 1.0;
	rot_axis[1] = 0.0;
	rot_axis[2] = 0.0;

	test1[0] = 0.0;
	test1[1] = 15.0;
	test1[2] = 0.0;

        DNA_start1[0] = 0.0;
        DNA_start1[1] = 2.695;
        DNA_start1[2] = -2.238;
        DNA_start2[0] = 1.087;
        DNA_start2[1] = -2.857;
        DNA_start2[2] = 2.041;
        DNA_angle = 360.0/10.5;

	sprintf(vizname,"zoomin_movie.%s.xtc",label);

	XDRFILE *xfp_out = xdrfile_open(vizname, "w");

	Write_frame(xfp_out,frame,time,num_atoms,box,prec,coord);
	frame_ctr++;
	for (frame = 1; frame <= num_frames; frame++) {

		rnap_angle = ((double) first_rnap[frame][2])/105.0;
		rnap_angle *= 360.0;

		Rotate_around_arbitrary_axis(rot_axis,rnap_angle,
			DNA_start1,last_DNA1);
		Rotate_around_arbitrary_axis(rot_axis,rnap_angle,
			DNA_start2,last_DNA2);

		for (i = 1; i <= num_DNA/2; i++) {
			coord[i][0] = DNA_start1[0] + ((double) i);
			Rotate_around_arbitrary_axis(rot_axis,
			   (DNA_angle * (1.0 + sigma[frame][i - 1])),
			    last_DNA1,test2);
			coord[i][1] = test2[1];
			coord[i][2] = test2[2];
			last_DNA1[1] = coord[i][1];
			last_DNA1[2] = coord[i][2];
			diameter_adj = sin((0.5984* (1.0 + sigma[frame][i - 1]))/2);
			diameter_adj = 1.029/diameter_adj;
			diameter_adj /= 3.5;
			coord[i][1] *= diameter_adj;
			coord[i][2] *= diameter_adj;
			coord[i][1] *= width_scaling;
			coord[i][2] *= width_scaling;

			j = num_DNA + 1;
			j -= i;
			coord[j][0] = DNA_start2[0] + ((double) i);
			Rotate_around_arbitrary_axis(rot_axis,
			   (DNA_angle * (1.0 + sigma[frame][i - 1])),
			    last_DNA2,test2);
			coord[j][1] = test2[1];
			coord[j][2] = test2[2];
			last_DNA2[1] = coord[j][1];
			last_DNA2[2] = coord[j][2];
			coord[j][1] *= diameter_adj;
			coord[j][2] *= diameter_adj;
			coord[j][1] *= width_scaling;
			coord[j][2] *= width_scaling;
		}

		for (i = start; i <= stop; i++) {
			coord[i][0] = ((double) RNAP_pos[frame][i - num_DNA]);
			coord[i + num_rnap][0] = coord[i][0];
			rnap_angle = ((double) RNAP_rot[frame][i - num_DNA])/105.0;
			rnap_angle *= 360.0;
			Rotate_around_arbitrary_axis(rot_axis,rnap_angle,test1,test2);
			coord[i][1] = 0.0;
			coord[i][2] = 0.0;
			coord[i + num_rnap][1] = test2[1];
			coord[i + num_rnap][2] = test2[2];
			if (coord[i][0] < -1.0) {
				if (frame < 4) {
					continue;
				}	
				if (RNAP_pos[frame - 1][i - num_DNA] > 0) {				// Move in a rectangle to
					coord[i][0] = ((double) num_DNA) + 100.0;			// allow smoothing of RNAP
					coord[i + num_rnap][0] = ((double) num_DNA) + 100.0;		// movement without finishing
					coord[i][1] = 0.0;						// RNAPs careening backwards
					coord[i + num_rnap][1] = 0.0;					// to tbe start site
				}
				else if (RNAP_pos[frame - 2][i - num_DNA] > 0) {
					coord[i][0] = ((double) num_DNA) + 100.0;
					coord[i + num_rnap][0] = ((double) num_DNA) + 100.0;
					coord[i][1] = 9999.99;
					coord[i + num_rnap][1] = 9999.9;
				}
				else if (RNAP_pos[frame - 3][i - num_DNA] > 0) {
					coord[i][0] = -100.0;
					coord[i + num_rnap][0] = -100.0;
					coord[i][1] = 9999.99;
					coord[i + num_rnap][1] = 9999.9;
				}
				else {
					continue;
				}
			} 
		}
		Write_frame(xfp_out,frame,time,num_atoms,box,prec,coord);
		frame_ctr++;
	}
	xdrfile_close(xfp_out);

	printf("FINSHED WRITING %d FRAMES\n", frame_ctr);

	return;
}

int main (int argc, char *argv[])
{
	int i,j;
	int start,stop;
	int num_rnap;
	int num_atoms;
	int num_frames;
	int num_DNA;
	double coord[MAX_ATOMS][3];
	double transloc[3];
	double rot_angle=0.0;
	char label[200],name[500];

	strcpy(label,argv[1]);

	num_DNA = strtol(argv[3],NULL,10);

	Make_initial_structure_and_bonds(label,num_DNA);

	sprintf(name,"zoomin_movie.%s.pdb",label);

	num_atoms = Read_pdb(name,coord);
	printf("FOUND %d ATOMS IN START PDB...\n",num_atoms);

	printf("...%d OF WHICH ARE DNA ATOMS\n",num_DNA*2);

	num_rnap = (num_atoms - (num_DNA*2))/2;
	printf("ASSUMING %d TWO-BEAD RNAPs FOR ROTATION\n",num_rnap);

	width_scaling = 1.3;
	printf("DNA WIDTH SCALING: %f\n",width_scaling);

	memset(sigma,0,sizeof(sigma));

	for (i = 1; i < MAX_FRAMES; i++) {
		for (j = 1; j < MAX_RNAP; j++) {
			RNAP_pos[i][j] = -100;
			RNAP_rot[i][j] = 0;
		}
	}

	strcpy(name, argv[2]);
	num_frames = Read_RNAP_posfile(name);
	printf("FOuND %d FRAMES IN POSITION FILE\n",num_frames);

	start = num_DNA*2 + 1;
	stop = num_DNA*2 + num_rnap;

	Make_xtc_file(num_atoms,coord,rot_angle,transloc,
		      num_frames,start,stop,num_DNA*2,num_rnap,
		      label);

	Process_rnap_and_supercoiling_changes_for_tcl(name,num_DNA);

	Customize_tcl_script(label,num_DNA);
	
	return(0);
}
	
