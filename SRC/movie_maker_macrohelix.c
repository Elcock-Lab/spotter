#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/xdrfile.h"
#include "INCL/xdrfile_xtc.h"
#include "INCL/vector_utilities.h"
#include "INCL/tagalong_tcl_template.h"

#define MAX_ATOMS 50000
#define MAX_FRAMES 10000
#define PI 3.141592654

#pragma GCC diagnostic ignored "-Wstrict-overflow"

int num_bonds=0,bonds[10000][2],num_dna=0,num_rnap=0,num_RNA=0,num_ribo=0;
int tx_unit[2], num_snapshots, num_atoms;
int ribo_info[3];
int num_frames,tx_end;
int rna_key[MAX_ATOMS];
int num_rna,max_ribo, max_rnap;
int frame_rank[10000][MAX_ATOMS];
double coord[MAX_ATOMS][3];
double sim_time[MAX_FRAMES];

int num_proteins;
int proteins[MAX_FRAMES];
int promoter_state[MAX_FRAMES];
int protein_offset,promoter_offset,ribostart_offset,ribostop_offset;

int num_genes;
int tsl_start[30],tsl_stop[30];
char gene_label[30][5];

double start_bp_radius,radius_incr_per_ring,bp_drop_per_step;
double rnap_pos[10000][3];
double rnap_rot[10000];

int num_RNAP_check;
int RNAP_BT_checklist[5000];

typedef struct {
	double time;
	int num_rna;
	int num_nascent;
	int num_mature;
	int num_ribo;
	int num_proteins[30];
	int num_new_ribo;
	int new_ribo_index1[50];
	int new_ribo_index2[50];
	int new_ribo_gene[50];
} frame_form;
frame_form frame_info[MAX_FRAMES];

char name_tag[100];

int Compare(const void *a, const void *b)
{
        return(*(int*)b - *(int*)a);
}

int Get_starting_trajectory_information(char *name, int *num_rna,
					int *max_rnap, int *max_ribo,
                                        double *sim_time, int *rna_key,
					int *tx_end,
					int frame_rank[][MAX_ATOMS])
{
	int i;
	int identity,num,sum=0,previous=0;
	int last=0;
	int frame,max_frame=0;
	int match,num_list=0;
	int rna_list[1000];
	int temp_rank[1000] = {0},rank_ctr=0;
	int num_addl_genes=0;
	double time;
	char line[1000];

	FILE *fp_init;
	fp_init = fopen(name,"r");
	while (fgets(line,1000,fp_init) != NULL) {
		sscanf(line,"%d %lf %d %d",&frame,&time,&identity,&num);
		if (identity > 0) {
		sscanf(line,"%d %lf %d %d %d",
		  &frame,&time,&identity,&num,&last);
			if (last > *tx_end) {
				*tx_end = last;
			}
		}
		if (identity == -1) {
			sscanf(line,"%d %lf %d %d %d",&frame,&time,
				&identity,&num_proteins,&num);
			proteins[frame] = num_proteins;
			promoter_state[frame] = num;
			frame_info[frame].num_proteins[1] = num_proteins;
			frame_info[frame].time = time;
			char *str = line;
			if (strlen(line) > 40) {
				str += 38;
				sscanf(str,"%d",&num_addl_genes);
			}
			if (num_addl_genes) {
				for (i = 1; i <= num_addl_genes; i++) {
				  str += 7;
				  sscanf(str,"%d",
					&frame_info[frame].num_proteins[i+1]);
				}
			}
		}	
	}
	fclose (fp_init);
	printf("DONE WITH FIRST PASS\n");
	
	FILE *fp;
	fp = fopen(name, "r");
	temp_rank[0] = 1000000000;
	while (fgets(line,1000,fp) != NULL) {
		char *str = line;
		sscanf(line,"%d %lf %d %d",&frame,&time,&identity,&num);
		if (frame > max_frame) {
			max_frame = frame;
		}
		if (frame != previous && frame > 0) {
			qsort(temp_rank, rank_ctr + 1, sizeof(int), Compare);
			for (i = 1; i <= rank_ctr; i++) {
				frame_rank[frame][temp_rank[i]] = i;
			}
			rank_ctr = 0;
			memset(temp_rank,0,sizeof(temp_rank));
			temp_rank[0] = 1000000000;
		}	
		sim_time[frame] = time;
		if (identity == 0 && num > *max_rnap) {
			*max_rnap = num;
		}
		if (identity > 0) {
			frame_info[frame].num_rna++;
			match = 0;
			for (i = 1; i <= num_list; i++) {
				if (identity == rna_list[i]) {
					match++;
					break;
				}
			}
			if (match == 0) {
				num_list++;
				rna_list[num_list] = identity;
				rna_key[identity] = num_list;
			}
			str += 31;
			sscanf(str,"%d",&num);
			if (num == *tx_end) {
				rank_ctr++;
				temp_rank[rank_ctr] = identity;
				frame_info[frame].num_mature++;
			}
			str += 14;
			sscanf(str,"%d",&num);
			frame_info[frame].num_ribo += num;
			if (frame == previous) {
				sum += num;
			}
			else {
				sum = num;
			}
			if (sum > *max_ribo) {
				*max_ribo = sum;
			}
			previous = frame;
		}
	}
	fclose(fp);
	qsort(temp_rank, rank_ctr + 1, sizeof(int), Compare);
	for (i = 1; i <= rank_ctr; i++) {
		frame_rank[frame][temp_rank[i]] = i;
	}
	
	*num_rna = num_list;
	return(max_frame);
}		

void Read_input_for_gene_info(char *name)
{
	int i,which,where,found_tsl_region=0;
	char line[1000], junk[500],check[50];

	num_genes = 0;
	FILE *fp;
	fp = fopen(name, "r");
	while (fgets(line,1000,fp) != NULL) {
		sprintf(check,"%.7s",line);
		if (strcmp(check,"num_RBS")== 0) {
			found_tsl_region++;
		}
		if (!(found_tsl_region)) {
			continue;
		}
		if (strcmp(check,"RBS_loc") == 0) {
			num_genes++;
			sscanf(line, "%s %d %d", junk, &which, &where);
			tsl_start[which] = where;
		}
		if (strcmp(check,"tsl_sto") == 0) {
			sscanf(line, "%s %d %d", junk, &which, &where);
			tsl_stop[which] = where;
		}
	}
	fclose(fp);
	for (i = 0; i < 10; i++) {
		sprintf(gene_label[i+1],"C%d",i);
		sprintf(gene_label[i+11],"O%d",i);
	}
	sprintf(gene_label[21],"NE");
	sprintf(gene_label[22],"NZ");
	sprintf(gene_label[23],"NH");
	return;
}

void Build_helical_track()
{
    int i,j;
    double k;
    double rot_axis[3];
    double radius_incr_per_degree;
    double curr_radial_displace,curr_rotation,spiral_rotation=0.0;
    double dist,test1[3],test2[3],test3[3],prev_dist,rot_dist_adj;

    memset(rnap_rot,0,sizeof(rnap_rot));
    if (start_bp_radius <= 10000) {
        radius_incr_per_degree = radius_incr_per_ring/360.0;
        curr_radial_displace = start_bp_radius;
        coord[1][2] = -1.0 * curr_radial_displace;
        coord[tx_unit[1]*2][2] = -1.0 * curr_radial_displace;
        test1[0] = 0.0;
        test1[1] = 0.0;
        rot_axis[0] = 0.0;
        rot_axis[1] = 1.0;
        rot_axis[2] = 0.0;
        curr_rotation = 0.0;
        rnap_rot[1] = curr_rotation;
        for (i = 2; i <= tx_unit[1]; i++) {
                test1[2] = -1.0 * curr_radial_displace;
                test2[0] = 0.0;
                test2[1] = -1.0 * bp_drop_per_step;
                prev_dist = 0.0;
                for (k = 0.01; k <= 10.0; k += 0.01) {
                        rot_dist_adj = k * radius_incr_per_degree;
                        test2[2] = test1[2] - rot_dist_adj;
                        Rotate_around_arbitrary_axis(
			  rot_axis,k,test2,test3);
                        dist = Find_interatom_distance(test3,test1);
                        if (dist >1.0 && prev_dist < 1.0) {
                                curr_rotation += k;
                                curr_radial_displace += rot_dist_adj;
                                break;
                        }
                        prev_dist = dist;
                }
                rnap_rot[i] = curr_rotation;
                test3[0] = coord[i][0];
                test3[0] -= ((double) (i - 1));
                test3[1] = coord[i][1];
                test3[2] = coord[i][2] - curr_radial_displace ;
                Rotate_around_arbitrary_axis(
		   rot_axis,curr_rotation,test3,coord[i]);
                coord[i][1] -= (((double) (i - 1)) * bp_drop_per_step);
                test3[0] = coord[tx_unit[1] * 2 - (i - 1)][0];
                test3[0] -= ((double) (i - 1));
                test3[1] = coord[tx_unit[1] * 2 - (i - 1)][1];
                test3[2] = coord[tx_unit[1] * 2 - (i - 1)][2] -
			   curr_radial_displace;
                Rotate_around_arbitrary_axis(rot_axis,curr_rotation,
                             test3,coord[tx_unit[1] * 2 - (i - 1)]);
                coord[tx_unit[1] * 2 - (i - 1)][1] -=
                   (((double) (i - 1)) * bp_drop_per_step);
        }

        test2[0] = 1.0;
        test2[1] = 0.0;
        test2[2] = 0.0;
        test3[0] = coord[tx_unit[1]][0];
        test3[1] = 0.0;
        test3[2] = coord[tx_unit[1]][2];
printf("XZ POS OF END WAS: %f %f\n",test3[0],test3[2]);
        spiral_rotation = Find_angle_between_vectors(test2,test3);
        spiral_rotation /= (2 * PI);
        spiral_rotation *= 360.0;
	if (test3[2] < 0.0) {
		spiral_rotation *= -1.0;
	}
printf("WILL ROTATE %f\n",spiral_rotation);
        for (i = 1; i <= tx_unit[1] * 2; i ++) {
                for (j = 0 ; j < 3; j++) {
                        test1[j] = coord[i][j];
                }
                Rotate_around_arbitrary_axis(
		   rot_axis,spiral_rotation,test1,coord[i]);
        }
    }
printf("XZ POS AFTER ROTATION IS: %f %f\n",coord[tx_unit[1]][0],coord[tx_unit[1]][2]);
       for (i = 1; i <= tx_unit[1]; i ++) {
                for (j = 0 ; j < 3; j++) {
                        rnap_pos[i][j] = (coord[i][j] +
                                          coord[tx_unit[1] * 2 -
					  (i - 1)][j]) / 2.0;
                }
                rnap_rot[i] += spiral_rotation;
        }
}


int Assign_coordinates_to_DNA_template(int *tx_unit, double coord[][3],
				       int bonds[][2], int num_rna)
{
        int i,j,ctr=0;
	int offset;
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
	for (i = 2; i <= 30; i++) {
		bonds[(i - 1) + (2 * (tx_unit[1] - 1))][0] =
			(i - 1) + promoter_offset;
		bonds[(i - 1) + (2 * (tx_unit[1] - 1))][1] =
			i + promoter_offset;
		bonds[29 + (i - 1) + (2 * (tx_unit[1] - 1))][0] =
			30 + (i - 1) + promoter_offset;
		bonds[29 + (i - 1) + (2 * (tx_unit[1] - 1))][1] =
			30 + i + promoter_offset;
		bonds[58 + (i - 1) + (2 * (tx_unit[1] - 1))][0] =
			60 + (i - 1) + promoter_offset;
		bonds[58 + (i - 1) + (2 * (tx_unit[1] - 1))][1] =
			60 + i + promoter_offset;
		bonds[58 + 29 + (i - 1) + (2 * (tx_unit[1] - 1))][0] =
			90 + (i - 1) + promoter_offset;
		bonds[58 + 29 + (i - 1) + (2 * (tx_unit[1] - 1))][1] =
			90 + i + promoter_offset;
	}
	offset = 2 * (tx_unit[1] - 1) + 2 * 29 + 2 * 29;
	for (i = 1; i <= num_rna; i++) {
		bonds[offset + i][0] = 2 * tx_unit[1] + ((i * 2) - 1);
		bonds[offset + i][1] = 2 * tx_unit[1] + (i * 2);
	}
	offset += num_rna;
	ctr = 0;
	for (i = 1; i <= num_rna * num_genes * 2; i += 2) {
		ctr++;
		bonds[offset + ctr][0] = ribostart_offset + i;
		bonds[offset + ctr][1] = ribostart_offset + i + 1;
	}

        return(ctr);
}


void Send_atoms_offscreen(int start, int stop, double coord[][3])
{
	int i,j;
	for (i = start; i <= stop; i++) {
		for (j = 0; j < 3; j++) {
			coord[i][j] = -999.999;
		}
	}
	return;
}

void Assign_RNAP_positions(char *line, int num_rna, int *tx_unit,
			   double coord[][3], int *rnap_assignments,
			   int *last_rnap_use, int *num_last_rnap,
                           int *current_rnap_use, int *num_current_rnap,
			   int max_rnap)
{
	int i,j,num_rnap=0;
	int junk1,junk2;
	int offset,index[200],rnap_slot[200];
	int bp[100];
	double time;
	
	offset = 2 * tx_unit[1] + 2 * num_rna;
	sscanf(line,"%d %lf %d %d",&junk1,&time,&junk2,&num_rnap);
	char *str = line;
	str += 31;
	num_RNAP_check = num_rnap;
	for (i = 1; i <= num_rnap; i++) {
		sscanf(str,"%d %d",&index[i],&bp[i]);
		RNAP_BT_checklist[i] = bp[i];
		if (rnap_assignments[index[i]] != 0) {
			rnap_slot[i] = rnap_assignments[index[i]];
		}
		else {
			for (j = 1; j <= max_rnap; j++) {
				if (last_rnap_use[j] == 0 &&
				    current_rnap_use[j] == 0) {
					rnap_slot[i] = j;
					rnap_assignments[index[i]] = j;
					break;
				}
			}
		}
		current_rnap_use[rnap_slot[i]]++;
		str += 14;
	}
	
	for (i = 1; i <= num_rnap; i++) {
		for (j = 0; j < 3; j++) {
			coord[offset + rnap_slot[i]][j] =
				rnap_pos[bp[i]][j];
		}
	}
	return;
}

void Assign_RNA_and_ribosome_positions(char *line, int num_rna,
				       int *tx_unit, double coord[][3],
                                       int max_rnap, int *last_ribo_use,
				       int *num_last_ribo,
				       int *current_ribo_use,
				       int *num_current_ribo,
                                       int *ribo_assignments,
				       int *rna_key, int max_ribo,
                                       int frame_rank[][MAX_ATOMS],
				       int frame)
{
	int i,j,k,num_ribo=0;
	int junk1,begin,end,label,connect;
	int offset,index[1000],ribo_slot[1000];
	int offset2,bt_offset;
	int bp[1000],closest;
	int gene_exists,gene_start,gene_stop=0;
	double time;
	double x,y,z,test1[3],test2[3];
	
	offset = 2*tx_unit[1];
	sscanf(line,"%d %lf %d %d %d %d %d",
	       &junk1,&time,&label,&begin,&end,&connect,&num_ribo);
	char *str = line;
	str += 52;
	for (i = 1; i <= num_ribo; i++) {
		sscanf(str,"%d %d",&index[i],&bp[i]);
		if (ribo_assignments[index[i]] != 0) {
			ribo_slot[i] = ribo_assignments[index[i]];
		}
		else {
			for (j = 1; j <= max_ribo; j++) {
				if (last_ribo_use[j] == 0 &&
				    current_ribo_use[j] == 0) {
					ribo_slot[i] = j;
					ribo_assignments[index[i]] = j;
					frame_info[frame].num_new_ribo++;
					frame_info[frame].new_ribo_index1[
					 frame_info[frame].num_new_ribo] =
					  (2 * tx_unit[1]) + (2 * num_rna) +
					   max_rnap + j; 
					frame_info[frame].new_ribo_index2[
					 frame_info[frame].num_new_ribo] =
					  (2 * tx_unit[1]) + (2 * num_rna) +
					   max_rnap + max_ribo + j;
					closest = 999999;
					for (k = 1; k <= num_genes; k++) {
					   if (abs(bp[i] - tsl_start[k]) <
					       closest) {
						 closest = abs(bp[i] -
							   tsl_start[k]);
						 frame_info[frame].new_ribo_gene[
						  frame_info[frame].num_new_ribo] = k;
					   }
					}
					break;
				}
			}
		}
		current_ribo_use[ribo_slot[i]]++;
		str += 14;
	}

	if (connect == 0) {								
		x = coord[tx_unit[1]][0] + (45.0 * frame_rank[junk1][label]);
		z = 0.0;
	}	
	else {
		bt_offset = 9999;
printf("CHECKING %d RNAPS: ",num_RNAP_check);
		for (i = 1; i <= num_RNAP_check; i++) {
printf("RNAP %d AT %d; ",i,RNAP_BT_checklist[i]);
		    if (RNAP_BT_checklist[i] > end + 1) {				// For post-transloc state
				continue;
		    }
		    bt_offset = (bt_offset < end - RNAP_BT_checklist[i] ? 
				 bt_offset : end - RNAP_BT_checklist[i]);
		}


printf("\n");
printf("FRAME IS %d; END WAS %d; BT OFFSET IS %d; END REASSIGNED TO %d\n",frame,end,bt_offset,end-bt_offset);

		end -= bt_offset;
		x = rnap_pos[end][0];
		z = rnap_pos[end][2];
		if (bt_offset == -1) {							
			end--;
		}
	}

	y = rnap_pos[end][1];
	offset += ((rna_key[label] - 1) *  2);		

	offset2 = ribostart_offset +
		(rna_key[label] - 1) * 2 * num_genes;

	coord[offset + 1][0] = x;							// RNA now two points with bond
	coord[offset + 1][1] = y + ((double) (end - begin));				// ...made above
	coord[offset + 1][2] = z;
	coord[offset + 2][0] = x;
	coord[offset + 2][1] = y;
	coord[offset + 2][2] = z;

	for (i = 1; i <= num_genes; i++) {
		gene_exists = 0;
		for (j = tsl_start[i]; j <= tsl_stop[i]; j++) {
			if (j >= begin && j <= end) {
				gene_exists++;
				gene_start = j;
				break;
			}
		}
		if (gene_exists) {
			for (j = tsl_stop[i]; j >= tsl_start[i]; j--) {
				if (j >= begin && j <= end) {
					gene_stop = j;
					break;
				}
			}
			coord[offset2 + (i - 1) * 2 + 1][0] = x;
			coord[offset2 + (i - 1) * 2 + 1][1] = 
			 y + ((double) (end - gene_start));
			coord[offset2 + (i - 1) * 2 + 1][2] = z;
			coord[offset2 + (i - 1) * 2 + 2][0] = x;
			coord[offset2 + (i - 1) * 2 + 2][1] = 
			 y + ((double) (end - gene_stop));
			coord[offset2 + (i - 1) * 2 + 2][2] = z;
		}
	}
	
	offset = (2 * tx_unit[1]) + (2 * num_rna) + max_rnap;
	test2[0] = 0.0;
	test2[1] = 1.0;
	test2[2] = 0.0;
	for (i = 1; i <= num_ribo; i++) {
		test1[0] = 10.0;
		test1[1] = y + ((double) (end - bp[i]));
		test1[2] = 0.0;
		if (connect) {
		  Rotate_around_arbitrary_axis(test2,rnap_rot[end],test1,
			coord[offset + ribo_slot[i]]);
		}
		else {
			for (j = 0; j < 3; j++) {
				coord[offset + ribo_slot[i]][j] = test1[j];
			}
		}
		coord[offset + ribo_slot[i]][0] += x;
		coord[offset + ribo_slot[i]][2] += z;
		test1[0] = -7.5;
		test1[1] = y + ((double) (end - bp[i]));
		test1[2] = 0.0;
		if (connect) {
		  Rotate_around_arbitrary_axis(test2,rnap_rot[end],test1,
			coord[offset + max_ribo + ribo_slot[i]]);
		}
		else {
			for (j = 0; j < 3; j++) {
				coord[offset + max_ribo + ribo_slot[i]][j] =
				  test1[j];
			}
		}
		coord[offset + max_ribo + ribo_slot[i]][0] += x;
		coord[offset + max_ribo + ribo_slot[i]][2] += z;
	}
	return;
}

void Assign_protein_positions(int frame)
{
	int i;
	int offset;

	offset = protein_offset;
	offset += num_proteins;
	if (!(promoter_state[frame])) {
		offset += 60;
	}
	for (i = 1; i <= 30; i++) {
		coord[i + offset][0] = coord[i][0];
		coord[i + offset][1] = coord[i][1];
		coord[i + offset][2] = coord[i][2];
		coord[i + offset + 30][0] = coord[i+(2*tx_unit[1]-30)][0];
		coord[i + offset + 30][1] = coord[i+(2*tx_unit[1]-30)][1];
		coord[i + offset + 30][2] = coord[i+(2*tx_unit[1]-30)][2];
	}
	return;
}
		


void Write_frame(XDRFILE *xfp_out, int frame, float time, int num_atoms,
		 matrix box, float prec, double coord[][3])
{
	int i,j,crap;
	rvec xrvec[MAX_ATOMS];

	for (i = 1; i <= num_atoms; i++) {
		for (j = 0; j < 3; j++) {
			xrvec[i-1][j] = coord[i][j]/10.0;
		}
	}
	crap = write_xtc(xfp_out, num_atoms, frame, time, box, xrvec, prec);
	crap++;
	return;
}
	
void Write_pdb(int *tx_unit,int num_rna, int max_rnap, int max_ribo,
	       double coord[][3])
{
	int i;
	int ctr,which_gene;
	char struct_nsme[200];
	
	sprintf(struct_nsme,"tx_tsl.phosDNA_march.two_lobe_ribo.proteins."
			    "macrohelix.multigene.%s.pdb",
			    name_tag);
	FILE *fp;
	fp = fopen(struct_nsme,"w");

	for (i = 1; i <= 2*tx_unit[1]; i++) {
		fprintf(fp, "ATOM%7d  P   DNA A%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,1,coord[i][0],coord[i][1],coord[i][2]);
	}

	for (i = 2*tx_unit[1] + 1; i <= 2 * num_rna + 2 * tx_unit[1]; i++) {
		fprintf(fp, "ATOM%7d  N   DNA B%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,2,coord[i][0],coord[i][1],coord[i][2]);
	}


	for (i = (2 * num_rna) + (2 * tx_unit[1]) + 1;
	     i <= (2 * num_rna) + (2 * tx_unit[1]) + max_rnap;
             i++) {
		fprintf(fp, "ATOM%7d  C   DNA C%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,3,coord[i][0],coord[i][1],coord[i][2]);
	}


	for (i = (2 * num_rna) + (2 * tx_unit[1]) + max_rnap + 1;
             i <= (2 * num_rna) +
		  (2 * tx_unit[1]) + max_rnap + max_ribo + max_ribo;
	     i++) {
	      if (i <= (2 * num_rna) +
		  (2 * tx_unit[1]) + max_rnap + max_ribo) {
		fprintf(fp, "ATOM%7d  O   DNA D%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,4,coord[i][0],coord[i][1],coord[i][2]);
	      }
	      else {
		fprintf(fp, "ATOM%7d  S   DNA D%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,4,coord[i][0],coord[i][1],coord[i][2]);
	      }
	}

	for (i = protein_offset + 1;
	     i <= protein_offset + num_proteins;
	     i++) {
		fprintf(fp, "ATOM%7d  CA  DNA C%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,3,coord[i][0],coord[i][1],coord[i][2]);
	}

	for (i = promoter_offset + 1; i <= promoter_offset + 60; i++) {
		fprintf(fp, "ATOM%7d  CB  DNA C%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,3,coord[i][0],coord[i][1],coord[i][2]);
	}

	for (i = promoter_offset + 61; i <= promoter_offset + 120; i++) {
		fprintf(fp, "ATOM%7d  CG  DNA C%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,3,coord[i][0],coord[i][1],coord[i][2]);
	}

	ctr = 0;
	for (i = ribostart_offset + 1; i <= num_atoms; i++) {
		ctr++;
		which_gene = ctr % (num_genes * 2);
		if (which_gene == 0) {
			which_gene = num_genes * 2;
		}
		which_gene++;
		which_gene /= 2;
		fprintf(fp, "ATOM%7d  %s  DNA C%4d%12.3f%8.3f%8.3f  "
			    "1.00  0.00\n",
                            i,gene_label[which_gene],3,coord[i][0],
			    coord[i][1],coord[i][2]);
	}
		
	fclose(fp);
	return;
}

void Assign_positions(char *name, int num_atoms, int *tx_unit,
		      double coord[][3], int *rna_key, int num_rna,
		      int max_rnap, int max_ribo,
		      int frame_rank[][MAX_ATOMS])
{
	int frame,previous = 0;
	int frame_ctr=0;
	int identity,num;
	float time;
	float prec = 1000.0;
	matrix box;
	char line[1000];
	int ribo_assignments[MAX_ATOMS] = {0};
	int rnap_assignments[MAX_ATOMS] = {0};
	int num_current_rnap, current_rnap_use[200] = {0};
	int num_last_rnap, last_rnap_use[200] = {0};
	int num_current_ribo, current_ribo_use[1000] = {0};
	int num_last_ribo, last_ribo_use[1000] = {0};
	char x_name[200];

	sprintf(x_name,"tx_tsl.phosDNA_march.two_lobe_ribo.proteins."
		       "macrohelix.multigene.%s.xtc",
			name_tag);
	XDRFILE *xfp_out = xdrfile_open(x_name, "w");

	FILE *fp;
	fp = fopen(name,"r");
	while (fgets(line,1000,fp) != NULL) {
		sscanf(line,"%d %f %d %d",&frame,&time,&identity,&num);
		if (frame > 1 && frame != previous) {
			if (frame == 2) {
				Write_pdb(tx_unit,num_rna,max_rnap,
					  max_ribo,coord);
			}
			Write_frame(xfp_out,frame,time,num_atoms,box,
				    prec,coord);
			frame_ctr++;
			Send_atoms_offscreen(2 * tx_unit[1] + 1, num_atoms,
					     coord);
			memcpy(last_rnap_use,current_rnap_use,
				(max_rnap + 1) * sizeof(current_rnap_use[0]));
			num_last_rnap = num_current_rnap;
			num_current_rnap = 0;
			memset(current_rnap_use,0,sizeof(current_rnap_use));
			memcpy(last_ribo_use,current_ribo_use,
				(max_ribo + 1) * sizeof(current_ribo_use[0]));
			num_last_ribo = num_current_ribo;
			num_current_ribo = 0;
			memset(current_ribo_use,0,sizeof(current_ribo_use));
		}
		if (identity == 0) {
		   if (num > 0) {
			Assign_RNAP_positions(line,num_rna,tx_unit,coord,
					      rnap_assignments,
					      last_rnap_use,&num_last_rnap,
					      current_rnap_use,
                                              &num_current_rnap, max_rnap);
		   }
		   else {
			num_RNAP_check = 0;
		   }
		}
		if (identity >= 1) {
			Assign_RNA_and_ribosome_positions(
				line,num_rna,tx_unit,coord,max_rnap,
				last_ribo_use,&num_last_ribo,
				current_ribo_use,&num_current_ribo,
                                ribo_assignments,rna_key,max_ribo,
				frame_rank,frame);
		}
		if (identity == -1) {
			Assign_protein_positions(frame);
		}
		previous = frame;
	}
	fclose(fp);
	xdrfile_close(xfp_out);

	printf("FINSHED WRITING %d FRAMES\n", frame_ctr);

	return;
}


void Write_bonds(int num_bonds, int num_atoms, int bonds[][2])
{
        int i;
        char name[200];

        sprintf(name,"tx_tsl.phosDNA_march.two_lobe_ribo.proteins."
		     "macrohelix.multigene.%s.psf",
		     name_tag);
        FILE *fp;
        fp = fopen(name,"w");
        fprintf(fp,"PSF whole-genome\n\n       1 !NTITLE\n REMARKS "
		   "psf for tx/transl snapshot\n\n");
        fprintf(fp,"%8d !NATOM\n", num_atoms);
        for (i = 1; i <= num_atoms; i++) {
                fprintf(fp,"%8d%2s%5d    DNA  N    NH3   -0.300000       "
			   "14.0070           0\n",
			   i,"U",1);
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

void Write_frame_info_for_tcl_script()
{
	int i,j,time_int;
	int min,sec,hundredth_int;
	double hundredth;
	FILE *fp;

	fp = fopen("temp_traj_fn.txt","w");
	fprintf(fp,"proc load_traj_line_array {} {\n");
	fprintf(fp,"  global lines\n");
	fprintf(fp,"  global num_lines\n\n");

	for (i = 1; i <= num_frames; i++) {
		time_int = (int) (frame_info[i].time + 0.0001);
		min = time_int/60;
		sec = time_int % 60;
		hundredth = frame_info[i].time - ((double) time_int);
		hundredth_int = (int) ((hundredth+0.00001)*100);
		fprintf(fp,"  set lines(%d) \"%8d%8d%8d%8d%8d%8d%8d%8d%8d",
			i,i,min,sec,hundredth_int,frame_info[i].num_rna,
			frame_info[i].num_rna - frame_info[i].num_mature,
			frame_info[i].num_mature,
			frame_info[i].num_ribo, num_genes);
		for (j = 1; j <= num_genes; j++) {
			fprintf(fp,"%8d",frame_info[i].num_proteins[j]);
		}
		fprintf(fp,"%8d",frame_info[i].num_new_ribo);
		for (j = 1; j <= frame_info[i].num_new_ribo; j++) {
			fprintf(fp,"%8d%8d%8d",
				frame_info[i].new_ribo_index1[j],
				frame_info[i].new_ribo_index2[j],
				frame_info[i].new_ribo_gene[j]);
		}
		fprintf(fp,"\"\n");
	}

	fprintf(fp,"  set num_lines %d\n}\n",num_frames);
	fclose(fp);

	return;
}

void Make_custom_tcl_script()
{
	int i,syst_ret;
	int max_mature = 0;
	double x_offset,y_offset;
	double x_min = 99999.9,x_max = -999999.9;
	char command[1000],tcl_template[1000000];

	y_offset = -1.0*bp_drop_per_step*((double) tx_unit[1]);
	for (i = 1; i <= tx_unit[1]; i++) {
		if (rnap_pos[i][0] < x_min) {
			x_min = rnap_pos[i][0];
		}
		if (rnap_pos[i][0] > x_max) {
			x_max = rnap_pos[i][0];
		}
	}
	for (i = 1; i <= num_frames; i++) {
		if (frame_info[i].num_mature > max_mature) {
			max_mature = frame_info[i].num_mature;
		}
	}
	x_max += ((double) (max_mature * 45));
	x_offset = (((x_min+x_max)/2.0) - 270.0) * -1.0;

	Assign_tcl_template_string(tcl_template);

	FILE *fp;
	fp = fopen("template.tcl","w");
	fprintf(fp,"%s\n",tcl_template);
	fclose(fp);

	sprintf(command,"sed \'s/LEFT_EDGE/%.1f/g\' %s > temp_rep.ribo.tcl",
		x_offset,"template.tcl");
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }


	sprintf(command,"sed -i \'s/Y_START/%.1f/g\' temp_rep.ribo.tcl",
		y_offset);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"sed -i \'s/MIN_INDEX/%d/g\' temp_rep.ribo.tcl",
		(2 * tx_unit[1]) + (2 * num_rna) + max_rnap + 1);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"sed -i \'s/MAX_INDEX/%d/g\' temp_rep.ribo.tcl",
		(2 * tx_unit[1]) + (2 * num_rna) + max_rnap + 2 * max_ribo);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"sed -i \'s/TRAJIC/%s/g\' temp_rep.ribo.tcl",
		name_tag);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"cat temp_rep.ribo.tcl temp_traj_fn.txt >"
			" custom_script_for_%s.tcl",
			name_tag);
	syst_ret = system(command); 
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"mkdir MOVIE_MATERIALS.%s", name_tag);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	sprintf(command,"mv tx_tsl.phosDNA_march.two_lobe_ribo.proteins."
			"macrohelix.multigene.%s* MOVIE_MATERIALS.%s; "
			"mv custom_script_for_%s.tcl MOVIE_MATERIALS.%s ;"
			"rm temp_rep.ribo.tcl; rm temp_traj_fn.txt",
			name_tag,name_tag,name_tag,name_tag);
	syst_ret = system(command);
        if (syst_ret == -1) {
                printf("FAILED TO EXECUTE %s\n",command);
        }

	return;
}

int main (int argc, char *argv[])
{
	int i;
	char sim_input_name[200], name[200];

	strcpy(name,argv[1]);
	strcpy(sim_input_name,argv[2]);
	start_bp_radius = strtod(argv[3],NULL);
	radius_incr_per_ring = strtod(argv[4],NULL);
	bp_drop_per_step = strtod(argv[5],NULL);
	strcpy(name_tag,argv[6]);

	memset(frame_info,0,sizeof(frame_info));

	Read_input_for_gene_info(sim_input_name);
	printf("TRANSCRIPT CONTAINS %d GENES\n",num_genes);
	printf("START AND STOP POSITIONS:\n");
	for (i = 1; i <= num_genes; i++) {
		printf("GENE%3d START:%7d; STOP: %7d\n",
			i, tsl_start[i],tsl_stop[i]);
	}

	num_frames = Get_starting_trajectory_information(
				name,&num_rna,&max_rnap,&max_ribo,
                                sim_time,rna_key,&tx_end,frame_rank);

	printf("TRAJ INFO:\nNUMBER OF FRAMES: %d\nNUMBER OF RNA: %d\n"
		"MAX RNAP: %d\nMAX RIBOSOMES: %d\n",
                num_frames,num_rna,max_rnap,max_ribo);

	printf("TRANSCRIPTION BOUNDARIES: 1 %d\n",tx_end);
	for (i = 1; i <= num_genes; i++) {
		printf("%d PROTEINS MADE DURING TRAJECTORY\n",
			frame_info[num_frames].num_proteins[i]);
	}
        printf("KEY:\n");
        for (i = 1; i <= 5; i++) {
                printf("%10d%10d\n",i,rna_key[i]);
        }

	tx_unit[0] = 1;
	tx_unit[1] = tx_end;

	num_proteins = 0;
	memset(proteins,0,sizeof(proteins));

	protein_offset = (2 * tx_unit[1]) + max_rnap +
			 (2 * num_rna) + max_ribo + max_ribo;
	promoter_offset = protein_offset + num_proteins;
	ribostart_offset = promoter_offset + 120;
	num_atoms = ribostart_offset + (2 * num_genes * num_rna);

	num_dna = Assign_coordinates_to_DNA_template(
			tx_unit,coord,bonds,num_rna);

	Build_helical_track();

	num_bonds = 2 * (tx_unit[1] - 1) + 2 * 29 + 2 * 29 +
		    (1 + num_genes) * num_rna;

	Write_bonds(num_bonds,num_atoms,bonds);

	Assign_positions(name,num_atoms,tx_unit,coord,rna_key,num_rna,
			 max_rnap, max_ribo,frame_rank);


	Write_frame_info_for_tcl_script();

	Make_custom_tcl_script();

	return(0);
}
	
