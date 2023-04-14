#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define MAX_ATOMS 50000
#define MAX_FRAMES 10000
#define PI 3.141592654


void Cross_product(double *test1, double *test2, double *cross)
{

        cross[0] = (test1[1] * test2[2]) - (test1[2] * test2[1]);               
        cross[1] = -1 * ((test1[0] * test2[2]) - (test1[2] * test2[0]));       
        cross[2] = (test1[0] * test2[1]) - (test1[1] * test2[0]);         

        return;
}

double Dot_product(double *test1, double *test2)
{
        int i;
        double dot;

        dot = 0.0;
        for (i = 0; i < 3; i++) {
                 dot += (test1[i] * test2[i]);
        }
        return dot;
}

double Find_angle_between_vectors(double *test1, double *test2)
{
        double mag1, mag2, dot, angle;

        dot = Dot_product(test1,test2);
        mag1 = sqrt((test1[0] * test1[0]) + (test1[1] * test1[1]) +
		    (test1[2] * test1[2]));
        mag2 = sqrt((test2[0] * test2[0]) + (test2[1] * test2[1]) +
		    (test2[2] * test2[2]));
        angle = acos(dot / (mag1 * mag2));
        return angle;
}


double Find_interatom_distance (double *atom1, double *atom2)
{
        int i;
        double dist;
        dist=0.0;
        for (i=0; i < 3; i++) {
                dist += (pow((atom1[i] - atom2[i]), 2));
        }
        dist = sqrt(dist);
        return dist;
}

void Translate(double *pos, double *move, double *moved)
{
        int i;
        for (i = 0; i < 3; i++) {
                moved[i] = pos[i] + move[i];
        }
        return;
}

void Rotate_around_arbitrary_axis (double *axis, double angle,
				   double *orig, double *moved)
{
        int i,j;
        double x,y,z,rotation_matrix[3][3], sin_theta, cos_theta, magnitude;

        angle /= 360;
        angle *= (2 * PI);
        sin_theta = sin(angle);
        cos_theta = cos(angle);
        magnitude = sqrt((axis[0] * axis[0]) + (axis[1] * axis[1]) +
			 (axis[2] * axis[2]));
        x = axis[0]/magnitude;
        y = axis[1]/magnitude;
        z = axis[2]/magnitude;
        rotation_matrix[0][0] = cos_theta + (x * x * (1 - cos_theta));
        rotation_matrix[0][1] = (x * y * (1 - cos_theta)) - (z * sin_theta);
        rotation_matrix[0][2] = (x * z * (1 - cos_theta)) + (y * sin_theta);
        rotation_matrix[1][0] = (y * x * (1 - cos_theta)) + (z * sin_theta);
        rotation_matrix[1][1] = cos_theta + (y * y * (1 - cos_theta));
        rotation_matrix[1][2] = (y * z * (1 - cos_theta)) - (x * sin_theta);
        rotation_matrix[2][0] = (z * x * (1 - cos_theta)) - (y * sin_theta);
        rotation_matrix[2][1] = (z * y * (1 - cos_theta)) + (x * sin_theta);
        rotation_matrix[2][2] = cos_theta + (z * z * (1 - cos_theta));

        for (i = 0; i < 3; i++) {
                moved[i] = 0.0;
                for (j = 0; j < 3; j++) {
                        moved[i] += (orig[j] * rotation_matrix[i][j]);
                }
        }
}

	
