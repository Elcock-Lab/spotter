#ifndef VEC_UTIL_H
#define VEC_UTIL_H

void Cross_product(double *test1, double *test2, double *cross);

double Dot_product(double *test1, double *test2);

double Find_angle_between_vectors(double *test1, double *test2);

double Find_interatom_distance (double *atom1, double *atom2);

void Translate(double *pos, double *move, double *moved);

void Rotate_around_arbitrary_axis (double *axis, double angle,
				   double *orig, double *moved);

#endif	
