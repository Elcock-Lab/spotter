#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "INCL/RNAP_rate_generator.h"
#include "INCL/sequence_utilities.h"
#include "INCL/read_utilities.h"

double kp,kp2,kp3,ke1,ke2,ke3;
double kN,kc,kd;
double kf,kb;
double p2,p3;
double p01,p12,p10,p23,p21;
double P03,Ptot,Pmin,Vmax;
double q0,q1,q2,q3;
double k0,k1,k2,k3;

double sum,tot_dwell,super_tot_dwell;
double tau_return,tau_repeat,q_repeat;
double loop_add;

int num_trials;
double x[500],y[500];

typedef struct {
	double dt;
	double kf;
	double kb;
	double kp;
	double kf_off;
	double kb_off;
} rate_form;

rate_form f_rates[5000000];

void Read_parameters(char *name)
{
        char label[100],line[200];
	double value;
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
		sscanf(line,"%s %lf",label,&value);
		if (strcmp(label,"kp2") == 0) {
			kp2 = value;
		}
		if (strcmp(label,"kp3") == 0) {
			kp3 = value;
		}
		if (strcmp(label,"ke1") == 0) {
			ke1 = value;
		}
		if (strcmp(label,"ke2") == 0) {
			ke2 = value;
		}
		if (strcmp(label,"ke3") == 0) {
			ke3 = value;
		}
		if (strcmp(label,"kN") == 0) {
			kN = value;
		}
		if (strcmp(label,"kc") == 0) {
			kc = value;
		}
		if (strcmp(label,"kd") == 0) {
			kd = value;
		}
	}
	fclose(fp);
}

void Assign_default_parameters()
{
	kp2 = 0.0;
	kp3 = 0.0;
	ke1 = 1.77;
	ke2 = 0.241;
	ke3 = 0.01;
	kN  = 250000000;
	kc  = 500;
	kd  = 1000;
	
	return;
}	

int Read_rate_file(char *name, rate_form *rate)
{
	int ctr = 0;
	char line[200];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
		ctr++;
		char *str = line;
		str += 46;
		sscanf(str,"%lf %lf %lf %lf", &rate[ctr].kf,
			&rate[ctr].kb,&rate[ctr].kf_off,
			&rate[ctr].kb_off);
	}
	fclose(fp);
	return(ctr);
}

int Read_dwell_file(char *name, rate_form *rate)
{
	int ctr = 0,junk;
	char line[200];
        FILE *fp;
        fp = fopen(name,"r");
        while (fgets(line,200,fp) != NULL) {
		ctr++;
		sscanf(line,"%d %lf",&junk,&rate[ctr].dt);
	}
	fclose(fp);
	return(ctr);
}
	

double Inter_or_extrapolate(int num_trials, double target)
{
        int i;
        int glb_index=0,lub_index=0;
        double glb=0.0,lub=10000000.0;
        double min=10000000.0,max=0.0;
        double m,rate;
        double x1=0,y1=0,x2=0,y2=0;

        for (i = 1; i <= num_trials; i++) {
                if (y[i] < target &&
                    y[i] > glb) {
                        glb = y[i];
                        glb_index = i;
                }
                if (y[i] > target &&
                    y[i] < lub) {
                        lub = y[i];
                        lub_index = i;
                }
        }
        if (glb_index != 0) {
                if (lub_index != 0) {
                        x1 = x[glb_index];
                        x2 = x[lub_index];
                        y1 = y[glb_index];
                        y2 = y[lub_index];
                }
                else {
                        x2 = x[glb_index];
                        y2 = y[glb_index];
                        for (i = 1; i <= num_trials; i++) {
                                if (i == glb_index) {
                                        continue;
                                }
                                if (y[i] > max) {
                                        max = y[i];
                                        x1 = x[i];
                                        y1 = y[i];
                                }
                        }
                }
        }
        else {
                x1 = x[lub_index];
                y1 = y[lub_index];
                for (i = 1; i <= num_trials; i++) {
                        if (i == lub_index) {
                                continue;
                        }
                        if (y[i] < min) {
                                min = y[i];
                                x2 = x[i];
                                y2 = y[i];
                        }
                }
        }
        m = (y2-y1)/(x2-x1);
        rate = (target - (y1-(m*x1)))/m;
        return(rate);
}

double Find_dwell_time(double kp)
{
	int i;

	p01 = kf/(kf+kp);
	p12 = kN/(kN+kb);							// These two kept
	p10 = 1 - p12;								// here for clarity
	P03 = (p01*p12*p23)/(1-(p10*p01)-(p12*p21));
	Ptot = 1- P03;
	Pmin = kp/(kf+kp);

	Vmax = ((kf+kp)*kc)/(kf+kp+kc);

	q0 = 1 - Ptot;
	q1 = Ptot*(1-(p2/(1-Ptot)));
	q2 = Ptot*((p2*(1-p3))/(1-Ptot));
	q3 = Ptot*((p2*p3)/(1-Ptot));
	q_repeat = 1 - (q0 + q1);

	k0 = 1/Vmax;
	k1 = (1/((1-Ptot+(Ptot*p2))*(ke1+kp2))) + (1/Vmax);
	k2 = (1/(ke2+kp3)) + 
	     (1/(((1-Ptot+Ptot*p2)/(1+Ptot-Ptot*p2))*(ke1+kp2))) +
	     (1/Vmax);
	k3 = k2+(1/ke3);

	sum = q0+q1+q2+q3;
	tot_dwell = k0*q0 + k1*q1 + k2*q2 + k3*q3;

	tau_repeat = tot_dwell + tau_return;
	loop_add = tau_repeat;
	super_tot_dwell = tot_dwell;
	for (i = 1; i <= 5; i++) {
		loop_add *= q_repeat;
		super_tot_dwell += loop_add;
	}
	return(super_tot_dwell);
}

int Find_offpathway_rate(double target)
{
	int i;
	int num_trials;
	double error;

	if ((1/kf)+(1/kc) > target) {
		kp = 0.0000001;
		return(0);
	}

	kp = 1;
	y[1] = Find_dwell_time(kp);
	x[1] = 1;

	kp = 5000;
	y[2] = Find_dwell_time(kp);
	x[2] = 5000;

	num_trials = 2;
	for (i = 3; i <= 100; i++) {
		kp = Inter_or_extrapolate(num_trials,target);
		if (kp < 0.0) {
			kp = 0.000001;
		}
		y[i] = Find_dwell_time(kp);
		error = fabs(y[i] - target);
		if (error < 0.0001) {
			break;
		}
		x[i] = kp;
		num_trials++;
	}
	return(num_trials);
}

void Write_final_rates(char *name, rate_form *rates, int num_bp)
{
	int i;
	FILE *fp;
	fp = fopen(name,"w");
	for (i = 1; i <= num_bp; i++) {
		fprintf(fp,"%10d%20.8f%20.8f%20.8f%20.8f%20.8f%20.8f\n",i,
			rates[i].dt,rates[i].kf,rates[i].kb,
			rates[i].kp,rates[i].kf_off,rates[i].kb_off);
	}
	fclose(fp);
	return;
}

void Calculate_final_rates(char *param_name, char *rate_name,
			   char *dwell_name)
{
	int i,j;
	int num_bp;
	double target;

	if (strcmp(
	       "dekker_derived_rnap_paramaters.kc_500.no_advanced_pauses",
	       param_name) != 0) {
		Read_parameters(param_name);
	}
	else {
		Assign_default_parameters();
	}

	p2 = kp2/(ke1+kp2);							// Need not be recalculated
	p3 = kp3/(ke2+kp3);			

	p23 = kc/(kc+kd);			
	p21 = 1 - p23;				

	num_bp = Read_rate_file(rate_name,f_rates);
	printf("FOUND %d ENTRIES IN RATE FILE\n",num_bp);

	num_bp = Read_dwell_file(dwell_name,f_rates);
	printf("FOUND %d ENTRIES IN DWELL FILE\n",num_bp);
	
	for (i = 1; i <= num_bp; i++) {
		tau_return = 0.0;
		for (j = 4; j >= 1; j--) {
			if ((i - j) > 0) {
				tau_return += f_rates[(i - j)].dt;
			}
		} 
		kf = f_rates[i].kf;
		kb = f_rates[i].kb;
		target = f_rates[i].dt;
		num_trials = Find_offpathway_rate(target);
		f_rates[i].kp = kp;
	}

	Write_final_rates("full_rnap_rate_set.seq_based",
			  f_rates,num_bp);

	return;
}	











	


