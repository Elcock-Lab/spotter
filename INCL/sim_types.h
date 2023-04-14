#ifndef SIM_TYPES_H
#define SIM_TYPES_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


typedef struct {
        int start;
        int stop;
        float k_init;
	float mean_rate;
} gene_form;

typedef struct {
        int object_type;
        int object_index;
        int gene_index;
        int reaction_id;
} reaction_form;

typedef struct {
	double net;		// Total dwelltime at position
	double mu_f;		// On-pathway forward translocation
	double mu_b;		// On-pathway reverse translocation
	double mu_p;		// Entry into off-pathway state
	double mu_f_off;	// Off-pathway forward transloc
	double mu_b_off;	// Off-pathway reverse transloc
	double mu_net_pre;	// Net rate for pretransloc rxns
	double mu_net_off;	// Net rate for offpath rxns
	double p_on_f;		// Rate-dependent prob of f transloc
	double p_off_f;		// Prob of for (v. rev) offpath step
} rnap_rate_form;		// Other rates invariant by position...

typedef struct {
	double codon;		// Expected dwell based on codon (Fluitt et al.)
	double exp;		// Riboseq-observed dwelltime
	double mu_f;		// Forward translocation
	double mu_add;		// Decoding step
} ribo_rate_form;

typedef struct {
        int promo_id;
        double add_time;
} promo_form;

typedef struct {		// Priority queue for reactions
        int reaction;
        double reaction_time;
} tree_form;

typedef struct {
        int RNAP_bound;
        int RNAP_state;
        int RNAP_rot;
        float sigma;
} kymo_form;

typedef struct {
        double last_p1_mu;
        double p1_time_adj;
        double last_e1_mu;
        double e1_time_adj;
        int checkpos;
        double checkstart;
        double last_resist_mu;
        double last_assist_mu;
        double last_f_mu;
        double f_time_adj;
} last_rate_form;

typedef struct {
	float start_tx;
	float end_tx;
	float start_deg;
	float end_deg;
	int ribo_start_nascent[25];
	int ribo_end_nascent[25];
	int ribo_start_mature[25];
	int ribo_end_mature[25];
} rna_log_form;

#endif

