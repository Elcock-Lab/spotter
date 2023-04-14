#ifndef GENDEF_H
#define GENDEF_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/inotify.h>
#include <unistd.h>

#define EVENT_SIZE (sizeof (struct inotify_event))
#define BUF_LEN (32 * (EVENT_SIZE + 16) )

#define MAX_TRAJ 2						
#define MAX_GENE_LEN 2500
#define MAX_RNAP_PER_TRAJ 10000
#define MAX_RNAP 2500
#define MAX_PROMO_CYCLES 1000
#define MAX_LOADING_ATTEMPTS 20000
#define MAX_RNA_PER_SNAP 1000
#define MAX_FISH_TIMES 50000
#define MAX_TEN_MIN_WINDOWS 10
#define MAX_FISH_PROBES 5
#define MAX_GENES 20
#define MAX_PROMO 33
#define PROMO_INFO 15								// Dimensions for object arrays:
#define RNA_INFO 150								// Adjust if necessary
#define TSL_INFO 10
#define RNAP_INFO 15
#define RIBO_INFO 15
#define MAX_TOPO 2500
#define TOPO_INFO 15
#define MAX_ACTIONS 50000
#define MAX_STEPS 1000000
#define GENOME 4639375
#define ORI 3923883
#define PI 3.1415927

#define ON 1
#define OFF 0
#define KILLED -1
#define NO_GENE 0
#define RNAP_ID 1
#define RNA_ID 2
#define RIBO_ID 3
#define HYBRID 2

#define POST_TRANSLOC 0
#define PRE_TRANSLOC 1								// Note: these are used for both
#define OFF_P1 2								// ribosomes and RNAP (STATE index)
#define OFF_P2 3
#define OFF_P3 4
#define CLOSED 5

#define INDEPENDENT 1								// Identity of "cause" for arrival
#define TRAILING_RNAP 2								// in the pretranslocated state
#define TRAILING_RIBO 3								// the "auto_forward" might be used
#define AUTO_FORWARD 4								// to prevent a cascade of pauses
#define UNKNOWN 5								// upstream of a paused RNAP

#define FORWARD 0
#define REVERSE 1

#define UPSTREAM 0
#define DOWNSTREAM 1
#define NEITHER 2
#define MIDDLE 3
#define NONE 4

#define BOTH 2

#define ASSIST_TORQUE_RELAX 7
#define RESIST_TORQUE_RELAX 8

#define IMPLICIT_RNAP_POS 1
#define IMPLICIT_EXTANT_RNA 2
#define EXPLICIT_RIBOSOMES 3
#define NO_RIBO_RNA_ONLY 4
#define FLAT_RESISTANCE 5

#define LEFT 1
#define RIGHT 2

#define TOPO_IA 0
#define GYRASE 1

#endif
