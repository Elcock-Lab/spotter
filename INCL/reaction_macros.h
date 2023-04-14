#ifndef REACTION_MACROS_H
#define REACTION_MACROS_H

#define MAX_REACTIONS 500000

#define PARENT(x) ((x)/(2))
#define FIRST_CHILD(x) ((2)*(x))
#define SECOND_CHILD(x) (((2)*(x))+(1))
#define MIN_TIME(x,y) (((x)<(y)) ? (x) : (y))
#define MIN_CHILD(x,y) (((x)<(y)) ? (0) : (1))
#define ORI_ABS(x) ((abs((x) - (ORI))))
#define ORI_REV(x) ((GENOME) - (ORI_ABS((x))))
#define ORI_DIST(x) ((MIN_TIME((ORI_ABS((x))),(ORI_REV((x))))))
#define ORI_MIN(x,y) ((MIN_TIME((ORI_DIST(x)),(ORI_DIST((y))))))
#define MAX_VAL(x,y) (((x)>(y)) ? (x) : (y))
#define PCONVERT(x) (((x)+(1))/(2))

#define PROMO_OFF_RXN(x) ((x) + (50))						// Macros for picking out the proper
#define RNAP_LOAD_RXN(x) ((x) + (250))						// "reserved" locations on the reaction
#define RNAP_MOVE_RXN(x) ((x) + (2500)) 					// list; ribosome locations also 
#define INIT_DEGR_RXN(x) ((x) + (5000)) 					// account for multiple genes on unit;
#define CONT_DEGR_RXN(x) ((x) + (15000))					// locations to be changed/added
#define RIBO_LOAD_RXN(x,y) ((((y)-(1))*(10000))+(25000)+(x))			// as necessary
#define RIBO_MOVE_RXN(x) ((x) + (200000)) 
#define CLOSED_TO_OPEN_RXN(x) ((x) + (220000))
#define UNBIND_RNAP_RXN(x) ((x) + (222500))
#define CSAT_TERM_RXN(x) ((x) + (225000))
#define ADD_SEQ_RXN(x) ((x) + (360000))
#define ADD_FISH_RXN(x) ((x) + (420000))
#define SNAPSHOT_RXN(x) ((x) + (245200))
#define UPDATE_RXN(x) ((x) + (245300))
#define ADD_PROMO_RXN(x) ((x) + (249000))

#define RIBO_DECODE_RXN(x) ((x) + (250000))
#define RNAP_ON_F_RXN(x) ((x) + (270000))
#define RNAP_ON_B_RXN(x) ((x) + (272500))
#define RNAP_ADD_RXN(x) ((x) + (275000))
#define RNAP_PAUSE_RXN(x) ((x) + (277500))
#define RNAP_EXIT_RXN(x) ((x) + (285000))
#define RNAP_OFF_F_RXN(x) ((x) + (292500))
#define RNAP_OFF_B_RXN(x) ((x) + (295000))

#define PROMO_KILL_RXN 299000

#define TOPO_I_BINDING_RXN(x) ((x) + (300000))
#define GYRASE_BINDING_RXN(x) ((x) + (305000))

#define NASCENT_PROTEIN_COUNT_RXN(x) ((x) + (319000))

#define TOPO_CLEAVAGE_RXN(x) ((x) + (320000))
#define TOPO_SWITCH_STATE_RXN(x) ((x) + (325000))
#define TOPO_UNBINDING_RXN(x) ((x) + (330000))

#define RNAP_TERM_RXN(x) ((x) + (335000))

#define RNAP_SURVEILLANCE_RXN(x) ((x) + (337500))
#define RNAP_EVAL_RXN(x) ((x) + (340000))

#define RNAP_RELAX_RXN(x) ((x) + (350000))

#define KYMO_RXN(x) ((x) + (480000))
										// Reaction stuctures define:
#define PROMOTER 1								//   Participating objects--
#define RNA 2									//   each object is indexed
#define RNAP 3									//   within its class
#define RIBOSOME 4
#define TOPOISOMERASE 5
										//   And reaction identifiers:
#define PROMO_ON 1 								//     [1] switch promoter on
#define PROMO_OFF 2								//     [2] switch promoter off
#define LOAD_RNAP 3								//     [3] load RNAP
#define LOAD_RIBO 4								//     [4] load ribosome
#define MOVE_RNAP 5								//     [5] advance RNAP 1 nt
#define ABORT_RNAP 6								//     [6] abort tx (not implemented)
#define START_RNAP_PAUSE 7							//     [7] add extra tx pause (not impl)
#define END_RNAP_PAUSE 8							//     [8] end extra tx pause (not impl)
#define MOVE_RIBO 9								//     [9] advance ribosome 3 nt
#define ABORT_RIBO 10								//     [10] abort tsl (not implemented)
#define START_RIBO_PAUSE 11							//     [11] add extra tsl pause (not impl)
#define END_RIBO_PAUSE 12							//     [12] end extra tsl pause (not impl)
#define START_DEGRADE 13							//     [13] start DNA degradation
#define CONT_DEGRADE 14								//     [14] degrade RNA 1 nt
#define COLLISION_ABORT 15							//     [15] abort tsl in CSAT
#define ACTIVATE_RNAP 16							//     [16] convert RNAP->RNAP(open)
#define REMOVE_RNAP 17								//     [17] Remove promoter-bound RNAP
#define RNAP_F_TRANSL 18							//     [18] On-pathway forward transloc
#define RNAP_B_TRANSL 19							//     [19] On-pathway reverse transloc
#define RNAP_NT_ADD 20								//     [20] Add NT to mRNA
#define RNAP_ENTER_PAUSE 21							//     [21] Entry into P1 pause state
#define RNAP_EXIT_PAUSE 22							//     [22] Departure from P1
#define RNAP_TO_P2 23								//     [23] Entry into P2 pause state
#define RNAP_FROM_P2 24								//     [24] Departure from P2
#define RNAP_TO_P3 25								//     [25] Entry into P3 pause state
#define RNAP_FROM_P3 26								//     [26] Departure from P3
#define RNAP_OFF_FOR 27								//     [27] Off-pathway for transloc
#define RNAP_OFF_REV 28								//     [28] Off-pathway rev transloc
#define RIBO_TRANSLOC 29							//     [29] Translocate ribosome
#define RIBO_DECODE 30								//     [30] Decode and extend polypep

										//   obj_index: ID for finding w/in obj.
										//   Pseudo-reactions for scheduling:
#define ADD_SEQ_DATA 31								//     [31] Collection of seq data
#define ADD_FISH_DATA 32							//     [32] Collection of FISH data
#define TAKE_SNAPSHOT 33							//     [33] Position snapshot
#define UPDATE_CONCENTRATIONS 34						//     [34] Adjustments to propensities
#define ADD_PROMOTER 35								//     [35] Adds promoter on replication

#define KILL_PROMOTERS 36							//     [36] Kills promoter at specified
										//	    time

#define TOPO_I_BINDING 37
#define GYRASE_BINDING 38

#define LOG_NASCENT_PROTEINS 44

#define TOPOISOMERASE_CLEAVAGE 46
#define TOPOISOMERASE_STATE_CHANGE 47
#define TOPOISOMERASE_UNBINDING 48

#define RNAP_TERMINATION 49

#define RNAP_SURVEILLANCE 50
#define RNAP_EVALUATION 51

#define RNAP_RELAX 52

#define KYMOGRAPH 53

#endif
