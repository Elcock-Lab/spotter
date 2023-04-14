PROGRAM	 = spotter
PRESRC	 = spotter.master io generate_trajectory reaction_manager \
	   transcription translation supercoiling replication \
	   regulation topo utilities
INCLPRE	 = gen_def sim_types reaction_macros reaction_manager \
	   io reaction_manager regulation replication \
	   supercoiling topo traj transcription translation \
	   utilities

RNG     ?= GSL_MT

ifeq (${RNG},C_RNG)
          SRCEND  = .C_RNG.c
          INCLEND = ${SRCEND:.c=.h}
          IPATH   =
          LPATH   =
          GLIB1   =
          GLIB2   =
endif

SRCEND  ?= .c
SRC      = ${patsubst %,%${SRCEND},${PRESRC}}
SRCDIR   = SRC

OVERDIR  = OBJ
OBJPRE	 = ${SRC:.c=.o}
OBJDIR	 = ${patsubst %,%/SIMRUNNER,${OVERDIR}}
OBJ	 = ${patsubst %,${OBJDIR}/%,${OBJPRE}}

INCLEND ?= .h
INCLDIR	 = INCL
INCL	 = ${patsubst %,${INCLDIR}/%${INCLEND},${INCLPRE}}


DIRPRE   = ${OBJDIR:SIMRUNNER=PRESIM}

PROGTX   = rnap_rate_file_generator
SRCTX    = sequence_utilities.c read_utilities.c RNAP_rate_generator.c \
           hybridization_energies.seq.c NET_seq_dwell_assign.c \
           calculate_pause_entry_rates.c
OBJTX_P  = ${SRCTX:.c=.o}
OBJTX    = ${patsubst %,${DIRPRE}/%,${OBJTX_P}}
INCPRETX = sequence_utilities.h read_utilities.h \
           RNAP_rate_generator.h
INCLTX   = ${patsubst %,${INCLDIR}/%,${INCPRETX}}

PROGTSL  = ribosome_rate_file_generator
SRCTSL   = sequence_utilities.c read_utilities.c ribosome_rate_generator.c \
           decode_based_dwell.c codon_info.c
OBJTSL_P = ${SRCTSL:.c=.o}
OBJTSL   = ${patsubst %,${DIRPRE}/%,${OBJTSL_P}}
INPRETSL = sequence_utilities.h read_utilities.h ribosome_rate_generator.h
INCLTSL  = ${patsubst %,${INCLDIR}/%,${INPRETSL}}


DIRPOST  = ${OBJDIR:SIMRUNNER=POSTSIM}

PROGMV   = trajectory_movie_maker
SRCMV    = movie_maker_macrohelix.c vector_utilities.c \
           xdrfile.c xdrfile_xtc.c
OBJMV_P  = ${SRCMV:.c=.o}
OBJMV    = ${patsubst %,${DIRPOST}/%,${OBJMV_P}}
INPREMV  = vector_utilities.h tagalong_tcl_template.h \
           xdrfile.h xdrfile_xtc.h
INCLMV   = ${patsubst %,${INCLDIR}/%,${INPREMV}}

PROGZMV  = zoomin_movie_maker
SRCZMV   = zoomin_trajectory_maker.c zoomin_structures.c \
           rnap_supercoiling_traj_for_tcl.c vector_utilities.c\
           vector_utilities.c xdrfile.c xdrfile_xtc.c
OBJZMV_P = ${SRCZMV:.c=.o}
OBJZMV   = ${patsubst %,${DIRPOST}/%,${OBJZMV_P}}
INPREZMV = vector_utilities.h zoomin_trajectory.h \
           zoomin_tagalong_tcl.h
INCLZMV  = ${patsubst %,${INCLDIR}/%,${INPREZMV}}

PROGKYMO = kymograph_maker
PREKYMO  = kymograph_maker.c
SRCKYMO  = ${patsubst %,${SRCDIR}/%,${PREKYMO}}

PROGRNA  = rna_trajectory_plotter
PRERNA   = rna_trajectory_plotter.c
SRCRNA   = ${patsubst %,${SRCDIR}/%,${PRERNA}}

CC	 = gcc
GFLAGS 	 = -g
OFLAGS	 = -O3
WFLAG1	 = -Wall
WFLAG2	 =
WFLAGS	 = ${WFLAG1} ${WFLAG2}

CFLAGS	 = ${WFLAGS} ${GFLAGS} ${OFLAGS}

MLIB	 = -lm

#IPATH	= -I/usr/local/include							# Default GSL incl path: if different,
										# change and uncomment or give on 
										# command line

#LPATH	= -L/usr/local/lib							# Default GSL lib path: if different
										# change and uncomment or give on
										# command line

ifneq (${RNG},C_RNG)
 ifneq (presim,${findstring presim,${MAKECMDGOALS}})
  ifneq (postsim,${findstring postsim,${MAKECMDGOALS}})
   ifneq (clean,${findstring clean,${MAKECMDGOALS}})

    ${info COMPILE REQUIRES GSL...}

    ifdef IPATH
     ${info ADDITONAL GSL INCLUDE PATH PROVIDED...}
     ${info ...WILL USE ${IPATH}...}
     ${info ...AND WILL NOT SEARCH FOR ADDITIONAL PATHS}
     ISTAT = GIVEN
    endif
    ifndef IPATH
     ${info LOOKING FOR GSL INCLUDE PATH...}

     LOCINCL  = ${shell locate gsl_rng.h}

     #LOCINCL   = /usr/include/gsl/gsl_rng.h

     GINCFNAL  = ${patsubst %/gsl/gsl_rng.h,%,${LOCINCL}}

     UCHECKINC = ${findstring usr,${GINCFNAL}}

     ifeq (,${UCHECKINC})
        ifneq (,${GINCFNAL})
          IPATH ?= ${patsubst %,-I%, ${GINCFNAL}}
        endif
     endif
    endif

   ifdef LPATH
     ${info ADDITONAL GSL LIBRARY PATH PROVIDED...}
     ${info ...WILL USE ${LPATH}...}
     ${info ...AND WILL NOT SEARCH FOR ADDITIONAL PATHS}
     LSTAT = GIVEN
    endif
   ifndef LPATH
     ${info LOOKING FOR GSL LIBRARY PATH...}
     LOCLIB   = ${shell locate libgsl.so}

     #LOCLIB    = /home/LAB/BIN/ccp4-7.1/lib/libgsl.so
     #LOCLIB    = /home/LAB/BIN/ccp4-7.1/lib/libgsl.so /usr/lib64/libgsl.so

     FINDER    = ${wildcard ${FINDCHK}.*}
     PRECHECK  = ${foreach FINDCHK, ${LOCLIB},${FINDER}}
     CHECKLIB  = ${filter-out ${PRECHECK}, ${LOCLIB}}

     ifneq (libgsl.so.,${findstring libgsl.so.,${CHECKLIB}})
      GLIBFNAL  = ${patsubst %/libgsl.so,%,${CHECKLIB}}
     endif

     UCHECKLIB = ${findstring usr,${GLIBFNAL}}

     ifeq (,${UCHECKLIB})
       ifneq (,${GLIBFNAL})
          LPATH ?= ${patsubst %,-L%, ${GLIBFNAL}}
       endif
     endif
    endif

   endif
  endif
 endif
endif


IPATH   ?= 

CURRDIR  = ${shell pwd}
IPATH2   = -I${CURRDIR}

LPATH   ?=
GLIB1   ?= -lgsl
GLIB2   ?= -lgslcblas

LDLIB 	 = ${MLIB} ${LPATH} ${GLIB1} ${GLIB2}

all: ${PROGRAM} ${PROGTX} ${PROGTSL} ${PROGMV} ${PROGZMV} ${PROGKYMO} ${PROGRNA}

simulation_runner: ${PROGRAM}

presim: ${PROGTX} ${PROGTSL}

postsim: ${PROGMV} ${PROGZMV} ${PROGKYMO} ${PROGRNA}


gslstat:
ifneq (${RNG},C_RNG)
ifneq (GIVEN,${findstring GIVEN,${ISTAT}})
	@echo FOUND THESE LOCATIONS FOR GSL INCLUDE: ${GINCFNAL}
ifneq (,${UCHECKINC})
	@echo ${UCHECKINC} ON STANDARD SEARCH PATH
	@echo NO ADDITONAL PATHS ADDED
endif
ifeq (,${UCHECKINC})
	@echo COULD NOT FIND GSL INCLUDE ON STANDARD SEARCH PATH...
ifneq (,${GINCFNAL})
	@echo ...BUT DID FIND IT IN THESE LOCATIONS: ${GINCFNAL}
endif
ifeq (,${GINCFNAL})
	@echo ...OR ANY ANOTHER LOCATION...WILL TRY COMPILING ANYWAY...
endif
endif
endif
ifneq (GIVEN,${findstring GIVEN,${LSTAT}})
	@echo FOUND THESE LOCATIONS FOR GSL LIBRARY: ${GLIBFNAL}
ifneq (,${UCHECKLIB})
	@echo ${UCHECKLIB} ON STANDARD SEARCH PATH
	@echo NO ADDITONAL PATHS ADDED
endif
ifeq (,${UCHECKLIB})
	@echo COULD NOT FIND GSL LIBRBARY ON STANDARD SEARCH PATH...
ifneq (,${GLIBFNAL})
	@echo ...BUT DID FIND IT IN THESE LOCATIONS: ${GLIBFNAL}
endif
ifeq (,${GLIBFNAL})
	@echo ...OR ANY ANOTHER LOCATION...WILL TRY COMPILING ANYWAY...
endif
endif
endif
endif

${PROGRAM} : ${OBJ}
	${CC}  -o $@ $^ ${CFLAGS} ${LDLIB}

${OBJDIR}/%.o: ${SRCDIR}/%.c gslstat ${INCL} 
	${CC} ${IPATH} ${IPATH2} -c -o $@ $< ${CFLAGS} ${LDLIB}


${PROGTX} : ${OBJTX}
	${CC}  -o $@ $^ ${CFLAGS} ${MLIB}

${PROGTSL} : ${OBJTSL}
	${CC}  -o $@ $^ ${CFLAGS} ${MLIB}

${DIRPRE}/%.o : ${SRCDIR}/%.c ${INCLTX} ${INCLTSL}
	${CC} ${IPATH} ${IPATH2} -c -o $@ $< ${CFLAGS} ${MLIB}


${PROGMV} : ${OBJMV}
	${CC}  -o $@ $^ ${CFLAGS} ${MLIB}

${PROGZMV} : ${OBJZMV}
	${CC}  -o $@ $^ ${CFLAGS} ${MLIB}

${PROGKYMO} : ${SRCKYMO}
	${CC} ${SRCKYMO} -o ${PROGKYMO} ${CFLAGS} ${MLIB} 

${PROGRNA} : ${SRCRNA}
	${CC} ${SRCRNA} -o ${PROGRNA} ${CFLAGS} ${MLIB} 

${DIRPOST}/%.o : ${SRCDIR}/%.c ${INCLMV} ${INCLZMV}
	${CC} ${IPATH} ${IPATH2} -c -o $@ $< ${CFLAGS} ${MLIB}


${OBJDIR}: | ${OVERDIR}
${DIRPRE}: | ${OVERDIR}
${DIRPOST}: | ${OVERDIR}

${OVERDIR} :
	mkdir ${OVERDIR}

${OBJ}: | ${OBJDIR}

${OBJDIR} :
	mkdir ${OBJDIR}

${OBJTX}: | ${DIRPRE}
${OBJTSL}: | ${DIRPRE}

${DIRPRE} :
	mkdir ${DIRPRE}

${OBJMV}: | ${DIRPOST}
${OBJZMV}: | ${DIRPOST}

${DIRPOST} :
	mkdir ${DIRPOST}

ifeq (${RNG},C_RNG)
            ADDOBJ    = ${OBJ:.C_RNG.o=.o}
endif

ADDOBJ   ?= ${OBJ:.o=.C_RNG.o}
CLEANOUT  = ${PROGRAM} ${PROGTX} ${PROGTSL} ${PROGMV} ${PROGZMV} \
            ${PROGKYMO} ${OBJ} ${OBJTX} ${OBJTSL} ${OBJMV} ${OBJZMV} \
            ${ADDOBJ} ${PROGRNA}

.PHONY: clean

clean:
	rm -f ${CLEANOUT}

