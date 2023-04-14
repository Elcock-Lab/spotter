<img align="right" width="500" height="302" src="https://github.com/Elcock-Lab/spotter/blob/main/IMG/220526_spotter_overview.png">

# _spotter_
<br/>

_spotter_ (**s**imulator of **p**rokaryotic **o**peron **t**ranscriptional and **t**ranslational **e**longation **r**eactions) is an integrated simulation model offering highly-detailed representations of prokaryotic transcription, translation, and DNA supercoiling. _spotter_ is designed to incorporate sequencing data and allows users to explore the interplay of transcription and translation.

<br/>

## CITING _spotter_

If you use _spotter_ in your work, please cite:

Hacker, W.C., Elcock, A.H. (2023) _spotter_: A single-nucleotide resolution stochastic simulation model of supercoiling-mediated transcription and translation in prokaryotes. (preprint) bioRxiv.

<br/>

## GETTING STARTED: A TEST RUN WITH THE _E. coli alaS_ OPERON

The test run below uses the _alaS_ sequence (provided in the EXAMPLES folder), but if a sequence
for a particular transcription unit is not available, simulations can proceed without one; as
described below, users have a variety of options for providing sequences, dwell time data, and
parameters defining alternative models of transcription and translation.

To set up the run, leave the spotter directory and create a directory for the trial run...

`mkdir ALA_S_TEST`

...change into that directory, and copy the provided _alaS_ sequence and two example files into it:

`cp /home/spotter/EXAMPLES/STARTUP/alaS.seq .`
	
`cp /home/spotter/EXAMPLES/STARTUP/consolidated.alaS.inp .`

`cp /home/spotter/EXAMPLES/STARTUP/alaS.gene .`

(note that `/home/spotter/` should be replaced by the location of the spotter folder on your system)

Once these files are in place, you can start making the files needed to run a simulation.

<br/>

**First, use the sequence file to generate a position-specific list of reaction rates for RNAPs at each DNA template position:**

  - `/home/spotter/PRECOMPILED/rnap_rate_file_generator  -sequence alaS.seq -tx_start 13 -tx_stop 2769`

	Several options are available for calculating RNAP reaction rates (these are detailed in [Making position-specific RNAP reaction rate files](#making-position-specific-rnap-reaction-rate-files) below). Here, we use default settings where dwell times are calculated using a NET-seq-derived energy function and RNAPs travel at a mean rate of 30 bp/s and are subject to short elemental pauses but do not enter longer-lived pauses or backtrack.

<br/>

**Rename the output RNAP rate file:**

  - `mv full_rnap_rate_set.seq_based rnap_rates.alaS`

<br/>

**Next, use the aligned sequence file made in generating RNAP rates to make a position-specific list of reaction rates for ribosomes at each mRNA template position:**

 - `/home/spotter/PRECOMPILED/ribosome_rate_file_generator  -sequence aligned_sequence.txt  -genefile alaS.gene`
 
 	Like RNAP reaction rates, ribosome rates can be calculated using a variety of options (detailed in [Making-position-specific ribosome reaction rate files](#making-position-specific-ribosome-reaction-rate-files) below). Here, we use default settings where ribosomes are assumed to move at a mean rate of 10 codons/s and their dwell times are determined by a tRNA competition model.
 
 <br/>
 
 **Rename the output ribosome rate file...**
 
  - `mv ribo_rates.seq.full_info ribo_rates.alaS`

<br/>

**...then modify the simulation input file (`consolidated.alaS.inp`) so that it lists the locations of the position-specific rate files just made:**

Change the lines to match the following:

 - <code>rnap_dwell_file&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;rnap_rates.alaS</code>

	(located under "Position-dependent reaction rates" in the "Transcription" section in the simulation input file)

 - <code>ribo_dwell_file&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;ribo_rates.alaS</code>


	(located under "Position-dependent reaction rates" in the "Translation" section in the simulation input file)

	The simulation input file used is "annotated." Annotations describe options available for general simulation parameters, transcription and supercoiling models, ribosome-RNAP and inter-RNAP interactions, and rates for position-independent reaction rates. Here, the input file sets up a 15-minute simulation with a simple model of RNAP activation, a single-step model of transcription (no elemental pauses, advanced pauses or backtracking), a ribosome initiaion rate of 0.12 s<sup>-1</sup>, and a promoter that toggles between "on" and "off" states. Many of these parameters can be adjusted without needing to remake the position-dependent rate files; note, though, that if you change position-**independent** rates for RNAP pause entry and escape in the simulation input file, you will need to remake position-**dependent** RNAP rate files to match these changes.

<br/>

**Now you're ready to run a simulation:**

 - `/home/spotter/PRECOMPILED/spotter.startup  -singlefile consolidated.alaS.inp -seed 17123 > simoutput`
 
 	**NOTE:** The above uses a version of the simulation runner compiled using a static version of the GSL library and so has a reasonable chance of running out of the box. If it's not working, you can also try a version compiled with the built-in C random number generator: `/home/spotter/PRECOMPILED/spotter.C_RNG`. See [Compiling _spotter_](#compiling-spotter) for further information.
 
 	There are several options available in running simulations (see [Running _spotter_ simulations](#running-spotter-simulations) for details). Input files can be given as modules (with separate files for transcription, translation, etc.) but are here collected in a single file.

	Simulation output files are described in [Running _spotter_ simulations: Simulation output](#simulation-output). Files include logs of transcriptional and translational activity and data for making plots and movies representing the trajectory.

<br/>

**If you have Visual Molecular Dynamics (VMD) installed on your desktop, you can now make a movie of the simulation trajectory. First, make VMD-readable files representing the trajectory:**

 - `/home/spotter/PRECOMPILED/trajectory_movie_maker  tx_tsl_traj.for_xtc.txt  consolidated.alaS.inp  200  50  0.1  alaS`

<br/>

**Next, copy the directory made by the movie-maker (`MOVIE_MATERIALS.alaS`) to your desktop. Start VMD by opening the `.pdb` file in this directory (if your version autoloads a ruler or other "molecules," don't open VMD first). In the Tk Console (available from the dropdown "Extensions" menu), enter...**

 - `source custom_script_for_alaS.tcl`
 
 <br/>
 
 **...followed by:**
 
  - `start_tx_tsl_viewer 1 white`

	Scroll in and out briefly to center text. The movie can played with "play" button in the lower right-hand corner of the display window. Because the fractional availability of the promoter is 1/3, RNAPs may not initiate at the beginning of the simulation; you can try other random seeds for comparison. See [Post-simulation visualization: Making movies from trajectories](#post-simulation-visualization-making-movies-from-trajectories) for additional details.

<br/>

## COMPILING _spotter_


All of the programs needed for running _spotter_ can be compiled using the
makefile included in this directory.

*REQUIRED:* gcc compiler (standard with Linux)

*HIGHLY RECOMMENDED:* GNU Scientific Library (GSL)

 ### Compiling with GSL
-------------------

If you have access to the GNU Scientific Library on your machine, you can
compile a version of the simulation code that uses the GSL implementation
of the Mersenne Twister random number generator.

To compile using the GSL RNG, enter (in this directory):

	make
or:

	make all

The makefile in this directory first looks for GSL on the machine it's being run
on and reports what it finds--so even if you're not sure about the existence
and/or location of GSL it's worth trying with the basic command to start.

 - If the GSL files are in standard locations (`/usr/local/include` or `/usr/include`
for the include path and `/usr/local/lib` or `/usr/lib` for the library linkage),
you should be good to go: the compiler will be able to include/link to files in
these locations.

 - If GSL library or header files exist but are in nonstandard locations, there's a
good chance you'll still be able to compile with the basic command and an
unaltered makefile: if found, paths to these locations will be passed
automtically to the compiler.

 - If running the command above fails to locate the correct filepaths but you know
where the files are located (or if you have a standard GSL installation but 
would prefer to use a version in a different location), enter:

	`make IPATH=-I/XXX/XXX LPATH=-L/YYY/YYY`

to have the compiler use one or both of the alternative include (`IPATH`) or
library (`LPATH`) locations: `/XXX/XXX` and `YYY/YYY` are the filepaths for these
locations, respectively (note that you need to keep the `-I/` and `-L/` prefixes).

Alternatively, you can modify the makefile and direct the compiler to your
GSL include and library locations by uncommenting and modifying lines 92 and 95.

--Note that if you are running on a cluster where GCC and GSL are loaded as
modules into your environment, the make command may not be able to locate the
appropriate filepaths but will still compile appropriately. JUST BE SURE TO LOAD
THESE MODULES (BOTH GCC and GSL, using "module load" or whatever is appropriate
to your system), or the program will not compile.

 ***If you don't have GSL on your system but want to run with its RNG:***

It's definitely worth trying the command above in case there's a rogue version
of GSL available on your system. If it's clear that GSL isn't available on
your system, you can install it.

You can download GSL following instructions at:

https://www.gnu.org/software/gsl/

With root access you can install to `/usr/` so that files are in standard
locations. Without root access you can install to an alternative directory
(instructions at https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/).
`-IPATH` and `-LPATH` should be modified to match the nonstandard location.

**COMPILATION OUTPUT:**

The compilation will produce executables for preparing simulation input files,
running simulations, and visualizing simulation output.

The files created are:

 - `rnap_rate_file_generator`

	Used to make required DNA-position-specific RNAP reaction rates. Discussed in [Preparing for simulations: Making position-specific RNAP reaction rate files](#making-position-specific-rnap-reaction-rate-files).

 - `ribosome_rate_file_generator`

	Used to make required mRNA-position-specific ribosome reaction rates. Discussed in [Preparing for simulations: Making position-specific ribosome reaction rate files](#making-position-specific-ribosome-reaction-rate-files).

 - `spotter`
 
 	Used to run simulations of transcription and translation. Discussed in [Running _spotter_ simulations](#running-spotter-simulations).

 - `trajectory_movie_maker`
 
 	Used to make operon-level movies of simulation trajectories. Discussed in [Post-simulation visualization: Making movies from trajectories](#post-simulation-visualization-making-movies-from-trajectories).

 - `zoomin_movie_maker`
 
 	Used to make transcription-only movies of a portion of the DNA template to illustrate RNAP rotation and local supercoiling. Discussed in [Post-simulation visualization: Making "zoomed-in" movies](#post-simulation-visualization-making-zoomed-in-movies).

 - `kymograph_maker`
 
 	Used for making 2-D plots of changing supercoiling density and RNAP positions as a function of time. Discussed in [Post simulation visualization: Making kymographs](#post-simulation-visualization-making-kymographs).

- `rna_trajectory_plotter`

	Used for making 3-D plots of production, degradation, and translation of individual RNAs made during the simulation. Discussed in [Post simulation visualization: Making RNA plots](#post-simulation-visualization-making-plots-of-simulation-rnas)

**NOTES ON COMPILING/RUNNING WITH THE GSL RNG:**

--If you're compiling with an older version of GCC, the compiler may issue
  warnings about missing braces around initializers--these are not a problem

--If you run simulation code in locations different from the location where it
  was compiled (e.g., compute nodes when the program was compiled on a head
  node), the *SAME* GSL library used in the compilation should be available
  in the location where the simulation is being run. The program *may* run with
  a version of GSL different from that with which it was compiled, but if you're
  on a system where different versions of GSL are available you'll want to keep
  track of the version you compiled with and use that in subsequent runs.


 ### Compiling WITHOUT GSL
----------------------

If you do not have access to the GNU Scientific Library, you can compile a
version of the simulation code that uses the built-in C random generator.

NB: THE BUILT-IN RNG IS NOT A HIGH-QUALITY RANDOM NUMBER GENERATOR! While
this may be a great way to get up and running quickly and simulation results
will be qualitatively very similar to those produced with a simulation-quality
RNG, you will want to install/recompile with GSL if you're going to be working
closely with simulation data.

To compile using the built-in RNG, enter (in the folder with the makefile):

	make RNG=C_RNG

or:

	make all RNG=C_RNG

Running with C_RNG option will produce the same files as are created in the
GSL-based compilation. Except for the simulation runner (spotter) all other
files (which do not depend on the RNG) are identical.



**Cleaning the directory**


To remove old executables and object files, enter (in the makefile directory):

	make clean

Note: this will remove ALL object files (from both GSL and C_RNG compilations,
if for some reason you've compiled both ways).


**Other targets**


If you only need to update the simulation code proper, you can use:

	make simulation_runner

Similarly, if you're only updating the reaction rate generators:

	make presim

And if you're only changing post-simulation code:

	make postsim

The last two options bypass the search for GSL and so can be faster on machines
where this search takes a while and  the simulation runner doesn't need to be
recompiled.


### Pre-compiled binaries
----------------------

Precompiled versions of the simulation code are located in the "PRECOMPILED"
folder in this directory. Note that the GSL-based version of _spotter_ (compiled
with a static link to GSL) may not work with an arbitrary set-up and the
program will need to be recompiled on your system.

The following have a good chance of being functional on an arbitrary
system:

	rnap_rate_file_generator

	ribosome_rate_file_generator

	spotter.C_builtin_RNG

	trajectory_movie_maker

	zoomin_movie_maker

	kymograph_maker
	
	rna_trajectory_plotter

There are three pre-compiled GSL-dependent versions of the simulation runner:

	spotter.startup

	spotter.static_build2

	spotter.dynamic_build

The first version was compiled with static linkage to the GSL library
(GCC 4.4.7 using archive files from GSL 1.15 on a machine using GLIBC 2.12 and
CentOS 6.9). The binary worked on all of the machines (older and newer) on 
which we attempted to run it, but the test set was limited and the binary may
not work in all environments.

The second version (spotter.static_build2) was compiled with the same GSL
archive files but with GCC 5.4 on a machine using GLIBC 2.17 and the CentOS
7.9.2009 operating system. It may work on machines with newer/updated C
libraries but will not work on older machines on which it does not have access
to these libraries.

The third version (spotter.dynamic_build) was compiled with dynamic linkage to
shared GSL library files (version 2.3) using GCC 5.4 on a machine with the
CentOS 7.9.2009 operating system. There's a good chance it won't work with
many setups (it will definitely not work if GSL has not been installed), but
it may be worth trying on a newer machine.

This folder also includes a GSL-independent version of the simulation runner:

	spotter.C_builtin_RNG

This version of the simulation runner is the most likely to run without
recompilation but uses a lower-quality RNG (see [Compiling without GSL](#compiling-without-gsl) above)--if you
have access to GSL, you'll want to recompile in order to take advantage of
the much better RNG it offers.

All pre- and post-simulation programs were compiled on a machine with the
CentOS 7.9.2009 OS using GCC 5.4.
<br/>
<br/>
<br/>


  
## PREPARING FOR SIMULATIONS
 
 


All simulations require a set of input files to run: a file listing position-specific
rates for reactions affecting RNAPs (described in the next section), a file listing
position-specific rates for reactions affecting ribosomes (described in [Making
position-specific ribosome reaction rate files](#making-position-specific-ribosome-reaction-rate-files)), and a group of files defining basic
simulation parameters and the rates of position-independent reactions (described in
[Preparing additional simulation input files](#preparing-additional-simulation-input-files)).

### Making position-specific RNAP reaction rate files

Simulations require a file detailing the rates at which RNAPs undergo position-specific
(translocation and pause-entry) reactions at each position along the DNA template.

The required rate file can be made using a single command, which can be run with a variety
of user-selected options. These rates can be calculated from a real transcribed sequence
from any organism provided by the user; alternatively, users can create random or uniform
sequences if they are more interested in general phenomena (e.g., the effect of isolated
pauses or pause clusters). Similarly, users can supply dwell times for their template or
or can choose to have these dwell times calculated for a sequence they provide.

#### Commands and flags
---------------------------------------

**To run the RNAP file-generating program, use:**

	./rnap_rate_file_generator

The options available in the program are accessed with flags following the basic command.
The flags are as follows:


  |Flag|Options|Comments|
  |----|-------|--------|
  |`-sequence`|file name / random / flat|Required. If a file name is given, the file should contain the sequence of interest in FASTA format (initial comment line OK). If no sequence is provided, a script will generate either a random sequence (random) or a uniform sequence (flat) in which translocation and pause-entry rates are identical at all sites.|
  |`-tx_start`|integer|Required. Gives the position of the TSS relative to the first NT in the sequence provided. Since reaction rates depend on the upstream sequence, ideally the file provided will include ~15 bp up- and downstream of the transcription unit. If these are not provided, they are filled in randomly as necessary. If a sequence is not provided, the start site position should be given as "1".|
  |`-tx_stop`|integer|Required. Gives the position of the end of the transcription unit in the sequence file provided (see note above on adding bp downstream of the stop). If a sequence is not provided, the stop position should be the length of the transcription unit desired (matching the dwell time file, if provided).|
  |`-reverse`|file name|Optional. If included, the sequence used is the reverse complement of that in the listed file. NOTE: the tx_start and tx_stop values provided should be those of the ORIGINAL, UNREVERSED sequence. These will be automatically changed in reversing the sequence (start and stop can be in either order). If you need a reverse complement, use the same file name for both the the -sequence and the `-reverse` flags.|
  |`-dwell`|file name|Optional. If included, the file listed is used to assign to dwell times at each position. If not included, the NET-seq energy function is used to determine dwell times based on the sequence that has been provided or generated.|
  |`-rate`|floating point number|Optional. If provided, dwell times are normalized so that the mean elongation rate over the template is the value given in nt/s. If a dwell time file is provided, the dwell times in it will NOT be normalized (i.e., they will be absolute dwell times in seconds) if this flag is not used. If no dwell time file is given, NET-seq dwell times are normalized to 30 nt/s when the -rate flag is not used.|
  |`-seed`|integer|Optional. If provided, it is used to seed the random sequence assignment when "random" is selected for -sequence. If a random sequence is requested without a seed, seeding is based on system time.|
  |`-repeat`|A/T/C/G|Optional. If provided, gives the base that will be repeated in making the uniform sequence. Defaults to "A" if the flag is not used.|
  |`-param`|file name|Optional. Gives the name of the file to be used in calculating elemental-pause entry reactions. Defaults to a file in which advanced pause entry rates are set to zero.|
  |`-scrunch`|n/a|Optional. If flag is included, rates of translocation for the intial 9 positions in the template are calculated assuming that DNA is "scrunched" into an RNAP that remains stationary at the promoter. If flag is not included, rates are  calculated as if a complete RNA/DNA hybrid existed at these positions|
  |`-slide`|n/a|Optional. If flag is included, the upstream bubble edge is assumed to move downstream immediately on the forward translocation step, producing a constant bubble size. If flag is not included, the upstream bubble edge advances downstream on nucleotide addition.|
						
<br/>

#### Examples
---------






 - **Example of a run where a user is interested in simulating the _E. coli gapA_ gene and wants to use the real sequence (padded by 15 bp at either end) and NET-seq derived dwell times in a simulation where transcription occurs at 50 nt/s over the gene:**

  	`./rnap_rate_file_generator -sequence gapA.seq  -tx_start 16 -tx_stop 1011 -rate 50.0`


 - **Example of a run where a user has a sequence and a dwell time file:**

  	`./rnap_rate_file_generator -sequence gapA.seq  -tx_start 16 -tx_stop 1011 -dwell gapA.dwell.txt`

 - **Example of a run where a user has a dwell time but is uninterested in sequence:**

  	`./rnap_rate_file_generator -sequence flat  -tx_start 1  -tx_stop 500  -repeat T  -dwell dwellfile.txt`

<br/>

--------------------------------------------------------------
#### Notes
--------------------------------------------------------------

Some details on files supplied to the script:

 - **SEQUENCE FILES**

All FASTA files (with or without header lines) will work. Sequences should be at least
as long as the transcription unit of interest, but ideally they will also contain 12-15
nt of "padding" upstream of the TSS and downstream of the transcriptional stop site. The
position of the transcription unit proper can then be designated with the `-tx_start` and
`-tx_stop` flags described above. In the absence of these flanking regions, sequences will
be filled out randomly; rates will be affected only in the first 11 and last two positions
on the template.

By default, it is assumed that the sequence provided is that of the coding strand. If
necessary, the sequence of the template strand can be provided instead and treated with
the `-reverse` flag above. As noted above, start and stop positions must be those of the
ORIGINAL, UNREVERSED strand.

 - **DWELL TIME FILES**

Dwell time files should list in each line the template position and (separated by one or
 more spaces) the dwell time. If the `-rate` flag is not used for normalization, dwell
times listed in the file are assumed to be ABSOLUTE (in seconds). If normalization
is used, dwell times are treated as relative and are adjusted to match the selected
elongation rate.

The number of entries in the dwell time file MUST MATCH the length of the transcription
unit as given by the start and stop flags (since the sequence can be longer than the
unit proper [see above] these are the key values defining its length). If no sequence
is provided (i.e., `-sequence` is set to "random" or "flat"), `-tx_start` and `-tx_stop`
ahould be set to 1 and N, respectively, where N is the number of entries in the
dwell time file.

 - **PARAMETER FILES**

These contain rates for position-independent reactions including advanced pause entry,
pause escape, nucleotide association, and nucleotide addition. All values listed in
the given files are required; the effect of different rates can be explored by 
modifying values in the template RNAP_invariant_rates.txt (in the "INPUT_FILES" 
directory) and including the modified file with the `-param` flag.




**FILE CHECKLIST:**

For random or flat sequences, the file generator is self-contained; if rate files
are to be derived from DNA sequence and/or lists of dwell times, these need to be
supplied, as should parameters for position-invariant reactions if users do not wish 
to use default rates:

	-sequence file (optional)
	-dwell time file (optional)
	-position-invariant reaction parameters (optional)
<br/>



### Making position-specific ribosome reaction rate files


Simulations require a file detailing the rates at which ribosomes undergo addition and
translocation reactions at each position along mRNA.

NOTE: THIS FILE NEEDS TO BE MADE FOR ALL SIMULATIONS, INCLUDING THOSE
THAT DO NOT MODEL TRANSLATION. A simple dummy file can be made for systems without
translation using commands below.

Like the RNAP rate file, the ribosome rate file can be made with a single command.
It is most straightforward to generate the ribosome file AFTER making the RNAP rate
file--the sequence file made in making the RNAP file (`aligned_sequence.txt`) can be used
directly.

Details about files required for the program are listed at the bottom of this README;
in every case users will need to supply a sequence and a list of genes. The rate file
can incorporate user-supplied dwell times (derived from ribosomal profiling or custom
made) and/or tRNA competition-driven decoding times. Files for simpler schemes with
uniform addition rates can also be made.

#### Commands and flags
---------------------------------------

**To run the ribosome rate file generator, use:**

	./ribosome_rate_file_generator

As with the RNAP generator, available options are accessed with flags following the basic
command. The flags are as follows:


  |Flag|Options|Comments|
  |----|-------|--------|
  |`-sequence`|file name|Required. The simplest option is to use the "aligned_sequence.txt" file made in running the RNAP file-making script. If a different file is used, users should verify that it starts at the TSS and ends at the transcription start site.|
  |`-genefile`|file name|Required. This is a list of of all genesincluded in the transcription unit (the file format is listed below) with their start position (first NT of the first codon) and stop position (last NT of the last/stop codon) relative to the TSS.|
  |`-dwell`|file name / flat / use_decode|Optional. If a file name is given, that file will contain relative dwell times for each *NT* (not codon) and so match the length of operon. If "flat" is selected, identical relative dwell times will be assigned at all positions. If "use_decode" is selected, dwell times proportional to decoding times for each codon will be assigned; this is the default setting.|
  |`-decode`|file name / flat|Optional. If a file name is given, that file will list decoding times for all codons in *milliseconds*. If the flag is not supplied, the script defaults to decoding times derived from a tRNA competition model for E. coli. If "flat" is selected, the same decoding is used for all codons; the default is 100 ms but changeable with the `-dt` flag.|
  |`-dt`|floating point number|Optional. If the "flat" option has been selected with the `-decode` flag, sets the fixed decoding time for all codons. The value given is in milliseconds.|

<br/>

#### Examples
---------

 - **Example of a run where a user has both a dwell time file and a codon-specific list of decoding times:**

  	`./ribosome_rate_file_generator  -sequence aligned_sequence.txt  -genefile gene.list  -dwell ribo_seq.gapA  -decode codon_add_times`


 -  **Example of a run where a user has a list of decoding time but not a dwell time file--in this case case dwell times will be assigned based on decoding times:**

  	`./ribosome_rate_file_generator  -sequence aligned_sequence.txt  -genefile gene.list -decode codon_add_times`


 - **Example of a run where a user has a dwell time file but not codon-based addition times--in this case decoding times will default to 100 ms at all codons:**

  	`./ribosome_rate_file_generator  -sequence aligned_sequence.txt  -genefile gene.list -dwell ribo_seq.gapA`


 - **Example of a run where a user wants a uniform rate of translation of 12.5 codons/s--the decoding time is set to 80 milliseconds:**

  	`./ribosome_rate_file_generator  -sequence aligned_sequence.txt  -genefile gene.list -dwell use_decode  -decode flat  -dt 80.0`

<br/>

--------------------------------------
#### Notes
--------------------------------------

Some details about files required for the ribosome rate generator:

 - **SEQUENCE FILES**

All FASTA files (with or without header lines) will work. If the `aligned_sequence.txt`
output file from the RNAP file-making script is used, it will automatically be trimmed
and aligned. If any other files are used, they MUST run from the TSS of the transcription
unit of interest to the last transcribed position in the unit and match the length of
any dwell time file supplied. 

 - **DWELL TIME FILES**

Dwell time files should list in each line RNA template position and (separated by one or
more spaces) the dwell time. The template position should be given relative to the TSS,
so that the position 1 in the file correspons to the first NT transcribed in the 
transcription unit. NOTE: The file MUST contain dwell times for all *nucleotides*, not
all codons on the template mRNA, and so match the length of the transcription unit.
Only the *relative* values of the times given matter for simulations: in another input
file used at simulation run-time, the mean rate of elongation on each gene is selected
and the times will be scaled accordingly.

 - **GENE INFORMATION FILE**

This file is required. It should contain one line for each gene included on the
transcription unit of interest. Each line should contain a name or unique identifier for
that gene, the start position (first NT of first codon relative to the TSS), and the
stop position (last NT of last codon relative to the TSS) of the gene, each separated by
one or more spaces. 


 - **DECODING INFORMATION FILES**

The default file used is fluitt_based_transloc_decoding_split.txt. If another file is
supplied, it should follow the same format: each line should contain: 1) the three-
letter codon identifier; 2) the fraction of total time spent in decoding; and 3) the
total decoding time in milliseconds. The order in which the codons appear in the file
will not affect its readability.


**FILE CHECKLIST:**


The following files are required for the rate generator to work:

	sequence file
	gene file

Additional files cen be supplied if desired:

	dwell time file
	decoding time file

 <br/>


### Preparing additional simulation input files
--------------------------------------------

#### Addtional files required to run simulations
------------------------------

In addition to the files made in following these instructions mentioned above, the
simulation runner also accepts (and in some cases requires) additional input files:

 - **a SYSTEM settings input file (REQUIRED)**

	  This file allows users to set basic simulation parameters, including the
	  length of the simulation and the frequency with which sequencing data are
	  collected, and to switch various output options on or off.

	  A template system settings file can be bound in EXAMPLES:

		sys.template.inp

 - **TRANSCRIPTIONAL settings input file (REQUIRED)**


	  This file allows users to set RNAP binding and activation rates, choose among
	  pausing and backtracking models, and define the rates of position-independent
	  transcription reactions.
	  
	  ***-->THIS FILE SHOULD LIST THE POSITION-SPECIFIC RNAP RATE FILE MADE IN THE STEPS ABOVE IN THE LINE SPECIFYING THE***
	  `rnap_dwell_file`
	  ***TO BE USED IN THE SIMULATION***
	  			

	  A template transcriptional input file can be bound in EXAMPLES:

		tx.template.inp

 - **a TRANSLATIONAL settings input file (REQUIRED)**

	  This file allows users to set gene boundaries and mean ribosomal loading and
	  translation rates. It is also used to set RNA lifetimes and the rate at which
	  the degradosome advances along RNAs.

	  ***-->THIS FILE SHOULD LIST THE POSITION-SPECIFIC RIBOSOMAL RATE FILE MADE IN THE STEPS ABOVE IN THE LINE SPECIFYING THE***
	  `ribo_dwell_file`
	  ***TO BE USED IN THE SIMULATION***


	  An example of a transcriptional input file can be bound in EXAMPLES:

		tsl.template.inp

 - **a SUPERCOILING settings input file (OPTIONAL)**

	  This file allows users to define the kind of template (chromosomal or plasmid)
	  on which transcription occurs and to set parameters governing topoisomerase I
	  and gyrase activity and the response of RNAPs to supercoiling.

	  ***-->IF A SUPERCOILING FILE IS NOT PROVIDED, THE SIMULATION IS ASSUMED TO BE A "NO SUPERCOILING" SIMULATION***



	  A template supercoiling input file can be bound in EXAMPLES:

		super.template.inp

 - **a GENE-REGULATION input file (OPTIONAL)**

	  This file allows users to set parameters affecting promoter availability for
	  the gene or operon of interest. If no gene-regulation file is provided the
	  promoter is assumed to remain on throughout the simulation.

	  A template regulation input file can be bound in EXAMPLES:

		reg.template.inp

----
   #### Notes
   --------------------------------

 - Fuller explanations of the settings available in each of the input files
	  listed above can be found in the "annotated" versions of each of these files
	  in the EXAMPLES folder

 - The input files above have broken into categories (system, transcription,
	  etc.) for clarity and to make modular changes to simulation easier. If users
	  want to combine all modules in a single input file, they can do so by cat-ing
	  the individual template files and running with the ` -singlefile ` flag (discussed below in RUNNING SIMULATIONS)

<br/>

## RUNNING _spotter_ SIMULATIONS


### Commands and flags
------------------------

To run simulations, use:

        ./spotter

The options available in the program are accessed with flags following the basic command.
The flags (which can be given in any order) are as follows:


  |Flag|Options/values|Comments|
  |----|--------------|--------|
  |`-sysfile`|system settings input file|Required. Described above.|
  |`-txfile`|transcriptional settings input file|Required. Described above.|
  |`-tslfile`|translational settings input file|Required. Described above.|
  |`-superfile`|supercoiling settings input file|Optional. Described above.|
  |`-regfile`|regulatory settings input file|Optional. Described above.|
  |`-singlefile`|consolidated input file|Optional. Use if you want to consolidate the modular input files above and supply a single input file; the modular files are *not* required if the singlefile option is selected. An example of a consolidated input file can be found in the EXAMPLES folder.|
  |`-seed`|integer|Optional. Sets the seed for the random number generator used in simulations (using either the GSL RNG or the native C RNG). If a seed is not provided, it will be set 1234.|
  |`-traj`|integer|Optional. Numerical label for the simulation attached to output files. Defaults to 1.
  |`-shut`|time (s)|Optional. Time at which promoter is shut off in simulations, preventing further transcription. Defaults to no shutoff.|

<br/>

**CHECKS BEFORE RUNNING SIMULATIONS:**


 Make sure that all required input files have been copied to the location where
    you're running simulations. These include:

 - **All module input files (system, transcriptional, and translational files, plus
	  supercoiling and regulatory files if you'll be using them)--or, if you're using the `-singlefile` option, a consolidated input file including parameters for all of these modules**


 - **The sequence-specific files for RNAP and ribosome reactions made using the
	  rate generators above**

<br/>

### Examples
---------

 - **A run where the user wants to use all of the modular input files made for the alaS gene:**

	`./spotter  -sysfile alaS.sys.inp  -txfile alaS.tx.inp  -tslfile alaS.tsl.inp  -superfile alaS.super.inp  -regfile alaS.reg.inp -seed 2413 > simoutput`

 - **A run where input files for alaS have been consolidated into one file:**

	`./spotter  -singlefile  alaS.consolidated.inp  -seed 11981 > simoutput`


NOTE that in the examples above output is sent to "simoutput"--when simulations are run
in "quiet mode" on-screen output is modest (~1000 lines for a frequently-initiated gene
on a constantly-on promoter) but it's usually best to redirect it to reduce visual
clutter.

<br/>

------------------
### Simulation output
------------------

At the end of a simulation run, a number of files will be produced. These include:

 - **A summary file:**
   

	NAME_summary_NUMBER.txt

		where the NAME is the name assigned to the system in the sys.inp file
		and NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file records parameters used in the simulation, as well as values
characterizing RNA and protein production, translational efficiency and
translational elongation rates


 - **A time course of protein production:**

	protein_log_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file lists the number of proteins present in the simulation at ten-second
intervals


 - **A time course of RNA production:**

	probe_based_RNA_abundance_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file reports the abundance of RNAs in the simulation at ten-second
intervals; abundances are reported for two positions


 - **An RNAP log:**

	tx_rate_info_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file gives information about each of the RNAPs that were activated during
the simulation. Each line represents an RNAP. The first column lists the
elongation rate of the RNAP over its entire path (the entire transcription
unit or the portion that the RNAP completed if the simulation ended before it
reached the unit end). The seventh column lists the position of the RNAP at the
close of the simulation (will be the length of the unit if it finished before
the end of the simulation); the eighth column lists the time at which the RNAP
began active transcription; the ninth column lists the time at which the RNAP
finished transcribing the first gene on the transcript; and the tenth column
lists the time at which reached the end of the transcription unit.
Columns 2-6 are for simulations in which promoters were shutoff at a user-
determined timepoint (or for comparison with such simulations); column 2 lists
the RNAP's position at the time of shutoff (or the equivalent if the promoter
was not shut off--value is 0 if initiation occurred after this time); column 3
is its elongation rate before the shutoff time or its equivalent. For
columns, 4-6 see source code.


 - **An RNA log:**
  
	rna_info_log_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file lists for each RNA that appeared at any point in the simulation
(in order): the index of the RNA, the time at which its 5' end appeared, the
time at which it matured, the time at which its degradation began, and the time
at which it was completely degraded. The next column lists the number of
ribosomes that initiated on its first gene in its nascent state, the number of
ribosomes that completed translation of its first gene in its nascent state, the
number of ribosomes that initiated translation of the first gene on the mature
mRNA, and the number of ribosomes that compeleted translation of the first gene
on the mature mRNA. If the transcription unit contains more than one gene,
successive blocks of four columns will give the same ribosomal information for
gene 2, gene 3, etc.


 - **A snapshot of simulation RNAs:**

	rna_ribosome_snapshot_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file provides a snapshot of all extant RNAs at the close of the simulation.
Reported (in order) for each RNA are: the index of the promoter from which the
RNA originated, the index of the RNA, the current 5' end of the RNA, the current
3' end of the RNA, the status of the RNA (0=mature,1=nascent), the number of
ribosomes currently translating the RNA, and the positions of each of those
ribosomes.


 - **A snapshot of the template DNA:**

	dna_rnap_snapshot_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file provides a snapshot of the DNA template at the close of the simulation.
Reported (in order) are the index of the template, the number of RNAPs currently
engaged with the template, and the positions of each of these RNAPs.



 - **A "riboseq" file, giving ribosome densities as a function of RNA template position:**

	pseudo_riboseq_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file lists the relative ribosomal occupany of each of the available
positions on the simulation mRNA sequence. The format of each line of the
file is "position  reads", where the number of reads is determined by
ribosome positional checks carrried out at a frequency specified in the
system settings input file.
	

 - **A "RNAPseq" file, giving RNAP densities as a function of DNA template position:**

	pseudo_RNAPseq_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file lists the relative occupany of RNAPs at each of the available
positions on the simulation DNA sequence. The format of each line of the
file is "position  reads", where the number of reads is determined by
RNAP positional checks carrried out at a frequency specified in the
system settings input file.

	
 - **A transcriptional and translational trajectory file, for movie-making:**

	tx_tsl_traj.for_xtc.txt

	This file provides information about the positions of RNAPs and ribosomes
throughout the simulation. These can be used to make movies representing
the simulation trajectory (See [Making movies from trajectories](#post-simulation-visualization-making-movies-from-trajectories) below)


 - **A file for making kymographs (if requested):**

	kymo_dump_rotate

	This file provides information about the position and rotation of RNAPs,
as well as the local supercoiling of the DNA template. It can be used to
make kymographs (see [Making kymographs](#post-simulation-visualization-making-kymographs) below) or for making "zoomed-in"
movies showing RNAP rotation and inter-RNAP supercoiling
(see [Making "zoomed-in" movies](#post-simulation-visualization-making-zoomed-in-movies) below)


 - **A file for making RNAP dwell-time distributions (if requested):**

	dwt_dump_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file lists the set of dwell times in windows of the size (in bp) specified
in the system settings input file (if dwell times have been requested).


 - **A file for making single-molecule traces of RNAP progress (if requested):**

	single_mol_trace_NUMBER.txt

		where NUMBER is the numerical label for the trajectory assigned with
		the " -traj " flag

	This file is used for making RNAP traces for comparison with traces made from
single-molecule experiments. The file lists in each line simulation time, RNAP
template position, and the 3' end of the mRNA produced by the RNAP--plotting
the last two values can be an effective way to visualize RNAP backtracking (in
the backtracked state, the 3' RNA end will remain at the farthest downstream
position an RNAP has achieved even as the enzyme itself moves upsstream).
"one RNAP run" and "single molecule trace" both need to be selected in the
system settings input file if you want to generate a trace.

<br/>

----------------
### Troubleshooting
----------------

The simulation runner has been written to reduce demand on the stack, but if
simulations make it through some set-up steps but run into a segmentation
fault before the simulation proper begins, you can try removing limits on
stack size with:

        ulimit -s unlimited

Current stack size can be checked with `ulimit -a`. Simulations seems to run
fine with default (8 MB) sizes but it might be helpful to remove the limit.

Simulations have been set so that the maximum number of RNAs present at any
one time in the simulation is capped at 10000, the maximum of RNAs made over
the whole trajectory is 100000, and the maximum number of proteins produced
over the entire simulation is 2000000. Overshooting these limits can cause
run-time problems. These limites are hardcoded (in SRC/generate_trajcetory.c)
but can be changed without effect on simulations (at the price of additional
RAM usage).









<br/>







## POST-SIMULATION VISUALIZATION: MAKING MOVIES FROM TRAJECTORIES




After a simulation run is complete, you can use output files from that run to visualize
the simulation trajectory. Movies made using the following commands represent all
transcription, translation, and RNA degradation that took place during the simulation
in a "whole operon" view: all active RNAPs and ribosomes, as well as all extant RNAs,
appear in each simulation frame.


<p align="center">
	<img width="800" height="530" src="https://github.com/Elcock-Lab/spotter/blob/main/IMG/220526_movie_example.png">
</p>

**Snapshot from a simulation of the _E. coli alaS_ operon**
<br/>
<br/>
<br/>



*REQUIREMENTS:* VMD (Visual Molecular Dynamics) molecular graphics viewer

### Note on VMD
------------

Making and watching movies does *NOT* require familiarity with VMD. The scripts described
below for use on VMD are essentially "self-activating" so that users will only need to
download and install the VMD software in order to watch trajectories.

If you need to download VMD, it is available at:

https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

The easiest way to watch movies is to install VMD on your desktop and copy folders from
machines on which simulations have been run to your workstation.

For a Windows machine, download version 1.9.3 for:

	Windows OpenGL (32-bit Intel x86) 

For a Mac, download version 1.9.3 for:

	MacOS X OpenGL (32-bit Intel x86) 

Installation is straightforward from the site listed above.

***IF YOU ARE USING AN OLDER VERSION OF VMD (1.9.2 or older)***, you may need to change the
  format of the psf files generated by the trajectory maker.

  If VMD can't read the psf made by trajectory_movie_maker, use the `psf-converter` script
  (located in SRC in the main directory--use `chmod +x psf_converter` to make it executable
  if it's not already) to change to the older format:

	./psf_converter  psf_file_name

  and replace the new-format psf with the identically-named old-format psf file by
  copying the output `converted.psf` file into the appropriately-named file.


### Making trajectory files for visualization
------------------------------------------

The movie-maker converts the simulation trajectory into a visualizable form, assigning spatial
positions to RNAPs, ribosomes, and mRNAs as a function of simulation time.

The key file for making a movie of a particular simulation trajectory is:

	tx_tsl_traj.for_xtc.txt

To make the files required for the VMD movie, run:

	./trajectory_movie_maker  trajectory_file  tsl_input_file  200  50  0.1  traj_label

where:

	trajectory_file			...is the simulation output trajectory file
					(tx_tsl_traj.for_xtc.txt or whatever you've renamed it)

	tsl_input_file			...is the translation input file you supplied for
					the simulation

	traj_label			...is a name you assign to the trajectory to be
					visualized

The additional numbers (200, 50, 0.1) define the shape of the DNA template included in the movie.
To fit in a movie frame, the template is represented as a widening helical track: the values
represent, respectively, the radius of the starting helix, the increase in this radius per
360-degree circuit, and the vertical drop per step (all in bp). The values given should work
well for most operons, but can be adjusted by trial and error for a better on-screen fit
as needed.
			
The output of the movie maker is a directory called:

	MOVIE_MATERIALS.traj_label

where the label is the one given above. The directory contains the files you need to make your
movie. These are:

 - a pdb file, with the initial structure of the system
 - a psf file, with bonds needed to plot the system
 - an xtc file, with the simulation trajectory
 - a tcl file, with a custom script for displaying the system in VMD.

Each of the files will be labeled with the unique label you assigned in running the program.

***IMPORTANT NOTE:*** Each simulation run--even of the same gene/transcription unit--will create different
pdb, psf, xtc, and tcl files. You CANNOT mix and match files. If you want to look at different
runs of the same system, you need to give each run its own "traj_label". If you want to look at
several runs of dusB-Fis, for example, you might want to name them dusB_Fis_1, dusb_Fis_2, etc.

To proceed to movie making, copy the MOVIE_MATERIALS directory to your local workstation.


### Visualizing the trajectory in VMD
----------------------------------

Once you have the directory on your workstation, do the following:

 - **Open the pdb file with VMD.**

	NOTE: If your implementation of VMD automatically creates a representation
(e.g., a ruler) so that the loaded pdb file will not be molecule 0, start your VMD session by opening
the pdb with VMD (that is, do not open VMD, then open the pdb from your session). If your version does not
autoload any representations at start-up, the order is inconsequential: you can either load your pdb
into an already-running session of VMD or open VMD with the opening of your pdb.

 - **Once the structure is loaded, open the VMD Tk console**


	The Tk console can be found in the drop-down menu nder "Extensions" in the VMD main window.

 - **Check your current location in the Tk console**


	If you are not already in the directory where your pdb is located, change to that directory.
You can check your current directory in the Tk console by entering "pwd" in the Tk console.
If the directory listed is not where your pdb file is located, change to that directory with:

                cd FILE_PATH/

	where `FILE_PATH` is the complete name of the folder containing the `.pdb` file


 - **Enter the following command:**

        source custom_script_for_NAME.tcl

	where NAME is the label mentioned above (which is also included in the directory name).


 - **Then enter:**


        start_tx_tsl_viewer NUMBER COLOR_NAME

	where NUMBER is "1" (use different ribosome colors in a polycistronic transcript) or "0" (do not use
different ribosome colors) and COLOR_NAME is the background color you would like (choose from the list
of names included in the VMD ColorID list)


Once the viewing script has been run, the trajectory can be played by clicking on the "play" button
at the lower right of the VMD "Main" window.

You can center the image by clicking on the visualization window and zooming in and out; simulation
info should remain centered and legible irrespective of zoom-in.

<br/>


## POST-SIMULATION VISUALIZATION: MAKING "ZOOMED-IN" MOVIES


In addition to the "big picture" movies designed to capture transcriptional dynamics
at the operon level (see documentation on whole-operon movies in this folder), users
can also use simulation output to make movies focused on a smaller section of the
DNA template. These movies are limited to a 500 bp section of the template and do not
include RNA or ribosomes but are well-suited for visualizing RNAP rotation and pausing
and local changes in supercoiling density.


*REQUIREMENTS:* VMD (Visual Molecular Dynamics)

*See* [Making movies from trajectories: Note on VMD](#note-on-vmd) *for further information*

### Making trajectory files for visualization
------------------------------------------

Like the whole-operon movie maker, the "zoomed-in" movie maker converts the
simulation trajectory into a visualizable form through the assignment of spatial
coordinates.

The key file for making zoomed-in movies of a trajectory is:

	kymo_dump_rotate

This file is produced in simulation in which users select the "make_kymograph"
option in simulation input files. It contains information about RNAP template
position, state, and extent of rotation, as well as base pair-resolution
information about supercoiling density. By default, the kymograph is set to
cover a 50-second portion of simulation time at 0.01-second resolution over
a template span of 500 bp. The 5000 resulting frames are a reasonable upper for
visualization, and the 0.01-second resolution is a reasonable minimum for capturing
rotation. The zoom-in movie maker currently includes a hard-coded assumption
that the portion of the template covered in the kymograph file is 500 bp.

To make the files required for a VMD movie, run:

	./zoomin_movie_maker traj_label kymo_file

where:

	traj_label			...is a unique label for assign to the 
					trajectory to be visualized

	kymo_file			...is the simulation output kymograph file
					(kymo_dump_rotate or whatever you've renamed it)

The output of the movie maker is a directory called:

	ZOOMIN_MOVIE_MATERIALS.traj_label

where "traj_label" is the label mentioned above. The directory will contain (like a 
whole-operon movie directory) four files:

 - a pdb file, with the initial structure of the system
 - a psf file, with bonds needed to plot the system
 - an xtc file, with the simulation trajectory
 - a tcl file, with a custom script for displaying the system in VMD.

Each of these files will be labeled with the "traj_label" identifier.

To proceed with moviemaking, copy this directory to your local workstation.


### Visualizing the trajectory in VMD
----------------------------------

Once you have the directory on your local workstation, do the following:

 - **Open the pdb file with VMD.**


	NOTE: If your implementation of VMD automatically creates a representation (e.g., a ruler)
so that the loaded pdb file will not be molecule 0, start your VMD session by opening
the pdb with VMD (that is, do not open VMD, then open the pdb from your session). If your
version does not autoload any representations at start-up, the order is inconsequential:
you can either load your pdb into an already-running session of VMD or open VMD with the 
opening of your pdb.

 - **Once the structure is loaded, open the VMD Tk console**

	The Tk console can be found in the drop-down menu nder "Extensions" in the VMD main window.

 - **Check your current location in the Tk console**

	If you are not already in the directory where your pdb is located, change to that directory.
You can check your current directory in the Tk console by entering "pwd" in the Tk console.
If the directory listed is not where your pdb file is located, change to that directory with:

                cd FILE_PATH/

	where "FILE_PATH" is the complete name of the folder containing the .pdb file


 - **Enter the following command:**


	`source zoomin_viewer.NAME.tcl`

	where NAME should be replaced by the "traj_label" you assigned above.

 - **To start the viewer, type:**

	`start_zoomin_movie `

	The movie can be played with the "play" buttom at lower right in the VMD main window.

In the trajectory, DNA coloration is an indicator of supercoiling density (sigma), with
white indicating relaxed DNA, blue representing negatively-supercoiled (underwound) DNA,
and red representing positively-supercoiled (overwound) DNA. RNAPs appear with two lobes
in order to make their rotation apparent.

<br/>

## POST-SIMULATION VISUALIZATION: MAKING KYMOGRAPHS


Simulation output can be used to make kymographs that illustrate changes in RNAP position
and the evolution of supercoiling density in a 2-D plot:



<p align="center">
	<img width="538" height="480" src="https://github.com/Elcock-Lab/spotter/blob/main/IMG/example_kymo.png">
</p>

Although the kymograph maker will provide data that can be used for plotting in a number of
graphing programs, the short custom script it generates is written specifically for R.

Users without R can download both R and the (free) RStudio Desktop package from links at:

	https://www.rstudio.com/products/rstudio/download/#download


### Making datafiles for the kymograph
-----------------------------------

The file--written out at the end of a simulation run--needed to make a kymograph is:

	kymo_dump_rotate

This file is produced in simulations in which users select the "make_kymograph"
option in simulation input files. It contains information about RNAP template
position, state, and extent of rotation, as well as base pair-resolution
information about supercoiling density. By default, reporting is set to
cover a 50-second portion of simulation time at 0.01-second resolution over
a template span of 500 bp. These values can be varied in input files but 
kymographs including significantly more than 500 bp or 50 s will result in
larger files and reduced clarity in kymograph plots.

To generate the materials you'll need to plot a kymograph, run:

	./kymograph_maker  kymo_file  500  5  10 0.01 50.01  traj_label

where:

        kymo_file                       ...is the simulation output kymograph file
                                        (kymo_dump_rotate or whatever you've renamed it)

        traj_label                      ...is a unique label you assign to the
                                        trajectory to be visualized

The numbers between these values determine the coarseness of the kymograph representation
(reducing the number of data points from 2500000 at default resolution). The numbers
(500, 5, 10, 0.01, 50.01) set the number of DNA bp, the number of bp to be included per
plotted point, the number of frames to be included per plotted point, the time between
frames in the kymo_file, and time at which to stop plotting data, respectively. The 
values listed work well for default settings in simulations.

Output of the kymograph is a directory,

	KYMOGRAPH_INFO.traj_label
	
where "traj_label" is the trajectory label used above. The directory contains three files:

 - a file with RNAP positional information
 - a file with supercoiling density information
 - a short R script for plotting information from the other two files

Each of these files is labeled with the "traj_label" identifier.


### Plotting the kymograph in R
----------------------------

To make a kymograph from the simulation trajectory, copy the the KYMOGRAPH_INFO folder to
your desktop.

**With RStudio:**

In an open session of RStudio, load the R script in the folder you've copied over

	custom_script_for_X.R		 (where "X" is the traj_label above)

using "Open File" from the dropdown menu under "File"

From the console, type

	CTL + Shift + Enter

to execute all lines of the script. You can also execute script lines one at a time by
placing the cursor on the first line and hitting CTL+Enter, which executes the current
line and advances to the next line.

**Without RStudio:**

The script can be run in an R console without RStudio by entering:

	source("C:/XXX/XXX/KYMO_INFO.name/custom_script_for_name.R",echo=TRUE)

where "XXX" etc. represents the file path and "name" is your trajectory label.

If you have any difficulties executing commands, check that you're in the right directory
in the console with:

	getwd()

and change to the proper directory with

	setwd("C:/XXX/XXX/KYMO_INFO.name/")

in order to read files from the correct location.

	
Modifications to the representation can be made by changing values in the ggplot command:

 - min and max sigma values can be changed with the `limits` term in
	  `scale_fill_gradient2` to compress or widen the scale

 - axes can be changed by altering "limits" terms in `scale_x_continuous` (for
	  the template position axis) and `scale_y_continuous` (for the time axis)


<br/>

## POST-SIMULATION VISUALIZATION: MAKING PLOTS OF SIMULATION RNAS

The transcription and translation trajectory file produced during a simulation can
also be used to create plots illustrating the production, translational usage, and
degradation of individual mRNA molecules.

Files for making three-dimensional plots of RNA production and translation can be
made with:

	./rna_trajectory_plotter tx_tsl_traj.for_xtc.txt 10 RNA_TRAJ_NAME

where the first argument, "tx_tsl_traj.for_xtc.txt", is the trajectory file mentioned
above used in making operon-scale movies, the second argument is the number of RNAs
to be plotted in the 3-D rendering (here, 10; values >30 van be difficult to visualize),
and the third argument is the name of the directory that will contain files necessary
for R visualization.

To display the plot, copy the directory made on executing the command above to your workstation,
load the "rna_and_ribosome_plotter.R" script and in RStudio run it. Note that the lines
designed to identify and reset the current working directory may not work on all systems. In this
case, the lines of the script can be executed one by one with the current directory set
manually at the point indicated by the comments included in the script. Once all commands have
been executed, the plot can be displayed by entering:

	multiRNA_3d_plot

Each RNA displayed will be displayed as a bold line with length at a given time indicated on
the z-axis; the position of translating ribosomes on that mRNA as a function of time is indicated by the
the thinner lines (one line for each ribosome) is also indicated by z-axis position. 
Production of mRNAs is illustrated in the increasing z-axis position; RNAs are mature
in the flat region of the RNA plot; and the progress of degradation is indicated by decreasing
RNA length in the period following.

<br/>

## EXTERNAL CONTRIBUTIONS

The programs for visualization of simulation trajectories included in this repository rely on the xtc
file format developed by Erik Lindahl and David van der Spoel. Files enabling the use of this format
from these authors have been included without modification. They are (in the INCL directory):

	xdrfile.h
	xdr_file_xtc.h

and (in the SRC directory):

	xdrfile.c
	xdrfile_xtc.c

Additional information can be found in the comments provided by the authors at the top of each
of the listed files.




