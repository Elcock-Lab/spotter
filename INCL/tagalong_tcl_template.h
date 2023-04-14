#ifndef TAGALONG_H
#define TAGALONG_H

void Assign_tcl_template_string(char *tcl_str)
{

sprintf(tcl_str,
"\n\
proc load_bonds_and_traj {} {\n\
	global traj_id;\n\
\n\
	set loaded_up [molinfo num]\n\
	if {$loaded_up == 0} {\n\
		set pdb_filename [format \"tx_tsl.phosDNA_march.two_lobe_ribo.proteins.macrohelix.multigene.%%s.pdb\" $traj_id]\n\
		mol new $pdb_filename type {pdb} first 0 last -1 step 1 waitfor 1\n\
		animate style Loop\n\
	}\n\
	set psf_filename [format \"tx_tsl.phosDNA_march.two_lobe_ribo.proteins.macrohelix.multigene.%%s.psf\" $traj_id]\n\
	set xtc_filename [format \"tx_tsl.phosDNA_march.two_lobe_ribo.proteins.macrohelix.multigene.%%s.xtc\" $traj_id]\n\
	mol addfile $psf_filename type {psf} first 0 last -1 step 1 waitfor all 0\n\
	animate style Loop\n\
	mol addfile $xtc_filename type {xtc} first 0 last -1 step 1 waitfor all 0\n\
	animate style Loop\n\
}\n\
\n\
proc set_initial_view {} {\n\
	global back_color;\n\
	global gene_color;\n\
\n\
	for {set i 0} {$i < 10} {incr i} {\n\
		set j [expr {$i + 1}];\n\
		set k [expr {$i + 11}];	\n\
		set gene_label($j) [format \"C%%d\" $i]\n\
		set gene_label($k) [format \"O%%d\" $i]\n\
	}\n\
	set gene_label(21) \"NE\" \n\
	set gene_label(22) \"NZ\" \n\
	set gene_label(23) \"NH\"\n\
	set gene_color(0) 20\n\
	set gene_color(1) 1\n\
	set gene_color(2) 0\n\
	set gene_color(3) 3\n\
	set gene_color(4) 11\n\
	set gene_color(5) 12\n\
	set gene_color(6) 6\n\
	set gene_color(7) 4\n\
	set gene_color(8) 15\n\
	set gene_color(9) 27\n\
\n\
	display projection Orthographic\n\
	axes location Off\n\
	color Display Background $back_color\n\
\n\
	mol modselect 0 0 name P\n\
	mol modstyle 0 0 CPK 2.2 5.3 12.0 12.0\n\
	mol modcolor 0 0 Beta\n\
	mol addrep 0\n\
	mol modselect 1 0 name N\n\
	mol modstyle 1 0 CPK 2.2 16.3 12.0 12.0\n\
	mol modcolor 1 0 Beta\n\
	mol addrep 0\n\
	mol modselect 2 0 name CB\n\
	mol modstyle 2 0 CPK 3.0 5.3 12.0 12.0\n\
	mol modcolor 2 0 ColorID 7\n\
	mol addrep 0\n\
	mol modselect 3 0 name CG\n\
	mol modstyle 3 0 CPK 3.0 5.3 12.0 12.0\n\
	mol modcolor 3 0 ColorID 1\n\
	mol addrep 0\n\
	mol modselect 4 0 name CD\n\
	mol modstyle 4 0 CPK 2.9 8.5 12.0 12.0\n\
	mol modcolor 4 0 ColorID 7\n\
	mol addrep 0\n\
	mol modselect 5 0 name CE\n\
	mol modstyle 5 0 CPK 2.9 8.5 12.0 12.0\n\
	mol modcolor 5 0 ColorID 1\n\
	mol addrep 0\n\
	mol modselect 6 0 name C\n\
	mol modstyle 6 0 VDW 10.0 42\n\
	mol modcolor 6 0 ColorID 23\n\
	mol modmaterial 6 0 BrushedMetal\n\
	material change ambient Ghost 0.15\n\
	material change diffuse Ghost 0.39\n\
	material change specular Ghost 0.34\n\
	material change shininess Ghost 0.3\n\
	material change opacity Ghost 1.0\n\
	color change rgb 17 1.0 0.0 0.03\n\
	mol addrep 0\n\
	mol modselect 7 0 name O\n\
	mol modstyle 7 0 VDW 11.5 42\n\
	mol modcolor 7 0 ColorID 17\n\
	mol modmaterial 7 0 Ghost\n\
	mol addrep 0\n\
	mol modselect 8 0 name S\n\
	mol modstyle 8 0 VDW 6.0 42\n\
	mol modcolor 8 0 ColorID 17\n\
	mol modmaterial 8 0 Ghost\n\
	mol addrep 0\n\
	mol modselect 9 0 name CA\n\
	mol modstyle 9 0 VDW 1.5 42\n\
	mol modcolor 9 0 ColorID 11\n\
	mol modmaterial 9 0 Opaque\n\
	for {set i 1} {$i < 24} {incr i} {\n\
		mol addrep 0\n\
		set rep_index [expr {$i + 9}]\n\
		mol modselect $rep_index 0 name $gene_label($i)\n\
		mol modstyle $rep_index 0 CPK 2.9 16.5 12.0 12.0\n\
		mol modcolor $rep_index 0 ColorID $gene_color([expr {($i - 1) %% 10}])\n\
	}\n\
		\n\
}\n\
	\n\
\n\
proc add_ribocolor_reps {} {\n\
	global num_genes;\n\
	global gene_color;\n\
	global y_offset;\n\
	\n\
	set curr_rep [expr {[molinfo 0 get numreps] - 1}]\n\
	for {set i 1} {$i <= $num_genes} {incr i} {\n\
		mol addrep 0\n\
		incr curr_rep\n\
		mol modselect $curr_rep 0 name O and floor(mass+0.001) == [expr {1000 + $i}]\n\
		mol modstyle $curr_rep 0 VDW 12.0 42\n\
		mol modcolor $curr_rep 0 ColorID $gene_color([expr {($i - 1) %% 10}])\n\
		mol modmaterial $curr_rep 0 Ghost\n\
		mol selupdate $curr_rep 0 1\n\
		mol addrep 0\n\
		incr curr_rep\n\
		mol modselect $curr_rep 0 name S and floor(mass+0.001) == [expr {1000 + $i}]\n\
		mol modstyle $curr_rep 0 VDW 6.5 42\n\
		mol modcolor $curr_rep 0 ColorID $gene_color([expr {($i - 1) %% 10}])\n\
		mol modmaterial $curr_rep 0 Ghost\n\
		mol selupdate $curr_rep 0 1\n\
	}\n\
\n\
}\n\
		\n\
proc remove_frames_and_adjust {frame_cut time_adj} {\n\
  global frame_offset\n\
  global time_reset\n\
  global frame_cut_made\n\
\n\
  animate delete beg 0 end $frame_cut skip 0 0\n\
  set time_reset $time_adj\n\
  set frame_offset $frame_cut\n\
  set frame_cut_made 1\n\
}\n\
\n\
proc start_tx_tsl_viewer {ribocolor backgr} {\n\
  global vmd_frame;\n\
  global min;\n\
  global sec;\n\
  global hdth;\n\
  global curr_rna;\n\
  global nascent_rna;\n\
  global mature_rna;\n\
  global ribo;\n\
  global protein;\n\
  global time_loc;\n\
  global ribo_loc;\n\
  global prot_loc;\n\
  global curr_rna_loc;\n\
  global nascent_rna_loc;\n\
  global mature_rna_loc;\n\
  global zoom_adjust;\n\
  global y_offset\n\
  global left_edge\n\
  global frame_offset\n\
  global frame_cut_made\n\
  global time_reset\n\
  global num_genes\n\
  global num_ribo_adds;\n\
  global ribo_adjust;\n\
  global ribo_id;\n\
  global back_color;\n\
  global traj_id;\n\
  global prev_frame;\n\
  global lines;\n\
  global num_lines;\n\
\n\
  set back_color $backgr\n\
  set traj_id \"TRAJIC\"\n\
  load_bonds_and_traj\n\
  set_initial_view\n\
\n\
  load_traj_line_array\n\
  for {set h 1} {$h <= $num_lines} {incr h} {\n\
	if {$lines($h) eq \"\"} continue\n\
        set entries [regexp -inline -all -- {\\S+} $lines($h)]\n\
        set first [lindex $entries 0]\n\
	set min($first) [lindex $entries 1]\n\
	set sec($first) [lindex $entries 2]\n\
	set hdth($first) [lindex $entries 3]\n\
	set curr_rna($first) [lindex $entries 4]\n\
	set nascent_rna($first) [lindex $entries 5]\n\
	set mature_rna($first) [lindex $entries 6]\n\
	set ribo($first) [lindex $entries 7]\n\
#	set protein($first) [lindex $entries 8]\n\
	set num_genes [lindex $entries 8]\n\
	set last_entry [expr {8 + $num_genes}]\n\
	for {set i 9} {$i <= $last_entry} {incr i} {\n\
		set protein($first,[expr {$i - 8}]) [lindex $entries $i]\n\
	}\n\
	incr last_entry\n\
	set num_ribo_adds($first) [lindex $entries $last_entry]\n\
	incr last_entry\n\
	for {set i 1} {$i <= $num_ribo_adds($first)} {incr i} {\n\
		set ribo_adjust($first,$i,1) [lindex $entries $last_entry]\n\
		incr last_entry\n\
		set ribo_adjust($first,$i,2) [lindex $entries $last_entry]\n\
		incr last_entry\n\
		set ribo_adjust($first,$i,3) [lindex $entries $last_entry]\n\
		incr last_entry\n\
	}\n\
   }\n\
\n\
\n\
  for {set i MIN_INDEX} {$i <= MAX_INDEX} {incr i} {\n\
 	set j [expr {$i - 1}]\n\
	set ribo_id($i) [atomselect 0 \"index $j\"]\n\
	$ribo_id($i) global\n\
  }\n\
\n\
  set time_reset 0\n\
  set frame_offset 0\n\
  set frame_cut_made 0\n\
  set left_edge LEFT_EDGE\n\
  set y_offset Y_START\n\
  set zoom_adjust 1.0\n\
  set min(0) 0\n\
  set sec(0) 0\n\
  set hdth(0) 0\n\
  set curr_rna(0) 0\n\
  set nascent_rna(0) 0\n\
  set mature_rna(0) 0\n\
  set ribo(0) 0\n\
  for {set i 1} {$i <= $num_genes} {incr i} {\n\
    set protein(0,$i) 0\n\
  }\n\
\n\
  if {$ribocolor > 0 && $num_genes > 1} {\n\
    add_ribocolor_reps;\n\
  }\n\
\n\
  if {$back_color eq \"black\"} {\n\
  	graphics 0 color 9\n\
  } else {\n\
  	graphics 0 color 16\n\
  }\n\
  set y_loc [expr $y_offset - 50.0]\n\
  graphics 0 text [list [expr {-1.0 * $left_edge}] $y_loc 0] \"time:\" size 1.5 thickness 1.5\n\
  graphics 0 text [list [expr {-1.0 * $left_edge + 270.0}] $y_loc 0] \"current RNA:\" size 1.5 thickness 1.5\n\
  set y_loc [expr $y_offset - 70.0]\n\
  graphics 0 text [list [expr {-1.0 * $left_edge}] $y_loc 0] \"ribosomes on:\" size 1.5 thickness 1.5\n\
  graphics 0 text [list [expr {-1.0 * $left_edge + 270.0}] $y_loc 0] \"nascent RNA:\" size 1.5 thickness 1.5\n\
  set y_loc [expr $y_offset - 90.0]\n\
  if {$num_genes == 1} {\n\
    graphics 0 text [list [expr {-1.0 * $left_edge}] $y_loc 0] \"proteins made:\" size 1.5 thickness 1.5\n\
  } else {\n\
    graphics 0 text [list [expr {-1.0 * $left_edge}] $y_loc 0] \"proteins made (1):\" size 1.5 thickness 1.5\n\
  }\n\
  graphics 0 text [list [expr {-1.0 * $left_edge + 270.0}] $y_loc 0] \"mature RNA:\" size 1.5 thickness 1.5\n\
  for {set i 2} {$i <= $num_genes} {incr i} {\n\
    set y_loc [expr {$y_loc - 20.0}];\n\
    set prot_label [format \"proteins made (%%d)\" $i]\n\
    graphics 0 text [list [expr {-1.0 * $left_edge}] $y_loc 0] $prot_label size 1.5 thickness 1.5\n\
  }\n\
  set y_loc [expr $y_offset - 50.0]\n\
  set hh [format \"%%02s:%%02s.%%.1s\" 0 0 0]\n\
  set time_loc [list [expr {-1.0 * $left_edge + 162.0}] $y_loc 0]\n\
  graphics 0 text $time_loc $hh size 1.5 thickness 1.5\n\
  set y_loc [expr $y_offset - 70.0]\n\
  set ii [format \"%%8s\" 0]\n\
  set ribo_loc [list [expr {-1.0 * $left_edge + 143.0}] $y_loc 0]\n\
  graphics 0 text $ribo_loc $ii size 1.5 thickness 1.5\n\
  for {set i 1} {$i <= $num_genes} {incr i} {\n\
    set y_loc [expr {($y_offset - 90.0) - (($i - 1) * 20.0)}]\n\
    set prot_loc($i) [list [expr {-1.0 * $left_edge + 143.0}] $y_loc 0]\n\
    graphics 0 text $prot_loc($i) $ii size 1.5 thickness 1.5\n\
  }\n\
  set y_loc [expr $y_offset - 50.0]\n\
  set curr_rna_loc [list [expr {-1.0 * $left_edge + 413.0}] $y_loc 0]\n\
  graphics 0 text $curr_rna_loc $ii size 1.5 thickness 1.5\n\
  set y_loc [expr $y_offset - 70.0]\n\
  set nascent_rna_loc [list [expr {-1.0 * $left_edge + 413.0}] $y_loc 0]\n\
  graphics 0 text $nascent_rna_loc $ii size 1.5 thickness 1.5\n\
  set y_loc [expr $y_offset - 90.0]\n\
  set mature_rna_loc [list [expr {-1.0 * $left_edge + 413.0}] $y_loc 0]\n\
  graphics 0 text $mature_rna_loc $ii size 1.5 thickness 1.5\n\
\n\
  set prev_frame 0\n\
\n\
   trace add variable ::vmd_frame([molinfo top]) write frame_changer	\n\
   trace add variable ::vmd_logfile write zoom_changer	\n\
}\n\
\n\
\n\
proc disabletrace {} {\n\
  global vmd_frame;\n\
  trace remove variable ::vmd_frame([molinfo top]) write shitter\n\
  trace remove variable ::vmd_logfile write shitter\n\
}\n\
\n\
proc reposition_labels {} {\n\
	global y_offset\n\
	global zoom_adjust\n\
  	global time_loc;\n\
  	global ribo_loc;\n\
  	global prot_loc;\n\
  	global curr_rna_loc;\n\
  	global nascent_rna_loc;\n\
  	global mature_rna_loc;\n\
	global left_edge;\n\
	global num_genes\n\
\n\
  	set y_loc [expr $y_offset - (50.0 * $zoom_adjust)]\n\
	set x_adj [expr {-270.0 * ($zoom_adjust - 1.0)}]\n\
	graphics 0 replace 1\n\
  	graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj}] $y_loc 0] \"time:\" size 1.5 thickness 1.5\n\
	graphics 0 replace 2\n\
  	graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj + (270.0 * $zoom_adjust)}] $y_loc 0] \"current RNA:\" size 1.5 thickness 1.5\n\
  	set y_loc [expr $y_offset - (70.0 * $zoom_adjust)]\n\
	graphics 0 replace 3\n\
  	graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj}] $y_loc 0] \"ribosomes on:\" size 1.5 thickness 1.5\n\
	graphics 0 replace 4\n\
  	graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj + (270.0 * $zoom_adjust)}] $y_loc 0] \"nascent RNA:\" size 1.5 thickness 1.5\n\
  	set y_loc [expr $y_offset - (90.0 * $zoom_adjust)]\n\
	graphics 0 replace 5\n\
	if {$num_genes < 2} {\n\
  	   graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj}] $y_loc 0] \"proteins made:\" size 1.5 thickness 1.5\n\
	} else {\n\
  	   graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj}] $y_loc 0] \"proteins made (1):\" size 1.5 thickness 1.5\n\
	}\n\
	graphics 0 replace 6\n\
  	graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj + (270.0 * $zoom_adjust)}] $y_loc 0] \"mature RNA:\" size 1.5 thickness 1.5\n\
  	for {set i 2} {$i <= $num_genes} {incr i} {\n\
	   set replace_index [expr {6 + ($i - 1)}]\n\
	   graphics 0 replace $replace_index\n\
	   set y_loc [expr $y_offset - ((90.0 + (($i - 1) * 20.0)) * $zoom_adjust)]\n\
	   set prot_label [format \"proteins made (%%d):\" $i]\n\
  	   graphics 0 text [list [expr {-1.0 * $left_edge + $x_adj}] $y_loc 0] $prot_label size 1.5 thickness 1.5\n\
	}\n\
	  \n\
\n\
	set time_loc [list [expr {-1.0 * $left_edge + $x_adj + (162.0 * $zoom_adjust)}]  [expr {$y_offset - (50.0 * $zoom_adjust)}] 0]\n\
	set ribo_loc [list [expr {-1.0 * $left_edge + $x_adj + (143.0 * $zoom_adjust)}]  [expr {$y_offset - (70.0 * $zoom_adjust)}] 0]\n\
	for {set i 1} {$i <= $num_genes} {incr i} {\n\
	   set prot_loc($i) [list [expr {-1.0 * $left_edge + $x_adj + (143.0 * $zoom_adjust)}]  [expr {$y_offset - ((90.0 + (($i - 1) * 20.0))  * $zoom_adjust)}] 0]\n\
	}\n\
	set curr_rna_loc [list [expr {-1.0 * $left_edge + $x_adj + (413.0 * $zoom_adjust)}]  [expr {$y_offset - (50.0 * $zoom_adjust)}] 0]\n\
	set nascent_rna_loc [list [expr {-1.0 * $left_edge + $x_adj + (413.0 * $zoom_adjust)}]  [expr {$y_offset - (70.0 * $zoom_adjust)}] 0]\n\
	set mature_rna_loc [list [expr {-1.0 * $left_edge+ $x_adj + (413.0 * $zoom_adjust)}]  [expr {$y_offset - (90.0 * $zoom_adjust)}] 0]\n\
}\n\
\n\
proc update_production_numbers {stuff stuff2 stuff3 stuff4 stuff5 stuff6 stuff7 stuff8 shite} {\n\
	upvar $stuff name\n\
	upvar $stuff2 name2\n\
	upvar $stuff3 name3\n\
	upvar $stuff4 name4\n\
	upvar $stuff5 name5\n\
	upvar $stuff6 name6\n\
	upvar $stuff7 name7\n\
	upvar $stuff8 name8\n\
  	global time_loc;\n\
  	global ribo_loc;\n\
  	global prot_loc;\n\
  	global curr_rna_loc;\n\
  	global nascent_rna_loc;\n\
  	global mature_rna_loc;\n\
	global zoom_adjust;\n\
	global y_offset;\n\
	global left_edge;\n\
	global frame_offset\n\
	global time_reset\n\
	global num_genes\n\
\n\
	set hh [format \"%%02s:%%02s.%%.1s\" $name($shite) $name2($shite) $name3($shite)]\n\
	if {$time_reset==1} {\n\
		set adj_frame [expr {$shite - ($frame_offset - 1)}];\n\
		set hh [format \"%%02s:%%02s.%%.1s\" $name($adj_frame) $name2($adj_frame) $name3($adj_frame)];\n\
	}\n\
	set replace_index [expr {7 + ($num_genes - 1)}]\n\
	graphics 0 replace $replace_index\n\
	incr replace_index\n\
	graphics 0 text $time_loc $hh size 1.5 thickness 1.5\n\
	set ii [format \"%%8s\" $name4($shite)]\n\
	graphics 0 replace $replace_index\n\
	incr replace_index\n\
	graphics 0 text $ribo_loc $ii size 1.5 thickness 1.5\n\
	for {set i 1} {$i <= $num_genes} {incr i} {\n\
	   set ii [format \"%%8s\" $name5($shite,$i)]\n\
	   graphics 0 replace $replace_index\n\
	   incr replace_index\n\
	   graphics 0 text $prot_loc($i) $ii size 1.5 thickness 1.5\n\
	}\n\
	set ii [format \"%%8s\" $name6($shite)]\n\
	graphics 0 replace $replace_index\n\
	incr replace_index\n\
	graphics 0 text $curr_rna_loc $ii size 1.5 thickness 1.5\n\
	set ii [format \"%%8s\" $name7($shite)]\n\
	graphics 0 replace $replace_index\n\
	incr replace_index\n\
	graphics 0 text $nascent_rna_loc $ii size 1.5 thickness 1.5\n\
	set ii [format \"%%8s\" $name8($shite)]\n\
	graphics 0 replace $replace_index\n\
	graphics 0 text $mature_rna_loc $ii size 1.5 thickness 1.5\n\
}\n\
\n\
proc frame_changer { name element op } {\n\
  global vmd_frame;\n\
  global min;\n\
  global sec;\n\
  global hdth;\n\
  global curr_rna;\n\
  global nascent_rna;\n\
  global mature_rna;\n\
  global ribo;\n\
  global protein;\n\
  global time_loc;\n\
  global ribo_loc;\n\
  global prot_loc;\n\
  global curr_rna_loc;\n\
  global nascent_rna_loc;\n\
  global mature_rna_loc;\n\
  global left_edge;\n\
  global frame_offset\n\
  global frame_cut_made;\n\
  global num_genes\n\
  global num_ribo_adds;\n\
  global ribo_adjust;\n\
  global ribo_id;\n\
  global prev_frame;\n\
\n\
  set finder [expr {$vmd_frame([molinfo top]) + $frame_offset}]\n\
  set checkster [expr {$prev_frame + 1}]\n\
  if {$finder != $checkster} {\n\
   for {set h 1} {$h < $finder} {incr h} { \n\
    for {set i 1} {$i <= $num_ribo_adds($h)} {incr i} {\n\
	$ribo_id($ribo_adjust($h,$i,1)) set mass [expr {1000.0 + $ribo_adjust($h,$i,3)}] \n\
	$ribo_id($ribo_adjust($h,$i,2)) set mass [expr {1000.0 + $ribo_adjust($h,$i,3)}]\n\
    }\n\
   }\n\
  }\n\
  for {set i 1} {$i <= $num_ribo_adds($finder)} {incr i} {\n\
	$ribo_id($ribo_adjust($finder,$i,1)) set mass [expr {1000.0 + $ribo_adjust($finder,$i,3)}] \n\
	$ribo_id($ribo_adjust($finder,$i,2)) set mass [expr {1000.0 + $ribo_adjust($finder,$i,3)}]\n\
  }\n\
  set prev_frame $finder \n\
  update_production_numbers min sec hdth ribo protein curr_rna nascent_rna mature_rna $finder\n\
}\n\
\n\
proc zoom_changer { name element op } {\n\
	global vmd_logfile;\n\
	global zoom_adjust;\n\
  	global vmd_frame;\n\
  	global min;\n\
  	global sec;\n\
  	global hdth;\n\
  	global curr_rna;\n\
  	global nascent_rna;\n\
  	global mature_rna;\n\
  	global ribo;\n\
  	global protein;\n\
  	global time_loc;\n\
  	global ribo_loc;\n\
  	global prot_loc;\n\
  	global curr_rna_loc;\n\
  	global nascent_rna_loc;\n\
  	global mature_rna_loc;\n\
	global y_offset\n\
	global left_edge;\n\
	global frame_offset;\n\
	\n\
	set wordz [regexp -inline -all -- {\\S+} $vmd_logfile]\n\
	set word [lindex $wordz 0]\n\
	if {$word==\"scale\" || $word==\"display\"} {\n\
		set hep [molinfo 0 get scale_matrix];\n\
		set hep2 [lindex $hep 0 0 0];\n\
		set zoom_adjust [expr {0.00573 / $hep2}];\n\
		reposition_labels;\n\
  		set finder [expr {$vmd_frame([molinfo top]) + $frame_offset}]\n\
  		update_production_numbers min sec hdth ribo protein curr_rna nascent_rna mature_rna $finder\n\
	}\n\
\n\
}");

}

#endif
