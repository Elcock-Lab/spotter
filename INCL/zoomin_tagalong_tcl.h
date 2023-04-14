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
                set pdb_filename [format \"zoomin_movie.%%s.pdb\" $traj_id]\n\
                mol new $pdb_filename type {pdb} first 0 last -1 step 1 waitfor 1\n\
                animate style Loop\n\
        }\n\
        set psf_filename [format \"zoomin_movie.%%s.psf\" $traj_id]\n\
        set xtc_filename [format \"zoomin_movie.%%s.xtc\" $traj_id]\n\
        mol addfile $psf_filename type {psf} first 0 last -1 step 1 waitfor all 0\n\
        animate style Loop\n\
        mol addfile $xtc_filename type {xtc} first 0 last -1 step 1 waitfor all 0\n\
        animate style Loop\n\
}\n\
\n\
\n\
proc set_initial_view {} {\n\
	display projection Orthographic\n\
	axes location Off\n\
\n\
	mol modselect 0 0 index < 1000 \n\
	mol modstyle 0 0 CPK 1.9 2.2 12.0 12.0\n\
	mol representation CPK 1.9 2.2 12.0 12.0\n\
	mol modcolor 0 0 ColorID 0\n\
	for {set i 1} {$i < 20} {incr i} {\n\
		set p [molinfo 0 get numreps]\n\
		set q [expr $i + 1]\n\
		if {$q > $p} {\n\
			mol addrep 0 \n\
		} else {\n\
			mol modrep $i 0\n\
		}\n\
		mol modcolor $i 0 ColorID $i\n\
		mol modselect $i 0 index > 10000\n\
	}\n\
	mol addrep 0\n\
	mol modselect 20 0 index 1000 \n\
	mol modstyle 20 0 VDW 12.2 22\n\
	mol modcolor 20 0 ColorID 20\n\
	mol representation VDW 12.2 22\n\
	set rnap_index 1001\n\
	for {set i 21} {$i < 31} {incr i} {\n\
		mol addrep 0\n\
		mol modcolor $i 0 ColorID $i\n\
		set rnap_top [expr $rnap_index + 15]\n\
		mol modselect $i 0 index $rnap_index or index $rnap_top\n\
		mol smoothrep 0 $i 5\n\
		incr rnap_index\n\
	}\n\
}\n\
	\n\
\n\
proc start_zoomin_movie {} {\n\
  global vmd_frame;\n\
  global start;\n\
  global stop;\n\
  global coloring;\n\
  global RNAP_changed;\n\
  global RNAP_color;\n\
  global traj_id;\n\
  global lines;\n\
  global lines2;\n\
  global numlines;\n\
  global numlines2;\n\
\n\
  set traj_id \"TRAJIC\"\n\
  load_bonds_and_traj\n\
  set_initial_view\n\
\n\
  load_super_traj_array\n\
\n\
  for {set h 1} {$h <= $numlines} {incr h} {\n\
        if {$lines($h) eq \"\"} continue\n\
        set entries [regexp -inline -all -- {\\S+} $lines($h)]\n\
        set first [lindex $entries 0]\n\
        set second [lindex $entries 1]\n\
        set third [lindex $entries 2]\n\
        set fourth [lindex $entries 3]\n\
        set fifth [lindex $entries 4]\n\
	set start($first,$second) $third\n\
	set stop($first,$second) $fourth\n\
	set coloring($first,$second) $fifth\n\
  }\n\
\n\
  load_rnap_traj_array\n\
\n\
  for {set h 1} {$h <= $numlines2} {incr h} {\n\
        if {$lines2($h) eq \"\"} continue\n\
        set entries [regexp -inline -all -- {\\S+} $lines2($h)]\n\
        set first [lindex $entries 0]\n\
        set second [lindex $entries 1]\n\
        set third [lindex $entries 2]\n\
        set fourth [lindex $entries 3]\n\
	set RNAP_changed($first,$second) $third\n\
	set RNAP_color($first,$second) $fourth\n\
   }\n\
\n\
   trace add variable ::vmd_frame([molinfo top]) write frame_changer	\n\
}\n\
\n\
\n\
proc disabletrace {} {\n\
  global vmd_frame;\n\
  trace remove variable ::vmd_frame([molinfo top]) write shitter\n\
}\n\
\n\
\n\
proc make_color_changes {stuff stuff2 stuff3 stuff4 stuff5 shite} {\n\
\n\
	upvar $stuff name\n\
	upvar $stuff2 name2\n\
	upvar $stuff3 name3\n\
	upvar $stuff4 name4\n\
	upvar $stuff5 name5\n\
\n\
	set ctr 0\n\
\n\
	foreach n [array names name \"$shite,*\"] {\n\
\n\
#		puts \"START: $name($n); STOP: $name2($n); COLOR: $name3($n)\"\n\
\n\
		set opp_low [expr 1001 - $name2($n)]\n\
		set opp_high [expr 1001 - $name($n)]\n\
\n\
#		mol modselect $ctr 0 index >= $name($n) and index <= $name2($n)\n\
		mol modselect $ctr 0 (index >= $name($n) and index <= $name2($n)) or (index >= $opp_low and index <= $opp_high)\n\
\n\
		if {$name3($n) > 0} {\n\
			set b 1\n\
#			set r [expr $name3($n) * -1.0]\n\
			set r [expr 1.0 - $name3($n)]\n\
			set g $r\n\
		} else {\n\
			set r 1\n\
#			set b $name3($n)\n\
#			set b [expr $name3($n) - 1.0]\n\
			set b [expr 1.0 + $name3($n)]\n\
			set g $b\n\
		}\n\
\n\
#		puts \"R: $r; G: $g; B: $b\"\n\
\n\
		color change rgb $ctr $r $g $b\n\
		incr ctr\n\
\n\
	}\n\
	for {set i $ctr} {$i < 20} {incr i} {\n\
		mol modselect $i 0 index > 10000\n\
	}\n\
	foreach n [array names name4 \"$shite,*\"] {\n\
#		set b 0\n\
#		set r 1\n\
		set b 1\n\
		set r 0\n\
#		if {$name5($n) == 4} {\n\
#			set g 0\n\
#		} elseif {$name5($n) == 3} {\n\
#			set g 0.59\n\
#		} elseif {$name5($n) == 2} {\n\
#			set g 1\n\
#		} else {\n\
#			set r 0\n\
#			set g 0.75\n\
#		}\n\
		set which_rnap [expr $name4($n) + 19]\n\
		color change rgb $which_rnap $r $g $b\n\
	}\n\
\n\
}\n\
\n\
proc frame_changer { name element op } {\n\
  global vmd_frame;\n\
  global start;\n\
  global stop;\n\
  global coloring;\n\
  global RNAP_changed;\n\
  global RNAP_color;\n\
\n\
  set finder $vmd_frame([molinfo top])\n\
  make_color_changes start stop coloring RNAP_changed RNAP_color $finder\n\
}\n\
\n");
}

#endif
