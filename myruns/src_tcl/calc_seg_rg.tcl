# Compute segmental Rg 
mol load psf ../init_files/L.psf ;# you can also use a PDB
set psfval [atomselect top "all"]
set seg [ lsort -unique [ $psfval get segname ] ]

mol addfile traj_npt_main.trr waitfor all 
set num_frames [ molinfo top  get numframes ]
set num_min    0 

# Write outputs
set fch [open "chaindetails_psf.dat" w]
puts $fch " chain/segmentID \t nmons "

set frg [open "rg_allsegs.dat" w]
puts $frg "segmentID \t n \t m \t m-n \t Rg"

# Main code
foreach s $seg {
    set res [lsort -unique [[atomselect top "segname $s"] get resid]]
    set Ntot [llength $res]
    puts $fch "$s \t $Ntot"
    for {set n 1} {$n < $Ntot} {inc n} {
	for {set m [ expr $n + 1 ] } {$m <= $Ntot } { inc m } { 
	    set sel [atomselect top "resid $n to $m and segname $s" ]
	    set separation [  expr abs ( $m - $n) ]
	    for {set i $num_min } { $i <= $num_frames-1  } { incr i } {
		$sel frame $i
		$sel update
		set d [measure rgyr $sel ]
		lappend N_$separation $d
	    } ;#end for {set i $num_min}...
	    set v1 N_$separation
	    set v2 [ expr $$v1 ]
	    set dmean [ expr [vecmean $v2 ] ]
	    puts $frg "$s \t $n \t $m \t $separation \t $dmean"
	    unset v1 
	    unset v2
	    unset dmean
	    unset d
	    unset separation 
	    unset sel
	} ;# end for {set m [expr $n+1]}...
	puts "finished residue $n of segment $s"
    } ;#end for {set n 1}...
    puts "finished all ($Ntot) residues of segment $s"
    unset res
    unset Ntot
} ;#end foreach seg
mol delete top
close $fch

set fp [open "RgvsN.dat" w]
puts $fp "N \t Rgmean \t Rgvar"
for { set l 1 } { $l < $Ntot } { incr l } {
    set var  N_$l
    set var2 [ expr $$var]
    set mean  [ expr [vecmean $var2] ]
    set c [expr $Ntot-$l] 
    set stdev [ expr [ vecstddev $var2]/$c ] ;# set stdev22 [llength $var2]
    puts $fp " $l $mean $stdev "
}
close $fp
