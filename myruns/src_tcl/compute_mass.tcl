package require topotools
mol new MYB_switchgrass_nch_1.pdb
set selection [atomselect top all]
set sum 0
foreach mass [$selection get mass] {
	set sum [expr {$sum + $mass}]
}	
puts "MW: $sum"
exit
