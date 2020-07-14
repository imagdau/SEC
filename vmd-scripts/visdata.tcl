if { $argc != 1 } {
  puts "Provide a file to view."
  exit
} else {
  set fname [ lindex $argv 0 ]
}

topo readlammpsdata $fname

mol delrep all top
mol representation PaperChain 1.000000 10.000000
mol addrep 0
mol selection not name 1 and not name 2
mol representation Lines 1.0
mol addrep 0
mol selection name 3
mol representation CPK 1.0 0.0 18.0 18.0
mol addrep 0

color Display Background white
display projection Orthographic
display depthcue off
axes location off

display nearclip set   0.0
display farclip set  100.0
pbc box -center com



