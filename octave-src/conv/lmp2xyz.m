#extracts the atom coordinates and lattice vectors form Lammps data file 
function [atoms,coords,latvec] = lmp2xyz(lmp)

  atoms = mass2atom(sort(lmp.Masses)(lmp.Atoms(:,3),2));
  coords = lmp.Atoms(:,5:7);
  boxdim = [lmp.xlo,lmp.xhi,lmp.ylo,lmp.yhi,lmp.zlo,lmp.zhi,lmp.xy,lmp.xz,lmp.yz];
  latvec = boxdim2latvec(boxdim);
  
endfunction