#computes lattice vectors from Lammps box definition
function latvec = boxdim2latvec(boxdim)

  latvec = zeros(3,3);
  latvec(1,1) = boxdim(2)-boxdim(1);
  latvec(2,2) = boxdim(4)-boxdim(3);
  latvec(3,3) = boxdim(6)-boxdim(5);
  latvec(1,2) = boxdim(7);
  latvec(1,3) = boxdim(8);
  latvec(2,3) = boxdim(9);

endfunction
