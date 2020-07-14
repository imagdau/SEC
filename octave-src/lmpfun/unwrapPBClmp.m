#unwraps Lammps data file according to image flags
function datout = unwrapPBClmp(datin)

  datout = datin;
  mat = [datin.xhi-datin.xlo, 0.0, 0.0;
        datin.xy, datin.yhi-datin.ylo, 0.0;
        datin.xz, datin.yz, datin.zhi-datin.zlo];
  pbc = datin.Atoms(:,8:10);
  datout.Atoms(:,5:7) += pbc*mat;
  datout.Atoms(:,8:10) = 0;
  
endfunction