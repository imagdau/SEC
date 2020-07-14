#writes Lammps data file to fname
function writelmp(fname,data)

  fid = fopen(fname,"w");
  fprintf(fid,"LAMMPS data file Created by Octave\n\n");
  
  #header
  fprintf(fid,"%8d atoms\n",data.atoms);
  fprintf(fid,"%8d bonds\n",data.bonds);
  fprintf(fid,"%8d angles\n",data.angles);
  fprintf(fid,"%8d dihedrals\n",data.dihedrals);
  fprintf(fid,"%8d impropers\n",data.impropers);
  fprintf(fid,"\n");
  fprintf(fid,"%8d atom types\n",data.atom_types);
  fprintf(fid,"%8d bond types\n",data.bond_types);
  fprintf(fid,"%8d angle types\n",data.angle_types);
  fprintf(fid,"%8d dihedral types\n",data.dihedral_types);
  fprintf(fid,"%8d improper types\n",data.improper_types);
  fprintf(fid,"\n");  

  #box
  fprintf(fid,"%12.6f%12.6f xlo xhi\n",data.xlo,data.xhi);
  fprintf(fid,"%12.6f%12.6f ylo yhi\n",data.ylo,data.yhi);
  fprintf(fid,"%12.6f%12.6f zlo zhi\n",data.zlo,data.zhi);
  if (isfield(data,'xy')&isfield(data,'xz')&isfield(data,'yz')),
    fprintf(fid,"%12.6f%12.6f%12.6f xy xz yz\n",data.xy,data.xz,data.yz);
  endif
  
  #coefficients
  fprintf(fid,"\nMasses\n\n");
  fprintf(fid,"%8d%12.3f\n",data.Masses');
  fprintf(fid,"\nPair Coeffs\n\n");
  fprintf(fid,"%8d%12.4f%12.7f\n",data.PairCoeffs');
  if isfield(data,"BondCoeffs"),
    fprintf(fid,"\nBond Coeffs\n\n");
    fprintf(fid,"%8d%12.4f%12.4f\n",data.BondCoeffs');
  endif
  if isfield(data,"AngleCoeffs"),
    fprintf(fid,"\nAngle Coeffs\n\n");
    fprintf(fid,"%8d%12.4f%12.4f\n",data.AngleCoeffs');
  endif
  if isfield(data,"DihedralCoeffs"),
    fprintf(fid,"\nDihedral Coeffs\n\n");
    fmt = repmat("%12.4f",1,size(data.DihedralCoeffs,2)-1);
    fprintf(fid,["%8d",fmt,"\n"],data.DihedralCoeffs');
  endif
  if isfield(data,"ImproperCoeffs"),
    fprintf(fid,"\nImproper Coeffs\n\n");
    fmt = repmat("%12.4f",1,size(data.ImproperCoeffs,2)-3);
    fprintf(fid,["%8d",fmt,"%8d%8d\n"],data.ImproperCoeffs');
  endif
  
  #topology
  fprintf(fid,"\nAtoms\n\n");
  if size(data.Atoms,2) == 7,
    afmt = "%8d%8d%8d%15.7f%20.10g%20.10g%20.10g\n";
  endif
  if size(data.Atoms,2) == 10,
    afmt = "%8d%8d%8d%15.7f%20.10g%20.10g%20.10g%8d%8d%8d\n";
  endif
  fprintf(fid,afmt,data.Atoms');
  if isfield(data,"Velocities");
    fprintf(fid,"\nVelocities\n\n");
    fprintf(fid,"%8d%20.10g%20.10g%20.10g\n",data.Velocities');
  endif
  if isfield(data,"Bonds"),
    fprintf(fid,"\nBonds\n\n");
    fprintf(fid,"%8d%8d%8d%8d\n",data.Bonds');
  endif
  if isfield(data,"Angles"),
    fprintf(fid,"\nAngles\n\n");
    fprintf(fid,"%8d%8d%8d%8d%8d\n",data.Angles');
  endif
  if isfield(data,"Dihedrals"),
    fprintf(fid,"\nDihedrals\n\n");
    fprintf(fid,"%8d%8d%8d%8d%8d%8d\n",data.Dihedrals');
  endif
  if isfield(data,"Impropers"),
    fprintf(fid,"\nImpropers\n\n");
    fprintf(fid,"%8d%8d%8d%8d%8d%8d\n",data.Impropers');
  endif
  fprintf(fid,"\n");
  
  fclose(fid);

endfunction