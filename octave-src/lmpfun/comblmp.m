# #combines two Lammps data files datin1 and datin2 into one datout
function datout = comblmp(datin1,datin2)

  #this code assumes indices in general run from 1:N
  #first part: concatenate the two data files
  datout.atoms     = datin1.atoms + datin2.atoms;
  datout.bonds     = datin1.bonds + datin2.bonds;
  datout.angles    = datin1.angles + datin2.angles;
  datout.dihedrals = datin1.dihedrals + datin2.dihedrals;
  datout.impropers = datin1.impropers + datin2.impropers;
  
  datout.atom_types     = datin1.atom_types + datin2.atom_types;
  datout.bond_types     = datin1.bond_types + datin2.bond_types;
  datout.angle_types    = datin1.angle_types + datin2.angle_types;
  datout.dihedral_types = datin1.dihedral_types + datin2.dihedral_types;
  datout.improper_types = datin1.improper_types + datin2.improper_types;
  
  datout.Masses         = [datin1.Masses;datin2.Masses+[datin1.atom_types,0]];
  datout.PairCoeffs     = [datin1.PairCoeffs;datin2.PairCoeffs+[datin1.atom_types,0,0]];

  bc1 = zeros(0,3); bc2 = zeros(0,3);
  if isfield(datin1,'BondCoeffs'),
    bc1 = datin1.BondCoeffs;
  endif
  if isfield(datin2,'BondCoeffs'),
    bc2 = datin2.BondCoeffs;
  endif
  datout.BondCoeffs = [bc1;bc2+[datin1.bond_types,0,0]];
  
  ac1 = zeros(0,3); ac2 = zeros(0,3);
  if isfield(datin1,'AngleCoeffs'),
    ac1 = datin1.AngleCoeffs;
  endif
  if isfield(datin2,'AngleCoeffs'),
    ac2 = datin2.AngleCoeffs;
  endif
  datout.AngleCoeffs = [ac1;ac2+[datin1.angle_types,0,0]];

  dc1 = zeros(0,5); dc2 = zeros(0,5);
  if isfield(datin1,'DihedralCoeffs'),
    dc1 = datin1.DihedralCoeffs;
  endif
  if isfield(datin2,'DihedralCoeffs'),
    dc2 = datin2.DihedralCoeffs;
  endif
  datout.DihedralCoeffs = [dc1;dc2+[datin1.dihedral_types,0,0,0,0]];

  ic1 = zeros(0,4); ic2 = zeros(0,4);
  if isfield(datin1,'ImproperCoeffs'),
    ic1 = datin1.ImproperCoeffs;
  endif
  if isfield(datin2,'ImproperCoeffs'),
    ic2 = datin2.ImproperCoeffs;
  endif
  datout.ImproperCoeffs = [ic1;ic2+[datin1.improper_types,0,0,0]];
    
  ### very important! connectFUN uses this, do not remove!
  if isempty(datin1.Atoms),
    plusmol = 0;
  else
    plusmol = max(datin1.Atoms(:,2));
  endif;
  datout.Atoms = [datin1.Atoms;datin2.Atoms+[datin1.atoms,plusmol,datin1.atom_types,zeros(1,size(datin2.Atoms,2)-3)]];

  bd1 = zeros(0,4); bd2 = zeros(0,4);
  if isfield(datin1,'Bonds'),
    bd1 = datin1.Bonds;
  endif
  if isfield(datin2,'Bonds'),
    bd2 = datin2.Bonds;
  endif  
  datout.Bonds = [bd1;bd2+[datin1.bonds,datin1.bond_types,datin1.atoms,datin1.atoms]];
  
  an1 = zeros(0,5); an2 = zeros(0,5);
  if isfield(datin1,'Angles'),
    an1 = datin1.Angles;
  endif
  if isfield(datin2,'Angles'),
    an2 = datin2.Angles;
  endif  
  datout.Angles = [an1;an2+[datin1.angles,datin1.angle_types,datin1.atoms,datin1.atoms,datin1.atoms]];
  
  di1 = zeros(0,6); di2 = zeros(0,6);
  if isfield(datin1,'Dihedrals'),
    di1 = datin1.Dihedrals;
  endif
  if isfield(datin2,'Dihedrals'),
    di2 = datin2.Dihedrals;
  endif  
  datout.Dihedrals = [di1;di2+[datin1.dihedrals,datin1.dihedral_types,datin1.atoms,datin1.atoms,datin1.atoms,datin1.atoms]];  
  
  im1 = zeros(0,6); im2 = zeros(0,6);
  if isfield(datin1,'Impropers'),
    im1 = datin1.Impropers;
  endif
  if isfield(datin2,'Impropers'),
    im2 = datin2.Impropers;
  endif  
  datout.Impropers = [im1;im2+[datin1.impropers,datin1.improper_types,datin1.atoms,datin1.atoms,datin1.atoms,datin1.atoms]];

  datout.xlo = min(datout.Atoms(:,5));
  datout.ylo = min(datout.Atoms(:,6));
  datout.zlo = min(datout.Atoms(:,7));

  datout.xhi = max(datout.Atoms(:,5));
  datout.yhi = max(datout.Atoms(:,6));
  datout.zhi = max(datout.Atoms(:,7));

  datout = reducelmp(datout);
  
endfunction
