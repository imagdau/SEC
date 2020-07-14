#adds N molecules to original datin Lammps data file
#molecules have random orientations and ensure a spherical exclusion volume controlled by minDist
function [datout,datcol] = addSolvent(datin,solvent,N,minDist=1.0)

  mat = [datin.xhi-datin.xlo, 0.0, 0.0;
        datin.xy, datin.yhi-datin.ylo, 0.0;
        datin.xz, datin.yz, datin.zhi-datin.zlo];

  datin.Atoms(:,5:7) -= mean(datin.Atoms(:,5:7));
  coords = datin.Atoms(:,5:7);
  internal = coords*inv(mat);
  
  SolvCoords = solvent.Atoms(:,5:7);
  SolvCM = mean(SolvCoords,1);
  SolvCoords -= SolvCM;
  SolvRG = max(sqrt(sum(SolvCoords.^2,2)));
  SolvInternal = SolvCoords*inv(mat);

  N1 = size(internal,1);
  N2 = size(SolvInternal,1);
  NewCoords = [];
  
  cnt = 0;
  while cnt < N,
    ROT = randRot();
    TRS = rand([1,3])-0.5;
    SolvNewCoords = SolvCoords*ROT+TRS*mat;
    
    vecs = permute(internal,[1,3,2])-permute(SolvNewCoords*inv(mat),[3,1,2]);
    pbc = floor(vecs+0.5);
    dists = sqrt(sum((reshape(vecs-pbc,[N1*N2,3])*mat).^2,2));
    
    if min(dists) >= minDist,
      internal = [internal;SolvNewCoords*inv(mat)];
      NewCoords = cat(3,NewCoords,SolvNewCoords);
      cnt++;
      N1+=N2;
      printf('%d/%d\n',cnt,N);
      fflush(stdout);
    endif
  endwhile
  
  datout = datin;
  solvent.Atoms(:,8:10) = 0;
  for i=1:N,
    solvent.Atoms(:,5:7) = NewCoords(:,:,i);
    datout = comblmp(datout,solvent); #WARNING: velocities are discarder and the box modified
    datcol{i} = comblmp(datin,solvent);
    
    #change box back
    datcol{i}.xlo = datin.xlo; datcol{i}.xhi = datin.xhi;
    datcol{i}.ylo = datin.ylo; datcol{i}.yhi = datin.yhi;
    datcol{i}.zlo = datin.zlo; datcol{i}.zhi = datin.zhi;
    datcol{i}.xy = datin.xy; datcol{i}.xz = datin.xz; datcol{i}.yz = datin.yz;
  endfor
  
  #change box back
  datout.xlo = datin.xlo; datout.xhi = datin.xhi;
  datout.ylo = datin.ylo; datout.yhi = datin.yhi;
  datout.zlo = datin.zlo; datout.zhi = datin.zhi;
  datout.xy = datin.xy; datout.xz = datin.xz; datout.yz = datin.yz;
  
endfunction
