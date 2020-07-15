clear all;
crash_dumps_octave_core(0);
warning off;
addpath(genpath("octave-src"));

minDist = 2.3;
N = 3;

### reading Lammps data file for the the equilibrated polymer
if exist('example/input/seed.data.oct','file'),
  load('example/input/seed.data.oct');
else,
  ### caching for large files
  seed = readlmp('example/input/seed.data',true);
  save('example/input/seed.data.oct','seed');
endif
seed = unwrapPBClmp(seed);

### reading Lammps data file for the ion
ion = readlmp('example/input/Li.data',true);

### adding ion to the equilibrium structure
[lmp,lmpcol] = addSolvent(seed,ion,N,minDist);

### writing Lammps data and cif files with all ion positions for illustration
[atoms,coords,latvec] = lmp2xyz(lmp);
writelmp('example/MD-sims/illus.data',lmp);
writecif('example/MD-sims/illus.cif',atoms,coords,latvec);

### writing Lammps data files to be used as input for MD simulations
for j=1:N,
  disp(j);
  fflush(stdout);
  writelmp(sprintf('example/MD-sims/md-%04d.data',j),lmpcol{j});
endfor

### example Lammps command, pass correct variables to md.input
j = 1;
atomstr = sprintf('%s ',mass2atom(lmpcol{j}.Masses(:,2)){:});
litp = find(strcmp(mass2atom(lmpcol{j}.Masses(:,2)),'Li'));
pair_types = 1:lmpcol{j}.atom_types;
pair_types(pair_types==litp) = [];
pair_types = [repmat(litp,size(pair_types));pair_types](:);
ptps = sprintf('%d ',pair_types);
jobname = sprintf("md-%04d",j);
printf('run with: mpirun -np 4 lmp_mpi -var fid %s -var astr "%s" -var litp %d -var ptps "%s" -var seed ${RANDOM} -in md.input\n',jobname,atomstr,litp,ptps);
system('cp lammps-scripts/md.input example/MD-sims/md.input');

