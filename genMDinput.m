clear all;
crash_dumps_octave_core(0);
warning off;
addpath(genpath("octave-src"));

######## input ########
minDist = 2.3;
N = 100;
inputPath = "example/input/";
outputPath = "example/MD-sims/";
#######################

### reading Lammps data file for the the equilibrated polymer
if exist([outputPath,'seed.data.oct'],'file'),
  load([outputPath,'seed.data.oct']);
else,
  ### caching for large files
  seed = readlmp([inputPath,'seed.data'],true);
  save([outputPath,'seed.data.oct'],'seed');
endif
seed = unwrapPBClmp(seed);

### reading Lammps data file for the ion
ion = readlmp([inputPath,'Li.data'],true);

### adding ion to the equilibrium structure
[lmp,lmpcol] = addSolvent(seed,ion,N,minDist);

### writing Lammps data and cif files with all ion positions for illustration
[atoms,coords,latvec] = lmp2xyz(lmp);
writelmp([outputPath,'illus.data'],lmp);
writecif([outputPath,'illus.cif'],atoms,coords,latvec);

### writing Lammps data files to be used as input for MD simulations
for j=1:N,
  disp(j);
  fflush(stdout);
  fname = sprintf([outputPath,'md-%04d.data'],j);
  writelmp(fname,lmpcol{j});
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

