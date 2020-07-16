clear all;
crash_dumps_octave_core(0);
warning off;
addpath(genpath("octave-src"));

######## input ########
inputPath = "example/MD-sims/";
outputPath = "example/SEC-ML/";
#######################

### reading the *.rdf, *.msd and *.eng files
vecs = readfvecs(inputPath);

### reading seed.data to build the atom types based on LJ params
load([inputPath,'seed.data.oct']);
vecs.types = [seed.PairCoeffs(:,2:end),seed.Masses(:,2:end)];
atoms = mass2atom(vecs.types(:,3));

### naming the atoms based on atom name, second index to differentiate like-atoms
buf = histc(vecs.types(:,3),unique(vecs.types(:,3)));
idx = [];
for i=1:numel(buf),
  idx = [idx,1:buf(i)];
endfor
vecs.atomnames = strcat(mass2atom(vecs.types(:,3)),cellstr(num2str(idx'))');

### saving vecs
save('-v7',[outputPath,'vecs.oct'],'vecs');

