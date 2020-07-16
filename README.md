The code supplied here can be used to classify and characterize solvation einvironments in the condensed phase, based on all-atom Lammps simulations.
Details of this work are published at zzzzz. The analysis for ProDOT-2Hex neat, presented in the paper, is uploaded here in the example directory.

The project consists of three components:


1) generating single-ion trajectories, starting from the an equilibrium simulations using <genMDinput.m>
  
  The code takes as input <seed.data> and <Li.data> (example/input) which are Lammps topology files.
  
  The code generates <md-*.data> files where the ion (in this case Li) is initialized at random wihtin the equilibrium <seed> simulation box.
  
  The trajectories are carried out with Lammps, using the template input script supplied in <lammps-scripts>. Each run generates *.rdf *.eng and *.msd files which correspond to radial distribution functions, binding energy and mean square displacement.
  
  A simple vmd script is supplied (vmd-scripts) to visualize the Lammps data files using: vmd -e visdata.tcl -args illus.data.
  
  
2) reading the output trajectories and packinging the data into objects <vecs> using <anlMDoutput.m>
  
  The code reads the Lammps output from <example/MD-sims> and creates object <vecs> which is written to <example/SEC-ML>
  
  
3) classifying the solvation environments (SEs) based on rdf feature vectors using <classifySEs.py>
  
  The code reads the <vecs> object and uses UMAP to embed the data in a 2D latent space.
  
  The code then labels the clusters and obtains a classification of SEs.
  
  The code clusters the rdfs, msds, engs based on SE membership and creates useful plots in <example/SEC-ML>
  
  
