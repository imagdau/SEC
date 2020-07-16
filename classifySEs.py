import sys
sys.path.append('python-src/');

import warnings
warnings.filterwarnings('ignore');

import embed as em
import seplot as sp
import joblib

import importlib
importlib.reload(em)
importlib.reload(sp)

import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import umap
import hdbscan

from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

import glob

#====== input =======
sigLi = 1.40
epsLi = 0.40
inputPath = "example/SEC-ML/"
outputPath = "example/SEC-ML/"
Nnbs = 400
Dmin = 0.1
Nclst = 200
Pcut = 0.35
#====================

### reading vecs
d = em.loadVecs(inputPath+'vecs.oct');

### computing type dependent Rcut
sigAt = d['types'][:,1];
epsAt = d['types'][:,0];
Rmin = (2*(epsAt*sigAt**12+epsLi*sigLi**12)/(epsAt*sigAt**6+epsLi*sigLi**6))**(1/6);
DR = Rmin*(1-2**(-1/6));
Rcut = Rmin+3.0*DR;

### embedding the feature vectors in latent space and classifying solvation environments
d = em.solvEnvGlob(d, Nnbs, Dmin, Nclst, Pcut, Rcut, 'allrdfs');
sp.plotEmbeddingGlob(d, outputPath+'gembedTabs.png', 'trjlabs', 'tcmap');
sp.plotEmbeddingGlob(d, outputPath+'gembedGabs.png', 'glabs_train', 'gcmap');

### constructing binding energy histrograms according to SE membership and plotting
d = em.potEngHis(d, 'glabs_train');
sp.plotEngHis(d, outputPath+'histPotEng.png', 'gcmap');

### plotting rdfs and cdfs according to SE membership
sp.plotDfClasif(d, outputPath+'rdfsGabs.png', 'rdf', 'glabs_train', 'gcmap', 1.0);
sp.plotDfClasif(d, outputPath+'cdfsGabs.png', 'cdf', 'glabs_train', 'gcmap', 1.0);

### plotting trajectories (mean square displacement or binding eng) based on SE membership
sp.plotTrajs(d, outputPath+'allpotsGabs.png', [10,10], 'allpots', 'glabs_train', 'gcmap');
sp.plotTrajs(d, outputPath+'allmsdsGabs.png', [10,10], 'allmsds', 'glabs_train', 'gcmap');

