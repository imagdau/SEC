import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import umap
import hdbscan
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

#loads the vecs object which contains rdfs, engs, msds
def loadVecs(fname):
  vecs = sio.loadmat(fname, struct_as_record=True)['vecs'];
  d = dict();
  d['R'] = vecs['R'][0,0][:,0];
  d['Nbins'] = vecs['Nbins'][0,0][0,0].astype(int);
  d['Ntypes'] = vecs['Ntypes'][0,0][0,0].astype(int);
  d['Nsteps'] = vecs['Nsteps'][0,0][0,0].astype(int);
  d['Nfvecs'] = vecs['Nfvecs'][0,0][0,0].astype(int);
  d['Ntrajs'] = vecs['Ntrajs'][0,0][0,0].astype(int);
  d['trjlabs'] = vecs['trjlabs'][0,0][0,:].astype(int)-1;
  d['tcmap'] = matplotlib.cm.rainbow(np.arange(d['Ntrajs'])/(d['Ntrajs']-1));
  d['types'] = vecs['types'][0,0];
  d['atnams'] = np.hstack(vecs['atomnames'][0,0][0,:]);
  d['allrdfs'] = vecs['allrdfs'][0,0].transpose();
  d['allcdfs'] = vecs['allcdfs'][0,0].transpose();
  d['allpots'] = vecs['allpots'][0,0].transpose();
  d['allmsds'] = vecs['allmsds'][0,0].transpose();
  d['maxrdfs'] = vecs['maxrdfs'][0,0][:,0];
  d['maxcdfs'] = vecs['maxcdfs'][0,0][:,0];
  return d

#embeds the feature vectors in 2D and clusters the solvation environments
def solvEnvGlob(d, umap_num_nbs, umap_min_dst, hdb_clst_size, hdb_prob_cut, Rmax, keyname='allrdfs'):
  embed = True;
  if ('umap_num_nbs_glob' in d) & ('umap_min_dst_glob' in d) & ('keyname' in d) & ('gRmax' in d):
    if (d['umap_num_nbs_glob'] == umap_num_nbs) & (d['umap_min_dst_glob'] == umap_min_dst) & (d['keyname'] == keyname) & (sum(~(d['gRmax'] == Rmax)) == 0):
      embed = False;
  if embed == True:
    d['umap_num_nbs_glob'] = umap_num_nbs;
    d['umap_min_dst_glob'] = umap_min_dst;
    d['keyname'] = keyname;
    d['gRmax'] = Rmax;
    mask = np.tile(d['R'],d['Ntypes'])<=np.tile(d['gRmax'],(d['Nbins'],1)).transpose().reshape(-1);
    data = d[keyname][:,mask];
    d['Rmask'] = mask;
    d['gembedding'] = umap.UMAP(n_neighbors=umap_num_nbs, min_dist=umap_min_dst, verbose=True).fit_transform(data);
  d['hdb_clst_size_glob'] = hdb_clst_size;
  cluster = hdbscan.HDBSCAN(min_cluster_size=hdb_clst_size, gen_min_span_tree=False, allow_single_cluster=True).fit(d['gembedding']);
  d['glabs'] = cluster.labels_;
  d['gprbs'] = cluster.probabilities_;
  d['glabs_train'] = np.copy(d['glabs']);
  d['glabs_train'][d['gprbs'] <  hdb_prob_cut] = -1;
  d['gNclst'] = np.max(d['glabs'])+1;
  d['gcmap'] = matplotlib.cm.rainbow(np.arange(d['gNclst'])/(d['gNclst']-1));
  print('Data clustered: '+str(sum(d['glabs']>=0))+', Number of clusters: '+str(d['gNclst']));
  return d

#builds binding energy histograms based on SEC membership
def potEngHis(d, keyname = 'glabs', nbin = 200):
  d['allcounts'], d['binedges'] = np.histogram(d['allpots'].flatten(), nbin);
  d['bincenters'] = (d['binedges'][1:] + d['binedges'][:-1])/2;
  d['gcounts'], _ = np.histogram(d['allpots'][d[keyname]>=0,:].flatten(), d['binedges']);
  d['lcounts'] = dict();
  for i in range(np.max(d[keyname])+1):
    d['lcounts'][i], _ = np.histogram(d['allpots'][d[keyname]==i,:].flatten(), d['binedges']);
  return d

