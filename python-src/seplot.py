import scipy.io as sio
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import umap
import hdbscan

#plots the solvation environments in latent space
def plotEmbeddingGlob(d, figname, keyname='glabs', cmapname='gcmap', sz = 2.0, axs = []):
  mask = (d[keyname]>=0);
  plt.scatter(d['gembedding'][~mask,0], d['gembedding'][~mask,1], s=sz, c=[0.7,0.7,0.7], alpha=0.1);
  plt.scatter(d['gembedding'][ mask,0], d['gembedding'][ mask,1], s=sz, c=d[cmapname][d[keyname][ mask],:], alpha=0.3);
  plt.xticks([]);
  plt.yticks([]);
  plt.title('global');
  if axs:
    plt.xlim(axs[0],axs[1]);
    plt.ylim(axs[2],axs[3]);
  plt.savefig(figname, dpi=500);
  plt.close();
  return

#plots the binding energy histograms
def plotEngHis(d, figname, colmap = 'gcmap'):
  plt.plot(d['bincenters'],d['allcounts'],color='black',label='all data');
  plt.plot(d['bincenters'],d['gcounts'],color='gray',label='clustered');
  for i in range(len(d['lcounts'])):
    plt.plot(d['bincenters'],d['lcounts'][i],color=d[colmap][i,:],label=('se'+str(i+1)),alpha=0.7);
  plt.xlabel('Energy (Kcal/mol)');
  plt.ylabel('Counts');
  plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=6, mode="expand", borderaxespad=0.)
  plt.tight_layout();
  plt.savefig(figname, dpi=500);
  plt.close();
  return

#plots the radial distribution functions (rdfs) and cummulative distribution functions (cdfs)
def plotDfClasif(d, figname, df, keyname = 'glabs', colmap = 'gcmap', fct = 1.0):
  clst = np.unique(d[keyname]);
  clst = np.delete(clst,np.where(clst==-1)[0]);
  Nclst = len(clst);
  fig, axs = plt.subplots(d['Ntypes'], Nclst, figsize=(Nclst*2.25,d['Ntypes']*1.5), sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0});
  if (axs.ndim==1):
    axs = axs.reshape(-1,1);
  dfs = d['all'+df+'s'].transpose(1,0).reshape(d['Nbins'],d['Ntypes'],d['Nfvecs']*d['Ntrajs'],order='F').transpose(0,2,1);
  R = d['R'];
  for i in range(Nclst):
    maskSE = (d[keyname] == clst[i]);
    for j in range(d['Ntypes']):
      maskRc = (R <= d['gRmax'][j]);
      data = dfs[:,maskSE,j];
      y_mean = np.mean(data,axis=1);
      y_std = np.std(data,axis=1);
      axs[j,i].plot(R[maskRc], data[maskRc,:], linewidth=0.02, color=d[colmap][clst[i],:], alpha=0.7);
      axs[j,i].fill_between(R[maskRc], np.amax([y_mean[maskRc]-fct*y_std[maskRc], np.zeros(sum(maskRc))],axis=0), y_mean[maskRc]+fct*y_std[maskRc], color=d[colmap][clst[i],:], alpha=0.3, linewidth=0.0);
    axs[0,i].set_title('se'+str(i+1)+' ('+str(sum(maskSE)*100/d[keyname].shape[0])+')');
    axs[d['Ntypes']-1,i].set_xlabel('R (A)');
  for j in range(d['Ntypes']):
    axs[j,0].set_ylabel(df);
    axs[j,Nclst-1].annotate(d['atnams'][j], xy=(1.05,0.5), xycoords='axes fraction');
  plt.tight_layout();
  plt.savefig(figname, dpi=500);
  plt.close();
  return

#plots all trajectories (msd or eng) labeled according to SE membership
def plotTrajs(d, figname, layout, keyname, labs='glabs', cmap='gcmap'):
  fig, axs = plt.subplots(layout[1], layout[0], figsize=(layout[0]*3.0,layout[1]*2.0));  
  for i in range(d['Ntrajs']):
    print(str(i+1)+'/'+str(d['Ntrajs']));
    [hidx,widx] = divmod(i,layout[0]);
    for j in range(i*d['Nfvecs'],(i+1)*d['Nfvecs']):
      tm = np.arange(j*d['Nsteps'],(j+1)*d['Nsteps'])-i*d['Nfvecs']*d['Nsteps'];
      if (d[labs][j]>=0):
        cl = d[cmap][d[labs][j],:];
      else:
        cl = [0.5,0.5,0.5,1.0];
      axs[hidx,widx].plot(tm, d[keyname][j,:], color=cl, linewidth=0.1);
  plt.tight_layout();
  plt.savefig(figname, dpi=500);
  plt.close();
  return

