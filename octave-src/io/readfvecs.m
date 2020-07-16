#reads *.rdf, *.msd and *.eng files and compiles all the inforamtion in vecs 
function vecs = readfvecs(fdir)

  fvecs = [fdir,'/vecs.cache.oct'];
  
  if exist(fvecs),
    tic;
    load('-v7',fvecs);
    toc;
  else,
  
    tic;
    listnames = glob([fdir,'/*.rdf']);
    [fdirs,fbases,fexts] = cellfun('fileparts',listnames,'UniformOutput',false);
    vecs.trjlabs = [];
    vecs.allrdfs = [];
    vecs.allcdfs = [];
    vecs.allmxyz = [];
    vecs.allmsds = [];
    vecs.allkins = [];
    vecs.allpots = [];
    vecs.Ntrajs = numel(fbases);
    
    for i = 1:vecs.Ntrajs,
    
      fpathrdf = [fdir,'/',fbases{i},'.rdf']
      fpathmsd = [fdir,'/',fbases{i},'.msd']
      fpatheng = [fdir,'/',fbases{i},'.eng']
      fflush(stdout);
      
      [vecs.R,rdfs,cdfs] = readrdf(fpathrdf);
      [vecs.Nbins,vecs.Ntypes,vecs.Nfvecs] = size(rdfs);
      vecs.trjlabs = [vecs.trjlabs,repmat(i,1,vecs.Nfvecs)];
      vecs.allrdfs = [vecs.allrdfs,reshape(rdfs,vecs.Nbins*vecs.Ntypes,vecs.Nfvecs)];
      vecs.allcdfs = [vecs.allcdfs,reshape(cdfs,vecs.Nbins*vecs.Ntypes,vecs.Nfvecs)];
      
      if (i==1),
        vecs.maxrdfs = max(max(rdfs,[],1),[],3)';
        vecs.maxcdfs = max(max(cdfs,[],1),[],3)';
      else,
        vecs.maxrdfs = max(vecs.maxrdfs,max(max(rdfs,[],1),[],3)');
        vecs.maxcdfs = max(vecs.maxcdfs,max(max(cdfs,[],1),[],3)');
      endif
      
      buf = load(fpathmsd);
      Nstepsall = size(buf,1)-1;
      vecs.Nsteps = Nstepsall/vecs.Nfvecs;
      msd = reshape(buf(2:end,:),vecs.Nsteps,vecs.Nfvecs,5);
      msdstart = reshape(buf(1:end-1,:),vecs.Nsteps,vecs.Nfvecs,5)(1,:,:);
      msd = cat(1,msdstart,msd);

      vecs.T = msd(:,1,1);
      vecs.allmxyz = cat(2,vecs.allmxyz,msd(:,:,2:4));
      vecs.allmsds = cat(2,vecs.allmsds,msd(:,:,5));
      
      buf = load(fpatheng);
      eng = reshape(buf(2:end,:),vecs.Nsteps,vecs.Nfvecs,3);
      engstart = reshape(buf(1:end-1,:),vecs.Nsteps,vecs.Nfvecs,3)(1,:,:);
      eng = cat(1,engstart,eng);
      vecs.allkins = cat(2,vecs.allkins,eng(:,:,2));
      vecs.allpots = cat(2,vecs.allpots,eng(:,:,3));
      
    endfor

    vecs.Nsteps++; 
    save('-v7',fvecs,'vecs');
    toc;
    
  endif
  
endfunction
