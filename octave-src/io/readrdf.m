#reads md*.rdf data files
function [R,rdfs,cdfs] = readrdf(fname)
  
  fid = fopen(fname,'r');
  
  #skip comment lines
  fskipl(fid,3);

  R = [];
  rdfs = [];
  cdfs = [];
  while !feof(fid),
    #read TimeStep and RdfBins
    buf = fgetl(fid);
    aux = str2num(buf);
    TimeStep = aux(1);
    RdfBins = aux(2);
    data = fgetcnk(fid,RdfBins,true);
    if isempty(R),
      R = data(:,2);
    endif
    rdfs = cat(3,rdfs,data(:,3:2:end));
    cdfs = cat(3,cdfs,data(:,4:2:end));
  endwhile
  
  fclose(fid);
  
endfunction
