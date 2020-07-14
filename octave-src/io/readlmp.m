#reads Lammps data file from fname
function data = readlmp(fname,prt=false)
  
  fid = fopen(fname,'r');
  
  strtrim(fgetl(fid));
  strtrim(fgetl(fid));
  
  #header
  header = fgetblk(fid);
  if numfields(header)==2,
    counts = header.('2');
    types = header.('3');
  else
    counts = header.('2');
    header = fgetblk(fid);
    types = header.('3');
  endif
  
  for i=1:size(counts,1),
    data.(counts{i,2}) = str2num(counts{i,1});
  endfor
  for i=1:size(types,1),
    name = [types{i,2},'_',types{i,3}];
    data.(name) = str2num(types{i,1});
  endfor
  
  #box
  boxdef = fgetblk(fid);
  lohi = boxdef.('4');
  for i=1:size(lohi,1),
    data.(lohi{i,3}) = str2num(lohi{i,1});
    data.(lohi{i,4}) = str2num(lohi{i,2});
  endfor
  
  if isfield(boxdef,'6')
    tilt = boxdef.('6');
    data.(tilt{1,4}) = str2num(tilt{1,1});
    data.(tilt{1,5}) = str2num(tilt{1,2});
    data.(tilt{1,6}) = str2num(tilt{1,3});
  endif
  
  #data
  while !feof(fid),
    label = [struct2cell(fgetblk(fid)){:}{:}];
    label = strsplit(label,"#"){1};
    if prt,
      label
      fflush(stdout);
    endif
    #assumes only numbers, cannot parse dihedral style
    vals = struct2cell(fgetblk(fid)){:};
    data.(label) = cellfun("str2num",vals);
  endwhile
  
  fclose(fid);

endfunction