#reads a block of data from Lammps data file
function data = fgetblk(fid)

  data = [];
  ln = strtrim(fgetl(fid));
  while !isempty(ln),
    buf = strsplit(ln);
    cnt = num2str(size(buf,2));
    if isfield(data,cnt),
      data.(cnt) = [data.(cnt); buf];
    else
      data.(cnt) = buf;
    endif
    if feof(fid),
      break;
    endif
    ln = strtrim(fgetl(fid));
  endwhile

endfunction