#reads a chunk of cnt lines from a file
function data = fgetcnk(fid,cnt,num=false)

  data = [];

  for i=1:cnt,
    ln = strtrim(fgetl(fid));
    if num,
      data = [data; str2num(ln)];
    else,
      data = [data; strsplit(ln)];
    endif
  endfor

endfunction