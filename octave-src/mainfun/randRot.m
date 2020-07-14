#generates random rotation matrix
function Rot = randRot()

  u = rand(3,1)-0.5;
  u = u./sqrt(dot(u,u));
  do
    v = rand(3,1)-0.5;
    v = v./sqrt(dot(v,v));
  until sqrt(dot(u,v)) > 0.9,
  v -= dot(u,v)*u;
  v = v./sqrt(dot(v,v));
  w = cross(u,v);
  Rot = [u,v,w];
  
endfunction
