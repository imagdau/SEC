#computes box lenghts and angles from lattice vectors
function lenang = latvec2lenang(latvec)

  mat = latvec';
  a = mat(1,:);
  b = mat(2,:);
  c = mat(3,:);
  length_a = sqrt(sum(a.^2));
  length_b = sqrt(sum(b.^2));
  length_c = sqrt(sum(c.^2));
  cos_alpha = b*c'/(length_b*length_c);
  cos_beta  = a*c'/(length_a*length_c);
  cos_gamma = a*b'/(length_a*length_b);
  angle_alpha = acos(cos_alpha)*180.0/pi;
  angle_beta  = acos(cos_beta )*180.0/pi;
  angle_gamma = acos(cos_gamma)*180.0/pi;
  lenang = [length_a,length_b,length_c,angle_alpha,angle_beta,angle_gamma,cos_alpha,cos_beta,cos_gamma];
  
endfunction
