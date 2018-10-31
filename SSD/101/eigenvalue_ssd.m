1;
clear all
more off;
format

# coordination number
z = 3;
# The number of shells
M = 4;
# hopping amplitude
t = 1;
# chemical potential
mu_list = [0, 0.5, 1, 1.5];
# sinsoidal deformation
a = 0;

# -----
N = M;
m = N+1;
param.t = t;
#param.mu = mu;
#param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
#param;

fid = fopen("ev_ssd.txt", "w");
fprintf(fid, "# z=%d, M=%d\n", z, M);
fprintf(fid, "# Datalist: i, E_i(OBC, direct), E_i(OBC, via block), E_i(SSD, direct), E_i(SSD, via block)\n\n");
for mu = mu_list
  a = 0;
  param.mu = mu;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
  ev_full
  ev_block

  n_eigen_full = prod(size(E));
  n_eigen_block = prod(size(E_Array));
  if n_eigen_full != n_eigen_block
    disp("Error!")
    break
  endif
  
  E0 = E;
  E0_Array = E_Array;
  a = 2;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
  ev_full
  ev_block

  fprintf(fid, "# mu = %.15e\n", mu);
  for j=1:n_eigen_full
    fprintf(fid, "%d %.15e %.15e %.15e %.15e\n",
	    j, E0(j), E0_Array(j), E(j), E_Array(j) );
  endfor
  fprintf(fid, "\n\n");
endfor
fclose(fid);
clear fid
