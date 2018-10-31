1;
clear all
more off
H_BL
quantity

z = 3
m = 7
t = 1
mu = 0.001
%%a = 2

fid = fopen("density.txt", "w");
fid_2 = fopen("density_2.txt", "w");
for a = 0:4
  a
  [H Lattice] = H_bc(z, m, t, mu, a);
  %%[H Lattice] = H_sc(z, m, t, mu, a);
  tic;
  [V E] = eig(H);
  toc
  E = diag(E);
  
  %% all sites
  n_eigen = prod(size(E));
  fprintf(fid, "#a = %f\n", a);
  for j=1:n_eigen
    rho = Density_GS(j, V, E);
    fprintf(fid, "%d: %.15e\n", j, rho);
  endfor
  fprintf(fid, "\n\n");
  
  %% each layer
  fprintf(fid_2, "#a = %f\n", a);
  for j=1:m
    x = Lattice{j}(1);
    rho = Density_GS(x, V, E);

    if j==m
      x2 = Lattice{1}(1);
    else
      x2 = Lattice{j+1}(1);
    endif
    C = Corr_GS(x, x2, V, E);

    fprintf(fid_2, "%d: %.15e %.15e\n", j, rho, C);
  endfor
  fprintf(fid_2, "\n\n");

endfor

fclose(fid);
fclose(fid_2);
