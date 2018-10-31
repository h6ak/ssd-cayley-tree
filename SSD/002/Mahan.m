1;
clear all;
more off;
quantity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G.D.Mahan, Phys.Rev.B.63, 155110 (2001)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = H_Mahan (z, m, t)
  x = (-t)*sqrt(z-1)*ones(m-1,1);
  x(1) = (-t)*sqrt(z);
  H = spdiags(x, -1, m, m);
  H = H + H';
endfunction

%% Parameters
z = 3
m = 7
t = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirtag = ["Data/Mahan/"];
mkdir(dirtag);
ftag = ["_" int2str(z) "-" int2str(m) "-sc.txt"];

fname = [dirtag "energy" ftag];
f_energy = fopen(fname, "w");

fname = [dirtag "dens" ftag];
f_dens = fopen(fname, "w");

fname = [dirtag "BondCorr" ftag];
f_bc = fopen(fname, "w");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = H_Mahan (z, m, t);
full(H);

[V E] = eig(H);
E = diag(E);

%% function handle
Dens = @(x) Density_GS(x, V, E);
Corr = @(x, y) Corr_GS(x, y, V, E);

%%%% Energy
fid = f_energy;
temp = [1:prod(size(E))];
fprintf(fid, "%d %.15e\n", [temp; E']);
fprintf(fid, "\n\n");
clear temp

%%%% Density
temp = [1:prod(size(E))];
rho = arrayfun(Dens, temp);

fid = f_dens;
fprintf(fid, "%d %.15e\n", [temp; rho]);
fprintf(fid, "\n\n");
clear temp

%%%% BondCorr
temp = [1:(prod(size(E))-1)];
temp2 = [2:prod(size(E))];
C = arrayfun(Corr, temp, temp2);

fid = f_bc;
fprintf(fid, "%d %.15e\n", [temp; C]);
fprintf(fid, "\n\n");
clear temp temp2

fclose(f_energy);
fclose(f_dens);
fclose(f_bc);

