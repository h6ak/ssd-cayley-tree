1;
clear all
more off
CayleyTree
Hamiltonian
quantity

global z = 2
global m = 30
global t = 1
global mu = 0
global opt = 'sc'

%% OUTPUT FILE
function fid = OpenOutput(fname)
  global z;
  global m;
  global t;
  global mu;
  global opt;

  fid = fopen(fname, "w");
  fprintf(fid, "## Parameters ##\n");
  fprintf(fid, "# z = %d\n", z);
  fprintf(fid, "# m = %d\n", m);
  fprintf(fid, "# t = %.15e\n", t);
  fprintf(fid, "# mu = %.15e\n", mu);
  fprintf(fid, "# opt = %s\n", opt);
  fprintf(fid, "#\n#\n", opt);
endfunction

info_date = datevec(now);
dirtag = [int2str(z) "-" int2str(m) "-" opt]
dirtag = ["Data/" dirtag "/"];
mkdir(dirtag);
ftag = ["_mu-" num2str(mu)];
ftag = [ftag ".txt"];

fname = [dirtag "energy" ftag];
f_energy = OpenOutput(fname);

fname = [dirtag "dens_all" ftag];
f_dens1 = OpenOutput(fname);
fname = [dirtag "dens_path" ftag];
f_dens2 = OpenOutput(fname);

fname = [dirtag "BondCorr_all" ftag];
f_bc1 = OpenOutput(fname);
fname = [dirtag "BondCorr_path" ftag];
f_bc2 = OpenOutput(fname);

%% make Cayley Tree
[Sites Bonds] = DefCayleyTree(z, m, opt);

for a = 0:4
  a
  H = H_Hopping(Sites, Bonds, opt, t, mu, a);
  tic;
  [V E] = eig(H);
  toc
  E = diag(E);
  
  %% function handle
  Dens = @(x) Density_GS(x, V, E);
  Corr = @(x, y) Corr_GS(x, y, V, E);

  %%%% Energy
  fid = f_energy;
  temp = [1:prod(size(E))];
  fprintf(fid, "#a = %f\n", a);
  fprintf(fid, "%d %.15e\n", [temp; E']);
  fprintf(fid, "\n\n");
  clear temp

  %%%% Density
  %% all sites
  temp = [1:prod(size(E))];
  rho = arrayfun(Dens, temp);

  fid = f_dens1;
  fprintf(fid, "#a = %f\n", a);
  fprintf(fid, "%d %.15e\n", [temp; rho]);
  fprintf(fid, "\n\n");
  clear temp
  
  %% each layer
  temp = [1:size(Sites)(2)];
  temp2 = arrayfun(@(x)Sites{x}(1), temp);
  rho = rho(temp2);

  fid = f_dens2;
  fprintf(fid, "#a = %f\n", a);
  fprintf(fid, "%d %.15e\n", [temp; rho]);
  fprintf(fid, "\n\n");
  clear temp temp2 rho

  %%%% BondCorr
  %% all bonds
  BondMat = cell2mat(Bonds);
  temp = [1:size(BondMat)(2)];
  C = arrayfun(Corr, BondMat(1,:), BondMat(2,:));

  fid = f_bc1;
  fprintf(fid, "#a = %f\n", a);
  fprintf(fid, "%d %.15e\n", [temp; C]);
  fprintf(fid, "\n\n");
  clear BondMat temp C

  %% each layer
  temp = [1:size(Bonds)(2)];
  x1 = arrayfun(@(x)Bonds{1, x}(1), temp);
  x2 = arrayfun(@(x)Bonds{2, x}(1), temp);
  C = arrayfun(Corr, x1, x2);
  
  fid = f_bc2;
  fprintf(fid, "#a = %f\n", a);
  fprintf(fid, "%d %.15e\n", [temp; C]);
  fprintf(fid, "\n\n");
  clear temp x1 x2 C
endfor

fclose("all");
clear f_energy;
clear f_dens1 f_dens2;
clear f_bc1 f_bc2;
