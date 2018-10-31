%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uniform NN hopping model on Cayley Tree
%% ------------
%% Diagonalization modules are packed in directory "NNH"
%% Add calculation of correlation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1;
clear all
more off

addpath("../NNH")
##Diag
Dens
Corr
Exact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonalization of Block Hamiltonian defined on the Cayley Tree
%% [Ref] T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 3
N = 10
param.t = ones(1,N);
param.mu = 0.2;
##param.f = @(s) 1;
param

## File Open
dirtag = [int2str(z) "-" int2str(N) "-sc"]
dirtag = ["Data/" dirtag "/"];
mkdir("Data");
mkdir(dirtag);
ftag = ["_mu-" num2str(param.mu)];
ftag = [ftag ".txt"];

fname = [dirtag "dens_path" ftag];
fid = fopen(fname, "w");

fname = [dirtag "Corr_path" ftag];
fid2 = fopen(fname, "w");
## --------------

##RCell = AllSites(z,N);
RCell{1} = [];
for j = 1:N
  RCell{j+1} = [RCell{j} 1];
endfor

exrho = Rho_exact(z, param.mu);
exbonde = E_exact(z, param.mu)*2/z;

for a = 0:4
  a
  tic;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
  fprintf(fid, "# a=%d\n", a);
  [Psi_All E_All] = Diag_All(z, N, param);
  
  for j = 1:prod(size(RCell))
    rho = Dens_GS(N, Psi_All, E_All, RCell{j});
    fprintf(fid, "%d %.15e\n", j-1, rho, exrho);
  endfor
  fprintf(fid, "\n\n");
  
  for j = 2: prod(size(RCell))
    C = Corr_GS(N, Psi_All, E_All, RCell{1}, RCell{j});
    fprintf(fid2, "%d %.15e\n", j-1, C);
  endfor
  fprintf(fid2, "\n\n");

  toc;
endfor

fclose("all");
clear fid fid2;
