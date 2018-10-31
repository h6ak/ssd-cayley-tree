%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Various Deformations for uniform NN hopping model on Cayley Tree
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
N = 20
param.t = 1;
param.mu = 1;
##param.f = @(s) 1;
param

fcell = cell(3,2);
fcell{1,1} = @(s) 1;
fcell{1,2} = "OBC";
fcell{2,1} = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^2;
fcell{2,2} = "SSD";
fcell{3,1} = @(s) 2^(-s);
fcell{3,2} = "exp";

## File Open
dirtag = [int2str(z) "-" int2str(N) "-sc"]
dirtag = ["Data/" dirtag "/"];
mkdir(dirtag);
ftag = ["_mu-" num2str(param.mu)];
ftag = [ftag ".txt"];

fname = [dirtag "dens_path" ftag];
fid = fopen(fname, "w");

fname = [dirtag "BondCorr_path" ftag];
fid2 = fopen(fname, "w");
## --------------

##RCell = AllSites(z,N);
RCell{1} = [];
for j = 1:N
  RCell{j+1} = [RCell{j} 1];
endfor

exrho = Rho_exact(z, param.mu);
exbonde = E_exact(z, param.mu)*2/z;

[amax tmp]= size(fcell);
for a = 1:amax
  a
  tic;
  param.f = fcell{a,1};
  fprintf(fid, "# %s\n", fcell{a,2});
  [Psi_All E_All] = Diag_All(z, N, param);
  
  for j = 1:prod(size(RCell))
    rho = Dens_GS(N, Psi_All, E_All, RCell{j});
    fprintf(fid, "%d %.15e %.15e\n", j-1, rho, exrho);
  endfor
  fprintf(fid, "\n\n");
  
  for j = 1:( prod(size(RCell))-1 )
    bonde = (-param.t)*Corr_GS(N, Psi_All, E_All, RCell{j}, RCell{j+1});
    fprintf(fid2, "%d %.15e %.15e\n", j-1, bonde, exbonde);
  endfor
  fprintf(fid2, "\n\n");

  toc;
endfor

fclose("all");
clear fid fid2;
