1;
clear all
more off

##Diag
Dens
Corr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonalization of Block Hamiltonian defined on the Cayley Tree
%% [Ref] T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 3
N = 50
param.t = 1;
param.mu = 0.5;
##param.f = @(s) 1;
param

## File Open
dirtag = [int2str(z) "-" int2str(N) "-sc"]
dirtag = ["TestData/" dirtag "/"];
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

for a = 0:4
  a
  tic;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
  fprintf(fid, "# a=%d\n", a);
  [Psi_All E_All] = Diag_All(z, N, param);
  
  Sum_rho = 0;
  for j = 1:prod(size(RCell))
    rho = Dens_GS(N, Psi_All, E_All, RCell{j});
    
    Sum_rho += rho;
    fprintf(fid, "%d %.15e %.15e\n", j-1, rho, Sum_rho/j);
  endfor
  fprintf(fid, "\n\n");

  Sum_C = 0;
  for j = 1:( prod(size(RCell))-1 )
    C = Corr_GS(N, Psi_All, E_All, RCell{j}, RCell{j+1});
    Sum_C += C;
    fprintf(fid2, "%d %.15e %.15e\n", j-1, C, Sum_C/j);
  endfor
  fprintf(fid2, "\n\n");

  toc;
endfor

fclose("all");
clear fid fid2;
