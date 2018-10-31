1;
clear all
more off

Diag
Dens
Corr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonalization of Block Hamiltonian defined on the Cayley Tree
%% [Ref] T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 5
N = 10
a = 2
param.t = 1;
param.mu = 0;
param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
param

##RCell = AllSites(z,N);
RCell{1} = [];
for j = 1:N
  RCell{j+1} = [RCell{j} 1];
endfor

## File Open
dirtag = ["TestData2/"];
mkdir(dirtag);
ftag = [int2str(z) "-" int2str(N) "-sc"];
ftag = [ftag ".txt"];

fname = [dirtag ftag];
fid = fopen(fname, "w");

#---------------
A = abs(2*param.t*sqrt(z-1)) - 1e-5;
mu = [-A:0.1:A];
size(mu)
count = 0;
fprintf(fid, "# param.mu, C, rho1, rho2, EG\n\n");
tic;
for param.mu = mu
  [Psi_All E_All] = Diag_All(z, N, param);

  rho1 = Dens_GS(N, Psi_All, E_All, RCell{1});
  rho2 = Dens_GS(N, Psi_All, E_All, RCell{2});
  C = Corr_GS(N, Psi_All, E_All, RCell{1}, RCell{2});
  EG = -param.t*C -param.mu*(rho1 + rho2)/2;
  fprintf(fid, "%.15e %.15e %.15e %.15e %.15e\n", param.mu, C, rho1, rho2, EG);
  
  count += 1;
  if mod(count,10)==0
    count
    toc;
    tic;
  endif
endfor

fclose("all");
clear fid;
