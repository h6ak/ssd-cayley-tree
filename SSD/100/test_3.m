%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate Bulk Quantities
%% Compare Results with Corresponding Exact Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1;
clear all
more off

Diag
Dens
Corr
Exact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonalization of Block Hamiltonian defined on the Cayley Tree
%% [Ref] T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 3
N = 50
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
dirtag = ["TestData3/"];
mkdir(dirtag);
ftag = [int2str(z) "-" int2str(N) "-sc"];
ftag = [ftag ".txt"];

fname = [dirtag ftag];
fid = fopen(fname, "w");

#---------------
A = abs(2*param.t*sqrt(z-1)) - 1e-5;
mu = [-A:0.01:A];
size(mu)
count = 0;
fprintf(fid, "# param.mu, -tC, rho1, E_exact*2/z, rho_exact\n\n");
tic;
for param.mu = mu
  [Psi_All E_All] = Diag_All(z, N, param);
  
  C = Corr_GS(N, Psi_All, E_All, RCell{1}, RCell{2});
  C_exact = E_exact(z, param.mu) * 2/z;

  rho = Dens_GS(N, Psi_All, E_All, RCell{1});
  rho_exact = Rho_exact(z, param.mu);
  fprintf(fid, "%.15e %.15e %.15e %.15e %.15e\n", param.mu, (-param.t)*C, rho,
	  C_exact, rho_exact);
  
  count += 1;
  if mod(count,10)==0
    count
    toc;
    tic;
  endif
endfor

fclose("all");
clear fid;
