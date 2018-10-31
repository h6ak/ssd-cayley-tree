1;
clear all
more off

##Diag
Dens
Corr
Exact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diagonalization of Block Hamiltonian defined on the Cayley Tree
%% [Ref] T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = input("z? :");
N = input("N? :");
param.t = 1;
param.mu = input("mu? :");
##param.f = @(s) 1;
param

## File Open
dirtag = [int2str(z) "-" int2str(N) "-sc"]
dirtag = ["TestData4/" dirtag "/"];
mkdir(dirtag);
ftag = ["_mu-" num2str(param.mu)];
ftag = [ftag "-v4.txt"];

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

exact_dens = Rho_exact(z, param.mu);
exact_bonde = (2/z)*E_exact(z, param.mu);

##a_list = [1.5:0.1:2.5];
##a_list = [0 1 a_list 3 4 5 6]
a_list = [0:6]

for FileID = [fid fid2 stdout]
  for j = 1:prod(size(a_list))
    fprintf(FileID, "#[%d] a=%g ", j-1, a_list(j));
    fprintf(FileID, "\n\n");
  endfor
endfor

for a = a_list
  a
  tic;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;
  fprintf(fid, "#a=%d\n", a);
  fprintf(fid2, "#a=%d\n", a);
  [Psi_All E_All] = Diag_All(z, N, param);
  
  for j = 1:prod(size(RCell))
    dens = Dens_GS(N, Psi_All, E_All, RCell{j});
    fprintf(fid, "%d %.15e %.15e %.15e\n",
	    j-1, dens, exact_dens, dens-exact_dens);
  endfor
  fprintf(fid, "\n\n");
  
  for j = 1:( prod(size(RCell))-1 )
    bonde = (-param.t)*Corr_GS(N, Psi_All, E_All, RCell{j}, RCell{j+1});
    fprintf(fid2, "%d %.15e %.15e %.15e\n",
	    j-1, bonde, exact_bonde, bonde-exact_bonde);
  endfor
  fprintf(fid2, "\n\n");

  toc;
endfor

fclose("all");
clear fid fid2;
