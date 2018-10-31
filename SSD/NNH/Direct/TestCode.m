clear all
more off

z = 3
M = 6
[Sites, Bonds] = DefCayleyTree(z, M, "sc");

param.mu = 0;
param.t = ones(1,M);
param.f = @(s) 1;
#param.f = @(s) sin()^2;

fname = ["Data/energy_" num2str(z) "-"  num2str(M) ".txt"];
fid = fopen(fname, "w");

for dlt = 0:0.1:0.9

  param.mu = zeros(1,M+1);
  param.dlt = dlt;
  for j=1:prod(size(param.mu))
    param.mu(j) += param.dlt * (-1)^j;
  endfor
  param
  
  H = HoppingModel(Sites, Bonds, param);
  E = eig(H);
  
  fprintf(fid, "# z = %d, M = %d\n", z, M);
  for j = 1:prod(size(E))
    fprintf(fid, "%d %.15f\n", j, E(j));
  endfor
  fprintf(fid, "\n\n");
endfor

fclose("all");
