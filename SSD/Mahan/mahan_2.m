1;
clear all;

function rho = DOS(z, E)
  A2 = 4*(z-1);
  rho = z*sqrt(A2-E*E)/(2*pi*(z*z-E*E));
endfunction

function rho = rhoG(z, mu)
  A = 2*sqrt(z-1);
  f = @(E) DOS(z, E);
  rho = quad(f, -A, mu);
endfunction

function rho = EG(z, mu)
  A = 2*sqrt(z-1);
  f = @(E) E*DOS(z, E);
  rho = quad(f, -A, mu);
endfunction

z = 5;
rhoG_ = @(mu) rhoG(z, mu);
EG_ = @(mu) EG(z, mu);

A = 2*sqrt(z-1);
A_ = A - 1e-12;
mu = [-A_:0.1:A_];
rho = arrayfun(rhoG_, mu);
E = arrayfun(EG_, mu);

fname = ["EG_" "z-" int2str(z) "_integral.txt"]
fid = fopen(fname, "w");
fprintf(fid, "%.15e %.15e %.15e\n", [mu; rho; E]);
fclose("all");
clear fid;
