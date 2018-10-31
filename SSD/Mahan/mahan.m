1;
clear all;

function y = MyAtan(x)
  if x < 0
    y = atan(x) + pi;
  else
    y = atan(x);
  endif
endfunction

function rho = rhoG(z, mu)
  A = 2*sqrt(z-1);
  theta = acos(-mu/A);

  rho = z*theta - (z-2)*MyAtan( z*tan(theta)/(z-2) );
  rho = rho/(2*pi);
endfunction

function E = EG(z, mu)
  A2 = 4*(z-1);
  E = sqrt(A2-mu.^2) - (z-2)*atan( sqrt(A2-mu.^2)/(z-2) );
  E = -z*E/(2*pi);
endfunction

z = 6;
rhoG_ = @(mu) rhoG(z, mu);
EG_ = @(mu) EG(z, mu);

A = 2*sqrt(z-1);
A_ = A - 1e-12;
mu = [-A_:0.01:A_];
rho = arrayfun(rhoG_, mu);
E = arrayfun(EG_, mu);
fname = ["EG_" "z-" int2str(z) ".txt"]
fid = fopen(fname, "w");
fprintf(fid, "%.15e %.15e %.15e\n", [mu; rho; E]);
fclose("all");
clear fid;
