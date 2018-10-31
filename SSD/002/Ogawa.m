1;
clear all
more off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% T.Ogawa, Prog.Theor.Phys.54, 1028 (1975) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% m: #layer except for root

function H = H_Ogawa(z, m, t, mu)
  H1_ = @(N) H1(N, z, t, mu);
  H2_ = @(N) H2(N, z, t, mu);
  sig = @(s) z*(z-1)^(s-1);

  H = H1_(m+1);
  for j = 1:z-1
    H = blkdiag(H, H2_(m));
  endfor

  for j = 1:m-1
    for k = 1:(z-2)*sig(j)
      H = blkdiag(H, H2_(m-j));
    endfor
  endfor
endfunction

function H = H1(N, z, t, mu)
  x = -mu*ones(N,1);
  y = -t*sqrt(z-1)*ones(N,1);
  H = spdiags([y x y], [-1:1], N, N);
  H(1,2) = -t*sqrt(z);
  H(2,1) = H(1,2);
endfunction

function H = H2(N, z, t, mu)
  x = -mu*ones(N,1);
  y = -t*sqrt(z-1)*ones(N,1);
  H = spdiags([y x y], [-1:1], N, N);
endfunction


z = 3
m = 7
t = 1
mu = 0
%%%%%
dirtag = ["Data/Ogawa/"];
mkdir(dirtag);
ftag = ["_" int2str(z) "-" int2str(m) "-sc.txt"];

fname = [dirtag "energy" ftag];
f_energy = fopen(fname, "w");
%%%%%

H = H_Ogawa(z, m-1, t, mu);
E = eig(H);

fid = f_energy;
temp = [1:prod(size(E))];
fprintf(fid, "%d %.15e\n", [temp; E']);
fprintf(fid, "\n\n");
clear temp

fclose(f_energy);

