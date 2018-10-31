%% R : coordinates
1;
Diag;

function rho = Dens_GS(z, N, param, R)
  rho = 0;
  rho += Dens_A(z, N, param, R);
  rho += Dens_B1(z, N, param, R);
  for l=1:(N-1)
    rho += Dens_B2(z, N, l, param, R);
  endfor
endfunction

function rho = Dens_A(z, N, param, R)
  [Psi E] = Diag_A(z, N, param);
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  rho = 0;
  s = min( prod(size(R)), N);
  for j = idx
    rho += abs( Psi(s+1,j) )^2;
  endfor
endfunction

function rho = Dens_B1(z, N, param, R)
  [PsiCell E] = Diag_B1(z, N, param);
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  rho = 0;
  s = min( prod(size(R)), N);
  if s==0
    return
  endif
  nmax = size(PsiCell)(1);
  r = R(1);
  for n = 1: nmax
    for j = idx
      rho += abs( PsiCell{n,r}(s+1,j) )^2;
    endfor
  endfor
endfunction

function rho = Dens_B2(z, N, l, param, R)
  [PsiCell E] = Diag_B2(z, N, l, param);
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  rho = 0;
  s = min( prod(size(R)), N);
  if s<=l
    return
  endif
  nmax = size(PsiCell)(1);
  r = R(l+1);
  for n = 1: nmax
    for j = idx
      rho += abs( PsiCell{n,r}(s+1,j) )^2;
    endfor
  endfor
endfunction
