%% R : coordinates
1;
Diag;

function rho = Dens_GS(N, Psi_All, E_All, R)
  rho = 0;
  rho += Dens_A(N, Psi_All{1}, E_All{1}, R);
  rho += Dens_B1(N, Psi_All{2}, E_All{2}, R);
  for l=1:(N-1)
    rho += Dens_B2(N, Psi_All{2+l}, E_All{2+l}, l, R);
  endfor
endfunction

function rho = Dens_A(N, Psi, E, R)
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  rho = 0;
  s = min( prod(size(R)), N);
  for j = idx
    rho += abs( Psi(s+1,j) )^2;
  endfor
endfunction

function rho = Dens_B1(N, PsiCell, E, R)
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

function rho = Dens_B2(N, PsiCell, E, l, R)
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
