%% R : coordinates
1;
Diag;

function C = Corr_GS(N, Psi_All, E_All, R1, R2)
  C = 0;
  C += Corr_A(N, Psi_All{1}, E_All{1}, R1, R2);
  C += Corr_B1(N, Psi_All{2}, E_All{2}, R1, R2);
  for l=1:N-1
    C += Corr_B2(N, Psi_All{2+l}, E_All{2+l}, l, R1, R2);
  endfor
endfunction

function C = Corr_A(N, Psi, E, R1, R2)
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  C = 0;
  s1 = min( prod(size(R1)), N);
  s2 = min( prod(size(R2)), N);
  for j = idx
    Z = [ Psi(s1+1,j) Psi(s2+1,j) ];
    X = real(Z);
    Y = imag(Z);
    C += 2*( prod(X) + prod(Y) );
  endfor
endfunction

function C = Corr_B1(N, PsiCell, E, R1, R2)
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  C = 0;
  s1 = min( prod(size(R1)), N);
  s2 = min( prod(size(R2)), N);
  if s1 == 0 || s2 == 0
    return
  endif
  nmax = size(PsiCell)(1);
  r1 = R1(1);
  r2 = R2(1);
  for n = 1: nmax
    for j = idx
      Z = [ PsiCell{n,r1}(s1+1,j) PsiCell{n,r2}(s2+1,j) ];
      X = real(Z);
      Y = imag(Z);
      C += 2*( prod(X) + prod(Y) );
    endfor
  endfor
endfunction

function C = Corr_B2(N, PsiCell, E, l, R1, R2)
  idx = find(E<0);
  idx = reshape(idx, 1, []);

  C = 0;
  s1 = min( prod(size(R1)), N);
  s2 = min( prod(size(R2)), N);
  if s1 <= l || s2 <= l
    return
  endif
  nmax = size(PsiCell)(1);  
  r1 = R1(l+1);
  r2 = R2(l+1);
  for n = 1: nmax
    for j = idx
      Z = [ PsiCell{n,r1}(s1+1,j) PsiCell{n,r2}(s2+1,j) ];
      X = real(Z);
      Y = imag(Z);
      C += 2*( prod(X) + prod(Y) );
    endfor
  endfor
endfunction
