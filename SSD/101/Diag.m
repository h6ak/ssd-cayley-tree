1;

Blocked_H

%% All
function [Psi_All E_All] = Diag_All(z, N, param)
  Psi_All = cell(1,N+1);
  E_All = cell(1,N+1);

  [Psi_All{1} E_All{1}] = Diag_A(z, N, param);
  [Psi_All{2} E_All{2}] = Diag_B1(z, N, param);
  for l=1:(N-1)
    [Psi_All{2+l} E_All{2+l}] = Diag_B2(z, N, l, param);
  endfor
endfunction

%% (R_0, 0)
function [Psi E] = Diag_A(z, N, param)
  H = H_A(z, N, param);
  [Psi_ E] = eig(H);
  E = diag(E);

  nf = @(s) 1/sqrt( NShell(z, s) );
  
  Psi = zeros(N+1, size(Psi_)(2));
  for s = 0:N
    Psi(s+1,:) = nf(s)*Psi_(s+1,:);
  endfor
endfunction

%% (R_0, n)
function [PsiCell E] = Diag_B1(z, N, param)
  H = H_B1(z, N, param);
  [Psi_ E] = eig(H);
  E = diag(E);

  nmax = z-1;
  rmax = z;
  PsiCell = cell(nmax, rmax);
  nf = @(s) 1/sqrt( NShell(z, s) );

  for r = 1:rmax
    for n = 1:nmax
      PsiCell{n, r} = zeros(N+1, size(Psi_)(2));
      for s = 1:N
	PsiCell{n, r}(s+1,:) = nf(s)*exp(2*pi*1i*n*r/z)*Psi_(s,:);
      endfor
    endfor
  endfor
endfunction

%% (R_l, n)
function [PsiCell E] = Diag_B2(z, N, l, param)
  H = H_B2(z, N, l, param);
  [Psi_ E] = eig(H);
  E = diag(E);

  nmax = z-2;
  rmax = z-1;
  PsiCell = cell(nmax, rmax);
  nf = @(s) sqrt( NShell(z, l) / NShell(z, l+s) );

  for r = 1:rmax
    for n = 1:nmax
      PsiCell{n, r} = zeros(N+1, size(Psi_)(2));
      for s = 1:N-l
	PsiCell{n, r}(l+s+1,:) = nf(s)*exp(2*pi*1i*n*r/(z-1))*Psi_(s,:);
      endfor
    endfor
  endfor
endfunction
