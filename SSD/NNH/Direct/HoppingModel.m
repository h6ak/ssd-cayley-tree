% Make Hamiltonian for
% Nearest-Neighbor Hopping Model on Cayley Tree
% We assume local energy scale is isotropic.
%
% ### needed variables ###
% Sites (cell)
% Bonds (cell)
% param (structure):
% -- t (1D array): hopping amplitude between j-th & (j+1)-th shells
% -- mu (1D array): chemical potential on each shell
% -- f (function bandle): modulation function
% ####
function H = HoppingModel(Sites, Bonds, param)
  NSites = 0;
  for j=1:prod(size(Sites))
    NSites += prod(size(Sites{j}));
  endfor
  H = sparse(NSites, NSites);

  % The number of shells
  M = prod(size(Sites)) - 1;

  % diagonal elements
  for j = 0:M
    for k = 1:prod(size(Sites{j+1}))
      x = Sites{j+1}(k);
      H(x, x) = -param.mu(j+1) * param.f(j);
    endfor
  endfor

  % off-diagonal elements
  for j=1:M
    for k = 1:size(Bonds{j})(2)
      x = Bonds{j}(1, k);
      y = Bonds{j}(2, k);
      H(x, y) = -param.t(j) * param.f(j-0.5);
      H(y, x) = H(x, y);
    endfor
  endfor
endfunction
