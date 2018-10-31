%% bond-centered tree
%% central bond with shells
%% [1, 2] : the 0-th shell 
function [Sites Bonds] = CayleyTree_bc(z, m)
  %% Sites
  Sites = cell(0);
  Sites{1} = [1 2];
  for j = 2:(m+1)
    x_min = Sites{j-1}(end);
    x_size = prod(size(Sites{j-1}));
    Sites{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor
  
  %% Make a list of Nearest-Neighbor Sites (toward the center)
  %% of the sites on the (j-1)-th shell
  NN = cell(0);
  temp = ones(1,z-1);
  for j=1:m
    NN{j} = kron(Sites{j}, temp);
  endfor

  %% Bonds
  Bonds = cell(0);
  Bonds{1} = [1;2];
  for j=2:(m+1)
    Bonds{j} = [NN{j-1};Sites{j}];
  endfor
endfunction
