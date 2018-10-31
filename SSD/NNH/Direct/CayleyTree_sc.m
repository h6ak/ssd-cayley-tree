%% site-centered tree
%% central site with shells
%% [1] : the 0-th shell 
function [Sites Bonds] = CayleyTree_sc(z, m)
  %% Sites
  Sites = cell(0);
  Sites{1} = [1];
  Sites{2} = [2: z+1];
  for j = 3:(m+1)
    x_min = Sites{j-1}(end);
    x_size = prod(size(Sites{j-1}));
    Sites{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor
  
  %% Make a list of Nearest-Neighbor Sites (toward the center)
  %% of the sites on the (j-1)-th shell
  NN = cell(0);
  NN{1} = ones(1,z);
  temp = ones(1,z-1);
  for j=2:m
    NN{j} = kron(Sites{j}, temp);
  endfor

  %% Bonds
  Bonds = cell(0);
  for j=1:m
    Bonds{j} = [NN{j};Sites{j+1}];
  endfor
endfunction
