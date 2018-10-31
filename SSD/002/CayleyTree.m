1;
function [Sites Bonds] = DefCayleyTree(z, m, opt)
  switch opt
    case 'bc'
      [Sites Bonds] = DefCayleyTree_bc(z, m);
    case 'sc'
      [Sites Bonds] = DefCayleyTree_sc(z, m);
    otherwise
      Sites = cell{0};
      Bonds = cell{0};
  endswitch
endfunction


%% bond-centered 
%% z legs, m layers
function [Sites Bonds] = DefCayleyTree_bc(z, m)
  %% Sites
  Sites = cell(0);
  Sites{1} = [1 2];
  for j = 2:m
    x_min = Sites{j-1}(end);
    x_size = prod(size(Sites{j-1}));
    Sites{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor

  %% Bonds
  Bonds = cell(0);
  Bonds{1, 1} = [1];
  Bonds{2, 1} = [2];
  ones_ = ones(1, z-1);
  for j = 2:m
    Bonds{1, j} = kron( Sites{j-1}, ones_ );
    Bonds{2, j} = Sites{j};
  endfor
endfunction

%% site-centered 
%% z legs, m layers
function [Sites Bonds] = DefCayleyTree_sc(z, m)
  %% Sites
  Sites = cell(0);
  Sites{1} = [1];
  Sites{2} = [2: z+1];
  for j = 3:m
    x_min = Sites{j-1}(end);
    x_size = prod(size(Sites{j-1}));
    Sites{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor

  %% Bonds
  Bonds = cell(0);
  Bonds{1, 1} = repmat([1],1,z);
  Bonds{2, 1} = [2:z+1];
  ones_ = ones(1, z-1);
  for j = 2:m-1
    Bonds{1, j} = kron( Sites{j}, ones_ );
    Bonds{2, j} = Sites{j+1};
  endfor
endfunction
