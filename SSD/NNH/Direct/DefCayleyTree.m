%% z: coordination number, m: the number of shells
%%
%% Sites (cell): Sites{j} denotes sites on the (j-1) th shell
%% Bonds (cell):
%% Bonds{j} is a 2 rows matrix,
%% where each column denotes edges of a bond
function [Sites Bonds] = DefCayleyTree(z, m, opt)
  switch opt
    case 'sc'
      [Sites Bonds] = CayleyTree_sc(z, m);
    case 'bc'
      [Sites Bonds] = CayleyTree_bc(z, m);
    otherwise
      Sites = cell{0};
      Bonds = cell{0};
  endswitch
endfunction
