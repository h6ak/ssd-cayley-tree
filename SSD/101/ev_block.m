1;
clear functions
Diag

[Psi_All E_All] = Diag_All(z, N, param);

E_Array = [];
n_cell = prod(size(E_All));
for j=1:n_cell
  if j==1
    E_Array = [E_Array;E_All{j}];
  elseif j==2
    E_Array = [E_Array; repmat(E_All{j}, z-1, 1)];
  else
    n_rep = (z-2)*z*(z-1)^(j-3);
    E_Array = [E_Array; repmat(E_All{j}, n_rep, 1)];
  endif
endfor
E_Array = sort(E_Array);

