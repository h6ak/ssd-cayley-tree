1;
%% E: diagonal matrix
function rho = Density_GS(x, V, E)
  idx = find(E<0,1,'last');
  V_ = V(x,1:idx);
  rho = dot(V_, V_);
endfunction

%% E: diagonal matrix
function C = Corr_GS(x1, x2, V, E)
  idx = find(E<0,1,'last');
  V1 = V(x1,1:idx);
  V2 = V(x2,1:idx);
  C = dot(V1, V2);
endfunction
