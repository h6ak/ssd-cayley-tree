1;

function y = MyAtan(x)
  if x < 0
    y = atan(x) + pi;
  else
    y = atan(x);
  endif
endfunction

function rho = Rho_exact(z, mu)
  A = 2*sqrt(z-1);
  theta = acos(-mu/A);

  rho = z*theta - (z-2)*MyAtan( z*tan(theta)/(z-2) );
  rho = rho/(2*pi);
endfunction

function E = E_exact(z, mu)
  A2 = 4*(z-1);
  E = sqrt(A2-mu.^2) - (z-2)*atan( sqrt(A2-mu.^2)/(z-2) );
  E = -z*E/(2*pi);
endfunction

