function H = make_H(z, N, param)
  f = param.f;
  mu2 = param.mu2;
  a = sqrt(z/(z-1));
  %% diagonal
  x = [0:N]';
  x = arrayfun(f, x);
  x = (- mu2)* x;

  %% off-diagonal
  y = [0:(N-1)]' + 0.5;
  y = arrayfun(f, y);
  y(1) = (-a) * y(1);
  y_ = [y;0];
  y = [0;y];
  
  H = spdiags([y_ x y], [-1:1], N+1, N+1);
endfunction
