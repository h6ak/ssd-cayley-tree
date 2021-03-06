1;
%% Blocked Hamiltonian

function H = H_A(z, N, param)
  t = param.t;
  mu = param.mu;

  %% diagonal
  x = [0:N]';
  x = arrayfun(param.f, x);
  x = -mu * x;

  %% off-diagonal
  y = [0:(N-1)]' + 0.5;
  y = arrayfun(param.f, y);
  temp = y(1);
  y = -t*sqrt(z-1)*y;
  y(1) = -t*sqrt(z)*temp;
  y_ = [y;0];
  y = [0;y];
  
  H = spdiags([y_ x y], [-1:1], N+1, N+1);
endfunction

function H = H_B1(z, N, param)
  t = param.t;
  mu = param.mu;

  %% diagonal
  x = [1:N]';
  x = arrayfun(param.f, x);
  x = -mu * x;

  %% off-diagonal
  y = [1:(N-1)]' + 0.5;
  y = arrayfun(param.f, y);
  y = -t*sqrt(z-1)*y;
  y_ = [y;0];
  y = [0;y];
  
  H = spdiags([y_ x y], [-1:1], N, N);
endfunction

function H = H_B2(z, N, l, param)
  t = param.t;
  mu = param.mu;

  %% diagonal
  x = [(l+1):N]';
  x = arrayfun(param.f, x);
  x = -mu * x;

  %% off-diagonal
  y = [(l+1):(N-1)]' + 0.5;
  y = arrayfun(param.f, y);
  y = -t*sqrt(z-1)*y;
  y_ = [y;0];
  y = [0;y];
  
  H = spdiags([y_ x y], [-1:1], N-l, N-l);
endfunction
