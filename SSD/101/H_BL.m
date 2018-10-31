1;
%%%%%%% Parameters %%%%%%%%
%% z, m can be array
function CheckNSite(z, m)
  sz = prod(size(z));
  sm = prod(size(m));
  z_ = reshape(z,[],1);
  z_ = repmat(z_, 1, sm)
  m_ = reshape(m,1,[]);
  m_ = repmat(m_, sz, 1)
  N_bc = 2.*( (z_-1).^m_ -1 )./(z_-2)
  N_sc = (z_.*(z_-1).^(m_-1) -2)./(z_-2)
endfunction


%%%%%% Hamiltonian for bond-centered Bethe lattice %%%%%%%
function [H, Lattice] = H_bc(z, m, t, mu, a)
  %% Lattice
  Lattice = cell(0);
  Lattice{1} = [1 2];
  for j = 2:m
    x_min = Lattice{j-1}(end);
    x_size = prod(size(Lattice{j-1}));
    Lattice{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor

  %% Hamiltonian for one-particle subspace
  sinf = @(x, p) sin( pi*x/(2*m) )^p;

  %% off-diagonal
  ones_ = ones(1, z-1);
  x = [1];
  y = [2];
  z = [(-t*sinf(m, a))];
  
  for j = 1:m-1
    x_ = kron( Lattice{j}, ones_ );
    y_ = Lattice{j+1};
    z_ = repmat([(-t*sinf(m-j, a))], size(y_));

    x = [x x_];
    y = [y y_];
    z = [z z_];
  endfor
  clear x_ y_ z_
  H = sparse([x y], [y x], [z z]);
  
  %% diagonal
  x = zeros(1,0);
  z = x;
  for j = 1:m
    x_ = Lattice{j};
    z_ = repmat([(-mu*sinf(m-j+0.5, a))], size(x_));
    
    x = [x x_];
    z = [z z_];
  endfor
  H += sparse(x, x, z);

endfunction

%%%%%% Hamiltonian for site-centered Bethe lattice %%%%%%%
function [H, Lattice] = H_sc(z, m, t, mu, a)
  %% Lattice
  Lattice = cell(0);
  Lattice{1} = [1];
  Lattice{2} = [2: z+1];
    
  for j = 3:m
    x_min = Lattice{j-1}(end);
    x_size = prod(size(Lattice{j-1}));
    Lattice{j} = [(x_min + 1) : (x_min + (z-1)*x_size)];
  endfor

  %% Hamiltonian for one-particle subspace
  sinf = @(x, p) sin( pi*x/(2*m-1) )^p;

  %% off-diagonal
  ones_ = ones(1, z-1);
  x = repmat([1],1,z);
  y = [2:z+1];
  z = repmat([(-t*sinf(m-1, a))], 1, z);
  for j = 2:m-1
    x_ = kron( Lattice{j}, ones_ );
    y_ = Lattice{j+1};
    z_ = repmat([(-t*sinf(m-j, a))], size(y_));

    x = [x x_];
    y = [y y_];
    z = [z z_];
  endfor
  clear x_ y_ z_
  H = sparse([x y], [y x], [z z]);

  %% diagonal
  x = zeros(1,0);
  z = x;
  for j = 1:m
    x_ = Lattice{j};
    z_ = repmat([(-mu*sinf(m-j+0.5, a))], size(x_));
    
    x = [x x_];
    z = [z z_];
  endfor
  H += sparse(x, x, z);

endfunction

%%E=eig(H);
