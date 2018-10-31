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


%%%%%% Hamiltonian for Cayley Tree %%%%%%%
function H = H_Hopping(Sites, Bonds, opt, t, mu, a)
  switch opt
    case 'bc'
      NChain = 2*size(Sites)(2);
    case 'sc'
      NChain = 2*size(Sites)(2)-1;
    otherwise
      H=[];
      return
  endswitch


  %% sinusoidal function 
  sinf = @(x, p) sin( pi*(x-0.5)/NChain )^p;

  %% off-diagonal elements
  x = zeros(1,0);
  y = zeros(1,0);
  f = zeros(1,0);
  jmax = size(Bonds)(2);
  for j = 1:jmax
    k = jmax-j+1;

    x_ = Bonds{1,k};
    y_ = Bonds{2,k};
    f_ = -t * sinf(j+0.5, a);
    f_ = repmat(f_, size(x_));
    
    x = [x x_];
    y = [y y_];
    f = [f f_];
  endfor
  H = sparse([x y], [y x], [f f]);
  
  %% diagonal elements
  x = zeros(1,0);
  f = zeros(1,0);
  clear y y_
  jmax = size(Sites)(2);
  for j = 1:jmax
    k = jmax-j+1;

    x_ = Sites{k};
    f_ = -mu * sinf(j, a);
    f_ = repmat(f_, size(x_));

    x = [x x_];
    f = [f f_];
  endfor
  H += sparse(x, x, f);
endfunction


