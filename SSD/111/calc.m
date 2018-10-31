function calc(z, N, a)
  param.mu = 0;
  param.mu2 = 0;
  param.f = @(s) sin( pi*(N-s+0.5)/(2*N+1) )^a;

  ## File Open
  dirtag = ["Data/z" int2str(z) "/"];
  mkdir(dirtag);
  ftag = [int2str(N) "-" num2str(a) "-"];
  fname = [dirtag ftag "n.txt"];
  fid_n = fopen(fname, "w");
  fname = [dirtag ftag "e.txt"];
  fid_e = fopen(fname, "w");

  mu2_max = 1.999;
  mu2_num = 101;
  mu2_step = 2*mu2_max/mu2_num;
  mu2_list = [-mu2_max : mu2_step : mu2_max];

  fprintf(fid_n, "# mu2 n0_GS n_exact n_diff\n");
  fprintf(fid_e, "# mu2 e0_GS e_exact e_diff\n");
  
  for mu2 = mu2_list
    param.mu2 = mu2;
    param.mu = mu2*sqrt(z-1);
    H = make_H(z, N, param);

    [V E] = eig(H);
    E = diag(E);
    idx = find(E<0);

    n0_GS = 0;
    e0_GS = 0;
    for j=1:prod(size(idx))
      n0_GS += V(1,j)*V(1,j);
      e0_GS += -2*V(1,j)*V(2,j)/sqrt(z);
    endfor
    
    n_exact = Rho_exact(z, param.mu);
    n_diff = abs(n0_GS - n_exact);
    fprintf(fid_n, "%.5f %.15f %.15f %.15f\n", mu2, n0_GS, n_exact, n_diff);
    
    e_exact = E_exact(z, param.mu)*2/z;
    e_diff = abs(e0_GS - e_exact);
    fprintf(fid_e, "%.5f %.15f %.15f %.15f\n", mu2, e0_GS, e_exact, e_diff);
  endfor

  fclose("all");
  clear fid_n, fid_e;
endfunction
