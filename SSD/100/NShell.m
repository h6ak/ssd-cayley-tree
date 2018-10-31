%% #Sites on the s-th shell
function y = NShell(z, s)
  if s == 0
    y = 1;
  elseif s > 0
    y = z*(z-1)^(s-1);
  else
    y = 0;
  endif
endfunction
