function y = NBranch (z, s)
  if s<0
    y = 0;
  elseif s==0
    y = z;
  else
    y = z-1;
  endif
endfunction
