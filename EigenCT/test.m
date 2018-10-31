function y = test(n, th)
  y1 = cos((n+2)*th) * sin( (n-1)*th ) ...
       - sin( (n+1)*th ) * cos( n*th );
  y2 = -sin(2*th)*cos(th);
  y = y1 -y2;
endfunction
