1;

function y = norm_V_ (a, n, th)
  y = n*(1-a*a)*cos(2*th) + (n-1)*(1-a*a)^2 /2 + (n+1)/2;
endfunction

## OK ##
function y = norm_V_0 (a, n, th)
  y = (a*sin(th))^2 + sin(2*th)^2;
  for j=1:n-1
    tmp = sin( (j+2)*th ) + (1-a^2)*sin(j*th);
    y += tmp*tmp;
  endfor
endfunction

## OK ##
function y = norm_V_1 (a, n, th)
  y1 = (a*a - 1)*sin(th)^2 + (n+1)/2;
  y2 = (n-1)*cos(2*th) ...
       - cos( (n+2)*th ) * sin( (n-1)*th ) /(2*sin(th));
  y3 = (n-1)/2 ...
       - cos( n*th ) * sin( (n-1)*th ) / (2*sin(th));
  y = y1 + y2*(1-a*a) + y3*(1-a*a)^2;
endfunction

function y = norm_V_2 (a, n, th)
  y1 = (a*a - 1)*sin(th)^2 + (n+1)/2;
  y2 = (n-1)*cos(2*th) + cos(th)**2;
  y3 = (n-1)/2;
  y = y1 + y2*(1-a*a) + y3*(1-a*a)^2;
endfunction
