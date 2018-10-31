function [E A] = EG_Mahan (z, t, mu)
  A = -2*t*sqrt(z-1);
  E = - (z/pi) *...
	( sqrt(A*A-mu*mu) - (z-2)*atan( sqrt(A*A-mu*mu)/(z-2) ) );
  A = abs(A);
endfunction
