## Eigenvalue Problem for Hamiltonian Matrix A
##
clear all
more off
norm_V

## parameters
z = input("Input: z = ");
N = input("Input: N = ");
z; N;
a = sqrt(z/(z-1))

## make Hamiltonian matrix A
## Be careful that A is a (N+1) x (N+1) matrix
y = ones(N, 1);
y(1) = a;
y1 = [y; 0];
y2 = [0; y];
A = spdiags([y1 y2], [-1 1], N+1, N+1);

## numerical solution
[V E] = eig(A);
E = diag(E);

## analytical
theta = acos(E/2);
func = @(x, n) sin( (n+1)*x ) + (1-a*a)*sin( (n-1)*x );
applied = arrayfun(@(x) func(x, N+1), theta);   # should be a zero vector

theta_ = theta(1);
V_ = zeros(N+1);
V_(1,:) = arrayfun(@(x) a*sin(x), theta);
V_(2,:) = arrayfun(@(x) sin(2*x), theta);
for j = 3:N+1
  V_(j,:) = arrayfun(@(x) func(x, j-1), theta);
endfor

F = diag( V_.'*V_ );
for j = 1:N+1
  V_(:,j) = V_(:,j)/sqrt(F(j));
endfor

##V.'*V_
F2 = arrayfun(@(x) norm_V_(a, N+1, x), theta);
F-F2
				    
