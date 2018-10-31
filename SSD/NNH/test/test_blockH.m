%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1;
clear all
more off

addpath("..")
Blocked_H

z = 3
N = 5
param.t = ones(1,N);
param.mu = 0.2;
param.f = @(s) 1;
param

H1 = H_A(z, N, param);
H2 = H_B1(z, N, param);
H3 = H_B2(z, N, 1, param);
