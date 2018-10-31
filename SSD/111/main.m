%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uniform NN hopping model on Cayley Tree
%% ------------
%% Diagonalization modules are packed in directory "NNH"
%% Add calculation of correlation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1;
clear all
more off
format long E

addpath("../NNH")
Exact

z = [3 4 5 10 50 100];
N = [10 20 50];
##a = [0:5];
a = 10;

for z_ = z
  for N_ = N
    for a_ = a
      printf("[%d, %d, %d]\n", z_, N_, a_);
      calc(z_, N_, a_);
    endfor
  endfor
endfor
