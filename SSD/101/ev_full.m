1;
clear functions
H_BL

[H Lattice] = H_sc(z, m, t, mu, a);
[V E] = eig(H);
E = diag(E);

