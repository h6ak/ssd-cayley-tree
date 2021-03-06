1;
clear all;

z = 3;
m = 3;

R0 = 300/2;
R1 = R0-20;
CX = R0;
CY = CX;
FontSize = 4;

eps_out = "CayleyTree.eps"
copyfile("template2.txt", eps_out);
fid = fopen(eps_out, "a");
%%fprintf(fid, "%d %d circle\n", CX, CY);

Theta = cell(1, m);
SiteR = cell(1, m);
R = zeros(1, m);
NShell = @(s) z*(z-1)^(s-1);

%% Theta, R
N = NShell(m);
for l=1:NShell(m)
  Theta{m}(l) = 2*pi*l/N;
  R(m) = R1;
endfor

for j=m-1:-1:1
  N = NShell(j);
  for l=1:N
    id_start = 1+(z-1)*(l-1);
    id_end = id_start + z - 2;
    Theta{j}(l) = mean( Theta{j+1}(id_start:id_end) );
    R(j) = R1*j/m;
  endfor
endfor

Theta0 = pi/2 - Theta{1}(1);
%%Theta0 = - Theta{1}(1);
for j=1:m
  Theta{j} = Theta{j} + Theta0;
endfor

%%SiteR
for j=1:z
  SiteR{1}{j} = int2str(j);
endfor
for j=2:m
  N = NShell(j);
  for l=1:N
    SiteR{j}{l} = [SiteR{j-1}{ceil(l/(z-1))} int2str(mod(l+1,z-1)+1)];
  endfor
endfor

%% Lines
X0 = CX;
Y0 = CY;
for k=1:z
  X1 = R(1)*cos(Theta{1}(k)) + CX;
  Y1 = R(1)*sin(Theta{1}(k)) + CY;
  fprintf(fid, "%d %d %d %d p\n", X0, Y0, X1, Y1);
endfor

for j=1:m-1
  for l=1:size(Theta{j})(2)
    X0 = R(j)*cos(Theta{j}(l)) + CX;
    Y0 = R(j)*sin(Theta{j}(l)) + CY;
    id_start = 1+(z-1)*(l-1);
    for k=1:(z-1)
      X1 = R(j+1)*cos(Theta{j+1}(id_start + k -1)) + CX;
      Y1 = R(j+1)*sin(Theta{j+1}(id_start + k -1)) + CY;
      fprintf(fid, "%d %d %d %d p\n", X0, Y0, X1, Y1);
    endfor
  endfor
endfor

%% Ponits
X = CX;
Y = CY;
fprintf(fid, "%d %d c\n", X, Y);
fprintf(fid, "%d %d moveto\n", X-FontSize, Y-FontSize);
fprintf(fid, "(0) show\n\n");

for j=1:m
  for l=1:size(Theta{j})(2)
    X = R(j)*cos(Theta{j}(l)) + CX;
    Y = R(j)*sin(Theta{j}(l)) + CY;
    fprintf(fid, "%d %d c\n", X, Y);

    fprintf(fid, "%d %d moveto\n", X-j*FontSize, Y-FontSize);
    fprintf(fid, "(%s) show\n\n", SiteR{j}{l});
  endfor
endfor

fprintf(fid, "\n%%%%EOF");
fclose("all");
clear fid
