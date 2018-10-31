#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
typedef unsigned int uint;

const double pi = acos(-1e0);

int sitesInShell(uint z, uint s){
  if(s==0) return 1;
  else {
    return z*pow(z-1, s-1);
  }
}

void cayleyTree 
(ostream &os, uint z, uint m, uint s,
 double x0, double y0, double r, double t, double dr){
  if(z<3 || m==0) return;

  double x, y;
  x = x0 + r*cos(t);
  y = y0 + r*sin(t);
  os << s << " " << x << " " << y << endl;
  //os << s << " " << r << " " << t << endl;

  if(s==0){
    r += dr;
    t = pi/2e0;
    double dt = 2e0*pi/z;
    for(int i=0; i<z; i++){
      cayleyTree(os, z, m, s+1, x0, y0, r, t, dr);
      t += dt;
    }
  }
  else if(s<m){
    r += dr;
    t -= pi*(z-1) / sitesInShell(z,s+1);
    double dt = 2e0*pi / sitesInShell(z,s+1);
    for(int i=0; i<z-1; i++){
      cayleyTree(os, z, m, s+1, x0, y0, r, t, dr);
      t += dt;
    }
  }
}

int main (void){
  uint z = 4, m = 3;
  ofstream ofs("data.txt");
  cayleyTree(ofs, z, m, 0, 0, 0, 0, 0, 1);
  return 0;
}
