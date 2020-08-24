#include<iostream>
//#include <vector>
#include <time.h>
#include "rings.h"  // Long,
#include "combinatorics.h"   //
using namespace std;
const int verbose = 2;
/*
Solve {X(z)}^{k+l} = {Z(z)}^l, Z(z)=(X(z)-1)/z  for X = c0 + c1 z + c2 z^2 + ...
 */
FormalPowerSeries<Rational<Long> > CatalanKernel(unsigned int k, unsigned int l) {
    if (verbose > 1) cout << "Solve x^" << k+l << " = ((x-1)/z)^" << l << " for x\n";
    FormalPowerSeries<Rational<Long> > X, Z;
    unsigned int maxD = X.getMaxD();
    X[0] = 1, Z[0] = X[1] = 1; Z[1] = 0;
    for (unsigned int N=1; N<maxD; N++) {
        Z[N] = X[N+1] = (polynomialPowerCoef(X,  N, k+l) - // [z^N] X(z)^{k+l}
                         polynomialPowerCoef(Z, N, l))/l;     // [z^N] Z(z)^l
        Z[N+1] = 0;
    }
    if (verbose > 1) cout << "K(" << k << "," << l << ")(z)\t= " << X << "\n";
    return X;
}

// Y(z) = X(-z);
template<typename T>
FormalPowerSeries<T> Fourier(FormalPowerSeries<T> X, unsigned int b) {
    FormalPowerSeries<T> Y;
    for (unsigned int i=0; i<X.readCoef().size(); i++) Y[i] = X[i] * (1-2*(int(i/b)%2));
    return Y;
}

// Y(z) = X(w^[a/b] z)
template<typename T>
FormalPowerSeries<Cyclotomic<T> > Fourier(FormalPowerSeries<T> X, unsigned int a, unsigned int b, unsigned int base) {
    FormalPowerSeries<Cyclotomic<T> > Y(base);
    for (unsigned int i=0; i<X.readCoef().size(); i++) Y[i][int(a*(i/b))%base] = X[i];
    return Y;
}

// Y(z) = X(z)
template<typename T>
FormalPowerSeries<T> FourierInv(FormalPowerSeries<Cyclotomic<T> > X) {
    FormalPowerSeries<T> Y;
    for (unsigned int i=0; i<X.readCoef().size(); i++) Y[i] = X[i][0];
    return Y;
}

FormalPowerSeries<Rational<Long> > ExtendedCatalan(FormalPowerSeries<Rational<Long> > X, unsigned int t, unsigned int m) {
    FormalPowerSeries<Rational<Long> > One, D, Yp, Yn;
    One[0] = 1;
    D[0] = -1; D[1] = 1; // D(z) = z-1

    FormalPowerSeries<Cyclotomic<Rational<Long> > > Yu(m), Yd(m);
    //Yu.setMaxD(X.getMaxD());

    Yp = D * X + One;       if (verbose > 1) cout << "\tD(C(z))\t= " << Yp << "\n";
    if (m%2 == 0) {
        Yn = Fourier(X,t);  if (verbose > 1) cout << "\tC(-z)\t= " << Yn << "\n";
        Yn = D * Yn + One;  if (verbose > 1) cout << "\tD(C(-z))\t= " << Yn << "\n";
        Yp *= Yn;           if (verbose > 1) cout << "\tD(C(z))D(C(-z))\t= " << Yp << "\n";
    }
    for (unsigned int h = 1; h<=(m-1)/2; h++) {
        Yu = Fourier(X, h, t, m);
        Yd = Fourier(X, m-h, t, m);
        if (verbose > 1) cout << "\t\tC"<< h << "\t= " << Yu << "\n\t\tC"<< m-h << "\t= " << Yd << "\n";
        Yu = Fourier(D,0,1,m) * Yu + Fourier(One,0,1,m);
        Yd = Fourier(D,0,1,m) * Yd + Fourier(One,0,1,m);
        if (verbose > 1) cout << "\t\tD(C"<< h << ")\t= " << Yu << "\n\t\tD(C"<< m-h << ")\t= " << Yd << "\n";
        Yu *= Yd;              if (verbose > 1) cout << "\tD(C"<< h << ")*D(C"<< m - h << ")\t= " << Yu << "\n";
        Yp *= FourierInv(Yu);  if (verbose > 1) cout << "\tMul\t= " << Yp << "\n";
    }
    Yp = (Yp - One)/D;
    return Yp;
}

void setMaxDegree(unsigned int s) {
    FormalPowerSeries<Rational<Long> > x;
    x.setMaxD(s);
    FormalPowerSeries<Cyclotomic<Rational<Long> > > xc;
    xc.setMaxD(s);
}

int main(){
    unsigned int k = 1, l = 2, m = 3, maxDegree = 12;
    cerr << "Maximum degree of the formal power series = "; cin >> maxDegree;
    setMaxDegree(maxDegree);
    cerr << "To compute the extended Catalan numbers C(k,l) and C(mk, ml), we set\n";
    cerr << "k = "; cin >> k;
    cerr << "l = "; cin >> l;
    cerr << "m = "; cin >> m;

    FormalPowerSeries<Rational<Long> > K_kl, C_kl,C_mkml, C_mCkl;
    clock_t start = clock();

    //for (l = 1; l <= 4; l++) for (k = 1; k <= 4; k++)
    {
    K_kl = CatalanKernel(k, l);
    C_kl = ExtendedCatalan(K_kl, 1, l);   if (verbose>0) cout << "C(" << k << "," << l << ")(z)\t= " << C_kl << "\n";
    //for (m = 2; m<=3; m++)
    {
    C_mkml = ExtendedCatalan(K_kl, 1, m * l); if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mkml << " from the Catalan Kernel K(" << k << "," << l << ")(z)\n";
    C_mCkl = ExtendedCatalan(C_kl, l, m);  if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mCkl << " from the Catalan numbers C(" << k << "," << l << ")(z)\n";
    }}
    cout << "Elaspsed time: " << (double)(clock() - start)/CLOCKS_PER_SEC;
    return 0;
}

