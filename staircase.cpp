#include<iostream>
//#include <vector>
#include <time.h>
#include "rings.h"  // Long,
#include "combinatorics.h"   //
using namespace std;
const int verbose = 1;
/*
Solve {X(z)}^{k+l} = {Z(z)}^l, Z(z)=(X(z)-1)/z  for X = c0 + c1 z + c2 z^2 + ...
 */
FormalPowerSeries<Rational<Long> > CatalanKernel(unsigned int k, unsigned int l) {
    if (verbose > 0) cout << "Solve x^" << k+l << " = ((x-1)/z)^" << l << " for x\n";
    FormalPowerSeries<Rational<Long> > X, Z;
    unsigned int maxD = X.getMaxD();
    X[0] = 1, Z[0] = X[1] = 1; Z[1] = 0;
    for (unsigned int N=1; N<maxD; N++) {
        Z[N] = X[N+1] = (polynomialPowerCoef(X,  N, k+l) - // [z^N] X(z)^{k+l}
                         polynomialPowerCoef(Z, N, l))/l;     // [z^N] Z(z)^l
        Z[N+1] = 0;
    }
    if (verbose > 0) cout << X << "\n";
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
    FormalPowerSeries<Rational<Long> > Yp, Yn, D, One;
    FormalPowerSeries<Cyclotomic<Rational<Long> > > Yu(m), Yd(m); //, Mul(m);
    D[0] = -1; D[1] = 1; // D(z) = z-1
    One[0] = 1;
    Yp = D * X + One;       if (verbose > 1) cout << "\tD(C(z))\t= " << Yp << "\n";
    if (m%2 == 0) {
        Yn = Fourier(X,t);  if (verbose > 1) cout << "\tC(-z)\t= " << Yn << "\n";
        Yn = D * Yn + One;  if (verbose > 1) cout << "\tD(C(-z))\t= " << Yn << "\n";
        Yp *= Yn;           if (verbose > 1) cout << "\tD(C(z))D(C(-z))\t= " << Yp << "\n";
    }
    // Mul = Fourier(Yp, 0, 1, m);
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

int main(){
    clock_t start = clock();
    int s = 1, t = 2, m = 3;
    for (t = 1; t <= 4; t++) for (s = 1; s <= 4; s++)
    {
    FormalPowerSeries<Rational<Long> > K = CatalanKernel(s, t);          if (verbose>0) cout << "K(" << s << "," << t << ")(z)\t= " << K << "\n";
    FormalPowerSeries<Rational<Long> > C = ExtendedCatalan(K, 1, t);  if (verbose>0) cout << "C(" << s << "," << t << ")(z)\t= " << C << "\n";
    for (m = 2; m<=3; m++)
    {
    FormalPowerSeries<Rational<Long> > D = ExtendedCatalan(K, 1, m * t);  if(verbose>0) cout << "C(" << m * s << "," << m * t << ")(z)\t= " << D << " derived from the Catalan Kernel K(" << s << "," << t << ")(z)\n";
    FormalPowerSeries<Rational<Long> > E = ExtendedCatalan(C, t, m);             if(verbose>0) cout << "C(" << m * s << "," << m * t << ")(z)\t= " << E << " derived from the Catalan numbers C(" << s << "," << t << ")(z)\n";
    }}
    cout << "Elaspsed time: " << (double)(clock() - start)/CLOCKS_PER_SEC;
    return 0;
}

