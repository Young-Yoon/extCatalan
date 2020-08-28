#include<iostream>
//#include <vector>
#include <time.h>
#include "rings.h"  // Long,
#include "combinatorics.h"   //
using namespace std;
int verbose = 2;
/*
Solve {X(z)}^{k+l} = {Z(z)}^l, Z(z)=(X(z)-1)/z  for X = c0 + c1 z + c2 z^2 + ...
 */
FormalPowerSeries<Rational<Long> > CatalanKernel(unsigned int k, unsigned int l) {
    if (verbose > 0) cout << "Solve x^" << k+l << " = ((x-1)/z)^" << l << " for x\n";
    FormalPowerSeries<Rational<Long> > X, Z;
    Rational<Long> kx, kz;
    unsigned int maxD = X.getMaxD();
    X[0] = 1, Z[0] = X[1] = 1; Z[1] = 0;
    for (unsigned int N=1; N<maxD; N++) {
        kx = polynomialPowerCoef(X,  N, k+l); // [z^N] X(z)^{k+l}
        kz = polynomialPowerCoef(Z, N, l);     // [z^N] Z(z)^l
        Z[N] = X[N+1] = (kx - kz) / l;
        if (verbose > 1) cout << "\ta" << N+1 << "= (" << kx << "-" << kz << ")/" << l << "=" << Z[N] << "\n";
        Z[N+1] = 0;
    }
    if (verbose > 0) cout << "A(" << k << "," << l << ")(z)\t= " << X << "\n";
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

FormalPowerSeries<Rational<Long> > ExtendedCatalan(FormalPowerSeries<Rational<Long> > X, unsigned int l, unsigned int m) {
    FormalPowerSeries<Rational<Long> > One, D, Yp, Yn;
    One[0] = 1;
    D[0] = -1; D[1] = 1; // D(z) = z-1

    FormalPowerSeries<Cyclotomic<Rational<Long> > > Yu(m), Yd(m), Mul(m);
    //Yu.setMaxD(X.getMaxD());

    Yp = D * X + One;       if (verbose > 1) cout << "\tD(C0)\t= " << Yp << "\n";
    if (m%2 == 0) {
        Yn = Fourier(X,l);  if (verbose > 1) cout << "\t\tC" << m*l/2 << "\t= " << Yn << "\n";
        Yn = D * Yn + One;  if (verbose > 1) cout << "\t\tD(C" << m*l/2 << ")\t= " << Yn << "\n";
        Yp *= Yn;           if (verbose > 1) cout << "\tD(C0)*D(C" << m*l/2 << ")\t= " << Yp << "\n";
    }
    Mul = Fourier(Yp,0,1,m);
    for (unsigned int h = 1; h<=(m-1)/2; h++) {
        Yu = Fourier(X, h, l, m);
        Yd = Fourier(X, m*l-h, l, m);
        if (verbose > 1) cout << "\t\tC"<< h << "\t= " << Yu << "\n\t\tC"<< m-h << "\t= " << Yd << "\n";
        Yu = Fourier(D,0,1,m) * Yu + Fourier(One,0,1,m);
        Yd = Fourier(D,0,1,m) * Yd + Fourier(One,0,1,m);
        if (verbose > 1) cout << "\t\tD(C"<< h << ")\t= " << Yu << "\n\t\tD(C"<< m-h << ")\t= " << Yd << "\n";
        Yu *= Yd;        if (verbose > 1) cout << "\tD(C"<< h << ")*D(C"<< m - h << ")\t= " << Yu << "\n";
        Mul*= Yu;        if (verbose > 1) cout << "\tMul\t= " << Mul << "\n";
    }
    Yp = FourierInv(Mul);
    Yp = (Yp - One)/D;
    return Yp;
}

FormalPowerSeries<Rational<Long> > ElementarySymmetric(FormalPowerSeries<Rational<Long> > X, unsigned int t) {

}


void setMaxDegree(unsigned int s) {
    FormalPowerSeries<Rational<Long> > x;
    x.setMaxD(s);
    FormalPowerSeries<Cyclotomic<Rational<Long> > > xc;
    xc.setMaxD(s);
}

void inputParam(unsigned int &deg, unsigned int &k, unsigned int &l, unsigned int &m, string app) {
    cerr << "Maximum degree of the formal power series = "; cin >> deg;
    setMaxDegree(deg);
    cerr << app << "\n";
    cerr << "k = "; cin >> k;
    cerr << "l = "; cin >> l;
    cerr << "m = "; cin >> m;
}

unsigned int gk = 4, gl = 4, gm = 4, gDeg = 20;

void showKernelEquivalence(int mode) {
    unsigned int k = 1, l = 2, m = 3, deg = gDeg;
    FormalPowerSeries<Rational<Long> > A_kl, A_mkml;
    clock_t start;
    string app = "To see the equivalence of the Catalan kernels A(k,l) and A(mk, ml), we set";
    if (mode == 1) {
        inputParam(deg, k, l, m, app);
        verboseCombin = 1;
        start = clock();
        A_kl = CatalanKernel(k, l);
        A_mkml = CatalanKernel(m * k, m * l);
        verboseCombin = 0;
    } else {
        setMaxDegree(deg);
        start = clock();
        for (l = 1; l <= gl; l++)
            for (k = 1; k <= gk; k++) {
                if (gcd(max(k,l), min(k,l)) > 1) continue;
                A_kl = CatalanKernel(k, l);
                for (m = 2; m<=gm; m++) A_mkml = CatalanKernel(m * k, m * l);
            }
    }
    cout << "Elaspsed time: " << (double)(clock() - start)/CLOCKS_PER_SEC << "\n";
}

void showCatalanEquality(int mode) {
    unsigned int k = 1, l = 2, m = 3, deg = gDeg;
    FormalPowerSeries<Rational<Long> > A_kl, C_kl,C_mkml, C_mCkl;
    clock_t start;
    string app = "To compute the extended Catalan numbers C(k,l) and C(mk, ml), we set\n";
    if (mode == 1) {
        inputParam(deg, k, l, m, app);
        start = clock();
        A_kl = CatalanKernel(k, l);
        C_kl = ExtendedCatalan(A_kl, 1, l);   if (verbose>0) cout << "C(" << k << "," << l << ")(z)\t= " << C_kl << "\n";
        C_mkml = ExtendedCatalan(A_kl, 1, m * l); if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mkml << " from the Catalan Kernel K(" << k << "," << l << ")(z)\n";
        C_mCkl = ExtendedCatalan(C_kl, l, m);  if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mCkl << " from the Catalan numbers C(" << k << "," << l << ")(z)\n";
    } else {
        verbose = 1;
        setMaxDegree(deg);
        start = clock();
        for (l = 1; l <= gl; l++) {
            for (k = 1; k <= gk; k++) {
                if (gcd(max(k,l), min(k,l)) > 1) continue;
                A_kl = CatalanKernel(k, l);
                C_kl = ExtendedCatalan(A_kl, 1, l);   if (verbose>0) cout << "C(" << k << "," << l << ")(z)\t= " << C_kl << "\n";
                for (m = 2; m<=gm; m++) {
                    C_mkml = ExtendedCatalan(A_kl, 1, m * l); if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mkml << " from the Catalan Kernel K(" << k << "," << l << ")(z)\n";
                    C_mCkl = ExtendedCatalan(C_kl, l, m);  if(verbose>0) cout << "C(" << m * k << "," << m * l << ")(z)\t= " << C_mCkl << " from the Catalan numbers C(" << k << "," << l << ")(z)\n";
                }
            }
        }
    }
    cout << "Elaspsed time: " << (double)(clock() - start)/CLOCKS_PER_SEC << "\n";
}
int main(){
    int mode = 1, app = 1;
    cerr << "Mode 1: individual case with details.\nMode 2: simple cases\n";
    cerr << "Select mode = "; cin >> mode;
    cerr << "App 1: Kernel Equivalence\nApp 2: Catalan Equality\nApp 3: Elementary Symmetric Polynomial\n";
    cerr << "Select app = "; cin >> app;

    if (app == 1) showKernelEquivalence(mode);
    else if (app == 2) showCatalanEquality(mode);

    return 0;
}

