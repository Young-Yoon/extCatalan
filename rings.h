//
// Created by euler.lee on 20. 8. 18..
//

#ifndef CPP_RINGS_H
#define CPP_RINGS_H

#include <iostream>
#include <string>
#include <iomanip>  //setw, setfill
#include <vector>

using namespace std;

const int verboseRing = 1;

class Long {
private:
    vector<long long> V;  // V[0]: Least S. digit --> V.back(): Most S. digit
    int sign;        // {-1, 0, 1}  sign=0 ==> V becomes empty
    static long long BASE, Digits;

    void unique();  // Unique representation rule: 1) V.back() != 0: 2) 0<= V[i] < BASE

    static Long addV(const Long& A, const Long& B);
    static Long subV(const Long& A, const Long& B);
public:
    Long() { sign=0, V.assign(1,0); }  // Default value: 0
    template<typename T> Long(T x) {
        sign = (x == T(0) ? 0 : (x>T(0) ? 1 : -1) );
        V.assign(1, sign>=0 ? x : -x);
        unique();
    }  // explicit: restrict to "long long"

    ~Long() { V.clear(); }

    const Long mag() const;

    Long& operator+= (const Long& A);
    Long& operator-= (const Long& A);
    Long& operator*= (const Long& A);
    Long& operator/= (const Long& A);
    Long& operator%= (const Long& A);
    Long& operator<<= (const unsigned int& s);
    Long& operator>>= (const unsigned int& s);
    Long& operator^= (const unsigned int& s);
    const Long operator- () const { Long A = *this; A.sign=-sign; return A; } // minus "copy"
    const Long operator+ (const Long& A) const { return Long(*this) += A; }
    const Long operator- (const Long& A) const { return Long(*this) -= A; }
    const Long operator* (const Long& A) const { return Long(*this) *= A; }
    const Long operator/ (const Long& A) const { return Long(*this) /= A; }
    const Long operator% (const Long& A) const { return Long(*this) %= A; }
    const Long operator<< (const unsigned int& s) const { return Long(*this) <<= s; }
    const Long operator>> (const unsigned int& s) const { return Long(*this) >>= s; }
    const Long operator^ (const unsigned int& s) const { return Long(*this) ^= s; }

    bool operator==(const Long& A) const {
        if (sign!=A.sign || V.size()!=A.V.size()) return false;
        for (unsigned int i=0; i<V.size(); i++) if (V[i]!=A.V[i]) return false;
        return true;
    }
    bool operator!=(const Long& A) const { return !(*this==A); }
    bool operator<(const Long& A) const {
        // 1:(1,-1), 2:(1,0), 3:(-1,0), 4:(11,1), 5:(-11,-1), 6:(2,1), 7:(-2,-1), 8:(1,1), 9:(-1,-1), 10:(0,0)
        if (*this == A) return false;  // Done 8,9,10
        if (sign != A.sign) return (sign<A.sign);  // Done 1,2,3
        if (V.size() != A.V.size()) return (V.size() < A.V.size() ? sign==1 : sign==-1);  // Done 4,5
        unsigned int i=V.size();
        while (i>0 && V[i-1]==A.V[i-1]) i--;
        return (V[i-1]<A.V[i-1] ? sign==1 : sign==-1);  // Done 6,7
    }
    bool operator<=(const Long& A) const { return (*this<A) || (*this==A); }
    bool operator>(const Long& A) const { return !(*this<=A); }
    bool operator>=(const Long& A) const { return !(*this<A); }

    friend ostream& operator<< (ostream& os, Long p);
};

template<typename T>
T gcd(T a, T b){  // a > b >= 0
    if (b==T(0)) return a;
    else return gcd(b, a%b);
}

const Long abs(const Long & x);

template<typename T=Long>
class Rational {  // p/q
private:
    T p, q; // q>0
    void simple() {
        // if (q==0) return;  // Error: divided by zero
        if (p== T(0)) q = 1;
        else {
            T c = (abs(p)>=abs(q) ? gcd(abs(p),abs(q)) : gcd(abs(q),abs(p)));
            p /= c, q /= c;
        }
        if (q<T(0)) p=-p, q=-q;  // Should be q>0
    }
public:
    Rational() { Rational<T>(0, 1); }
    template<typename U> Rational(U _p=0, U _q=1):p(_p), q(_q) { simple(); }

    Rational<T>& operator+= (const Rational<T> & A) {
        Rational<T> sum(p*A.q + q*A.p, q*A.q);
        sum.simple();
        return *this = sum;
    }
    Rational<T>& operator-= (const Rational<T> & A) {
        //Rational<T> B = A;  // C2678: -A
        //*this += (-B);
        return *this += (-A);
    }
    Rational<T>& operator*= (const Rational<T> & A) {
        Rational<T> mul(p*A.p, q*A.q);
        mul.simple();
        return *this = mul;
    }
    Rational<T>& operator/= (const Rational<T> & A) {
        Rational<T> div(p*A.q, q*A.p);
        div.simple();
        return *this = div;
    }
    Rational<T>& operator^= (const unsigned int& s) {
        Rational<T> mul = Rational<T>(1,1);
        Rational<T> base = *this;
        unsigned int p = s;
        while (p > 0) {
            if (p & 1) mul *= base;
            base *= base;
            p >>= 1;
        }
        return *this = mul;
    }

    const Rational<T> operator- () const { Rational<T> A = *this; A.p=-p; return A; }  // minus
    const Rational<T> operator+ (const Rational<T>& A) const { return Rational<T> (*this) += A; }
    const Rational<T> operator- (const Rational<T>& A) const { return Rational<T> (*this) -= A; }
    const Rational<T> operator* (const Rational<T>& A) const { return Rational<T> (*this) *= A; }
    const Rational<T> operator/ (const Rational<T>& A) const { return Rational<T> (*this) /= A; }
    const Rational<T> operator^ (const unsigned int& s) const { return Rational<T> (*this) ^= s; }

    bool operator==(const Rational<T>& A) const { return (p==A.p) && (q==A.q); }
    bool operator!=(const Rational<T>& A) const { return !(*this==A); }
    bool operator<(const Rational<T>& A) const { return (p*A.q < A.p*q); }
    bool operator<=(const Rational<T>& A) const { return (*this<A) || (*this==A); }
    bool operator>(const Rational<T>& A) const { return !(*this<=A); }
    bool operator>=(const Rational<T>& A) const { return !(*this<A); }

    template <typename U> friend ostream & operator <<(ostream&, const Rational<U>&);
};

template<typename T=Rational<Long>, char V='z'>
class Polynomial {
protected:
    vector<T> coef;
    unsigned int dimN;
    string var;
    static Polynomial<T,V> Zero, One;

    virtual void unique();
    void resetCoef() {
        coef.erase(coef.begin()+1, coef.end());
        coef[0] *= T(0);
    }
public:
    Polynomial<T, V>(unsigned int N=0):dimN(N) { coef.assign(1, T(0)); var = V; }
    Polynomial<T, V>(const T& x, unsigned int N=0):dimN(N) {coef.assign(1, x); var = V; }
    Polynomial<T, V>(const vector<T>& x, unsigned int N=0):dimN(N) { coef.assign(x.begin(), x.end()); var = V;}
    ~Polynomial<T, V>() { coef.clear(); }

    const vector<T>& readCoef() const { return coef; }   // C2662: const Body & const Return needed
    T& operator[] (unsigned int t) {
        while (coef.size() <= t) {
            T x = coef[0];
            x *= T(0);
            coef.push_back(x);
        }
        return coef.at(t);
    }

    void setAsZero() { resetCoef(); }

    Polynomial<T,V>& operator+= (const Polynomial<T,V> & x);
    Polynomial<T,V>& operator-= (const Polynomial<T,V> & x) { return *this += (-x); }
    Polynomial<T,V>& operator*= (const T & x);
    Polynomial<T,V>& operator*= (const Polynomial<T,V> & x);
    Polynomial<T,V>& operator/= (const Polynomial<T,V> & x);
    Polynomial<T,V>& operator%= (const Polynomial<T,V> & x);
    Polynomial<T,V>& operator<<= (const unsigned int t);
    Polynomial<T,V>& operator>>= (const unsigned int t);

    const Polynomial<T,V> operator - () const;
    const Polynomial<T,V> operator + (const Polynomial<T,V> &x) const { return Polynomial<T,V> (*this) += x;}
    const Polynomial<T,V> operator - (const Polynomial<T,V> &x) const { return Polynomial<T,V> (*this) -= x;}
    const Polynomial<T,V> operator * (const T &x) const { return Polynomial<T,V> (*this) *= x;}
    const Polynomial<T,V> operator * (const Polynomial<T,V> &x) const { return Polynomial<T,V> (*this) *= x;}
    const Polynomial<T,V> operator / (const Polynomial<T,V> &x) const { return Polynomial<T,V> (*this) /= x;}
    const Polynomial<T,V> operator % (const Polynomial<T,V> &x) const { return Polynomial<T,V> (*this) %= x;}
    const Polynomial<T,V> operator << (const unsigned int t) const { return Polynomial<T,V> (*this) <<= t;}
    const Polynomial<T,V> operator >> (const unsigned int t) const { return Polynomial<T,V> (*this) >>= t;}

    bool operator==(const Polynomial<T,V> &x) const {
        if (coef.size()!=x.coef.size()) return false;
        for (unsigned int i=0; i<coef.size(); i++) if (coef[i]!=x.coef[i]) return false;
        return true;
    }
    bool operator!=(const Polynomial<T,V> &x) const { return !(*this==x); }

    template<typename U, char W> friend ostream & operator <<(ostream&os, const Polynomial<U,W>& p);
};
template<typename T, char V> Polynomial<T,V> Polynomial<T,V>::Zero;
template<typename T, char V> Polynomial<T,V> Polynomial<T,V>::One(T(1));


template<typename T, char V>
void Polynomial<T,V>::unique() {
    while (coef.size()>1 && coef.back()==T(0))  coef.pop_back(); // remove Zero highest coefficients
}

template<typename T, char V>
const Polynomial<T,V> Polynomial<T,V>::operator- () const {
    Polynomial<T,V> A(*this);
    for (unsigned int i=0; i<A.coef.size(); i++)
        A.coef[i] = -A.coef[i];
    return A;
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator+= (const Polynomial<T,V> &x) {
    //Polynomial<T> sum = *this;
    for (unsigned int i=0; i<x.coef.size(); i++)
        (*this)[i] += x.coef[i];
    unique();
    return *this;  // C4172: return sum: return value disappears
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator*= (const T &x) {
    if (x==T(0)) resetCoef(); //*this = Polynomial<T,V>();
    else
        for (unsigned int i=0; i<coef.size(); i++) coef[i] *= x;
    return *this;
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator<<= (const unsigned int t) {
    if (*this != Polynomial<T,V>()) {
        T zero(coef[0]); zero*=T(0);
        vector<T> zeroCoefs(t, zero);
        coef.insert(coef.begin(), zeroCoefs.begin(), zeroCoefs.end());
    }
    unique();
    return *this;
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator>>= (const unsigned int t) {
    if (t < coef.size()) coef.erase(coef.begin(), coef.begin()+t);
    else *this = Polynomial<T,V>();
    return *this;
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator*= (const Polynomial<T,V> &x) {
    Polynomial<T,V> Zero = *this;  Zero.resetCoef();
    if (*this == Zero || x == Zero) return *this = Zero; // if x==0 skip unique()
    Polynomial<T,V> A = Zero;
    for (unsigned int i = 0; i < coef.size(); i++)
        A += ((x * coef[i]) << i);
    // A.unique(); // Only Cyclic.unique() needed
    return *this = A;
}


template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator/= (const Polynomial<T,V> &x) {
    Polynomial<T,V> Zero(*this);   Zero.resetCoef();
    if (x.coef.back() == T(0) || (coef.size() > 1 && coef.back()==T(0))) {cerr<<"Zero Highest Coeff"; return *this = Zero;}
    Polynomial<T,V> Q(Zero), R(*this);
    while (R != Zero && R.coef.size() >= x.coef.size()) {
        unsigned int s = R.coef.size() - x.coef.size();
        Q[s] = R.coef.back() / x.coef.back();
        R -= (x * Q[s])<<s;
    }
    return *this = Q;
}

template<typename T, char V>
Polynomial<T,V>& Polynomial<T,V>::operator%= (const Polynomial<T,V> &x) {
    Polynomial<T,V> Zero(*this);   Zero.resetCoef();
    if (x == Zero) { cerr<<"Divided by zero"; return *this = Zero; }
    Polynomial<T,V> y = *this / x;
    y *= x;
    if (y!=Zero) *this -= y;  // Avoid Loops
    return *this;
}

template<typename T=Rational<Long> >
class FormalPowerSeries: public Polynomial<T> {
protected:
    static unsigned int maxD;  // c_0 + c_1 z + ... + c_maxD z^maxD
    virtual void unique();
    using Polynomial<T>::coef;   // Use Parent Class Members
    using Polynomial<T>::var;

public:
    FormalPowerSeries<T>() : Polynomial<T>() { }
    FormalPowerSeries<T>(const T& x) : Polynomial<T>(x) { }
    FormalPowerSeries<T>(vector<T> x) : Polynomial<T>(x) { }
    ~FormalPowerSeries<T>() { };

    const unsigned int getMaxD() const { return maxD; }
    using Polynomial<T>::readCoef;

    FormalPowerSeries<T>& operator+= (const FormalPowerSeries<T> & x) { Polynomial<T>::operator+=(x); return *this; }
    FormalPowerSeries<T>& operator-= (const FormalPowerSeries<T> & x) { Polynomial<T>::operator+=(-x); return *this; } //FormalPowerSeries<T> y = x; *this += (-y); return *this; }
    FormalPowerSeries<T>& operator*= (const T & x) { Polynomial<T>::operator*=(x); return *this; }
    FormalPowerSeries<T>& operator*= (const FormalPowerSeries<T> & x) { Polynomial<T>::operator*=(x); return *this; }
    FormalPowerSeries<T>& operator/= (const FormalPowerSeries<T> & x);
    FormalPowerSeries<T>& operator<<= (const unsigned int t) { Polynomial<T>::operator<<=(t); return *this; }
    FormalPowerSeries<T>& operator>>= (const unsigned int t) { Polynomial<T>::operator>>=(t); return *this; }

    const FormalPowerSeries<T> operator - () const { //FormalPowerSeries<T> A(*this); A.Polynomial<T>::operator-(); return A; }
        FormalPowerSeries<T> A = *this;
        for (unsigned int i = 0; i < coef.size(); i++)
            A.coef[i] = -coef[i];
        return A;
    }
    const FormalPowerSeries<T> operator + (const FormalPowerSeries<T> &x) const { return FormalPowerSeries<T> (*this) += x;}
    const FormalPowerSeries<T> operator - (const FormalPowerSeries<T> &x) const { return FormalPowerSeries<T> (*this) -= x;}
    const FormalPowerSeries<T> operator * (const T &x) const { return FormalPowerSeries<T> (*this) *= x;}
    const FormalPowerSeries<T> operator * (const FormalPowerSeries<T> &x) const { //return FormalPowerSeries<T> (*this) *= x;}
        FormalPowerSeries<T> A = *this;
        A *= x;
        return A;
    }
    const FormalPowerSeries<T> operator / (const FormalPowerSeries<T> &x) const { return FormalPowerSeries<T> (*this) /= x;}
    const FormalPowerSeries<T> operator << (const unsigned int t) const { return FormalPowerSeries<T> (*this) <<= t;}
    const FormalPowerSeries<T> operator >> (const unsigned int t) const { return FormalPowerSeries<T> (*this) >>= t;}

    bool operator==(const FormalPowerSeries<T> &x) const {
        if (coef.size()!=x.coef.size()) return false;
        for (unsigned int i=0; i<coef.size(); i++) if (coef[i]!=x.coef[i]) return false;
        return true;
    }
    bool operator!=(const FormalPowerSeries<T> &x) const { return !(*this==x); }
};
template<typename T> unsigned int FormalPowerSeries<T>::maxD = 7; //11

template<typename T>
void FormalPowerSeries<T>::unique() {
    if (coef.size()>maxD+1) coef.erase(coef.begin()+maxD+1, coef.end());  // Ignore highest coefficients
    Polynomial<T>::unique();
}

template<typename T>
FormalPowerSeries<T>& FormalPowerSeries<T>::operator/= (const FormalPowerSeries<T> &x) {
    FormalPowerSeries<T> Zero(*this);   Zero.resetCoef();
    FormalPowerSeries<T> Q(Zero), A(Zero), R = *this;
    unsigned int s = 0;
    while (s<x.coef.size() && x.coef[s] == Zero.coef[0]) s++;
    if (s==x.coef.size()) { cerr<<"div error"; return *this = Zero; }
    A = x >> s;
    while (R != Zero && s <= maxD) {
        Q[s] = R[0] / A.coef[0];
        R -= A * Q[s];
        R >>= 1;
        s++;
    }
    return *this = Q;
}

template<typename T=Rational<Long>, char V='w'>
class Cyclotomic:public Polynomial<T,V> {
private:
    using Polynomial<T,V>::dimN;
    virtual void unique();
    using Polynomial<T,V>::coef;   // Use Parent Class Members
    using Polynomial<T,V>::var;
    static vector<Cyclotomic<T,V> > Phi;   // Lists of the cyclotomic polynomials Phi(n)
    static vector<unsigned int> phiIdx;
    static Cyclotomic<T,V> Zero, One;
    unsigned int get_phiIdx(unsigned int n);
public:
    Cyclotomic<T,V> (unsigned int N=1) : Polynomial<T,V>(N) { if (N>1) var+= to_string(N); }
    Cyclotomic<T,V> (T x, unsigned int N=1) : Polynomial<T,V>(x,N) { if (N>1) var+= to_string(N); }
    Cyclotomic<T,V> (const Cyclotomic<T,V>& x) : Polynomial<T,V>(x.coef, x.dimN) { if (x.dimN>1) var+= to_string(x.dimN); }

    using Polynomial<T,V>::readCoef;

    Cyclotomic<T,V>& operator+= (const Cyclotomic<T,V> & x) { Polynomial<T, V>::operator+=(x);  return *this; }
    Cyclotomic<T,V>& operator-= (const Cyclotomic<T,V> & x) { Polynomial<T, V>::operator+=(-x); return *this; }

    Cyclotomic<T,V>& operator*= (const T & x) { Polynomial<T,V>::operator*=(x); unique(); return *this; }
    Cyclotomic<T,V>& operator<<= (const unsigned int t) { Polynomial<T,V>::operator<<=(t); return *this; }
    Cyclotomic<T,V>& operator>>= (const unsigned int t) { Polynomial<T,V>::operator>>=(t); return *this; }
    Cyclotomic<T,V>& operator*= (const Cyclotomic<T,V> &x) { Polynomial<T,V>::operator*=(x); return *this; }
    Cyclotomic<T,V>& operator/= (const Cyclotomic<T,V> & x);

    const Cyclotomic<T,V> operator- () const {
        Cyclotomic<T,V> A(*this);
        for (unsigned int i=0; i<A.coef.size(); i++) A.coef[i] = -coef[i];
        return A;
    }
    
    const Cyclotomic<T,V> operator+ (const Cyclotomic<T,V> &x) const { return Cyclotomic<T,V>(*this) += x; }
    const Cyclotomic<T,V> operator- (const Cyclotomic<T,V> &x) const { return Cyclotomic<T,V>(*this) -= x; }

    const Cyclotomic<T,V> operator* (const T &x) const { return Cyclotomic<T,V> (*this) *= x;}
    const Cyclotomic<T,V> operator* (const Cyclotomic<T,V> &x) const { return Cyclotomic<T,V>(*this) *= x; }
    const Cyclotomic<T,V> operator/ (const Cyclotomic<T,V> &x) const { return Cyclotomic<T,V>(*this) /= x; }

    const Cyclotomic<T,V> operator << (const unsigned int t) const { return Cyclotomic<T,V> (*this) <<= t;}
    const Cyclotomic<T,V> operator >> (const unsigned int t) const { return Cyclotomic<T,V> (*this) >>= t;}

    bool operator==(const Cyclotomic<T,V> &x) const { return Polynomial<T,V>::operator==(x); }
    bool operator!=(const Cyclotomic<T,V> &x) const { return !(*this==x); }
};

template<typename T, char V> vector<Cyclotomic<T,V> > Cyclotomic<T,V>::Phi; // Phi(n) for 1<n<dimN, n|dimN
template<typename T, char V> vector<unsigned int> Cyclotomic<T,V>::phiIdx;
template<typename T, char V> Cyclotomic<T,V> Cyclotomic<T,V>::Zero;

// Phi={x2,x4,x3} => phiIdx={0,1,3,2,0} => xn=Phi[phiIdx[n-1]-1] : get_phiIdx(n) returns phiIdx[n-1]
template<typename T, char V> unsigned int Cyclotomic<T,V>::get_phiIdx(unsigned int n) {
    while (n > phiIdx.size()) phiIdx.push_back(0);     // phiIdx[n] = 0 if Phi[n] not ready.
    if (phiIdx[n - 1] > 0) return phiIdx[n - 1];       // Phi(n) = Phi[phiIdx[n-1]-1]
    Cyclotomic<T, V> x(n);
    for (unsigned int i = 0; i < n; i++) x[i] = T(1);  // (x^n)/(x-1)=1+x+x^2+x^{n-1}
    unsigned int p = 2;
    while (p * p <= n && n % p != 0) p++;              // check if n is prime
    if (p * p <= n)                                    // Phi(n) = (x^n - 1)/Prod_{d|n, d<n} Phi(d)
        for (unsigned int d = 2; d <= n / 2; d++)
            if (n % d == 0)
                x.Polynomial<T, V>::operator/=(Phi[get_phiIdx(d) - 1]);
    Phi.push_back(x);
    phiIdx[n - 1] = Phi.size();
    if (verboseRing > 0) {
        cout << "For Phi(" << n << ") phiIdx.size=" << phiIdx.size() << " Phi.size=" << Phi.size() << "\n";
        for (unsigned int i = 0; i < phiIdx.size(); i++) {
            cout << "Phi(" << i+1 << ") in list(" << (int)phiIdx[i]-1 << ")";
            if (phiIdx[i]>0) cout << "=" << Phi[phiIdx[i]-1];
            cout << "\n";
        }
    }
    return phiIdx[n-1];
}

template<typename T, char V>
void Cyclotomic<T, V>::unique() {
    Polynomial<T,V>::unique();  // Remove the highest order zeroes first.
    Polynomial<T,V>::operator%=(Phi[get_phiIdx(dimN)-1]);
}

template<typename T>
ostream& operator << (ostream& os, const Rational<T>& A) {
    os<<A.p;
    if (A.q!=0 && A.q!=1) std::cout<<"/"<<A.q;
    return os;
}

template<typename T, char V>
ostream& operator << (ostream& os, const Polynomial<T, V>& p) {
    os << p.coef[0];
    for (unsigned int j = 1; j < p.coef.size(); j++) {
        //if (p.coef[j] > 0) std::cout << " +";
        if (p.coef[j] != T(0)) {
            std::cout << "+(" << p.coef[j] << ")" << p.var;
            if (j > 1) std::cout << "^" << j;
        }
    }
    return os;
}

#endif //CPP_RINGS_H
