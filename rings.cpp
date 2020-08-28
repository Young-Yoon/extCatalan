//
// Created by euler.lee on 20. 8. 18..
//
#include "rings.h"

long long Long::BASE = 10000000, Long::Digits = 7;

void Long::unique() {  // V >= 0
    while (V.size()>1 && V.back()==0)  V.pop_back(); // 1) V.back() != 0 : remove Leading Zero
    long long carry = 0;
    for (unsigned int i=0; i<V.size(); i++) {
        carry += V[i];
        V[i] = carry % BASE;
        carry /= BASE;
    }
    while (carry > 0) {
        V.push_back(carry%BASE);
        carry /= BASE;
    }
    if (V.size()==1 && V[0] == 0) sign=0;  // when Zero
}

const Long Long::mag() const {
    Long mag = *this;
    if (mag.sign==-1) mag.sign=1;
    return mag;
}

const Long abs(const Long & x) {
    return x.mag();
}

Long Long::addV(const Long& A, const Long& B){  // return A+B for A,B>0
    Long sum = (A.V.size()>=B.V.size() ? A : B);
    const Long *small = (A.V.size()>=B.V.size() ? &B : &A);
    for (unsigned int i = 0; i < small->V.size(); i++) sum.V[i] += small->V[i];
    sum.unique();
    return sum;
}

Long Long::subV(const Long& A, const Long& B){  // return A-B for A>B>0
    Long dif;
    unsigned int md = A.V.size();  // md = argmax_d A.V[d-1]!=B.V[d-1]
    if (A.V.size()==B.V.size()) {  // md < A.V.size --> md==0 || A.V[md-1] > B.V[md-1]
        while (md>0 && A.V[md-1] == B.V[md-1]) md--;
    } // else --> md == A.V.size
    // if (md==0) return Long();
    dif.V.resize(md);
    for (unsigned int j=0; j<md; j++) dif.V[j] = A.V[j]+BASE-1;
    dif.V[0] += 1;
    dif.V[md-1] -= BASE;
    for (unsigned int j=0; j<(md==A.V.size() ? B.V.size() : md); j++)
        dif.V[j] -= B.V[j];
    dif.unique();
    return dif;
}

Long& Long::operator+= (const Long& A) {
    Long sum = *this;
    if (sign * A.sign >= 0) {  // 1:(+,+) 2:(-,-), 3:(+,0), 4:(-,0), 5:(0,0)
        sum = addV(mag(), A.mag());
        sum.sign = sign + A.sign;
        if (sum.sign>1) sum.sign = 1;
        else if(sum.sign< -1) sum.sign = -1;
    } else {  // 5+(-3), 3+(-5)
        if (mag() == A.mag()) sum = Long(0);
        else if (mag() > A.mag()) {
            sum = subV(mag(), A.mag());
            sum.sign = sign;
        } else {
            sum = subV(A.mag(), mag());
            sum.sign = A.sign;
        }
    }
    return *this = sum;
}

Long& Long::operator-= (const Long& A) {
    return *this += (-A);
}

Long& Long::operator<<= (const unsigned int& s) {
    if (*this != Long()) {
        vector<int> padZero(s, 0);
        V.insert(V.begin(), padZero.begin(), padZero.end() );
    }
    return *this;
}

Long& Long::operator>>= (const unsigned int& s) {
    if (s<V.size()) {
        V.erase(V.begin(), V.begin()+s);
    } else *this=Long();
    return *this;
}

Long& Long::operator^= (const unsigned int& s) {
    Long mul = Long(1);
    Long base = *this;
    unsigned int p = s;
    while (p > 0) {
        if (p & 1) mul *= base;
        base *= base;
        p >>= 1;
    }
    return *this = mul;
}

Long& Long::operator*= (const Long& A) {
    sign *= A.sign;
    if (sign==0) {*this=Long(); return *this;}

    Long sum = Long(0);
    for (unsigned int i = 0; i < A.V.size(); i++) {
        if (A.V[i]==0) continue;
        Long mul = *this;
        for (unsigned int j = 0; j < mul.V.size(); j++)
            mul.V[j] *= A.V[i];
        //mul.unique();
        mul <<= i;
        sum += mul;  // with unique();
    }
    sum.sign = sign;
    return *this = sum;
}

Long& Long::operator/= (const Long& C) {
    if (C==Long()) return *this;  // Error: divided by zero
    Long Q, sum = Long();
    Long R = mag(), A = C.mag();
    long long q = 1, a = A.V.back(), dig = 1;
    // b = log10 BASE, lB = log10 *this = b*nB+rB, lA = log10 A = b*nA+rA,
    // if nA >0 : la = log10 a = b, ld = log10 dig = b-rA;
    // if nA==0 : la = rA, ld = 0
    // lq = log10 q = b+rB-la+ld = b+rB-rA
    // lQ = log10 Q = lq+b*(nB-1-nA)
    // log10 A*Q = lA+lQ = b*nB+rB
    if (A.V.size()>1) {
        a *= BASE;
        a += A.V[A.V.size() - 2];
        dig = BASE;
        while (a >= BASE) {
            a /= 10;
            dig /= 10;
        }
        a += 1;
    } // BASE/10 + 1 <= a <= BASE
    while (R.V.size()>A.V.size() && q>0) {
        q = (R.V.back() * BASE + R.V[R.V.size() - 2]) / a * dig;
        Q = Long(q)<<(R.V.size()-1-A.V.size());
        sum+=Q;
        Q *= A;
        R -= Q;
    } // --> R.V.size()==A.V.size() || q==0
    if (R >= A) {
        q = R.V.back() / a * dig;
        Q = Long(q)<<(R.V.size()-A.V.size());
        sum+=Q;
        Q *= A;
        R -= Q;
    }
    while ( R >= A ) {
        sum += 1;
        R -= A;
    }
    if (sum.sign != 0) sum.sign = sign * C.sign;  // if sum is not zero
    return *this = sum;
}

Long& Long::operator%= (const Long& A) {
    if (A==Long()) return *this;  // Error: divided by zero
    if (A==Long(1)) return *this = Long(0);  // Simple & avoid infinte loops
    Long B = (*this / A) * A;
    *this -= B;
    return *this;
}

bool Long::operator==(const Long& A) const {
    if (sign!=A.sign || V.size()!=A.V.size()) return false;
    for (unsigned int i=0; i<V.size(); i++) if (V[i]!=A.V[i]) return false;
    return true;
}

bool Long::operator<(const Long& A) const {
    // 1:(1,-1), 2:(1,0), 3:(-1,0), 4:(11,1), 5:(-11,-1), 6:(2,1), 7:(-2,-1), 8:(1,1), 9:(-1,-1), 10:(0,0)
    if (*this == A) return false;  // Done 8,9,10
    if (sign != A.sign) return (sign<A.sign);  // Done 1,2,3
    if (V.size() != A.V.size()) return (V.size() < A.V.size() ? sign==1 : sign==-1);  // Done 4,5
    unsigned int i=V.size();
    while (i>0 && V[i-1]==A.V[i-1]) i--;
    return (V[i-1]<A.V[i-1] ? sign==1 : sign==-1);  // Done 6,7
}

// Print types of numbers
ostream& operator<< (ostream& os, Long p) {
    if (p.sign == 0) os<<"0";
    else {
        os << (p.sign == -1 ? "-":"") << p.V.back();
        for (unsigned int i = p.V.size()-1; i > 0; i--)
            os << setfill('0') << setw(p.Digits) << p.V[i-1];
    }
    return os;
}
