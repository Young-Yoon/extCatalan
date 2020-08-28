//
// Created by euler.lee on 20. 8. 20..
//
#ifndef CATALAN_COMBINATORICS_H
#define CATALAN_COMBINATORICS_H

#include<iostream>
#include <vector>
#include "rings.h"  // typename T can be Long, Rational, Polynomial, Cyclotomic.

using namespace std;

int verboseCombin = 0;  // 0: Print-off   1: Print details of computation

/* In computing [z^N] {x(z)}^H,
 * this function gives "ans = C * P", where
 *   x(z) = c0 + c1 z + c2 z^2 + ...
 *   H = sum_{i=0}^{n} hi,
 *   N = sum_{i=0}^{n} i*hi,
 *   C = combin{H}{h0, h1,..., hn},
 *   P = Prod_{i=0}^{n} ci^hi.
 */
template<typename T>
T prodCoef(const vector<T>& c, const vector<unsigned int>& h) {
    T C = T(1), P = T(1);
    unsigned int H = 0;
    for (unsigned int i = 0; i < h.size(); i++) H += h[i];
    for (unsigned int i = 0; i < h.size(); i++) {
        for (unsigned int j = 0; j < h[i]; j++) {
            C *= H;
            C /= (j + 1);
            H--;
        }
        P *= c[i] ^ h[i];
    }
    if (verboseCombin > 0) cout << C << "*\t";
    return C * P;
}

/* Integer Weak Partition of N with H elements
 * N = sum_i i*hist[i-1]
 * For example, N=5, H=4
 *        hist  [0 1 2 3 4 5]
 * 5 = 5+0+0+0  [3 0 0 0 0 1]
 *   = 4+1+0+0  [2 1 0 0 1 0]
 *   = 3+2+0+0  [2 0 1 1 0 0]
 *   = 3+1+1+0  [1 2 0 1 0 0]
 *   = 2+2+1+0  [1 1 2 0 0 0]
 *   = 2+1+1+1  [0 3 1 0 0 0]
 */
template<typename T>
T integerWeakPartition(unsigned int N, unsigned int C, unsigned int H, vector<unsigned int> * const hist,
                    T (*valFcn)(const vector<T>&, const vector<unsigned int>&), const vector<T>& coef ) {
    T value = T(0);
    if (N == 0 || C == 1) {
        hist->at(1) = N;
        hist->front() -= N;
        value += valFcn(coef, *hist);
        if (verboseCombin > 0) {
            for (unsigned int i = 0; i < hist->size(); i++) { cout << "(" << coef[i] << "z^" << i << ")^" << hist->at(i) << "\t"; }
            cout << "= " << value << "\n";
        }
        hist->at(1) = 0;
        hist->front() += N;
        return value;
    }
    for (int h = N/C; h >= 0; --h) {
        hist->at(C) = (unsigned int)h;
        hist->front() -= h;
        int nextN = N-C*h;
        int nextC = nextN < C ? nextN : C-1;
        int nextH = H-h;
        if (nextN > nextC * nextH) {
            hist->front() += h;
            hist->at(C) = 0;
            break;
        }
        for (unsigned int i=nextC+1; i<C; i++) hist->at(i) = 0;
        value += integerWeakPartition(nextN, nextC, nextH, hist, valFcn, coef);
        hist->front() += h;
        hist->at(C) = 0;
    }
    return value;
}

/* For  x(z) = c0 + c1 z + c2 z^2 + ...
 * return [z^N] {x(z)}^H
 */
template<typename T>
T polynomialPowerCoef(const FormalPowerSeries<T> &x, unsigned int N, unsigned int H) {
    vector<unsigned int> *hist = new vector<unsigned int>;
    hist->assign(N+1,0);
    hist->front() = H;
    const vector<T> Coef = x.readCoef();    // C2662

    if (verboseCombin > 0) {
        cout << "[z^" << N << "]\t(" << Coef[0] << " +\t";
        for (unsigned int i = 1; i <= N; i++) cout << Coef[i] << "z^" << i << "+\t";
        cout << "..)^" << H << "\n";
    }
    return integerWeakPartition(N, N, H, hist, prodCoef, Coef);
}
/* Example : N=5, H=4, x(z) = 1 + z + 2z^2 + 5z^3 + 14z^4 + 42z^5 + ..
[z^5]	(1 +	    1z^1+	    2z^2+	    5z^3+	    14z^4+	    42z^5+	..)^4
4*	    (1z^0)^3	(1z^1)^0	(2z^2)^0	(5z^3)^0	(14z^4)^0	(42z^5)^1	= 168
12*	    (1z^0)^2	(1z^1)^1	(2z^2)^0	(5z^3)^0	(14z^4)^1	(42z^5)^0	= 168
12*	    (1z^0)^2	(1z^1)^0	(2z^2)^1	(5z^3)^1	(14z^4)^0	(42z^5)^0	= 120
12*	    (1z^0)^1	(1z^1)^2	(2z^2)^0	(5z^3)^1	(14z^4)^0	(42z^5)^0	= 60
12*	    (1z^0)^1	(1z^1)^1	(2z^2)^2	(5z^3)^0	(14z^4)^0	(42z^5)^0	= 48
4*	    (1z^0)^0	(1z^1)^3	(2z^2)^1	(5z^3)^0	(14z^4)^0	(42z^5)^0	= 8
 */

#endif //CATALAN_COMBINATORICS_H


/* Select H elements from N with lower limit M
 * For example, S_k,l = A1, A2, A3, ... , Al
 *
 * 1 2 3  [1 1 1 0 0]
 * 1 2 4  [1 1 0 1 0]
 * 1 2 5  [
 * 1 3 4
 * 1 3 5
 * 1 4 5
 * 2 3 4
 * 2 3 5
 * 3 4 5
 *
template<typename T>
T Combin(unsigned int N, unsigned int H, unsigned int C, vector<unsigned int> * const hist, vector<vector<int> > *visit, vector<vector<T> > *val) {
    if (H == 0 || M == N) {}


    for (int h = 1; h >= 0; --h) {
        hist->at(C) = (unsigned int)h;
        hist->front() -= h;
// for (unsigned int i=0; i < hist->size(); i++) cout << hist->at(i) << "\t"; cout << "N " << N << " C " << C << " h " << h << "\n";
        int nextN = N-C*h;
        int nextC = nextN < C ? nextN : C-1;
        int nextH = H-h;
        if (nextN > nextC * nextH) {
            hist->front() += h;
            hist->at(C) = 0;
            break;
        }
        for (unsigned int i=nextC+1; i<C; i++) hist->at(i) = 0;
        printIntWeakPartitionH(nextN, nextC, nextH, hist, visit, val, valFcn, pcCoef);
        hist->front() += h;
        hist->at(C) = 0;
    }
    visit->at(N-1)[C-1] ++;

    return value;
}

if (N == 0 || C == 1) {
hist->at(1) = N;
hist->front() -= N;
T value = valFcn(pcCoef, *hist);
unsigned int sum = 0, h = 0;
for (unsigned int i=0; i < hist->size(); i++) {
cout << hist->at(i) << "\t";
h += hist->at(i);
sum += i*hist->at(i);
if (sum>0 && i<sum && visit->at(sum-1)[i]==0)
val->at(sum-1)[i] += value;
//hist->at(i) = 0;
}
cout << value << "\n";
if (N>0) visit->at(N-1)[C-1] ++;
hist->front() += N;
hist->at(1) = 0;
return;
}


}
*/