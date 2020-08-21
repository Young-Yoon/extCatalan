The Extended Catalan Numbers
---
### Definition 
Consider the enumeration of lattice paths in 2D plane with steps *E*=(1,0) and *N*=(0,1).</br>
Let <img align="center" height="16" src="https://latex.codecogs.com/svg.latex?\small&space;c_n^{k,l}" title="catalan" /> 
denote **the extended Catalan numbers** which equals the number of lattice paths staying weakly below 
<img align="center" height="16" src="https://latex.codecogs.com/svg.latex?\small&space;(E^lN^k)^{\lfloor\frac{n}{l}\rfloor}E^{n-l\lfloor\frac{n}{l}\rfloor}" title="path"/>.

For example, <img align="center" height="16" src="https://latex.codecogs.com/svg.latex?\small&space;c_9^{2,3}" title="catalan" /> equals 
the number of lattice paths from (0,0) to (9,6) with the possible paths as shown in bold lines of the following figure:

![Lattice](./c9_2_3.png)

Let <img align="center" height="16" src="https://latex.codecogs.com/svg.latex?\small&space;C^{k,l}(z)" title="catalanGF" /> 
denote the generating function for the extended Catalan numbers.

We present three equations.

### The Catalan Kernel  
Let <img align="center" height="16" src="https://latex.codecogs.com/svg.latex?\small&space;K^{k,l}(z)\in\mathbb{Q}(\zeta_l)[[z]]"/>
denote **the Catalan kernel** which is the primitve root of the equation 
<img align="center" height="40" src="https://latex.codecogs.com/svg.latex?\small&space;x^{k+l}=\left(\frac{x-1}{z}\right)^{l}"/>.

Using the function `CatalanKernel` in `staircase.cpp`, one can compute the Catalan kernels by
```
unsigned int k, l;
FormalPowerSeries<Rational<Long> > K_kl = CatalanKernel(k, l); 
std::cout << K << "\n";
```

### Relation between the extended Catalan numbers and the Catalan kernel  

<img align="center" height="50" src="https://latex.codecogs.com/svg.latex?\large&space;C^{k,l}(z)=\frac{1}{1-z}\left(1-\prod_{j=0}^{l-1}\left(1-(1-z)K^{k,l}(e^{i\frac{2j\pi}{l}}z)\right)\right)"/>.

Using the function `ExtendedCatalan` in `staircase.cpp`, one can compute the generating function for the extended Catalan numbers derived from the Catalan kernel as:
```
unsigned int k, l, m;
FormalPowerSeries<Rational<Long> > C_kl = ExtendedCatalan(K_kl, 1, l); 
FormalPowerSeries<Rational<Long> > C_mkml = ExtendedCatalan(K_kl, 1, m*l); 
std::cout << C_kl << "\n";
```

### Relation between the extended Catalan numbers and the simple Catalan numbers  

<img align="center" height="90" src="https://latex.codecogs.com/svg.latex?\large&space;C^{mk,ml}(z)=\frac{1}{1-z}\left(1-\prod_{h=0}^{m-1}\left(1-(1-z)\sum_{i=0}^{l-1}e^{-i\frac{2hi\pi}{ml}}\widetilde{C}_i^{k,l}(e^{i\frac{2h\pi}{ml}}z)\right)\right)=\frac{1}{1-z}\left(1-\prod_{h=0}^{m-1}\sum_{d\geq0}e^{i\frac{2hd\pi}{m}}\sum_{i=0}^{l-1}c^{k,l}_{ld+i}z^{ld+i}\right)"/>.

Using the function `ExtendedCatalan` in `staircase.cpp`, one can compute the generating function for the extended Catalan numbers derived from from the simple Catalan numbers as: 
```
unsigned int k, l, m;
FormalPowerSeries<Rational<Long> > C_mCkl = ExtendedCatalan(C_kl, l, m);
std::cout << C_mCkl << "\n";
```
