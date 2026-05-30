# Q3D LOM Analysis — From Maxwell Capacitance Matrix to $E_C$ for a QPD
If one wants to compute E_C for QPD, there are two ways (in both cases, there is some ambiguity due to bare qubit VS dressed qubit):
- Use HFSS: Use the qubit-like HFSS eigenmode (not f_ND) and the plasma frequency equation and the estimated E_J by AB formula (bare)
- Use Q3D: use L_J by AB formula and then capacitance matrix. You need to have 4 nodes and the ground ring: the inductor half of the resonator (floating, Q=0), the capacitor part of the resonator (floating, Q=0), the upper half of pad and the lower half of pad. Then you will compute the Schur-reduced Maxwell matrix:

$$
\left(\begin{array}{l}{{Q_{1}}}\\ {{Q_{2}}}\\ {{0}}\end{array}\right)=\left(\begin{array}{l l l l}{{C_{11}}}&{{C_{12}}}&{{C_{13}}}&{{C_{14}}}\\ {{C_{21}}}&{{C_{22}}}&{{C_{23}}}&{{C_{24}}}\\ {{C_{31}}}&{{C_{32}}}&{{C_{33}}}&{{C_{34}}}\\ {{C_{41}}}&{{C_{42}}}&{{C_{43}}}&{{C_{44}}}\end{array}\right)\left({V_{3}}\right)
$$

- Node 1: Qubit Pad 1
- Node 2: Qubit Pad 2
- Node 3: Floating Metal A
- Node 4: Floating Metal B

Because nodes 3 and 4 are completely floating, they have no path to ground to shed or gain electrons. Therefore, their boundary conditions are $Q_3=Q_4=0$ 

The simplification will lead to 

$$
C_{\Sigma}=\frac{C_{11}^{\prime}C_{22}^{\prime}-\left(C_{12}^{\prime}\right)^{2}}{C_{11}^{\prime}+C_{22}^{\prime}+2C_{12}^{\prime}}
$$

$$
C_{11}^{\prime}=C_{11}-\frac{C_{13}^{2}C_{44}-2C_{13}C_{14}C_{34}+C_{14}^{2}C_{33}}{C_{33}C_{44}-C_{34}^{2}}
$$

$$
C_{22}^{\prime}=C_{22}-\frac{C_{23}^{2}C_{44}-2C_{23}C_{24}C_{34}+C_{24}^{2}C_{33}}{C_{33}C_{44}-C_{34}^{2}}
$$

$$
C_{12}^{\prime}=C_{12}-\frac{C_{13}C_{23}C_{44}-C_{14}C_{23}C_{34}-C_{13}C_{24}C_{34}+C_{14}C_{24}C_{33}}{C_{33}C_{44}-C_{34}}
$$
