# Course Project

[TOC]

# Q.1:

## Solution:

### (a) Use Composite Trapezoidal Rule:

$$
\begin{eqnarray}
h=\frac{b-a}{n}&,&r_j=a+jh\\
f(r)&=&\rho\upsilon2\pi r\\
\int_{0}^{R}f(r)\mathrm{d}r&=&\frac{h}{2}\left [f(a)+2\sum_{j=1}^{n-1}f(r_j)+f(b)  \right ] 
\end{eqnarray}
$$



### (b) Use Composite Simpson's Rule:

$$
\int_{0}^{R}f(r)\mathrm{d}r=\frac{h}{3}\left [f(a)+2\sum_{j=1}^{(n/2)-1}f(r_{2j})+4\sum_{j=1}^{n/2}f(r_{2j-1})+f(b)  \right ] 
$$



### (c) Select a method of your choice that you think is more accurate than the two methods above:

​	Composite Simpson Rule is more accurate because the remainder term is $o(h^4)$ while that of the trapezoidal method is $o(h^2)$.		

## Python Code Output:

Composite Trapezoidal Rule= 0.1790632414322496 kg/s
Composite Simpson Rule= 0.18621350631182 kg/s

------

# Q.2 :

## Solution:

$$
C=1.83;t_0=1,k(t_0)=1,\epsilon(t_0)=0.2176
$$



### (a) Use the Euler method:

$$
\left\{  

\begin{array}{**lr**}  

 k(t_{i+1})=k(t_i)+h(-\epsilon(t_i))\\

\epsilon(t_{i+1})=\epsilon(t_i)+h(-C\frac{\epsilon(t_i)^2}{k(t_i)})

\end{array}  

\right.
$$

### (b) Use the modified Euler method:

$$
\left\{ \begin{eqnarray}
\bar{k}(t_{i+1}) & = & k(t_i)+h(-\epsilon (t_i))\\
\bar{\epsilon}(t_{i+1})&=&\epsilon(t_i)+h(-C\frac{\epsilon (t_i)^2}{k(t_i)} )\\
 {k}(t_{i+1})& = &k(t_i)+\frac{h}{2}(-\epsilon (t_i)-\bar{\epsilon}(t_{i+1}) )\\
\epsilon(t_{i+1})&=&\epsilon(t_i)+\frac{h}{2}\left ( -C\frac{\epsilon (t_i)^2}{k(t_i)} -
C\frac{\bar{\epsilon } (t_{i+1})^2}{\bar{k} (t_{i+1})} \right ) 

\end{eqnarray}\right.
$$

### (c) Use the 4th order Runge-Kutta method:

$$
\left\{\begin{eqnarray}
K_{k,1}&=&-\epsilon
\\K_{\epsilon ,1}&=&-C\frac{\epsilon ^2}{k} \\
K_{k,2}&=&-\epsilon+h\frac{K_{k,1}}{2} 
\\K_{\epsilon ,2}&=&-C\frac{\epsilon ^2}{k}+h\frac{K_{\epsilon,1}}{2}\\
K_{k,3}&=&-\epsilon+h\frac{K_{k,2}}{2}
\\K_{\epsilon ,3}&=&-C\frac{\epsilon ^2}{k}+h\frac{K_{\epsilon,2}}{2}\\
K_{k,4}&=&-\epsilon+hK_{k,3}
\\K_{\epsilon ,4}&=&-C\frac{\epsilon ^2}{k}+hK_{\epsilon,3}\\
\end{eqnarray}\right.
$$

### Comment:

1. RK4 Method converges the fastest following the Modified Euler Method, the Euler Method converges the slowest.
2. RK4 Method and Modified Euler Method has almost the same accuracy,while the Euler Method has the worst accuracy.

## Python Code Output:

Results of $k$ at $t=0.5s$ ,with different number of intervals dividing [1,5]

<img src="image-20250103113451689.png" alt="image-20250103113451689" style="zoom:33%;" />

|                 | Euler method | modified Euler method | Runge-Kutta 4 method |
| :-------------: | :----------: | :-------------------: | :------------------: |
| $k(t),t=5;N=40$ |   0.51674    |        0.51943        |       0.51939        |



# Q.3:

## Solution:

Let $y_+=\frac{u_\tau y}{\nu}, U_+=\frac{U}{u_\tau}$,
$$
f(u_\tau)=-y_{+}+U_{+}+e^{-\kappa B}\left[e^{\kappa U_{+}}-1-\kappa U_{+}-\frac{1}{2!}\left(\kappa U_{+}\right)^{2}-\frac{1}{3!}\left(\kappa U_{+}\right)^{3}-\frac{1}{4!}\left(\kappa U_{+}\right)^{4}\right]=0
$$

$$
\tau_{wall}={u_\tau}^2\rho
$$


### (a) Use bisection method

​	Since $f(0.1)=3.0524259021078516e+36>0,f(2)=-1318.893031144245<0$, let $a=0.1,b=2,mid=\frac{b+a}{2}$, if $f(a)f(mid)<0$, $b=mid$, else $a=mid$, repeat until $toler=\vert b-a\vert<10^{-3}$

### (b) Use Newton’s method

​	Take initial root approximation $p_0=1$, use $p_n=p_{n-1}-\frac{f(p_{n-1})}{f'(p_{n-1})}$ to find the approximation root. $f'(x_0)=\frac{f(x_{0}+h)-f(x_0)}{h},h=0.0001$

### (c) Try to find a fixed point iteration formula that converges

$$
f(U_+)=-y_{+}+U_{+}+e^{-\kappa B}\left[e^{\kappa U_{+}}-1-\kappa U_{+}-\frac{1}{2!}\left(\kappa U_{+}\right)^{2}-\frac{1}{3!}\left(\kappa U_{+}\right)^{3}-\frac{1}{4!}\left(\kappa U_{+}\right)^{4}\right]=0\\
U_+=U_+-\frac{f(U_+)}{U_+^4}
$$

is a fixed point iteration formula that converges.

## Python Code Output:

1. Output of Bisection Method
   Root of u_tau= **0.9975830078125001**
   Wall Shear Stress(bisection method)= **1.2439648218452932<img src="image-20250103150630698.png" alt="image-20250103150630698" style="zoom:33%;" />**

2. Output of the Newton's Method

   Root of u_tau= **0.9976727452597954**
   Wall Shear Stress(Newton's method)= **1.2441886332927707**<img src="image-20250103150637473.png" alt="image-20250103150637473" style="zoom:33%;" />

3. Output of a fixed point iteration that converges
   Wall Shear Stress(Fixed Point Iteration)= **1.3218985298910881** 
   Iterations= **321**
   Tolerance= **0.0009987081191376035**

# Q.4:

## Solution:

First use the Householder transformation transfer matrix A into a tridiagonal matrix and than use the Gram-Schmidt process to decompose A into QR,than let A=RQ, do the QR decomposition again until the sum of the non-diagonal element of A is no bigger than 1e-10.Than the diagonal elements are the eigenvalue of A.	

## Python Code Output:

Eigenvalues: [132.62787533  52.4423     -11.54113078  -3.52904455]
Number of positive eigenvalues: 2
Number of negative eigenvalues: 2

# Q.5:

## Solution:

$$
y_{+}=U_{+}+e^{-\kappa B}\left[e^{\kappa U_{+}}-1-\kappa U_{+}-\frac{1}{2!}\left(\kappa U_{+}\right)^{2}-\frac{1}{3!}\left(\kappa U_{+}\right)^{3}-\frac{1}{4!}\left(\kappa U_{+}\right)^{4}\right]\label{14}
$$

Use the fsolve function to calculate the value of $U_+$ at $y_+ = 1,590$.Use these value as an interval to plot equation $\ref{14}$ 



## Python Code Output:

1. Draw the $y_+$ vs $U_+$,Plotting the $y_+$ range using the log10 scale.

   <img src="clipboard-5921566-5921569.png" alt="clipboard" style="zoom:33%;" />

2. Construct piecewise linear interpolation, and draw on a figure: <img src="image-20250117201257892.png" alt="image-20250117201257892" style="zoom:%;" />
3. For each subinterval you choose, use the Monomial polynomials up to order 3 to construct least-squares approximation and plot the answers in a figure: The interval I choose is [1, 5, 25,100,300, 590]<img src="image-20250117201354205.png" alt="image-20250117201354205" style="zoom: 15%;" />
4. Repeat the previous step by using the Legendre polynomials up to order 3:<img src="image-20250117201418131.png" alt="image-20250117201418131" style="zoom:15%;" />

# Q.6:

## Solution:

### (a) Use N equally spaced points in [−1,1] to interpolate the Runge function.

Use Newton interpolation method to do the interpolation, use N point to divide the interval.

<img src="截屏2025-01-03 16.31.34.png" alt="截屏2025-01-03 16.31.34" style="zoom:33%;" />

### (b) Use roots of Legendre polynomial of degree (N−1), to interpolate the Runge function

 Use scipy and numpy to compute the roots of Legendre polynomial of degree (N−1) and do (a) again.

### (c) Plot

| N    |                             (c)                              |
| :--- | :----------------------------------------------------------: |
| N=15 | <img src="image-20250102102235703.png" alt="image-20250102102235703" style="zoom:33%;" /> |
| N=16 | <img src="image-20250102102358540.png" alt="image-20250102102358540" style="zoom:33%;" /> |
| N=17 | <img src="image-20250102102420548.png" alt="image-20250102102420548" style="zoom:33%;" /> |
| N=18 | <img src="image-20250102102437734.png" alt="image-20250102102437734" style="zoom:33%;" /> |
| N=19 | <img src="image-20250102102457806.png" alt="image-20250102102457806" style="zoom:33%;" /> |
| N=20 | <img src="image-20250102102517847.png" alt="image-20250102102517847" style="zoom:33%;" /> |

### (d) Use least square approximation to approximate the Runge function

​	Since the Legendre polynomials are orthogonal polynomials on the interval [-1, 1] with respect to the power function 1, which constitute a set of standard orthogonal bases of the polynomial space, the least squares approximation of the Runge function can be obtained by projecting the Runge function onto this set of bases.

$P_i$  is $i$-th order Legendre polynomials.$R(x)=\frac{1}{1+25x^2}$.The least square approximation of $R(x)$:
$$
\begin{eqnarray}
L(x)&=&\frac{\langle{P_0,R(x)}\rangle}{\langle{P_0,P_0}\rangle} P_0+
\frac{\langle{P_1,R(x)}\rangle}{\langle{P_1,P_1}\rangle} P_1+\cdots+
\frac{\langle{P_{N-1},R(x)}\rangle}{\langle{P_{N-1},P_{N-1}}\rangle} P_{N-1}\\
\langle{P_{N-1},R(x)}\rangle&=&\int_{-1}^{1}P_{N-1}R(x)\mathrm{d}x \label{least_square}
\end{eqnarray}
$$
Use sympy to calculate the integral $\eqref{least_square}$ analytically.

### (e) Plot

| N    |                             (e)                              |
| :--- | :----------------------------------------------------------: |
| N=15 | <img src="image-20250102103031040.png" alt="image-20250102103031040" style="zoom:33%;" /> |
| N=16 | <img src="image-20250102103122637.png" alt="image-20250102103122637" style="zoom:33%;" /> |
| N=17 | <img src="image-20250102103215418.png" alt="image-20250102103215418" style="zoom:33%;" /> |
| N=18 | <img src="image-20250102103419075.png" alt="image-20250102103419075" style="zoom:33%;" /> |
| N=19 | <img src="image-20250102103450800.png" alt="image-20250102103450800" style="zoom:33%;" /> |
| N=20 | <img src="image-20250102103523684.png" alt="image-20250102103523684" style="zoom:33%;" /> |

### (f) Now we use a reduced Gauss-Legendre quadrature and do Q.6d again.

​	Use Gauss-Legendre quadrature to compute the integral $\ref{least_square}$ numerically.

### (g) Plot

| N    |                             (g)                              |
| :--- | :----------------------------------------------------------: |
| N=15 | <img src="image-20250102103048265.png" alt="image-20250102103048265" style="zoom:33%;" /> |
| N=16 | <img src="image-20250102103145218.png" alt="image-20250102103145218" style="zoom:33%;" /> |
| N=17 | <img src="image-20250102103233065.png" alt="image-20250102103233065" style="zoom:33%;" /> |
| N=18 | <img src="image-20250102103432993.png" alt="image-20250102103432993" style="zoom:33%;" /> |
| N=19 | <img src="image-20250102103503319.png" alt="image-20250102103503319" style="zoom:33%;" /> |
| N=20 | <img src="image-20250102103534822.png" alt="image-20250102103534822" style="zoom:33%;" /> |

​	As we can see from the plot the result of  Legendre polynomials to interpolate a function is the same as using polynomials and reduced quadrature to approximate a function. 

