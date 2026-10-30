---
title: "Numerov's method"
#date:
#excerpt_separator: "<!--more-->"
last_modified_at: 2026-10-14T16:20:02-05:00
categories:
  - Blog
tags:
  - Post Formats
  - math
  - physics
toc: true   # Table-Of-Content
toc_label: "Table of Contents"
toc_icon: "cog"
# link: https://github.com


---

# Numerov's Method

References: 

- [1] J. L. M. Quiroz González, D. Thompson: _"Getting started with Numerov's method"_, Computers in Physics, Vol. 11 Num. 5 pag. 514–515 (1997), DOI [10.1063/1.168593](https://doi.org/10.1063/1.168593)

## TL;DR

Equation to be solved numerically

$$
\boxed{ \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = U(x) + V(x)\,y} \tag{1.1}
$$

**Assumptions**: we suppose the following:

- to have a set of $N+1$ values $x_0, ..., x_N$ EQUALLY SPACED by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N-1$
- to know all the $U_n$ and $V_n$
- defining $y_n := y(x_n)$, to know $y_0$ and $y_0'$


$$
\boxed{
y_{n+1} = \frac{
    2y_n - y_{n-1} + \dfrac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big)
  }{
    1 - \dfrac{h^2 \,V_{n+1}}{12}
  } + O(h^6)} \tag{2.8}
$$

For the first step

$$
\boxed{
y_1 = \tfrac{
    y_0\left(1 - \tfrac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \tfrac{V_2h^2}{12}\right) \\
   + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) - \tfrac{h^4\,V_2}{36}\big(F_0 + 2U_1\big)
  }{
    1 - \tfrac{V_1h^2}{4} + \tfrac{V_1V_2h^4}{18}
  } + O(h^5)} \tag{3.8}
$$

The methods has a global error of $O(h^4)$

## 1. Introduction

Numerov's method is an efficient algorithm for solving second-order differential equations of the following form:

$$
\boxed{ \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = U(x) + V(x)\,y} \tag{1.1}
$$

Numerov's method relies on two properties of Eq.(1.1):

1. the first derivative $dy/dx$ is absent;
2. $U(x) + V(x)\, y$ is linear wrt $y$.


Examples in physics of this type of equations:

- 1D time-independent Schrödinger equation. From $[-\frac{\hbar^2}{2m}\tfrac{\mathrm{d}^2}{\mathrm{d}x^2} + V(x)]\psi = E\psi$:

$$
 \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} = \frac{2m}{\hbar^2}\,\big[V(x) - E\big]\,\psi \tag{1.2}
$$


- Undamped forced harmonic oscillator ($x$ plays the role of time):

$$
m\, \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = f_0\cos{(\omega x)} - k\,y \tag{1.3}
$$

- Poisson equation with spherically symmetric distribution in SI units:

$$
\nabla^2\phi = \frac{1}{r} \frac{\mathrm{d}^2(r\phi)}{\mathrm{d}r^2} = -\frac{\rho(r)}{\varepsilon_0} \tag{1.4}
$$

Cast into the form of Eq.(1.1):

| equation                    | $y$                    | $U(x)$                            | $V(x)$ of Eq.(1.1)             |
| --------------------------- | ---------------------- | --------------------------------- | ------------------------------ |
| (1.2) Schrödinger           | $\psi$                 | $0$                               | $\frac{2m}{\hbar^2}[V(x) - E]$ |
| (1.3) forced oscillator     | $y$                    | $\frac{f_0}{m}\cos\omega x$       | $-\frac{k}{m}$                 |
| (1.4) Poisson, spherical    | $\varphi \equiv r\phi$ | $-\frac{r\rho(r)}{\varepsilon_0}$ | $0$                            |




## 2. Derivation

We now want to apply the central difference scheme to Eq.(1.1).

**Assumption**: We suppose to have a set of $N+1$ values $x_0, ..., x_N$ EQUALLY SPACED by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N-1$, and we define $y_n := y(x_n)$.

We can then apply the Taylor expansion to $y_{n\pm1}$, forward ($+$) and backward ($-$):

$$
y_{n\pm1} = y(x_n \pm h) = 
    y_n \pm h\,y'_n + \frac{h^2}{2!}\,y''_n \pm \frac{h^3}{3!}\,y'''_n  + 
    \frac{h^4}{4!}\,y^{(4)}_n \pm \frac{h^5}{5!}\,y^{(5)}_n + O(h^6)
 \tag{2.1}
$$

Add the two expansions. All odd powers of $h$ cancel (symmetric stencil), so the first neglected term is $O(h^6)$:

$$
\Rightarrow \quad
y_{n+1} + y_{n-1} = 2\left(y_n + \frac{h^2}{2!}\,y''_n + \frac{h^4}{4!}\,y^{(4)}_n \right) + O(h^6) \tag{2.2a}
$$

$$
\Rightarrow \quad
y_{n+1} = - y_{n-1} + 2\,y_n + h^2 \,y''_n + \frac{h^4}{12}\,y^{(4)}_n + O(h^6) \tag{2.2b}
$$


Defining the right-hand side of Eq.(1.1), which is linear in $y$, as

$$
F := U(x) + V(x)\,y \; , \tag{2.3a}
$$

and on the grid $F_n := F(x_n) = U_n + V_n\,y_n$, with $U_n := U(x_n)$ and $V_n := V(x_n)$.

$$
(1.1) \; \Rightarrow \;  \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = F \; , \tag{2.3b}
$$


$$
\Rightarrow \quad
\left. \frac{\mathrm{d}^2y}{\mathrm{d}x^2}\right|_n = F_n \; , \qquad
\left. \frac{\mathrm{d}^4y}{\mathrm{d}x^4}\right|_n = \left. \frac{\mathrm{d}^2}{\mathrm{d}x^2}\!\left( \frac{\mathrm{d}^2y}{\mathrm{d}x^2}\right)\right|_n = \left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n
$$



$$
(2.2b) \; \Rightarrow \quad  y_{n+1} = - y_{n-1} + 2\,y_n + h^2F_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + O(h^6) \tag{2.4}
$$

Unfortunately, Eq.(2.4) is not yet a usable scheme, we cannot compute $\mathrm{d}^2F/\mathrm{d}x^2$ analytically from its definition:

$$
\begin{align*}
\frac{\mathrm{d}^2F}{\mathrm{d}x^2}  
    &= \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x}
        \left[U(x) + V(x)\, y(x) \right] \\[8pt]
    &=  \frac{\mathrm{d}}{\mathrm{d}x}
        \left[U' + V'\,y + V \, y'\right] \\[8pt]
    &= U'' + V''\,y + 2V'\,y' + V\,y'' \\[8pt]
    &= U'' + V''\,y + 2V'\,y' + V\,F
\end{align*}
$$

because that requires derivatives of $U$ and $V$ and the slope $y'$, which the scheme does not carry.


However, $F = \mathrm{d}^2y/\mathrm{d}x^2$ is itself smooth, so Eq.(2.1)–(2.2b) hold with $y \to F$:

$$
\Rightarrow \quad
F_{n+1} - 2 F_n + F_{n-1} = h^2\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^4F}{\mathrm{d}x^4}\right|_n + O(h^6) \tag{2.5a}
$$

$$
\Rightarrow \quad 
\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n = \frac{F_{n+1} - 2F_n + F_{n-1}}{h^2} + O(h^2) \tag{2.5b}
$$


$$
\Rightarrow \quad 
\frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n = \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \tag{2.5c}
$$

We then substitute Eq.(2.5c) into Eq.(2.4):

$$
\begin{aligned}
(2.5c) \; &\mathrm{into} \; (2.4) \\[7pt]
\; \Rightarrow \quad 
y_{n+1} &= 2\,y_n - y_{n-1} + h^2\,F_n + \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + (12 - 2)\,F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} \tag{2.6}
$$

We didn't finish yet, because Eq.(2.6) is implicit: $F_{n+1} = U_{n+1} + V_{n+1}\,y_{n+1}$ contains the unknown $y_{n+1}$.
Splitting $F_{n+1}$ according to its definition Eq.(2.3a), and collecting the $y_{n+1}$ terms on the left:

$$
\begin{align*}
\Rightarrow \quad 
y_{n+1} &= 2y_n - y_{n-1} + \frac{h^2}{12}\left[U_{n+1} + V_{n+1}\,y_{n+1} + 10F_n + F_{n-1}\right] + O(h^6) \\[8pt]
 &= 2y_n - y_{n-1} + \frac{h^2}{12}\left[U_{n+1} + 10F_n + F_{n-1}\right]+ \frac{V_{n+1}h^2}{12}\,y_{n+1} + O(h^6)
\end{align*}
$$


$$
\Rightarrow \quad 
\left(1 - \frac{V_{n+1}h^2}{12}\right) y_{n+1} =  2y_n - y_{n-1} + \frac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big) + O(h^6) \tag{2.7}
$$

Divide by $1 - V_{n+1}h^2/12 = 1 + O(h^2)$, which is non-zero for $V_{n+1}h^2 \neq 12$ and leaves the remainder $O(h^6)$:

$$
\boxed{
y_{n+1} = \frac{
    2y_n - y_{n-1} + \dfrac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big)
  }{
    1 - \dfrac{h^2 \,V_{n+1}}{12}
  } + O(h^6)} \tag{2.8}
$$

- The right-hand side needs only $y_n$, $y_{n-1}$ (hence $F_n$, $F_{n-1}$, already computed in the previous steps) and $U_{n+1}$, $V_{n+1}$: one evaluation of $U$ and $V$ per step, for a local error $O(h^6)$.
- Linearity in $y$ is used only in this step: for a general $F(x, y)$, Eq.(2.6) is a nonlinear equation for $y_{n+1}$ (root finding at every step).
- Eq. (2.8) is a two-step recursion: it needs $y_0$ and $y_1$ to start, which is the subject of Section 3 (the paper's Eqs. (8)–(15)).



### 2.1 Leading local error term

We can make one step further to get the leading error term of Eq.(2.8).
Expanding Eq.(2.1) and Eq.(2.5a) two orders further,

$$
\begin{align*} \tag{2.9a}
y_{n\pm1} = y(x_n \pm h) = y_n &\pm   h\,y'_n + \frac{h^2}{2!}\,y''_n 
    \pm \frac{h^3}{3!}\,y'''_n  + \frac{h^4}{4!}\,y^{(4)}_n \\[8pt] 
    &\pm \frac{h^5}{5!}\,y^{(5)}_n + \frac{h^6}{6!}\,y^{(6)}_n 
    \pm \frac{h^7}{7!}\,y^{(7)}_n + O(h^8)
\end{align*}
$$

$$
\Rightarrow \quad
y_{n+1} + y_{n-1} = 2\left(y_n + \frac{h^2}{2!}\,y''_n + \frac{h^4}{4!}\,y^{(4)}_n +
\frac{h^6}{6!}\,y^{(6)}_n \right) + O(h^8) \tag{2.9b}
$$

$$
\Rightarrow \quad
y_{n+1} - 2\,y_n + y_{n-1} = h^2 \,y''_n + \frac{h^4}{12}\,y^{(4)}_n 
 + \frac{h^6}{360}\,y^{(6)}_n + O(h^8) \tag{2.9c}
$$

using the definition of $F$ in Eq.(2.3b),

$$
\Rightarrow \quad
y_{n+1} - 2y_n + y_{n-1} = h^2F_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + \frac{h^6}{360}\,y^{(6)}_n + O(h^8)
 \tag{2.10}
$$

rederiving Eq.(2.5c) from Eq.(2.5a) but keeping one step further,

$$
(2.5a) \; : \quad 
F_{n+1} - 2 F_n + F_{n-1} = h^2\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^4F}{\mathrm{d}x^4}\right|_n + O(h^6)
$$

$$
\Rightarrow \quad 
\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n = \frac{F_{n+1} - 2F_n + F_{n-1}}{h^2} - 
\frac{h^4}{12}\left. \frac{\mathrm{d}^4F}{\mathrm{d}x^4}\right|_n+ O(h^4)
\tag{2.11a}
$$


$$
\begin{align*}
\Rightarrow \quad 
\frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n &= \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) - \frac{h^4}{144}\left. \frac{\mathrm{d}^4F}{\mathrm{d}x^4}\right|_n + O(h^8) \\[10pt]
&= \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) - \frac{h^6}{144}\,y^{(6)}_n + O(h^8)
\tag{2.11b}
\end{align*}
$$


and comparing Eq.(2.11b) with (2.9c) we get that the $O(h^6)$ in Eq.(2.6) and Eq.(2.8) is:

$$
O(h^6)\; \mathrm{in \; Eq.(2.6-2.8)} = \left(\frac{1}{360} - \frac{1}{144}\right) h^6\,y^{(6)}_n + O(h^8) = -\frac{h^6}{240}\,y^{(6)}_n + O(h^8)
$$

## 2.2 Global Error

At first glance, one could say that the global error of Numerov's method is $O(h^5)$:
- to reach a fixed $x_n$, you take $n$ steps;
- if each step adds an error $\tau \propto h^6$,
- and if the errors just pile up over all the $N \propto 1/h$ points,
- the total error is $N \cdot \tau  \propto  (1/h) \cdot h^6  =  h^5$

This counting rule works for one-step methods for first-order ODEs, like Euler or Runge–Kutta. There the step error lands directly on the value $y$, so the global error is the local error divided by $h$.


However, Numerov is a two-step recursion for a second-order equation. Its local error doesn't land on $y_n$, but on the second difference, the discrete $y''$.
It's like an error in acceleration rather than in position, and such an error in acceleration gets integrated twice before it shows up in position.
In other words, one factor of h is lost per "integration": $h^6$ becomes $h^5$ for the slope, then $h^4$ for the value.


Redoing the original argument to get Numerov's method global error s $O(h^4)$ :
- to reach a fixed $x_n$, you take $n$ steps;
- if each step adds an error $\tau \propto h^6$,
- and if the errors influence quadratically over all the $N \propto 1/h$ points
- the total error is $N^2 \cdot \tau  \propto  (1/h^2) \cdot h^6  =  h^4$


## 3. Initial conditions

Eq.(2.8) is a two-step recursion: to compute $y_{n+1}$ we need the two previous values $y_n$ and $y_{n-1}$. Tha means, to start the algorithm we need both $y_0$ and $y_1$.

**Assumption**: As initial conditions we know the value $y_0 = y(x_0)$, which gives $F_0 = U_0 + V_0\,y_0$, and the slope $y'_0 = y'(x_0)$. Without loss of generality, $x_0 = 0$.

How accurate must $y_1$ be? To compute $y_2$ at $O(h^6)$ with Eq.(2.8), we would need $y_1$ with the same accuracy. However, it is sufficient to compute $y_1$ with a precision $O(h^5)$: $y_1$ is computed only once, and its error would be still smaller than the global one.

The best we can do from scratch, without any analytical differentiation, is to Taylor-expand $y_1$ about $x_0$, using $y'' = F$, $y''' = F'$ and $y^{(4)} = F''$:

$$
\Rightarrow \quad
y_1 = y_0 + h\,y'_0 + \frac{h^2}{2}\,F_0 + \frac{h^3}{3!}\,F'_0 + \frac{h^4}{4!}\,F''_0 + O(h^5) \tag{3.1}
$$

The problem is then: how to estimate $F'_0$ and $F_0''$? 


### 3.1 A naive first step

A first attempt is to replace $F'_0$ with the forward difference and to drop the $F_0''$ term by including it into the error (which goes from $O(h^5)$ to $O(h^4)$):

$$
F'_0 = \frac{F_1 - F_0}{h} + O(h) \tag{3.2}
$$

$$
\begin{align*}
(3.2) \; \mathrm{into} \; &(3.1)\; \mathrm{+ \;removing \;}F_0'' \ \\[8pt]
\Rightarrow \quad 
y_1 &= y_0 + h\,y'_0 + \frac{h^2}{2}\,F_0 + \frac{h^3}{6}\left[\frac{F_1 - F_0}{h} + O(h)\right] + 
    \cancel{\frac{h^4}{4!}\,F''_0} + O(h^{\cancel{5}\,4}) \\[8pt]
    &= y_0 + h\,y'_0 + \frac{h^2}{2}\,F_0 + \frac{h^2}{6}\big(F_1 - F_0\big) + O(h^4) +
    + O(h^4) \\[8pt]
    &= y_0 + h\,y'_0 + \frac{h^2}{6}\big(F_1 + 2F_0\big) + O(h^4)
\end{align*}
$$

As in Eq.(2.6), $F_1 = U_1 + V_1\,y_1$ contains the unknown $y_1$, and we collect it on the left:

$$
\begin{align*}
\Rightarrow \quad 
y_1 &= y_0 +  h\,y_0' + \frac{h^2}{6}(U_1 + V_1\,y_1 + 2 F_0)+ O(h^4) \\[8pt]
 &= y_0 +  h\,y_0' + \frac{h^2}{6}(U_1 + 2 F_0 ) + \frac{V_1\,h^2}{6}\,y_1 + O(h^4)
\end{align*}
$$


$$
\Rightarrow \quad 
\left(1 - \frac{V_1\,h^2}{6}\right) y_1 =  y_0 +  h\,y_0' + \frac{h^2}{6}\left(U_1 + 2F_0 \right) + O(h^4)
$$

$$
\Rightarrow \quad
y_1 = \frac{
    y_0 + h\,y'_0 + \dfrac{h^2}{6}\big(U_1 + 2F_0\big)
  }{
    1 - \dfrac{h^2\,V_1}{6}
  } + O(h^4) \tag{3.3}
$$



Eq.(3.3) is one order short: $y_1$ at $O(h^4)$ costs a factor $h$ in the global error, which becomes $O(h^3)$.
In the literature, three main alternatives can be found:

1. compute $F''$ analytically,which is not always practical;
2. compute $y_1$ from Eq.(3.1) without $F_0'$ and $F_0''$, get a first estimate of $y_2$ with Eq.(2.8), use it to estimate $F_0''$, and iterate until $y_1$ and $y_2$ stabilise;
3. perform the first step with a different self-starting algorithm (e.g. Runge–Kutta), then switch to Numerov.

None of them is necessary.


### 3.2 The $y_1$ starting formula

In [1], the authors find a better approach.

Suppose for a moment that we somehow obtained an estimation for $y_2$, hence also for $F_2=U_2 + V_2 \, y_2$. Using $F_0$, $F_1$ and $F_2$ we can estimate, via e.g. forward difference, $F_0''$.
That means, we look for $y_1$ in the form:

$$
y_1 = y_0 + h\,y'_0 + h^2\big(a\,F_0 + b\,F_1 + c\,F_2\big) \tag{3.4a}
$$

The constants $a$, $b$, $c$ are fixed by requiring that Eq.(3.4) agrees with Eq.(3.1). 

We expand via Taylor series $F_1$ and $F_2$ about $x_0$,

$$
\begin{align*}
F_1 &= F(x_0+h) = F_0 + h\,F'_0 + \frac{h^2}{2}\,F''_0 + O(h^3) \\[8pt]
F_2 &= F(x_0 + 2h) = F_0 + 2h\,F'_0 + 2h^2\,F''_0 + O(h^3)
\end{align*}
$$

which means that

$$
\begin{align*}
\Rightarrow \quad
h^2&\big(a\,F_0 + b\,F_1 + c\,F_2\big) = \\[6pt]
&= a\,h^2 F_0 + b\,h^2\left[F_0 + h\,F'_0 + \tfrac{h^2}{2}\,F''_0 + O(h^3)\right] +  \\[6pt]
&\quad\quad + c\,h^2\left[F_0 + 2h\,F'_0 + 2h^2\,F''_0 + O(h^3)\right]\\[6pt]
&=(a + b + c)\,h^2\,F_0 + (b + 2c)\,h^3\,F'_0 + \left(\tfrac{b}{2} + 2c\right) h^4\,F''_0 + O(h^5)
\end{align*}
$$

This gives from Eq.(3.1)

$$
\begin{align*}\tag{3.4b}
\Rightarrow \quad 
y_1 &= y_0 + h\,y'_0 + (a + b + c)\,h^2\,F_0  \\[6pt] 
&\quad\quad + (b + 2c)\,h^3\,F'_0 + \left(\tfrac{b}{2} + 2c\right) h^4\,F''_0 + O(h^5)

\end{align*}
$$

and comparing term by term Eq.(3.4b) with Eq.(3.1):

$$
\begin{cases}
a + b + c = \tfrac{1}{2} \\[6pt]
b + 2c = \tfrac{1}{6} \\[6pt]
\tfrac{b}{2} + 2c = \tfrac{1}{24}
\end{cases}
\; \Rightarrow \;
\tfrac{b}{2} = \tfrac{1}{6} - \tfrac{1}{24} = \tfrac{1}{8}
$$

$$
\quad \Rightarrow \quad
a = \tfrac{7}{24} \; ,\quad b = \tfrac{1}{4} \; , \quad c = -\tfrac{1}{24}
$$

$$
\Rightarrow \quad
y_1 = y_0 + h\,y'_0 + \frac{h^2}{24}\big(7F_0 + 6F_1 - F_2\big) + O(h^5) \tag{3.5}
$$

$F_2$ needs $y_2$, which is given by the standard Numerov step Eq.(2.8) with $n = 1$:

$$
y_2 = \frac{
    2y_1 - y_0 + \dfrac{h^2}{12}\big(U_2 + 10F_1 + F_0\big)
  }{
    1 - \dfrac{h^2\,V_2}{12}
  } + O(h^6) \tag{3.6}
$$

Eqs.(3.5) and (3.6) are two linear equations in the two unknowns $y_1$ and $y_2$. Splitting $F_1 = U_1 + V_1\,y_1$ and $F_2 = U_2 + V_2\,y_2$, we get

$$
\begin{align*}
&(3.5) \; \Rightarrow \;
 y_1 = y_0 + h\,y'_0 + \tfrac{h^2}{24}\big(7F_0 + 6U_1 + 6V_1\,y_1 - U_2 - V_2\,y_2\big) \\[6pt]
&\quad\Rightarrow \;
\left(1 - \tfrac{V_1h^2}{4}\right) y_1 + \tfrac{V_2h^2}{24}\,y_2 = y_0 + h\,y'_0 + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \\[16pt]
&(3.6) \; \Rightarrow \;
 \left(1 - \tfrac{V_2h^2}{12}\right) y_2 = 2y_1 - y_0 + \tfrac{h^2}{12}\big(U_2 + 10U_1 + 10V_1\,y_1 + F_0\big) \\[8pt]
&\quad \Rightarrow \;
-\left(2 + \tfrac{5V_1h^2}{6}\right) y_1 + \left(1 - \tfrac{V_2h^2}{12}\right) y_2 = -y_0 + \tfrac{h^2}{12}\big(F_0 + 10U_1 + U_2\big)
\end{align*}
$$

i.e. we have to solve the linear system

$$
\begin{cases}
a_{11}\,y_1 + a_{12}\,y_2 = b_1 \\[4pt]
a_{21}\,y_1 + a_{22}\,y_2 = b_2
\end{cases} \tag{3.7}
$$

$$
\begin{align*}
a_{11} &= 1 - \tfrac{V_1h^2}{4} \, , & a_{12} &= \tfrac{V_2h^2}{24} \, , & b_1 &= y_0 + h\,y'_0 + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \, , \\[8pt]
a_{21} &= -2 - \tfrac{5V_1h^2}{6} \, , & a_{22} &= 1 - \tfrac{V_2h^2}{12} \, , & b_2 &= -y_0 + \tfrac{h^2}{12}\big(F_0 + 10U_1 + U_2\big) \, .
\end{align*}
$$

Solving with Cramer's rule:

$$
y_1 = \frac{a_{22}\,b_1 - a_{12}\,b_2}{a_{11}\,a_{22} - a_{12}\,a_{21}}
$$

where:

$$
\begin{align*}
a_{22}\,&b_1 - a_{12}\,b_2 =\\[6pt]
    &= \left(1 - \tfrac{V_2h^2}{12}\right) \left[
            y_0 + h\,y'_0 + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big)
        \right]  \\[6pt]
        &\quad\quad 
        - \tfrac{V_2h^2}{24} \left[
        -y_0 + \tfrac{h^2}{12}\big(F_0 + 10U_1 + U_2\big)
        \right]\\[6pt]
    &= y_0\left(1 - \tfrac{V_2h^2}{12} + \tfrac{V_2h^2}{24}\right) 
        + h\,y'_0\left(1 - \tfrac{V_2h^2}{12}\right) 
        + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \\[8pt]
        &\quad\quad - \tfrac{V_2h^4}{288}\Big[\big(7F_0 + 6U_1 - \cancel{U_2}\big) 
        + \big(F_0 + 10U_1 + \cancel{U_2}\big)\Big] \\[8pt]
    &= y_0\left(1 - \tfrac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \tfrac{V_2h^2}{12}\right) \\[8pt]
    &\quad\quad + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) - \tfrac{V_2h^4}{36}\big(F_0 + 2U_1\big) \\[16pt]
a_{11}\,&a_{22} - a_{12}\,a_{21} = \\[6pt]
    &= \left(1 - \tfrac{V_1h^2}{4}\right)\left(1 - \tfrac{V_2h^2}{12}\right) + \tfrac{V_2h^2}{24}\left(2 + \tfrac{5V_1h^2}{6}\right) \\[8pt]
    &= 1 - \tfrac{V_1h^2}{4} + V_1V_2h^4\left(\tfrac{1}{48} + \tfrac{5}{144}\right) \\[6pt]
    &= 1 - \tfrac{V_1h^2}{4} + \tfrac{V_1V_2h^4}{18}
\end{align*}
$$

We get then the final expression:

$$
\boxed{
y_1 = \tfrac{
    y_0\left(1 - \tfrac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \tfrac{V_2h^2}{12}\right) \\
   + \tfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) - \tfrac{h^4\,V_2}{36}\big(F_0 + 2U_1\big)
  }{
    1 - \tfrac{V_1h^2}{4} + \tfrac{V_1V_2h^4}{18}
  } + O(h^5)} \tag{3.8}
$$

Eq.(3.8) gives $y_1$ at $O(h^5)$ from the initial conditions $y_0$, $y_0'$ and from $U$, $V$ at $x_1$ and $x_2$ only: no iteration and no analytical derivatives are needed.

Then, $y_2$ follows from Eq.(3.6), and Eq.(2.8) continues for $n \geq 2$.


### 3.3 Numerical check

Test problem: $V(x) = -(1 + x)$ and $U(x) = y'' - V\,y$ chosen so that the exact solution is $y(x) = e^{x/2}\cos(2x) + x^2$, integrated on $x \in [0, 2]$ with Eq.(2.8). Observed orders $p$ of the errors $\propto h^p$, for $h = 0.02,\ 0.01,\ 0.005$:

| starting value $y_1$           | error on $y_1$ | global error at $x = 2$ |
| ------------------------------ | -------------- | ----------------------- |
| paper's Eq.(10), as printed    | $h^2$          | $h^1$                   |
| Eq.(3.3), sign corrected       | $h^4$          | $h^3$                   |
| Eq.(3.8)                       | $h^5$          | $h^4$                   |

