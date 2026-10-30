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

- J. L. M. Quiroz González, D. Thompson: _"Getting started with Numerov’s method"_, Computers in Physics, Vol. 11 Num. 5 pag. 514–515 (1997), DOI [10.1063/1.168593](https://doi.org/10.1063/1.168593)


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
F := U(x) + V(x)\,y \; , \tag{2.3}
$$

and on the grid $F_n := F(x_n) = U_n + V_n\,y_n$, with $U_n := U(x_n)$ and $V_n := V(x_n)$.

$$
(1.1) \; \Rightarrow \;  \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = F \; ,
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

because that requires derivatives of $U$ and $V$ and the slope $$y'$$, which the scheme does not carry.


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
(2.5c) \; \mathrm{into} \; (2.4)
\; \Rightarrow \quad 
y_{n+1} &= 2\,y_n - y_{n-1} + h^2\,F_n + \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + (12 - 2)\,F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} \tag{2.6}
$$

We didn't finish yet, because Eq.(2.6) is implicit: $F_{n+1} = U_{n+1} + V_{n+1}\,y_{n+1}$ contains the unknown $y_{n+1}$.
Splitting $F_{n+1}$ according to its definition Eq.(2.3), and collecting the $y_{n+1}$ terms on the left:

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

---

## Remark: leading local error (not in the paper)

Expand Eq.(2.1) and Eq.(2.5a) one order further, and use $\mathrm{d}^4F/\mathrm{d}x^4 = y^{(6)}$:

$$
y_{n+1} - 2y_n + y_{n-1} = h^2F_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + \frac{h^6}{360}\,y^{(6)}_n + O(h^8)
$$

$$
\frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n = \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) - \frac{h^6}{144}\,y^{(6)}_n + O(h^8)
$$

Hence the $O(h^6)$ in Eq.(2.6) and Eq.(2.8) is

$$
\left(\frac{1}{360} - \frac{1}{144}\right) h^6\,y^{(6)}_n + O(h^8) = -\frac{h^6}{240}\,y^{(6)}_n + O(h^8)
$$


---

## 3. Initial conditions

Eq.(2.8) is a two-step recursion: to compute $y_{n+1}$ we need the two previous values $y_n$ and $y_{n-1}$. The first step it can perform is $n=1$, i.e. to start we need both $y_0$ and $y_1$.

**Assumption**: As initial conditions we know the value $y_0 = y(x_0)$, which gives $F_0 = U_0 + V_0\,y_0$, and the slope $$y'_0 = y'(x_0)$$. Without loss of generality, $x_0 = 0$.

How accurate must $y_1$ be? To compute $y_2$ at $O(h^6)$ with Eq.(2.8) we would need $y_1$ with the same accuracy. However, it is sufficient to compute $y_1$ at $O(h^5)$: $y_1$ is computed only once, and an error $\delta$ on it propagates through the recursion growing at most linearly, $\sim N\,\delta = O(h^4)$ after $N \sim 1/h$ steps, which is not worse than the global error of the method.

**Note**: the paper states that the global error of Numerov's method is $O(h^5)$. It is $O(h^4)$ (Numerov's method is a fourth-order method): the local errors $O(h^6)$ of $N \sim 1/h$ steps accumulate through the second-order recursion as $N^2 \cdot O(h^6) = O(h^4)$, see the numerical check in Section 3.4. The conclusion that $y_1$ at $O(h^5)$ is sufficient still holds.
{: .notice--warning}


### 3.1 Taylor expansion of $y_1$

The best we can do from scratch, without any analytical differentiation, is to Taylor-expand $y_1$ about $x_0$, using $$y'' = F$$, $$y''' = F'$$ and $$y^{(4)} = F''$$:

$$
y_1 = y_0 + h\,y'_0 + \frac{h^2}{2}\,F_0 + \frac{h^3}{3!}\,F'_0 + \frac{h^4}{4!}\,F''_0 + O(h^5) \tag{3.1}
$$

The problem is then how to estimate $$F'_0$$ and $$F''_0$$, which (as for Eq.(2.4)) we cannot compute analytically.


### 3.2 A naive first step

A first attempt is to replace $$F'_0$$ with the forward difference and to drop the $$F''_0$$ term:

$$
F'_0 = \frac{F_1 - F_0}{h} + O(h) \tag{3.2}
$$

$$
\begin{align*}
(3.2) \; \mathrm{into} \; (3.1) \; \Rightarrow \quad
y_1 &= y_0 + h\,y'_0 + \frac{h^2}{2}\,F_0 + \frac{h^2}{6}\big(F_1 - F_0\big) + O(h^4) \\[8pt]
    &= y_0 + h\,y'_0 + \frac{h^2}{6}\big(F_1 + 2F_0\big) + O(h^4)
\end{align*}
$$

where both the $\frac{h^3}{6}\,O(h)$ error of Eq.(3.2) and the dropped $$\frac{h^4}{24}F''_0$$ term are $O(h^4)$.
As in Eq.(2.7), $F_1 = U_1 + V_1\,y_1$ contains the unknown $y_1$, and we collect it on the left:

$$
\Rightarrow \quad
y_1 = \frac{
    y_0 + h\,y'_0 + \dfrac{h^2}{6}\big(U_1 + 2F_0\big)
  }{
    1 - \dfrac{h^2\,V_1}{6}
  } + O(h^4) \tag{3.3}
$$

**Typo in the paper**: Eq.(10) of the paper prints $(U_1 - 2F_0)$ in the numerator; the correct sign is $(U_1 + 2F_0)$. With the minus sign, $y_1$ is only $O(h^2)$ accurate and the global error drops to $O(h)$ (Section 3.4).
{: .notice--warning}

Eq.(3.3) is one order short: $y_1$ at $O(h^4)$ costs a factor $h$ in the global error, which becomes $O(h^3)$.
The paper mentions three alternatives found in the literature:

- compute $$F''$$ analytically, as in Eq.(2.4), which is not always practical;
- compute $y_1$ from Eq.(3.1) without $$F'_0$$ and $$F''_0$$, get a first estimate of $y_2$ with Eq.(2.8), use it to estimate $$F''_0$$, and iterate until $y_1$ and $y_2$ stabilise;
- perform the first step with a different self-starting algorithm (e.g. Runge–Kutta), then switch to Numerov.

None of them is necessary.


### 3.3 The paper's starting formula

Suppose for a moment that we know $y_2$, hence $F_2$. Using $F_0$, $F_1$ and $F_2$ we can estimate $$F''_0$$, so we look for $y_1$ in the form:

$$
y_1 = y_0 + h\,y'_0 + h^2\big(a\,F_0 + b\,F_1 + c\,F_2\big) \tag{3.4}
$$

The constants $a$, $b$, $c$ are fixed by requiring that Eq.(3.4) agrees with Eq.(3.1). Taylor-expanding $F_1$ and $F_2$ about $x_0 = 0$:

$$
F_1 = F_0 + h\,F'_0 + \frac{h^2}{2}\,F''_0 + O(h^3) \; , \qquad
F_2 = F_0 + 2h\,F'_0 + 2h^2\,F''_0 + O(h^3)
$$

$$
\Rightarrow \quad
h^2\big(a\,F_0 + b\,F_1 + c\,F_2\big) = (a + b + c)\,h^2\,F_0 + (b + 2c)\,h^3\,F'_0 + \left(\frac{b}{2} + 2c\right) h^4\,F''_0 + O(h^5)
$$

and comparing term by term with Eq.(3.1):

$$
\begin{cases}
a + b + c = \dfrac{1}{2} \\[8pt]
b + 2c = \dfrac{1}{6} \\[8pt]
\dfrac{b}{2} + 2c = \dfrac{1}{24}
\end{cases}
\quad \Rightarrow \quad
\frac{b}{2} = \frac{1}{6} - \frac{1}{24} = \frac{1}{8}
\quad \Rightarrow \quad
b = \frac{1}{4} \; , \quad c = -\frac{1}{24} \; , \quad a = \frac{7}{24}
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

Eqs.(3.5) and (3.6) are two linear equations in the two unknowns $y_1$ and $y_2$. Splitting $F_1 = U_1 + V_1\,y_1$ and $F_2 = U_2 + V_2\,y_2$:

$$
\begin{align*}
(3.5) \; \Rightarrow \quad
& y_1 = y_0 + h\,y'_0 + \frac{h^2}{24}\big(7F_0 + 6U_1 + 6V_1\,y_1 - U_2 - V_2\,y_2\big) \\[8pt]
& \Rightarrow \quad
\left(1 - \frac{V_1h^2}{4}\right) y_1 + \frac{V_2h^2}{24}\,y_2 = y_0 + h\,y'_0 + \frac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \\[16pt]
(3.6) \; \Rightarrow \quad
& \left(1 - \frac{V_2h^2}{12}\right) y_2 = 2y_1 - y_0 + \frac{h^2}{12}\big(U_2 + 10U_1 + 10V_1\,y_1 + F_0\big) \\[8pt]
& \Rightarrow \quad
-\left(2 + \frac{5V_1h^2}{6}\right) y_1 + \left(1 - \frac{V_2h^2}{12}\right) y_2 = -y_0 + \frac{h^2}{12}\big(F_0 + 10U_1 + U_2\big)
\end{align*}
$$

i.e. the linear system

$$
\begin{cases}
a_{11}\,y_1 + a_{12}\,y_2 = b_1 \\[4pt]
a_{21}\,y_1 + a_{22}\,y_2 = b_2
\end{cases} \tag{3.7}
$$

$$
\begin{align*}
a_{11} &= 1 - \frac{V_1h^2}{4} \; , & a_{12} &= \frac{V_2h^2}{24} \; , & b_1 &= y_0 + h\,y'_0 + \frac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \; , \\[8pt]
a_{21} &= -2 - \frac{5V_1h^2}{6} \; , & a_{22} &= 1 - \frac{V_2h^2}{12} \; , & b_2 &= -y_0 + \frac{h^2}{12}\big(F_0 + 10U_1 + U_2\big) \; .
\end{align*}
$$

Solving with Cramer's rule, $y_1 = (a_{22}\,b_1 - a_{12}\,b_2)/(a_{11}\,a_{22} - a_{12}\,a_{21})$, where:

$$
\begin{align*}
a_{22}\,b_1 - a_{12}\,b_2
    &= y_0\left(1 - \frac{V_2h^2}{12} + \frac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \frac{V_2h^2}{12}\right) + \frac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) \\[8pt]
    &\quad - \frac{V_2h^4}{288}\Big[\big(7F_0 + 6U_1 - U_2\big) + \big(F_0 + 10U_1 + U_2\big)\Big] \\[8pt]
    &= y_0\left(1 - \frac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \frac{V_2h^2}{12}\right) \\[8pt]
    &\quad + \frac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) - \frac{V_2h^4}{36}\big(F_0 + 2U_1\big) \\[16pt]
a_{11}\,a_{22} - a_{12}\,a_{21}
    &= \left(1 - \frac{V_1h^2}{4}\right)\left(1 - \frac{V_2h^2}{12}\right) + \frac{V_2h^2}{24}\left(2 + \frac{5V_1h^2}{6}\right) \\[8pt]
    &= 1 - \frac{V_1h^2}{4} + V_1V_2h^4\left(\frac{1}{48} + \frac{5}{144}\right)
     = 1 - \frac{V_1h^2}{4} + \frac{V_1V_2h^4}{18}
\end{align*}
$$

$$
\boxed{
y_1 = \frac{
    \begin{aligned}
    & y_0\left(1 - \dfrac{V_2h^2}{24}\right) + h\,y'_0\left(1 - \dfrac{V_2h^2}{12}\right) \\
    & + \dfrac{h^2}{24}\big(7F_0 + 6U_1 - U_2\big) - \dfrac{h^4\,V_2}{36}\big(F_0 + 2U_1\big)
    \end{aligned}
  }{
    1 - \dfrac{V_1h^2}{4} + \dfrac{V_1V_2h^4}{18}
  } + O(h^5)} \tag{3.8}
$$

- Eq.(3.8) gives $y_1$ at $O(h^5)$ from the initial conditions $y_0$, $$y'_0$$ and from $U$, $V$ at $x_1$ and $x_2$ only: no iteration and no analytical derivatives are needed.
- Then $y_2$ follows from Eq.(3.6), and Eq.(2.8) continues for $n \geq 2$.


### 3.4 Numerical check (not in the paper)

Test problem: $V(x) = -(1 + x)$ and $$U(x) = y'' - V\,y$$ chosen so that the exact solution is $y(x) = e^{x/2}\cos(2x) + x^2$, integrated on $x \in [0, 2]$ with Eq.(2.8). Observed orders $p$ of the errors $\propto h^p$, for $h = 0.02,\ 0.01,\ 0.005$:

| starting value $y_1$           | error on $y_1$ | global error at $x = 2$ |
| ------------------------------ | -------------- | ----------------------- |
| paper's Eq.(10), as printed    | $h^2$          | $h^1$                   |
| Eq.(3.3), sign corrected       | $h^4$          | $h^3$                   |
| Eq.(3.8)                       | $h^5$          | $h^4$                   |

