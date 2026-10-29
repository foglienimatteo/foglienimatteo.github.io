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

- J. L. M. Quiroz González, D. Thompson: _"Getting started with Numerov’s method"_, Computer in Physics, Vol. 11 Num. 5 pag. 514–515 (1997), DOI [10.1063/1.168593](https://doi.org/10.1063/1.168593)


## 1. Introduction

Numerov's method is an efficient algorithm for solving second-order differential equations of the following form:

$$
\boxed{ \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = U(x) + V(x)\,y} \tag{1.1}
$$

Numerov's method relies on two properties of (1):

1. the first derivative $dy/dx$ is absent;
2. $U(x) + V(x)\, y$ is linear wrt $y$.


Examples in physics of this type of equations:

- 1D time-independent Schrödinger equation. From $[-\frac{\hbar^2}{2m}\tfrac{d^2}{\mathrm{d}x^2} + V(x)]\psi = E\psi$:

$$
 \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2} = \frac{2m}{\hbar^2}\,\big[V(x) - E\big]\,\psi \tag{1.2}
$$


- Undamped forced harmonic oscillator ($x$ plays the role of time):

$$
m\, \frac{\mathrm{d}^2y}{\mathrm{d}x^2} = f_0\cos{(\omega x - k\,y)} \tag{1.3}
$$

- Poisson equation with spherically symmetric distribution in SI units:

$$
\nabla^2\phi = \frac{1}{r} \frac{\mathrm{d}^2(r\phi)}{\mathrm{d}r^2} = -\frac{\rho(r)}{\varepsilon_0} \tag{1.4}
$$

Cast into the form (1):

| equation                    | $y$                    | $U(x)$                            | $V(x)$ of Eq. (1)              |
| --------------------------- | ---------------------- | --------------------------------- | ------------------------------ |
| (2) Schrödinger             | $\psi$                 | $0$                               | $\frac{2m}{\hbar^2}[V(x) - E]$ |
| (3) forced oscillator       | $y$                    | $\frac{f_0}{m}\cos\omega x$       | $-\frac{k}{m}$                 |
| Poisson, spherical symmetry | $\varphi \equiv r\phi$ | $-\frac{r\rho(r)}{\varepsilon_0}$ | $0$                            |




## 2. Derivation

We now want to apply the central difference scheme to the Numerov's method of Eq.(1).

**Assumption**: We suppose to have a set of $N$ values $x_0, ..., x_N$ EQUALLY SPACED by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N-1$, and we define $y_n := y(x_n)$.

We can then apply the Taylor expansion to $y_{n\pm1}$, forward ($+$) and backward ($-$). Adding the two expansions, all odd powers of $h$ cancel (symmetric stencil), so the first neglected term is $O(h^8)$:

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

Unfortunately, Eq.(2.4) is not yet a usable scheme, we cannot compute $\mathrm{d}^2F/\mathrm{d}x^2$ analytically from its definiion:

$$
\begin{align*}
\frac{\mathrm{d}^2F}{\mathrm{d}x^2}  
    &= \frac{\mathrm{d}}{\mathrm{d}x} \frac{\mathrm{d}}{\mathrm{d}x}
        \left[U(x) + V(x)\, y(x) \right] \\[8pt]
    &=  \frac{\mathrm{d}}{\mathrm{d}x}
        \left[U' + V'\,y + V \, y'\right] \\[8pt]
    &= U'' + V''\,y + 2V'\,y' + V\,y'' \\[8pt]
    &= U'' + V''\,y + 2V'\,y' + V\,F \\[8pt]
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
(2.5c) \; \mathrm{into} \; (2.4)
\; \Rightarrow \quad 
y_{n+1} &= 2\,y_n - y_{n-1} + h^2\,F_n + \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + (12 - 2)\,F_n + F_{n-1}\big) + O(h^6) \\[8pt]
        &= 2\,y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} \tag{2.6}
$$

We didn't finish yet, because Eq.(2.6) is implicit: $F_{n+1} = U_{n+1} + V_{n+1}\,y_{n+1}$ contains the unknown $y_{n+1}$.
Spliting $F_{n+1}$ according to its definition Eq.(2.3), and collect the $y_{n+1}$ terms on the left:

$$
\begin{align*}
\Rightarrow \quad 
y_{n+1} &= 2y_n - y_{n-1} + \frac{h^2}{12}\left[U_{n+1} + V_{n+1}\,y_{n+1} + 10F_n + F_{n-1}\right]+ + O(h^6) \\[8pt]
 &= 2y_n - y_{n-1} + \frac{h^2}{12}\left[U_{n+1} + 10F_n + F_{n-1}\right]+ \frac{V_{n+1}h^2}{12}\,y_{n+1} + O(h^6) \\[8pt]
\end{align*}
$$


$$
\begin{aligned}\tag{2.7}
\Rightarrow \quad 
\left(1 - \frac{V_{n+1}h^2}{12}\right) y_{n+1} =  2y_n - y_{n-1} + \frac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} 
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
- Linearity in $y$ is used only in this step: for a general $F(x, y)$, Eq. (2.7) is a nonlinear equation for $y_{n+1}$ (root finding at every step).
- Eq. (2.8) is a two-step recursion: it needs $y_0$ and $y_1$ to start, which is the subject of the paper's Eqs. (8)–(15).

---

## Remark: leading local error (not in the paper)

Keep the $h^6$ terms of (b) and (d), and use $d^4F/dx^4 = y^{(6)}$:

$$
y_{n+1} - 2y_n + y_{n-1} = h^2F_n + \frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n + \frac{h^6}{360}\,y^{(6)}_n + O(h^8)
$$

$$
\frac{h^4}{12}\left. \frac{\mathrm{d}^2F}{\mathrm{d}x^2}\right|_n = \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) - \frac{h^6}{144}\,y^{(6)}_n + O(h^8)
$$

Hence the $O(h^6)$ in (g) and (7) is

$$
\left(\frac{1}{360} - \frac{1}{144}\right) h^6\,y^{(6)}_n + O(h^8) = -\frac{h^6}{240}\,y^{(6)}_n + O(h^8)
$$

