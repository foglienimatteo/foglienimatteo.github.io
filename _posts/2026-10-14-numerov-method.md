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

# Numerov's method — step-by-step derivation of Eqs. (1)–(7)

Source: J. L. M. Quiroz González and D. Thompson (1997), *Getting started with Numerov's method*,
Computers in Physics 11, 514–515, DOI: [10.1063/1.168593](https://doi.org/10.1063/1.168593).

Notation is the paper's. Equation numbers (1)–(7) are the paper's; intermediate steps are labelled (a)–(i).

- uniform grid of step $h$: $x_{n\pm1} = x_n \pm h$;
- values on the grid: $y_n = y(x_n)$, $U_n = U(x_n)$, $V_n = V(x_n)$, $F_n = F(x_n)$;
- derivatives at $x_n$: $y'_n$, $y''_n$, $y'''_n$, $y^{(k)}_n$ stand for $d^ky/dx^k$ evaluated at $x = x_n$
  (the paper writes the same thing as $\left.\frac{d^2F}{dx^2}\right|_n$ in Eq. (6), and with primes, $y'_0$, $F''_0$, from Eq. (8) on).

---

## 0. Basics of numerical analysis


### Taylor's theorem

**Taylor's theorem** Given a function ${\displaystyle f:\mathbb{R} \to \mathbb{R} }$ that is ${N \geq 1}$ times differentiable at the point ${a\in \mathbb {R} }$, there exists a function ${\displaystyle g_{N}:\mathbb {R} \to \mathbb {R} }$ such that

$$ 
f(x) = \sum_{k=0}^{N}{\frac {f^{(k)}(a)}{k!}}(x-a)^{k} + g_{N}(x)(x-a)^{k} \; , \quad \quad \lim _{x\to a}g_{N}(x)=0.
$$

The polynomial appearing in Taylor's theorem is the ${\textstyle {\boldsymbol {N}}}$-th order Taylor polynomial of the function ${\displaystyle f}$ at the point ${\displaystyle a}$

$$
P_{N}(x) = f(a)+f'(a)(x-a)+{\frac {f''(a)}{2!}}(x-a)^{2}+\cdots +{\frac {f^{(N)}(a)}{N!}}(x-a)^{N}
$$

Using the little-o notation, the statement in Taylor's theorem reads as

$$
f(x) = \sum_{k=0}^{N}{\frac {f^{(k)}(a)}{k!}}(x-a)^{k} + o(\lvert x-a \rvert^N)
$$

We can rewrite the theorem in a slightly different notation, with the substitutions $x \rightarrow x+h$ and $a \rightarrow x$

$$
\begin{align*}
f(x+h)  &= \sum_{k=0}^{N}{\frac {h^k}{k!}}f^{(k)}(x) + o(h^N) \\[15pt]
        &= f(x) + o(1) \\[10pt]
        &= f(x) + h \, f'(x) + o(h) \\[10pt]
        &= f(x) + h \, f'(x) + \frac{h^2}{2} f''(x) + o(h^2) \\[10pt]
        &= ...
\end{align*}
$$


### First order finite differences

The derivative of a function $f$ at a point $x$ is defined by the limit

$$
f'(x) := \lim_{h \rightarrow 0^{+}} \frac{f(x+h)-f(x)}{h}
$$

We define the forward difference ${\Delta_{h}[f]}$, backward difference ${\nabla_{h}[f]}$ and central difference ${\delta_{h}[f]}$ as:

$$
\begin{align*}
    \mathrm{Forward  \; difference:} \quad &\Delta_{h}[f](x) := f(x+h)-f(x) \\[10pt]
    \mathrm{Backward \; difference:} \quad &\nabla_{h}[f](x) := f(x)-f(x-h) \\[10pt]
    \mathrm{Central  \; difference:} \quad &\delta_{h}[f](x) := f\left(x+\frac{h}{2}\right)-f\left(x-\frac{h}{2}\right) \\[15pt]

    \Rightarrow \quad \nabla_{h}[f](x) &= \Delta_{h}[f](x-h) \\[15pt]
    \Rightarrow \quad \delta_{h}[f](x) &= \Delta_{\tfrac {h}{2}}[f](x) + \nabla_{\tfrac {h}{2}}[f](x) \\[15pt]
\end{align*}
$$

We can approximate the derivative of a function (assuming is twice differentiable, three times for centered-difference) using each of the previous three differences. Applying Taylor's theorem, we can get the error of each approximation:

$$
\begin{align*}
    \frac{\Delta_h[f](x)}{h} - f'(x) 
        &= \frac{f(x+h) - f(x)}{h} - f'(x)  \\[10pt]
        &= \frac{1}{h}\bigg\{
            \left[\cancel{f(x)} + h \, f'(x) + o(h) \right] - \cancel{f(x)}
            \bigg\} - f'(x)  \\[10pt]
        &= o(h) \\[15pt]    
    \frac{\nabla_h[f](x)}{h} - f'(x) 
        &= \frac{f(x) - f(x-h)}{h} - f'(x)  \\[10pt]
        &= \frac{1}{h}\bigg\{
                \cancel{f(x)} - \left[\cancel{f(x)} - h \, f'(x) + o(h)\right]\bigg\} - f'(x)  \\[10pt]
        &= o(h) \\[15pt]
    \frac{\delta_h[f](x)}{h} - f'(x) 
        &= \frac{1}{h}\left[f\left(x+\frac{h}{2}\right) - f\left(x-\frac{h}{2}\right)\right] - f'(x)  \\[10pt]
        &= \frac{1}{h}\bigg\{
            \left[\cancel{f(x)} + \frac{h}{2} \, f'(x) + \cancel{\frac{h^2}{8} f''(x)} + o(h^2)\right]-\\[15pt]
            &\quad\quad\;
            \left[\cancel{f(x)} - \frac{h}{2} \, f'(x) + \cancel{\frac{h^2}{8} f''(x)} + o(h^2)\right]
                \bigg\} - f'(x)  \\[10pt]
        &= \frac{1}{h} \left[ 2\,\frac{h}{2} \, f'(x) + o(h^2)\right] - f'(x)\\[15pt]
        &= o(h^2) \\[20pt]
\end{align*}
$$


$$
\begin{align*}
 \quad f'(x) &= \frac{f(x+h) - f(x)}{h} + o(h) \\[10pt]
\Rightarrow \quad f'(x) &= \frac{f(x) - f(x-h)}{h} + o(h) \\[10pt]
 \quad f'(x) &= \frac{1}{h}\left[f\left(x+\frac{h}{2}\right) - f\left(x-\frac{h}{2}\right)\right] + o(h^2)
\end{align*}
$$

Rewritten using the index notation $x_n$ and $f_n = f(x_n)$ and assuming that all the elements are equally spaced by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N$

$$
\begin{align*}
 \quad f_{n}' &= \frac{f_{n+1} - f_n}{h} + o(h) \\[10pt]
\Rightarrow \quad f_{n}' &= \frac{f_{n} - f_{n-1}}{h} + o(h) \\[10pt]
\end{align*}
$$



### Second order fainite differences

In an analogous way, one can obtain finite difference approximations to higher order derivatives and differential operators. 
For example, by using the above central difference formula for $f'(x + ⁠h/2⁠)$ and $f'(x − ⁠h/2⁠)$ and applying a central difference formula for the derivative of $f'$ at ${x}$, we obtain the central difference approximation of the second derivative of ${f}$ :


$$
f'\left(x+\frac{h}{2}\right) = \frac{f(x+h) - f(x)}{h} + o(h^2)\\[12pt]
f'\left(x-\frac{h}{2}\right) = \frac{f(x) - f(x-h)}{h} + o(h^2) \\[18pt]
$$

$$
\begin{align*}
\Rightarrow f''(x) &= 
    \frac{1}{h}\left[f'\left(x+\frac{h}{2}\right) - f'\left(x-\frac{h}{2}\right)\right] + o(h^2) \\[15pt]
    &=\frac{1}{h}\left[
        \frac{f(x+h) - f(x)}{h} - \frac{f(x) - f(x-h)}{h}
    \right] + o(h^2) \\[15pt]
    &= \frac{f(x+h) - 2\,f(x) + f(x-h)}{h^2} 
     + o(h^2) \\[15pt]
\end{align*}
$$
## 1. Class of equations — Eqs. (1)–(3)

Numerov's method is an efficient algorithm for solving second-order differential equations of the following form:

$$
\frac{d^2y}{dx^2} = U(x) + V(x)\,y \tag{1}
$$

Numerov's method relies on two properties of (1):

1. the first derivative $dy/dx$ is absent;
2. $U(x) + V(x)y$ is linear wrt $y$.

The paper calls $U(x) + V(x)y$ the "LHS" of Eq. (1); as Eq. (1) is printed, it is the right-hand side.

Examples given in the paper:

- 1D time-independent Schrödinger equation. From $-\frac{\hbar^2}{2m}\psi'' + V(x)\psi = E\psi$:

$$
\frac{d^2\psi}{dx^2} = \frac{2m}{\hbar^2}\,\big[V(x) - E\big]\,\psi \tag{2}
$$

> The paper prints $2m/\hbar$: typo, the prefactor is $2m/\hbar^2$.

- Undamped forced harmonic oscillator ($x$ plays the role of time):

$$
m\,\frac{d^2y}{dx^2} = f_0\cos\omega x - k\,y \tag{3}
$$

Cast into the form (1):

| equation                    | $y$                    | $U(x)$                            | $V(x)$ of Eq. (1)              |
| --------------------------- | ---------------------- | --------------------------------- | ------------------------------ |
| (2) Schrödinger             | $\psi$                 | $0$                               | $\frac{2m}{\hbar^2}[V(x) - E]$ |
| (3) forced oscillator       | $y$                    | $\frac{f_0}{m}\cos\omega x$       | $-\frac{k}{m}$                 |
| Poisson, spherical symmetry | $\varphi \equiv r\phi$ | $-\frac{r\rho(r)}{\varepsilon_0}$ | $0$                            |

In the Schrödinger row, $V(x)$ inside the bracket is the potential energy (the paper uses the same symbol).
The Poisson case is only mentioned in words in the paper; with $x \to r$ it follows from
$\nabla^2\phi = \frac{1}{r}\frac{d^2(r\phi)}{dr^2} = -\frac{\rho(r)}{\varepsilon_0}$ (SI units).

---

## 1. Centered difference — Eq. (4)

Taylor-expand $y$ about $x_n$, forward ($+$) and backward ($-$):

$$
\begin{aligned}
y_{n\pm1} = y(x_n \pm h) = {} & y_n \pm h\,y'_n + \frac{h^2}{2!}\,y''_n \pm \frac{h^3}{3!}\,y'''_n \\
& + \frac{h^4}{4!}\,y^{(4)}_n \pm \frac{h^5}{5!}\,y^{(5)}_n + \frac{h^6}{6!}\,y^{(6)}_n + O(h^7)
\end{aligned} \tag{a}
$$

Add the two expansions. All odd powers of $h$ cancel (symmetric stencil), so the first neglected term is $O(h^8)$:

$$
y_{n+1} + y_{n-1} = 2y_n + 2\left(\frac{h^2}{2!}\,y''_n + \frac{h^4}{4!}\,y^{(4)}_n + \frac{h^6}{6!}\,y^{(6)}_n\right) + O(h^8) \tag{b}
$$

Move $2y_n$ to the left and absorb the $h^6$ term into $O(h^6)$. This is the paper's centered-difference equation (derivatives at $x_n$):

$$
y_{n+1} - 2y_n + y_{n-1} \cong 2\left(\frac{h^2}{2}\frac{d^2y}{dx^2} + \frac{h^4}{4!}\frac{d^4y}{dx^4} + O(h^6)\right) \tag{4}
$$

Since the $h^3$ and $h^5$ terms cancel, truncating after the $h^4$ term costs only $O(h^6)$.

---

## 2. Insert the ODE — Eqs. (5)–(6)

Denote the linear-in-$y$ side of (1) by

$$
F = U(x) + V(x)\,y \tag{5}
$$

so that on the grid $F_n = U_n + V_n y_n$. Eq. (1) reads $d^2y/dx^2 = F$ for all $x$; evaluating it at $x_n$, and differentiating it twice, gives the two derivatives appearing in (4):

$$
\left.\frac{d^2y}{dx^2}\right|_n = F_n, \qquad
\left.\frac{d^4y}{dx^4}\right|_n = \left.\frac{d^2}{dx^2}\!\left(\frac{d^2y}{dx^2}\right)\right|_n = \left.\frac{d^2F}{dx^2}\right|_n \tag{c}
$$

Insert (c) into (4), with the coefficients

$$
2\cdot\frac{h^2}{2} = h^2, \qquad 2\cdot\frac{h^4}{4!} = \frac{2h^4}{24} = \frac{h^4}{12}, \qquad 2\,O(h^6) = O(h^6),
$$

and solve for $y_{n+1}$:

$$
y_{n+1} = 2y_n - y_{n-1} + h^2F_n + \frac{h^4}{12}\left.\frac{d^2F}{dx^2}\right|_n + O(h^6) \tag{6}
$$

Eq. (6) is not yet a usable scheme: computing $d^2F/dx^2$ analytically (with $y'' = F$),

$$
\frac{d^2F}{dx^2} = U'' + V''\,y + 2V'\,y' + V\,F ,
$$

requires derivatives of $U$ and $V$ and the slope $y'$, which the scheme does not carry.

---

## 3. Second difference for $d^2F/dx^2$

$F = d^2y/dx^2$ is itself smooth, so (b)–(4) hold with $y \to F$:

$$
F_{n+1} - 2F_n + F_{n-1} = h^2\left.\frac{d^2F}{dx^2}\right|_n + \frac{h^4}{12}\left.\frac{d^4F}{dx^4}\right|_n + O(h^6) \tag{d}
$$

Divide by $h^2$:

$$
\left.\frac{d^2F}{dx^2}\right|_n = \frac{F_{n+1} - 2F_n + F_{n-1}}{h^2} + O(h^2) \tag{e}
$$

Multiply by $h^4/12$. The $O(h^2)$ error becomes $O(h^6)$, i.e. the same order as the remainder already present in (6), so no accuracy is lost:

$$
\frac{h^4}{12}\left.\frac{d^2F}{dx^2}\right|_n = \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \tag{f}
$$

---

## 4. Substitute (f) into (6)

$$
\begin{aligned}
y_{n+1} &= 2y_n - y_{n-1} + h^2F_n + \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) + O(h^6) \\
        &= 2y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + (12 - 2)\,F_n + F_{n-1}\big) + O(h^6) \\
        &= 2y_n - y_{n-1} + \frac{h^2}{12}\big(F_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} \tag{g}
$$

Eq. (g) is implicit: $F_{n+1} = U_{n+1} + V_{n+1}y_{n+1}$ contains the unknown $y_{n+1}$.

---

## 5. Use linearity in $y$ — Eq. (7)

Split $F_{n+1}$ according to (5):

$$
\begin{aligned}
y_{n+1} = {} & 2y_n - y_{n-1} + \frac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big) \\
& + \frac{V_{n+1}h^2}{12}\,y_{n+1} + O(h^6)
\end{aligned} \tag{h}
$$

Collect the $y_{n+1}$ terms on the left:

$$
\begin{aligned}
\left(1 - \frac{V_{n+1}h^2}{12}\right) y_{n+1} = {} & 2y_n - y_{n-1} \\
& + \frac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big) + O(h^6)
\end{aligned} \tag{i}
$$

Divide by $1 - V_{n+1}h^2/12 = 1 + O(h^2)$, which is non-zero for $V_{n+1}h^2 \neq 12$ and leaves the remainder $O(h^6)$:

$$
y_{n+1} = \frac{2y_n - y_{n-1} + \dfrac{h^2}{12}\big(U_{n+1} + 10F_n + F_{n-1}\big)}{\left(1 - \dfrac{V_{n+1}h^2}{12}\right)} + O(h^6) \tag{7}
$$

- The right-hand side needs only $y_n$, $y_{n-1}$ (hence $F_n$, $F_{n-1}$, already computed in the previous steps) and $U_{n+1}$, $V_{n+1}$: one evaluation of $U$ and $V$ per step, for a local error $O(h^6)$.
- Linearity in $y$ is used only in this step: for a general $F(x, y)$, Eq. (g) is a nonlinear equation for $y_{n+1}$ (root finding at every step).
- Eq. (7) is a two-step recursion: it needs $y_0$ and $y_1$ to start, which is the subject of the paper's Eqs. (8)–(15).

---

## Remark: leading local error (not in the paper)

Keep the $h^6$ terms of (b) and (d), and use $d^4F/dx^4 = y^{(6)}$:

$$
y_{n+1} - 2y_n + y_{n-1} = h^2F_n + \frac{h^4}{12}\left.\frac{d^2F}{dx^2}\right|_n + \frac{h^6}{360}\,y^{(6)}_n + O(h^8)
$$

$$
\frac{h^4}{12}\left.\frac{d^2F}{dx^2}\right|_n = \frac{h^2}{12}\big(F_{n+1} - 2F_n + F_{n-1}\big) - \frac{h^6}{144}\,y^{(6)}_n + O(h^8)
$$

Hence the $O(h^6)$ in (g) and (7) is

$$
\left(\frac{1}{360} - \frac{1}{144}\right) h^6\,y^{(6)}_n + O(h^8) = -\frac{h^6}{240}\,y^{(6)}_n + O(h^8)
$$

