# Numerov's method — step-by-step derivation of Eqs. (1)–(7)

Source: J. L. M. Quiroz González and D. Thompson, *Getting started with Numerov's method*,
Computers in Physics **11**(5), 514–515 (1997).

Notation is the paper's. Equation numbers (1)–(7) are the paper's; intermediate steps are labelled (a)–(i).

- uniform grid of step $h$: $x_{n\pm1} = x_n \pm h$;
- values on the grid: $y_n = y(x_n)$, $U_n = U(x_n)$, $V_n = V(x_n)$, $F_n = F(x_n)$;
- derivatives at $x_n$: $y'_n$, $y''_n$, $y'''_n$, $y^{(k)}_n$ stand for $d^ky/dx^k$ evaluated at $x = x_n$
  (the paper writes the same thing as $\left.\frac{d^2F}{dx^2}\right|_n$ in Eq. (6), and with primes, $y'_0$, $F''_0$, from Eq. (8) on).

---

## 0. Class of equations — Eqs. (1)–(3)

$$
\frac{d^2y}{dx^2} = U(x) + V(x)\,y \tag{1}
$$

Numerov's method relies on two properties of (1):

1. the first derivative $dy/dx$ is absent;
2. $U(x) + V(x)y$ is linear in $y$.

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

| equation | $y$ | $U(x)$ | $V(x)$ of Eq. (1) |
|---|---|---|---|
| (2) Schrödinger | $\psi$ | $0$ | $\frac{2m}{\hbar^2}[V(x) - E]$ |
| (3) forced oscillator | $y$ | $\frac{f_0}{m}\cos\omega x$ | $-\frac{k}{m}$ |
| Poisson, spherical symmetry | $\varphi \equiv r\phi$ | $-\frac{r\rho(r)}{\varepsilon_0}$ | $0$ |

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
