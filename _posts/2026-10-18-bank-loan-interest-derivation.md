---
title: "Bank interest functions: analytical derivation"
#date:
#excerpt_separator: "<!--more-->"
last_modified_at: 2026-10-18T16:20:02-05:00
categories:
  - Blog
tags:
  - Post Formats
  - math
toc: true   # Table-Of-Content
toc_label: "Table of Contents"
toc_icon: "cog"
# link: https://github.com


---

Notation:

- $D$ : loan amount
- $R$ : monthly installment
- $I$ : nominal annual interest rate
- $i = \frac{I}{12}$ : monthly interest rate, 30/360 convention
- $N$ : total number of installments
- $T = D + Z$ : total paid (loan + interests)
- $Z$ : total interests paid

For the $n$-th installment:

- $d_n$ : debt left after the $n$-th installment ($d_0=D$)
- $p_n = R - i\,d_{n-1}$ : principal part of the $n$-th installment
- $z_n = i\,d_{n-1}$ : interest part of the $n$-th installment

Here a summary of all the equations derived in this article:

$$
\boxed{N = \lceil n^* \rceil}
\quad \quad \quad
\boxed{\;n^* = \frac{\ln{(R)} - \ln{(R-i\,D)}}{\ln{(i+1)}}} 
$$

$$
R_n = \begin{cases}
R & \quad \forall \, n=1,...,N-1 \\[8pt]
R_{\text{last}} \le R & \quad n={N}
\end{cases}
\quad \quad \quad
\boxed{\;R_\mathrm{last} = \frac{R}{i}(i+1) + \left(D - \frac{R}{i}\right) (i+1)^N}
$$

$$
\boxed{
\; R = \frac{i \, (1 + i)^{N}}{(1 + i)^{N} - 1}D
= \frac{i}{1 - (1 + i)^{-N}} D
}
\quad \quad \quad
\boxed{\;
D = \frac{(1 + i)^{N} - 1}{i \, (1 + i)^{N}} R 
= \frac{1 - (1 + i)^{-N}}{i}\, R
}
$$

$$d_n = D\,(1+i)^n - R\,\frac{(1+i)^n - 1}{i}$$
$$p_n = (R - iD)\,(i+1)^{n-1}, \qquad z_n = R - (R - iD)\,(i+1)^{n-1}\,$$

$$
\boxed{
\;T =(N - 1)\,R + R_\mathrm{last} \approx n^* R
}
\quad \quad \quad 
\boxed{
\;Z = (N - 1)\,R + R_\mathrm{last} - D \approx n^* R - D
}
\quad \quad \quad 
P_\mathrm{tot} = D
$$


### Check: amortization condition

The debt decreases iff:

$$
R > i\,D \quad\Longleftrightarrow\quad 12R > I\,D
$$

If $R \le iD$, then $d_n \ge D$ for all $n$ and the loan is never repaid.

### Debt left to be payed after $n$ installments

Each month the bank adds the interest on the current balance, then subtracts the installment:

$$
\tag{1}
\begin{cases}
d_0 &= D \\[10pt]
d_n &= d_{n-1} - (R - i\,d_{n-1}) & \quad\quad \forall n > 0\\[4pt]
\end{cases}
$$

With simple algebra:

$$
\begin{aligned}
\Rightarrow \quad d_n &= d_{n-1} - (R - i\,d_{n-1}) \quad\quad \forall n > 0 \\[6pt]
    & \quad\quad q := 1 + i\\[6pt]
    &= q\,d_{n-1} - R  \\[6pt]
    &= q\,(q\,d_{n-2} - R) -R = q^2 d_{n-2} - R(1 + q) \\[6pt]
    &= q^2\,(q\,d_{n-3} - R) - R(1 + q) = q^3 d_{n-3} - R(1 + q + q^2) \\[6pt]
    &= \dots \\[6pt]
    &= q^n d_0 - R (1 + q + q^2 + ... + q^{n-1}) \\[6pt]
    & \quad\quad d_0 := D\\[6pt]
    &= q^n D - R \sum_{k=0}^{n-1} q^k \\[6pt]
    &= q^n D - R\,\frac{q^n - 1}{q - 1} \\[6pt]
    & \quad\quad q = 1 + i\\[6pt]
    &= D\,(1+i)^n - R\,\frac{(1+i)^n - 1}{i}
\end{aligned}
$$

$$
\boxed{\,d_n = D\,(1+i)^n - R\,\frac{(1+i)^n - 1}{i}}
$$

From this equation, we can easily derive the relation between the total number of installments $N$, the loan $D$, the interest rate $i$ and the installment $R$ by imposing the condition $d_N=0$ (i.e. the debt must be repaid at the end of last installment). 

$$
\begin{align*}
d_N = 0 
\quad&\iff\quad 0 = \frac{R}{i} + \left(D - \frac{R}{i}\right) (i+1)^N \\[10pt]
&\iff\quad \frac{R-i\,D}{i}(i+1)^N = \frac{R}{i} \\[10pt]
&\iff\quad \boxed{(i+1)^N = \frac{R}{R-i\,D}}
\end{align*}
$$

We'll derive all the expressions later.

**Equivalent derivation from fixed-point theorem**

The fixed point of $d \mapsto q d - R$ is $d^* = R/i$, and the deviation from it grows geometrically:

$$
\begin{aligned}
d_{n-1} \mapsto d_{n} &= q d_n - R \\[8pt]
d_{n-1} = d_{n} &\iff d_n=d^*:=\frac{R}{q-1}=\frac{R}{i} \\[8pt]
d_{n} - d^* &= d_{n} - \frac{R}{i} \\[6pt]
        &= q \, d_{n-1} - R - \frac{R}{i}  \\[6pt] 
        &= q \left( d_{n-1} - \frac{1}{\cancel{i+1}}\frac{\cancel{(i+1)}R}{i}\right) \\[8pt]
        &= q \left(d_{n-1} - d^*\right) \\[10pt]
        &= q^2 \left(d_{n-2} - d^*\right) \\[8pt]
        &= ... \\[10pt]
        &= q^n \left(d_{0} - d^*\right) \\[8pt]
\end{aligned}
$$

$$
\;\Longrightarrow\;
\boxed{d_n  - d^* = q^n \left(d_{0} - d^*\right)}
$$


$$
\;\Longrightarrow\;
\boxed{d_n = \frac{R}{i} + \left(D - \frac{R}{i}\right) (i+1)^n}
$$

### Total number of installments

To get the number of installments needed to repay a debt, we can set $d_n = 0$ and get the value $n^*$:

$$
\begin{align*}
d_n = 0 
\quad&\iff\quad 0 = \frac{R}{i} + \left(D - \frac{R}{i}\right) (i+1)^n \\[10pt]
&\iff\quad \frac{R-i\,D}{i}(i+1)^n = \frac{R}{i} \\[10pt]
&\iff\quad (i+1)^n = \frac{R}{R-i\,D} \\[10pt]
&\iff\quad n \, \ln{(i+1)} = \ln{\left(\frac{R}{R-i\,D}\right)} \\[10pt]
&\iff\quad n = n^* :=\frac{\ln{(R)} - \ln{(R-i\,D)}}{\ln{(i+1)}} \\[10pt]
\end{align*}
$$

In the previous expression, $n^{\*}$ is real-valued => the (integer) number of installments is $N = \lceil n^{\*} \rceil$.
However, the last installment, that pays the remaining debt plus one month of interest on it, is smaller than (or equal to) the standard $R$:

$$
\boxed{\,n^* = \frac{\ln{(R)} - \ln{(Ri\,D)}}{\ln{(i+1)}}}
$$
$$
\mathrm{Total \; number \; of \; installments:} \quad \boxed{N = \lceil n^* \rceil}
$$

$$
R_n = \begin{cases}
R & \quad \forall \, n=1,...,N-1 \\[8pt]
R_{\text{last}} \le R & \quad n={N}
\end{cases}
$$

$$
\boxed{
R_\mathrm{last} := (i+1)\,d_{N-1}
= \frac{R}{i}(i+1) + \left(D - \frac{R}{i}\right) (i+1)^N}
$$

### Installment and Loan from the other terms

Recalling the equation we derived by setting $d_N=0$, we can get all the inverses:

$$
\begin{align*}
d_N = 0 
\quad&\iff\quad 0 = \frac{R}{i} + \left(D - \frac{R}{i}\right) (i+1)^N \\[10pt]
\quad&\iff\quad R + (iD - R) (i+1)^N = 0 \\[10pt]
\quad&\iff\quad \left[1-(i+1)^N\right]R + i(i+1)^N D = 0 \\[10pt]
\end{align*}
$$

**Installment from $D$, $I$, $N$**


$$
\boxed{
R = \frac{i \, (1 + i)^{N}}{(1 + i)^{N} - 1}D 
= \frac{i}{1 - (1 + i)^{-N}} D
}
$$

**Loan from $R$, $I$, $N$**

$$
\boxed{
D = \frac{(1 + i)^{N} - 1}{i \, (1 + i)^{N}} R
= \frac{1 - (1 + i)^{-N}}{i}\, R
}
$$

**Interest rate from $D$, $R$, $N$: no closed form**

Substitute $x = 1/q$ into the present-value condition:

$$
g(x) \equiv R\,(x + x^2 + \dots + x^N) - D = 0
$$

This is a polynomial of degree $N$ in $x$. For $N \ge 5$ there is no general solution in radicals (Abel–Ruffini), and $N$ is usually in the hundreds. The same equation defines the APR / internal rate of return.

Existence and uniqueness of a positive rate: $g$ is strictly increasing on $x > 0$, with $g(0) = -D < 0$ and $g(x) \to \infty$. So there is exactly one positive root $x^*$. It satisfies $x^* < 1$, i.e. $i > 0$, iff

$$
g(1) = R\,N - D > 0
$$

Numerically, in terms of $i$: the present value

$$
\mathrm{PV}(i) = R\,\frac{1 - (1 + i)^{-N}}{i}, \qquad
\lim_{i \to 0} \mathrm{PV}(i) = R\,N
$$

is strictly decreasing in $i$. Bisection on $\mathrm{PV}(i) = D$ therefore converges to the unique root, in $\lceil \log_2(\text{width}/\text{tol}) \rceil \approx 43$ steps for width $1$ and tolerance $10^{-13}$.

For small $i$, evaluate $1 - (1 + i)^{-N}$ as `-expm1(-N*log1p(i))` to avoid cancellation.

### Principal part and Interest part

Recalling the definitions of the principal part $p_n$ and interest part $z_n$ of each installment

$$
\begin{align*}
z_n &:= i\,d_{n-1}               \\[10pt]
p_n &:= R - i\,d_{n-1} = d_{n-1} - d_n
\end{align*}
$$

we can get their expressions by inserting the fixed-point form of $d_{n-1}$ into $p_n$:

$$
\begin{align*}
p_n &= R - i\,d_{n-1} \\[6pt]
    &=\cancel{R} - i\left[\cancel{\frac{R}{i}} + \left(D - \frac{R}{i}\right) q^{n-1}\right] \\[8pt]
    &= (R-iD)\,q^{n-1} \\[10pt]
z_n &= i\,d_{n-1} \\[6pt]
    &= R - (R - iD) q^{n-1}
\end{align*}
$$

$$
\boxed{\,p_n = (R - iD)\,(i+1)^{n-1}, \qquad z_n = R - (R - iD)\,(i+1)^{n-1}\,}
$$

The principal part grows geometrically, so the interest part shrinks.

### Total paid and total paid interest

We can get the total paied $T$ and the total interest paid $Z$ easily (the total principal paid $P_\mathrm{tot}$ is trivial):

$$
\boxed{
T =(N - 1)\,R + R_\mathrm{last} \approx n^* R}
$$

$$
\boxed{
Z = (N - 1)\,R + R_\mathrm{last} - D \approx n^* R - D
}
$$

$$
P_\mathrm{tot} = D
$$

Treating $N$ as continuous gives $Z \approx R\,n^* - D$, with an error below $R$.


## Function definitions

### Analytical definitions

Note: for numerical reasons, it's better to rewrite some terms as follows and use the numerical functions $\mathrm{log1p}(x)$ and $\mathrm{expm1}(x)$:

$$
\begin{align*}
\ln(1+i) &= \mathrm{log1p}(i) \\[6pt]
1-(1+i)^{-N} &= - \left[(1+i)^{-N} -1 \right] = - \left[e^{-N\,\ln(1+i)} -1 \right] = -\mathrm{expm1}(-N \, \mathrm{log1p}(i))\\[10pt]
\mathrm{log1p}(x) :&= \ln(1+x) \\[6pt]
\mathrm{expm1}(x) :&= e^x -1
\end{align*}
$$
