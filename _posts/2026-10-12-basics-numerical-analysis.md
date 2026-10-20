---
title: "Basics of Numerical Analysis"
#date:
#excerpt_separator: "<!--more-->"
last_modified_at: 2026-10-12T16:20:02-05:00
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

# Basics of numerical analysis

We summarise the important equations derived in this article, for a function  ${\displaystyle f:\mathbb{R} \to \mathbb{R} }$ that is ${N \geq 3}$ times differentiable:


$$
f'(x) = \frac{\delta_h[f](x)}{h} + o(h^2) \; , \quad \quad \delta_h[f](x) := f(x+\tfrac{h}{2}) - f(x-\tfrac{h}{2})  \tag{2.4a}
$$

$$
f''(x) = \frac{\delta^2_h[f](x)}{h^2} + o(h^2) \; , \quad \quad \delta^2_h[f](x) = f(x+h) - 2\,f(x) + f(x-h) \tag{3.1a} 
$$

Rewritten using the index notation $x_n$ and $f_n = f(x_n)$ and assuming that all the elements are equally spaced by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N$, for the first derivative:


$$
f_{n}' = \begin{cases}\tag{2.5a-c}
\displaystyle \frac{f_{n+1} - f_{n-1}}{2h} + O(h^2) \\[10pt]
\displaystyle \frac{f_{n+1} - f_n}{h}  + o(h)  \\[10pt]
\displaystyle \frac{f_{n} - f_{n-1}}{h}  + o(h)  \\[10pt]
\end{cases}
$$

and for the second:

$$
f_{n}'' = \begin{cases}\tag{3.2a-c}
\displaystyle \frac{f_{n+1} - 2f_n + f_{n-1}}{h} + o(h^2) \\[10pt]
\displaystyle \frac{f_{n+2} - 2f_{n+1} + f_{n} }{h} + o(h)  \\[10pt]
\displaystyle \frac{f_{n} - 2f_{n-1} + f_{n-2} }{h} + o(h)  \\[10pt]
\end{cases}
$$




## Taylor's theorem

**Taylor's theorem** Given a function ${\displaystyle f:\mathbb{R} \to \mathbb{R} }$ that is ${N \geq 1}$ times differentiable at the point ${a\in \mathbb {R} }$, there exists a function ${\displaystyle g_{N}:\mathbb {R} \to \mathbb {R} }$ such that

$$ 
f(x) = \sum_{k=0}^{N}{\frac {f^{(k)}(a)}{k!}}(x-a)^{k} + g_{N}(x)(x-a)^{k} \; , \quad \quad \lim _{x\to a}g_{N}(x)=0.
\tag{1.1}
$$

The polynomial appearing in Taylor's theorem is the ${\displaystyle {\boldsymbol {N}}}$-th order Taylor polynomial of the function ${\displaystyle f}$ at the point ${\displaystyle a}$

$$
P_{N}(x) = f(a)+f'(a)(x-a)+{\frac {f''(a)}{2!}}(x-a)^{2}+\cdots +{\frac {f^{(N)}(a)}{N!}}(x-a)^{N}
\tag{1.2}
$$

Using the little-$o$ and big-$O$ notation, the statement in Taylor's theorem reads as

$$
f(x) = \sum_{k=0}^{N}{\frac {f^{(k)}(a)}{k!}}(x-a)^{k} + o(\lvert x-a \rvert^N)
\tag{1.3a}
$$

$$
f(x) = \sum_{k=0}^{N}{\frac {f^{(k)}(a)}{k!}}(x-a)^{k} + O(\lvert x-a \rvert^{N+1})
\tag{1.3b}
$$

We can rewrite the theorem in a slightly different notation, with the substitutions $x \rightarrow x+h$ and $a \rightarrow x$

$$
\begin{align*} \tag{1.4a-b}
f(x+h)  &= \sum_{k=0}^{N}{\frac {h^k}{k!}}f^{(k)}(x) + o(h^N) & 
f(x+h)  &= \sum_{k=0}^{N}{\frac {h^k}{k!}}f^{(k)}(x) + O(h^{N+1})\\[15pt]
        &= f(x) + o(1) &
        &= f(x) + O(h) \\[7pt]
        &= f(x) + h \, f'(x) + o(h) &
        &= f(x) + h \, f'(x) + O(h^2) \\[7pt]
        &= f(x) + h \, f'(x) + \frac{h^2}{2} f''(x) + o(h^2) & 
        &= f(x) + h \, f'(x) + \frac{h^2}{2} f''(x) + O(h^3) \\[7pt]
        &= ... &
        &= ...
\end{align*}
$$

Rewritten using the index notation $x_n$ and $f_n = f(x_n)$ and assuming that all the elements are equally spaced by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N$

$$
\begin{align*} \tag{1.5}
f_{n+1}  &= \sum_{k=0}^{N}{\frac {h^k}{k!}}f^{(k)}_n + o(h^N) \\[15pt]
        &= \sum_{k=0}^{N}{\frac {h^k}{k!}}f^{(k)}_n + O(h^{N+1}) \\[15pt]
        &= f_n + h \, f_n' + \frac{h^2}{2} f_n'' + \frac{h^3}{6} f_n''' + ... \\[7pt]
\end{align*}
$$



## First order finite differences

The derivative of a function $f$ at a point $x$ is defined by the limit

$$
f'(x) := \lim_{h \rightarrow 0^{+}} \frac{f(x+h)-f(x)}{h}
\tag{2.1}
$$

We define the forward difference ${\Delta_{h}[f]}$, backward difference ${\nabla_{h}[f]}$ and central difference ${\delta_{h}[f]}$ as:

$$
\begin{align}
    \mathrm{Forward  \; difference:} \quad &\Delta_{h}[f](x) := f(x+h)-f(x) \tag{2.2a}\\[10pt]
    \mathrm{Backward \; difference:} \quad &\nabla_{h}[f](x) := f(x)-f(x-h) \tag{2.2b}\\[10pt]
    \mathrm{Central  \; difference:} \quad &\delta_{h}[f](x) := f\left(x+\frac{h}{2}\right)-f\left(x-\frac{h}{2}\right) \tag{2.2c}\\[15pt]
\end{align}
$$
$$
\begin{align}
    \Rightarrow \quad \nabla_{h}[f](x) &= \Delta_{h}[f](x-h) \tag{2.3a}\\[15pt]
    \Rightarrow \quad \delta_{h}[f](x) &= \Delta_{\tfrac {h}{2}}[f](x) + \nabla_{\tfrac {h}{2}}[f](x) \tag{2.3b}\\[15pt]
\end{align}
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

To summarize:

$$
\begin{align}
f'(x) &= \frac{\delta_h[f](x)}{h} + o(h^2) \; , \quad \quad \delta_h[f](x) := f(x+\tfrac{h}{2}) - f(x-\tfrac{h}{2}) \tag{2.4a} \\[10pt]
f'(x) &= \frac{\Delta_h[f](x)}{h} + o(h) \; , \quad \quad \Delta_h[f](x) := f(x+h) - f(x) \tag{2.4b}\\[10pt]
f'(x) &= \frac{\nabla_h[f](x)}{h} + o(h) \; , \quad \quad \nabla_h[f](x) := f(x) - f(x-h) \tag{2.4c}\\[10pt]
\end{align}
$$

Rewritten using the index notation $x_n$ and $f_n = f(x_n)$ and assuming that all the elements are equally spaced by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N$

$$
\begin{align}
f_{n}' &= \frac{f_{n+1} - f_n}{h} + o(h) \tag{2.5b}\\[10pt]
f_{n}' &= \frac{f_{n} - f_{n-1}}{h} + o(h) \tag{2.5c} \\[10pt]
\end{align}
$$

Unfortunately, the central difference has (and will for every odd $n$) $h$ multiplied by the non-integer $1/2$. This is often a problem because it amounts to changing the interval of discretization, but we can remedy by applying the Taylor theorem of Eq.(1.4):

$$
(1.4) \; : \begin{cases}
f(x + h) = f(x) + h \, f'(x) + \frac{h^2}{2} \, f''(x) + \frac{h^3}{6} \, f'''(x) + ... \\[7pt]
f(x - h) = f(x) - h \, f'(x) + \frac{h^2}{2} \, f''(x) - \frac{h^3}{6} \, f'''(x) + ...\\[10pt]
\end{cases}
$$

which for the sum leads to:

$$
\, \Rightarrow \; 
\frac{f(x + h) + f(x-h)}{2} =  f(x) + \frac{h^2}{2} \, f''(x) + ... 
$$

$$
\tag{2.6}
\Rightarrow \quad f(x) = \frac{f(x+h) + f(x-h)}{h} + O(h^2)
$$


and for the difference leads to:

$$
\, \Rightarrow \; 
\frac{f(x + h) - f(x-h)}{2} =  h\,f'(x) + \frac{h^3}{6} \, f'''(x) + ... 
$$

$$
\tag{2.7}
\Rightarrow \quad f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
$$

$$
\tag{2.5a}
\Rightarrow \quad f_{n}' = \frac{f_{n+1} - f_{n-1}}{2h} + O(h^2)
$$


## Second order finite differences

In an analogous way, one can obtain finite difference approximations to higher order derivatives and differential operators. 
For example, by using the above central difference formula for $f'(x + ⁠h/2⁠)$ and $f'(x − ⁠h/2⁠)$ and applying a central difference formula for the derivative of $f'$ at ${x}$, we obtain the central difference approximation of ${f''}$:

$$
f'\left(x+\frac{h}{2}\right) = \frac{f(x+h) - f(x)}{h} + o(h^2) \quad , \quad \quad
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

Analogously for backward and forward differences:

$$
\begin{align*}
f''(x) &= 
    \frac{1}{h}\left[f'(x+h) - f'(x)\right] + o(h) \\[15pt]
    &=\frac{1}{h}\left[
        \frac{f(x+2h) - f(x+h)}{h} - \frac{f(x+h) - f(x)}{h}
    \right] + o(h) \\[15pt]
    &= \frac{f(x+2h) - 2\,f(x+h) + f(x)}{h^2} 
     + o(h) \\[15pt]
f''(x) &= 
    \frac{1}{h}\left[f'(x) - f'(x-h)\right] + o(h) \\[15pt]
    &=\frac{1}{h}\left[
        \frac{f(x) - f(x-h)}{h} - \frac{f(x-h) - f(x-2h)}{h}
    \right] + o(h) \\[15pt]
    &= \frac{f(x) - 2\,f(x-h) + f(x-2h)}{h^2} 
     + o(h) \\[15pt]
\end{align*}
$$

Note that we can also rewrite it in the following way:

$$
\begin{align*}
\Delta^2_h[f](x) &:=  \Delta_h[\Delta_h[f]](x)  \\[10pt]
    &= \Delta_h[g_h](x) \quad , \quad \quad g_h(x):=f(x+h)-f(x) \\[10pt]
    &= g_h(x+h)-g_h(x) \\[10pt]
    &=\left[ f(x+2h) - f(x+h)\right] - \left[ f(x+h) - f(x)\right] \\[10pt]
    &=f(x+2h) - 2\,f(x+h) + f(x) \\[10pt]
\Rightarrow \quad f''(x) &= \frac{\Delta^2_h[f](x)}{h^2} + o(h) \\[15pt]

\nabla^2_h[f](x) &:=  \nabla_h[\nabla_h[f]](x)  \\[10pt]
    &= \nabla_h[g_h](x) \quad , \quad \quad g_h(x):=f(x)-f(x-h) \\[10pt]
    &= g_h(x)-g_h(x-h) \\[10pt]
    &=\left[ f(x) - f(x-h)\right] - \left[ f(x-h) - f(x-2h)\right] \\[10pt]
    &=f(x) - 2\,f(x-h) + f(x-2h) \\[10pt]
\Rightarrow \quad f''(x) &= \frac{\nabla^2_h[f](x)}{h^2} + o(h) \\[15pt]

\delta^2_h[f](x) &:=  \delta_h[\delta_h[f]](x)  \\[10pt]
    &= \delta_h[g_h](x) \quad , \quad \quad g_h(x):=f(x+\tfrac{h}{2})-f(x-\tfrac{h}{2}) \\[10pt]
    &= g_h(x+\tfrac{h}{2})-g_h(x-\tfrac{h}{2}) \\[10pt]
    &=\left[ f(x+h) - f(x)\right] - \left[ f(x) - f(x-h)\right] \\[10pt]
    &=f(x+h) - 2\,f(x) + f(x-h) \\[10pt]
\Rightarrow \quad f''(x) &= \frac{\delta^2_h[f](x)}{h^2} + o(h^2) \\[15pt]
\end{align*}
$$

So, to summarize:

$$
\begin{align}
f''(x) &= \frac{\delta^2_h[f](x)}{h^2} + o(h^2) \; , \quad \quad \delta^2_h[f](x) = f(x+h) - 2\,f(x) + f(x-h) \tag{3.1a}\\[8pt]
f''(x) &= \frac{\Delta^2_h[f](x)}{h^2} + o(h) \; , \quad \quad \Delta^2_h[f](x) = f(x+2h) - 2\,f(x+h) + f(x)    \tag{3.1b}\\[8pt]
f''(x) &= \frac{\nabla^2_h[f](x)}{h^2} + o(h) \; , \quad \quad \nabla^2_h[f](x) = f(x) - 2\,f(x-h) + f(x-2h)    \tag{3.1c}\\[8pt] 
\end{align}
$$


Rewritten using the index notation $x_n$ and $f_n = f(x_n)$ and assuming that all the elements are equally spaced by $h=x_{n+1}-x_{n}\, ,\; \forall\, n=0,...,N$

$$
\begin{align}
f_{n}'' &= \frac{f_{n+1} - 2f_n + f_{n-1}}{h} + o(h^2) \tag{3.2a}\\[10pt]
f_{n}'' &= \frac{f_{n+2} - 2f_{n+1} + f_{n} }{h} + o(h) \tag{3.2b} \\[10pt]
f_{n}'' &= \frac{f_{n} - 2f_{n-1} + f_{n-2} }{h} + o(h) \tag{3.2b} \\[10pt]
\end{align}
$$


## Higher order finite differences

More generally, the $n$-th order forward, backward, and central differences are given by the following expressions, which we will not derive here:

$$
\begin{align}
\frac{\mathrm{d}^{n}f}{\mathrm{d}x^n}(x) = \frac{\delta_h^n[f](x)}{h^n}+o(h^2) 
\; , \quad \quad 
\delta_h^n[f](x) &= \sum_{j=0}^n (-1)^{j}  \, \binom{n}{j} \, f(x+(\tfrac{n}{2}-j)h) \tag{4.1a}\\[15pt]
\frac{\mathrm{d}^{n}f}{\mathrm{d}x^n}(x) = \frac{\Delta_h^n[f](x)}{h^n}+o(h) 
\; , \quad \quad 
\Delta_h^n[f](x) &= \sum_{j=0}^n (-1)^{n-j}\, \binom{n}{j} \, f(x+jh) \tag{4.1b}\\[10pt]
\frac{\mathrm{d}^{n}f}{\mathrm{d}x^n}(x) = \frac{\nabla_h^n[f](x)}{h^n}+o(h) 
\; , \quad \quad 
\nabla_h^n[f](x) &= \sum_{j=0}^n (-1)^{j}  \, \binom{n}{j} \, f(x-jh) \tag{4.1c}\\[15pt]
\binom{n}{j} :&= \frac{n!}{j!(n-j)!} \notag
\end{align}
$$

