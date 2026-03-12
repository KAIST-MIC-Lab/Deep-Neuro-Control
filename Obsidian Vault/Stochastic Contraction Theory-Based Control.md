### 0.1 Contributions
- Energy-based effective space approach in contraction based-control design is extended for the stochastic systems.
- Stochastic constraints in controller, such as input noise, are handled.

## 1 Conventional Stochastic Contraction Theory-Based Control
### 1.1 System and Control Input Declaration 
Consider we have the following stochastic system as follows:
$$
\der\mv{x}
=
\left(
\mv{f}(\mv{x},t)+\mm{B}(\mv{x},t)\mv{u]}
\right)
\der t 
+
\mm{G}_c(\mv{x},t) \der\mathscr{W}(t)
,
$$
where $\mv{f},\mm{B}$ are system functions and $\mm{G_c}$ denotes bounded external disturbance matrix, $\mathscr{W}$ denotes $w$-dimensional Wiener process.
We use the following feedback control input 
$$
\mv{u}
=
\mv{u}_d
-
\mm{R}^{-1} \mm{B}^\top \mm{M}\mv{e}
,
$$
where $\mm{R}$ denotes weighting matrix and $\mm{M}$ denotes contraction metric which satisfies $\underline{m}\mm{I}\preceq\mm{M}\preceq\overline{m}\mm{I}$. 

---
### 1.2 Stochastic Contraction Theory
Assume that the following inequalities hold for uniqueness of the solution:
$$
\exists L\in\R_{\ge0},
\quad \text{s.t.}
\norm{\mv{f}(\mv{x},t) - \mv{f}(\mv{x}',t)}
+
\norm{
\mm{B}(\mv{x},t)\mv{u}
-
\mm{B}(\mv{x}',t)\mv{u}
}
+
\norm{
\mm{G}_c(\mv{x},t)-\mm{G}_c(\mv{x}',t)
}_F
\le
L\norm{\mv{x}-\mv{x}'}
$$
and
$$
\exists \overline{L}\in\R_{\ge0},
\quad \text{s.t.}
\norm{\mv{f}(\mv{x},t)}^2
+
\norm{
\mm{B}(\mv{x},t)\mv{u}(\mv{x})
}^2
+
\norm{
\mm{G}_c(\mv{x},t)
}_F^2
\le
\overline{L}\norm{1+\norm{\mv{x}}^2}
$$
In addition, assume that the contraction metric satisfies the following Lipsitz condition:
$$
\norm{
\pptfrac{{\mm{M}(\mv{x},t)}}{x_i} - \pptfrac{{\mm{M}(\mv{x}',t)}}{x'_i}
}
\le L_m \norm{\mv{x}-\mv{x}'}
$$


Let say what we have two solutions of the given stochastic system as $\mv{\xi}_0$ and  $\mv{\xi}_1$. 


$$
\alpha_s = 
L_m (\overline{g}_0^2+\overline{g}_1^2)(\alpha_G+\tfrac{1}{2})
$$

$$
2ab \le \alpha_G^{-1} a^2 + \alpha_G b^2
$$
$$ 
a = \sqrt{2\overline{m}}, b =\sqrt{L_m}\norm{\partial_\mu x}
$$

Path-length integral
$$
V_{s\mathscr{l}} = \int_{\mv{x}_d}^{\mv{x}} \delta\mv{q}^\top\mm{M}\delta\mv{q}
$$

$$
\ddtt\mm{M}
+
2\mysym
\left(
\mm{M} \pptfrac{\mv{f}}{\mv{x}}
\right)
\preceq
-2\alpha\mm{M}
-
\beta\mm{I}
,
$$

Then the trajectory error is 
$$
\mathbb{E}
\left[\norm{\mv{x}(t)-\mv{x}_d(t)}^2 \right]
\le
\tfrac{\mathbb{E}(V_{xl}(0))}{\underline{m}}\exp(-2\alpha t)
+ \tfrac{C_C}{2\alpha} \chi
$$

$$ 
C_C = \overline{g}_c^2 ( 2\alpha_G^{-1}+1)
$$


$$
\mathbb{P}
\left[
\norm{\mv{x}(t)-\mv{x}_d(t)}
\ge \epsilon
\right]
\le
\tfrac{1}{\epsilon^2}
\left(
\tfrac{\mathbb{E}[V_{s\ell}(0)]}{\underline{m}}
\exp(-2\alpha t)
+
\tfrac{C_C}{2\alpha}\chi
\right)
$$
