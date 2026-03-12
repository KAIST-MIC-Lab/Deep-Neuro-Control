





---
## 1 Introduction
- About contraction theory
- Conventional approach of contraction theory based control
- Limitations of the existing methods
- 

## 2 Review of Conventional Contraction based Control
We have contraction constraint which is BMI, as follows:
$$
\ddtt\mm{M}+\mathbb{A}^\top\mm{M}+\mm{M}\mathbb{A}
-2\mm{M}\mm{B}\mm{R}^{-1}\mm{B}^\top\mm{M}
\preceq
-2\alpha\mm{M}
.
$$
---
Additionally, we have input saturation constraint to prevent the norm of the control input from exceeding its limit $\overline{u}$ as follows:
$$
\norm{\mv{u}-\mv{u}_d}^2 \le \overline{u}^2
.
$$
---
Conventionally, we minimize the condition number of the contraction metric $\mm{M}$ according to the ### theorem.


---
## 3 Proposed Method
By pre-and post-multiplying $\mv{e}$, we obtain
$$
-\mv{y}^\top (2\mm{B}\mm{R}^{-1}\mm{B}^\top) \mv{y}
+
e^\top (\tfrac{1}{T_s} \mm{I}+2\mathbb{A}^\top+2\alpha\mm{I})
\mv{y}
-
\left(
\tfrac{\mv{e}^\top\mm{M}_\rm{pre}\mv{e}}{T_s}
\right)
\le
0
$$

and     
$$
\mv{y}(\mm{B}\mm{R}^{-1}\mm{R}^{-1}\mm{B}^\top)\mv{y}
+
(-2\mv{u}_d^\top\mm{R}^{-1}\mm{B}^\top)\mv{y}
+
(\mv{u}_d^\top\mv{u}_d-\overline{u}^2)
\le
0
.
$$
Each inequalities represents the contraction condition and input saturation constraint, respectively. They can be illustrated as ellipses whose feasible domain is outside and inside, respectively. 

In conclusion, following the conventional optimization formulation, we can formulate the following optimization problem:
$$
\begin{aligned}
\min_{\mv{y}}
\mv{y}
\end{aligned}
$$



---
From the first equation, we have the following two scenarios in two-dimensional case.
![[scene_one.jpg]]
In the first scenario, 





Therefore, the optimization problem should be formulated to minimize the steady-state tracking error while satisfying the two constraints as follows:
$$
\begin{aligned}
\min_{\mm{M},\mv{y}} \norm{\mm{M}-\mm{I}}^2
\\
\text{subject to}
\begin{cases}
&
y=\mm{M}\mv{e}
\\
&
\mv{y}^\top (-2\mm{B}\mm{R}^{-1}\mm{B}^\top) \mv{y}
+
e^\top (\tfrac{1}{\der t} \mm{I}+2\mathbb{A}+2\alpha\mm{I})
\mv{y}
+
\left(
-
\tfrac{\mv{e}^\top\mm{M}_\rm{pre}\mv{e}}{\der t}
\right)
\preceq
0
\\
&
\mv{y}(\mm{B}\mm{R}^{-1}\mm{R}^{-1}\mm{B}^\top)\mv{y}
+
(-2\mv{u}_d^\top\mm{R}^{-1}\mm{B}^\top)\mv{y}
+
(\mv{u}_d^\top\mv{u}_d-\overline{u}^2)
\preceq
0
.
\end{cases}
\end{aligned}
$$
The objective function $\norm{\mm{M}-\mm{I}}$ is alternatively defined to minimize the condition number of the contraction metric $\mm{M}$ indirectly. ~~Note that the optimization problem is quadratic constrained quadratic problem(QCQP).~~
