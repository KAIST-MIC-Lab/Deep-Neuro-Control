- Author:
	- Myeongseok Ryu
	- Kyunghwan Choi
	- Sesun You

> This research includes:
> - Contraction theory-based control using LMI approach
> 	- disturbances are considered
> 	- effects of saturation is considered, as well
> - Input saturation handling
> 	- The handling process is integrated in LMI process
> 	- Sector condition and s-procedure are used
> - Numerical validation is demonstrated

## 1 Preliminaries




> [!theorem]
> Consider a nonlinear system of the form $\ddtt\mv{x} = \mv{f}(\mv{x},\mv{x},t)$ and the virtual system $\ddtt\mv{q} = \mv{f}(\mv{q},\mv{x},t)$ are contracting with respect to $\mv{q}$. If a particular solution of virtual system verifies a smooth specific property, then all solutions of the original system verifies this property exponentially.


## 2 System Declaration
Let say we have the following system,
$$
\ddtt{\mv x} = \mv f(\mv x)+\mm B(\mv x)\mysat(\mv u)+\mv{d}
,
$$

^f0d1d1

where .... (omitted here).

The control objective is to design a tracking controller to track some desired trajectory $\mv{x}_d$. In this study, we are assuming that desired tacking controller $\mv{u}_d$ is already given. In other words, the following desired system is given as follows:
$$
\ddtt{\mv x}_d 
= 
\mv f(\mv x_d)
+
\mm B(\mv x_d)
\mv u_d
,
$$

^32b756

where ... (omitted here). 
Note that the desired controller is designed not to exceed the input limits with respect to the zero-disturbances, *i.e.,* $\mv{d}=\mv{0}$.
## 3 Controller Design

### 3.1 LMI Formulation of Contraction Theory
Using the perturbed system [[#^f0d1d1]] and the desired system [[#^32b756]], the error dynamics can be obtained as follows:
$$
\ddtt\mv{x}-\ddtt\mv{x}_d
=
\mv{f}-\mv{f}_d
+
\mv{g}\mysat(\mv{u})-\mv{g}_d\mv{u}_d
+
\mv{d}
.
$$

^20228c

The state-dependent coefficient (SDC) formulation [] is used to rewrite the error dynamics [[#^20228c]] linearly as follows:
$$
\begin{aligned}
\ddtt
(
\mv{x} - \mv{x}_d
)
=&
\mv{f} - \mv{f}_d
+ \mm{B} \mv{u}
- \mm{B}_d \mv{u}_d
+ \mv{d}
\\
=&
\underbrace
{
\mv{f} - \mv{f}_d
+
(
\mm{B}
- \mm{B}_d
) \mv{u}_d
}_{= \mathbb{A}(\mv{x},\mv{x}^*) (\mv{x} - \mv{x}_d)}
- \mm{B}
(
\mysat(\mv{u}) - \mv{u}_d
)
+ \mv{d}
\\
=&
\mathbb{A} (\mv{x} - \mv{x}_d)
- \mm{B} (\mv{u} - \mv{u}_d)
+ \mm{B}
\underbrace{
(\mysat(\mv{u}) - \mv{u})
}_{=:\\mv{phi}(\mv{u})}
+ \mv{d}
,
\end{aligned}
$$

^cdaa47

The control law $\mv{u}$ is defined with contraction metric $mm{M}$ as follows:
$$
\mv{u} = \mv{u}_d - \mm{R}^{-1}\mm{B}^\top\mm{M}(\mv{x}-\mv{x}_d)
,
$$
Then, [[#^cdaa47]] is represented as follows:
$$
\begin{aligned}
\ddtt
(
\mv{x} - \mv{x}_d
)
=&
\left(
\mathbb{A} - \mm{B} \mm{R}^{-1} \mm{B} ^\top \mm{M}
\right)
(\mv{x} - \mv{x}_d)
+ \mm{B}\mv{\phi}(\mv{u})
+
\mv{d}
,
\end{aligned}
$$

According to [], the partial contraction theory is used. First, the virtual system is defined as follows:
$$
\ddtt\mv{q} 
= 
\ddtt\mv{x}_d
+
\left(
\mathbb{A} - \mm{B} \mm{R}^{-1} \mm{B} ^\top \mm{M}
\right)
(\mv{q}-\mv{x}_d)
+ \mm{B}\mv{\phi}(\mv{u})
+
\mv{d}_\mu
,
$$

^9d8748

where $\mu\in[0,1]$ and $\mv{d}_\mu:=\mu\mv{d}$. The virtual system is parameterized by $\mu$, which has $\mv{q}(\mu=0,t)=\mv\xi_0$ and $\mv{q}(\mu=1,t)=\mv\xi_1$.

Note the the feasible solutions of the virtual $\mv{q}$  system in [[#^9d8748]] are $\mv{q}=\mv{x}$ and $\mv{q}=\mv{x}_d$. That is, the two solutions will converged, if the virtual system is contracting.

Assume that the contraction metric $\mm{M}$ satisfies the contraction condition presented in []. Then the following theorem can be used to ...
> [!Theorem]
> If the system [[#^9d8748]] is contracting with contraction metric $\mm{M}$ and decay rate $\alpha$, *i.e.,* $\mm{M}\mathbb{A}+\mathbb{A}\mm{M}-\mm{M}\mm{B}\mm{R}^{-1}\mm{B}^\top\mm{M}\le-2\alpha\mm{M}$, and $\Theta\mv{d}\in\mathcal{L}_\infty$ and $\exists\underline{m},\overline{m}\in\R_{>0}$, and $\exists\overline{d}\in\R_{\ge0}$, subject to $\overline{d}=\sup_{\mv{x},t}\norm{\mv{d}}$, and $\underline{m}\mm{I}\le \mm{M}\le\overline{m}\mm{I}$.
> Then, we have the following relation: 
> $$
> \norm{\xi_1(t)-\xi_0(t)}
> \le
> 
> $$
> .....

`\begin{proof}`
According the [], we have $\mm{M}=\mm{\Theta}^\top\mm{\Theta}$. The 

`\end{proof}`


### 3.2 Saturation Handling
In this study, we use *sector condition* to formulate the constraint into inequality condition. First, let us consider following saturation constraint:
$$

$$
Then, the following constraint can be obtained:
$$

$$

## 4 Numerical Validation



## 5 Conclusion