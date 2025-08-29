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
Consider a nonlinear system as follows:
$$
\ddtt\mv{x}
=
\mv{f}(\mv{x},t)
.
$$

^ffa982

The differential dynamics can be obtained using partial derivative as follows:
$$
\ddtt\delta\mv{x} = \pptfrac{\mv{f}}{\mv{x}}\delta\mv{x}
.
$$

The system [[#^ffa982]] is contracting according to the following theorem.
> [!theorem] Contraction Condition
> The system [[#^ffa982]] is contracting, if there exists the contraction metric $\mm{M}:=\mm{M}(\mv{x},t)$ satisfying 
> $$
> \ddtt\mm{M}+2\mysym(\mm{M}\pptfrac{\mv{f}}{\mv{x}}) \le -2\alpha\mm{M}
> $$
> with an exponential rate $\alpha\in\R_{>0}$.

^4b076f

`\begin{proof}`
Consider the following Lyapunov function $V:=V(\mv{x},\delta\mv{x},t)=\tfrac{1}{2}\delta\mv{x}^\top\delta\mv{x}$. The time derivative of $V$ is
$$
\ddtt V = 
\delta\mv{x}^\top 
\left(\ddtt\mm{M}+2\mysym(\mm{M}\pptfrac{\mv{f}}{\mv{x}}) \right)
\delta\mv{x}
\le
-2\alpha
\delta\mv{x}^\top 
\mm{M}
\delta\mv{x}
=
-\alpha V
.
$$
Therefore, $\norm{\delta\mv{x}}\to 0$ as $t\to 0$ with an exponential rate $\alpha$.
`\end{proof}`
Now, consider a perturbed system as follows:
$$
\ddtt\mv{x}
=
\mv{f}(\mv{x},t)
+ \mv{d}(\mv{x},t)
.
$$

^aa0442

> [!theorem] Partial Contraction
> Consider a nonlinear system of the form $\ddtt\mv{x} = \mv{f}(\mv{x},\mv{x},t)$ and the virtual system $\ddtt\mv{q} = \mv{f}(\mv{q},\mv{x},t)$ are contracting with respect to $\mv{q}$. If a particular solution of virtual system verifies a smooth specific property, then all solutions of the original system verifies this property exponentially.

Assume that the contraction metric $\mm{M}$ satisfies the contraction condition presented in  [[#^4b076f]] . Then the following theorem can be used to ...
 
> [!Theorem] Contraction (Perturbed System)
> If the system [[#^aa0442]]  is contracting with contraction metric $\mm{M}$ and decay rate $\alpha$, *i.e.,* $\mm{M}\mathbb{A}+\mathbb{A}\mm{M}-\mm{M}\mm{B}\mm{R}^{-1}\mm{B}^\top\mm{M}\le-2\alpha\mm{M}$, and $\Theta\mv{d}\in\mathcal{L}_\infty$ and $\exists\underline{m},\overline{m}\in\R_{>0}$, and $\exists\overline{d}\in\R_{\ge0}$, subject to $\overline{d}=\sup_{\mv{x},t}\norm{\mv{d}}$, and $\underline{m}\mm{I}\le \mm{M}\le\overline{m}\mm{I}$.
> Then, we have the following relation: 
> $$
> \norm{\xi_1(t)-\xi_0(t)}
> \le
> 
> $$
> .....

^c3ed0a

`\begin{proof}`
According the [[#^4b076f]], we have $\mm{M}=\mm{\Theta}^\top\mm{\Theta}$.  Then, the following relation can be obtained:
$$
\begin{aligned}
\ddtt \norm{\mm\Theta\pptfrac{\mv{q}}{\mu}}
=&
\tfrac{1}{2}
\left(
\pptfrac{\mv{q}}{\mu}^\top
\mm M
\pptfrac{\mv{q}}{\mu}
\right)^{-\tfrac{1}{2}}
\left[
\pptfrac{\mv{q}}{\mu}^\top
\ddtt \mm M
\pptfrac{\mv{q}}{\mu}
+
\pptfrac{\mv{q}}{\mu}^\top
2\mysym(\mm\Theta\pptfrac{\mv{f}}{\mv{q}})
\pptfrac{\mv{q}}{\mu}
+
\pptfrac{\mv{q}}{\mu}^\top
2
\mm M \mv{d}
\right]
\\
\le&
\tfrac{1}{2}
\left(
\pptfrac{\mv{q}}{\mu}^\top
\mm M
\pptfrac{\mv{q}}{\mu}
\right)^{-\tfrac{1}{2}}
\left[
-2\alpha
\pptfrac{\mv{q}}{\mu}^\top
\mm M
\pptfrac{\mv{q}}{\mu}
+
\pptfrac{\mv{q}}{\mu}^\top
2\mm M \mv{d}
\right]
\\
\le &
-\alpha
\norm{\mm\Theta\pptfrac{\mv{q}}{\mu}}
+
\norm{\mm\Theta\pptfrac{\mv{q}}{\mu}}^{-1}
\norm{\mm\Theta\pptfrac{\mv{q}}{\mu}}
\norm{\mm\Theta \mv{d}}
\\
\le &
-\alpha
\norm{\mm\Theta\pptfrac{\mv{q}}{\mu}}
+
\norm{\mm\Theta \mv{d}}
.
\end{aligned}
$$
By integrating 

`\end{proof}`
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
}_{=:\mv{\phi}(\mv{u})}
+ \mv{d}
,
\end{aligned}
$$

^cdaa47

The control law $\mv{u}$ is defined with contraction metric $\mm{M}$ as follows:
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

If the perturbed system [[#^9d8748]] is contracting with respect to [[#^4b076f]],  Accordingly to [[#^c3ed0a]], the steady error of [[#^9d8748]] (*i.e.,* $\norm{\mv{x}-\mv{x}_d}$) can be obtained as follows:
$$
\norm{\mv{x}-\mv{x}_d}
\le 
...
+
\sup(\mm{B}\mv{\phi}+\mv{d})
\tfrac{1}{\alpha}
\sqrt{
\tfrac{\overline{m}}{\underline{m}}
}
(1-\exp(-\alpha t))
.
$$
Therefore, saturation should be handled and $\mm{M}$ should be carefully selected such that $\tfrac{\overline{m}}{\underline{m}}$.
Then, we have the following optimization problem:
$$
\ddtt\mm{M}+2\mysym(\mm{M}\mathbb{A})+\mm{M}\mm{B}\mm{R}^{-1}\mm{B}^\top\mm{M}\le -2\alpha\mm{M}
$$

^92bbd5

To reformulate [[#^92bbd5]] into LMI from BMI, we pre-and post-multiply $\mm{W}:=\mm{M}^{-1}$ to [[#^92bbd5]]. 

### 3.2 Saturation Handling
We consider the following input saturation constraint:
$$
\norm{\mv{u}}
=
\norm{\mv{u}_d-\mm{R}^{-1}\mm{B}^\top\mm{M}(\mv{x}-\mv{x}_d)}
\le
\overline{u}
.
$$
However, ...
$$
\left(\mv{u}_d-\mm{R}^{-1}\mm{B}^\top\mm{M}(\mv{x}-\mv{x}_d\right)^\top
\left(\mv{u}_d-\mm{R}^{-1}\mm{B}^\top\mm{M}(\mv{x}-\mv{x}_d\right)
\le
\overline{u}^2
.
$$
Therefore, ...
$$
\begin{bmatrix}
\overline{u}^2 & \left(\mv{u}_d - \mm{R}^{-1}\mm{B}^\top\mm{M}(\mv{x}-\mv{x}_d\right)^\top \\
\star & \mm{I}
\end{bmatrix}
\succeq
0

$$


Multiply $\mm{T}=\mydiag(\mm{I},\mm{W}(\mm{R}^{-1}\mm{B}^{\top})^{-1})$ 
$$
\begin{bmatrix}
\overline{u}^2 & \left(\mm{W}(\mm{R}^{-1}\mm{B}^\top)^{-1}\mv{u}_d - (\mv{x}-\mv{x}_d\right)^\top \\
\star & 
\underbrace{
\mm{W}(\mm{R}^{-1}\mm{B}^\top)^{-1}(\mm{R}^{-1}\mm{B}^\top)\mm{W}
}_{=:\mm{X}}
\end{bmatrix}
\succeq
0
$$
Then the constraint is ...
$$
\begin{aligned}
\begin{bmatrix}
\overline{u}^2 & \left(\mm{W}(\mm{R}^{-1}\mm{B}^\top)^{-1}\mv{u}_d - (\mv{x}-\mv{x}_d\right)^\top \\
\star & \mm{X}
\end{bmatrix}
\succeq
0
,\quad
&
\mm{X}
\succeq
\mm{W}(\mm{R}^{-1}\mm{B}^\top)^{-1}(\mm{R}^{-1}\mm{B}^\top)\mm{W}

\end{aligned}
$$

## 4 Numerical Validation



|                       | C1     | C2              | C3              |
| --------------------- | ------ | --------------- | --------------- |
| Contraction Condition | yes    | yes             | yes             |
| Saturation Handling   | in LMI | penalty (large) | penalty (small) |
|                       |        |                 |                 |
Lorenz system
$$
\ddtt
\begin{pmatrix}
x\\y\\z
\end{pmatrix}
=
\begin{pmatrix}
\sigma(y-x)
\\
x(\rho-z)-y
\\
xy-\beta z
\end{pmatrix}
+
\begin{pmatrix}
u_x\\u_y\\u_z
\end{pmatrix}
+
\begin{pmatrix}
d_x\\d_y\\d_z
\end{pmatrix}
$$
SDC
$$
\mathbb{A}
=
\begin{bmatrix}
\sigma&-\sigma&0
\\
\rho-\tfrac{1}{2}(z+z_d) & -1 & -\tfrac{1}{2}(x+x_d)
\\
\tfrac{1}{2}(y+y_d) & \tfrac{1}{2}(x+x_d) & -\beta
\end{bmatrix}
$$

### 4.1 Scenario 1: Contraction Property
- 본 시나리오에서는 contracting system의 거동을 알아본다
- 아래와 같은 거동을 보이길 기대한다
	- C1은 saturation 방지가 LMI에서 가능
	- C2와 C3는 penalty를 이용해서 간접적으로 함
		- 즉, penalty로 인한 control input에 bias가 존재
			- C2가 더 작은 control input을 사용할 것임
	- 이를 통해 C1는 bias없이 제어 성능이 더 좋다는 것을 보임
- 이를 위해 아래와 같은 시나리오를 선정한다
	- Contraction property에 집중하기 위하여  포화는 걸리지 않도록한다
	- 계단 desired input을 사용
	- 외란을 지속적으로 적용
### 4.2 Scenario 2: Saturation Handling
- 본 시나리오는 포화 해결에 집중
- 아래와 같은 거동을 보이길 기대
	- C1은 LMI 내부에 saturation handling이 포함
		- 제약조건을 범하지 않는 contraction condition을 만족하는 contraction metric 사용
		- 물론 saturation을 범하지 않음으로 성능이 낮아질수있음
	- C2는 penalty를 이용함으로 saturation violation 빈도가 낮음
		- 성능이 C1보다 약간 좋을 수도 있음
	- C3는 penalty가 작기 때문에 satuation이 크게 발생할 수있음
		- 이로인해 control input이 과도하게 커짐
- 


## 5 Conclusion
In this study, we 
