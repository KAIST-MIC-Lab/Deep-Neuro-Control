## System Declaration
Let say we have the following system,
$$
\ddtt{\mv x} = \mv f(\mv x)+\mv g(\mv x)\mysat(\mv u)+\mv{d}
,
$$
where .... (omitted here).

The control objective is to design a tracking controller to track some desired trajectory $\mv{x}_d$. In this study, we are assuming that desired tacking controller $\mv{u}_d$ is already given. In other words, the following desired system is given as follows:
$$
\ddtt{\mv x}_d 
= 
\mv f(\mv x_d)
+
\mv g(\mv x_d)
\mv u_d
,
$$
where ... (omitted here). Note that the desired controller is designed not to exceed the input limits with respect to the zero-disturbances, *i.e.,* $\mv{d}=\mv{0}$.


