# Deep-Neuro Control with Contraction Theory

## Description
 
The objectives of this research are as follows:
- Mathematical stability analysis of the controller and estimator with deep neural networks using the contraction theory.
- Development of the controller and estimator with deep neural networks using the contraction theory.

## Tables of Contents

You can find draft [./manuscript.pdf](./manuscript.pdf) and [./doc/document.pdf](./doc/main.pdf).

```
├── README.md                   // you are here
├── manuscript.tex              // template of manuscript
├── dding_template              // pre-defined template
├── figures                     // general figures
├── public                      // do not touch
├── localRefs.bib               // bibtex file
├── docs                        // documentation
└── src                         // source code
    ├── script_simulation       // script sim. example
    │   ├── figures             //             figures
    │   └── results             //             results
    └── simulink_simulation     // simulink sim. example
        ├── figures             //               figures
        └── results             //               results
```

## How to start

> If you have questions, please, please, let the author knows.

Download this remote repository on your local device.
```
git clone https://gitlab.com/dding_friends/dding_research
```
Initialize submodule to download [Template](https://gitlab.com/dding_friends/dding_template) repository.
> This command will download very recent version of [Template](https://gitlab.com/dding_friends/dding_template)
```
git submodule init
git submodule update
```

I provide you the keywords that you can google what you need for what you want to do.

| Keywords 	    | Descriptions 	|
|---------	    |--------------	|
| branch       	| Want to      	|

<!-- 
Talbe generator
https://www.tablesgenerator.com/markdown_tables
https://insight.infograb.net/blog/2025/01/22/gitlab-wiki/
-->

## Example research simulation and manuscript

> This simulation in `./src ` provides simple feedback control example.
> Please, read `src/script_simulation/main.m` and `src/simulink_simulation/main_slx.m` for script and Simulink simulations, respectively.
> When you run those scripts, the results and figures shall be saved in the directories named `results` and `figures` (of course, you need to check setting in the scripts.).

Consider following control-affine system represented as
$$
    \dot {\boldsymbol{x}} = \boldsymbol{A}\boldsymbol{x} + \boldsymbol{B}\boldsymbol{u}
$$
where $\boldsymbol{x}=[x_1;x_2]\in\mathbb{R}^2$ denotes state variables, $\boldsymbol{u}=[u_1;u_2]\in\mathbb{R}^2$ denotes control input variables, and $\boldsymbol{A}=[0,1;-2,-3]\in\mathbb{R}^{2\times2}$ and $\boldsymbol{B}=I_2\in\mathbb{R}^{2\times2}$ denote system and input matrices, respectively.
The initial condition are $\boldsymbol x\vert_{t=0} = [0;0]^\top$ and $\boldsymbol u\vert_{t=0} = [0;0]^\top$.

Suppose that we have smooth reference trajectory of $\boldsymbol{x}$ as follows:
$$
    \boldsymbol{r} 
    = 
    \begin{bmatrix}
        r_1 \\
        r_2
    \end{bmatrix}
    =
    \begin{bmatrix}
        \sin(t)\\
        \cos(t)
    \end{bmatrix}
    \in\mathbb{R}^2
$$
A feedback control law can be designed as follows:
$$
    \boldsymbol{u} = -\boldsymbol{K}\boldsymbol{e}
$$
where $\boldsymbol{e}\triangleq \boldsymbol{x}-\boldsymbol{r}$ denotes tracking error and $\boldsymbol{K}=\text{diag}([-2,-3])\in\mathbb{R}^{2\times 2}$ denotes control gain matrix.

Then the simulation results are plotted like below.

<img src="src/script_simulation/figures/9-Feb-2025_22-32-23/Fig1.png" width="350">
<img src="src/script_simulation/figures/9-Feb-2025_22-32-23/Fig2.png" width="350">
<img src="src/script_simulation/figures/9-Feb-2025_22-32-23/Fig3.png" width="350">
<img src="src/script_simulation/figures/9-Feb-2025_22-32-23/Fig4.png" width="350">

## Authors

- Ryu Myeongseok [@DDingR](https://gitlab.com/DDingR)
- You Sesun [@yousesun95](https://gitlab.com/yousesun95)
- Choi Kyunhhwan 

## License

