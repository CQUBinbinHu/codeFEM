# Robust shape and topology optimization considering geometric uncertainties using the sensitivity analysis within level set method

**Abstract**

> By applying the sensitivity analysis with respect to the level set function, the descent direction of the level set function is obtained. 

## Introduction

## Robust optimization model with implicit method description

**Robust optimization  formulation**

Based on the work [^7] of (Zhang 2017), the robust shape and topology optimization problem are formulated as:
$$
\min _{\Omega} {\mu [f(\Omega(\xi),u(\xi))] + \beta_0 \sigma [f(\Omega(\xi),u(\xi))] }
$$
*s.t.*
$$
\mu[g_i(\Omega(\xi),u(\xi))] + \beta_i \sigma[g_i(\Omega(\xi),u(\xi))] \le 0 \\
\int_\Omega H(\phi) \le \overline V
$$
where the state variable $u(\xi)$ is governed by the state equations:
$$
-\nabla \cdot (Ae(u(\xi))) = f \ \ \ in\ \Omega(\xi) \\
(Ae(u(\xi)))\cdot n = g \ \ \  on \ \Gamma_N(\xi) \\
u = 0 \ \ \ on \ \Gamma_D(\xi)
$$
Employing the **implicit method** [^4], the optimization problem can be described as:
$$
\min _{\phi} {\mu [f(\varphi(\xi),u(\xi))] + \beta_0 \sigma [f(\varphi(\xi),u(\xi))] }
$$
s.t.
$$
\mu[g_i(\varphi(\xi),u(\xi))] + \beta_i \sigma[g_i(\varphi(\xi),u(\xi))] \le 0 \\
\int_D H(\phi) \le \overline V
$$
with the variation form of state equation constraint:
$$
a(u,v,\phi) = L(v,\phi), \ \forall v\in U
$$
By applying the augmented Lagrangian method, the robust optimization problem is converted into:
$$
\min \hat L(\varphi,\lambda,Y,r)
$$
where,
$$
\hat L(\phi,\lambda,Y,r) = \mu [f(\varphi(\xi),u(\xi))] + \beta_0 \sigma [f(\varphi(\xi),u(\xi))] \\
+ \sum_i \lambda_i (h_i(\phi,Y)+ Y_i^2) \\
+ r \sum_i (h_i(\phi,Y) + Y_i^2)^2 \\
h_i =\mu[g_i(\varphi(\xi),u(\xi))] + \beta_i \sigma[g_i(\varphi(\xi),u(\xi))] \\
$$
Here, the $Y$ is the vector of slack variables. The penalty parameter $r$ is problem oriented. For sake of convenience, equation (8) is equivalent to
$$
\hat L(\phi,\lambda,r) = \mu [f(\varphi(\xi),u(\xi))] + \beta_0 \sigma [f(\varphi(\xi),u(\xi))] \\
+ \sum_i \lambda_i (\alpha_i) 
+ r \sum_i (\alpha_i)^2 \\
\alpha_i = \max\{ h_i(\phi),-\frac{\lambda_i}{2 r} \} \\
$$
Then the **PCE** is used for the discretization of each response with the following formulas :
$$
\mu[g] \approx \sum_j^q w^j g(\xi^{(j)})
$$

$$
\sigma[g] \approx \sqrt{\sum_i\frac{1}{\gamma_i}(\sum_j w^{j}g(\xi^{(j)})\Psi_i(\xi^{(j)}))^2}
$$

substitute (6)(7) into (5), and then take the sensitivity analysis,
$$
\frac{\part \hat L}{\part \phi} = \sum w^j\frac{\part f(\xi^j)}{\part \phi} \\

+ \beta_0 \frac{1}{\sigma(f)} \sum \sum \frac{1}{\gamma_i} w^j \Psi_i \frac{\part f(\xi^j)}{\part \phi} \\

+ \sum \lambda_k (\sum w^j\frac{\part g_k(\xi^j)}{\part \phi}+ \frac{1}{\sigma(g_k)} \sum \sum \frac{1}{\gamma_i} w^j \Psi_i \frac{\part g_k(\xi^j)}{\part \phi}) \\

+ 2r_k\sum (\mu[g_k] + \beta_i \sigma[g_k] + Y^2_i)(\sum w^j\frac{\part g_k(\xi^j)}{\part \phi}+ \frac{1}{\sigma(g_k)} \sum \sum \frac{1}{\gamma_i} w^j \Psi_i \frac{\part g_k(\xi^j)}{\part \phi})
$$
the obtained formula (11) can be represented as:
$$
\frac{\part \hat L}{\part \phi} = \sum_{j=0}^k A_j \frac{\part F_j}{\part \phi}
$$
where the response function $F_j$ can be generally expressed as:
$$
F = \int_D G(u,x)H(\varphi) dx + \int_D M(u,x)\delta(\varphi) \lVert\nabla\varphi\rVert dx
$$
with the constraint 
$$
-\nabla \cdot (Ae(u)) = f \ \ \ in\ \Omega(\varphi) \\
(Ae(u))\cdot n = g \ \ \  on \ \Gamma_N \\
u = 0 \ \ \ on \ \Gamma_D
$$


## Sensitivity analysis with respect to the level set function

In order to apply a gradient method to the minimization of the objective function (8), we need to analyze the variation of $L$ with respect to $\phi$. Also, we have the equation (12) which implies that we only need to take analysis on each $F_j$.

Following the techniques of Allaire [^2] , the Lagangian is introduced:
$$
L(\phi,v,q) = \int_D G(u,x) H(\phi) dx 
+ \int_D M(u,x)\delta(\phi) |\nabla \phi| dx \\
+ \int_D Ae(v)\cdot e(q) H(\phi) dx
- \int_D q\cdot f H(\phi) dx \\
- \int_D q\cdot g \delta(\phi) |\nabla \phi| dx \\
- \int_D q\cdot A(e(v)) \cdot n \delta(\phi) |\nabla \phi| dx
- \int_D v\cdot A(e(q)) \cdot n \delta(\phi) |\nabla \phi| dx
$$
According to the Lagrangian stationarity,
$$
\frac{\part L}{\part q}(\phi,u,p)(\eta) = 0
$$

$$
\frac{\part L}{\part v}(\phi,u,p)(\eta) = 0
$$

the proposed formula (17)  induce the objective function and state equation and boundary conditions. While the formula (18) gives the adjoint state equation with boundary conditions.

Notes that, the formula (19) includes two types of integral concerning the level set function. Such that, we only need to analysis the **Frechet derivative** of the following two types:
$$
\hat G(\phi) = \int_D G(x)H(\phi) dx
$$

$$
\hat M(\phi) = \int_D M(x)\delta(\phi) |\nabla \phi| dx
$$

After analysis, we can get:
$$
\frac{\part \hat G}{\part \phi}(\eta) = \int_D G(x) \delta (\phi) \eta dx
$$

$$
\frac{\part \hat M}{\part \phi}(\eta) = \int_D - \frac{\nabla M \cdot \nabla \phi}{|\nabla \phi|} \delta (\phi) \eta dx
$$

Thus, by substituting Eqs. (24)(25) into (19), the derivative of $L(\phi,v,q)$ along  direction $\eta$ can be represented as:
$$
\frac{\part L}{\part \phi}(\phi,u,p)(\eta) = \int_D K^* \delta(\phi) \eta dx
$$
where 
$$
K^* =  P - \frac{\nabla Q \cdot \nabla \phi}{|\nabla \phi|} \\
P =G(u,x) + Ae(v) \cdot e(q) - q \cdot f \\
Q =  M(u,x) 
- q\cdot g  
- q\cdot A(e(v)) \cdot n 
- v\cdot A(e(q)) \cdot n 
$$
Taking $\eta(x) = -K^*(x)$ as a descent direction.

## Algorithm

### Robust Optimization Problem within geometry uncertainties

According to the previous analysis, an algorithm for solving the robust optimization problem with the augmented Lagrangian method is proposed at the following.

1. Set the initial information $\phi^{(0)}, \lambda^{(0)},Y^{(0)},r^{(0)} $

2. Minimize $\hat L(\phi,\lambda^{(k)} ,Y,r)$ from initial $\phi^{(k)}$ , find the $\phi^{(k+1)}$, $Y^{(k+1)}$.

3. Check for convergence of $\lambda^{(k)}$, and $\phi^{(k)}$. 

   **If** it is convergent, take $\phi^* = \phi^{(k)}$ and **stop**. **Else**, the procedure continue to **next step**.

4. Set $\lambda_j^{(k+1)} = \lambda_j^{(k)} + 2 r^{(k)}h_j(\phi^{(k)})$

5. Set $r^{(k+1)} = cr^{(k)}$

6. If $r^{(k+1)}>r_{max}$, then set $r^{k+1} = r_{max}$

7. Set k = k+1 and return to **step 2.**

During the proposed algorithm, the **Step 2**  involves the minimization process which is described at the following :

1. Take $N$ samples of input parameters $\xi^{(j)}, j=1,...,N$ 

   where $\xi_i= [ \xi_i^{(j)},...,\xi_i^{(j)} ]^T$

2. For each $\xi^{(j)}$, calculate the perturbation based on K-L basis:

   $\delta\phi^{(j)} = A \mathcal{F} (\sum\lambda_i\xi_i^{(j)}\Psi_i) $

   and the **Sample** Level Set functions:

   $\psi^{(j)} = \phi + \delta\phi^{(j)}$

3. For each Sample, solve the state equation (3) and get the displacement $u^{(j)}$ 

4. Calculate each response $g_k(\psi^{(j)},u^{(j)})$

5. Get the descent direction $\eta = -K^*$

6. Do one-dimension search along the descent direction get the updating step

   $\Delta\phi = \epsilon*\eta$

7. Update the level set function through

   $\phi^* = \phi + \Delta \phi$

Notes that the **Step-1**, **Step-5** need to be specially handled.

## Compliance Optimization with the only volume constraint

In order to verify the feasibility of our method, we first solve a simpler problem, that is, the optimization within compliance objective function under volume constraint.

$$
J(\phi) = l \int_D H(\phi) dx + \int_D fu_y \delta(\phi) |\nabla \phi| dx \\
 =  \int_D (l + Ae(u)\cdot e(u)) H(\phi) dx
$$
where the Lagrange multiplier $l = 100.$
$$
J'(\eta) = l \int_D \delta (\phi) \eta dx + \\
\int_D - \frac{\nabla M \cdot \nabla \phi}{|\nabla \phi|} \delta (\phi) \eta dx
$$




















## References

------

[^1]: Rao SS. Engineering Optimization: Theory and Practice (4th edn). Wiley: New York, 2009.
[^2]: Allaire, Grégoire, François Jouve, and Anca-Maria Toader. “Structural Optimization Using Sensitivity Analysis and a Level-Set Method.” *Journal of Computational Physics* 194, no. 1 (February 10, 2004): 363–93. https://doi.org/10.1016/j.jcp.2003.09.032.
[^3]: Sethian JA. Level Set Methods and Fast Marching Methods: Evolving Interfaces in Computational Geometry, Fluid Mechanics, Computer Vision, and Materials Science. Cambridge University Press: Cambridge, 1999.
[^4]: Wang, Michael Yu, Xiaoming Wang, and Dongming Guo. “A Level Set Method for Structural Topology Optimization.” *Computer Methods in Applied Mechanics and Engineering* 192, no. 1 (January 3, 2003): 227–46. https://doi.org/10.1016/S0045-7825(02)00559-5.
[^5]: Huang, S. P., S. T. Quek, and K. K. Phoon. “Convergence Study of the Truncated Karhunen–Loeve Expansion for Simulation of Stochastic Processes.” *International Journal for Numerical Methods in Engineering* 52, no. 9 (2001): 1029–43. https://doi.org/10.1002/nme.255.
[^6]: Hughes, T. J. R., J. A. Cottrell, and Y. Bazilevs. “Isogeometric Analysis: CAD, Finite Elements, NURBS, Exact Geometry and Mesh Refinement.” *Computer Methods in Applied Mechanics and Engineering* 194, no. 39 (October 1, 2005): 4135–95. https://doi.org/10.1016/j.cma.2004.10.008.
[^7]: Zhang, Wenbo, and Zhan Kang. “Robust Shape and Topology Optimization Considering Geometric Uncertainties with Stochastic Level Set Perturbation: ROBUST TOPOLOGY OPTIMIZATION CONSIDERING GEOMETRIC UNCERTAINTIES.” *International Journal for Numerical Methods in Engineering* 110, no. 1 (April 6, 2017): 31–56. https://doi.org/10.1002/nme.5344.

