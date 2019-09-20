# CodeDoc

> This is a document for the c++ codes for solving shape optimization problems.

## Finite Element Method within Isoparametric transformation

### Elastic Matrix

$$
\mathbf{D} = D_0 \begin{bmatrix}
1 & \mu_0 & 0 \\
\mu_0 & 1 & 0 \\
0 & 0 & \frac{1-\mu_0}{2}
\end{bmatrix}
$$

with the parameters:
$$
D_0 = \frac{E_0}{1-\mu_0^2}
$$
where,
$$
E_0 = E \\
\mu_0 = \mu
$$
for the plain stress problem;
$$
E_0 = \frac{E}{1-\mu^2} \\
\mu_0 = \frac{\mu}{1-\mu}
$$
for the plain strain problem.

### Isoparametric Shape Function

#### for linear function**

$$
N_i = \frac{1}{4} (1 + \hat\xi_i)(1 + \hat\eta_i)
$$
where $\hat\xi_i = \xi_i \xi ,\  \hat\eta_i = \eta_i \eta$

with the nodes : $(\xi_1,\eta_1),(\xi_2,\eta_1),(\xi_2,\eta_2),(\xi_1,\eta_2)$

#### for quadratic function

$$
N_5 = \frac{1}{2} (1-\xi^2)(1-\eta) \\
N_6 = \frac{1}{2} (1-\eta^2)(1+\xi) \\
N_7 = \frac{1}{2} (1-\xi^2  )(1 +\eta) \\
N_8 = \frac{1}{2} (1-\eta^2)(1-\xi) \\
N_1 = \hat N_1 - \frac{1}{2}N_5 - \frac{1}{2}N_8 \\
N_2 = \hat N_2 - \frac{1}{2}N_5 - \frac{1}{2}N_6 \\
N_3 = \hat N_3 - \frac{1}{2}N_6 - \frac{1}{2}N_7 \\
N_4 = \hat N_4 - \frac{1}{2}N_7 - \frac{1}{2}N_8
$$
