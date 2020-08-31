## Abstract
The WENO (Weighted Essentially Non-Oscillatory) methods were developed after the ENO (Essentially Non-Oscillatory) methods in 1990s in order to solve hyperbolic PDEs with discontinuities and shocks. In this project, we will introduce the derivation of WENO scheme developed by Liu, Chan, and Osher in 1994, and the following scheme developed by Jiang and Shu in 1996 as well as their order of accuracy. We will use $r = 3$ as an example to show how the scheme is carried out. Then, we will briefly discuss WENO's linear stability analysis. A 1 dimensional numerical example will be presented in the end as well.

## Introduction
First developed by Harten, Engquist, Osher, and Chakravarthy in 1987, the Essentially Non-Oscillatory (ENO) methods are part of high-resolution schemes in finding the numerical solution of partial differential equations. Since hyperbolic PDEs may contain singularities in their solutions, the ENO methods are desirable because they are able to provide high-order accuracy when the function is smooth while avoiding Gibbs phenomenon at discontinuities. The ENO schemes could achieve such goals by reconstructing a piece-wise polynomial of the solution from its cell average. In other words, the procedure uses the smoothest stencil to approximate fluxes at cell boundaries to avoid oscillations near shocks.

As researchers later pointed out, the ENO methods, although uniformly high order accurate, can have certain drawbacks. Its reconstruction procedure from cell averages could be complicated and computationally expensive. In addition, its adaptive stencil could change from perturbations and might not be in the regions where solution is smooth. In 1994, a new version of ENO schemes was proposed by Liu, Osher, and Chan. This new modification is called weighted ENO (WENO). The main idea is that, unlike ENO that chooses the smoothest stencil for reconstruction, WENO will use a convex combination of all candidate stencils. A weight will be assigned to each candidate, which represents how much this stencil contributes to the final solution. The scheme will remain essentially non-oscillatory as ENO, while adding one more order of accuracy.

Two years later, an improved version of WENO method was proposed by Jiang and Shu. The previous WENO scheme requires weights on candidate stencils, which depends on their relative smoothness. Thus, following the previous WENO scheme, Jiang and Shu devised a new way of evaluating the smoothness of a stencil. This new smoothness indicator achieves the highest possible order of accuracy (fifth order) when $r = 3$. We will now present the derivations of two WENO schemes.



## WENO, version 1
### Problem set-up
Recall the hyperbolic conservation law $\bf{u}_t + div \bf{f}(\bf{u}) = 0$. In 1 dimension, the conservation law becomes
$$u_t + f(u)_x = 0.$$

For a uniform discretization in space $\cdots < x_0 < x_1 < x_2 < \cdots$, $\Delta x = x_{i+1} - x_i$. We use $I_i = [x_{i - 1/2}, x_{i+1/2}]$ to represent i-th cells. The cell average can be written as

$$\bar{u}_i =\frac{1}{\Delta x} \int_{I_i}u(x, t)dx
= \frac{1}{\Delta x} \int^{x_{i+1/2}}_{x_{i-1/2}}u(x, t)dx. $$

Integrating the conservation law over each cell, we obtain

$$\frac{d}{dt}\bar{u}_i(t) = - \frac{1}{\Delta x}(f(u(x_{i+1/2}, t)) - f(u(x_{i-1/2}, t))).$$

In order to evaluate $f(u(x, t))$ at each cell boundaries $x_{i+1/2}$, we can reconstruct a piecewise polynomial $R(x) = \{R_j(x)\}$ approximating $u(x, t)$ given the cell averages $\bar{u} = \{\bar{u}_j\}$. For a general flux, there are multiple methods used for flux-splitting, including Engquist-Osher, Godunov, Roe with entropy fix, and Lax-Friedrichs etc.

Putting the time discretization in a more compact form, we can write
$$
\begin{align*}
    \frac{d}{dt}\bar{u}_i(t) &\approx L_i(\bar{u})\\
    L_i(\bar{u}) &= -\frac{1}{\Delta x}[\hat{f}(R_i(x_{i+1/2}), R_{i+1}(x_{i+1/2}) - \\ &\hat{f}(R_{i-1}(x_{x_{i-1/2}}), R_i(x_{i-1/2}))],
\end{align*}
where $\hat{f}(R_i(x_{i+1/2}), R_{i+1}(x_{i+1/2})$ approximates $f(u(x_{i+1/2}, t))$, $\hat{f}(R_{i-1}(x_{x_{i-1/2}}), R_i(x_{i-1/2}))$ approximates $f(u(x_{i-1/2}, t))$.
$$

### Reconstruction
We denote the r candidate stencils by $S_k$, $k = 0, 1, \cdots, r-1$, where

$$S_k = \{I_{i+k-r+1}, I_{i+k-r+2}, \cdots, I_{i+k}\}$$

For example, when $r = 3$, the candidate stencils would be
$$
\begin{align*}
    S_0 &= \{I_{i - 2}, I_{i-1}, I_i\} = \{x_{i-5/2}, x_{i-3/2}, x_{i-1/2}, x_{i+1/2}\}\\
    S_1 &= \{I_{i - 1}, I_{i}, I_{i+1}\} = \{x_{i-3/2}, x_{i-1/2}, x_{i+1/2}, x_{i+3/2}\}\\
    S_2 &= \{I_{i}, I_{i+1}, I_{i+2}\} = \{x_{i-1/2}, x_{i+1/2}, x_{i+3/2}, x_{i+5/2}\}\\
\end{align*}
$$
The three candidate stencils are illustrated below. Each color represents a different stencil.

![Candidate Stencils](images/stencils.PNG)

Figure 1: Candidate Stencils.

We would then obtain the interpolation polynomials on each of the stencils such that

$$
\bar{u}_i = \frac{1}{\Delta x_i} \int_{I_i} p_j (x)dx.
$$

Again, using $r = 3$ as an example, we have three corresponding quadratic polynomials,
$$
\begin{align*}
    p_i(x) &= \frac{\bar{u}_i - 2\bar{u}_{i-1} + \bar{u}_{i-2}}{2 \Delta x^2} (x - x_{i-1})^2 + \frac{\bar{u}_i - \bar{u}_{i-2}}{2 \Delta x}(x - x_{i-1}) + \bar{u}_{i-1} - \frac{\bar{u}_i - 2\bar{u}_{i-1} + \bar{u}_{i-2}}{24}\\
    p_{i+1}(x) &= \frac{\bar{u}_{i+1} - 2\bar{u}_{i} + \bar{u}_{i-1}}{2 \Delta x^2} (x - x_{i})^2 + \frac{\bar{u}_{i+1} - \bar{u}_{i-1}}{2 \Delta x}(x - x_{i}) + \bar{u}_{i} - \frac{\bar{u}_{i+1} - 2\bar{u}_{i} + \bar{u}_{i-1}}{24}\\
    p_{i+2}(x) &= \frac{\bar{u}_{i+2} - 2\bar{u}_{i+1} + \bar{u}_{i}}{2 \Delta x^2} (x - x_{i+1})^2 + \frac{\bar{u}_{i+2} - \bar{u}_{i}}{2 \Delta x}(x - x_{i+1}) + \bar{u}_{i+1} - \frac{\bar{u}_{i+2} - 2\bar{u}_{i+1} + \bar{u}_{i}}{24}.\\
\end{align*}
$$
Let the approximated value $u^{k}_{i+1/2} = p_k(x_{i+1/2})$, the expressions above can be simplified to
$$
\begin{align*}
    u^0_{i+1/2} &= \frac{1}{3} \bar{u}_{i-2} - \frac{7}{6}\bar{u}_{i-1} + \frac{11}{6}\bar{u}_{i}\\
    u^1_{i+1/2} &= -\frac{1}{6} \bar{u}_{i-1} + \frac{5}{6}\bar{u}_{i} + \frac{1}{3}\bar{u}_{i+1}\\
    u^3_{i+1/2} &= \frac{1}{3} \bar{u}_{i} + \frac{5}{6}\bar{u}_{i+1} - \frac{1}{6}\bar{u}_{i+2}\\
\end{align*}
$$
Note that this approximation is third order accurate, namely, $u^k_{i+1/2} - u(x_{i+1/2}) = \mathcal{O}(\Delta x^3)$. However, we are not satisfied with interpolation on a single stencil, we want to achieve more by using all of them together. In other words, we aim to find a convex combination of the polynomials that follows the conservation law and it is essentially non-oscillatory. Since the hyperbolic PDEs may contain shocks and discontinuities, we should use the smoother part of the region to form our solution. The smoother the stencil is, the heavier weight it will have on the combination. The difference between two versions of the WENO methods is largely based on how the smoothness indicator is defined.

### First Smoothness Measurement
The indicator of smoothness (IS, or $\beta$'s in some literatures) on $j$-th stencil proposed by Liu, Osher, and Chan is defined to be
$$
IS_j = \sum^{r-1}_{l=1} \left ( \sum^l_{k=1} (\triangle^{r-1}[u_{j-r+k}])^2 \right)/l,
$$

where $\triangle[u]$ comes from a table of differences of $\{\bar{u}_i\}$ on $S_i$ and satisfies $\triangle[u_l] = u_{l+1} - u_l, \triangle^k[u_l] = \triangle^{k-1}[u_{l+1}] - \triangle^{k-1}[u_l]$. This way, the smoothness indicator becomes the summation of all averages of squared same order differences. When $r = 3$,
$$IS_j= \frac{(\triangle[u_{j-2}])^2 + (\triangle[u_{j-1}])^2}{2} + (\triangle^2[u_{j-2}])^2.$$
When the true solution $u(x,t)$ is discontinuous on $S_j$, $IS_j \approx \mathcal{O}(1)$; if $u(x, t)$ is continuous on $S_j$, $IS_j \approx \mathcal{O}(\Delta x^2)$.

The convex combination $R_j(x)$ that comes from $r$ interpolating polynomials ($p_k$ obtained before) on $r$ stencils is of the form
$$
\begin{align*}
    R_j(x) &= \sum^{r-1}_{k=0} \omega_k p_{j+k}(x)\\
    \omega_k &= \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l},
\end{align*}
$$
where $\sum^{r-1}_{k=0} \omega_k = 1, k = 0, 1, \cdots, r-1$ and $\alpha^j_k >0$. The ENO property is satisfied if the coefficients $\frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} = \mathcal{O}(1)$ if the stencil $S_{j+k}$ is in the smooth regions, and $\frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} \leq \mathcal{O}(\Delta x^r)$ if the stencil is in the discontinuous region of the solution.

Let us consider the coefficients defined as
$$
\alpha^j_k = \frac{C^j_k}{(\epsilon+ IS_{j+k})^r},
$$
where $\epsilon = (10^{-5} \sim 10^{-7})$ is a small positive number added here to ensure that the denominator will not become zero because $IS$ could potentially be zero. The coefficients $C_k^r$'s can be found by computing the weights of polynomials on each stencil to form the polynomial across all candidate stencils, ie.
$$
Q(x_{i+1/2}) = \sum^{r-1}_{k=0}C_k^r p_k(x_{i+1/2}),
$$
where $Q(x)$ is the interpolating polynomial on the largest stencil.

For example, when $r = 3$, we have 3 candidate stencils and the largest stencil that includes all contains 5 cells $\{I_{i-2}, I_{i-1}, I_{i}, I_{i+1}, I_{i+2}\}$. $Q(x)$, evaluated at $x_{i+1/2}$, has the expression
$$
Q(x_{i+1/2}) = \frac{1}{30} \bar{u}_{i-2} - \frac{13}{60}\bar{u}_{i-1} + \frac{47}{60} \bar{u}_i + \frac{9}{20} \bar{u}_{i+1} - \frac{1}{20}\bar{u}_{i+2}
$$

Combined with previous expressions for interpolation polynomials $p_k$'s we found before, it is easy to see that the coefficients are given by $C_0^3 = \frac{1}{10}, C_1^3 = \frac{6}{10}, C_2^3 = \frac{3}{10}$. When $r = 2$, $C_0^2 = \frac{1}{3}, C_1^2 = \frac{2}{3}$. The optimal weight for $C_k^r$ leads to one order of improvement in accuracy.  Now it is clear to see that
$$
\begin{align*}
    \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} &= \mathcal{O}(1), \text{when stencil is smooth} \\
    \frac{\alpha^j_k}{\sum^{r-1}_{l=0}\alpha^j_l} &\leq \max(\mathcal{O}(\epsilon^r), \mathcal{O}((\Delta x)^r)), \text{when stencil meets discontinuities}.
\end{align*}
$$
Thus, these weights $\{\alpha^j_k\}^{r-1}_{k=0}$ satisfies the ENO properties.

Going back to our example when $r = 3$, the reconstructed solution is expressed as
$$
R_j (x) = \omega_0 p_j(x) + \omega_1 p_{j+1}(x) + \omega_2 p_{j+2}(x),
$$
where $p_k(x)$'s are calculated above and the weights are given by
$$
\omega_0 = \frac{\alpha_0^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j}, \omega_1 = \frac{\alpha_1^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j},
\omega_2 = \frac{\alpha_2^j}{\alpha_0^j + \alpha_1^j + \alpha_2^j}
$$
and
$$
\alpha_0^j = \frac{1}{10(\epsilon + IS_j)^3},
\alpha_1^j = \frac{6}{10(\epsilon + IS_{j+1})^3},
\alpha_2^j = \frac{3}{10(\epsilon + IS_{j+2})^3}.
$$

The interpolating method described above satisfies
$$
\begin{align*}
    &u(x, t) = R_j(x) + \mathcal{O}((\Delta x)^r)\\
    &u(x_j^*, t) = R_j(x_j^*) + \mathcal{O}((\Delta x)^{r+1}),
\end{align*}
$$
where $x_j^*$'s are chosen to be points at the cell boundaries. For a general upwind method, when $f'(R(x)) > 0$,
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(R_j(x_{j+1/2})) - f(R_{j-1}(x_{j-1/2}))].
$$
The points chosen are $x_j^* = x_{j+1/2}, x_{j-1}^* = x_{j - 1/2}$. Following the reconstruction scheme presented above, we obtain
$$
\begin{align*}
    &R_i(x_{j+1/2}) - u(x_{j+1/2}, t) = \mathcal{O}(\Delta x^{r+1})\\
    &R_{j-1}(x_{j-1/2}) - u(x_{j-1/2}, t) = \mathcal{O}(\Delta x^{r+1}).
\end{align*}
$$
Hence,
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(u(x_{j+1/2}, t)) - f(u(x_{j-1/2}, t))] + \mathcal{O}(\Delta x^{r+1}).
$$
Similarly, for $f'(R(x)) < 0$, we have
$$
L_j(\bar{u}) = -\frac{1}{\Delta x}[f(R_{j+1}(x_{j+1/2})) - f(R_{j}(x_{j-1/2}))].
$$
Choosing $x_j^* = x_{j-1/2}, x_{j+1}^* = x_{j+1/2}$, we achieve the same order of accuracy. Therefore, in teh smooth regions, the spatial discretization $L$ is able to approximate $\frac{d \bar{u}}{dt}$ to the $r+1$-th order.




## WENO, version 2
Jiang and Shu later presented a new smoothness measurement for the candidate stencils and subsequently a different set of weights for the ENO scheme. Total variation is a good indicator for smoothness, hence, they proposed to minimize total variation of the approximation. Given $p_k(x)$ the interpolation polynomial on stencil $S_k$, the new smoothness indicator is defined as

$$
IS_k = \sum^{r-1}_{l=1} \int^{x_{j+1/2}}_{x_{j-1/2}} (\Delta x)^{2l-1}\left( \frac{d^l}{dx^l}p_k(x) \right)^2 dx.
$$

This new definition is a sum of the $L^2$ norms of all the derivatives of $p_k(x)$ over cell $I_i$. $(\Delta x)^{2l-1}$ is added in order to remove the $\Delta x$'s that arise while taking derivatives of the polynomials.
When $r = 3$, we obtain
$$
\begin{align*}
    IS_0 &= \frac{13}{12}(\bar{u}_{i-2} - 2\bar{u}_{i-1} + \bar{u}_i)^2 + \frac{1}{4}(2 \bar{u}_{i-2} - 4\bar{u}_{i-1} + \bar{u}_i)^2\\
    IS_1 &= \frac{13}{12}(\bar{u}_{i-1} - 2\bar{u}_{i} + \bar{u}_{i+1})^2 + \frac{1}{4}(3 \bar{u}_{i-1} - \bar{u}_{i+1})^2\\
    IS_2 &= \frac{13}{12}(\bar{u}_{i} - 2\bar{u}_{i+1} + \bar{u}_{i+2})^2 + \frac{1}{4}(\bar{u}_{i} - 4\bar{u}_{i+1} + \bar{u}_{i+2})^2
\end{align*}
$$
The weights follow previous form
$$
\omega_k = \frac{\alpha^j_k}{\sum^{r-1}_{l=0} \alpha^j_l}, \alpha^j_k = \frac{C_k^j}{(\epsilon + IS_k)^r}.
$$

In the smooth regions, using Taylor expansion, we can get
$$
\begin{align*}
    IS_0 &= \frac{13}{12}(f''(u_i)\Delta x^2)^2 + \frac{1}{4}(2f'(u_i) \Delta x- \frac{2}{3}f'''(u_i)\Delta x^3)^2 + \mathcal{O}(\Delta x^6)\\
    IS_1 &= \frac{13}{12}(f''(u_i)\Delta x^2)^2 + \frac{1}{4}(2f'(u_i) \Delta x- \frac{1}{3}f'''(u_i)\Delta x^3)^2 + \mathcal{O}(\Delta x^6)\\
    IS_2 &= \frac{13}{12}(f''(u_i)\Delta x^2)^2 + \frac{1}{4}(2f'(u_i) \Delta x- \frac{2}{3}f'''(u_i)\Delta x^3)^2 + \mathcal{O}(\Delta x^6)\\
\end{align*}
$$
Therefore, these new weights give the optimal order ($2r-1 = 5$) for $r = 3$. This WENO method is usually refered to as WENO-5, or WENO-JS.

## Time discretization
For the time discretization scheme, both versions of WENO methods are combined with high order TVD Runge-Kutta developed by Shu and Osher. With the spatial operator $L(u)$, the third order Runge-Kutta is
$$
\begin{align*}
    &u^{(1)} = u^n + \Delta t L(u^n)\\
    &u^{(2)} = \frac{3}{4} u^n + \frac{1}{4} u^{(1)} + \frac{1}{4}\Delta t L(u^{(1)})\\
    &u^{n+1} = \frac{1}{3}u^n + \frac{2}{3}u^{(2)} + \frac{2}{3}\Delta t L(u^{(2)}).
\end{align*}
$$
The fourth order Runge-Kutta scheme is given by
$$
\begin{align*}
    &u^{(1)} = u^n + \frac{1}{2} \Delta t L(u^n)\\
    &u^{(2)} = u^n + \frac{1}{2} \Delta t L(u^{(1)})\\
    &u^{(3)} = u^n + \Delta t L(u^{(2)})\\
    &u^{n+1} = \frac{1}{3}(-u^n + u^{(1)} + 2u^{(2)} + u^{(3)}) + \frac{1}{6}\Delta t L(u^{(3)}).
\end{align*}
$$

## Stability analysis
In this section we briefly discuss the stability region of fifth order WENO method using von Neumann analysis. The discrete Fourier series in space gives
$$
\begin{align*}
    u_j(t) = \sum^{N/2}_{m = -N/2} \hat{u}_m(t) e^{i \omega_m j \Delta x}\\
    u_j(t) = \hat{u}_m(t)e^{i j \theta_m}, \theta_m = \omega_m \Delta x.
\end{align*}
$$
For a linear spatial operator $L$, we also have
$$
L(u_{j-r}, \cdots, u_{j+s}) = z(\theta_m) u_j.
$$
Combining with a time discretization with step size $\Delta t$. Let $u^n_j = u_j(t^n)$, we obtain
$$
u^{n+1}_j = g(-\sigma z(\theta_m)) u^n_j.
$$
The method is stable if $|g(-\sigma z(\theta))| \leq 1$ for $\theta \in [0, 2\pi]$. Normally, we can use a more relaxed condition, $|g(-\sigma z(\theta))| \leq 1 + C\Delta t$.

Now we can linearize WENO method on advection equation.
$$
u_t + u_x = 0.
$$
The space operator is
$$
L = \hat{f}_{j+1/2} - \hat{f}_{j-1/2}.
$$
The flux reconstruction and weights are provided as example above. In the smooth regions, Taylor expansion of the smoothness indicators gives
$$
IS_k = \frac{13}{12}((\Delta x)^2 u''_j)^2 + \frac{1}{4}(2 \Delta x u'_j + c_k (\Delta x)^3 u'''_j)^2 + \mathcal{O}((\Delta x) ^6),
$$
where $c_0 = -\frac{2}{3}, c_1 = \frac{1}{3}, c_2 = -\frac{2}{3}$, $u'_j$ is the first derivative of $u(x,t)$ evaluated at $x_j$. Then, we have
$$
IS_k = (u' \Delta x)^2(1 + \mathcal{O}((\Delta x)^2)), k = 0, 1, 2,
$$
where $u'_j \neq 0$, and
$$
IS_k = \frac{13}{12}(u'' (\Delta x)^2)^2(1 + \mathcal{O}(\Delta x)^2), k = 0, 1, 2.
$$
Then we obtain the flux
$$
\hat{f}_{j+1/2} = \frac{1}{30}u_{j-2} - \frac{13}{60}u_{j-1} + \frac{47}{60}u_j + \frac{9}{20}u_{j+1} - \frac{1}{20}u_{j+2},
$$
so the corresponding linearized operator for fifth order WENO becomes
$$
L = -\frac{1}{30}u_{j-3} + \frac{1}{4}u_{j-2} - u_{j-1} + \frac{1}{3}u_j + \frac{1}{2}u_{j+1} - \frac{1}{20}u_{j+2}.
$$
Then we found $z(\theta)$ to be
$$
z(\theta) = \frac{16}{15}\sin^6{\frac{\theta}{2}} + i\left(\frac{1}{6}\sin{2\theta} - \frac{4}{3} \sin{\theta} - \frac{16}{15}\sin^5{\frac{\theta}{2}}\cos{\frac{\theta}{2}} \right).
$$
The discrete spectrum $S$ of a spatial discretization scheme is the set of eigenvalues $S = \{-z(\theta_m): \theta_m \in \{0, \Delta \theta, 2\Delta \theta, \cdots, 2\pi \}, \Delta \theta = 2\pi \Delta x\}$. We can now obtain the plot of spectrum for linearized WENO method on advection equation. The figure below shows the boundary of stability domain of forward Euler together with WENO's spectrum.

![Spectrum of Linearized WENO](images/WENOspectrum.PNG

Figure2: Spectrum of Linearized WENO.

The absolute stability condition can be expressed as
$$
(1 - \sigma Re(z_m))^2 + (\sigma Im(z_m))^2 \leq 1,
$$
which means that linearized WENO combined with Forward Euler would be stable if the above condition is satisfied. Analysis on other time discretization and WENO can be done in a similar way.



## Numerical Example
In this section we run a numerical example through PyClaw that uses fifth order WENO methods to solve one dimensional Euler equations for inviscid, compressible flow. The Woodward-Colella Blast Wave problem describes a one dimensional tube filled with ideal gas. Two infinitely thin layers were put into the tube, thus separating the gas into three sections. The middle section is the largest and has the lowest pressure. The problem is formulated as
$$
\begin{align*}
    \rho_t + (\rho u)_x &= 0\\
    (\rho u)_t + (\rho u^2 + p)_x &= 0\\
    E_t + (u(E + p))_x &= 0,
\end{align*}
$$
where $p = \rho(\gamma - 1)e$ is the pressure of the ideal gas, and $e$ is the internal energy. Since the gas is separated and both ends are sealed, we expect to see two shocks in the solution and they will approach the center due to the difference in pressure. The coefficients used in this problems are $\gamma = 1.4, \rho_l = 1, \rho_r = \frac{1}{8}, p_l = 1, p_r =0.1$ and we want to know the solution between $t = 0$ and $t = 0.04$. The sharpClawSolver from PyClaw uses high order WENO methods. Density and energy at a given time step are shown below.

![t0](images/t=0.PNG)

![t3](images/t3.PNG)

![t8](images/t8.PNG)

![t10](images/t10.PNG)


## Further Improvements
The scheme proposed by Jiang and Shu lays the foundation for future modifications for this class of weighted essentially non-oscillatory methods. Some refinements are done to the scheme in order to achieve the same order of accuracy near critical points as well. One of the improvements is Mapped WENO scheme that increases the accuracy of the weights. Some mapping functions are given by
$$
g_k(w) = \frac{w(\bar{w}_k + \bar{w}_k^2) - 3\bar{w}_kw + w^2}{\bar{w}_k^2 + \bar{w}(1 - 2\bar{w}_k)}, \bar{w}_k \in (0, 1)
$$
for $k = 0, 1, 2$. The functions are monotonically increasing and have $g_k(0) = 0, g_k(1) = 1, g_k(\bar{w}_k) = \bar{w}_k, g'_k(\bar{w}_k) = 0$, and $g''_k(\bar{w}_k) = 0$. Then the approximation of the weights can be written as
$$
\alpha^*_k = g_k(w_k),
$$
where $w_k$'s are the original weights and Jiang and Shu's scheme. This modified scheme is normally refered to as WENO-M.

Besides adding a mapping function, another smoothness indicator of higher order is devised and the new scheme (WENO-Z) requires less computing power than WENO-M. The whole 5-points stencil is used and the new indicators, $\tau_5$, is defined as the absolute difference between two classical smoothness indicators $\beta_k$ at $x_i$,
$$
\tau_5 = |\beta_0 - \beta_2|,
$$
and the rest of the weights change accordingly.
$$
w_k^Z = \frac{\alpha_k^Z}{\sum^2_{l=0} \alpha^Z_l}, \alpha^Z_k = \frac{C_k}{\beta^Z_k} = C_k(1 + \frac{\tau_5}{\epsilon + \beta_k}).
$$

Other research has been conducted in the relevant field. For example,  the accuracy and performance are studied for WENO methods in two dimensional Cartesian meshes. WENO methods are also compared with Runge-Kutta discontinuous Galerkin methods. In the future, I would also like to read about how WENO scheme is applied to solve Fokker-Planck-Kolmogorov equations as nonlinear one dimensional system, which relates my project for last semester with this one.





## Bibliography
Borges, R., Carmona, M., Costa, B., Don, W.S.
An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws. Journal of Computational Physics 207 (2008), pp. 3191-3221.
https://www.sciencedirect.com/science/article/pii/S0021999107005232


Clawpack Development Team (2019),
Clawpack Version 5.6.0.
http://www.clawpack.org


Colella, P., Woodward, P.R.
The Piecewise Parabolic Method for gas-dynamical simulations. Journal of Computational Physics 54 (1984), pp. 174-201.
https://doi.org/10.1016/0021-9991(84)90143-8


Henrick, A.K., Aslam, T.D., Powers, J.M
Mapped weighted essentially non oscillatory scheme: Achieving optimal order near critical points. Journal of Computational Physics 207 (2005), pp.542-567
https://www3.nd.edu/~powers/paper.list/jcp.weno5m.05.pdf


Jiang, G.S., Shu, C.W.
Efficient Implementation of Weighted ENO Schemes. Journal of Computational Physics 206 (1996), pp. 202-228.
https://www.jstor.org


Ketcheson, D.T., Mandli, K.T., Ahmedia, A.J., Alghamdi, A., Quezada de Lunea, M., Parsani, M., Knepley, M.G., Emmett, M.
PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation Problems. SIAM Journal on Scientific Computing 34 (2012), pp. C210-C231.


Ketcheson, D.T., Parsani, M., LeVeque, R.J.,
High-order Wave Propagation Algorithms for Hyperbolic Systems. SIAM Journal on Scientific Computing 35 (2013), pp. A351-A377.



Liu, X.D., Osher, S., Chan, T.
Weighted Essentially Non-oscillatory Schemes. Journal of Computational Physics 115 (1994), pp.200-212.


Mandli, K.T., Ahmadia, A.J., Berger, M.J., Calhoun, D.A., George, D.L.,
Hadjimichael, Y., Ketcheson, D.I., Lemoine, G.I., LeVeque, R.J.,
Clawpack: building an open source ecosystem for solving hyperbolic PDEs. PeerJ Computer Science.
doi:10.7717/peerj-cs.68


Motamed, M., Macdonald, C.B., Ruuth, S.J.
On the Linear Stability of the Fifth-Order WENO Discretization. Simon Fraser University.
http://people.math.sfu.ca/~sruuth/weno5\stability.pdf


Shu, C.W.
Essentially Non-oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic Conservation Laws. NASA, ICASE Report, No. 97-65.
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19980007543.pdf


Wang, W.J., Feng, J.H., Xu, W.
The Numerical Solution of the TVD Runge-Kutta and WENO Scheme to the FPK Equations to Nonlinear System of One-Dimension. Applied and Computational Mathematics 5 (2016), pp. 160-164
http://article.sciencepublishinggroup.com/html


Zhang, R., Zhang, M.P., Shu, C.W.
On the order of accuracy and numerical performance of two classes of ginite volume WENO schemes. University of Science and Technology of China, 2009.
https://www.brown.edu/research/projects/scientific-computing


Zhang, Y.T., Shu, C.W.
ENO and WENO Schemes. Handbook of Numerical Analysis, Vol. 17. Ch. 5.
https://www3.nd.edu/~yzhang10/WENO\_ENO.pdf


Zhou, T., Li, Y.F., Shu, C.W.
Numerical Comparison of WENO Finite Volume and Runge-Kutta Discontinuous Galerkin Methods. Journal of Scientific Computing, 16 (2001).
https://link.springer.com/content/pdf
