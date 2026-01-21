![pendulum meme](/assets/double-pendulum-predictability-meme.png)

# Double pendulum — Theory and modelling

This document provides a theory-focused presentation of the planar double pendulum used in the project. It covers kinematics, energy expressions, a Lagrangian derivation of the equations of motion, common approximations (small-angle), damping options, numerical considerations, parameter estimation strategy, and a short comparison with projectile motion.

**Goal:** give clear, reproducible equations that match the implementation in `pendule.jl` and provide the background needed to calibrate and extend the numerical model.

---

## What is a double pendulum?

A double pendulum is a mechanical system of two pendula attached end-to-end: an upper pendulum (point mass $m_1$ on a rigid, massless rod of length $L_1$) with a second pendulum (mass $m_2$ on rod length $L_2$) attached to its end. It is often called a "compound pendulum" or — in dynamical contexts — a "chaotic pendulum" or "planar double pendulum".

Key properties:

- **Deterministic but often chaotic**: the motion follows Newton's laws (or an equivalent Lagrangian formulation), yet many initial conditions display sensitive dependence on initial conditions — small changes in angles or velocities may produce rapidly diverging trajectories.
  -- **Low-dimensional, rich dynamics**: the planar double pendulum has four state variables (two angles and two angular velocities) and can show periodic, quasi-periodic, or chaotic motion depending on energy and initial conditions.
- **Common uses**: demonstration of chaos and nonlinear dynamics, benchmarks for numerical integrators, and pedagogical examples in mechanics. It also appears in robotics and control for multi-link analyses.

Practical consequence for this project: because of this sensitivity, parameter estimation and short-term prediction must be treated carefully (short fitting windows, noise filtering, and uncertainty quantification via ensembles).

## 1. System definition and kinematics

Geometry and notation:

- Two rigid, massless rods of lengths $L_1, L_2$.
- Two point masses $m_1$ (upper) and $m_2$ (lower).
- Angles $\theta_1(t),\theta_2(t)$ measured from the downward vertical (so $\theta=0$ is the resting vertical configuration).

Cartesian positions (useful to convert tracking data to angles):

$$
\mathbf{r}_1=(x_1,y_1)=(L_1\sin\theta_1,\;L_1\cos\theta_1),
\qquad
\mathbf{r}_2=(x_2,y_2)=(L_1\sin\theta_1+L_2\sin\theta_2,\;L_1\cos\theta_1+L_2\cos\theta_2).
$$

Velocities (time derivatives):

$$
\dot{\mathbf{r}}_1=(L_1\cos\theta_1\,\dot\theta_1,\;-L_1\sin\theta_1\,\dot\theta_1),\qquad
\dot{\mathbf{r}}_2=(L_1\cos\theta_1\,\dot\theta_1+L_2\cos\theta_2\,\dot\theta_2,\;-L_1\sin\theta_1\,\dot\theta_1-L_2\sin\theta_2\,\dot\theta_2).
$$

These relations are used to compute kinetic energy and to transform tracked coordinates (pixels) into the angular variables used by the model.

## 2. Energies and Lagrangian

- **Kinetic energy:**

$$
T=\tfrac12 m_1\lVert\dot{\mathbf{r}}_1\rVert^2+\tfrac12 m_2\lVert\dot{\mathbf{r}}_2\rVert^2.
$$

Expanding the squared norms yields explicit terms in $\dot\theta_1^2$, $\dot\theta_2^2$ and the coupling $\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2)$.

- **Potential energy** (gravitational, reference at origin):

$$
V=m_1 g y_1 + m_2 g y_2 = (m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2.
$$

- **Lagrangian:** $\mathcal{L}=T-V$.

Derive the Euler–Lagrange equations for generalized coordinates $\theta_1,\theta_2$:

$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_i}\right)-\frac{\partial\mathcal{L}}{\partial\theta_i}=Q_i,\qquad i=1,2,
$$

where $Q_i$ are generalized non-conservative torques (zero for the conservative model).

### Hamilton's principle and Euler–Lagrange

These equations follow from Hamilton's principle (stationary action): the actual trajectory $q(t)$ between fixed times $t_1$ and $t_2$ extremizes the action

$$
S[q]=\int_{t_1}^{t_2}\mathcal{L}(q,\dot q,t)\,dt.
$$

Taking variations and integrating by parts yields the Euler–Lagrange equation

$$
\frac{\partial\mathcal{L}}{\partial q}-\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot q}\right)=0.
$$

If non-conservative generalized forces $Q_i$ act, they appear on the right-hand side.

Applied to $\theta_1,\theta_2$ this procedure (compute $T$ and $V$, form $\mathcal{L}=T-V$, compute derivatives) leads to the explicit angular-acceleration formulas given next. The Lagrangian approach is convenient because it automatically enforces kinematic constraints and uses generalized coordinates directly.

## 3. Explicit equations of motion (final form)

Let $\delta=\theta_2-\theta_1$ and define the practical denominator

$$
\mathrm{den}=m_1+m_2\sin^2\delta.
$$

Then the angular accelerations can be written (algebraically equivalent to the implementation in `pendule.jl`):

$$
\begin{aligned}
\ddot\theta_1 &= \frac{m_2\sin\delta\big(L_1\dot\theta_1^2\cos\delta + L_2\dot\theta_2^2\big) + m_2 g \sin\theta_2\cos\delta - (m_1 + m_2) g \sin\theta_1}{L_1\,\mathrm{den}},\\[6pt]
\ddot\theta_2 &= \frac{-m_2 L_2\dot\theta_2^2\sin\delta\cos\delta - (m_1 + m_2) L_1\dot\theta_1^2\sin\delta + (m_1 + m_2) g \sin\theta_1\cos\delta - (m_1 + m_2) g \sin\theta_2}{L_2\,\mathrm{den}}.
\end{aligned}
$$

These non-linear expressions couple angles and angular velocities and correspond to the formulas used in `equations_double_pendulum!` in `pendule.jl`.

### Algebraic comments

- The coupling appears through $\sin\delta$ and $\cos\delta$ factors and via squared angular-velocity products.
- The denominator arises when eliminating internal reaction forces and depends on masses and geometry.

## 4. Common approximations and limits

- **Small-angle linearization:** for small $\theta_i$ use $\sin\theta\approx\theta$, $\cos\theta\approx1$ and keep linear terms. This yields a linear 2×2 system useful for analyzing normal modes and natural frequencies but invalid for large amplitudes or chaotic regimes.
- **Single pendulum limit:** letting $m_2\to0$ (or $L_2\to0$) recovers the simple pendulum equation $\ddot\theta + \frac{g}{L}\sin\theta=0$.
- **Energy conservation:** in the conservative model $E=T+V$ is constant. Observed energy decay indicates non-conservative effects (damping, friction), motivating the inclusion of loss models and parameter estimation.

## 5. Damping and non-conservative forces

Two common damping models:

- **Viscous damping:** generalized torques $Q_i=-c_i\dot\theta_i$. In practice one can subtract a term proportional to $c_i\dot\theta_i/L_i$ from the expressions for $\ddot\theta_i$ (an approximation consistent with torque→angular-acceleration mapping used here).
- **Aerodynamic drag:** more complex torque laws depending on $\|\dot\theta\|$ or fluid-dynamic coefficients; typically nonlinear in velocity and possibly different for each mass.

Estimate $c_i$ (or an equivalent loss parameter) by matching the observed energy decay over a calibration window, or by optimizing simulated positions/angles against tracked data on a later window. Start with simple viscous damping.

## 6. Numerical integration and implementation notes

- The system is a set of four first-order ODEs for $[\theta_1,\omega_1,\theta_2,\omega_2]$ with right-hand side given by the formulas above. Use a robust non-stiff integrator (e.g. `Tsit5`) from `DifferentialEquations.jl` as in `pendule.jl`.
- Save the solution at the same sampling instants as the tracker (`saveat`) to directly compare simulated angles with observed angles.
- Chaotic divergence: even with a correct model, trajectories may diverge from measurements due to sensitivity to initial conditions. Use short calibration windows and evaluate prediction skill probabilistically (ensembles) rather than relying on a single long deterministic trajectory.

### Runge–Kutta methods — why they are appropriate

Runge–Kutta (RK) methods evaluate the right-hand side $f(t,y)$ multiple times per step to construct a higher-order update. Compared to explicit Euler, RK methods (e.g. RK4, embedded pairs like Tsit5) give much higher accuracy per step. Key points:

- **Order and accuracy:** a scheme of order $p$ has global error $O(h^p)$; higher order strongly reduces error for the same step size.
- **Adaptive step control:** embedded RK pairs provide stepwise error estimates and adapt the step size to meet tolerances.
- **Stability:** explicit RK methods are efficient for non-stiff problems (like this one); for stiff systems implicit methods are preferable.
- **Cost vs accuracy:** higher-order RK methods evaluate $f$ more times per step but permit much larger steps for a given accuracy, often reducing total cost.

For the double pendulum, `Tsit5` (an explicit 5/4 order pair) is a good default: accurate, efficient, and with built-in error control. For long-term energy preservation consider symplectic integrators, but for parameter estimation and short-term prediction RK adaptive methods are typically preferable.

## 7. Parameter estimation strategy (practical guidance)

1. **Geometry & initial angles:** estimate $L_i$ from calibration or image measurements (the project already uses measured lengths). Convert tracked pixel coordinates to metres before computing angles.
2. **Initial angular velocities:** estimate $\dot\theta$ by finite differences on early frames (central differences with noise filtering recommended).
3. **Masses & initial speeds:** optimize masses (and optionally small corrections to initial angular velocities) on a short early window where damping is minimal (typical windows: 0–0.5 s or 0–1 s to limit chaotic mismatch).
4. **Damping:** estimate damping coefficients from later windows where energy decay is visible, or by matching energy loss between frames. Begin with viscous damping.

Use objective functions based on squared angle errors or Cartesian position errors (convert simulated angles to positions to compare geometrically).

## 8. State-space formulation

Define the state vector

$$
y=[\theta_1,\omega_1,\theta_2,\omega_2]^\top
$$

and the vector field

$$
F(t,y)=\begin{bmatrix}\omega_1\\[4pt] f_1(\theta_1,\theta_2,\omega_1,\omega_2)\\[4pt] \omega_2\\[4pt] f_2(\theta_1,\theta_2,\omega_1,\omega_2)\end{bmatrix}
$$

where $f_1,f_2$ return the angular accelerations given in §3. At each integrator call the solver evaluates $F$ and advances the state; embedded RK methods evaluate $F$ multiple times per step to obtain higher-order accuracy.

## 9. Statistics and error metrics

Useful metrics to evaluate fit quality and divergence between measurements and simulation:

- **Coefficient of determination $R^2$:**

$$
R^2 = 1 - \frac{\sum_i (y_i-\hat y_i)^2}{\sum_i (y_i-\bar y)^2},
$$

where $\bar y$ is the mean of measurements.

- **RMSE (root mean squared error):**

$$
\mathrm{RMSE}=\sqrt{\frac{1}{N}\sum_{i=1}^N (y_i-\hat y_i)^2}.
$$

- **MAE (mean absolute error):**

$$
\mathrm{MAE}=\frac{1}{N}\sum_{i=1}^N |y_i-\hat y_i|.
$$

- **Energy statistics:** for observed energy series $E_i$ (or residuals $E_i-\hat E_i$) use the sample standard deviation

$$
\sigma_E=\sqrt{\frac{1}{N-1}\sum_{i=1}^N (E_i-\overline{E})^2}
$$

to quantify dispersion; measure mean energy decay rate by linear regression on $E(t)$.

- **Lyapunov exponent (practical estimate):** estimate the maximal exponent $\lambda$ by observing the exponential divergence of nearby trajectories: if $\delta(t)\approx C e^{\lambda t}$ then

$$
\lambda \approx \frac{1}{\Delta t}\left\langle\ln\frac{\delta(t+\Delta t)}{\delta(t)}\right\rangle.
$$

Practical comparisons:

- Compute RMSE/MAE on angles or on Cartesian positions (convert simulated angles to positions).
- Track energy difference $\Delta E(t)=E_{obs}(t)-E_{sim}(t)$ and report mean, standard deviation $\sigma_{\Delta E}$ and mean slope (energy loss rate).
- Measure trajectory distance $d(t)$ (e.g. Euclidean norm of mass positions) and plot $\ln d(t)$: a positive slope indicates exponential divergence (related to the Lyapunov exponent).
- Use short sliding windows and ensembles of perturbed initial conditions to obtain distributions of metrics (median, quartiles) rather than single scalars, due to sensitivity.

Normalisations (e.g. RMSE divided by data standard deviation) help compare series with different amplitudes.

## 10. Chaos, validation and uncertainty

- The double pendulum shows sensitivity to initial conditions: validate on multiple short windows and report ensemble statistics (error distributions) rather than single long-horizon fits.
- For predictions beyond the calibration window, provide confidence bands (ensemble forecasts) and avoid over-interpreting a single trajectory past the Lyapunov time.

## 11. Comparison with ballistic equations

- Projectile motion models a rigid body's center-of-mass translation. Without drag the canonical equations are

$$
x''(t)=0,\qquad y''(t)=-g,
$$

which are decoupled and admit analytic solutions. Adding drag couples components through $v=\|\mathbf v\|$ and requires numeric integration, but the structure remains simpler than pendulum dynamics.

- Pendulums use angular generalized coordinates and produce nonlinear, geometrically coupled ODEs with $\sin(\cdot)$ terms and products of angular velocities. The double pendulum's phase space can be chaotic; projectile motion is typically non-chaotic.

## 12. Practical tips and extensions

- Filter tracking noise before differentiating to estimate velocities.
- When optimizing masses, enforce realistic bounds (positive masses in plausible ranges) to avoid pathological fits.
- For long-term prediction consider adding stochastic forcing or using ensemble forecasts to capture model error and initial-condition uncertainty.
- Visualise both angle time series and Cartesian traces of the masses (use the kinematics in §1 to convert angles to positions).

---

## Appendix: symbolic Euler–Lagrange derivation

This appendix expands the kinetic and potential energies, computes the partial derivatives needed, and shows the Euler–Lagrange manipulations that lead to the explicit formulas in §3.

1) Kinematics recap:

$$
\mathbf{r}_1=(L_1\sin\theta_1,\;L_1\cos\theta_1),\qquad
\mathbf{r}_2=(L_1\sin\theta_1+L_2\sin\theta_2,\;L_1\cos\theta_1+L_2\cos\theta_2).
$$

2) Velocities squared (useful identities):

$$
\lVert\dot{\mathbf{r}}_1\rVert^2=L_1^2\dot\theta_1^2,\qquad
\lVert\dot{\mathbf{r}}_2\rVert^2=L_1^2\dot\theta_1^2+L_2^2\dot\theta_2^2+2L_1L_2\dot\theta_1\dot\theta_2\cos\delta.
$$

3) Kinetic and potential energy (expanded):

$$
T=\tfrac12 m_1 L_1^2\dot\theta_1^2 + \tfrac12 m_2\big(L_1^2\dot\theta_1^2 + L_2^2\dot\theta_2^2 + 2L_1L_2\dot\theta_1\dot\theta_2\cos\delta\big),
\qquad
V=(m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2,
$$

with $\delta=\theta_2-\theta_1$.

Combine terms to write the Lagrangian:

$$
\mathcal{L}=\tfrac12(m_1+m_2)L_1^2\dot\theta_1^2 + \tfrac12 m_2 L_2^2\dot\theta_2^2 + m_2 L_1 L_2\dot\theta_1\dot\theta_2\cos\delta - V.
$$

4) Partial derivatives with respect to velocities:

$$
\frac{\partial\mathcal{L}}{\partial\dot\theta_1}=(m_1+m_2)L_1^2\dot\theta_1 + m_2 L_1 L_2\dot\theta_2\cos\delta,\qquad
\frac{\partial\mathcal{L}}{\partial\dot\theta_2}=m_2 L_2^2\dot\theta_2 + m_2 L_1 L_2\dot\theta_1\cos\delta.
$$

5) Time derivatives (product rule, $\dot\delta=\dot\theta_2-\dot\theta_1$):

$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_1}\right)=(m_1+m_2)L_1^2\ddot\theta_1 + m_2 L_1 L_2\big(\ddot\theta_2\cos\delta - \dot\theta_2\sin\delta\,(\dot\theta_2-\dot\theta_1)\big),
$$

$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_2}\right)=m_2 L_2^2\ddot\theta_2 + m_2 L_1 L_2\big(\ddot\theta_1\cos\delta - \dot\theta_1\sin\delta\,(\dot\theta_1-\dot\theta_2)\big).
$$

6) Partial derivatives with respect to angles (note signs from $V$ and dependence of $\cos\delta$):

$$
\frac{\partial\mathcal{L}}{\partial\theta_1}=m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + (m_1+m_2)gL_1\sin\theta_1,
\qquad
\frac{\partial\mathcal{L}}{\partial\theta_2}=-m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + m_2 g L_2\sin\theta_2.
$$

7) Euler–Lagrange equations $\frac{d}{dt}(\partial_{\dot\theta_i}\mathcal{L})-\partial_{\theta_i}\mathcal{L}=0$ yield two scalar equations. Substituting the expressions above and rearranging terms (collecting $\ddot\theta_1,\ddot\theta_2$, gravitational and velocity nonlinear terms) gives the explicit forms in §3.

After algebraic simplification the pair can be written using the denominator

$$
\mathrm{den}=m_1 + m_2\sin^2\delta
$$

which yields the formulas shown in §3.

Remarks:

- The algebraic steps involve grouping terms containing $\sin\delta$ and $\cos\delta$ and factoring common coefficients; the derivatives here allow direct algebraic verification (or symbolic checking with a CAS).
