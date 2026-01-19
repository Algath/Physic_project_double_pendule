![pendulum meme](/assets/double-pendulum-predictability-meme.png)

# Double pendulum — Theory and modeling
This document gives a focused, theory‑heavy presentation of the planar double pendulum used in the project, including kinematics, energy expressions, a Lagrangian derivation of the equations of motion, common approximations (small‑angle), damping options, numerical considerations, parameter estimation strategy, and a short comparison with ballistic motion.

**Goal:** provide clear, reproducible equations that match the implementation in `pendule.jl` and background needed to calibrate and extend the numerical model.

---

## What is a double pendulum?

A double pendulum is a simple mechanical system consisting of two pendula attached end‑to‑end: a pendulum (mass m1 on a rod of length L1) with a second pendulum (mass m2 on rod L2) attached to its end. It is also commonly referred to as a "compound pendulum" or—in contexts emphasizing its dynamical behavior—the "chaotic pendulum" or "planar double pendulum".

Key properties:
- Deterministic but often chaotic: the motion is fully determined by Newton's laws (or equivalently by a Lagrangian formulation), yet for many initial conditions the system shows sensitive dependence on initial conditions — small changes in initial angles or velocities produce rapidly diverging trajectories.
- Low dimensional but rich behavior: the minimal planar double pendulum has four state variables (two angles and two angular velocities) and can produce periodic, quasi-periodic, or chaotic trajectories depending on energy and initial conditions.
- Common uses: demonstration of chaos and nonlinear dynamics, validation benchmarks for numerical integrators, and pedagogical examples in mechanics courses. It also appears in robotics and control contexts when multiple hinged links are analyzed.

Practical consequence for this project: because of the sensitivity, parameter estimation and short-term prediction must be handled with care (short fitting windows, careful noise filtering, and uncertainty quantification via ensembles).

## 1. System definition and kinematics

Geometry and notation:
- Two rigid, massless rods of lengths $L_1, L_2$.
- Two point masses $m_1$ (upper) and $m_2$ (lower).
- Angles $\theta_1(t),\theta_2(t)$ measured from the downward vertical (so $\theta=0$ is the rest vertical configuration).

Cartesian positions (useful to relate tracking data to angles):
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

These expressions are useful for computing kinetic energy and for transforming tracking coordinates (pixels) into the angular variables used by the model.

## 2. Energies and Lagrangian

- Kinetic energy:
$$
T=\tfrac12 m_1\lVert\dot{\mathbf{r}}_1\rVert^2+\tfrac12 m_2\lVert\dot{\mathbf{r}}_2\rVert^2.
$$
Expanding the squared norms gives explicit terms in $\dot\theta_1^2$, $\dot\theta_2^2$ and the coupling $\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2)$.

- Potential energy (gravitational, reference at origin):
$$
V=m_1 g y_1 + m_2 g y_2 = (m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2.
$$

- Lagrangian: $\mathcal{L} = T - V$.

Derive the Euler–Lagrange equations for generalized coordinates $\theta_1,\theta_2$:
$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_i}\right)-\frac{\partial\mathcal{L}}{\partial\theta_i}=Q_i,\qquad i=1,2
$$
where $Q_i$ are generalized non-conservative torques (zero for the conservative model).

### Hamilton's principle and the Euler–Lagrange derivation

The equations above follow from Hamilton's principle (principle of stationary action): the actual trajectory $q(t)$ of a mechanical system between fixed times $t_1$ and $t_2$ makes the action
$$
S[q]=\int_{t_1}^{t_2}\mathcal{L}(q,\dot q,t)\,dt
$$
stationary under all variations $\delta q(t)$ that vanish at the endpoints. The variation of the action is
$$
\delta S=\int_{t_1}^{t_2}\left(\frac{\partial\mathcal{L}}{\partial q}\delta q + \frac{\partial\mathcal{L}}{\partial\dot q}\delta\dot q\right)dt.
$$
Integrating the second term by parts and using $\delta q(t_1)=\delta q(t_2)=0$ gives
$$
\delta S=\int_{t_1}^{t_2}\left(\frac{\partial\mathcal{L}}{\partial q}-\frac{d}{dt}\frac{\partial\mathcal{L}}{\partial\dot q}\right)\delta q\,dt.
$$
Since $\delta q$ is arbitrary, the integrand must vanish, yielding the Euler–Lagrange equation for each generalized coordinate:
$$
\frac{\partial\mathcal{L}}{\partial q}-\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot q}\right)=0.
$$
If non-conservative generalized forces $Q_i$ act on the system, they appear on the right-hand side and recover the form used above.

Applied to $q=\theta_1,\theta_2$, this procedure (compute $T$ and $V$, form $\mathcal{L}=T-V$, compute the partial derivatives and time derivatives) yields the explicit angular-acceleration formulas given in §3. The Lagrangian approach is convenient because it automatically accounts for kinematic constraints and uses generalized coordinates directly.

## 3. Explicit equations of motion (final form)

Let $\delta=\theta_2-\theta_1$, and define the convenient denominator
$$\mathrm{den}=m_1+m_2\sin^2\delta.
$$
Then the angular accelerations can be written (algebraically equivalent to the implementation in `pendule.jl`):

$$
\begin{aligned}
\ddot\theta_1 &= \frac{m_2\sin\delta\big(L_1\dot\theta_1^2\cos\delta + L_2\dot\theta_2^2\big) + m_2 g \sin\theta_2\cos\delta - (m_1 + m_2) g \sin\theta_1}{L_1\,\mathrm{den}},\\[6pt]
\ddot\theta_2 &= \frac{-m_2 L_2\dot\theta_2^2\sin\delta\cos\delta - (m_1 + m_2) L_1\dot\theta_1^2\sin\delta + (m_1 + m_2) g \sin\theta_1\cos\delta - (m_1 + m_2) g \sin\theta_2}{L_2\,\mathrm{den}}.
\end{aligned}
$$

These expressions are non-linear and couple angles and angular velocities; they are the ones implemented in `equations_double_pendulum!` in `pendule.jl`.

### Comments on the algebra
- The coupling appears via $\sin\delta$ and $\cos\delta$ factors and via products of squared angular velocities.
- The denominator arises from eliminating internal constraint reaction forces and depends on geometry and masses.

## 4. Common approximations and limits

- Small‑angle linearization (for small $\theta_i$): replace $\sin\theta\approx\theta$, $\cos\theta\approx1$, keep only linear terms. This yields a linear 2×2 system with constant coefficients; useful to analyze normal modes and natural frequencies but invalid for large motions or chaotic regimes.

- Single pendulum limit: set $m_2\to0$ (or $L_2\to0$) and recover the simple pendulum equation
$\ddot\theta + \frac{g}{L}\sin\theta=0.$

- Energy conservation: in the conservative model $E=T+V$ remains constant. Observed energy decay in data indicates non-conservative forces (damping, friction), motivating the addition of damping terms and parameter estimation.

## 5. Damping and non-conservative forces

Two common damping models:

- Viscous (torque proportional to angular velocity): add generalized torques $Q_i=-c_i\dot\theta_i$. In practice one can subtract a term from $\ddot\theta_i$ of the form $c_i\dot\theta_i/L_i$ (as an approximation consistent with how torques map to angular accelerations in the chosen equations).

- Aerodynamic drag (more complex): a drag torque depending on $\dot\theta$ magnitude or on fluid dynamic coefficients; often non-linear in velocity and possibly different for each mass.

Fitting $c_i$ (or an equivalent loss model) can be done by matching the observed energy decay over a calibration window, or by optimizing simulated positions/angles against tracked data over a later time window.

## 6. Numerical integration and implementation notes

- The system is a set of 4 first‑order ODEs for $[\theta_1,\omega_1,\theta_2,\omega_2]$ with right-hand side given by the above formulas. Use a robust non-stiff integrator (e.g. Tsit5) from `DifferentialEquations.jl` as in `pendule.jl`.
- Save the solution at the same sampling instants as the tracker (`saveat`) to directly compare simulated angles with observed angles.
- Chaotic divergence: even with a perfect model, trajectories can diverge from measurements due to sensitivity to initial conditions. Use short windows for parameter fitting and evaluate prediction skill probabilistically (ensembles) rather than deterministically for long horizons.

## 7. Parameter estimation strategy (practical guidance)

- Step 1 — Geometry & initial angles: estimate lengths $L_i$ from calibration or image measurements (the project already uses measured lengths). Convert tracked pixel coordinates to meters before converting to angles.
- Step 2 — Initial angular velocities: estimate $\dot\theta$ from finite differences on the first frames (central difference with noise filtering helps).
- Step 3 — Masses & initial speeds: optimize masses (and possibly small corrections to initial angular velocities) on a short early window where damping is minimal (the project uses 0–0.5s or 0–1s windows to avoid chaotic mismatch).
- Step 4 — Damping: estimate damping coefficients by fitting on later windows where energy decay is visible (or by matching energy loss between frames). Try simple viscous damping first.
- Use objective functions based on squared angle errors or Cartesian position errors (convert simulated angles back to positions for a direct geometric comparison).

## 8. Chaos, validation and uncertainty

- The double pendulum exhibits sensitive dependence on initial conditions: validate results on multiple short windows and report skill metrics (mean error, RMSE, or distribution over ensembles with perturbed initial conditions).
- When predicting beyond the observed data window, provide confidence bands (ensemble forecast) and avoid over‑interpreting single-trajectory predictions beyond the Lyapunov time scale.

## 9. Comparison with ballistic equations

- Ballistics (projectile motion) typically models a single rigid body's center-of-mass translation. The canonical (no-drag) equations are
$$x''(t)=0,\qquad y''(t)=-g,$$
which are decoupled and admit closed-form solutions. Adding drag couples velocity components through $v=\|\mathbf v\|$ and requires numerical integration, but the state dimensionality remains small and the equations are simpler in structure than for pendulums.
- Pendulums use angular generalized coordinates and produce inherently nonlinear, geometrically coupled ODEs with terms in $\sin(\cdot)$ and products of angular velocities. The double pendulum's phase space is higher-dimensional and can be chaotic; ballistics is typically non‑chaotic.

## 10. Practical tips and extensions

- Filter tracking noise before differentiating to estimate velocities.
- When optimizing masses, keep realistic bounds (masses positive and in a plausible range) to avoid pathological fits.
- If long-term prediction is required, consider estimating a stochastic forcing term or using ensemble forecasts to represent model error and initial condition uncertainty.
- For visualization, compare both angle time series and Cartesian traces of mass positions (convert angles using the kinematics in §1).

---

## Appendix: symbolic Euler–Lagrange derivation

This appendix expands the kinetic and potential energies, computes the necessary partial derivatives, and shows the Euler–Lagrange manipulation that leads to the explicit angular-acceleration formulas in §3.

1) Kinematics recap (for convenience):
$$
\mathbf{r}_1=(L_1\sin\theta_1,\;L_1\cos\theta_1),\qquad
\mathbf{r}_2=(L_1\sin\theta_1+L_2\sin\theta_2,\;L_1\cos\theta_1+L_2\cos\theta_2).
$$

2) Velocities squared (useful identities):
$$
\lVert\dot{\mathbf{r}}_1\rVert^2=L_1^2\dot\theta_1^2,
\lVert\dot{\mathbf{r}}_2\rVert^2=L_1^2\dot\theta_1^2+L_2^2\dot\theta_2^2+2L_1L_2\dot\theta_1\dot\theta_2\cos(\theta_1-\theta_2).
$$

3) Kinetic and potential energy (expanded):
$$
T=\tfrac12 m_1 L_1^2\dot\theta_1^2 + \tfrac12 m_2\big(L_1^2\dot\theta_1^2 + L_2^2\dot\theta_2^2 + 2L_1L_2\dot\theta_1\dot\theta_2\cos\delta\big),
V=(m_1+m_2)gL_1\cos\theta_1 + m_2 g L_2\cos\theta_2,
$$
with $\delta=\theta_2-\theta_1$.

Combine terms in $T$ to write the Lagrangian $\mathcal{L}=T-V$ as:
$$
\mathcal{L}=\tfrac12(m_1+m_2)L_1^2\dot\theta_1^2 + \tfrac12 m_2 L_2^2\dot\theta_2^2 + m_2 L_1 L_2\dot\theta_1\dot\theta_2\cos\delta - V.
$$

4) Partial derivatives with respect to velocities:
$$
\frac{\partial\mathcal{L}}{\partial\dot\theta_1}=(m_1+m_2)L_1^2\dot\theta_1 + m_2 L_1 L_2\dot\theta_2\cos\delta,

\frac{\partial\mathcal{L}}{\partial\dot\theta_2}=m_2 L_2^2\dot\theta_2 + m_2 L_1 L_2\dot\theta_1\cos\delta.
$$

5) Time derivatives of the above (use product rule and $\dot\delta=\dot\theta_2-\dot\theta_1$):
$$
\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_1}\right)=(m_1+m_2)L_1^2\ddot\theta_1 + m_2 L_1 L_2\big(\ddot\theta_2\cos\delta - \dot\theta_2\sin\delta\,(\dot\theta_2-\dot\theta_1)\big),

\frac{d}{dt}\left(\frac{\partial\mathcal{L}}{\partial\dot\theta_2}\right)=m_2 L_2^2\ddot\theta_2 + m_2 L_1 L_2\big(\ddot\theta_1\cos\delta - \dot\theta_1\sin\delta\,(\dot\theta_1-\dot\theta_2)\big).
$$

6) Partial derivatives with respect to angles (note signs from $V$ and dependence of $\cos\delta$ on the angles):
$$
\frac{\partial\mathcal{L}}{\partial\theta_1}=m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + (m_1+m_2)gL_1\sin\theta_1,

\frac{\partial\mathcal{L}}{\partial\theta_2}=-m_2 L_1 L_2\dot\theta_1\dot\theta_2\sin\delta + m_2 g L_2\sin\theta_2.
$$

7) Euler–Lagrange equations $\frac{d}{dt}(\partial_{\dot\theta_i}\mathcal{L})-\partial_{\theta_i}\mathcal{L}=0$ give two scalar equations. Substituting the expressions from steps 5 and 6 and rearranging terms (grouping $\ddot\theta_1$, $\ddot\theta_2$, gravitational and nonlinear velocity terms) yields the explicit forms shown in §3.

After algebraic simplification the two equations can be written in the compact form used in the implementation; introducing the denominator
$$\mathrm{den}=m_1 + m_2\sin^2\delta$$
one obtains the pair (reproduced from §3):
$$
\begin{aligned}
\ddot\theta_1 &= \frac{m_2\sin\delta\big(L_1\dot\theta_1^2\cos\delta + L_2\dot\theta_2^2\big) + m_2 g \sin\theta_2\cos\delta - (m_1 + m_2) g \sin\theta_1}{L_1\,\mathrm{den}},\\
\ddot\theta_2 &= \frac{-m_2 L_2\dot\theta_2^2\sin\delta\cos\delta - (m_1 + m_2) L_1\dot\theta_1^2\sin\delta + (m_1 + m_2) g \sin\theta_1\cos\delta - (m_1 + m_2) g \sin\theta_2}{L_2\,\mathrm{den}}.
\end{aligned}
$$

Remarks:
- The intermediate steps involve combining terms that contain $\sin\delta$ or $\cos\delta$ and factoring common coefficients; the appendix here shows the key derivatives so the algebraic rearrangement is straightforward to verify symbolically (or with a CAS).
- If desired I can produce a fully expanded symbolic derivation (explicit polynomial in sines and cosines) using a CAS (e.g. SymPy or SymEngine) and include it as a supplementary file.
