# Code Walkthrough: `01_reproduce_fig3.py`

This document explains the baseline MBL simulation script that reproduces the three-disorder-strength imbalance curves.

---

### 1. The `quantum_operator` Template Pattern
```python
operator_dict = dict(H0=operator_list_0)
for i in range(L):
    operator_dict["n" + str(i)] = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]
H_dict = quantum_operator(operator_dict, basis=basis, **no_checks)
```
Instead of building a separate Hamiltonian for each disorder realization, the script builds one `quantum_operator` template. The clean part `H0` contains the hopping and interaction terms (fixed). The disorder potentials are injected as named parameters (`n0` through `n7`) at realization time. This saves significant overhead — the basis and operator structure are computed only once.

### 2. Disorder Injection per Realization
```python
params_dict = dict(H0=1.0)
for j in range(L):
    params_dict["n" + str(j)] = uniform(-w, w)
H = H_dict.tohamiltonian(params_dict)
```
Each realization draws 8 fresh random on-site energies $\epsilon_j \sim \text{Uniform}(-W, W)$ and assembles the full Hamiltonian in a single `tohamiltonian()` call. The `np.random.seed(42)` at the top ensures every run of the script produces identical results — essential for reproducibility.

### 3. Time Evolution via `exp_op`
```python
U_op = exp_op(H, a=-1j, start=t.min(), stop=t.max(), num=len(t), iterate=True)
psi_t = U_op.dot(psi_0)
obs_t = obs_vs_time(psi_t, t_grid, dict(I=I_op))
```
QuSpin's `exp_op` iterator computes $e^{-iHt}|\psi_0\rangle$ across a linearly spaced time grid using Krylov subspace methods. The `iterate=True` flag makes it a generator — each step is computed on-demand, minimizing memory usage. `obs_vs_time` extracts $\langle \mathcal{I}(t) \rangle$ at every time point without storing the full state history.

### 4. Bootstrap Error Estimation
```python
boot_means = np.array([
    I_data[choice(n_real, size=n_real)].mean(axis=0)
    for _ in range(n_boot)
])
I_error = np.std(boot_means, axis=0)
```
The 100 realization trajectories are stored in `I_data` (shape: `n_real × n_times`). The bootstrap procedure resamples rows with replacement 100 times, computing the mean trajectory each time. The standard deviation of these bootstrap means gives **pointwise error bars** at every time step — more statistically robust than the standard error of the mean for correlated, heavy-tailed disorder distributions.

### 5. The Initial State Assertion
```python
assert abs(I_op.expt_value(psi_0).real - 1.0) < 1e-10
```
This hard assertion verifies that the Néel state has exactly $\mathcal{I}(0) = 1.0$ before any realization runs. If the string indexing is wrong (e.g., string length mismatch or wrong site ordering), the assertion fails immediately rather than silently producing incorrect physics.
