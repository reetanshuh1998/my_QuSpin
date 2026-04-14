# Code Walkthrough: `08_filling_level_spacing.py`

This document explains the filling-fraction level spacing scan, which probes spectral statistics across four particle numbers at fixed $U=5$, $W=10$.

---

### 1. `eigvalsh()` — Eigenvalues Only
```python
E = H.eigvalsh()
```
Unlike entanglement entropy (which needs eigenvectors to compute $\rho_A$), the level spacing ratio requires **only the sorted eigenvalues**. QuSpin's `.eigvalsh()` skips the computationally expensive eigenvector computation entirely and returns only the eigenvalues in ascending order. For a dim=4900 matrix, this is approximately **3× faster** than a full `eigh()` call — which is critical when running 100 realizations at the heaviest filling.

### 2. Middle-Spectrum Filter
```python
cut = int(dim * 0.25)
E_mid = E[cut:-cut]
```
The code discards the lowest 25% and highest 25% of energy levels, retaining only the central 50%. This is essential for two reasons:
1. The density of states goes to zero at the band edges, causing artificial level clustering that biases $\langle r \rangle$ toward Poisson regardless of phase.
2. The MBL transition is energy-dependent (there is a many-body mobility edge). The bulk of the spectrum is the most representative and homogeneous region.

With 100 realizations and dim ranging from 64 to 4900, the number of mid-spectrum levels per realization ranges from ~32 (tiny, noisy) to ~2450 (reliable). This is why the $N=2$ data point will have larger error bars than $N=8$.

### 3. Vectorized $r$-Ratio Computation
```python
deltaE = np.diff(E_mid)
r_n = np.minimum(deltaE[:-1], deltaE[1:]) / np.maximum(deltaE[:-1], deltaE[1:])
return np.mean(r_n)
```
This is a fully vectorized implementation of the Oganesyan-Huse ratio (PRB 2007). `np.diff` computes all consecutive gaps in a single array operation. `np.minimum` and `np.maximum` then broadcast element-wise to compute $r_n = \min(\delta_n, \delta_{n+1}) / \max(\delta_n, \delta_{n+1})$ for all pairs simultaneously. The function returns a single scalar $\langle r \rangle$ per realization, which is then aggregated across all 100 realizations.

### 4. Per-Filling Data Aggregation
```python
r_data = np.array([
    run_realization_rstat(H_dict_obj, w_val, dim, i, n_real)
    for i in range(n_real)
])
r_avg = np.mean(r_data)
```
Each element of `r_data` is one $\langle r \rangle$ value from one disorder realization. The final plot uses `np.mean(r_data)` with bootstrap error bars. Unlike the time-resolved scripts, the result here is a single scalar per disorder value per filling — the final plot is therefore a clean bar/line chart rather than a time series.

### 5. Plot Design: Reference Lines + Dimension Annotation
```python
ax.axhline(y=0.5307, ...)  # GOE
ax.axhline(y=0.386, ...)   # Poisson
ax.annotate(f"dim={dims[i]}", ...)
```
The two horizontal reference lines are the analytically known values that define the two phases. Every data point is additionally annotated with its Hilbert space dimension, making it transparent to the reader exactly how much spectral information underlies each $\langle r \rangle$ value. This is important because the $N=2$ (dim=64) point is computed from far fewer levels than the $N=8$ (dim=4900) point, and the annotation lets the reader calibrate their confidence accordingly.
