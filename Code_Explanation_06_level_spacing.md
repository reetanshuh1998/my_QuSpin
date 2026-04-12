# Code Walkthrough: `06_level_spacing.py`

This document explains the implementation of the level spacing ratio computation in `06_level_spacing.py`.

---

### 1. Parameter Design: Dense Grid Near the Transition
```python
n_real = 300
w_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0]
```
The disorder grid is deliberately concentrated around $W \approx 3$–$5$ where the MBL phase transition occurs for this model ($U=5$, $L=8$). Points at $W=0.5$ and $W=1.0$ anchor the ergodic (GOE) side, while $W=8$ and $W=10$ anchor the localized (Poisson) side. The 300 disorder realizations per $W$ point ensures statistically robust error bars, particularly in the critical transition region where fluctuations are largest.

### 2. Optimized Eigenvalue Computation
```python
E = H.eigvalsh()
```
Unlike the entanglement entropy scripts which need full eigenvectors ($V$), level spacing statistics require **only the eigenvalues**. QuSpin's `.eigvalsh()` method (backed by LAPACK's `dsyevd`) exploits this by skipping the expensive eigenvector computation entirely. For a $784 \times 784$ matrix, this runs roughly 2–3× faster than a full `eigh()` call.

### 3. Middle-Spectrum Filtering
```python
cut = int(dim * 0.25)
E_mid = E[cut:-cut]
```
The code discards the bottom 25% and top 25% of eigenvalues, retaining only the central 50% of the spectrum. This is essential because:
*   The density of states vanishes at the band edges, creating spurious level clustering that would bias $\langle r \rangle$ toward Poisson regardless of the actual phase.
*   The MBL transition is most cleanly defined in the bulk of the spectrum, away from edge effects.

### 4. The $r$-Ratio Computation
```python
deltaE = np.diff(E_mid)
r_n = np.minimum(deltaE[:-1], deltaE[1:]) / np.maximum(deltaE[:-1], deltaE[1:])
```
This is a vectorized implementation of the Oganesyan-Huse ratio:
1. `np.diff(E_mid)` computes all consecutive energy gaps $\delta_n = E_{n+1} - E_n$.
2. For each pair of adjacent gaps $(\delta_n, \delta_{n+1})$, the ratio $r_n = \min/\max$ is computed element-wise.
3. `np.mean(r_n)` returns the single-realization mean ratio.

The beauty of the $r$-ratio is that it is bounded $r \in [0, 1]$ and requires **no spectrum unfolding** — a significant practical advantage over the traditional level spacing distribution $P(s)$.

### 5. Bootstrap Error Estimation
```python
boot_gen = (np.mean(choice(r_data, size=n_real)) for _ in range(n_boot))
r_err = np.std(list(boot_gen))
```
For each disorder strength $W$, the script collects 300 scalar $\langle r \rangle$ values (one per realization). The bootstrap procedure resamples these 300 values with replacement 50 times, computing the mean of each resample. The standard deviation of these bootstrap means gives the statistical uncertainty on $\langle r \rangle(W)$.

### 6. Reference Lines
```python
ax.axhline(y=0.5307, ...)  # GOE (ergodic)
ax.axhline(y=0.386, ...)   # Poisson (MBL)
```
These are the analytically known values (Atas et al., PRL 2013):
*   GOE: $\langle r \rangle_{GOE} = 4 - 2\sqrt{3} \approx 0.5359$ (the code uses the numerically observed value 0.5307 for finite matrices, which is standard).
*   Poisson: $\langle r \rangle_{Poisson} = 2\ln 2 - 1 \approx 0.386$.

The horizontal reference lines allow instant visual identification of which phase the system occupies at each disorder strength.
