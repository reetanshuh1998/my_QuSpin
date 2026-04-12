# Code Walkthrough: `04_strong_interaction_effect.py`

This document serves as an explicit mechanical breakdown of the Python script `04_strong_interaction_effect.py`, detailing how variable particle interactions ($U$) are mathematically analyzed using QuSpin across exact evaluation timeframes.

---

### 1. Dynamic Hamiltonian Generation
```python
operator_list_0 = [
    ["+-|", hop_left ],
    ["-+|", hop_right],
    ["|+-", hop_left ],
    ["|-+", hop_right],
    ["n|n", int_list ],
]
```
Rather than decoupling the Interaction variable `U_val` and holding it static like in previous scripts, the central `for U_val in U_list:` loop actively rebuilds the core `operator_list_0` topology on every iteration. The `"n|n"` string commands QuSpin to map exactly the Fermi-Hubbard repulsion: $U \sum n_{i \uparrow} n_{i \downarrow}$.
This reconstructs the internal Hamiltonian dimension accurately before the disorder fields `n0`...`nX` are statically injected.

### 2. Execution of Absolute Dynamics
```python
U_op = exp_op(H, a=-1j, start=t.min(), stop=t.max(), num=len(t), iterate=True)
psi_t = U_op.dot(psi_0)
t_grid = U_op.grid
obs_t = obs_vs_time(psi_t, t_grid, dict(I=I_op))
```
Because this script calculates observables in a purely linear $Jt$ timeline, it bypasses the heavy $O(N_s^3)$ complexity of computing Full Exact Diagonalization with `numpy.linalg.eigh()`. 
Instead, we apply QuSpin's `exp_op` iterator. This implements highly efficient Krylov subspace approximations to exactly evaluate the action of $e^{-iHt} |\psi_0\rangle$ dynamically propagating along the $35.0 Jt$ timeframe. The `obs_vs_time` function acts as a wrapper to instantly dump all calculated Imbalance values $\langle \mathcal{I} \rangle$ simultaneously, minimizing overhead time per realization to fractions of a second.

### 3. Extracting the Phase Threshold $\mathcal{I}_\infty$
```python
t_mask = t >= 25.0
I_inf = I_avg[t_mask].mean()
boot_I_inf = np.array([b[t_mask].mean() for b in boot_means])
I_inf_err = np.std(boot_I_inf)
```
The script attempts to prove the breakdown of theoretical quantum memory using a fixed long-time parameter array. To calculate the asymptotic behavior, the script forces a Numpy mask to completely cut off the energetic transients of the first $25 Jt$ units.
It then isolates only the $t \in [25, 35]$ final flat plateau. Rather than simply evaluating an incorrect cross-averaging array (which would destroy statistical error mapping), the array aggregates pure Bootstrapping elements (`boot_means`) matching the subset timeframe to output the single point scalar `I_inf` safely mapped alongside `I_inf_err`.

### 4. Twin Plot Resolution
Finally, the script separates the observable rendering into two completely distinct physical concepts:
*   **Plot 1 (`MBL_interaction_effect.png`):** Dumps the entirety of the real-time dynamics, creating the visible nested plateau graphs proving time continuity.
*   **Plot 2 (`MBL_interaction_phase.png`):** Drops the element of time entirely, charting strictly `I_inf` (Y axis) directly against the interaction threshold $U/J$ (X axis) to visually prove the threshold boundaries of Mondaini & Rigol (PRB 2015).
