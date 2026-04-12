# Code Walkthrough: `05_entanglement_entropy.py`

This document provides a detailed breakdown of how `05_entanglement_entropy.py` computes the half-chain entanglement entropy $S(t)$ using Full Exact Diagonalization and QuSpin's built-in entropy routines.

---

### 1. The Log-Spaced Time Grid
```python
t_eval = np.logspace(-1, 3, 50)
```
Unlike scripts `01` and `04` which use `np.linspace` (uniformly spaced time points), this script uses `np.logspace` to generate 50 time points logarithmically distributed from $Jt = 0.1$ to $Jt = 1000$.

This is **critical** because the observable $S(t) \sim \ln(t)$ is logarithmic. On a linear grid, the first decade ($0.1$ to $1.0$) would have very few points while the last decade ($100$ to $1000$) would have thousands — wasting computation on the saturated regime and under-sampling the growth regime. The log grid ensures equal resolution across all decades.

> **Historical Bug Note:** The earlier version of this script used QuSpin's `exp_op(iterate=True)`, which internally calls `np.linspace(start, stop, num)` regardless of the input array. This silently destroyed the log-spacing and produced an artificial curve. The fix was to switch to Full Exact Diagonalization, which evaluates $|\psi(t)\rangle$ at exactly the requested time points.

### 2. Full Exact Diagonalization (The Critical Fix)
```python
H_mat = H.toarray()
E, V = np.linalg.eigh(H_mat)
c = V.conj().T @ psi_0
```
The function `np.linalg.eigh` (backed by LAPACK's `dsyevd`) fully diagonalizes the $784 \times 784$ Hermitian Hamiltonian matrix in one shot. This yields:
*   `E`: All 784 eigenvalues (energy levels), sorted in ascending order.
*   `V`: The corresponding eigenvectors as columns of a unitary matrix.
*   `c`: The expansion coefficients $c_n = \langle n | \psi_0 \rangle$ of the initial Néel state in the eigenbasis.

Because $L=8$ with quarter-filling gives dimension $784$, the diagonalization takes only $\sim 0.1$–$0.5$ seconds per call — extremely fast for exact results.

### 3. Exact Time Evolution at Arbitrary Times
```python
for it, tt in enumerate(t_eval):
    psi_t = V @ (c * np.exp(-1j * E * tt))
```
This is the mathematical heart of the script. Given the eigendecomposition, the time-evolved state is computed analytically:
$$|\psi(t)\rangle = \sum_n c_n \, e^{-iE_n t} \, |n\rangle = V \cdot (\vec{c} \odot e^{-i\vec{E}t})$$
where $\odot$ denotes element-wise multiplication. 

This formula can evaluate $|\psi(t)\rangle$ at **any** time $t$ — whether linearly spaced, logarithmically spaced, or even a single isolated point. This is precisely why Full ED was necessary: `exp_op` could only step forward sequentially through a fixed linear grid.

### 4. Computing Entanglement Entropy via QuSpin
```python
subsys_A = (range(L // 2), range(L // 2))
# ...
ent = basis.ent_entropy(psi_t, sub_sys_A=subsys_A, density=False)
S_t[it] = ent['Sent_A']
```
QuSpin's `basis.ent_entropy()` performs the following internally:
1. Reshapes the state vector $|\psi(t)\rangle$ into a bipartite tensor $\Psi_{ij}$ where $i$ indexes the left half-chain (subsystem A) and $j$ indexes the right half-chain (subsystem B).
2. Computes the reduced density matrix $\rho_A = \text{Tr}_B(|\psi\rangle\langle\psi|)$ via a partial trace.
3. Diagonalizes $\rho_A$ to obtain its eigenvalues $\lambda_k$.
4. Returns the Von Neumann entropy: $S = -\sum_k \lambda_k \ln \lambda_k$.

The `sub_sys_A` argument is a tuple of two ranges because the basis is `spinful_fermion_basis_1d`: the first range specifies which **sites** belong to subsystem A for spin-up fermions, and the second range specifies the same for spin-down fermions. Using `range(L//2)` for both means subsystem A is the left half of the chain for both spin species.

The flag `density=False` tells QuSpin that `psi_t` is a pure state vector, not a density matrix.

### 5. The Interaction Sweep Loop
```python
for U_val in U_list:
    # ... rebuild operator_list_0 with new int_list ...
    H_dict_obj = quantum_operator(operator_dict, basis=basis, **no_checks)
```
Like in Script `04`, the Hamiltonian template is rebuilt inside the loop for each value of $U$. The hopping terms remain fixed; only the `"n|n"` (Hubbard repulsion) coupling list changes. The disorder operator slots (`n0` through `n7`) are injected fresh for each realization inside `run_realization_ent`.

### 6. Semi-Log Plotting
```python
ax.set_xscale('log')
```
The final plot uses a logarithmic x-axis. On this scale:
*   Logarithmic growth $S(t) \sim \alpha \ln(t)$ appears as a **straight line** with slope $\alpha$.
*   Saturation (Anderson case, $U=0$) appears as a **flat horizontal line**.
*   The clear visual separation between these two behaviors is the primary evidence for MBL presented in this plot.
