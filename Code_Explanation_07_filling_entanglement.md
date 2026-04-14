# Code Walkthrough: `07_filling_entanglement.py`

This document explains the filling-fraction entanglement entropy scan, which extends Script 05 across four particle numbers at fixed $U=5$, $W=10$.

---

### 1. Dynamic Basis Rebuilding
```python
for n_fill in filling_list:
    N_up = n_fill
    N_down = n_fill
    basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
```
Unlike previous scripts where the basis is built once outside the loop, this script rebuilds the `spinful_fermion_basis_1d` for each filling. This is necessary because the Hilbert space dimension changes dramatically with filling — from dim=64 ($N=2$) to dim=4900 ($N=8$). Each loop iteration operates in a completely different vector space.

### 2. Filling-Adaptive Initial State
```python
def build_initial_state(L, N_up, N_down, basis):
    s_up = ['0'] * L
    s_down = ['0'] * L
    even_sites = [0, 2, 4, 6]
    odd_sites  = [1, 3, 5, 7]
    for i in range(N_up):
        s_up[even_sites[i]] = '1'
    for i in range(N_down):
        s_down[odd_sites[i]] = '1'
```
Rather than a hard-coded string (which only works for a fixed filling), this helper function builds the initial state programmatically. It places up-spins on the first `N_up` even sites and down-spins on the first `N_down` odd sites. This gives the maximum spatial ordering for any filling fraction, ensuring all curves start from a comparably well-ordered (zero-entanglement) initial condition.

For example:
- $N=2$: `"10000000"` / `"01000000"` (1 particle each, far apart)
- $N=4$: `"10001000"` / `"00100010"` (2 particles each, alternating)
- $N=8$: `"10101010"` / `"01010101"` (4 particles each, fully alternating)

### 3. Full Exact Diagonalization at Each Filling
```python
H_mat = H.toarray()
E, V = np.linalg.eigh(H_mat)
c = V.conj().T @ psi_0
for it, tt in enumerate(t_eval):
    psi_t = V @ (c * np.exp(-1j * E * tt))
```
Full ED is used identically to Script 05 — this correctly handles the log-spaced time grid by evaluating at arbitrary $t$ values analytically. The computational cost scales as $O(\dim^3)$ for diagonalization, which is why $N=8$ (dim=4900) is much slower than $N=2$ (dim=64).

### 4. Hilbert-Space-Adaptive `subsys_A`
```python
subsys_A = (range(L // 2), range(L // 2))
```
This is defined inside the loop after the basis is created. It always specifies the left half of the chain (sites 0–3) for both spin species, regardless of filling. QuSpin's `ent_entropy` then correctly partitions the Hilbert space into the bipartition $A \otimes B$ and computes $S_A = -\text{Tr}(\rho_A \ln \rho_A)$.

### 5. Memory-Efficient Bootstrap via Generator
```python
boot_gen = (S_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
sq_gen = ((b - S_avg)**2 for b in boot_gen)
S_error = np.sqrt(sum(sq_gen) / n_boot)
```
This generator-based approach computes the bootstrap standard deviation without storing all 50 bootstrap sample arrays simultaneously. It is especially important at higher fillings where `S_data` has more rows and each row is longer. Memory usage remains bounded at $O(n\_real \times n\_times)$.
