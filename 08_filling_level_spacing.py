"""
Phase 8: Filling-Fraction Dependence of Level Spacing Statistics
-----------------------------------------------------------------
Fixed U=5.0, W=10.0, L=8. Scan N_total = 2, 4, 6, 8 fermions.
Computes the mean level spacing ratio <r> for each filling.

Expected:
  - N=2 (dilute):  Poisson-like (r ~ 0.386), localized
  - N=4 (quarter):  Near Poisson (deep MBL at W=10)
  - N=6:            May shift toward GOE if interactions overcome disorder
  - N=8 (half):     Richest many-body effects
"""

from quspin.operators import quantum_operator
from quspin.basis import spinful_fermion_basis_1d
import numpy as np
from numpy.random import uniform, choice
from time import time
import matplotlib.pyplot as plt

np.random.seed(42)

#####################################################################
# Parameters
#####################################################################
n_real = 100
n_boot = 50

L = 8
U_val = 5.0
J_hop = 1.0
w_val = 10.0

filling_list = [1, 2, 3, 4]  # N_up = N_down

#####################################################################
# Realization function
#####################################################################
def run_realization_rstat(H_dict, w, dim, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    
    # Only need eigenvalues for level spacing
    E = H.eigvalsh()
    
    # Filter to middle 50% of spectrum to avoid edge effects
    cut = int(dim * 0.25)
    E_mid = E[cut:-cut]
    
    # Level spacings
    deltaE = np.diff(E_mid)
    
    # Ratio r_n = min(delta_n, delta_{n+1}) / max(delta_n, delta_{n+1})
    r_n = np.minimum(deltaE[:-1], deltaE[1:]) / np.maximum(deltaE[:-1], deltaE[1:])
    
    elapsed = time() - ti
    if (i_real + 1) % 20 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
    return np.mean(r_n)

#####################################################################
# Main loop
#####################################################################
results = {}

for n_fill in filling_list:
    N_up = n_fill
    N_down = n_fill
    N_total = N_up + N_down
    
    print(f"\n=== N_total = {N_total} (N_up={N_up}, N_down={N_down}) ===")
    
    basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
    dim = basis.Ns
    print(f"  Hilbert space dimension: {dim}")
    
    # Build Hamiltonian template
    hop_right = [[-J_hop, i, i + 1] for i in range(L - 1)]
    hop_left  = [[ J_hop, i, i + 1] for i in range(L - 1)]
    int_list  = [[U_val, i, i] for i in range(L)]

    operator_list_0 = [
        ["+-|", hop_left ], ["-+|", hop_right],
        ["|+-", hop_left ], ["|-+", hop_right],
        ["n|n", int_list ]
    ]
    operator_dict = dict(H0=operator_list_0)
    for i in range(L):
        operator_dict["n" + str(i)] = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]

    no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
    H_dict_obj = quantum_operator(operator_dict, basis=basis, **no_checks)
    
    t0 = time()
    r_data = np.array([
        run_realization_rstat(H_dict_obj, w_val, dim, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")
    
    r_avg = np.mean(r_data)
    boot_gen = (np.mean(choice(r_data, size=n_real)) for _ in range(n_boot))
    r_err = np.std(list(boot_gen))
    
    results[N_total] = (r_avg, r_err, dim)

#####################################################################
# Plot
#####################################################################
fig, ax = plt.subplots(figsize=(8, 6))

N_totals = [2 * n for n in filling_list]
r_avgs = [results[N][0] for N in N_totals]
r_errs = [results[N][1] for N in N_totals]
dims   = [results[N][2] for N in N_totals]

ax.errorbar(N_totals, r_avgs, yerr=r_errs, fmt="-o", color="crimson",
            linewidth=2, markersize=10, capsize=5,
            label=rf"$U={U_val:.0f}$, $W={w_val:.0f}$")

# Annotate Hilbert space dimensions
for i, N in enumerate(N_totals):
    ax.annotate(f"dim={dims[i]}", (N, r_avgs[i]),
                textcoords="offset points", xytext=(10, 10),
                fontsize=9, color='gray')

# Reference lines
ax.axhline(y=0.5307, color='navy', linestyle='--', alpha=0.7,
           label=r'GOE $\approx 0.53$ (Ergodic)')
ax.axhline(y=0.386, color='darkgreen', linestyle=':', alpha=0.7,
           label=r'Poisson $\approx 0.386$ (MBL)')

ax.set_xlabel(r"Total particle number $N$", fontsize=16)
ax.set_ylabel(r"Level spacing ratio $\langle r \rangle$", fontsize=16)
ax.set_title(rf"Spectral Statistics vs Filling ($U={U_val:.0f}$, $W={w_val:.0f}$, $L={L}$)",
             fontsize=15)
ax.set_xticks(N_totals)
ax.set_xticklabels([f"{N}\n({N}/{2*L})" for N in N_totals])
ax.tick_params(labelsize=13)
ax.legend(loc="best", fontsize=12)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig("MBL_filling_level_spacing.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_filling_level_spacing.png")
plt.close(fig)
