"""
Phase 6: Level Spacing Statistics
---------------------------------
Analyzes the level spacing ratio r of the energy spectrum.
In the ergodic (thermal) phase, states exhibit level repulsion (Wigner-Dyson): r ~ 0.53
In the MBL (localized) phase, localized levels are uncorrelated (Poisson): r ~ 0.386

This requires full exact diagonalization across a range of disorders.
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
n_real = 100   # Increased to 100 for proper statistical representation
n_boot = 50

# Literature firmly requires L >= 8 to see the crossover sharply. 
# Anything smaller will be completely smeared out by finite-size effects.
L = 8         
N = L // 2
N_up = N // 2 + N % 2
N_down = N // 2

# Finer grid across the critical MBL transition region
w_list = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 11.0, 16.0]
J_hop = 1.0
U_val = 5.0

print(f"System: L={L}, U={U_val}, Realizations={n_real}")

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
dim = basis.Ns
print(f"Hilbert space dimension: {dim}")

#####################################################################
# Site-coupling lists
#####################################################################
hop_right = [[-J_hop, i, i + 1] for i in range(L - 1)]
hop_left  = [[ J_hop, i, i + 1] for i in range(L - 1)]
int_list  = [[U_val, i, i] for i in range(L)]

operator_list_0 = [
    ["+-|", hop_left ],
    ["-+|", hop_right],
    ["|+-", hop_left ],
    ["|-+", hop_right],
    ["n|n", int_list ],
]

operator_dict = dict(H0=operator_list_0)
for i in range(L):
    operator_dict["n" + str(i)] = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]

no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
H_dict_obj = quantum_operator(operator_dict, basis=basis, **no_checks)

#####################################################################
# Realization function
#####################################################################
def run_realization_rstat(H_dict, w, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    
    # Fast diagonalization: only need eigenvalues, not eigenvectors
    # eigvalsh is heavily optimized by LAPACK vs full eigh
    E = H.eigvalsh()
    
    # Filter to middle of the spectrum to avoid edge effects
    # Keep middle 50% of states
    cut = int(dim * 0.25)
    E_mid = E[cut:-cut]
    
    # Calculate level differences
    deltaE = np.diff(E_mid)
    
    # Calculate ratio r
    r_n = np.minimum(deltaE[:-1], deltaE[1:]) / np.maximum(deltaE[:-1], deltaE[1:])
    
    elapsed = time() - ti
    if (i_real + 1) % 5 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
        
    return np.mean(r_n)

#####################################################################
# Main loop
#####################################################################
r_avg_list = []
r_err_list = []

for w in w_list:
    print(f"\n=== Disorder w = {w:.1f} ===")
    t0 = time()
    
    r_data = np.array([
        run_realization_rstat(H_dict_obj, w, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")
    
    r_avg = np.mean(r_data)
    # Bootstrap error estimation
    boot_gen = (np.mean(choice(r_data, size=n_real)) for _ in range(n_boot))
    r_err = np.std(list(boot_gen))
    
    r_avg_list.append(r_avg)
    r_err_list.append(r_err)

#####################################################################
# Plot
#####################################################################
fig, ax = plt.subplots(figsize=(8, 6))

ax.errorbar(w_list, r_avg_list, yerr=r_err_list, fmt="-o", color="crimson",
            linewidth=2, markersize=8, capsize=4, label=f"Fermi-Hubbard (L={L})")

# Theoretical reference lines
ax.axhline(y=0.5307, color='navy', linestyle='--', alpha=0.7, 
           label=r'GOE $\approx 0.53$ (Ergodic)')
ax.axhline(y=0.386, color='darkgreen', linestyle=':', alpha=0.7, 
           label=r'Poisson $\approx 0.386$ (MBL)')

ax.set_xlabel(r"Disorder strength $W/J$", fontsize=16)
ax.set_ylabel(r"Level spacing ratio $\langle r \rangle$", fontsize=16)
ax.set_title("MBL Phase Transition: Spectral Statistics", fontsize=16)
ax.tick_params(labelsize=14)
ax.legend(loc="upper right", fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("MBL_level_spacing.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_level_spacing.png")
plt.close(fig)
