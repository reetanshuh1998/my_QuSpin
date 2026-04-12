"""
Phase 7: Filling-Fraction Dependence of MBL Entanglement Growth
----------------------------------------------------------------
Fixed U=5.0, W=10.0, L=8. Scan N_total = 2, 4, 6, 8 fermions.
Measures half-chain entanglement entropy S(t) for each filling.

Expected:
  - N=2 (dilute):  Anderson-like, S(t) saturates quickly
  - N=4 (quarter):  MBL, S(t) ~ alpha*ln(t), small alpha
  - N=6:            MBL, larger alpha (more interactions)
  - N=8 (half):     Maximum interactions, largest alpha or potential thermalization
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
n_real = 30
n_boot = 50

L = 8
U_val = 5.0
J_hop = 1.0
w_val = 10.0

# Filling scan: N_up = N_down = 1, 2, 3, 4
filling_list = [1, 2, 3, 4]  # N_up = N_down values

t_eval = np.logspace(-1, 3, 50)

#####################################################################
# Helper: build initial Neel-like product state for any filling
#####################################################################
def build_initial_state(L, N_up, N_down, basis):
    """
    Place up-spins on even sites (0, 2, 4, ...) and 
    down-spins on odd sites (1, 3, 5, ...).
    This gives a clean antiferromagnetic Neel pattern at half-filling,
    and a partial Neel for lower fillings.
    """
    s_up = ['0'] * L
    s_down = ['0'] * L
    
    even_sites = [0, 2, 4, 6]
    odd_sites  = [1, 3, 5, 7]
    
    for i in range(N_up):
        s_up[even_sites[i]] = '1'
    for i in range(N_down):
        s_down[odd_sites[i]] = '1'
    
    s_up_str = "".join(s_up)
    s_down_str = "".join(s_down)
    
    i_0 = basis.index(s_up_str, s_down_str)
    psi_0 = np.zeros(basis.Ns)
    psi_0[i_0] = 1.0
    
    print(f"  Initial state: |{s_up_str}>(up) x |{s_down_str}>(dn)")
    return psi_0

#####################################################################
# Realization function
#####################################################################
def run_realization_ent(H_dict, w, t_eval, psi_0, basis, subsys_A, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    
    # Full Exact Diagonalization
    H_mat = H.toarray()
    E, V = np.linalg.eigh(H_mat)
    c = V.conj().T @ psi_0
    
    S_t = np.zeros(len(t_eval))
    for it, tt in enumerate(t_eval):
        psi_t = V @ (c * np.exp(-1j * E * tt))
        ent = basis.ent_entropy(psi_t, sub_sys_A=subsys_A, density=False)
        S_t[it] = ent['Sent_A']
        
    elapsed = time() - ti
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
    return S_t

#####################################################################
# Main loop over filling fractions
#####################################################################
results = {}

for n_fill in filling_list:
    N_up = n_fill
    N_down = n_fill
    N_total = N_up + N_down
    
    print(f"\n=== N_total = {N_total} (N_up={N_up}, N_down={N_down}) ===")
    
    # Build basis for this filling
    basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
    print(f"  Hilbert space dimension: {basis.Ns}")
    
    # Build initial state
    psi_0 = build_initial_state(L, N_up, N_down, basis)
    
    # Subsystem A: left half-chain for both spin species
    subsys_A = (range(L // 2), range(L // 2))
    
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
    S_data = np.vstack([
        run_realization_ent(H_dict_obj, w_val, t_eval, psi_0, basis, subsys_A, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")
    
    S_avg = S_data.mean(axis=0)
    boot_gen = (S_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
    sq_gen = ((b - S_avg)**2 for b in boot_gen)
    S_error = np.sqrt(sum(sq_gen) / n_boot)
    
    results[N_total] = (S_avg, S_error, basis.Ns)

#####################################################################
# Plot
#####################################################################
fig, ax = plt.subplots(figsize=(9, 6))

colors = ['navy', 'forestgreen', 'darkorange', 'crimson']
markers = ['o', 's', '^', 'D']

for idx, n_fill in enumerate(filling_list):
    N_total = 2 * n_fill
    S_avg, S_error, dim = results[N_total]
    filling_frac = N_total / (2 * L)  # fraction of max (2L slots)
    ax.errorbar(t_eval, S_avg, S_error, marker=markers[idx], markersize=5,
                color=colors[idx],
                label=rf"$N={N_total}$ (dim={dim}, $\nu={filling_frac:.2f}$)",
                linewidth=1.5, elinewidth=0.8, capsize=2)

ax.set_xscale('log')
ax.set_xlabel(r"$Jt$", fontsize=18)
ax.set_ylabel(r"$S(t)$", fontsize=18)
ax.set_title(rf"Entanglement Growth vs Filling ($U={U_val:.0f}$, $W={w_val:.0f}$, $L={L}$)",
             fontsize=16)
ax.set_xlim([1e-1, 1e3])
ax.tick_params(labelsize=14)
ax.legend(loc="upper left", fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, which="both")

fig.tight_layout()
fig.savefig("MBL_filling_entanglement.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_filling_entanglement.png")
plt.close(fig)
