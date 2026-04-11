"""
Phase 5: Entanglement Entropy Growth
------------------------------------
Reproduces the canonical Bardarson, Pollmann, Moore (2012) setup.
We compute the Half-Chain von Neumann Entanglement Entropy S(t).
We fix a strong disorder (W=10.0) to ensure the system is localized.
We sweep the interaction strength U = 0.0, 1.0, 5.0.

- U = 0.0 (Anderson Localization): S(t) saturates rapidly.
- U > 0.0 (Many-Body Localization): S(t) grows logarithmically S(t) ~ ln(t).
"""

from quspin.operators import quantum_operator
from quspin.basis import spinful_fermion_basis_1d
from quspin.operators import exp_op
import numpy as np
from numpy.random import uniform, choice
from time import time
import matplotlib.pyplot as plt

#####################################################################
# Parameters
#####################################################################
n_real = 30
n_boot = 50

L = 8
N = L // 2
N_up = N // 2 + N % 2
N_down = N // 2

U_list = [0.0, 1.0, 5.0, 10.0]
J_hop = 1.0
w_val = 10.0  # Fixed deep in the localized regime

# Log-spaced time array to clearly see logarithmic growth
t_eval = np.logspace(-1, 3, 50)

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
print(f"Hilbert space dimension: {basis.Ns}")

#####################################################################
# Initial pure Fock state
#####################################################################
# As asked, we start from a pure product state. 
# We use the standard N\'eel state for N_up=N_down.
s_up   = "".join("10" for _ in range(N_up))
s_down = "".join("01" for _ in range(N_down))
i_0    = basis.index(s_up, s_down)
psi_0  = np.zeros(basis.Ns)
psi_0[i_0] = 1.0

subsys_A = (range(L // 2), range(L // 2))

#####################################################################
# Realization function
#####################################################################
def run_realization_ent(H_dict, w, t_eval, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    
    # Time evolution operator
    U_op = exp_op(H, a=-1j, start=t_eval.min(), stop=t_eval.max(), 
                  num=len(t_eval), iterate=True)
    
    S_t = np.zeros(len(t_eval))
    for it, psi_t in enumerate(U_op.dot(psi_0)):
        # Calculate entanglement entropy
        ent = basis.ent_entropy(psi_t, sub_sys_A=subsys_A, density=False)
        S_t[it] = ent['Sent_A']
        
    elapsed = time() - ti
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
    return S_t

#####################################################################
# Main loop
#####################################################################
results = {}

for U_val in U_list:
    print(f"\n=== Interaction U = {U_val:.1f} ===")
    t0 = time()
    
    # Build H_dict for this specific U
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
    
    S_data = np.vstack([
        run_realization_ent(H_dict_obj, w_val, t_eval, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")
    
    S_avg = S_data.mean(axis=0)
    boot_gen = (S_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
    sq_gen = ((b - S_avg)**2 for b in boot_gen)
    S_error = np.sqrt(sum(sq_gen) / n_boot)
    
    results[U_val] = (S_avg, S_error)

#####################################################################
# Plotting like Bardarson Fig 1
#####################################################################
fig, ax = plt.subplots(figsize=(9, 6))

colors = ['navy', 'forestgreen', 'darkorange', 'crimson']
markers = ['o', 's', '^', 'D']

for idx, U_val in enumerate(U_list):
    S_avg, S_error = results[U_val]
    ax.errorbar(t_eval, S_avg, S_error, marker=markers[idx], markersize=5,
                color=colors[idx], label=rf"$U/J = {U_val:.1f}$",
                linewidth=1.5, elinewidth=0.8, capsize=2)

ax.set_xscale('log')
ax.set_xlabel(r"$Jt$", fontsize=18)
ax.set_ylabel(r"$S$", fontsize=18)
ax.set_title(rf"Entanglement Growth (Fixed Disorder $W={w_val:.1f}$)", fontsize=16)

# Mimic the publication axis
ax.set_xlim([1e-1, 1e3])
ax.tick_params(labelsize=14)
ax.legend(loc="upper left", fontsize=14, framealpha=0.9)
ax.grid(True, alpha=0.3, which="both")

fig.tight_layout()
fig.savefig("MBL_entanglement_growth_U_sweep.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_entanglement_growth_U_sweep.png")
plt.close(fig)
