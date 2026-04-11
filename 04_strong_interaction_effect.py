"""
Phase 4: Effect of Strong Interactions on MBL
----------------------------------------------
Study how varying the on-site interaction U affects
Many-Body Localization in the 1D disordered Fermi-Hubbard model.

We fix the disorder strength w = 4.0 (near the MBL transition)
and vary U/J = 0.0, 2.0, 5.0, 10.0

Observables:
  1. Full I(t) curves for select U values
  2. Long-time imbalance I_inf = <I(t)>_{t in [25,35]} vs U

Parameters:
    L = 8, J = 1.0, w = 4.0
    n_real = 100 disorder realizations
"""

from quspin.operators import hamiltonian, exp_op, quantum_operator
from quspin.basis import spinful_fermion_basis_1d
from quspin.tools.measurements import obs_vs_time
import numpy as np
from numpy.random import uniform, choice
from time import time
import matplotlib.pyplot as plt

#####################################################################
# Parameters
#####################################################################
n_real = 100
n_boot = 100

L = 8
N = L // 2
N_up = N // 2 + N % 2
N_down = N // 2

w = 4.0            # fixed disorder strength (near MBL transition)
J_hop = 1.0        # hopping
U_list = [0.0, 2.0, 5.0, 10.0]  # interaction strengths to scan

start, stop, num = 0.0, 35.0, 101
t = np.linspace(start, stop, num=num, endpoint=True)

#####################################################################
# Build basis (same for all U)
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
print(f"Hilbert space dimension: {basis.Ns}")

#####################################################################
# Site-coupling lists (hopping is U-independent)
#####################################################################
hop_right = [[-J_hop, i, i + 1] for i in range(L - 1)]
hop_left  = [[ J_hop, i, i + 1] for i in range(L - 1)]
sublat_list = [[(-1.0)**i / N, i] for i in range(L)]

imbalance_list = [["n|", sublat_list], ["|n", sublat_list]]

#####################################################################
# Initial state
#####################################################################
s_up   = "".join("1000" for _ in range(N_up))
s_down = "".join("0010" for _ in range(N_down))
no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
i_0    = basis.index(s_up, s_down)
psi_0  = np.zeros(basis.Ns)
psi_0[i_0] = 1.0

# imbalance observable (U-independent)
I_op = hamiltonian(imbalance_list, [], basis=basis, **no_checks)

#####################################################################
# Realization function
#####################################################################
def run_realization(H_dict, I_op, psi_0, w, t, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    U_op = exp_op(H, a=-1j, start=t.min(), stop=t.max(),
                  num=len(t), iterate=True)
    psi_t = U_op.dot(psi_0)
    t_grid = U_op.grid
    obs_t = obs_vs_time(psi_t, t_grid, dict(I=I_op))
    print(f"    realization {i_real+1}/{n_real} done in {time()-ti:.2f} s")
    return obs_t["I"].real

#####################################################################
# Main loop over interaction strengths
#####################################################################
results = {}  # U -> (I_avg, I_error)

for U_val in U_list:
    print(f"\n=== U/J = {U_val:.1f} ===")

    # rebuild operator dict with this U
    int_list = [[U_val, i, i] for i in range(L)]
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

    H_dict = quantum_operator(operator_dict, basis=basis, **no_checks)

    t0 = time()
    I_data = np.vstack([
        run_realization(H_dict, I_op, psi_0, w, t, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")

    I_avg = I_data.mean(axis=0)
    bootstrap_gen = (
        I_data[choice(n_real, size=n_real)].mean(axis=0)
        for _ in range(n_boot)
    )
    sq_fluc_gen = ((bootstrap - I_avg)**2 for bootstrap in bootstrap_gen)
    I_error = np.sqrt(sum(sq_fluc_gen) / n_boot)

    results[U_val] = (I_avg, I_error)

#####################################################################
# Plot 1: Full I(t) curves for all U values
#####################################################################
fig1, ax1 = plt.subplots(figsize=(9, 5.5))
cmap = plt.cm.viridis
colors = [cmap(i / (len(U_list) - 1)) for i in range(len(U_list))]

for idx, U_val in enumerate(U_list):
    I_avg, I_error = results[U_val]
    ax1.errorbar(t, I_avg, I_error, marker=".", markersize=3,
                 color=colors[idx], label=f"U/J={U_val:.1f}",
                 linewidth=0.8, elinewidth=0.4, capsize=1)

ax1.set_xlabel(r"$Jt$", fontsize=18)
ax1.set_ylabel(r"$\mathcal{I}$", fontsize=18)
ax1.set_title(f"MBL: Interaction Effect (w={w:.1f}, L={L})", fontsize=16)
ax1.tick_params(labelsize=14)
ax1.legend(loc="best", fontsize=12)
ax1.set_xlim(0, 35)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("MBL_interaction_effect.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_interaction_effect.png")
plt.close(fig1)

#####################################################################
# Plot 2: Long-time imbalance vs U
#####################################################################
fig2, ax2 = plt.subplots(figsize=(7, 5))
t_mask = t >= 25.0  # long-time window [25, 35]

I_inf_vals = []
I_inf_errs = []
for U_val in U_list:
    I_avg, I_error = results[U_val]
    I_inf = I_avg[t_mask].mean()
    I_inf_err = I_error[t_mask].mean()  # approximate error
    I_inf_vals.append(I_inf)
    I_inf_errs.append(I_inf_err)

ax2.errorbar(U_list, I_inf_vals, I_inf_errs, marker="o", markersize=8,
             color='#d62728', linewidth=2, capsize=4,
             label=f"w = {w:.1f}")
ax2.set_xlabel(r"$U / J$", fontsize=18)
ax2.set_ylabel(r"$\mathcal{I}_\infty$", fontsize=18)
ax2.set_title(f"Long-Time Imbalance vs Interaction (L={L})", fontsize=16)
ax2.tick_params(labelsize=14)
ax2.legend(fontsize=14)
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("MBL_interaction_phase.png", dpi=150, bbox_inches="tight")
print("Saved: MBL_interaction_phase.png")
plt.close(fig2)
