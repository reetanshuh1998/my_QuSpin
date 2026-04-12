"""
Phase 2: Effect of Finite Temperature on MBL (Microcanonical Approach)
----------------------------------------------------------------------
Instead of the canonical thermal density matrix, we use Energy-Shell 
Sampling (Luitz et al., 2015). We map the thermal states to an energy density eps.

Protocol:
1. Start from the Néel Fock state |psi_Neel>
2. Diagonalize disordered H to get eigenvalues E_n and eigenstates |n>
3. Define eps_n = (E_n - E_min) / (E_max - E_min)
4. Select target energy density eps corresponding to an effective T.
5. Create a pure-state projection within an energy window around eps, 
   weighted by their overlap with the initial Neel state.
6. Evolve this pure state natively.

Parameters:
    L=8, J=1.0, U=5.0, w=4.0
"""
from quspin.operators import hamiltonian, quantum_operator
from quspin.basis import spinful_fermion_basis_1d
import numpy as np
from time import time
from numpy.random import uniform, choice
import matplotlib.pyplot as plt

np.random.seed(42)

n_real = 50
n_boot = 100
L = 8
N = L              # half-filling
N_up = L // 2      # = 4
N_down = L // 2    # = 4

w = 4.0
J_hop = 1.0
U_val = 5.0

# Energy density targets (representing different Effective Temperatures)
# 0.1 -> Cold (close to ground state)
# 0.5 -> Infinite T (mid-spectrum)
# 0.9 -> Negative T (highly energetic bounds)
eps_list = [0.1, 0.5, 0.9] 

start, stop, num = 0.0, 35.0, 71
t = np.linspace(start, stop, num=num, endpoint=True)

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
print(f"Hilbert space dimension: {basis.Ns}")

#####################################################################
# Structural Definitions
#####################################################################
hop_right = [[-J_hop, i, i + 1] for i in range(L - 1)]
hop_left  = [[ J_hop, i, i + 1] for i in range(L - 1)]
int_list  = [[U_val, i, i] for i in range(L)]
sublat_list = [[(-1.0)**i / N, i] for i in range(L)]

operator_list_0 = [
    ["+-|", hop_left ],
    ["-+|", hop_right],
    ["|+-", hop_left ],
    ["|-+", hop_right],
    ["n|n", int_list ],
]
imbalance_list = [["n|", sublat_list], ["|n", sublat_list]]

operator_dict = dict(H0=operator_list_0)
for i in range(L):
    operator_dict["n" + str(i)] = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]

no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
H_dict_template = quantum_operator(operator_dict, basis=basis, **no_checks)
I_op = hamiltonian(imbalance_list, [], basis=basis, **no_checks)
I_mat = I_op.toarray()

# Pure initial state (CDW Néel State — all particles on even sites)
s_up   = "10101010"
s_down = "10101010"
i_0    = basis.index(s_up, s_down)
psi_Neel = np.zeros(basis.Ns)
psi_Neel[i_0] = 1.0

# Safety trigger
assert abs(I_op.expt_value(psi_Neel).real - 1.0) < 1e-10, "Initial imbalance must be 1.0"

#####################################################################
# Microcanonical Shell Decomposition
#####################################################################
def compute_imbalance_shell(H_dict, I_mat, eps_target, w, t, i_real, n_real):
    ti = time()
    
    # 1. Map disordered H
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    
    # 2. Complete Exact Diagonalization
    E_n, V_n = np.linalg.eigh(H.toarray())
    
    # 3. Calculate overlap constants for the initial Neel State
    c_n = V_n.T.conj() @ psi_Neel
    
    E_min = E_n[0]
    E_max = E_n[-1]
    
    # 4. Filter projection spectrum
    E_target = E_min + eps_target * (E_max - E_min)
    delta_E = 0.05 * (E_max - E_min)
    shell_mask = np.abs(E_n - E_target) <= delta_E
    
    # Edge-case safety 
    if not np.any(shell_mask):
        idx = np.argmin(np.abs(E_n - E_target))
        shell_mask[idx] = True
        
    # Project Neel state strictly onto the isolated Energy Shell
    psi_eps_eig = np.zeros_like(E_n, dtype=np.complex128)
    psi_eps_eig[shell_mask] = c_n[shell_mask]
    
    # Formal Renormalization of probabilities in new spectrum
    norm = np.linalg.norm(psi_eps_eig)
    if norm < 1e-12:
        psi_eps_eig[shell_mask] = 1.0 / np.sqrt(np.sum(shell_mask))
    else:
        psi_eps_eig /= norm
        
    # 5. Native Quantum Dynamics I(t) = <psi(t)| I_eig |psi(t)>
    I_eig = V_n.T.conj() @ I_mat @ V_n
    A = np.outer(psi_eps_eig.conj(), psi_eps_eig) * I_eig
    dE = E_n[:, np.newaxis] - E_n[np.newaxis, :]
    
    I_t = np.zeros(len(t))
    for it, tt in enumerate(t):
        I_t[it] = np.sum(A * np.exp(-1j * dE * tt)).real
        
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {time()-ti:.2f} s")
    return I_t

#####################################################################
# Loop Phase
#####################################################################
results_eps = {}
for eps in eps_list:
    print(f"\n=== Energy Density eps = {eps:.2f} ===")
    t0 = time()
    I_data = np.vstack([
        compute_imbalance_shell(H_dict_template, I_mat, eps, w, t, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")

    I_avg = I_data.mean(axis=0)
    boot_gen = (I_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
    sq_gen = ((b - I_avg)**2 for b in boot_gen)
    I_error = np.sqrt(sum(sq_gen) / n_boot)

    results_eps[eps] = (I_avg, I_error)

#####################################################################
# Pure Neel state evaluation limit (Equivalent to unrestricted integration)
#####################################################################
print("\n=== Pure Neel state (Full Unrestricted Spectrum) ===")
# We simply calculate exact decomposition taking the complete mask
I_data_pure = np.zeros((n_real, len(t)))
for i in range(n_real):
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict_template.tohamiltonian(params_dict)
    
    E_n, V_n = np.linalg.eigh(H.toarray())
    c_n = V_n.T.conj() @ psi_Neel
    
    I_eig = V_n.T.conj() @ I_mat @ V_n
    A = np.outer(c_n.conj(), c_n) * I_eig
    dE = E_n[:, np.newaxis] - E_n[np.newaxis, :]
    for it, tt in enumerate(t):
        I_data_pure[i, it] = np.sum(A * np.exp(-1j * dE * tt)).real
        
I_avg_pure = I_data_pure.mean(axis=0)
boot_gen = (I_data_pure[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
I_err_pure = np.sqrt(sum(((b - I_avg_pure)**2 for b in boot_gen)) / n_boot)

#####################################################################
# Graphical Representation
#####################################################################
fig, ax = plt.subplots(figsize=(9, 6))

ax.errorbar(t, I_avg_pure, I_err_pure, marker=".", markersize=3,
            color='black', label=r"Pure Néel State",
            linewidth=1.2, elinewidth=0.4, capsize=1)

colors_th = ['#1f77b4', '#ff7f0e', '#d62728']
labels_th = ["Low T (eps=0.1)", "Infinite T (eps=0.5)", "Negative T (eps=0.9)"]

for idx, eps in enumerate(eps_list):
    I_avg, I_error = results_eps[eps]
    ax.errorbar(t, I_avg, I_error, marker=".", markersize=3,
                color=colors_th[idx],
                label=rf"{labels_th[idx]}",
                linewidth=0.8, elinewidth=0.4, capsize=1)

ax.set_xlabel(r"$Jt$", fontsize=18)
ax.set_ylabel(r"$\mathcal{I}$", fontsize=18)
ax.set_title(f"MBL: Thermal Bounds via Microcanonical Shells (w={w:.1f}, L={L})", fontsize=16)
ax.tick_params(labelsize=14)
ax.legend(loc="best", fontsize=11)
ax.set_xlim(0, 35)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("MBL_temperature_effect_shell.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_temperature_effect_shell.png")
plt.close(fig)
