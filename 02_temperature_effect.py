"""
Phase 2: Effect of Finite Temperature on MBL
----------------------------------------------
Study how finite-temperature initial states affect
Many-Body Localization in the 1D disordered Fermi-Hubbard model.

Instead of a pure Fock state, we prepare a thermal density matrix
at inverse temperature beta and evolve it. The imbalance is computed
as I(t) = Tr[rho(t) I_hat] using eigenstate decomposition.

We fix disorder w = 4.0 and vary beta*J = 0.1, 1.0, 10.0
(high, intermediate, low temperature). We also show the pure-state
(beta -> infinity) result for comparison.

Parameters:
    L = 6 (reduced for density matrix tractability)
    J = 1.0, U = 5.0
    n_real = 50 disorder realizations
"""

from quspin.operators import hamiltonian, quantum_operator
from quspin.basis import spinful_fermion_basis_1d
import numpy as np
from numpy.random import uniform, choice
from time import time
import matplotlib.pyplot as plt

#####################################################################
# Parameters
#####################################################################
n_real = 50    # fewer realizations (density matrix is more expensive)
n_boot = 100

L = 6          # smaller system for density matrix approach
N = L // 2
N_up = N // 2 + N % 2
N_down = N // 2

w = 4.0
J_hop = 1.0
U_val = 5.0
beta_list = [10.0, 1.0, 0.1]   # inverse temperatures (beta*J)

start, stop, num = 0.0, 35.0, 101
t = np.linspace(start, stop, num=num, endpoint=True)

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
print(f"Hilbert space dimension: {basis.Ns}")

#####################################################################
# Site-coupling lists
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

# pure-state initial state (for comparison)
s_up   = "".join("100" for _ in range(N_up))
s_down = "".join("010" for _ in range(N_down))
i_0    = basis.index(s_up, s_down)
psi_0_pure = np.zeros(basis.Ns)
psi_0_pure[i_0] = 1.0

#####################################################################
# Realization function (full eigendecomposition)
#####################################################################
def run_realization_thermal(H_dict, I_mat, beta, w, t, i_real, n_real, dim):
    """
    For a single disorder realization:
    1. Build H and diagonalize fully
    2. Construct thermal density matrix rho_0 = exp(-beta H) / Z
    3. Evolve: rho(t) = U(t) rho_0 U^dag(t) where U(t) = exp(-iHt)
    4. Measure I(t) = Tr[rho(t) I_hat]
    
    Using eigendecomposition:
       rho_0 in eigenbasis: rho_{mn} = delta_{mn} * exp(-beta E_n) / Z
       I(t) = sum_{m,n} rho_{mn}(0) * I_{nm} * exp(i(E_n - E_m)t)
            = sum_n p_n * I_{nn}   +   sum_{m!=n} off-diagonal (zero for diag rho_0)
    
    For diagonal rho_0, I(t) = sum_n p_n * <n|I|n> = const (no dynamics!)
    
    This is because a thermal state of the DISORDERED Hamiltonian is stationary.
    So we use the thermal state of the CLEAN Hamiltonian instead, which is NOT
    stationary under the disordered Hamiltonian.
    """
    ti = time()
    
    # Disordered Hamiltonian
    params_dict = dict(H0=1.0)
    disorder = [uniform(-w, w) for j in range(L)]
    for j in range(L):
        params_dict["n" + str(j)] = disorder[j]
    H = H_dict.tohamiltonian(params_dict)
    H_mat = H.toarray()
    
    # Diagonalize the disordered H
    E, V = np.linalg.eigh(H_mat)
    
    # Thermal weights from a REFERENCE (e.g. clean or specific) Hamiltonian
    # We use the Boltzmann weights based on the INITIAL (clean) Hamiltonian energies
    # But we project the thermal state onto the disordered eigenbasis for evolution
    
    # Actually, the physically meaningful approach:
    # Start from thermal state of the CLEAN Hamiltonian, evolve under DISORDERED H
    # This captures the quench from clean to disordered system at finite T
    
    return E, V

def compute_imbalance_thermal(H_dict, I_mat, beta, w, t, i_real, n_real, 
                                H_clean_E, H_clean_V):
    """
    Prepare thermal state of the clean Hamiltonian at temperature 1/beta.
    Evolve under the disordered Hamiltonian.
    Measure imbalance.
    """
    ti = time()
    dim = len(H_clean_E)
    
    # Build disordered Hamiltonian
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    H_mat = H.toarray()
    
    # Diagonalize disordered H
    E_dis, V_dis = np.linalg.eigh(H_mat)
    
    # Thermal density matrix of clean system (in clean eigenbasis)
    boltzmann = np.exp(-beta * H_clean_E)
    Z = boltzmann.sum()
    p_n = boltzmann / Z  # thermal weights
    
    # rho_0 in clean eigenbasis: diagonal with p_n
    # Transform to site basis: rho_0 = V_clean @ diag(p_n) @ V_clean^dag
    # Then express in disordered eigenbasis:
    #   rho_mn^(dis) = <m_dis|rho_0|n_dis>
    
    # Overlap matrix: S_mn = <m_dis | n_clean>
    S = V_dis.T.conj() @ H_clean_V  # shape (dim, dim)
    
    # rho in disordered eigenbasis:
    # rho_mn = sum_k p_k * S_{m,k} * S_{n,k}^*
    rho_dis = (S * p_n[np.newaxis, :]) @ S.T.conj()  # (dim, dim)
    
    # I operator in disordered eigenbasis
    I_dis = V_dis.T.conj() @ I_mat @ V_dis  # (dim, dim)
    
    # Time evolution: rho_mn(t) = rho_mn(0) * exp(-i(E_m - E_n)t)
    # I(t) = Tr[rho(t) I] = sum_{m,n} rho_mn(0) * I_nm * exp(-i(E_m-E_n)t)
    I_t = np.zeros(len(t))
    
    # Efficient vectorized computation
    dE = E_dis[:, np.newaxis] - E_dis[np.newaxis, :]  # E_m - E_n
    A = rho_dis * I_dis.T   # element-wise: rho_mn * I_nm, shape (dim, dim)
    
    for it, tt in enumerate(t):
        phase = np.exp(-1j * dE * tt)
        I_t[it] = np.sum(A * phase).real
    
    elapsed = time() - ti
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
    return I_t

#####################################################################
# Build clean Hamiltonian and diagonalize once
#####################################################################
print("Diagonalizing clean Hamiltonian...")
params_clean = dict(H0=1.0)
for j in range(L):
    params_clean["n" + str(j)] = 0.0  # no disorder
H_clean = H_dict_template.tohamiltonian(params_clean)
H_clean_mat = H_clean.toarray()
H_clean_E, H_clean_V = np.linalg.eigh(H_clean_mat)

#####################################################################
# Also compute pure-state result for comparison
#####################################################################
print("\n=== Pure state (beta -> inf, Fock state) ===")
from quspin.operators import exp_op
from quspin.tools.measurements import obs_vs_time

def run_pure_realization(H_dict, I_op_quspin, psi_0, w, t, i_real, n_real):
    ti = time()
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict.tohamiltonian(params_dict)
    U_op = exp_op(H, a=-1j, start=t.min(), stop=t.max(),
                  num=len(t), iterate=True)
    psi_t = U_op.dot(psi_0)
    t_grid = U_op.grid
    obs_t = obs_vs_time(psi_t, t_grid, dict(I=I_op_quspin))
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {time()-ti:.2f} s")
    return obs_t["I"].real

I_data_pure = np.vstack([
    run_pure_realization(H_dict_template, I_op, psi_0_pure, w, t, i, n_real)
    for i in range(n_real)
])
I_avg_pure = I_data_pure.mean(axis=0)
boot_gen = (I_data_pure[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
sq_gen = ((b - I_avg_pure)**2 for b in boot_gen)
I_err_pure = np.sqrt(sum(sq_gen) / n_boot)

#####################################################################
# Main loop over temperatures
#####################################################################
results_thermal = {}

for beta in beta_list:
    print(f"\n=== beta*J = {beta:.1f} (T/J = {1.0/beta:.2f}) ===")
    t0 = time()
    I_data = np.vstack([
        compute_imbalance_thermal(H_dict_template, I_mat, beta, w, t, i, n_real,
                                   H_clean_E, H_clean_V)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")

    I_avg = I_data.mean(axis=0)
    boot_gen = (I_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
    sq_gen = ((b - I_avg)**2 for b in boot_gen)
    I_error = np.sqrt(sum(sq_gen) / n_boot)

    results_thermal[beta] = (I_avg, I_error)

#####################################################################
# Plot
#####################################################################
fig, ax = plt.subplots(figsize=(9, 6))

# pure state
ax.errorbar(t, I_avg_pure, I_err_pure, marker=".", markersize=3,
            color='#1f77b4', label=r"Pure Fock state ($\beta\to\infty$)",
            linewidth=0.8, elinewidth=0.4, capsize=1)

# thermal states
colors_th = ['#ff7f0e', '#2ca02c', '#d62728']
for idx, beta in enumerate(beta_list):
    I_avg, I_error = results_thermal[beta]
    ax.errorbar(t, I_avg, I_error, marker=".", markersize=3,
                color=colors_th[idx],
                label=rf"$\beta J = {beta:.1f}$ (T/J={1.0/beta:.2f})",
                linewidth=0.8, elinewidth=0.4, capsize=1)

ax.set_xlabel(r"$Jt$", fontsize=18)
ax.set_ylabel(r"$\mathcal{I}$", fontsize=18)
ax.set_title(f"MBL: Temperature Effect (w={w:.1f}, U={U_val:.1f}, L={L})", fontsize=16)
ax.tick_params(labelsize=14)
ax.legend(loc="best", fontsize=11)
ax.set_xlim(0, 35)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("MBL_temperature_effect.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_temperature_effect.png")
plt.close(fig)
