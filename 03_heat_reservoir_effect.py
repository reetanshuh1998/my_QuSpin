"""
Phase 3: Effect of Heat Reservoir on MBL
-----------------------------------------
Study how coupling to an external heat bath destroys
Many-Body Localization in the 1D disordered Fermi-Hubbard model.

We model dephasing via the Lindblad master equation:
  d rho / dt = -i[H, rho] + gamma * sum_k (n_k rho n_k - 1/2 {n_k, rho})

Because n_k (the number operator) is diagonal in the real-space
Fock basis (the basis QuSpin uses by default), the dissipative part
is super-simple in this basis:
  Dissipator(rho)_{mn} = -Gamma_{mn} * rho_{mn}
where Gamma_{mn} = (gamma/2) * sum_k (m_k - n_k)^2.
Here m_k and n_k are the occupation numbers of state m and n on site/spin k.
(m_k - n_k)^2 is exactly the Hamming distance between the two Fock states!

We can thus simulate the full Lindblad equation in the site basis using
an ODE solver without needing any large 4D tensors or superoperators,
greatly reducing memory and computation time.

Parameters:
    L = 6, J = 1.0, U = 5.0, w = 10.0
    n_real = 30
"""

from quspin.operators import hamiltonian, quantum_operator
from quspin.basis import spinful_fermion_basis_1d
import numpy as np
from numpy.random import uniform, choice
from scipy.integrate import solve_ivp
from time import time
import matplotlib.pyplot as plt

np.random.seed(42)

#####################################################################
# Parameters
#####################################################################
n_real = 20   # Reduced: density matrix at half-filling (dim=400) is expensive
n_boot = 50

L = 6
N = L              # half-filling
N_up = L // 2      # = 3
N_down = L // 2    # = 3

w = 5.0
J_hop = 1.0
U_val = 5.0
gamma_list = [0.0, 0.1, 0.5, 1.0]

start, stop, num = 0.0, 35.0, 71
t_eval = np.linspace(start, stop, num=num, endpoint=True)

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
dim = basis.Ns
print(f"Hilbert space dimension: {dim}")

#####################################################################
# Initial state
#####################################################################
# CDW Néel state: all particles on even sites (doubly occupied)
s_up   = "101010"   # up-spins at sites 0, 2, 4
s_down = "101010"   # down-spins at sites 0, 2, 4
i_0    = basis.index(s_up, s_down)
psi_0  = np.zeros(dim)
psi_0[i_0] = 1.0
rho_0  = np.outer(psi_0, psi_0).astype(np.complex128)

#####################################################################
# Precompute Dissipator Decay Matrix Gamma_{mn}
#####################################################################
print("Precomputing dissipator matrix...")
no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
jump_diags = []
# Collect diagonals of all jump operators (n_i,up and n_i,down)
for i in range(L):
    n_up = hamiltonian([["n|", [[1.0, i]]]], [], basis=basis, **no_checks).toarray().diagonal()
    n_dn = hamiltonian([["|n", [[1.0, i]]]], [], basis=basis, **no_checks).toarray().diagonal()
    jump_diags.append(n_up)
    jump_diags.append(n_dn)

jump_diags = np.array(jump_diags)  # shape (2L, dim)
# Gamma_{mn} = 1/2 * sum_k (m_k - n_k)^2
# Using broadcasting:
m_k = jump_diags[:, :, np.newaxis]  # (2L, dim, 1)
n_k = jump_diags[:, np.newaxis, :]  # (2L, 1, dim)
Gamma_matrix = 0.5 * np.sum((m_k - n_k)**2, axis=0)  # shape (dim, dim)

#####################################################################
# Site-coupling lists for Hamiltonian
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

I_op = hamiltonian(imbalance_list, [], basis=basis, **no_checks)
I_mat = I_op.toarray()

# Sanity assertion
assert abs(I_op.expt_value(psi_0).real - 1.0) < 1e-10, "Initial imbalance must be 1.0"

#####################################################################
# Lindblad RHS function using split real/imag to please solve_ivp
#####################################################################
def lindblad_rhs(t_val, y, H_mat, Gamma_matrix_gamma, dim):
    n2 = dim**2
    rho = (y[:n2] + 1j * y[n2:]).reshape(dim, dim)
    
    # Commutator
    drho = -1j * (H_mat @ rho - rho @ H_mat)
    # Dissipator
    drho -= Gamma_matrix_gamma * rho
    
    drho_flat = drho.flatten()
    return np.concatenate([drho_flat.real, drho_flat.imag])

#####################################################################
# Realization function
#####################################################################
def run_realization(H_dict_obj, gamma, w, i_real, n_real):
    ti = time()
    
    # Disordered Hamiltonian
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    H = H_dict_obj.tohamiltonian(params_dict)
    H_mat = H.toarray()
    
    if gamma == 0:
        # Unitary evolution
        E, V = np.linalg.eigh(H_mat)
        rho_eig = V.T.conj() @ rho_0 @ V
        I_eig = V.T.conj() @ I_mat @ V
        dE = E[:, np.newaxis] - E[np.newaxis, :]
        A = rho_eig * I_eig.T
        
        I_t = np.zeros(len(t_eval))
        for it, tt in enumerate(t_eval):
            I_t[it] = np.sum(A * np.exp(-1j * dE * tt)).real
    else:
        # Lindblad evolution
        Gamma_matrix_gamma = gamma * Gamma_matrix
        rho0_flat = rho_0.flatten()
        y0 = np.concatenate([rho0_flat.real, rho0_flat.imag])
        
        sol = solve_ivp(
            lindblad_rhs,
            [t_eval[0], t_eval[-1]],
            y0,
            t_eval=t_eval,
            args=(H_mat, Gamma_matrix_gamma, dim),
            method='RK45',
            rtol=1e-6,
            atol=1e-8,
        )
        
        n2 = dim**2
        I_t = np.zeros(len(t_eval))
        for it in range(len(t_eval)):
            rho_t = (sol.y[:n2, it] + 1j * sol.y[n2:, it]).reshape(dim, dim)
            I_t[it] = np.trace(rho_t @ I_mat).real
            
    elapsed = time() - ti
    if (i_real + 1) % 10 == 0 or i_real == 0:
        print(f"    realization {i_real+1}/{n_real} done in {elapsed:.2f} s")
    return I_t

#####################################################################
# Main loop
#####################################################################
results = {}

for gamma in gamma_list:
    print(f"\n=== gamma = {gamma:.3f} ===")
    t0 = time()
    
    H_dict_obj = quantum_operator(operator_dict, basis=basis, **no_checks)
    I_data = np.vstack([
        run_realization(H_dict_obj, gamma, w, i, n_real)
        for i in range(n_real)
    ])
    print(f"  Completed in {time()-t0:.1f} s")
    
    I_avg = I_data.mean(axis=0)
    boot_gen = (I_data[choice(n_real, size=n_real)].mean(axis=0) for _ in range(n_boot))
    sq_gen = ((b - I_avg)**2 for b in boot_gen)
    I_error = np.sqrt(sum(sq_gen) / n_boot)
    
    results[gamma] = (I_avg, I_error)

#####################################################################
# Plot
#####################################################################
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

for idx, gamma in enumerate(gamma_list):
    I_avg, I_error = results[gamma]
    label = rf"$\gamma = {gamma:.2g}$" if gamma > 0 else r"$\gamma = 0$ (isolated)"
    
    # Linear scale subplot
    ax1.errorbar(t_eval, I_avg, I_error, marker=".", markersize=3,
                 color=colors[idx], label=label,
                 linewidth=0.8, elinewidth=0.4, capsize=1)
    
    # Log-log scale subplot (skip t=0)
    ax2.errorbar(t_eval[1:], I_avg[1:], I_error[1:], marker=".", markersize=3,
                 color=colors[idx], label=label,
                 linewidth=0.8, elinewidth=0.4, capsize=1)

ax1.set_xlabel(r"$Jt$", fontsize=16)
ax1.set_ylabel(r"$\mathcal{I}$", fontsize=16)
ax1.set_title("Linear Decay", fontsize=14)
ax1.tick_params(labelsize=12)
ax1.legend(loc="best", fontsize=11)
ax1.set_xlim(0, 35)
ax1.grid(True, alpha=0.3)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r"$Jt$ (log)", fontsize=16)
ax2.set_ylabel(r"$\mathcal{I}$ (log)", fontsize=16)
ax2.set_title("Subdiffusive Power-Law Decay", fontsize=14)
ax2.tick_params(labelsize=12)
ax2.legend(loc="best", fontsize=11)
ax2.grid(True, alpha=0.3)

fig.suptitle(f"MBL: Heat Reservoir Effect (w={w:.1f}, U={U_val:.1f}, L={L})", fontsize=16)
fig.tight_layout()
fig.savefig("MBL_heat_reservoir.png", dpi=150, bbox_inches="tight")
print("\nSaved: MBL_heat_reservoir.png")
plt.close(fig)
