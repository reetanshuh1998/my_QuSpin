"""
Phase 1: Reproduce Figure 3 from SciPost Phys. 7, 020 (2019)
--------------------------------------------------------------
Many-Body Localization in the 1D disordered Fermi-Hubbard model.

Measures the sublattice imbalance I(t) = (N_A - N_B) / N_tot
averaged over 100 disorder realizations for disorder strengths
w = 1.0, 4.0, 10.0.

Parameters:
    L = 8       (chain length)
    J = 1.0     (hopping)
    U = 5.0     (on-site interaction)
    N_up = 2, N_down = 2  (quarter filling)
    OBC (open boundary conditions)

Reference: Section 2.3 of the QuSpin paper.
"""

from quspin.operators import hamiltonian, exp_op, quantum_operator
from quspin.basis import spinful_fermion_basis_1d
from quspin.tools.measurements import obs_vs_time
import numpy as np
from numpy.random import uniform, choice
from time import time
import matplotlib.pyplot as plt

np.random.seed(42)

#####################################################################
# Parameters
#####################################################################
n_real = 100   # number of disorder realizations
n_boot = 100   # number of bootstrap samples for error estimation

L = 8          # system size
N = L // 2     # total number of particles
N_up = N // 2 + N % 2   # spin-up fermions
N_down = N // 2          # spin-down fermions

w_list = [1.0, 4.0, 10.0]  # disorder strengths
J = 1.0    # hopping strength
U = 5.0    # interaction strength

# time grid
start, stop, num = 0.0, 35.0, 101
t = np.linspace(start, stop, num=num, endpoint=True)

#####################################################################
# Build basis
#####################################################################
print("Building basis...")
basis = spinful_fermion_basis_1d(L, Nf=(N_up, N_down))
print(f"Hilbert space dimension: {basis.Ns}")

#####################################################################
# Build Hamiltonian structure (using quantum_operator for disorder)
#####################################################################
# hopping (OBC)
hop_right = [[-J, i, i + 1] for i in range(L - 1)]
hop_left  = [[ J, i, i + 1] for i in range(L - 1)]
# on-site interaction
int_list = [[U, i, i] for i in range(L)]
# sublattice imbalance observable
sublat_list = [[(-1.0)**i / N, i] for i in range(L)]

# operator lists for the clean Hamiltonian
operator_list_0 = [
    ["+-|", hop_left ],   # up hop left
    ["-+|", hop_right],   # up hop right
    ["|+-", hop_left ],   # down hop left
    ["|-+", hop_right],   # down hop right
    ["n|n", int_list ],   # on-site interaction
]
imbalance_list = [["n|", sublat_list], ["|n", sublat_list]]

# build operator dictionary: H0 + site-dependent disorder potentials
operator_dict = dict(H0=operator_list_0)
for i in range(L):
    operator_dict["n" + str(i)] = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]

#####################################################################
# Build operators
#####################################################################
no_checks = dict(check_pcon=False, check_symm=False, check_herm=False)
H_dict = quantum_operator(operator_dict, basis=basis, **no_checks)
I_op   = hamiltonian(imbalance_list, [], basis=basis, **no_checks)

# initial state: alternating spins on every other site
s_up   = "".join("1000" for _ in range(N_up))
s_down = "".join("0010" for _ in range(N_down))
i_0    = basis.index(s_up, s_down)
psi_0  = np.zeros(basis.Ns)
psi_0[i_0] = 1.0

# Sanity Check
assert abs(I_op.expt_value(psi_0).real - 1.0) < 1e-10, "Initial imbalance must be 1.0"

print(f"Initial state: |{s_up}>(x)|{s_down}>")

#####################################################################
# Function for a single disorder realization
#####################################################################
def run_realization(H_dict, I_op, psi_0, w, t, i_real):
    ti = time()
    # random on-site potential
    params_dict = dict(H0=1.0)
    for j in range(L):
        params_dict["n" + str(j)] = uniform(-w, w)
    # construct Hamiltonian for this realization
    H = H_dict.tohamiltonian(params_dict)
    # time evolution via matrix exponential
    U_op = exp_op(H, a=-1j, start=t.min(), stop=t.max(),
                  num=len(t), iterate=True)
    psi_t = U_op.dot(psi_0)
    # measure imbalance
    t_grid = U_op.grid
    obs_t = obs_vs_time(psi_t, t_grid, dict(I=I_op))
    print(f"  realization {i_real+1}/{n_real} completed in {time()-ti:.2f} s")
    return obs_t["I"].real

#####################################################################
# Main loop over disorder strengths
#####################################################################
fig, ax = plt.subplots(figsize=(8, 5))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for idx, w in enumerate(w_list):
    print(f"\n--- Disorder strength w = {w:.2f} ---")
    t0 = time()
    I_data = np.vstack([
        run_realization(H_dict, I_op, psi_0, w, t, i)
        for i in range(n_real)
    ])
    print(f"All realizations for w={w:.2f} completed in {time()-t0:.1f} s")

    # disorder average
    I_avg = I_data.mean(axis=0)

    # bootstrap error estimation
    bootstrap_gen = (
        I_data[choice(n_real, size=n_real)].mean(axis=0)
        for _ in range(n_boot)
    )
    sq_fluc_gen = ((bootstrap - I_avg)**2 for bootstrap in bootstrap_gen)
    I_error = np.sqrt(sum(sq_fluc_gen) / n_boot)

    # plot
    ax.errorbar(t, I_avg, I_error, marker=".", markersize=3,
                color=colors[idx], label=f"w={w:.2f}", linewidth=0.8,
                elinewidth=0.5, capsize=1)

#####################################################################
# Formatting
#####################################################################
ax.set_xlabel(r"$Jt$", fontsize=18)
ax.set_ylabel(r"$\mathcal{I}$", fontsize=18)
ax.tick_params(labelsize=14)
ax.legend(loc="best", fontsize=14)
ax.set_xlim(0, 35)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("fermion_MBL_fig3.png", dpi=150, bbox_inches="tight")
print("\nSaved: fermion_MBL_fig3.png")
plt.close()
