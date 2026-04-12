# Physics Interpretation: MBL Effective Temperature (Energy-Shell)

**Target Plot:** `MBL_temperature_effect_shell.png`
**Observable:** Sublattice Imbalance $\mathcal{I}(t)$
**Hamiltonian:** 1D Disordered Fermi-Hubbard Model ($W=4.0, U=5.0, J=1.0$)

---

### 1. The Necessity of the Microcanonical Energy Shell
In Many-Body Localization (MBL) simulations, directly simulating a classical "Thermal Density Matrix" ($\rho \propto e^{-\beta H}$) using the clean system collapses the initial imbalance to $\mathcal{I}(0) \approx 0$ due to translational symmetry. If $\mathcal{I}(0) = 0$, there is no sublattice order for the MBL system to preserve, rendering the time-evolution flat and physically uninformative.

To bypass this and correctly study temperature effects, this plot utilizes **Microcanonical Energy-Shell Sampling**. 
1. We start with the strict **Néel State** (where all particles occupy even sites), guaranteeing a maximal initial imbalance of $\mathcal{I}(0) = 1.0$.
2. We project this Néel state directly onto narrow energy "tiers" (shells) of the disordered Hamiltonian's eigenspectrum. 
3. By tuning the normalized energy density target ($\epsilon$), we can emulate different "effective temperatures" while retaining the necessary $\mathcal{I}(0) \approx 1$ starting condition.

### 2. Physical Meaning of the Energy Densities ($\epsilon$)
The plot separates the time-evolution of the imbalance into three distinct thermal shells:

*   **$\epsilon = 0.1$ (Low-Temperature proxy):** 
    This shell isolates eigenstates residing in the bottom $10\%$ of the energy spectrum. Because the density of states is low and energy is close to the ground state, it maps to a "cold" thermal limit.
*   **$\epsilon = 0.5$ (Infinite-Temperature proxy):** 
    This isolates the dead-center of the energy spectrum. In bounded quantum systems, the mid-spectrum contains the highest density of states, heavily packed with highly excited, chaotic states. This corresponds strictly to $T \to \infty$.
*   **$\epsilon = 0.9$ (Negative-Temperature proxy):** 
    This isolates the top $10\%$ of the spectrum. Because the density of states is strictly bounded and actively decreasing at the upper limit, the thermodynamic relation $1/T = \partial S / \partial E$ turns negative, marking highly-energetic "Negative Temperature" states.

### 3. Proof of the Many-Body Mobility Edge
The behavior shown in the plot correctly captures a fundamental element of MBL physics: **The Mobility Edge**. 

In the deep MBL phase ($W=4.0$), all three shells remain localized (none of the curves decay down to zero). However, the *saturation height* (the plateau $\mathcal{I}_{\infty}$) is not identical:
*   The mid-spectrum states ($\epsilon=0.5$) are the most fundamentally energetic and chaotic. They are the hardest to localize, meaning they experience the most thermal "leakage" before freezing. Consequently, the infinite-temperature curve settles at the lowest plateau.
*   The band-edge states ($\epsilon=0.1$ and $\epsilon=0.9$) have a much shorter localization length. They freeze rapidly and securely under the protective random field, resulting in significantly higher plateau saturation.

### Conclusion
This plot correctly demonstrates that while strong disorder ($W=4.0$) can localize the entire spectrum (resulting in non-zero $\mathcal{I}_{\infty}$ plateaus across the board), **higher effective temperatures weaken the localization effect**, explicitly lowering the survival plateau of the initial Néel state due to increased internal ergodic pressure.
