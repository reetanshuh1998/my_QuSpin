# Physics Interpretation: MBL Level Spacing Statistics

**Target Plot:** `MBL_level_spacing.png`
**Observable:** Mean Level Spacing Ratio $\langle r \rangle$
**Hamiltonian:** 1D Disordered Fermi-Hubbard Model ($U=5.0, J=1.0$)
**Core Literature Reference:** Oganesyan & Huse (PRB 2007), Atas et al. (PRL 2013), Luitz et al. (PRB 2015)

---

### 1. What the Level Spacing Ratio Measures
While observables like imbalance and entanglement entropy probe **dynamics** (how states evolve in time), the level spacing ratio $\langle r \rangle$ probes **spectral statistics** — the static structure of the energy spectrum itself.

For each pair of consecutive energy gaps $\delta_n = E_{n+1} - E_n$ and $\delta_{n+1} = E_{n+2} - E_{n+1}$, the ratio is defined as:
$$r_n = \frac{\min(\delta_n, \delta_{n+1})}{\max(\delta_n, \delta_{n+1})}$$

This ratio is bounded between 0 and 1, and its mean $\langle r \rangle$ reveals crucial information about the nature of the quantum states:

### 2. The Two Universal Limits

*   **GOE (Gaussian Orthogonal Ensemble), $\langle r \rangle \approx 0.53$:**
    In the ergodic/thermal phase, energy levels exhibit **level repulsion** — nearby levels actively push each other apart. This is the hallmark of quantum chaos and the Eigenstate Thermalization Hypothesis (ETH). The energy spectrum looks like the eigenvalues of a random matrix (Wigner-Dyson statistics).

*   **Poisson, $\langle r \rangle \approx 0.386$:**
    In the MBL/localized phase, eigenstates are localized in different regions of Hilbert space and do not "see" each other. Their energies are uncorrelated, like random numbers drawn independently from a distribution. There is no level repulsion — levels can cluster arbitrarily close together.

### 3. Reading the Plot: The MBL Phase Transition
The plot shows $\langle r \rangle$ as a function of disorder strength $W/J$:

*   **Weak disorder ($W \lesssim 2$):** The system is in the ergodic phase. The spectrum exhibits GOE statistics ($\langle r \rangle \approx 0.53$), confirming that eigenstates are delocalized and thermalize.
*   **Strong disorder ($W \gtrsim 6$):** The system is in the MBL phase. The spectrum shows Poisson statistics ($\langle r \rangle \approx 0.386$), confirming eigenstate localization.
*   **Critical region ($W \approx 3$–$5$):** The crossover between GOE and Poisson reveals the **MBL phase transition**. The curve smoothly interpolates between the two limits, with the sharpness of the transition increasing with system size $L$.

### 4. Why We Filter to the Middle of the Spectrum
The code only analyzes the middle 50% of the eigenvalue spectrum (discarding the top and bottom 25%). This is critical because:
*   The band edges always exhibit anomalous statistics (the density of states drops to zero, creating artificial level clustering).
*   The MBL transition is an energy-dependent phenomenon — the middle of the spectrum is where the transition is sharpest and most representative of the bulk behavior.

### 5. Finite-Size Effects
At $L=8$ (Hilbert space dimension 784), the transition curve is **broadened** compared to the thermodynamic limit. In larger systems ($L \geq 14$), the crossover would be much sharper, approaching a true phase transition. The broadened crossover at $L=8$ is expected and consistent with all published exact diagonalization studies at comparable sizes. The qualitative S-shaped interpolation from GOE to Poisson is the key physical signature, not the exact location of the crossing point.

### 6. Why $\langle r \rangle$ is Preferred Over $P(s)$
Older studies used the full level spacing distribution $P(s)$ (requiring "spectrum unfolding" to remove the local density of states). The ratio $r_n$ was introduced by Oganesyan & Huse (2007) specifically because it is **independent of the local density of states** and requires no unfolding. This makes it computationally simpler and statistically more robust — which is why it has become the standard MBL diagnostic in modern literature.
