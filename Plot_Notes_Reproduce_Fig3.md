# Physics Interpretation: Reproducing Fig. 3 — MBL Baseline

**Target Plot:** `fermion_MBL_fig3.png`
**Observable:** Sublattice Imbalance $\mathcal{I}(t)$
**Hamiltonian:** 1D Disordered Fermi-Hubbard Model ($U=5.0, J=1.0$)
**Core Literature Reference:** Schreiber et al. (Science 349, 842, 2015), Mondaini & Rigol (PRB 91, 2015)

---

### 1. What the Imbalance Measures
The sublattice imbalance is defined as:
$$\mathcal{I}(t) = \frac{1}{N} \sum_{i=0}^{L-1} (-1)^i \left( \langle n_{i\uparrow}(t) \rangle + \langle n_{i\downarrow}(t) \rangle \right)$$

Particles start on even sites only (Néel state), giving $\mathcal{I}(0) = 1.0$ exactly. In the ergodic (thermal) phase, particles spread uniformly and $\mathcal{I}(\infty) \to 0$. In the MBL phase, particles remain spatially localized and $\mathcal{I}(\infty)$ stays at a large positive plateau.

### 2. The Three Disorder Regimes

*   **$W = 1.0$ (Ergodic phase):**
    Weak disorder cannot prevent thermalization. Particles hop freely, collide via Hubbard-$U$ interactions, and rapidly redistribute. The imbalance decays toward zero on a timescale $Jt \sim$ a few units. This confirms eigenstate thermalization hypothesis (ETH): all local observables relax to thermal values.

*   **$W = 4.0$ (Near the MBL transition):**
    At intermediate disorder, the system sits close to the critical point $W_c \approx 3.5$–$5J$ (for $U=5J$ at this filling). The imbalance decays slowly and settles at a small but non-zero plateau. This crossover regime is the most sensitive to system parameters and shows the competition between hopping, interaction, and disorder.

*   **$W = 10.0$ (Deep MBL phase):**
    Strong disorder localizes all particles exponentially. The imbalance barely decays from its initial value of 1.0 and remains at a large plateau $\mathcal{I}_\infty \approx 0.5$–$0.8$. This is the hallmark of MBL: the system retains memory of its initial conditions indefinitely.

### 3. Why the Plot Is Consistent with Literature
The three-curve fan structure — highest plateau for strongest disorder, decay to zero for weakest disorder — directly reproduces the qualitative behavior of Figure 3 of Schreiber et al. (Science 2015), originally demonstrated in a cold-atom optical lattice experiment. The numerical protocol (Néel initial state, disorder-averaged dynamics, long-time plateau $\mathcal{I}_\infty$) is the standard benchmark in MBL numerical studies.

### 4. System Details
- **Quarter-filling** ($N_{up}=2, N_{down}=2$, $L=8$): Hilbert space dim = 784
- The imbalance normalization by $N = N_{up} + N_{down} = 4$ ensures $\mathcal{I}(0) = 1.0$ regardless of filling
- 100 disorder realizations with bootstrap error bars provide statistically reliable averages
