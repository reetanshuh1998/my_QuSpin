# Physics Interpretation: MBL Entanglement Entropy Growth

**Target Plot:** `MBL_entanglement_growth_U_sweep.png`
**Observable:** Half-Chain Von Neumann Entanglement Entropy $S(t)$
**Hamiltonian:** 1D Disordered Fermi-Hubbard Model ($W=10.0, J=1.0$)
**Core Literature Reference:** Bardarson, Pollmann & Moore (PRL 109, 017202, 2012), Žnidarič et al. (PRB 2008)

---

### 1. Why Entanglement Entropy is the Gold Standard for MBL
The Sublattice Imbalance $\mathcal{I}(t)$ (used in earlier scripts) measures *local* memory — whether particles stay on odd or even sites. Entanglement Entropy $S(t)$ probes something far deeper: it measures how much *quantum information* has leaked from one half of the chain to the other.

In an **ergodic** (thermalizing) system, entanglement spreads ballistically: $S(t) \sim t$, rapidly reaching the thermodynamic maximum ("Page value"). In contrast, in an **MBL** system, entanglement spreads incredibly slowly — *logarithmically*: $S(t) \sim \alpha \ln(t)$. This unbounded but ultra-slow growth is one of the most striking and unique signatures of Many-Body Localization, distinguishing it sharply from simple Anderson Localization.

### 2. The Three Regimes Shown in the Plot

*   **$U = 0$ (Anderson Localization — Navy Curve):** 
    Without interactions, the system is in the Anderson Localized phase. The initially unentangled product (Fock) state evolves, but because single-particle orbitals are exponentially localized, entanglement can only spread within a localization length $\xi$. Once the entanglement frontier reaches $\xi$ (which happens very quickly), $S(t)$ **saturates** at a small, finite value. There is no logarithmic growth — the entanglement simply freezes.

*   **$U = 1.0$ (Weak MBL — Green Curve):**
    Turning on weak interactions introduces dephasing between the localized orbitals (called "l-bits" in the MBL literature). This dephasing slowly generates entanglement across the chain. The result is the hallmark **logarithmic growth**: $S(t) \sim \alpha(U) \ln(Jt)$, with a small prefactor $\alpha$. The growth is unbounded in principle (for infinite systems) but extremely slow — it takes exponentially long times to entangle distant regions.

*   **$U = 5.0, 10.0$ (Stronger MBL — Orange, Red Curves):**
    Stronger interactions produce faster dephasing between the l-bits. The logarithmic growth persists, but with a **larger prefactor** $\alpha(U)$: the curves rise more steeply on the semi-log plot. Despite this faster entanglement growth, the system remains localized — the growth is still logarithmic, not linear. All curves plateau well below the ergodic Page value, confirming the area-law character of MBL eigenstates.

### 3. The Physical Origin of Logarithmic Growth
The logarithmic entanglement growth is explained by the **l-bit (localized bit) picture** of MBL (Huse, Nandkishore & Oganesyan, 2014):
*   In the MBL phase, the true eigenstates of the Hamiltonian can be described by a set of quasi-local integrals of motion (l-bits): $\tau_i^z$.
*   These l-bits interact via exponentially decaying couplings: $J_{ij} \sim e^{-|i-j|/\xi}$.
*   When we start from a product state (which is not an eigenstate of the l-bits), the weak inter-l-bit interactions slowly generate entanglement.
*   The time required to entangle l-bits separated by distance $r$ scales as $t \sim e^{r/\xi} / J_{ij}$, meaning distance grows as $r \sim \xi \ln(t)$.
*   Since entanglement entropy scales with the number of entangled l-bits, $S(t) \sim \ln(t)$.

### 4. Why the Log-Spaced Time Axis is Essential
Logarithmic growth ($S \sim \ln t$) appears as a **straight line** only when the x-axis is logarithmic. On a linear time axis, the growth would look deceptively flat and indistinguishable from saturation. The semi-log plot spanning $Jt \in [0.1, 1000]$ (four decades) is therefore essential to resolve and confirm the logarithmic law.

### 5. Finite-Size Saturation
At late times ($Jt \gtrsim 100$–$1000$), all curves eventually plateau. This is **not** a sign of thermalization — it is a finite-size effect. For $L=8$, the maximum possible entanglement entropy of the left half-chain is bounded by $S_{max} = \ln(\dim_A)$. Once this ceiling is reached, the entropy cannot grow further regardless of the dynamics. In the thermodynamic limit ($L \to \infty$), the logarithmic growth would continue indefinitely.

### 6. Comparison with Literature
The qualitative features of this plot — logarithmic growth for $U > 0$, rapid saturation for $U = 0$, increasing prefactor with $U$ — are fully consistent with the foundational results of Bardarson et al. (PRL 2012). The original paper studied spinless fermions at half-filling, while this simulation uses spinful fermions at quarter-filling. The qualitative signatures are model-independent; only the numerical values of the prefactors and saturation levels differ.
