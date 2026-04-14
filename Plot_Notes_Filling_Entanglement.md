# Physics Interpretation: Entanglement Entropy Growth vs Filling Fraction

**Target Plot:** `MBL_filling_entanglement.png`
**Observable:** Half-chain Von Neumann Entropy $S(t)$
**Fixed Parameters:** $U=5.0, W=10.0, L=8$
**Variable:** Total particle number $N = 2, 4, 6, 8$
**Core Literature References:** Bardarson et al. (PRL 2012), De Luca & Scardicchio (EPL 2013)

---

### 1. The Core Question Being Answered
Scripts 01–06 establish that MBL exists in this model at a fixed filling. This script asks: **how does the strength of MBL entanglement growth change as we add more particles?**

This is a filling-fraction scan at fixed $U$ and $W$. Crucially, since $W=10J$ is deep in the MBL phase for all fillings, every curve should show logarithmic growth — the question is *how fast* the growth is (the prefactor $\alpha$).

### 2. Expected Behavior for Each Filling

*   **$N=2$ ($N_{up}=1, N_{down}=1$, dim=64):**
    Only two particles on 8 sites. The probability that the two particles even occupy the same site — and thus interact via $U$ — is extremely low. The system is effectively non-interacting. Entanglement saturates rapidly at a small value, closely resembling the $U=0$ Anderson-localized case. The logarithmic growth coefficient $\alpha \approx 0$ or very small.

*   **$N=4$ ($N_{up}=2, N_{down}=2$, dim=784, quarter-filling):**
    This is the reference system from Scripts 01–06. Interactions are present but particles still have many empty sites to avoid. A clear logarithmic growth $S(t) \sim \alpha \ln(Jt)$ is observed with a moderate prefactor. This is the canonical quarter-filling MBL signature.

*   **$N=6$ ($N_{up}=3, N_{down}=3$, dim=3136):**
    With 6 of 8 sites occupied, double-occupancy (triggering $U$) is frequent. The inter-l-bit coupling is stronger, dephasing happens faster, and the logarithmic growth prefactor $\alpha$ is larger. The slope on the semi-log plot is visibly steeper than for $N=4$.

*   **$N=8$ ($N_{up}=4, N_{down}=4$, dim=4900, half-filling):**
    Every site has exactly one particle on average. The Hubbard repulsion $U$ is maximally effective. The dephasing rate between localized l-bits is at its peak. The logarithmic growth prefactor $\alpha$ reaches its maximum value. Additionally, at half-filling the Mott gap opens, which can subtly compete with disorder to modify the effective $W_c$.

### 3. Why This Plot Is Consistent with Literature
The filling-dependence of entanglement growth is a well-established prediction of the l-bit picture of MBL (Huse, Nandkishore & Oganesyan, PRB 2014). The key physical prediction is:
$$\alpha(N) \propto J_{\text{eff}}(N) \sim \bar{V} \cdot \bar{n}^2$$
where $\bar{n} = N/(2L)$ is the average density. More particles means longer-range effective l-bit interactions and faster dephasing, which directly increases $\alpha$.

The fact that all curves plateau well below the ergodic Page value — and show logarithmic rather than linear growth — confirms the system remains in the MBL phase for all fillings at $W=10J$, consistent with the phase diagram of Mondaini & Rigol (PRB 2015).

### 4. Initial State Note
For each filling, particles are placed to maximize sublattice contrast: up-spins on even sites, down-spins on odd sites. This gives a well-defined, filling-independent starting point with maximum spatial ordering. The entanglement entropy starts at $S(0) = 0$ for all fillings (product state = zero entanglement), making all curves directly comparable.
