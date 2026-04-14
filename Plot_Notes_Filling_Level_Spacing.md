# Physics Interpretation: Level Spacing Statistics vs Filling Fraction

**Target Plot:** `MBL_filling_level_spacing.png`
**Observable:** Mean Level Spacing Ratio $\langle r \rangle$
**Fixed Parameters:** $U=5.0, W=10.0, L=8$
**Variable:** Total particle number $N = 2, 4, 6, 8$
**Core Literature References:** Oganesyan & Huse (PRB 2007), Atas et al. (PRL 2013)

---

### 1. What This Plot Shows
This is the **spectral** complement to Script 07. Instead of measuring dynamics ($S(t)$), it probes the static energy spectrum at each filling. The mean level spacing ratio $\langle r \rangle$ acts as a filling-independent order parameter:
- $\langle r \rangle \approx 0.53$ (GOE) → ergodic/thermal phase
- $\langle r \rangle \approx 0.386$ (Poisson) → MBL/localized phase

With $W=10J$ fixed deep in the localized regime, all fillings should give $\langle r \rangle$ close to the Poisson value. The key physics is how *quickly* the system localizes (converges to Poisson) as particle number increases.

### 2. Expected Behavior Across Fillings

*   **$N=2$ (dim=64):**
    Few particles, rare interactions. The spectrum is dominated by single-particle Anderson localization, which is Poisson for any non-zero disorder in 1D. $\langle r \rangle$ should be closest to 0.386.

*   **$N=4$ (dim=784, quarter-filling):**
    Many-body interactions begin to play a role. At $W=10J$, the system is still deeply localized, so $\langle r \rangle$ remains near Poisson. This is the crosscheck of the quarter-filling level spacing already seen in Script 06.

*   **$N=6$ (dim=3136):**
    More interaction channels open up, enhancing many-body complexity. The level statistics may start to shift slightly toward GOE, as more inter-l-bit connections are active. However at $W=10J$, the system should remain in the MBL phase.

*   **$N=8$ (dim=4900, half-filling):**
    Maximum complexity. At half-filling, the Mott physics and disorder compete. Even here at $W=10J$, localization should dominate. However, this point is most susceptible to finite-size effects.

### 3. Why This Plot Is Consistent with Literature
At the fixed deep-MBL disorder $W=10J$ (well above the critical $W_c \approx 3.5$–$5J$), all fillings should show Poisson-like statistics, confirming the MBL phase is robust across filling fractions. The gradual shift of $\langle r \rangle$ toward GOE with increasing filling demonstrates the general principle: **more interactions → closer to the phase boundary** (even if not crossing it at $W=10J$).

This filling dependence of spectral statistics is consistent with the general many-body localization phase diagram discussed in Nandkishore & Huse (Annual Reviews 2015) and Altman & Vosk (Annual Reviews 2015).

### 4. Practical Value of This Plot
This plot serves as a **self-consistency check** for Script 07: if Script 07 shows logarithmic entanglement growth for all fillings (consistent with MBL), Script 08 should show all fillings giving Poisson statistics (also consistent with MBL). Agreement between both diagnostics — dynamic and spectral — provides strong evidence that the system is genuinely in the MBL phase across all probed filling fractions.
