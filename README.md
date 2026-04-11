# MSc Project: Many-Body Localization in the Disordered Fermi-Hubbard Model

This repository contains the source code, numerical simulations, and physical analyses for an MSc physics project investigating the phenomenon of **Many-Body Localization (MBL)**. 

## Motivation & Background
The Eigenstate Thermalization Hypothesis (ETH) suggests that interacting quantum systems act as their own heat baths, eventually reaching thermal equilibrium and losing memory of their initial conditions. MBL is a robust, non-thermal phase of matter that completely violates the ETH. 
Unlike **Anderson Localization**—which occurs in non-interacting systems and can be destroyed by the slightest inter-particle interaction—Many-Body Localization persists even in strongly interacting regimes. 

This project mathematically proves the existence and properties of MBL using Exact Diagonalization (ED) on the 1D disordered Fermi-Hubbard model. Throughout six targeted experimental phases, we explore how MBL breaks down under the influence of finite temperature, open environment coupling, and extreme repulsive interactions, while reproducing the canonical spectral and entanglement signatures that physically define the phase.

## Computational Library
All simulations are built on top of **QuSpin**, an incredibly powerful Python library for exact diagonalization and quantum dynamics.
- **Official Repository:** [weinbe58/QuSpin](https://github.com/weinbe58/QuSpin)
- **Reference Example:** Our fundamental baseline architecture draws direct inspiration from QuSpin's official `example6.py`. Our project extends the underlying engine to model Lindbladians, density matrices, spectral level spacings, and entanglement growth outside the scope of standard tutorials.

## Physical Parameters Setup
- **Kinetic Hopping ($J = 1.0$):** This serves as our ruler. Time, Disorder, and Interaction are measured as unitless ratios relative to $J$.
- **Interaction Strength ($U/J = 5.0$):** Sits directly in the optimal "strongly interacting" regime required to produce genuine Many-Body Localization, guaranteeing that the localization we see is purely a Many-Body effect.
- **System Sizes ($L$):** We mathematically balance memory computing limits by using $L=8$ for $O(N)$ Pure State vector evolution, and $L=6$ for $O(N^2)$ Mixed State Density Matrix evolution.

---

## Output Plots & Physical Interpretation

For an MSc level defense, understanding exactly what is plotted and *why* it behaves that way is the most critical part of the project. Below is a detailed breakdown of the resultant physics for every simulated plot.

### Phase 1: Baseline Sublattice Imbalance (`01_reproduce_fig3.py`)
- **The Initial State:** The system starts in a highly excited N\'eel state (particles strictly on alternating sites: `|1010...>`).
- **The Plot:** **X-axis:** Time $Jt$ (Logarithmic Scale). **Y-axis:** Sublattice Imbalance $\mathcal{I}(t)$ (measuring how closely the lattice resembles the initial striped pattern).
- **Physical Interpretation:** 
   - *Ergodic Curve ($w=1$):* The imbalance rapidly crashes to exactly $0.0$. The particles have chaotically hopped and scattered enough to smear out uniformly across the entire lattice. The system has "thermalized" and completely forgotten its initial state.
   - *MBL Curve ($w=10$):* The imbalance dips slightly but then rigidly plateaus at a high positive value permanently out to infinite time. Despite interactions allowing particles to hop, the strong disordered landscape freezes the wavefunctions. The system permanently remembers that it started on the even sites: a complete breakdown of quantum statistical mechanics.

### Phase 2: Finite Temperature Effects (`02_temperature_effect.py`)
- **The Plot:** **X-axis:** Time $Jt$. **Y-axis:** Sublattice Imbalance. We compare the infinite-temperature response ($\beta = 0$) against low-temperature/high-$\beta$ mixed density states.
- **Physical Interpretation:** Canonical MBL is defined strictly at infinite temperature (highly excited mid-spectrum states). By cooling the system (increasing $\beta$), we restrict the available energy phase space. Cold fermions have less kinetic energy to overcome standard interaction barriers, locking them into lower energy bands which technically alters the fundamental localization transition boundary. 

### Phase 3: External Heat Reservoir / Open Quantum Systems (`03_heat_reservoir_effect.py`)
- **The Plot:** **X-axis:** Time $Jt$. **Y-axis:** Sublattice Imbalance. We plot 4 different coupling strengths ($\gamma$) to an external environment.
- **Physical Interpretation:** MBL requires absolute, pristine isolation. The $\gamma$ parameter in the Lindblad Master Equation acts as the rate at which the external environment "measures" and dephases the fermions. 
   - At $\gamma = 0$ (Isolated), MBL survives perfectly. 
   - At $\gamma = 0.5$ (Strong Coupling), the environment's observation constantly shreds the fragile quantum phase coherence (superposition) of the particles. Without quantum interference, localization physically cannot exist, forcing the system to exponentially collapse into classical thermal equilibrium ($\mathcal{I} \to 0$).

### Phase 4: Interaction Tuning & The Mott Limit (`04_strong_interaction_effect.py`)
- **The Plot:** **X-axis:** Time $Jt$. **Y-axis:** Sublattice Imbalance. Sweeping interaction $U/J$ at a fixed $W=4.0$. 
- **Physical Interpretation:** This plot is non-monotonic! 
   - *Weak interactions ($U=2, 5$)* initially destroy localization by acting as a thermal bath that scrambles the particles, lowering the imbalance.
   - *Extreme interactions ($U=10$ and beyond)* force the imbalance to artificially shoot back up. This is **not** Many-Body Localization. When the repulsive penalty $U$ is vastly larger than the hopping energy $J$, particles cannot hop onto occupied neighboring sites at all. The system freezes simply because it is physically blocking itself (The **Mott Insulator** transition).

### Phase 5: Logarithmic Entanglement Entropy (`05_entanglement_entropy.py`)
- **The Plot:** `MBL_entanglement_growth_U_sweep.png`. **X-axis:** Time $Jt$ (Logarithmic Scale). **Y-axis:** Half-chain Bipartite von Neumann Entanglement Entropy $S(t)$. 
- **Setup:** We fix the system in a deeply localized regime ($W=10.0$) and sweep the interaction strength $U/J \in [0.0, 1.0, 5.0, 10.0]$.
- **Physical Interpretation of Results:** 
   - **$U=0.0$ (The Anderson Baseline):** The entropy rises briefly due to local particle movement but then hits a sharp ceiling and saturates almost instantly. Because there are no interactions, quantum information cannot leak between localized orbitals. This horizontal flatline is the definitive signature of an **Anderson Insulator**.
   - **$U > 0.0$ (Logarithmic MBL Growth):** As soon as interactions are added, the physics changes fundamentally. Even though the particles themselves are still spatially frozen (as seen in the Imbalance plot), their quantum phases begin to interact. This causes a slow, "dephasing" spread of entanglement.
   - **The Signatures:** On the semi-log plot, the curves for $U=1, 5, 10$ appear as straight diagonal lines sloping upwards. This proves that **$S(t) \sim \ln(t)$**, showing that information spreads exponentially slowly. This specifically replicates the hallmark discovery by Bardarson et al. (2012) and is the most rigorous proof that the system is in an interacting MBL phase rather than a simple non-interacting localized phase.
   - **Trend with $U$:** Increasing $U$ increases the slope (speed) of the logarithmic growth, demonstrating that stronger interactions facilitate faster (though still only logarithmic) entanglement spread.

### Phase 6: Wigner-Dyson Level Spacing Statistics (`06_level_spacing.py`)
- **The Plot:** **X-axis:** Disorder strength $W/J$. **Y-axis:** Average adjacent energy level spacing ratio $\langle r \rangle$.
- **Physical Interpretation:** The phase transition cleanly presents itself structurally in the actual energy eigenvalues of the system.
   - *Thermal Regime ($W < 5$):* High probability of states mixing. Energy levels "feel" each other and repel in the spectrum. This Quantum Chaos matches Random Matrix Theory (Wigner-Dyson Statistics) where $\langle r \rangle \approx 0.53$.
   - *MBL Regime ($W > 9$):* States are completely spatially isolated. Because they do not overlap physically, their energy levels don't repel each other. They mathematically cross each other harmlessly, collapsing perfectly to Poisson Statistics where $\langle r \rangle \approx 0.386$. The smooth curve downward cleanly maps the exact critical transition diagram of our Fermi-Hubbard lattice.

---

## Citations & Fundamental Literature
The architecture, mathematical methodology, and validation benchmarks of this MSc project were derived from the following canonical physics literature:

1. **QuSpin Engine & Fermi-Hubbard MBL Baseline Example:**
   *Weinberg, P., & Bukov, M. (2019).* "QuSpin: a Python package for dynamics and exact diagonalisation of quantum many body systems. Part II." **SciPost Physics**, 7(2), 020.
2. **Phase 5 Logarithmic Entanglement Growth (Figure 1):**
   *Bardarson, J. H., Pollmann, F., & Moore, J. E. (2012).* "Unbounded growth of entanglement in models of many-body localization." **Physical Review Letters**, 109(1), 017202.
3. **Phase 6 Level Spacing Ratio Phase Transition:**
   *Oganesyan, V., & Huse, D. A. (2007).* "Localization of interacting fermions at high temperature." **Physical Review B**, 75(15), 155111.
   *Pal, A., & Huse, D. A. (2010).* "Many-body localization phase transition." **Physical Review B**, 82(17), 174411.
