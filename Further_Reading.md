## Further Reading & Resources

### Papers

1.  Nagel & Schreckenberg (1992) — A Cellular Automaton Model for Freeway Traffic Journal de Physique I, 2, 2221–2229 https://ui.adsabs.harvard.edu/abs/1992JPhy1...2.2221N The original discrete traffic model. Monte Carlo simulations show a sharp phase transition from free flow to stop-and-go jams as vehicle density crosses a critical threshold — the clearest demonstration that jamming is a collective, density-driven phenomenon, not caused by a bottleneck.

---
2. Bando, Hasebe, Nakayama, Shibata & Sugiyama (1995) — Dynamical Model of Traffic Congestion and Numerical Simulation
Physical Review E, 51, 1035–1042
https://doi.org/10.1103/PhysRevE.51.1035
Introduces the Optimal Velocity (OV) model, a microscopic car-following model. Stability analysis reveals a Hopf bifurcation at a critical density: below it, uniform flow is stable; above it, small perturbations grow into oscillatory stop-and-go waves. The bifurcation parameter maps cleanly onto the fish growth-rate sensitivity.
---

3.  Kerner & Konhäuser (1994) — Structure and Parameters of Clusters in Traffic Flow Physical Review E, 50(1), 54–83 https://doi.org/10.1103/PhysRevE.50.54 Studies the nonlinear structure of jams (clusters) once they form. Derives explicit formulas for the density, velocity, and propagation speed inside and outside a jam as functions of the background density — i.e. the constants that characterise whether a jam will form and what it will look like.

---
4. Richards (1956) — Shock Waves on the Highway
Operations Research, 4(1), 42–51
https://doi.org/10.1287/opre.4.1.42
The third paper in the LWR trio (with Lighthill & Whitham 1955). Shows that the conservation PDE develops shock waves — sharp discontinuities in density — precisely when the flux-density curve is concave. The condition for shock formation is the direct traffic analogue of the jamming instability criterion in the fish model.
---

5.  Flynn, Kasimov, Nave, Rosales & Seibold (2009) — Self-Sustained Nonlinear Waves in Traffic Flow Physical Review E, 79, 056113 https://doi.org/10.1103/PhysRevE.79.056113 Proves analytically that the Payne–Whitham second-order traffic model admits jamitons — self-sustaining nonlinear traveling waves analogous to detonation waves in combustion. Derives existence conditions and wave speed from first principles. The closest mathematical paper to the fish "pile-up" wave the project is trying to find.

---
6. Seibold, Flynn, Kasimov & Rosales (2013) — Constructing Set-Valued Fundamental Diagrams from Jamiton Solutions
Networks and Heterogeneous Media, 8(3), 745–772
https://arxiv.org/abs/0809.2828
Shows how the multi-valued (unstable) region of the fundamental diagram — the flow-density curve — arises directly from jamiton solutions. The set-valued region is exactly the regime where jams spontaneously form, giving a precise mathematical characterisation of the unstable parameter space.
---

7.  Sugiyama et al. (2008) — Traffic Jams Without Bottlenecks: Experimental Evidence New Journal of Physics, 10, 033001 https://doi.org/10.1088/1367-2630/10/3/033001 The landmark ring-road experiment: 22 cars drive in a circle and, without any obstruction, a jam emerges spontaneously above a critical density of 22 vehicles. Provides empirical confirmation that the Hopf bifurcation predicted by OV models is real. This is the traffic equivalent of what the MIZER simulations are trying to demonstrate for fish.

---
8. Colombo (2002) — Hyperbolic Phase Transitions in Traffic Flow
SIAM Journal on Applied Mathematics, 63, 708–721
https://doi.org/10.1137/S0036139901393184
Formulates traffic as a system with two distinct phases — free flow (scalar conservation law) and congested flow (2×2 system) — coupled at a free boundary. The phase transition is mathematically sharp and rigorously solved. Directly relevant to understanding the fish model's transition between "smooth spectrum" and "jammed" states.
---

9.  Bifurcation Analysis of Macroscopic Traffic Flow Models Based on Driver Behavior (2023)\* ResearchGate preprint / journal https://www.researchgate.net/publication/376606973_Bifurcation_analysis_of_macroscopic_traffic_flow_models_based_on_driver_behavior A more recent paper that explicitly proves the existence of Hopf bifurcations and saddle-node bifurcations in macroscopic traffic flow models as sensitivity parameters vary. Directly bridges the traffic and dynamical-systems language — useful for understanding what type of instability the fish project is dealing with.

---
10. Datta, Delius & Law (2011) — A Stability Analysis of the Power-Law Steady State of Marine Size Spectra
Journal of Mathematical Biology, 63(4), 779–799traffic flow models as sensitivity parameters vary. Directly bridges the traffic and dynamical-systems language — useful for understanding what type of instability the fish project is dealing with.




11. Sheldon, Prakash & Sutcliffe (1972) — The Size Distribution of Particles in the Ocean
Limnology and Oceanography, 17(3), 327–340
https://doi.org/10.4319/lo.1972.17.3.0327
The empirical starting point for the whole field. Shows that roughly equal biomass exists at every logarithmic size class from bacteria to whales — the observation that motivated size-spectrum modelling.
---

12. Andersen & Beyer (2006) — Asymptotic Size Determines Species Abundance in the Marine Size Spectrum The American Naturalist, 168(1) https://doi.org/10.1086/504849 Derives, from first principles, how fish abundance scales with body mass using the predation kernel and physiology. The theoretical backbone of MIZER.

---
13. Hartvig, Andersen & Beyer (2011) — Food Webd Populations
Journal of Theoretical Biology, 272(1), 113–122
https://doi.org/10.1016/j.jtbi.2010.12.006
The paper that defines the MIZER framework in he semi-chemostat plankton dynamics, and theMcKendrick–von Foerster PDE as used in this project. Essential reading.
---

14. Benoît & Rochet (2004) — A Continuous Model of Biomass Size Spectra Governed by Predation and the Effects of Fishing Journal of Theoretical Biology, 226(1) https://doi.org/10.1016/j.jtbi.2003.08.007 Introduces the continuous PDE model for marinerecedes MIZER, and already examines how fishing deforms the spectrum. Good for understanding where the equations come from.

---
15. Scott, Blanchard & Andersen (2014) — mizer: An R Package for Multispecies, Trait-Based and Community Size Spectrum Ecological Modelling
Methods in Ecology and Evolution, 5(10)
https://doi.org/10.1111/2041-210X.12256
The primary software paper for MIZER. Read this before running any simulations — it explains how the equations are
discretised and which parameters control what.
---

16. Andersen & Pedersen (2010) — Fishing Destabilises the Biomass Flow in Marine Size Spectra Proceedings of the Royal Society B, 279(1727), 284–292 https://pmc.ncbi.nlm.nih.gov/articles/PMC3223683/ Shows that fishing amplifies oscillations in the size spectrum, with selective fishing (targeting a narrow size range) being most destabilising. Directly relevant to the question of whether fishing can relieve jams.

---
17. Brechner et al. (2021) — Spatial Drivers of Instability in Marine Size-Spectrum Ecosystems
Journal of Theoretical Biology, 516, 110631
https://doi.org/10.1016/j.jtbi.2021.110631
Extends size-spectrum instability into space, showing that spatial movement and fitness-taxis can produce Turing-like
patterns of abundance. A natural next step oncial jamming instability.
---

18. Zhang & Chen (2017) — Evaluating Fishing Effects on the Stability of Fish Communities Using a Size-Spectrum Model Fisheries Research https://www.sciencedirect.com/science/article/ Applied study using a size-spectrum model to test how varying fishing pressure affects community stability indicators. Good for seeing how the theoretical instability plays out in a practical management context.

---
19. Lighthill & Whitham (1955) — On Kinematic ong Rivers
Proceedings of the Royal Society A, 229(1178)
https://doi.org/10.1098/rspa.1955.0112
The original paper deriving the LWR conservation PDE — the exact traffic-flow equation your supervisor used as the analogy for the fish model. Reading the source makes the mathematical parallel to the McKendrick–von Foerster equation very clear.
---

20.https://rpubs.com/gustav/plankton-anchovy

### Online Tools & Simulations

-   **MIZER official documentation & tutorials** — <https://sizespectrum.org/mizer/>
-   **MIZER interactive course (fishing scenarios)** — <https://course.mizer.sizespectrum.org/use/fishing-scenarios.html>
-   **Community Size Spectrum Simulator (DTU Aqua)** — <http://oceanlife.dtuaqua.dk/cspectrum/>
-   **Phantom traffic jam simulator (Martin Treiber)** — <https://mtreiber.de/MicroApplet/RingRoad.html>
-   **MIT traffic & jamitons interactive** — <https://math.mit.edu/traffic/>
