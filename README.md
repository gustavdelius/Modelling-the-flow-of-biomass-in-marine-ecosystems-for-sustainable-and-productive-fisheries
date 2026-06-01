# Understanding of the Project

## Overview

The central question is whether **fishing can improve, rather than harm, the productivity of marine ecosystems**. Conventionally, fishing is viewed as purely detrimental. This project challenges that assumption by modelling fish population dynamics as a transport problem — an analogy that reveals a mechanism by which targeted fishing could actively relieve ecological "traffic jams" and boost overall ecosystem productivity.

---

## The Traffic Analogy

### Traffic flow as a PDE

To build intuition, consider cars on a motorway. Let $n(x, t)$ be the density of cars at position $x$ and time $t$. By conservation of cars, the rate of change of density equals the net flux into a given road segment:

$$\frac{\partial n}{\partial t} = -\frac{\partial J}{\partial x}$$

The **flux** $J(x) = V(x) \cdot n(x)$ is the product of the local car velocity $V(x)$ and the car density $n(x)$. The minus sign makes physical sense: if the flux decreases along $x$ (more cars entering a segment than leaving it), the density in that segment must increase.

### The jamming instability

The crucial non-linearity is that velocity $V$ depends on density $n$. When traffic becomes dense, drivers slow down. This creates a positive feedback loop:

> High density → lower velocity → cars arriving faster than they leave → even higher density

This is the **phantom traffic jam**: a smooth flow can spontaneously break down into a propagating jam from a single small perturbation, but only when density is above a critical threshold. Motorway agencies exploit this insight by imposing variable speed limits — counterintuitively *reducing* speed to prevent jams from forming in the first place.

---

## The Fish Population Equation

### From cars to fish

Now replace spatial position $x$ with body mass $w$, and cars with fish. Let $n(w, t)$ be the **number density of fish at body mass $w$**. Fish begin as tiny eggs and grow continuously to large adults, so biomass is transported from small $w$ to large $w$ — exactly like cars moving along a road.

The governing equation (known as the McKendrick–von Foerster equation) is:

$$\frac{\partial n}{\partial t} = -\frac{\partial}{\partial w}\big[g(w)\, n(w)\big] - \mu(w)\, n(w)$$

where:

| Symbol | Meaning |
|--------|---------|
| $g(w)$ | Growth rate of an individual fish of mass $w$ (the analogue of car velocity) |
| $\mu(w)$ | Size-dependent mortality rate (natural predation **plus** fishing mortality) |

The mortality term $\mu(w)\,n(w)$ means this is not a pure conservation equation — fish can leave the system by dying at any size, not only by reaching the largest size class.

### The non-linearity: growth depends on density

Fish grow exclusively by eating. The prey that determine $g(w)$ are smaller fish (at masses far to the *left* on the $w$-axis), so the dependence of $G$ on $n$ is **non-local** — unlike traffic, where a car's speed depends only on the density immediately nearby.

This produces the same feedback mechanism as traffic jams:

> High fish abundance at mass $w$ → prey depleted → $g(w)$ falls → fish spend longer at mass $w$ → local density increases further → prey depleted further

This self-reinforcing cycle is the **jamming instability in size space**. A small random perturbation in abundance at one size class can amplify into a full pile-up, provided the system is in the right parameter regime.

### The intervention: fishing as a speed limit

On a motorway, variable speed limits artificially reduce velocity to prevent pile-ups. We cannot impose speed limits on fish, but we can do something motorway agencies cannot: **selectively remove individuals**. Because fishing mortality is a component of $\mu(w)$, fishing is our primary control lever. The central hypothesis of this project is that targeted fishing at the right size class could relieve a jam, improving overall ecosystem productivity rather than diminishing it.

---

## The MIZER Model

Simulations will be run using [MIZER](https://sizespectrum.org/mizer/), an R package for size-spectrum modelling. MIZER computes $g(w)$ via a **predation kernel** that specifies the preferred prey sizes for a predator of mass $w$, then integrates that preference over the actual size spectrum to obtain the growth rate. This naturally captures the non-local, density-dependent structure of $g(w)$: increased predator abundance directly suppresses $g$ through prey depletion.

### Plankton resource dynamics

Fish also consume plankton. MIZER models the plankton resource $R(w)$ using **semi-chemostat dynamics**:

$$\frac{dR}{dt} = r(w)\,\big[c(w) - R(w)\big] - \mu_R(w)\, R(w)\, n(w)$$

where:

| Symbol | Meaning |
|--------|---------|
| $r(w)$ | Intrinsic replenishment rate of the resource |
| $c(w)$ | Carrying capacity (the abundance the resource reaches without predation) |
| $\mu_R(w)\,R(w)\,n(w)$ | Losses due to predation by fish |

Without predation, the resource grows toward $c(w)$. With predation, it settles at a lower equilibrium determined by the balance between replenishment and consumption.

### Why a low replenishment rate triggers jams

Setting $dR/dt = 0$ and solving for the equilibrium resource abundance gives approximately:

$$R^* \approx \frac{r \cdot c}{r + \mu_R \cdot n}$$

Two contrasting regimes emerge:

- **Large $r$:** $R^* \approx c$ regardless of how many fish are present. The resource replenishes quickly, so $g(w)$ is nearly independent of fish density — no feedback, no jam.
- **Small $r$:** $R^*$ is highly sensitive to $n$. A small increase in fish abundance sharply depletes the resource, strongly reducing $g(w)$, causing fish to accumulate at that size — triggering the jam.

A low replenishment rate therefore maximises the coupling between fish density and growth rate, placing the system in the unstable regime where jams can form spontaneously from small random fluctuations in abundance.

---


