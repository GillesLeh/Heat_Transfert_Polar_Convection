# Heat Transfert polar convection 

# Content

* [Objectives](#objectives)
* [About this repository](#about-this-repository)
* [Code description](#code-description)

# Objectives
Computer programmes are useful tools for understanding concepts that are not accessible to us. In this contribution I present the equivalence of a Matlab code in Python that allows to compute pressure, velocity, stress and strain fields in a two-dimensional system in cylindrical polar coordinates for a case of bottom-heated convection with an isothermal mantle. 

# About this repository
In this repository you can find:
- the main Python code [Essais_v9.py](Essais_v9.py)
- the Matlab code base [] ()
- the code that allows the differences between the two main codes to be perceived [comparateur3.py](comparateur3.py)

# Code description
* [Equation of state](#Equation-of-state)
* [Equation of Mass Conservation](#Equation-of-Mass-Conservation)
* [Constitutive Equations](#Constitutive-Equations)
* [Equation of Conservation of Linear Momentum](#Equation-of-Conservation-of-Linear-Momentum)
* [Equation of Energy Conservation](#Equation-of-Energy-Conservation)
* [Effective viscosities](#Effective-viscosities)

## Equation of state
Because of the mainly temperature and compositional differences the flow in the mantle is driven by local variation in density. The equation of state, in the form of a linearized equation (LEOS), allows calculation of the effective density times the gravitational acceleration of the system dependent on the temperature:

$\rho_{f}(T) = \rho_{0}(1-\alpha\Delta T + \beta\Delta P)$
```md
rhofg = rhog_v * (1 - alph * T + beta * P)
```
But we can ignore the pressure because our material is incompressible.

⤴️ [_back to content_](#content)

## Equation of Mass Conservation
We then use the continuity equation. In this equation mass is conserved and any change in density must be accompanied by an inflow or outflow of material, as the volume remains constant. 

The continuity equation:
$\frac{\delta\rho}{\delta t} + div(\rho\vec{v}) = 0,$

We can couple this continuity equation with the beta compressibility (1/bulk modulus):
$\beta = \frac{1}{\rho}(\frac{\delta\rho}{\delta p})$

by combining these two equations we obtain:
$0 = \frac{d P}{d t}-\frac{1}{\beta}\nabla V$

```md
dPdt = -1/bet * divV
```
∇V corresponds to the divergence of velocities calculated by the expression: 

$\nabla V = \xi_{\varphi\varphi} + \xi_{rr}$

in which $\xi_{rr}$ and $\xi_{\varphi\varphi}$ are the deviatoric strain rate in the radial and angular direction defined by:

$\xi_{rr} = \frac{\delta V_{r}}{\delta r}$

$\xi_{\varphi\varphi} = \frac{1}{r}\frac{\delta V_{\varphi}}{\delta \varphi}+\frac{V_{r}}{r}$

The $\xi_{\varphi r}$ represents the deviatoric tensor of viscous stresses:

$\xi_{\varphi r} = \xi_{r \varphi} = \frac{1}{2}(\frac{\delta V_{r}}{\delta \varphi} + \frac{\delta V_{\varphi}}{\delta r}-\frac{V_{\varphi}}{r})$

```md
Err = np.diff(Vr, axis=0)/np.diff(radr, axis=0)
Epp = (np.diff(Vp, axis=1)/np.diff(phip, axis=1) + avr(Vr)) / radn
Erp = ((np.diff(Vr[1:-1, :], axis=1)/np.diff(phir[1:-1, :], axis=1) - avr(Vp[:, 1:-1] )) / radc + np.diff(Vp[:, 1:-1], axis=0)/np.diff(radp[:, 1:-1], axis=0))/2
```

⤴️ [_back to code description_](#code-description)

## Constitutive Equations
the constitutive equations explain how a material deforms under certain pressures and stresses. The values $\tau_{rr}$ and $\tau_{\varphi\varphi}$ describe the normal stress in the radial and angular directions, while $\tau_{r \varphi}$ describes the shear stress. All three values are determined by these equations:

$\tau_{rr} = 2\eta(\dot{\xi}_{rr}-\frac{1}{3}\nabla V)$

$\tau_{\varphi\varphi} = 2\eta(\dot{\xi}_{\varphi\varphi}-\frac{1}{3}\nabla V)$

$\tau_{r \varphi} = 2\eta\dot{\xi} _{r \varphi}$

```md 
taurr = 2 * Eta * (Err - 1/3 * divV)
taupp = 2 * Eta * (Epp - 1/3 * divV)
taurp = 2 * Eta_rp * Erp
```

⤴️ [_back to code description_](#code-description)

## Equation of Conservation of Linear Momentum 
We use a force balance to understand the pressure, dynamic stresses, and gravity forces acting on the mantle. Since the mantle is highly viscous, we don't need to consider the effects of inertia. The Navier-Stokes equations are used to describe the behavior of the material, including its elasticity and viscosity. These equations help us understand the specific mechanical properties of the mantle's materials. By relating velocity and strain rate, we can understand how the mantle rock moves and deforms.

$\frac{dV_r}{dt}= \frac{1}{\rho}(\frac{\delta\sigma_{rr}}{\delta r} + \frac{1}{r}\frac{\delta\tau_{r\varphi}}{\delta\varphi} + \frac{\Delta\sigma}{r} -\rho fg)$

$\frac{dV_{\varphi}}{dt} = \frac{1}{\rho}(\frac{1}{r}\frac{\delta\sigma_{\varphi\varphi}}{\delta\varphi} + \frac{\delta\tau_{\varphi}}{\delta r} + 2\frac{\tau_{r\varphi}}{r})$

```md
dVrdt = 1/rho * (np.diff(Srr, axis=0)/np.diff(radn[:, 1:-1], axis=0) +
                (np.diff(taurp, axis=1)/np.diff(phic, axis=1) + deltaS) / radr[1:-1, 1:-1] -
                (rhofg[:-1, 1:-1] + rhofg[1:, 1:-1])/2)
dVpdt = 1/rho * (np.diff(taurp, axis=0)/np.diff(radc, axis=0) +
                (np.diff(Spp, axis=1)/np.diff(phin[1:-1, :], axis=1) + 2 * avr(taurp)) / radp[1:-1, 1:-1])
```
where $\rho$ corresponds to the mechanical density, $\rho f$ corresponds to the effective density calculated in the equation of state.

the variables $\delta\sigma_{rr}$ and $\delta\sigma_{\varphi\varphi}$ define the difference between cells of radial and angular total stress.

⤴️ [_back to code description_](#code-description)

## Equation of Energy Conservation
Using the energy equation, we can determine the temperature distribution on which density and rheology depend. The temperature equations allow calculation of temperature changes by advection and diffusion as well as the update of the model temperature.

$\frac{dT}{dt} = -V_{r}\frac{\delta T}{\delta r}-\frac{V_{\varphi}}{r}\frac{\delta T}{\delta\varphi}+\frac{1}{\rho C_{p}}(\frac{\delta}{\delta r}(\lambda(\frac{\delta T}{\delta r})+\frac{\lambda}{r}\frac{\delta T}{\delta r}+\frac{1}{r^2}\frac{\delta}{\delta\varphi}(\lambda(\frac{\delta T}{\delta\varphi}))$

```md
dTdt_1 = - np.maximum(0, Vr[1:-2, 1:-1]) * np.diff(T[:-1, 1:-1], axis=0) / drad 
dTdt_2 = - np.minimum(0, Vr[2:-1, 1:-1]) * np.diff(T[1:, 1:-1], axis=0) / drad 
dTdt_3 = - np.maximum(0, Vp[1:-1, 1:-2]) * np.diff(T[1:-1, :-1], axis=1) / dphi / radn[1:-1, 1:-1]
dTdt_4 = - np.minimum(0, Vp[1:-1, 2:-1]) * np.diff(T[1:-1, 1:], axis=1) / dphi / radn[1:-1, 1:-1]
dTdt_5 = (np.diff(lam * np.diff(T[:, 1:-1], axis=0) / drad, axis=0) / drad) / rhoCp
dTdt_6 = (lam * np.diff(avr(T[:, 1:-1]), axis=0) / drad / radn[1:-1, 1:-1]) / rhoCp
dTdt_7 = (np.diff(lam * np.diff(T[1:-1, :], axis=1) / dphi, axis=1) / dphi / radn[1:-1, 1:-1] ** 2) / rhoCp
dTdt = dTdt_1 + dTdt_2 + dTdt_3 + dTdt_4 + dTdt_5 + dTdt_6 + dTdt_7
```
⤴️ [_back to code description_](#code-description)

## Effective viscosities
The effective viscosity, $\eta$, which is used in the above equations can define several types of viscous flow. In the original Matlab code, $\eta L$ represented linear (Newtonian) viscous flow. This represents diffusion creep. In our case, I have added a non-linear viscous flow power law type with the effective viscosity depending on the strain rate. The effective viscosity for a power-law viscous fluid, termed here $\eta PL$, can be written as:

$\eta PL = \eta L (\frac{T_{II}}{T_{R}})\^{1-\eta}$

where $\eta L$ is the linear viscosity, $T_{R}$ is a constant reference stress, $\eta$ is the stress exponent, which for rocks is ≥1, and

$T_{II} = (0.5(T_{rr}^2+T{\varphi\varphi}^2)+(T_{r \varphi})^2)^{0.5}$

These two viscosities were averaged by a harmonic mean:

$\eta C = \frac{1}{\frac{1}{\eta L}+\frac{1}{\eta PL}}$

```md
tau2 = (0.5*(taurr**2+taupp**2)+(c2n(taurp)**2))**0.5

if n_exp > 1:
    Eta_pl_it = Eta_pl
    Eta_pl    = Eta*((tau2/s_ref)**(1-n_exp))
    Eta_pl    = np.exp(np.log(Eta_pl*rel+np.log(Eta_pl_it)*rel))
    Eta       = 1/(1/Eta_l + 1/Eta_pl)
```
⤴️ [_back to code description_](#code-description)
