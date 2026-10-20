---
title: "Supporting Text S1 — derivations behind Section 3.3"
subtitle: "draft, 2026-10-20, untracked working file"
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
---

*Notes (delete before pasting): checked against Ribe (2001) on 2026-10-20; equation numbers cited for Ribe are his. S1.1 is the shear endmember, the derivation the submitted paper's Text S1 gave for the $H^1$ scaling, written here as an endmember rather than as the result. S1.2 is the moment-balance route that Section 3.3 now rests on. Coordinates and signs match the extraction and Tests A and C ($s$ down-dip, $n$ across the slab, $K = d\theta/ds$). One difference from the submitted Text S1 is deliberate: the old S1 wrote the shear stress as $\eta\, dv_n/ds$ and evaluated $dv_n/ds$ as the derivative of the component $\mathbf{v}\cdot\hat{\mathbf{n}}$, which equals $\omega + Kv_s$ (Equation S2 below); the $Kv_s$ part is the rotation of the coordinate basis, not a strain, so the shear-mode stress is $\eta\omega$. This changes the coefficients in Equation S6 relative to the old Equation 8 but not the dependence on $\eta$, $H$, $K$ and $v_s$.*

---

## Text S1. Derivation of the in-slab shear force gradient

### S1.0 Coordinates

We write the slab in coordinates aligned with its midsurface: $s$ is arc length along the midsurface, positive down-dip; $n$ is the distance from the midsurface; $\theta(s)$ is the dip; and $K = d\theta/ds$ is the curvature. With $x$ to the right and $y$ upward, the unit vectors along and across the slab are $\hat{\mathbf{t}} = (\cos\theta, -\sin\theta)$ and $\hat{\mathbf{n}} = (\sin\theta, \cos\theta)$, so that $v_s = v_x\cos\theta - v_y\sin\theta$ and $v_n = v_x\sin\theta + v_y\cos\theta$. Because $\theta$ varies along the slab, the unit vectors rotate with $s$:

$$
\frac{d\hat{\mathbf{t}}}{ds} = -K\hat{\mathbf{n}}, \qquad \frac{d\hat{\mathbf{n}}}{ds} = K\hat{\mathbf{t}}. \qquad (S1)
$$

### S1.1 The through-thickness shear endmember ($dQ_s/ds \propto H$)

The rotation rate of the slab midsurface is the rate at which its tangent turns, $\omega = \hat{\mathbf{n}}\cdot\partial\mathbf{v}/\partial s$. Written in terms of the velocity components,

$$
\omega = \frac{\partial(\mathbf{v}\cdot\hat{\mathbf{n}})}{\partial s} - \mathbf{v}\cdot\frac{d\hat{\mathbf{n}}}{ds} = \frac{\partial v_n}{\partial s} - Kv_s , \qquad (S2)
$$

which is the quantity that appears in Equation 6 of the main text. The second term is there because the component $v_n$ is measured against a direction that itself rotates along the slab; it is not a deformation. In the same way, the rotation rate of a material line drawn across the slab is $-\partial v_s/\partial n$. The shear strain rate at the midsurface is the difference between the rotation of a cross-section and the rotation of the midsurface it was normal to,

$$
2\dot\varepsilon_{sn} = \frac{\partial v_s}{\partial n} + \omega , \qquad (S3)
$$

and the shear stress is $\tau_{sn} = 2\eta\,\dot\varepsilon_{sn}$. Equation S3, with $\tau_{sn} = 2\eta\,\dot\varepsilon_{sn}$, is the exact midsurface shear stress of a curved viscous sheet (Ribe, 2001, the expression for $\sigma_{zs}$ below his equation 2.1, with his $\omega = w_0' + Ku_0$; his $K$ is positive with the centre of curvature on the $+z$ side, ours on the $-n$ side, which is the only difference). A rigid rotation of a curved slab, $v_s = \Omega(R + n)$ with $R = 1/K$, gives $\partial v_s/\partial n = \Omega$ and $\omega = -\Omega$, so Equation S3 vanishes as it must.

The shear endmember is the assumption that cross-sections translate without rotating, $\partial v_s/\partial n = 0$. All of the midsurface rotation must then be accommodated by shearing the section:

$$
\tau_{sn} = \eta\omega, \qquad Q_s = \int \tau_{sn}\,dn = \eta H \omega . \qquad (S4)
$$

For a slab that moves down-dip through a geometry that is steady in the slab frame, the velocity is along the slab, $\mathbf{v} = v_s\hat{\mathbf{t}}$, and Equation S2 with Equation S1 gives

$$
\omega = v_s\,\hat{\mathbf{n}}\cdot\frac{d\hat{\mathbf{t}}}{ds} = -Kv_s . \qquad (S5)
$$

The shear force is then $Q_s = -\eta H K v_s$, and its gradient along the slab is

$$
\frac{dQ_s}{ds} = -\eta H\left(K\frac{dv_s}{ds} + v_s\frac{dK}{ds}\right). \qquad (S6)
$$

If $K$ and $v_s$ vary over a common length $L$ along the slab, $|dQ_s/ds| \sim 2\eta H K v_s/L$, i.e.

$$
\frac{dQ_s}{ds} \propto \eta H K v_c \qquad (S7)
$$

with $v_s \approx v_c$ and a prefactor absorbing $H/L$. This is the scaling of the submitted manuscript. It is linear in plate thickness because the shear stress is assumed uniform across the slab, so that doubling $H$ doubles the shear force at the same stress.

The models do not deform this way, and thin-sheet theory says they should not: for a sheet much thinner than the length over which it deforms, the along-slab velocity varies across the sheet as $-\omega n$ at leading order (Ribe, 2001, equation 4.6a), which is the Kirchhoff behaviour. The measured $\partial v_s/\partial n$ equals $-\omega$ to within a few per cent across the whole suite (Figure S-bend b), so the two terms of Equation S3 cancel, the constitutive shear stress is several times smaller than $\eta|\omega|$, and the shear force gradient scales as $H^3$ rather than $H$ across the $H = 70$, 100 and 130 km models (Figure S-varH). Section 3.3 therefore uses the bending endmember, derived next.

### S1.2 The bending endmember: moment balance and the thin-sheet moment

When cross-sections rotate with the midsurface ($\partial v_s/\partial n = -\omega$), Equation S3 gives $\dot\varepsilon_{sn} \approx 0$ at leading order, and the shear stress can no longer be obtained from the constitutive law. It is fixed instead by equilibrium.

**Moment balance.** The along-slab component of the momentum equation is

$$
\frac{\partial\sigma_{ss}}{\partial s} + \frac{\partial\sigma_{sn}}{\partial n} + \rho g_s = 0 \qquad (S8)
$$

up to terms in $K\sigma_{sn}$ that come from the curvature of the coordinates (Ribe, 2001, equation 2.1b); $\sigma_{ss}$ is the total along-slab normal stress (deviatoric stress minus pressure), $\sigma_{sn} = \tau_{sn}$, and $g_s$ is the down-dip component of gravity. Multiplying by $n$ and integrating across the slab from $n = -H/2$ to $+H/2$ (the curvature terms cancel exactly in this weighted integral):

$$
\int n\,\frac{\partial\sigma_{ss}}{\partial s}\,dn = \frac{dM}{ds}, \qquad M = \int n\,\sigma_{ss}\,dn ,
$$

the bending moment;

$$
\int n\,\frac{\partial\sigma_{sn}}{\partial n}\,dn = \Big[n\,\tau_{sn}\Big]_{-H/2}^{+H/2} - \int \tau_{sn}\,dn = C - Q_s , \qquad C = \frac{H}{2}\left[\tau_{sn}\!\left(+\tfrac{H}{2}\right) + \tau_{sn}\!\left(-\tfrac{H}{2}\right)\right],
$$

the couple exerted by the mantle shear tractions on the two faces; and

$$
\int n\,\rho g_s\,dn = 0
$$

when the density is uniform across the slab (the lithostatic part of the pressure balances the reference density, and the excess density is symmetric about the midsurface). Hence

$$
Q_s = \frac{dM}{ds} + C . \qquad (S9)
$$

Equation S9 is the exact torque balance of thin viscous sheet theory (Ribe, 2001, equation 2.4; Buffett, 2006, equation 5, where the face tractions are neglected). Retaining the curvature of the slab in the volume element, $(1 + Kn)\,dn$, adds a couple from the slab's own weight, $\tfrac{1}{12}KH^3\Delta\rho\,g_s$ (Ribe, 2001, equation 2.4; it appears with $H$ in Buffett's equation 5, a misprint, since his equation 7 and the dimensions require $H^3$); in the models it is 2 % of $Q_s$ at the median and 13 % at the 90th percentile, and it is neglected here. The measured partition of $Q_s$ in the models is $dM/ds \approx 0.86\,Q_s$ and $C \approx -0.12\,Q_s$ at the median (Figure S-bend c), so $Q_s$ is carried almost entirely by the moment gradient.

**Constitutive relation for the moment.** With cross-sections rotating with the midsurface, the along-slab velocity varies linearly across the slab, $v_s(s, n) = v_s(s, 0) - n\,\omega(s)$, and so does the along-slab strain rate:

$$
\dot\varepsilon_{ss} = \Delta + n\dot\kappa, \qquad \dot\kappa = -\frac{\partial\omega}{\partial s} = \frac{\partial^2 v_s}{\partial s\,\partial n}, \qquad (S10)
$$

where $\Delta$ is the stretching rate of the midsurface and $\dot\kappa$ is the bending rate, $\Omega = -\omega'$ in Ribe (2001). The rate of change of curvature of the material midsurface is $\Omega + K\Delta$ (Ribe, 2001, after equation 3.1); the $K\Delta$ part, curvature times stretching, is small for slabs that do not stretch, and the moment it induces is dropped here, as in Ribe's own thin-sheet equations (his 4.10b versus 4.11). The slab is incompressible and in plane strain, so $\dot\varepsilon_{nn} = -\dot\varepsilon_{ss}$, and across a thin sheet the normal stress $\sigma_{nn} = -p + 2\eta\,\dot\varepsilon_{nn}$ is small compared with the bending stresses, which fixes the pressure at $p = -2\eta\,\dot\varepsilon_{ss}$ (the linear-in-$n$ pressure $p = -2\eta\Omega n$ of Ribe, 2001, equation 4.6b). The along-slab normal stress is therefore

$$
\sigma_{ss} = -p + 2\eta\,\dot\varepsilon_{ss} = 4\eta\,\dot\varepsilon_{ss} \qquad (S11)
$$

(Buffett, 2006, equation 12). Substituting Equations S10 and S11 into the definition of $M$, the $\Delta$ term integrates to zero and

$$
M = 4\eta\dot\kappa \int n^2\,dn = \frac{\eta H^3\dot\kappa}{3}, \qquad (S12)
$$

which is equation 4.8b of Ribe (2001), equation 16 of Ribe (2010) and equation 13 of Buffett (2006), the latter with $\dot\kappa = u_0\,dK/ds$ for a steady shape translating at speed $u_0$. This is Equation 10 of the main text. In the models the measured $M$ is within 10 % of Equation S12 at the median (Figure S-bend d).

**The shear force gradient.** Differentiating Equation S9 and neglecting $C$ and along-slab variations in $\eta$ and $H$,

$$
\frac{dQ_s}{ds} \approx \frac{d^2M}{ds^2} = \frac{\eta H^3}{3}\,\frac{d^2\dot\kappa}{ds^2}. \qquad (S13)
$$

For a slab moving through a curved geometry, $\omega \approx -Kv_s$ (Equation S5; in Ribe, 2001, equation 6.7, $\omega$ is the rate of change of the dip following the material, which for a slab sliding along a fixed shape is $v_s\,d\theta/ds$), so $\dot\kappa \approx \partial(Kv_s)/\partial s$ and Equation S13 contains the third along-slab derivative of $Kv_s$. Representing the along-slab variation of $Kv_s$ by a single length $L_K$, and taking $v_s \approx v_c$,

$$
\frac{dQ_s}{ds} \approx \frac{\eta H^3 K v_c}{3L_K^3}, \qquad (S14)
$$

which is Equation 12 of the main text. Dimensionally, $\eta$ [Pa s] $\times$ $K$ [m$^{-1}$] $\times$ $v_c$ [m s$^{-1}$] $\times$ $H^3/L_K^3$ is a stress, as $dQ_s/ds$ must be. This is the bending-dominated stress balance of Ribe (2001, Section 3), in which an applied normal or tangential stress is carried by $\eta H^3\omega/L^3$. $L_K$ is fitted to the models (Figure S-LK; 200 km); it is about twice the local curvature length $K/|dK/ds|$ of the modelled slabs (60–140 km), and much shorter than the length over which $v_s$ varies ($\sim 4500$ km), so $\dot\kappa \approx v_s\,dK/ds$. In the models $\omega$ tracks $-Kv_s$ with a slope of 0.75 overall and 0.86–0.89 for the 100 and 130 km slabs, the difference being the slabs' normal motion. The factor 3 is inherited from Equation S12 and is kept so that $L_K$ retains its meaning as the length over which $Kv_s$ varies; it is not an independent fitted constant. Fitting the exponent of $K$ across the suite gives 0.65 (model bootstrap 0.4–0.7), so a single $L_K$ independent of curvature is the right description; the alternative $L_K \propto 1/K$, which would make the dependence $K^4$, is excluded.

The ratio of the two endmembers, $(\eta H^3 K v_c/3L_K^3)/(\eta H K v_c/L) = H^2 L/(3L_K^3)$, shows why they cannot be told apart at a single plate thickness: at $H = 100$ km, $L_K = 200$ km and $L$ of order a few thousand kilometres, the ratio is of order one, and either form can be fitted to the $H = 100$ km models with its own length. Only the thickness dependence separates them (Figure S-varH).
