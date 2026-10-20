---
title: "Section 3.3 — redraft v3"
subtitle: "2026-10-20, untracked working file"
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{amsmath}
---

*Notes (delete before pasting): checked against Ribe (2001) on 2026-10-20; equation numbers cited for Ribe are his. "Figure S-bend" = bending_supp.z300.0; "Figure S-varH" = H-scaling_supp; "Figure S-LK" = LK_fit; "Text S1" = the supplement derivation file (shear endmember in S1.1, moment-balance details in S1.2). Relative to v2: the kinematics are set out explicitly, $M$ and the moment balance are explained and their origin stated, the constitutive step for $M$ is given, the neglected terms are named with their measured size, and the placeholders are filled. The claim "an order of magnitude below $\eta\omega$" has been corrected to what the whole population shows. The steps $\omega \approx -Kv_s$ and $L_K$ versus the curvature length are now measured (2026-10-20).*

---

## 3.3 Scaling of the in-slab shear force gradient

Besides the mantle pressure difference, the main term balancing slab buoyancy is the gradient of the in-slab shear force, $dQ_s/ds$ (Figures 3–4). To know when slab dip can be read as a pressure constraint, we need a scaling for this term in quantities observable on Earth: plate thickness, slab viscosity, curvature and convergence rate.

**Kinematics.** We use coordinates aligned with the slab: $s$ is arc length along the slab midsurface, positive down-dip; $n$ is the distance from the midsurface; and $K = d\theta/ds$ is the curvature, with $\theta$ the dip. The shear force is the integral of the shear stress across the slab, $Q_s = \int \tau_{sn}\,dn$, and for a viscous slab $\tau_{sn} = 2\eta\,\dot\varepsilon_{sn}$. At the midsurface the shear strain rate is

$$
2\dot\varepsilon_{sn} = \frac{\partial v_s}{\partial n} + \omega, \qquad \omega = \frac{\partial v_n}{\partial s} - K v_s , \qquad (6)
$$

where the term $-Kv_s$ arises because the $s$ and $n$ directions rotate along the slab; Equation 6 is the exact expression for the shear strain rate at the midsurface of a curved sheet (Ribe, 2001). It has a simple reading: $\omega$ is the rotation rate of the slab midsurface, and $-\partial v_s/\partial n$ is the rotation rate of a material line drawn across the slab. The shear strain rate is the difference between the two, i.e., the rate at which cross-sections tilt relative to the midsurface they were normal to.

**Two endmembers.** If cross-sections translate without rotating ($\partial v_s/\partial n = 0$), all of the midsurface rotation must be taken up by shearing the section, and $\tau_{sn} = \eta\omega$. For a slab moving down-dip through a fixed curved geometry, $\omega \approx -Kv_s$, so $\tau_{sn} \approx -\eta K v_s$ is uniform across the slab, $Q_s \approx -\eta H K v_s$, and the shear force gradient is linear in plate thickness, $dQ_s/ds \sim \eta H K v_c/L$, with $L$ the length over which $Kv_s$ varies (Supporting Text S1.1). This is the through-thickness shear endmember. If instead cross-sections rotate with the midsurface ($\partial v_s/\partial n = -\omega$), the two terms in Equation 6 cancel at leading order. This is not an assumption of thin viscous sheet theory but a result of it: for a sheet much thinner than the length over which it deforms, the along-slab velocity varies across the sheet as $-\omega n$ at leading order (Ribe, 2001). The shear stress is then not set by the constitutive law at all: it is whatever equilibrium requires, and it must be found from the balance of moments across the slab. This is the bending endmember.

Our slabs follow the second behaviour. Across every model, the measured $\partial v_s/\partial n$ tracks $-\omega$ with a slope of 1.0 (Figure S-bend b), the shear stress inside the slab is several times smaller than the shear-endmember estimate $\eta|\omega|$ (typically a fifth of it), and the along-slab normal stress varies linearly across the slab, compressional on one face and tensional on the other (Figure S-bend a). That stress profile is what a bending slab carries.

**Moment balance.** The torque per unit trench length exerted by the along-slab normal stress about the midsurface is the bending moment,

$$
M = \int n\,\sigma_{ss}\,dn , \qquad (7)
$$

where $\sigma_{ss}$ is the total (pressure included) along-slab normal stress. Integrating the along-slab component of the momentum equation across the slab, weighted by $n$, gives (Supporting Text S1.2)

$$
Q_s = \frac{dM}{ds} + C, \qquad C = \frac{H}{2}\left[\tau_{sn}\!\left(+\tfrac{H}{2}\right) + \tau_{sn}\!\left(-\tfrac{H}{2}\right)\right], \qquad (8)
$$

where $C$ is the couple exerted by the mantle shear tractions on the two slab faces. Equation 8 is the torque balance of thin-sheet theory (Ribe, 2001, his equation 2.4; Buffett, 2006, who neglects the face tractions), and it is exact: a gradient in bending moment along the slab must be balanced by a shear force across it, just as in a beam. Two terms are dropped. The couple from the slab's own weight, which is non-zero only because the slab is curved, is $\tfrac{1}{12} K H^3 \Delta\rho\, g \sin\theta$ (Ribe, 2001) and amounts to 2 % of $Q_s$ in the models (median; 13 % at the 90th percentile). The face couple $C$ is measured directly and is small, $-0.12\,Q_s$ at the median, so $Q_s$ is carried almost entirely by $dM/ds$ (Figure S-bend c).

**Constitutive relation for $M$.** When cross-sections rotate with the midsurface, the along-slab strain rate varies linearly across the slab, $\dot\varepsilon_{ss} = n\dot\kappa$, where

$$
\dot\kappa = \frac{\partial^2 v_s}{\partial s\,\partial n} = -\frac{\partial\omega}{\partial s} \qquad (9)
$$

is the bending rate of the midsurface ($\Omega$ in Ribe, 2001), equal to the rate of change of curvature of the material midsurface when stretching is negligible. For an incompressible viscous slab in plane strain, the pressure adjusts so that the normal stress across the thin slab is small, and $\sigma_{ss} = 4\eta\,\dot\varepsilon_{ss}$ (Ribe, 2001; Buffett, 2006). Equation 7 then gives

$$
M = \frac{\eta H^3 \dot\kappa}{3} \qquad (10)
$$

(Ribe, 2001, his equation 4.8b; Buffett, 2006). This is the viscous analogue of the elastic flexural rigidity, and it is the only place plate thickness enters cubed. Equation 10 holds closely in the models, with the measured $M$ within 10 % of the thin-sheet value at the median (Figure S-bend d).

**The scaling.** Neglecting $C$ and along-slab variations in $\eta$ and $H$, Equations 8 and 10 give

$$
\frac{dQ_s}{ds} \approx \frac{d^2 M}{ds^2} = \frac{\eta H^3}{3}\,\frac{d^2\dot\kappa}{ds^2}. \qquad (11)
$$

For a slab moving through a curved geometry $\omega \approx -Kv_s$ (in the models $\omega$ tracks $-Kv_s$ with a slope of 0.75 overall and 0.86–0.89 for the 100 and 130 km slabs, the difference being the slabs' normal motion), so $\dot\kappa \approx \partial(Kv_s)/\partial s$ and Equation 11 involves the third along-slab derivative of $Kv_s$. That derivative cannot be measured reliably on Earth, and is dominated by noise even in the models. We therefore represent it with a single length scale $L_K$ over which $Kv_s$ varies along the slab, and take $v_s \approx v_c$:

$$
\frac{dQ_s}{ds} \approx \frac{\eta H^3 K v_c}{3 L_K^3}. \qquad (12)
$$

This is the bending resistance of thin-sheet theory, $\eta H^3\omega/L^3$ for a sheet rotating at rate $\omega$ over a length $L$ (Ribe, 2001; written as $\eta h^3 V/\ell_b^4$ with the bending length $\ell_b$ in Ribe, 2010), with the measured curvature standing in for one power of the length. $L_K$ is fitted to the models; we use $L_K = 200$ km (Supporting Figure S-LK), about twice the length over which curvature varies along the modelled slabs ($K/|dK/ds| \approx 60$–140 km at 300 km depth), while the along-slab velocity varies over thousands of kilometres, so $\dot\kappa$ is dominated by $v_s\,dK/ds$.

The $H^3$ dependence is what distinguishes bending from through-thickness shear, and it can only be tested by varying plate thickness: at a single thickness, a scaling linear in $H$ with a fitted length fits equally well. Across our $H = 70$, 100 and 130 km models, only the cubic form brings all thicknesses onto one trend (Supporting Figure S-varH). The thinnest slabs carry almost no shear force gradient, as Equation 12 predicts. The curvature dependence is also as written: fitting the exponent of $K$ across the suite gives 0.65 with a model-bootstrap range of 0.4–0.7, so a single length $L_K$, independent of curvature, is the right description.

Figure 5b compares the modelled $dQ_s/ds$ with Equation 12. The model means follow the 1:1 line closely ($r = 0.94$). Individual outputs scatter more, as expected when one length scale stands in for the along-slab structure of each slab: the scatter is the temporal variability of $dQ_s/ds$ at a fixed depth as bends migrate through it, not measurement noise. Overturned slabs follow the scaling to $\sim 15$ MPa and then saturate, which we attribute to the breakdown of the assumed geometry at very high curvature (Section 3.2).

Dividing Equation 12 by the slab buoyancy $B = \Delta\rho\, g H$ defines

$$
\Lambda = \frac{\eta H^3 K v_c}{3 L_K^3 B} = \frac{\eta H^2 K v_c}{3 L_K^3\,\Delta\rho\, g}, \qquad (13)
$$

the size of the in-slab stress term relative to buoyancy. Where $\Lambda < 0.1$, the modelled $\Delta P$ matches $B_{slab}$ to within 10 % of $B$ in 88 % of cases, and this agreement breaks down quickly above that value (Figure 5c). We therefore use $\Lambda < 0.1$ as the condition for reading slab dip as a measure of mantle pressure. Because $B$ itself scales with $H$, $\Lambda$ scales with $H^2$, so thin plates satisfy this condition more readily than thick ones (Section 4.2).

---

*Assumptions carried into Equation 12, with their measured sizes: (i) cross-sections rotate with the midsurface (slope of $-\partial v_s/\partial n$ on $\omega$: 1.04, $r = 0.99$); (ii) the face couple $C$ is neglected ($-0.12\,Q_s$ at the median); (iii) the body-force couple is neglected (2 % of $Q_s$); (iv) $\eta$ and $H$ are uniform along $s$; (v) the slab's normal velocity is small, so $\omega \approx -Kv_s$ (slope 0.75 overall, 0.86–0.89 for $H = 100$ and 130 km; $|v_n|/|v_s| \approx 0.3$ supplies the shortfall); (vi) stretching is small, so the stretching contribution to the moment (the $K\Delta$ term of Ribe, 2001, equation 4.10b, which he also drops in his thin-sheet equations) is neglected; (vii) $v_s \approx v_c$, and the $K\,dv_s/ds$ term of $\dot\kappa$ is dropped (the velocity length $|v_s/(dv_s/ds)| \approx 4500$ km is far longer than $L_K$); (viii) $L_K$ is a single constant: the fitted curvature exponent is 0.65 (0.4–0.7), not the 4 that $L_K \propto 1/K$ would give, and on Earth the headline fraction is 93–100 % for any $L_K$ between 150 and 300 km, with the transition viscosity moving by a factor of 1.3 across the bootstrap range 194–212 km.*
