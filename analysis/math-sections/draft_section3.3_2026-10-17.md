# Section 3.3 — redraft v2 (2026-10-17, untracked working file)

*Notes (delete before pasting): [brackets] = check or fill. "Figure S-bend" = bending_supp.z300.0; "Figure S-varH" (H¹ vs H³ across thickness) is not made yet. Text S1 needs rewriting along the same lines.*

---

## 3.3. Scaling of the in-slab shear force gradient

Besides the mantle pressure difference, the main term balancing slab buoyancy is the gradient of the in-slab shear force, dQ_s/ds (Figures 3–4). To know when slab dip can be read as a pressure constraint, we need a scaling for this term in quantities observable on Earth.

The shear force is the integral of the shear stress across the slab, Q_s = ∫ τ_sn dn. In a viscous slab, to leading order,

τ_sn ≈ η (∂v_s/∂n + ω),   ω = ∂v_n/∂s − K v_s,   (6)

where ω is the rotation rate of the slab midsurface. If material cross-sections translate without rotating (∂v_s/∂n = 0), the slab is sheared through its thickness and dQ_s/ds scales with H. If cross-sections instead remain normal to the midsurface, as in thin viscous sheet theory (e.g., Ribe, 2001, 2010; Buffett, 2006), then ∂v_s/∂n = −ω. The two terms then cancel, and the shear force is set by the bending moment. Our slabs follow the second behaviour in every model: cross-sections rotate with the midsurface, the shear stress within the slab is an order of magnitude below the through-thickness estimate ηω, and the along-slab stress varies linearly across the slab (Figure S-bend a, b).

The shear force therefore comes from moment balance across the slab:

Q_s = dM/ds + C,   M = ∫ n σ_ss dn,   (7)

where M is the bending moment and C is the couple exerted by mantle shear tractions on the slab faces. In the models, C is small and Q_s is carried almost entirely by dM/ds (Figure S-bend c). For a viscous slab in plane strain,

M = η H³ κ̇ / 3,   (8)

where κ̇ = ∂²v_s/∂s∂n is the bending rate. This relation holds closely in the models (Figure S-bend d). Neglecting C and along-slab variations in η and H, Equations 7 and 8 give

dQ_s/ds ≈ d²M/ds² = (η H³ / 3) d²κ̇/ds².   (9)

For a slab moving through a curved geometry, ω ≈ −K v_s, so κ̇ ≈ ∂(K v_s)/∂s, and Equation 9 depends on the third along-slab derivative of K v_s. That derivative cannot be measured reliably on Earth, and is dominated by noise even in the models. We therefore represent it with a single length scale L_K over which K v_s varies along the slab, and take v_s ≈ v_c:

dQ_s/ds ≈ η H³ K v_c / (3 L_K³).   (10)

L_K is fitted to the models; we use L_K = 200 km (Supporting Figure S[LK]), comparable to the length over which curvature varies along the modelled slabs.

The H³ dependence is what distinguishes bending from through-thickness shear, and it can only be tested by varying plate thickness. At a single thickness, a scaling linear in H with a fitted length fits equally well. Across our H = 70, 100 and 130 km models, only the cubic form brings all thicknesses onto one trend (Supporting Figure S-varH). The thinnest slabs carry almost no shear force gradient, as Equation 10 predicts.

Figure 5b compares the modelled dQ_s/ds with Equation 10. The model means follow the 1:1 line closely ([r = 0.94]). Individual outputs scatter more, as expected when one length scale stands in for the along-slab structure of each slab. Overturned slabs follow the scaling to ~15 MPa and then saturate, which we attribute to the breakdown of the assumed geometry at very high curvature (Section 3.2).

Dividing Equation 10 by the slab buoyancy B = ΔρgH defines

Λ = η H³ K v_c / (3 L_K³ B),   (11)

the size of the in-slab stress term relative to buoyancy. Where Λ < 0.1, the modelled ΔP matches B_slab to within 10 % of B in [~90 %] of cases, and this agreement breaks down quickly above that value (Figure 5c). We therefore use Λ < 0.1 as the condition for reading slab dip as a measure of mantle pressure. Because B itself scales with H, Λ scales with H², so thin plates satisfy this condition more readily than thick ones (Section 4.2).
