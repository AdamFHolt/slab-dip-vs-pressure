---
title: "Response to reviewers — GJI-26-0707"
subtitle: "Working draft, started 2026-10-22"
geometry: margin=1in
fontsize: 11pt
---

*Notes (delete before sending): replies are added as they are drafted. Reviewer comments are quoted in italics, our replies in plain text.*

---

## Reviewer 1

*(to be drafted)*

---

## Reviewer 2

### Comment 1 — thermal diffusion
> *The first is related to neglecting thermal diffusion. The assumptions around this were stated very clearly, but I feel like I would have benefited from a little more justification around this. Looking at the length of the model runs (~20 Myr) and the reference thermal diffusivity used (8.0e-7 m^2/s), slab temperatures could diffuse over length scales of ~20 km. This would definitely change the buoyancy force predicted from the models. For me a simple enough test to see the impact of ignoring thermal diffusion would be to re-run a single well behaved model with a non-zero diffusivity and see what this does to the force balance presented in Figure 4.*

We thank the reviewer for raising this. Our original description did not make clear that the slab buoyancy in our models is compositional, not thermal. The models are isothermal: temperature is a uniform 1573 K at the start and on every boundary. Density depends only on composition, with the lithosphere 50 kg m⁻³ denser than the mantle, and viscosity is set by composition and depth. There is therefore no temperature gradient for diffusion to act on. A rerun with κ = 8 × 10⁻⁷ m² s⁻¹ would give the same solution, and the Figure 4 force balance would be unchanged. We have revised Section 2.2 to say this explicitly and to write out the governing equations, which also responds to Reviewer 1.

The reviewer's underlying point still applies to Earth's slabs, where buoyancy is thermal and diffusion acts on it over the timescales considered. We do not think it changes our conclusions, for three reasons:

1. Diffusion redistributes the slab's thermal anomaly but does not remove it: heat gained by the slab is lost by the adjacent mantle, so the depth-integrated negative buoyancy, B, is conserved. What diffusion changes is how that buoyancy is distributed across the slab.
2. The force balance (Eq. 4) is diagnostic. B, ΔP and dQ_s/ds are all measured from the model state at the time of analysis. The balance does not depend on how that state was reached, so it would hold for a thermally diffusing slab, just with different values of each term.
3. The main effect of diffusion that matters for our argument is a change in effective slab thickness. Its importance is set by Λ through the H² dependence (Eq. 13). We test that dependence directly with the new variable-thickness models (H = 70–130 km), whose range is larger than the ~20 km diffusion length the reviewer estimates. In the Earth application, slab buoyancy and effective thickness come from a plate-cooling model of the subducting plate's age, so thermal structure is included there (Section 4).

**Change to the manuscript (Section 2.2):**

> The models are isothermal (T = 1573 K throughout), so slab buoyancy is purely compositional and thermal diffusion plays no role; the effect of slab thickness, which diffusion would modify on Earth, is tested directly with models of H = 70, 100 and 130 km (Section 3.3).

### Comment 2 — plate thickness

> *The second concern is that the models seem to use the same subducting plate thickness (100 km). However, I would imagine that if the thickness of the subducting plate were smaller, that this would change slab stresses, potentially increasing the slab stress term in the force balance equation? This would be important to note because, as the authors present in the discussion, subduction zones globally have difference ages and thicknesses, so this could impact where the model results can be extrapolated. I think running a couple models with a thinner/thicker subducting plate and showing how this influences the force balance would be great.*

We thank the reviewer for this suggestion, which we have followed. We have run 18 new models with subducting plates 70 and 130 km thick, alongside the original 100 km plates: slab viscosity ratios η′ = 250, 500 and 1000, each with free plates, a fixed subducting plate and a fixed overriding plate. They are analysed in exactly the same way as the original models, and all of the model figures now include them: the model evolution (Figure 2), the force balance (Figures 3 and 4) and the scaling of the in-slab shear force gradient (Figure 5), together with the corresponding supplementary figures.

The force balance holds at all three thicknesses (Figure 3). The new models also showed that our original scaling for the in-slab shear force gradient, which was linear in plate thickness, was not the right form: with one length scale it fits the 100 km models but not the thinner and thicker ones. The gradient instead scales with H³, as expected for viscous bending of the slab rather than shear across its thickness. We have rewritten Section 3.3 to derive the scaling from the moment balance of thin viscous sheet theory, with the derivation and the shear alternative in Supporting Text S1, and compare the two forms across the thickness models in a new supplementary figure (Figure S-varH). In answer to the reviewer's question, a thinner plate therefore has a smaller in-slab stress term, not a larger one. Λ is now defined with this scaling, so it varies as H² (relative to buoyancy, which is itself proportional to H), and thin plates satisfy Λ < 0.1 more readily than thick ones. We have recomputed the application to Earth's subduction zones (Section 4) with the new scaling. Because Earth's plates are mostly thinner than 100 km, many more segments now satisfy the criterion: at the reference slab viscosity, 98 % of segments have Λ < 0.1, compared with 89 % with the original scaling, and the criterion holds for most segments up to a slab viscosity of ~10²³ Pa s.

### Comment 3 — shallow slabs
> *My third concern is that the authors choose to extract information from their numerical models once the slab tip passes 600 km depth, which they then apply to subduction zones globally. However, not all subduction zones have slabs that reach 600 km depth (Cascadia, Scotia, Aleutians, Central America). I think this could be an issue because the numerical models used in this study have a high viscosity lower mantle at 660 km depth. As the slab tip approaches the lower mantle, the high viscosity layer could start to support the subducting plate, which I imagine would change the stress distribution within the slab and also the mantle wedge. At the aforementioned subduction zones, which only get down to depths of 300-400 km or so, this effect would not be present. Some additional discussion on this, and perhaps testing what happens if you re-run the same force balance calculation on your model output but once the slab reaches a depth of 400-500 km, would help to clarify this.*

We agree that this is an important consideration. We would note that at several of these margins tomographic models image the slab continuing to the base of the transition zone (e.g., the compilation of Goes et al., 2017, their Figure 2, which includes Cascadia, the Aleutians and Central America), although the deeper geometry is less certain at some of them, including Cascadia. Interaction with the lower mantle may therefore be relevant to more of these slabs than their upper-mantle extent suggests.

The reviewer is right that the slabs in our reference analysis are supported viscously by the lower mantle: they are all extracted after the tip has reached it. To test whether the force balance depends on that support, we repeated the calculation, at the same 300 km depth, for earlier outputs in which the slab tip has not yet reached the lower mantle (new Figure S-early). The forces balance in all of them. Where the tip lies more than about two plate thicknesses below the analysis depth, the results are the same as with lower-mantle support: ΔP tracks B_slab to within a few percent of B, and Λ < 0.1 still identifies the outputs where the two agree. The relationship between slab dip and mantle pressure therefore holds whether or not the slab is supported by the lower mantle.

The models do show one limit that bears on the reviewer's point. Where the slab tip lies less than about two plate thicknesses below the depth of analysis, the in-slab stresses become a significant part of the balance (~10 % of B), so ΔP no longer equals B_slab even though the forces still balance. Near its end the slab's shear force must fall to zero, so the in-slab stresses there are set by the distance to the tip rather than by the local curvature, and Λ does not always capture them. Slab dip therefore records the mantle pressure difference where the slab continues at least ~2 plate thicknesses below the depth at which it is measured, whether or not it is supported by the lower mantle. For Earth's plates (H ≈ 75 km) and dips measured at 180–320 km depth, this requires the slab to reach ~450–500 km, as the tomography suggests it does at the margins the reviewer lists. We state this condition in the revised Discussion.

**Change to the manuscript:** In Section 4 (Discussion), after the application to global subduction zones, we have added: "In the models, slab dip records the mantle pressure difference whether or not the slab is supported by the lower mantle, provided the slab continues at least ~2 plate thicknesses below the depth at which dip is measured (Figure S-early); closer to the slab tip, in-slab stresses take up part of the buoyancy. Several slabs that appear to end at 300–400 km depth, including Cascadia, the Aleutians and Central America, are imaged tomographically to the base of the transition zone (Goes et al., 2017)." *[Placement and wording to check against the revised Section 4.]*

**New supplementary figure (caption draft):** Figure S-early. Force balance before the slab reaches the lower mantle, at 300 km depth. Coloured symbols: model outputs taken before the slab tip reaches the lower mantle, coloured by the distance of the tip below the analysis depth in plate thicknesses (brown: less than 2; blue: more than 2); grey: the outputs of Figures 3 and 5. Symbol shape gives the plate configuration; open symbols are overturned slabs. (a) ΔP + σ_slab against B_slab (as Figure 3); (b) ΔP against B_slab (as Figure 5a). Fourteen early outputs lie outside the axes in (a) and four in (b).

**Reference added:** Goes, S., Agrusta, R., van Hunen, J., & Garel, F. (2017). Subduction-transition zone interaction: A review. *Geosphere*, 13(3), 644–664. https://doi.org/10.1130/GES01476.1



