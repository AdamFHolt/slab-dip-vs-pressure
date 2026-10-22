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

# Comment 1 — thermal diffusion
> *The first is related to neglecting thermal diffusion. The assumptions around this were stated very clearly, but I feel like I would have benefited from a little more justification around this. Looking at the length of the model runs (~20 Myr) and the reference thermal diffusivity used (8.0e-7 m^2/s), slab temperatures could diffuse over length scales of ~20 km. This would definitely change the buoyancy force predicted from the models. For me a simple enough test to see the impact of ignoring thermal diffusion would be to re-run a single well behaved model with a non-zero diffusivity and see what this does to the force balance presented in Figure 4.*

We thank the reviewer for raising this. Our original description did not make clear that the slab buoyancy in our models is compositional, not thermal. The models are isothermal: temperature is a uniform 1573 K at the start and on every boundary. Density depends only on composition, with the lithosphere 50 kg m⁻³ denser than the mantle, and viscosity is set by composition and depth. There is therefore no temperature gradient for diffusion to act on. A rerun with κ = 8 × 10⁻⁷ m² s⁻¹ would give the same solution, and the Figure 4 force balance would be unchanged. We have revised Section 2.2 to say this explicitly and to write out the governing equations, which also responds to Reviewer 1.

The reviewer's underlying point still applies to Earth's slabs, where buoyancy is thermal and diffusion acts on it over the timescales considered. We do not think it changes our conclusions, for three reasons:

1. Diffusion redistributes the slab's thermal anomaly but does not remove it: heat gained by the slab is lost by the adjacent mantle, so the depth-integrated negative buoyancy, B, is conserved. What diffusion changes is how that buoyancy is distributed across the slab.
2. The force balance (Eq. 4) is diagnostic. B, ΔP and dQ_s/ds are all measured from the model state at the time of analysis. The balance does not depend on how that state was reached, so it would hold for a thermally diffusing slab, just with different values of each term.
3. The main effect of diffusion that matters for our argument is a change in effective slab thickness. Its importance is set by Λ through the H² dependence (Eq. 13). We test that dependence directly with the new variable-thickness models (H = 70–130 km), whose range is larger than the ~20 km diffusion length the reviewer estimates. In the Earth application, slab buoyancy and effective thickness come from a plate-cooling model of the subducting plate's age, so thermal structure is included there (Section 4).

**Change to the manuscript (Section 2.2):**

> The models are isothermal (T = 1573 K throughout), so slab buoyancy is purely compositional and thermal diffusion plays no role; the effect of slab thickness, which diffusion would modify on Earth, is tested directly with models of H = 70, 100 and 130 km (Section 3.3).

### Comment 2 — plate thickness

*(to be drafted)*

### Comment 3 — shallow slabs
> *My third concern is that the authors choose to extract information from their numerical models once the slab tip passes 600 km depth, which they then apply to subduction zones globally. However, not all subduction zones have slabs that reach 600 km depth (Cascadia, Scotia, Aleutians, Central America). I think this could be an issue because the numerical models used in this study have a high viscosity lower mantle at 660 km depth. As the slab tip approaches the lower mantle, the high viscosity layer could start to support the subducting plate, which I imagine would change the stress distribution within the slab and also the mantle wedge. At the aforementioned subduction zones, which only get down to depths of 300-400 km or so, this effect would not be present. Some additional discussion on this, and perhaps testing what happens if you re-run the same force balance calculation on your model output but once the slab reaches a depth of 400-500 km, would help to clarify this.*



