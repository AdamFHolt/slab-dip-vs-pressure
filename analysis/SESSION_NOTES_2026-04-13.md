# Session Notes (2026-04-13)

## Summary
Discussed the dQ/ds scaling narrative for paper rewriting, clarified the role of
each panel in Fig. 6, updated plot_forces.no-ot.py and plot_forces.supp-shear.no-ot.py
to use the new L_eff scaling.

---

## dQ/ds recap: why L_eff and not separate terms

### The two-term structure (Eq. 8)
```
dQ/ds ≈ η H [ 2K·dvs/ds  +  vs·dK/ds ]
```
These two terms vary on independent length scales L_v and L_K. The old Eq. 8→9
step collapsed them into a common L — that step is unjustified.

### Why not evaluate each term separately?
Two compounding problems:
1. **dK/ds is not reliably computable** — it is a 3rd derivative of raw slab
   geometry. Noise dominates even from model output. Restricting to low KH
   (where thin-sheet is most valid) makes correlation *worse*, not better —
   confirming the failure is noise-limited, not thin-sheet breakdown.
2. **dvs/ds and L_v are unobservable on Earth** — no equivalent of slab-parallel
   velocity is available for real subduction zones.

### Why L_eff works
Both terms scale as ηHKvs/L, so their sum is proportional to ηHKvc regardless
of the individual L_v and L_K. L_eff absorbs whatever combination is operating.
Its stability across the suite (median 1624 km, IQR relatively tight for normal
models) means the ratio L_v:L_K is not wildly variable — the sum behaves as if
there is one effective scale. This collapses the problem to H, K, vc — all
estimable for real subduction zones.

---

## Why high-K points fail in Fig. 2 panel 1 (ηHKvc/L_eff vs measured)

High predicted values (x > ~10 MPa) fall below the 1:1 line — formula overpredicts.
These are essentially the overturned models. Two reinforcing reasons:

1. **L_eff calibrated on normal models only (1624 km)** — overturned models have
   larger actual L_eff (right tail in fig1), so applying 1624 km overpredicts them.
2. **Geometric breakdown at high K** — thin-sheet assumptions (uniform τ_sn across H,
   small metric corrections) become less valid at high curvature. The formula
   inherits this since it is derived from Eq. 8. Confirmed by fig1 panel 2: high-K
   points sit above the OLS line (L_eff grows with K at the high end), though not
   dramatically.

The high-K outliers are almost exclusively the overturned models — not high-K
normal models. This means the geometric breakdown and "not observed on Earth"
exclusion criteria are not independent: overturned models are both physically
unrealistic AND in the regime where the thin-sheet assumptions are least valid.

**For the paper:** the exclusion of overturned models is doubly justified — on
physical grounds (not observed on Earth) and on geometric grounds (thin-sheet
framework breaks down at their curvatures). Both arguments should appear.

---

## Fig. 6 (plot_forces.no-ot.py): panel roles and narrative

Three panels, now on one row:

### Panel 1: Full stress misfit vs K
- Plots (DP + σ_slab − B_slab) vs K for all models
- This uses the **numerically computed** force balance terms, not the scaling
- Normal models cluster near zero across all K — force balance closes
- Overturned (hollow) models sit off-zero even here — **the force balance
  framework itself struggles at high K**, independent of any scaling approximation
- Key implication: the problem with overturned models is not a failure of the
  L_eff scaling — it is more fundamental. Even with exact numerical terms, the
  balance doesn't close for these geometries.
- Also motivates why K matters: residual correlates with K

### Panel 2: DP misfit vs ηHKvc/L_eff
- Shows (DP − B_slab) as a function of the scaling parameter
- Application panel: once you know dQ/ds from the scaling, you can assess which
  slabs have DP close to B_slab (within ~10 MPa) and which don't
- The scatter band directly quantifies uncertainty in DP reconstruction for Earth

### Panel 3: dQ/ds vs ηHKvc/L_eff
- The calibration panel — direct test of the scaling law
- 1:1 line shown; r=0.87 (normal models)

### Narrative flow for paper (suggested order: 1 → 3 → 2)
1. Force balance closes numerically (panel 1) — and overturned cases fail even here
2. The shear term is well predicted by the scaling (panel 3) — calibration
3. Therefore we can assess DP vs B_slab for any slab using the scaling (panel 2)

---

## Script changes (2026-04-13)

### plot_forces.no-ot.py
- `coeff` changed from `100./1000.0` (= 0.1) to `100./1624.0` (= H[km]/L_eff[km])
  → implements ηHKvc/L_eff with H=100 km, L_eff=1624 km
- Layout changed from GridSpec(2,3) to GridSpec(1,3) — all panels on one row
- x-axis tightened to −5→40 MPa on panels 2 and 3 (was −15→65)
- Panels reordered: full stress misfit vs K | dQ/ds vs scaling | DP misfit vs scaling
  (narrative: balance closes → calibrate scaling → apply to DP)
- 1:1 line added to DP misfit panel
- Axis labels updated to `(ηHKV_C)/L_eff [MPa]`

### plot_forces.supp-shear.no-ot.py
- Same `coeff` change: `100./1624.0`
- Layout changed from GridSpec(2,3) to GridSpec(1,3) — dropped repeated combined-scaling
  panel (already in main figure); now 3 panels: dQ/ds vs K, vs Vc, vs η
- Removed unused imports for `plot_forcecomponent_dqds` and `plot_forcecomponent_dqds_overturned`

### plot_DPvsDP.color-points.no-ot.py (bottom row)
- Color variable updated from `0.1·ηKvc` to `(100/1624)·ηKvc` = `ηHKvc/L_eff`
  — both inline scatter and via `plot_BvsDP_scalingcolored` / `_overturned`
- colorbar label updated to `(ηHKV_C)/L_eff [MPa]`
- color_max2 kept at 20 MPa

### functions_plotting.py
- `plot_BvsDP_scalingcolored`: `0.1 *` → `(100./1624.) *`
- `plot_BvsDP_scalingcolored_overturned`: same

---

## Next session tasks (updated)

1. **Fill in [X–Y] km IQR** for L_eff from fig1 panel 1 (normal models only)

2. **Update Eq. 11 everywhere in the paper**:
   - Was: dQ/ds ≈ 0.1ηKvc
   - Now: dQ/ds ≈ ηHKvc / L_eff,  L_eff = 1624 km (excl. overturned)
   - Check: abstract, Section 3.3, Section 4.1, conclusions, Fig. 6 caption,
     Fig. 7 caption, any place 0.1ηKvc appears

3. **Reorder Fig. 6 panels** (if not already): 1 → 3 → 2 for logical narrative flow

4. **Update Earth application figures**:
   - Use dQ/ds / DP_anal = HKvc / (L_eff · Δρ g cosθ) with L_eff = 1624 km
   - Find existing Earth-application script (in analysis/observations/ ?)

5. **Update Fig. 7 (Earth application map)**:
   - Update to use ηHKvc/L_eff with L_eff = 1624 km

6. **Paper text — justification for excluding overturned models**:
   - Add that these models also fail the *exact* numerical force balance (panel 1),
     not just the scaling — this is independent evidence they are geometrically
     anomalous, not just empirical outliers
   - Do NOT say "thin-sheet assumption breaks down" without qualification —
     it is valid for normal models; overturned cases are a separate geometric regime
