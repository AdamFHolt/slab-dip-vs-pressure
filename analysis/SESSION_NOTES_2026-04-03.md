# Session Notes (2026-04-03)

## Summary
Rewrote dQ/ds diagnostic, reframed scaling law around L_eff, and worked out the
correct narrative for the manuscript derivation section.

---

## Key Result: L_eff framing

The scaling law is:

    dQ/ds  =  η H K v_c / L_eff

where L_eff is the effective length scale implied by the data:

    L_eff  =  η H K v_c / |dQ/ds_meas|

**Calibrated values:**
- L_eff = 1624 km  (median, N=119, excl. overturned slabs)  ← use for Earth application
- L_eff = 1963 km  (median, N=183, all models)

**Overturned models excluded** (not observed on Earth):
  η'=1000 free plates, η'=1000 fixed OP, η'=500 fixed OP

This replaces the previous 0.1ηKvc prefactor framing.

**Earth application form:**
    dQ/ds / DP_anal  =  H K v_c / (L_eff · Δρ g cos θ)

---

## Where the derivation fails and why

### The thin-sheet step (Eq. 6) is valid
The thin-sheet assumption is used at Eq. 6:
    dQ/ds ≈ H·dτ_sn/ds ≈ ηH·d²vn/ds²
This integrates τ_sn across H assuming it is uniform in the slab-normal direction.
It is valid when: H ≪ L (thin sheet), K²H² ≪ 1, dvs/dn ≈ 0.
This step gives Eq. 8 as the thin-sheet result.

### Eq. 8 has TWO independent length scales
    dQ/ds ≈ ηH [ 2K·dvs/ds  +  vs·dK/ds ]
               \_____________/  \__________/
               scales as K·vs/L_v   K·vs/L_K

The old Eq. 8→9 step collapsed L_v and L_K into a single common L — unjustified.

### Why direct evaluation of Eq. 8 fails — primary cause: numerical noise
Direct computation of ηH·vs·dK/ds from model outputs gives:
  - All same-sign:  N=122, r=0.71
  - KH < 0.2:       N=116, r=0.80
  - KH < 0.1:       N=100, r=0.22  ← much worse
  - KH < 0.05:      N=84,  r=0.10  ← no correlation

Counterintuitively, restricting to cases where KH ≪ 1 (where thin-sheet is most
valid) makes the prediction WORSE, not better. The reason: at low K, dQ/ds itself
is small and below the noise floor. The signal in the dK/ds prediction comes from
the higher-K cases where dQ/ds is large enough to see. At low K, numerical noise
in dK/ds dominates entirely.

**Conclusion: numerical noise in dK/ds is the primary limitation, not thin-sheet
breakdown.** The thin-sheet assumptions are valid for the retained models (KH < 0.2,
K² small) — they just don't help when the required derivatives can't be computed
reliably.

### Why dK/ds is noisy
dK/ds requires differentiating K, which is itself computed by differentiating a
Savitzky-Golay-smoothed slab contour. This is effectively a 3rd derivative of the
raw slab geometry, which amplifies noise substantially. Additionally, dK/ds changes
sign depending on whether the slab is increasing or decreasing in curvature at that
depth — producing wrong-sign predictions ~half the time.

### Revised manuscript narrative (replacing earlier "L_K ≈ H" story)
The correct statement is: Eq. 8 cannot be evaluated directly because dK/ds is not
reliably estimable — not because the thin-sheet approximation breaks down.

---

## Replacement manuscript text (Section 3.3, after Eq. 8)

Replace "To obtain a first-order scaling, we further assume that vs and K vary over
a common mantle length scale L..." through old Eq. 11 (and its justification) with:

---

Equation 8 contains two terms, involving the along-slab gradients of vs and K
respectively, which vary on different and independently evolving length scales.
Both are difficult to constrain reliably — model-derived dK/ds in particular
requires repeated numerical differentiation of the slab geometry, and direct
evaluation of Eq. 8 yields large overpredictions of dQ/ds. We tested whether
restricting to cases where thin-sheet assumptions are most strictly satisfied
(KH ≪ 1) improves the prediction, but the correlation actually decreases at lower
KH values — indicating that the failure is driven by numerical noise in dK/ds
rather than a breakdown of the thin-sheet approximation. At low curvature, dQ/ds
is sufficiently small that numerical noise in dK/ds dominates.

We therefore take an empirical approach. Both terms in Eq. 8 scale as ηHKvs/L,
where L represents the relevant along-slab length scale, giving:

    dQ/ds ∝ η H K vs                                                      (9)

Rather than estimating L analytically, we define the effective length scale
implied directly by the data:

    L_eff = η H K vc / dQ/ds                                              (10)

where we substitute vc ≈ vs as before. If L_eff is approximately constant across
our model suite, Eq. 9 can be written as a practical scaling law with a single
calibrated parameter:

    dQ/ds ≈ η H K vc / L_eff                                              (11)

Excluding models with rollover slab geometry not observed on Earth, we find a
median L_eff = 1624 km (N = 119; Fig. X), and Eq. 11 predicts measured dQ/ds
with r = 0.87 (Fig. 6c). The scatter in L_eff reflects the fact that Eq. 11
approximates the full two-term expression of Eq. 8 with a single empirical
constant — the residuals encode real differences in the along-slab variation of
vs and K between models that are not captured by η, H, K, and vc alone. The IQR
of L_eff ([X–Y] km; read from fig1 panel 1, normal models only) directly
quantifies the uncertainty on dQ/ds estimates when Eq. 11 is applied to Earth slabs.

---

**Do NOT say:**
- "H/L_eff ≈ 1/25 is consistent with expected H/L suppression" — circular
- "L_K ≈ H causes the failure" — likely a numerical artifact
- "the thin-sheet assumption breaks down" without qualification — it is valid for
  the retained models; the issue is that dK/ds cannot be computed reliably

---

## Does K belong in the scaling law?
Yes. Removing K (no-K formula ηHvc/α) degrades r from 0.87 → 0.69 (excl. overturned).
K must stay in the formula.

---

## Current figures (plots/tmp/)

### fig1: dqds_diag.fig1_Leff_distribution.{png,pdf}
- Panel 1: L_eff histogram — filled=normal (N=119, median=1624 km),
           outlined=overturned (N=64, median pulls overall to 1963 km)
- Panel 2: L_eff vs K scatter, OLS fit through origin on normal models only
           hollow = overturned
- Panel 3: per-model strip, ordered fixedSP/free/fixedOP × η'=50/250/375/500/1000
           dashed bars = overturned models; two median lines (normal + all)

### fig2: dqds_diag.fig2_Leff_prediction.{png,pdf}
- Panel 1: ηHKvc/L_eff vs measured — excl.overturned r=0.87, all r=0.86
- Panel 2: ηHvc/α vs measured (no K) — worse than panel 1
- Panel 3: ηH·vs·dK/ds vs measured — r=0.71 all, r=0.80 (KH<0.2), r=0.22 (KH<0.1)
           → dK/ds not a viable predictor; noise-limited

### Script
`plot_dqds_diagnostic.py` — self-contained, auto-loads text_files/TESTC/*.txt
OVERTURNED_MODELS and MODEL_ORDER defined at top of script.

---

## Next Session Tasks

1. **Fill in [X–Y] km IQR** for L_eff from fig1 panel 1 (normal models only)

2. **Update Eq. 11 everywhere in the paper**:
   - Was: dQ/ds ≈ 0.1ηKvc
   - Now: dQ/ds ≈ ηHKvc / L_eff,  L_eff = 1624 km (excl. overturned)
   - Check: abstract, Section 3.3, Section 4.1, conclusions, Fig. 6c caption,
     Fig. 7 caption, any place 0.1ηKvc appears

3. **Update Earth application figures**:
   - Use dQ/ds / DP_anal = HKvc / (L_eff · Δρ g cosθ) with L_eff = 1624 km
   - Find existing Earth-application script (in analysis/observations/ ?)
   - Inputs already in the obs scripts: H, K, vc, Δρ, θ per subduction zone

4. **Update Fig. 6c** (currently shows 0.1ηKvc vs measured dQ/ds):
   - Replot as ηHKvc/L_eff vs measured, with L_eff = 1624 km
   - Add 1:1 line; this becomes the key calibration figure

5. **Fig. 7 (Earth application map)**:
   - Currently uses 0.1ηKvc as the shear stress scaling parameter
   - Update to use ηHKvc/L_eff with L_eff = 1624 km
   - Note H needs to be estimated per subduction zone (plate model)
