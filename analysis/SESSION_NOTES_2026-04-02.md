# Session Notes (2026-04-02)

## Summary
Rewrote dQ/ds diagnostic entirely, consolidated to 2 publication-quality figures,
and reframed the scaling law around L_eff.

## Key Result: L_eff framing

The scaling law is now expressed as:

    dQ/ds  =  η H K v_c / L_eff

where L_eff is the effective length scale implied by the data:

    L_eff  =  η H K v_c / |dQ/ds_meas|

**Calibrated value: L_eff (median) = 1963 km  (N=183)**

This replaces the previous dimensionless-prefactor framing (e.g. ~0.035 × ηKvcH).

## Why L_eff is more useful than a prefactor
- More physical: has units of length, can be compared to slab geometry
- Directly portable to Earth application:
    dQ/ds / DP_anal  =  H K v_c / (L_eff · Δρ g cos θ)

## Where the derivation fails and why (important for manuscript)

### The thin-sheet step (Eq. 6) is fine
The thin-sheet assumption is used at Eq. 6:
    dQ/ds ≈ H·dτ_sn/ds ≈ ηH·d²vn/ds²
This integrates τ_sn across H assuming it is uniform in the slab-normal direction.
This step is valid and gives Eq. 8 as the thin-sheet result — it is NOT the failure point.

### Eq. 8 has TWO different length scales
    dQ/ds ≈ ηH [ 2K·dvs/ds  +  vs·dK/ds ]
               \_____________/  \__________/
               scales as K·vs/L_v   K·vs/L_K

The old derivation (Eq. 8 → 9) collapsed L_v and L_K into a single common L — an
unjustified assumption. These are two independent length scales.

### Why direct evaluation of Eq. 8 fails
Direct computation of ηH·vs·dK/ds from model outputs (col12) predicts dQ/ds with
r=0.75 but:
- Wrong sign ~half the timesteps (N drops from 247 → 122)
- Noisy magnitude (dK/ds = 3rd derivative of slab geometry after smoothing)
The failure is NOT a thin-sheet violation — Eq. 8 is the thin-sheet result.
It is a practical estimation problem: dK/ds is unreliable.

### Consequence for L_K ≈ H diagnosis (revision of earlier interpretation)
The previous claim that "L_K ≈ H violates thin-sheet assumption" was likely wrong.
The small L_K values seen in the models may be a numerical artifact of differentiating
already-smoothed K, producing spuriously large |dK/ds| → small L_K.
The physical curvature gradient length scale may be larger than H.

### Correct narrative for the manuscript
The failure is not theoretical (thin-sheet) — it is practical (dK/ds unreliable).
Eq. 8 cannot be evaluated directly for either models or real slabs.
Instead, note that both terms scale as ηHKvs/L, giving dQ/ds ∝ ηHKvs (Eq. 9),
and define L_eff empirically. L_eff ≈ constant is a non-trivial empirical finding.

### What NOT to say
- Do NOT say "H/L_eff ≈ 1/25 is consistent with expected H/L suppression" —
  this is circular (L_eff is defined from dQ/ds, so H/L_eff = dQ/ds/(ηKvc) by construction)
- Do NOT say the thin-sheet assumption breaks down at Eq. 8 → 9

## Does K belong in the scaling law?
Yes. Removing K (no-K formula ηHvc/α) degrades r from 0.86 → 0.75.
The middle panel of fig1 shows L_eff is only weakly proportional to K (r=0.51),
so K does not cancel cleanly — it belongs in the formula.

## Current figures (plots/tmp/)
- `dqds_diag.fig1_Leff_distribution.{png,pdf}`
    - Panel 1: L_eff histogram, median = 1963 km, N=183
    - Panel 2: L_eff vs K scatter with OLS-through-origin fit (slope = 1449×10⁶ m, r=0.51)
    - Panel 3: per-model strip sorted by median, IQR bars
- `dqds_diag.fig2_Leff_prediction.{png,pdf}`
    - Panel 1: ηHKvc/L_eff vs measured (all same-sign), r=0.86, N=247
    - Panel 2: ηHvc/α vs measured (no-K formula, α=1.47×10¹² m²), r=0.75, N=193
    - Panel 3: ηH·vs·dK/ds vs measured (direct Eq. 8 dominant term), r=0.75, N=122
      → good correlation but wrong sign ~half the time; not a viable predictor

## Script
`plot_dqds_diagnostic.py` — self-contained, auto-loads `text_files/TESTC/*.txt`

## Replacement manuscript text (Section 3.3, after Eq. 8)

Replace "To obtain a first-order scaling, we further assume that vs and K vary over
a common mantle length scale L..." through old Eq. 11 (and its justification) with:

---

Equation 8 contains two terms, involving the along-slab gradients of vs and K
respectively, which vary on different and independently evolving length scales.
Both are difficult to constrain reliably — model-derived dK/ds in particular is
noisy due to the repeated numerical differentiation required to compute it, and
direct evaluation of Eq. 8 yields large overpredictions of dQ/ds. We therefore
take an empirical approach. Both terms in Eq. 8 scale as ηHKvs/L, where L
represents the relevant along-slab length scale, giving:

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
of L_eff ([X–Y] km; compute from fig1 panel 1) directly quantifies the
uncertainty on dQ/ds estimates when Eq. 11 is applied to Earth slabs.

---

**Do NOT:**
- Say the thin-sheet assumption breaks down (it was already used to get Eq. 8)
- Say "L_K ≈ H" causes the failure (likely a numerical artifact)
- Say "H/L_eff ≈ 1/25 is consistent with expected suppression" (circular)

## Next Session Tasks
1. **Fill in [X–Y] km IQR** for L_eff from fig1 panel 1 histogram (normal models only)
2. **Update Eq. 11 coefficient** everywhere it appears in the paper (was 0.1ηKvc,
   now ηHKvc/L_eff with L_eff = 1624 km excl. overturned, 1963 km all models)
3. **Update Earth application figures**:
   - Use `dQ/ds / DP_anal = H K vc / (L_eff · Δρ g cos θ)` with L_eff = 1624 km
   - Existing Earth-application script exists — find it and update coefficient
   - Inputs per subduction zone: H, K, vc, Δρ, θ (already in the obs. scripts)
4. **Update Fig. 6c** (currently shows 0.1ηKvc vs measured) to show ηHKvc/L_eff
