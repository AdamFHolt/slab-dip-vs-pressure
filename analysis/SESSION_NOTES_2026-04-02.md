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
- Eq. 8 (thin-sheet theory) fails because L_K ≈ H violates thin-sheet assumption;
  L_eff framing sidesteps this without needing to explain the failure mode in detail.

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

## Script
`plot_dqds_diagnostic.py` — self-contained, auto-loads `text_files/TESTC/*.txt`

## Next Session Tasks
1. **Update paper text**:
   - Reframe Eq. 11 in terms of L_eff = 1963 km (not a dimensionless prefactor)
   - Justify: Eq. 8 fails (L_K ≈ H), L_eff empirically ~constant, K stays in
   - Note L_eff scatter (IQR) and what drives it
2. **Update Earth application figures**:
   - Use `dQ/ds / DP_anal = H K v_c / (L_eff · Δρ g cos θ)` with L_eff = 1963 km
   - Existing Earth-application script exists — find it and update coefficient
   - Inputs needed per subduction zone: H, K, v_c, Δρ, θ
