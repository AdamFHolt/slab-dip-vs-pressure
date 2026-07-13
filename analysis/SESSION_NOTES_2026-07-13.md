# Session Notes (2026-07-13) — Wiki's reply accepted Λ; log-x Fig. 6

Follow-up to `SESSION_NOTES_2026-07-10.md` (Royden review response).

## 1. Wiki's reply (email, 2026-07-13)

She **accepted** keeping ηHKv_c/L_eff ("I take your point on your parameter
actually working!") and the dimensionless Λ compromise, with one change:
**remove cos θ from the normalization**, i.e.

    Λ = (ηHKv_c/L_eff) / B,   B = Δρ g H     (not B_slab = B cos θ)

She suggests defining B = Δρ g H in the paper alongside the slab-normal
B_slab, for clarity. Zoom: Tues/Wed late afternoon (2026-07-14/15).

Notes on her version:
- Her algebra cancels H → Λ = ηKv_c/(L_eff Δρ g). Exact in the models
  (H = 100 km in both numerator and buoyancy) but NOT on Earth: scaling uses
  mechanical H_eff (0.9·T_m isotherm) while buoyancy uses the full
  plate+crust density integral (Text S2). Keep the general ratio form.
- Λ_wiki = Λ_old · cos(dip). Recoverable from the existing
  `observations/text_files/maps.*.txt` without re-extraction: file rows
  equal the filtered segment order of `data/segment_data.txt`, and the
  filter (dip_shall, age<250, vc, K present) is independent of the swept
  parameters, so per-row dip can be row-matched.
- Model-side Fig. 6c normalization: in models B is a constant (≈49 MPa),
  so normalizing by B leaves r identical to the unnormalized values
  (0.841/0.886) — the 0.917 quoted on 07-10 was specific to B_slab.
- **Λ implementation still NOT done** (Adam handling manuscript edits
  himself). The 07-10 plan §1(d) stands but amended: denominator B, not
  B_slab → thresholds (0.15/0.3) need re-deriving before touching the
  observation scripts. A comparison script was drafted at
  scratchpad `lambda_cos_test.py` (session tmp, not run).

## 2. Her follow-up: log x-axis for Fig. 6b/c — DONE

Her point: on linear axes the many small-scaling, small-misfit points pile
up at the origin and get covered, making the fit look worse than it is.
Requires |K| (she pre-approved: "should be OK?").

New **permanent, separate** script `plot_forces.no-ot.logx.py`
(commit `4b191028`, pushed; original `plot_forces.no-ot.py` untouched):
- Panels b/c: x = (ηH|K|V_C)/L_eff on log axis, 1e-2..1e2 MPa; gridlines
  at decades only (minor ticks kept); 1:1 line retained (curves on log-x);
  panel a (misfit vs K, linear) unchanged.
- Outputs `plots/DP-comparisons/compilations/*.no-ot.logx.png/pdf`
  (gitignored, like all plot outputs). Run:
  `python3 plot_forces.no-ot.logx.py 300000 10000 10000 1000`.
- Data ranges (tmin=3): normal models |scaling| 0.013–13 MPa (N=187),
  overturned 6–57 MPa (N=76). 63/187 normal points (34%) have K<0 and fold
  onto |K| — display-only change, r values unaffected.
- Preference: log version for main text (shows the 3.5-decade dynamic
  range, the flat noise-floor below ~1 MPa = sensor regime, and cleanly
  separates the overturned cluster). Caveats for caption: folded negative-K
  points lose sign-agreement information; deviations from 1:1 harder to
  judge by eye (r carries that). Linear version could go to supplement.

## 3. Open items (carried + amended)
1. Λ rework in observation scripts — now with B (no cos θ); re-derive
   thresholds ≈ equivalent to 5/10 MPa, then implement per 07-10 plan §1(d).
2. Manuscript: notation pass (τ_s, σ_n=τ_n+P, Q_s/Q_n), pressure
   derivation, Fig. 1 relabel, define B and B_slab — Adam editing docx.
3. Data Availability URL mismatch (Slabs-as-pressure-probes vs
   slab-dip-vs-pressure) — still unresolved.
4. `plots/tmp/` artifacts (test_royden_scaling.py, fig5_LM_split.py, ...)
   still in gitignored dir.
5. Zoom with Wiki Tues/Wed (2026-07-14/15) late afternoon.
