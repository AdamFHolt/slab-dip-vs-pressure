# Session Notes (2026-07-10) — Royden review response

Context: Wiki (Leigh Royden) returned manuscript edits on 2026-07-09
(`~/Downloads/royden-emails.pdf`, edited docx attached to her email) just
before planned GJI submission. This session evaluated her comments,
especially the proposed replacement of the scaling parameter, and prepared
the response (email draft: `~/Downloads/email_to_wiki.txt`).

Manuscript: "Mantle slabs as dynamic pressure sensors" (Holt, Goldberg,
Zhao, Royden) — `~/Downloads/manuscript.pdf`, `supplement.pdf`.

---

## 1. Her comments and our positions

### (a) Pressure term (P_wedge + P_subslab) — ACCEPTED
She argued you cannot drop P_avg = (P_subslab+P_wedge)/2 by gauge argument
alone: the in-slab pressure acts on the slab-normal (edge) faces and must be
carried through the derivation, after which the neglected term becomes one
that is "readily argued near zero."

Verification that she is right: the analysis code integrates the
**deviatoric** normal stress for the K·N term (`extract_properties.py:275`,
`slab_norm_devstress_cut`; the `_fullstress` variant exists in
`functions.py:convert_to_slabnorm_normstress` but is not what goes into
col 17). For a closed segment, a uniform pressure shift must exert zero net
force, which requires the total stress on the edge faces. With deviatoric N,
the term actually neglected in Eq. 3→4 is the gauge-invariant
`HK(P̄_slab − P_avg)` — mean in-slab pressure minus mean adjacent mantle
pressure — small because pressure is nearly continuous across a thin slab,
and empirically small because the full force balance closes (Fig. 3).

**Consequences: derivation text only. No data or figure recomputation.**
Cosmetic: Fig. 1 sketch labels (add P on edge faces), captions with renamed
symbols.

### (b) Notation — ACCEPTED
- τ_ss/τ_ns → τ_s (shear on slab-normal faces), σ_n = τ_n + P (total normal).
- Q → Q_s, N → Q_n (ripples into "dQ_s/ds" in text + figure axes).
- K vs dθ/ds: keeping K (standard in bending literature: Ribe 2001,
  Buffett 2006); define clearly at first use.

### (c) Scaling parameter replacement — REJECTED, with model evidence
Her proposal: replace ηHKv_c/L_eff with η·H·v_s·(d²θ/ds²) and use the
dimensionless ratio buoyancy/(in-slab stress), claiming one-step derivation,
dimensionless, non-empirical.

Key recognition: her parameter **is already in the paper** — it is the
v_s·dK/ds term of Eq. 8 (DK/Dt for a steadily sinking slab). Direct test
against the model suite, where measured dQ/ds is ground truth
(`plots/tmp/test_royden_scaling.py`; TESTC z=300 km, shear-dz10, ds10,
non-overturned models, N=220 timesteps, N=105 with dQ/ds>0):

| predictor of measured −dQ/ds     | r (all) | r (dQ/ds>0) |
|----------------------------------|---------|-------------|
| ηHKv_c/L_eff (Eq. 10, L=1497 km) | 0.841   | 0.886       |
| ηH·v_s·dK/ds (Royden)            | 0.176   | 0.324       |
| ηH·2K·dv_s/ds (other Eq.8 term)  | 0.634   | 0.748       |
| full Eq. 8 (both terms)          | 0.403   | 0.571       |

Magnitude/sign: Royden predictor median 23× too large; wrong sign for 81%
of dQ/ds>0 timesteps (implied prefactor median −19, IQR −42 to −7). Not a
sign-convention artifact: magnitude-only r(|pred|,|meas|) = 0.48 vs 0.85.
Implied L_eff by contrast is stable: median 1497 km, IQR 907–2571
(IQR/median = 1.11).

Root cause: dK/ds is a 3rd derivative of slab shape → noise-dominated even
in clean model output with Sav-Gol smoothing. On Earth it is worse: our K
needs one derivative of Slab2 catalogued dips; d²θ/ds² needs two.

### (d) Dimensionless compromise — PROPOSED to Wiki (NOT yet implemented)
Keep ηHKv_c/L_eff as the calibrated dimensional parameter (calibration,
Fig. 6b, L_eff = 1497 km all unchanged); state the **Earth criterion** as

    Λ ≡ (ηHKv_c/L_eff) / B_slab        (resisting / driving stresses)

with Λ < 0.15 replacing σ < 5 MPa, and Λ < 0.3 replacing 10 MPa.
Λ = col1/col0 of existing `observations/text_files/maps.*.txt` — **no
re-extraction needed**.

Numbers (from `data/segment_data.txt`, reference thermal params):

| η [Pa·s] | σ<5 MPa | Λ<0.15 | mean ΔP (Λ<0.15) | set difference vs σ<5 |
|----------|---------|--------|------------------|------------------------|
| 2e22     | 95%     | 95%    | 31.6 MPa         | identical set           |
| 4e22     | 74%     | 75%    | 32.6 MPa         | 55 shared, 4 lost, 5 gained |
| 8e22     | 42%     | 38%    | 32.8 MPa         | similar                 |

Λ distribution at 4e22: median 0.096, p95 0.29 (so Λ<0.3 ≈ the 95%/10 MPa
result). Model-side check: normalizing Fig. 6c axes by B_slab changes r
0.841→0.832 (all), 0.886→**0.917** (calibration branch) — model dips max
82°, so cosθ never pathological.

Implementation plan when Wiki agrees (deliberately NOT done yet):
- `observations/plot_final_map.py`: `scaling_thresh=5` → `lambda_thresh=0.15`
  (line 34); reorder segment loop (compute ΔP before edge criterion, then
  `lam = |stress_scaling|/DP`); ax4 y-label → Λ<0.15; ax3/ax4 stats calls.
- `observations/functions.py`: `stats_data_file`, `stats_DP` criterion →
  `|scaling|/DP < thresh`.
- `observations/plot_final_just-maps.py`: same loop change (line ~148).
- `observations/print_stats.py`: add Λ<0.15 / Λ<0.3 rows.
- `observations/plot_param_exploration.py`: `thresh=5` → 0.15 (3 sites).
- Maps keep MPa colorbars (hybrid); only the qualifying criterion changes.
- Manuscript: new Λ equation after Eq. 10; Figs. 7b/8a captions; abstract
  stats 74%→~75%, mean ΔP 30.8→~32.6 MPa.

### (e) Fig. 5 upper-vs-lower-mantle split (her Jul 8 email; she then said
"hold off") — CHECKED PRIVATELY; no figure change warranted
Slab tip depths extracted from every `csv_outputs/*/full.t.csv` (deepest
llith>0.5 point; `plots/tmp/scan_tip_depths.sh` →
`plots/tmp/slab_tip_depths.txt`, 590 model-timesteps). "Reached LM" =
tip ≥ 660 km at any t' ≤ t (cumulative). Analysis:
`plots/tmp/fig5_LM_split.py` → `plots/tmp/fig5_LM_split.png`.

| group                | N   | r     | mean(ΔP−B) | mean\|ΔP−B\| |
|----------------------|-----|-------|------------|--------------|
| normal, UM only      | 36  | 0.732 | −1.2 MPa   | 2.2 MPa      |
| normal, reached LM   | 184 | 0.875 | −0.9 MPa   | 2.3 MPa      |
| normal, combined     | 220 | 0.867 | −0.9 MPa   | 2.3 MPa      |
| overturned (excl.)   | 88  | 0.249 | −8.5 MPa   | 8.9 MPa      |

Every normal model reaches 660 km within its first ~2–4 analyzed timesteps,
so "UM-only" is a small early transient with the same 1:1 behavior. The
genuinely different group is overturned vs normal — already shown hollow.
Zoom answer if raised: split separates a transient, not a regime; sub-1:1
excursions at intermediate B_slab are high-ηHKv_c/L_eff timesteps (Fig. 5b
coloring), not LM contact.

---

## 2. Housekeeping done this session
- Fixed stale `coeff = 100./1000.0` → `100./1497.0` at
  `plot_forces.supp-shear.no-ot.py:115` (fed only the unplotted
  combined-scaling column; regenerated figure confirmed unchanged).
  Commit `38f5d1d5`, pushed.
- Committed regenerated observation plots, retired zips, archived
  `old/plot_final_map_lower-scaling-thresh.py`. Commit `10647a89`, pushed.
- Email draft to Wiki: `~/Downloads/email_to_wiki.txt` (concedes pressure
  derivation + notation; holds on scaling with the r numbers; proposes Λ;
  Fig. 5 deliberately omitted).

## 3. Open items
1. Await Wiki's reply / Zoom (proposed Tues/Wed) → then implement Λ rework
   per plan in §1(d).
2. Manuscript notation pass: τ_s, σ_n=τ_n+P, Q_s/Q_n; adopt her pressure
   derivation; Fig. 1 relabel.
3. Data Availability URL mismatch: manuscript cites
   `github.com/AdamFHolt/Slabs-as-pressure-probes`; actual repo is
   `slab-dip-vs-pressure`. Fix or create before submission.
4. `plots/tmp/` artifacts (test_royden_scaling.py, fig5_LM_split.py,
   scan_tip_depths.sh, slab_tip_depths.txt) are in a gitignored dir —
   decide whether to move somewhere tracked.
5. `observations/for_Tao/tmp` still untracked (informal notes).

## 4. Environment note
Observation map scripts (Basemap) fail with user-site numpy 2.2; run with:
`PYTHONNOUSERSITE=1 ~/miniconda3/envs/mantle-flow-modeling/bin/python`
(numpy 1.26.4). `print_stats.py` and model-side scripts run fine with
system python3.
