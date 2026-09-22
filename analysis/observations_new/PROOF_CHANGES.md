# GJI-26-0707 — changes to make when the proof comes back

Cause: the trench segmentation discarded the leftover <300 km at the end of every
PB2002 boundary string, leaving 421 km unsampled at the Japan/Izu-Bonin cusp (and
17% of global trench length). Fixed in this folder; the pipeline
has been re-run end to end. Details: `README_SEGMENTATION_FIX.md`.

**N goes from 80 to 85 segments. Mean ΔP is unchanged (31.7 MPa); the qualifying
fraction improves, 86.2% → 89.4%.**

## Status (2026-09-21)

- **Figures — DONE.** All new figures swapped into the Illustrator files.
  Nothing shifted but the data layers, as expected (every color scale, axis
  limit, basemap and legend in these scripts is hardcoded).
- **Text edits — outstanding.** The table below, plus the checks at the end.
- **The proof is not back yet.** Page and line numbers below are from the
  submitted manuscript and will not match the proof's pagination.

## Text edits

| page | ms line | current | change to |
|---|---|---|---|
| 3 | 13 | "for 86% of segments" | **89%** |
| 3 | 14–15 | "ΔP ranges from 11 to 63 MPa (mean ~32 MPa)" | **7 to 62 MPa** (mean ~32 unchanged) |
| 13 | 303 | "92% of segments have K < 0.002 km⁻¹" | **94%** — 80 of 85, was 74 of 80 (92.5%) |
| 18 | 432 | "segments of 500 km trench-parallel length" | **"segments of nominally 300 km"** — see note 1 |
| 20 | 478 | "Λ spans ~4 × 10⁻⁴ to 0.5, median of 0.05" | **~5 × 10⁻⁴ to 0.45** (median 0.05 unchanged) |
| 20 | 479 | "86.2% ... (69 of 80) have Λ < 0.1 and 95.0% have Λ < 0.2" | **89.4% (76 of 85) ... and 94.1%** |
| 20 | 485 | "the 11 excluded segments" | **9** (same places: N Tonga, Ryukyu near Taiwan, south-central Chile) |
| 21 | 498 | "range from 11 to 63 MPa, with a mean of 31.7 MPa" | **7 to 62 MPa**, mean 31.7 unchanged |
| 22 | Fig 7 caption | "over the 69 segments that satisfy Λ < 0.1" | **76** |
| 23 | 546 | "47.5% of segments still meet ... (η = 8 × 10²²)" | **50.6%** |
| 23 | 548 | "rises to almost all segments (95.0%)" (η = 2 × 10²²) | **94.1%** |
| 26 | 650 | "86% at our reference slab viscosity" | **89%** |
| 38 | Fig S3 caption | "K < 0.002 km⁻¹ (92% of segments)" | **94%** — same 80 of 85 |
| 44 | Fig S9 caption (not S8; corrected 2026-09-22 against the docx) | "over the 69 segments that satisfy Λ < 0.1" | **76** |

## Verified 2026-09-22 against `~/Downloads/Holt_etal_GJI_manuscript.docx`

Every row above was recomputed from `observations_new/text_files` and
`observations/text_files` (script kept in the session scratchpad; numbers agree with
the table to the quoted precision). The Fig S3 caption already reads 94% in the docx;
the main-text "92%" (Sec 3.2) does not. Two pre-existing inaccuracies, present with the
OLD data too and independent of the segmentation fix:

- Sec 4.3 "mean ΔP values between ~28 and 38 MPa (Figure 7b)": the fixed-mask field over
  the plotted window (1250–1450 °C, 80–135 km) is 28.5–40.1 MPa, both old and new; 40 MPa
  is reached only in the 1450 °C / 135 km corner. Suggest "~28 and 40 MPa".
- Sec 4.3 "change the mean ΔP value by < 25%": that corner is +27% (40.1/31.7); the α–κ
  window is −17%..+22% and the eclogite term −14%. Suggest "by up to ~25%" or "< 30%".

"at most 0.6 MPa" (Fig 7 and Fig S9 captions) still holds: new max |fixed − reselected|
is 0.53 MPa in the Fig 7b window, 0.21 in the S9a window, 0.52 for the no-crust contours.
Λ range is 5.5e-4..0.453 (write ~5 × 10⁻⁴ or ~6 × 10⁻⁴), median 0.050. Eclogite removal
is 4.4 MPa (was 4.35), still "4 MPa".

No change needed: mean ΔP 31.7 MPa (p. 24, line 599); the Fig 7b range "~28 to 38 MPa";
"~70% of segments" using 300 km curvature (p. 19, line 456 — now 73%); "exceed 60 MPa
beneath the Northeast Japan slab" (new max 62.4 MPa, still NE Japan, and the three
highest ΔP segments are now all Japan Trench); Data Availability repo name is correct.

## Figures to replace — done

All regenerated in `plots/`, all swapped in:

- **Fig 6** and **Fig 7** — `maps.slab4e+22...pdf`. Note the two right-column
  insets differ: the T-vs-plate-thickness contour panel (the "Fig 7b range
  ~28 to 38 MPa") is visually identical old vs new — frozen-mask averaging,
  max cell change 0.029 MPa against 2 MPa contour intervals, same as Fig S9.
  The %-vs-viscosity scatter beneath it **does** change and must be replaced:
  the reference marker rises 86.2% → 89.4%, the 8e22 point 47.5% → 50.6%, the
  2e22 point 95.0% → 94.1%. Do not treat the two insets as one unit.
- **Fig S3** — `maps_H-curvature-vc.pdf`
- **Fig S6** — `maps_dip-age.pdf`
- **Fig S7** — `maps_scaling-B.pdf`
- **Fig S8** — `just-maps.slab2e+22...pdf` and `just-maps.slab8e+22...pdf`

## Figures that do NOT need replacing — left as submitted

- **Fig S9** (mean ΔP as a function of parameters) — `DP-param-exploration.pdf`.
  Earlier drafts of this file listed it as a "Fig 7b sensitivity panel" and put it
  in the replace list. Both were wrong. It is Fig S9, and its artwork is unchanged:
  `DPmean_fixed` averages over a mask frozen at the reference parameters, so the
  segmentation fix only shifts the field by at most

  | panel | field range, submitted → new | max cell change |
  |---|---|---|
  | T vs plate thickness, crust = 7 km | 24.3–40.1 → 24.3–40.1 MPa | 0.029 MPa |
  | T vs plate thickness, crust = 0 km | 20.0–35.8 → 20.0–35.8 MPa | 0.018 MPa |
  | α vs κ | 24.4–44.4 → 24.4–44.5 MPa | 0.118 MPa |

  Contour levels are 1 MPa (dashed) and 2 MPa (solid), so the worst case moves a
  line by ~12% of a contour interval — less than its width. The two rendered
  versions are visually indistinguishable. Keep the submitted artwork.

  **But check the S9 caption.** The population it averages over went from 69 to
  76 segments. If the caption quotes that N, or says "the 69 segments that
  satisfy Λ < 0.1", the caption still needs 69 → 76 even though the figure does
  not get redrawn.

## Two notes

1. **The "500 km" on line 432 is wrong independently of this fix.** The code has
   used 300 km segments in every version (checked the git history), so that
   sentence has never matched the analysis. Worth correcting either way.
2. **Two new segments fill the cusp**, both qualifying: 142.22°E, 35.34°N
   (Λ = 0.031, ΔP = 61.1 MPa) and 142.00°E, 32.96°N (Λ = 0.053, ΔP = 49.2 MPa).
   The new ΔP minimum of 7.2 MPa is a new Mexican-trench segment (101.7°W, 17.0°N),
   which is what moves the lower end of the range from 11 to 7 MPa.

## Still to check in the proof

- **Fig S9 caption** — the artwork is unchanged, but the population it averages
  over went 69 → 76. If the caption quotes that N, it still needs editing.
- **Stray numbers** — grep for any other bare "92" (K < 0.002) or "69"
  (Λ < 0.1 population) beyond the rows in the table above.
- **Supplement figure numbering** — the labels above (S3, S6, S7, S8) come from
  the submitted manuscript. One entry in this file was already wrong: it had
  `DP-param-exploration.pdf` as a main-text "Fig 7b sensitivity panel" when it
  is Fig S9. Confirm the rest, in particular that the Fig S3 and Fig S8 caption
  rows point at the K/v_c/H_eff map and the 2e22/8e22 map pair respectively.
