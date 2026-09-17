# GJI-26-0707 — changes to make when the proof comes back

Cause: the trench segmentation discarded the leftover <300 km at the end of every
PB2002 boundary string, leaving 421 km unsampled at the Japan/Izu-Bonin cusp (and
17% of global trench length). Fixed in this folder; the pipeline
has been re-run end to end. Details: `README_SEGMENTATION_FIX.md`.

**N goes from 80 to 85 segments. Mean ΔP is unchanged (31.7 MPa); the qualifying
fraction improves, 86.2% → 89.4%.**

## Text edits

| page | ms line | current | change to |
|---|---|---|---|
| 3 | 13 | "for 86% of segments" | **89%** |
| 3 | 14–15 | "ΔP ranges from 11 to 63 MPa (mean ~32 MPa)" | **7 to 62 MPa** (mean ~32 unchanged) |
| 13 | 303 | "92% of segments have K < 0.002 km⁻¹" | **94%** |
| 18 | 432 | "segments of 500 km trench-parallel length" | **"segments of nominally 300 km"** — see note 1 |
| 20 | 478 | "Λ spans ~4 × 10⁻⁴ to 0.5, median of 0.05" | **~5 × 10⁻⁴ to 0.45** (median 0.05 unchanged) |
| 20 | 479 | "86.2% ... (69 of 80) have Λ < 0.1 and 95.0% have Λ < 0.2" | **89.4% (76 of 85) ... and 94.1%** |
| 20 | 485 | "the 11 excluded segments" | **9** (same places: N Tonga, Ryukyu near Taiwan, south-central Chile) |
| 21 | 498 | "range from 11 to 63 MPa, with a mean of 31.7 MPa" | **7 to 62 MPa**, mean 31.7 unchanged |
| 22 | Fig 7 caption | "over the 69 segments that satisfy Λ < 0.1" | **76** |
| 23 | 546 | "47.5% of segments still meet ... (η = 8 × 10²²)" | **50.6%** |
| 23 | 548 | "rises to almost all segments (95.0%)" (η = 2 × 10²²) | **94.1%** |
| 26 | 650 | "86% at our reference slab viscosity" | **89%** |
| 38 | Fig S3 caption | "K < 0.002 km⁻¹ (92% of segments)" | **94%** |
| 44 | Fig S8 caption | "over the 69 segments that satisfy Λ < 0.1" | **76** |

No change needed: mean ΔP 31.7 MPa (p. 24, line 599); the Fig 7b range "~28 to 38 MPa";
"~70% of segments" using 300 km curvature (p. 19, line 456 — now 73%); "exceed 60 MPa
beneath the Northeast Japan slab" (new max 62.4 MPa, still NE Japan, and the three
highest ΔP segments are now all Japan Trench); Data Availability repo name is correct.

## Figures to replace

All regenerated in `plots/`:

- **Fig 6** and **Fig 7** — `maps.slab4e+22...pdf`
- **Fig S3** — `maps_H-curvature-vc.pdf`
- **Fig S6** — `maps_dip-age.pdf`
- **Fig S7** — `maps_scaling-B.pdf`
- **Fig S8** — `just-maps.slab2e+22...pdf` and `just-maps.slab8e+22...pdf`
- **Fig 7b sensitivity panel** — `DP-param-exploration.pdf`

## Two notes

1. **The "500 km" on line 432 is wrong independently of this fix.** The code has
   used 300 km segments in every version (checked the git history), so that
   sentence has never matched the analysis. Worth correcting either way.
2. **Two new segments fill the cusp**, both qualifying: 142.22°E, 35.34°N
   (Λ = 0.031, ΔP = 61.1 MPa) and 142.00°E, 32.96°N (Λ = 0.053, ΔP = 49.2 MPa).
   The new ΔP minimum of 7.2 MPa is a new Mexican-trench segment (101.7°W, 17.0°N),
   which is what moves the lower end of the range from 11 to 7 MPa.
