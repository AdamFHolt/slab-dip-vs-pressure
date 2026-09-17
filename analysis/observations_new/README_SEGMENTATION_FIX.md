# observations_new — trench segmentation fix

Working copy of `observations/` with one change: how PB2002 boundary strings are
split into segments. `observations/` is left exactly as submitted.

## The problem

`extract_map_properties.py` segmented each PB2002 boundary string with

    num_segments = int(total_length // segment_length)   # 300 km

The floor division discarded the remainder at the end of **every** string — up to
300 km of trench per string end, 8,890 km of the 52,390 km of PB2002 subduction
boundary (17%) in total.

It is most visible at the Japan/Izu-Bonin cusp, because PB2002 breaks the trench
at the Boso triple junction (141.88degE, 34.21degN) and two discarded tails meet there:

| string | what it is | length | segments | discarded tail |
|---|---|---|---|---|
| `PA\OK` | Japan Trench, running S to the junction | 793 km | 2 | 193 km |
| `PS/PA` | Izu-Bonin, running N to the junction | 1128 km | 3 | 228 km |

That left 421 km of unsampled trench from ~35.8degN to ~32.2degN — a visible hole in
the DP and Lambda maps next to the reported DP maximum beneath NE Japan.

## The change

    if total_length < min_string_length:   # 150 km
        continue
    num_segments = max(1, int(round(total_length / segment_length)))
    seg_spacing  = total_length / num_segments

Segments are now spaced `total_length/num_segments`, so the whole string is
sampled. Segment length is nominal (~250-350 km instead of exactly 300 km);
it never enters any calculation — it only places the sample centers, so the
physics is untouched.

Verified: no boundary string loses segments (145 -> 174 centers). Seven strings
of 159-229 km get their first ever sample, and all seven are removed anyway by
the existing age/v_c/K quality filters, so `min_string_length` does not affect
any result.

## What changed in the numbers (reference eta = 4e22 Pa s)

| quantity | submitted | new |
|---|---|---|
| segments passing the quality filters | 80 | **85** |
| Lambda < 0.1 | 86.2% (69) | **89.4% (76)** |
| Lambda < 0.2 | 95.0% | **94.1%** |
| mean DP, Lambda < 0.1 | 31.71 MPa | **31.72 MPa** |
| DP range, Lambda < 0.1 | 10.6-62.7 MPa | **7.2-62.4 MPa** |
| Lambda < 0.1 at eta = 2e22 | 95.0% | **94.1%** |
| Lambda < 0.1 at eta = 8e22 | 47.5% | **50.6%** |

The two segments that fill the cusp:

| center | Lambda | DP |
|---|---|---|
| 142.22degE, 35.34degN (Japan Trench) | 0.031 — qualifies | 61.1 MPa |
| 142.00degE, 32.96degN (Izu-Bonin)    | 0.053 — qualifies | 49.2 MPa |

The new DP minimum, 7.2 MPa, is a new segment on the Mexican trench
(101.68degW, 17.01degN); the previous minimum was 10.6 MPa.

## Layout

- Scripts are copies of `observations/` as of commit 0f76801d.
- `data/Slab2` and `data/Muller2008` are **hardlinks** to the originals (no extra
  disk, inputs are read-only). `data/Lallemand` and `data/PB2002` are real copies.
- `data/segment_data.OLD.txt` is the submitted segment file, kept for diffing.
- `papers/`, `old/`, `for_Tao/` were not copied — reference material, not pipeline.

## How to reproduce

    cd observations_new
    # 1. segment file (~1.5 min)
    MPLBACKEND=Agg python3 extract_map_properties.py
    # 2. parameter sweep, 5838 files (~5 min on 24 cores; many_extractions.sh runs
    #    the same grid serially, ~2 h)
    # 3. figures and stats (~2.5 min)
    bash all_plots.sh

`python3` must be an interpreter with Basemap **and** cmcrameri — i.e.
`~/miniconda3/envs/mantle-flow-modeling/bin/python` *with* user site-packages
enabled. Do not set `PYTHONNOUSERSITE=1`: cmcrameri lives in `~/.local`.
