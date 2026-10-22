"""
Reviewer 2, comment 3: does each Earth segment's slab continue at least ~2 plate thicknesses below the
depth at which its curvature is measured?  In the models, dip records the mantle pressure difference
only there (Figure S-early).  For every segment of data/segment_data.txt this walks the same down-dip
profile as extract_map_properties.py (and four more at +-10, +-20 deg), extended to 2500 km from the trench (1000 km there, which
reaches the curvature depths but not the deepest slab at shallow-dipping margins), through the Slab2 depth
grids, takes the deepest Slab2 depth on it, and compares it with z_meas + 2 H_eff, where z_meas is
300 km when the deep curvature (280-320 km) is used and 200 km when only the shallow one (180-220 km)
is, and H_eff is the reference plate-cooling thickness.  Slab2 follows seismicity, so its deepest
depth is a lower bound on how far the slab extends.

Usage:  python3 slab2_depth_check.py      (writes text_files/slab2_depth_check.txt)
"""
import glob
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator
from functions import read_pb2002_boundaries, haversine, calculate_bearing, destination_point
from functions import make_strictly_ascending, compute_H_eff

TM, DIFFUSIVITY, PLATE_THICK = 1333.0, 8.044e-7, 88000.0      # reference values, as print_stats.py

grids = []
for f in glob.glob("data/Slab2/*_dep_*grd"):
    nc = Dataset(f, 'r')
    lon, lat, z = make_strictly_ascending(nc.variables['x'][:]), make_strictly_ascending(nc.variables['y'][:]), nc.variables['z'][:]
    nc.close()
    grids.append(dict(interp=RegularGridInterpolator((lat, lon), z, bounds_error=False, fill_value=np.nan),
                      lon0=lon.min(), lon1=lon.max(), lat0=lat.min(), lat1=lat.max(), name=f.split('/')[-1][:3]))

# segment centres and down-dip bearings, exactly as extract_map_properties.py
segs = []
for seg in read_pb2002_boundaries("data/PB2002/PB2002_boundaries.dig.txt"):
    if seg['symbol'] not in ["/", "\\"] or len(seg['coords']) < 2:
        continue
    coords = seg['coords']
    cumdist = [0]
    for i in range(1, len(coords)):
        cumdist.append(cumdist[-1] + haversine(*coords[i-1], *coords[i]))
    if cumdist[-1] < 150:
        continue
    n = max(1, int(round(cumdist[-1] / 300)))
    for i in range(n):
        c = i * cumdist[-1] / n + cumdist[-1] / n / 2.
        for j in range(len(cumdist) - 1):
            if cumdist[j] <= c <= cumdist[j+1]:
                fr = (c - cumdist[j]) / (cumdist[j+1] - cumdist[j])
                lonc = coords[j][0] + fr * (coords[j+1][0] - coords[j][0])
                latc = coords[j][1] + fr * (coords[j+1][1] - coords[j][1])
                b = calculate_bearing(coords[j][0], coords[j][1], coords[j+1][0], coords[j+1][1])
                segs.append((lonc, latc, (b - 90) % 360 if seg['symbol'] == "/" else (b + 90) % 360))

sd = np.genfromtxt("data/segment_data.txt", delimiter=',', dtype=str)
sd[sd == 'None'] = np.nan
sd = sd.astype(float)
assert len(segs) == len(sd) and np.allclose(np.array(segs)[:, :2], sd[:, :2]), 'segment order differs'

out = []
for (lonc, latc, perp), r in zip(segs, sd):
    shallow_K, deep_K, vc, age, dip_shall = r[2], r[3], r[4], r[5], r[9]
    K = deep_K if not np.isnan(deep_K) else shallow_K
    if not (not np.isnan(dip_shall) and age < 250 and not np.isnan(vc) and not np.isnan(K)):
        continue                       # the same selection as the Earth Lambda (N = 85)
    zmax, name = 0., ''
    # a fan of profiles, +-20 deg about down-dip: a single profile can run past the deepest part of a
    # slab whose strike curves with depth (Tonga at 18 S: 425 km on the profile, 590 km in the grid)
    for brg in perp + np.arange(-20, 21, 10):
        for d in np.linspace(0, 2500, 1250):
            lo, la = destination_point(lonc, latc, brg % 360, d)
            lo = lo % 360
            for g in grids:
                if g['lon0'] <= lo <= g['lon1'] and g['lat0'] <= la <= g['lat1']:
                    z = -float(g['interp']((la, lo)))
                    if np.isfinite(z) and z > zmax:
                        zmax, name = z, g['name']
    z_meas = 300. if not np.isnan(deep_K) else 200.
    H = compute_H_eff(age, TM, k=DIFFUSIVITY, plate_thick=PLATE_THICK) / 1e3
    out.append((lonc, latc, age, H, z_meas, zmax, (zmax - z_meas) / H, name))

out_arr = np.array([o[:7] for o in out])
np.savetxt('text_files/slab2_depth_check.txt', out_arr, fmt='%.3f',
           header='lon lat age_Ma H_eff_km z_meas_km slab2_max_depth_km (zmax-z_meas)/H_eff')
d = out_arr[:, 6]
print('segments: %d   z_meas 300: %d, 200: %d' % (len(d), np.sum(out_arr[:, 4] == 300), np.sum(out_arr[:, 4] == 200)))
for thr in (1, 1.5, 2, 2.5, 3):
    print('  slab extends >= %.1f H_eff below z_meas (Slab2 alone): %d of %d (%.0f %%)' % (thr, np.sum(d >= thr), len(d), 100*np.mean(d >= thr)))
print('\nsegments short of 2 H_eff (Slab2 zone, lon, lat, age, H_eff, z_meas, Slab2 max depth, distance in H):')
for o in sorted(out, key=lambda o: o[6]):
    if o[6] < 2:
        print('  %s  %7.1f %6.1f  %5.1f Ma  H %3.0f  z %3.0f  max %4.0f km  %.2f H' % (o[7], o[0], o[1], o[2], o[3], o[4], o[5], o[6]))
