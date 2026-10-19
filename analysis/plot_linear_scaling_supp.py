"""
Supplement figure: Fig. 5b redrawn twice, a) for a scaling LINEAR in plate
thickness, dQ_s/ds ~ eta*H*K*v_c/L, with the length L fitted to the H = 100 km
models, and b) for the cubic moment-balance form of Fig. 5b.  Per-thickness
trend lines (slope through the origin of the model means) show that the linear
form puts H = 130 km on twice the H = 100 km slope, while the cubic form brings
them onto one trend.  H = 70 km gets no line (no resolvable gradient).  Same model selection, trimming, marker
styles and panel extent as Fig. 5b (plot_forces.varH.py ... fig5 means inset
lamd Hsize strip): every output pale, model means with +-1 sigma whiskers,
marker size = H, overturned hollow.

Prints, for the caption, the fitted L, and the per-H slopes through the origin
of the model means under the linear scaling and under the cubic one (L_K).

Usage:  python3 plot_linear_scaling_supp.py [z_m] [L_km]
        z_m  analysis depth (default 300.e3); L_km overrides the fitted length
Writes plots/DP-comparisons/compilations/manuscript_figs/H-scaling_supp.z<d>...png/pdf
"""
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

from functions import nominal_plate_thickness

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})

analysis_depth = float(sys.argv[1]) if len(sys.argv) > 1 else 300.e3
L_override = float(sys.argv[2]) if len(sys.argv) > 2 else None
analysis_depth_dz, ds, dz = 10.e3, 10.e3, 1.e3
L_K = 200.                                    # km, the cubic scaling's length (Fig. 5b)
tactual_min = 11
tmin = tactual_min - 8
suffix = '.z%s.shear-dz%s.ds%s.prof-dz%skm' % (analysis_depth/1.e3, analysis_depth_dz/1.e3, ds/1.e3, dz/1.e3)

K_ind = 11; sn_ind = 6; vc_ind = 19; cov_ind = 38
mant_visc  = 2.5e20
cmyr_to_ms = 0.01/(365.25*24*3600)

# --- model selection: Fig. 5's (alleta, nodrho71, overturned hollow, fixed-SP 50 trimmed)
H100 = {'free':    {50: '2D_compositional_subd_lower-res_new_50plates',
                    250: '2D_compositional_subd_lower-res_new_250plates',
                    375: '2D_compositional_subd_lower-res_new_375plates',
                    500: '2D_compositional_subd_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_1000plates'},
        'FixedSP': {50: '2D_compositional_subd_lower-res_new_FixedSP_50plates',
                    250: '2D_compositional_subd_lower-res_new_FixedSP_250plates',
                    375: '2D_compositional_subd_lower-res_new_FixedSP_375plates',
                    500: '2D_compositional_subd_FixedSP_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_FixedSP_1000plates2'},
        'FixedOP': {50: '2D_compositional_subd_lower-res_new_FixedOP_50plates',
                    250: '2D_compositional_subd_lower-res_new_FixedOP_250plates',
                    375: '2D_compositional_subd_lower-res_new_FixedOP_375plates',
                    500: '2D_compositional_subd_FixedOP_lower-res_new',
                    1000: '2D_compositional_subd_lower-res_new_FixedOP_1000plates'}}

def name(bc, eta, H):
    if H == 100:
        return H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)

BCS  = ['free', 'FixedSP', 'FixedOP']
ETAS = [50, 250, 375, 500, 1000]
HS   = [70, 100, 130]
bc_label  = {'free': 'free', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^'}
colors    = {50: 'tan', 250: 'peru', 375: 'firebrick', 500: 'maroon', 1000: 'black'}
MEAN_S    = {70: 16, 100: 30, 130: 50}
INSET_S   = {70: 5, 100: 9, 130: 15}
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70),
              ('FixedOP', 375, 100)}
LAST_OUTPUT = {('FixedSP', 50, 100): 22}

def load(bc, eta, H):
    m = name(bc, eta, H)
    f = 'text_files/withT_varH/' + m + suffix + '.txt'
    if not os.path.exists(f):
        if not (eta in (50, 375) and H != 100):
            print('skipping %s (no withT_varH file)' % m)
        return None
    d = np.loadtxt(f)[tmin:]
    d = d[d[:, cov_ind] == 1]
    d = d[d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9)]
    Hn = nominal_plate_thickness(m)
    base = d[:, K_ind] * eta * mant_visc * d[:, vc_ind] * cmyr_to_ms * 1e-6      # eta*K*v_c  [MPa/km]
    return dict(lin=Hn * base,                                # eta*H*K*v_c, divide by L [km] for MPa
                cub=(Hn/L_K)**3 / 3 * base,                   # eta*H^3*K*v_c/(3 L_K^3)  [MPa]
                dqds=-d[:, sn_ind] / 1.e6, Hn=Hn)             # MPa, paper sign

data = {(bc, e, H): load(bc, e, H) for bc in BCS for e in ETAS for H in HS}
data = {k: v for k, v in data.items() if v is not None}
normal = {k: d for k, d in data.items() if k not in OVERTURNED}

# --- the linear length, fitted to the normal H = 100 km models (every output,
# slope through the origin, as the L_K fit): L = sum(lin*dqds)/sum(dqds^2)... no:
# dqds = lin/L  ->  1/L = sum(lin*dqds)/sum(lin^2)
sel = [d for k, d in normal.items() if k[2] == 100]
lin = np.concatenate([d['lin'] for d in sel]); dq = np.concatenate([d['dqds'] for d in sel])
L_fit = np.sum(lin**2) / np.sum(lin*dq)
L_lin = L_override if L_override else L_fit
print('linear length fitted to H = 100 km normal models: L = %.0f km (%d models, N = %d);  using L = %.0f km'
      % (L_fit, len(sel), len(lin), L_lin))
for d in data.values():
    d['scal'] = d['lin'] / L_lin

# --- per-H slopes of the model means (and every output), both scalings, for the caption
def slopes(label, key):
    print(label)
    for H in HS + ['all']:
        ks = [k for k in normal if H == 'all' or k[2] == H]
        xm = np.array([np.mean(normal[k][key]) for k in ks]); ym = np.array([np.mean(normal[k]['dqds']) for k in ks])
        x = np.concatenate([normal[k][key] for k in ks]); y = np.concatenate([normal[k]['dqds'] for k in ks])
        print('  H = %-4s %2d models, N = %3d: means slope %.2f (r %.2f, rms %.2f MPa);  outputs slope %.2f (r %.2f)'
              % (H, len(ks), len(x), np.sum(xm*ym)/np.sum(xm**2), pearsonr(xm, ym)[0] if len(xm) > 2 else np.nan,
                 np.sqrt(np.mean((ym - xm)**2)), np.sum(x*y)/np.sum(x**2), pearsonr(x, y)[0]))
slopes('linear scaling, L = %.0f km (slope = 1 at H = 100 by construction unless L overridden):' % L_lin, 'scal')
slopes('cubic scaling, L_K = %.0f km (Fig. 5b):' % L_K, 'cub')

# --- figure: two panels in Fig. 5b's style, a) the linear scaling, b) the cubic
# one (Fig. 5b's), same per-thickness trend lines in both; legend to the right.
# H = 70 km gets no line: its outputs show no slope at all (free fit
# -0.03 +- 0.10, r = -0.03; means scatter -2..+1 MPa), so a line through the
# origin would invent a trend the data do not hold
NO_LINE = {70}
lims = (-4, 16, -4, 16)
fig = plt.figure(figsize=(6.2, 2.6))
axA = fig.add_axes([0.075, 0.16, 0.33, 0.78])
axB = fig.add_axes([0.455, 0.16, 0.33, 0.78])

def panel(ax, key, xlab, lab):
    x0, x1, y0, y1 = lims
    for (bc, e, H), d in data.items():
        m, ot = bc_marker[bc], (bc, e, H) in OVERTURNED
        x, y = d[key], d['dqds']
        if ot:
            ax.scatter(x, y, s=5, facecolors='none', edgecolor=colors[e], marker=m, linewidth=0.3, alpha=0.35, zorder=2)
        else:
            ax.scatter(x, y, s=4, color=colors[e], marker=m, linewidth=0, alpha=0.45, zorder=2)
        xm, ym = np.mean(x), np.mean(y)
        if not (x0 <= xm <= x1 and y0 <= ym <= y1):
            continue
        ax.errorbar(xm, ym, xerr=np.std(x), yerr=np.std(y), fmt='none', ecolor=colors[e], elinewidth=0.4,
                    alpha=0.5, capsize=0, zorder=3)
        if ot:
            ax.scatter(xm, ym, s=MEAN_S[H], facecolors='none', edgecolor=colors[e], marker=m, linewidth=0.75,
                       zorder=4, clip_on=False)
        else:
            ax.scatter(xm, ym, s=MEAN_S[H], color=colors[e], edgecolor='black', marker=m, linewidth=0.4,
                       zorder=5, clip_on=False)
    ax.plot([x0, x1], [x0, x1], color='gray', linestyle='--', linewidth=1, zorder=1)
    ax.set_xticks([-4, 0, 4, 8, 12, 16]); ax.set_yticks([-4, 0, 4, 8, 12, 16])
    ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
    ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1, adjustable='box')
    # per-thickness trend: slope through the origin of the normal model means,
    # drawn behind the points and labelled, so the fan of slopes (or its
    # absence) is what the eye reads
    H_slope = {}
    for H in HS:
        ks = [k for k in normal if k[2] == H]
        xm = np.array([np.mean(normal[k][key]) for k in ks]); ym = np.array([np.mean(normal[k]['dqds']) for k in ks])
        H_slope[H] = sl = np.sum(xm*ym) / np.sum(xm**2)
        if H in NO_LINE:
            continue
        xe = min(x1, y1/sl) if sl > 0 else x1            # where the line leaves the axes
        ax.plot([0, xe], [0, sl*xe], color='steelblue', linewidth=1.1, zorder=1.5)
        ang = np.degrees(np.arctan(sl * (x1 - x0)/(y1 - y0)))
        if sl*xe >= y1:
            # leaves through the top: label above the axis (right-aligned near the corner)
            near = xe > x1 - 2
            ax.text(x1 if near else xe, y1 + 0.25, 'H = %d km' % H, color='steelblue', fontsize=7,
                    ha='right' if near else 'center', va='bottom', clip_on=False, zorder=10)
        else:
            # leaves through the right: label just below the line, along it
            ax.text(11.0, sl*11.0 - 0.35, 'H = %d km' % H, color='steelblue', fontsize=7, rotation=ang,
                    rotation_mode='anchor', ha='center', va='top', zorder=10)
    print('%s) %s: trend slopes through the origin, model means: ' % (lab, key)
          + ', '.join('H %d: %.2f' % (H, H_slope[H]) for H in HS))
    # misfit of the normal model means about the 1:1 line, lower right (the
    # points leave that corner empty): correlation, rms, and the H slopes
    xm = np.array([np.mean(d[key]) for k, d in normal.items()]); ym = np.array([np.mean(d['dqds']) for k, d in normal.items()])
    r_, rms_ = pearsonr(xm, ym)[0], np.sqrt(np.mean((ym - xm)**2))
    print('%s) model means: r %.2f, rms about 1:1 %.2f MPa' % (lab, r_, rms_))
    ax.text(0.97, 0.04, 'model means:\nr = %.2f, rms = %.1f MPa\nslopes H%s: %s' % (r_, rms_,
            ' / '.join('%d' % H for H in HS if H not in NO_LINE),
            ' / '.join('%.2f' % H_slope[H] for H in HS if H not in NO_LINE)),
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5, color='black', zorder=10,
            linespacing=1.25, bbox=dict(facecolor='white', edgecolor='none', pad=1.2))
    ax.set_xlabel(xlab)
    ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
            bbox=dict(facecolor='white', edgecolor='none', pad=1))

panel(axA, 'scal', r'$\eta H K v_c / L$   [MPa]', 'a')
panel(axB, 'cub', r'$\eta H^3 K v_c / (3 L_K^3)$   [MPa]', 'b')
axA.set_ylabel(r'$-dQ_{s}/ds$   [MPa]')
axB.set_yticklabels([])
ax = axB          # the legend hangs off panel b)

# legend: one vertical column just right of the panel, each group under a bold header
header = lambda t: Line2D([], [], linestyle='', marker='', label=t)
eta_h = [Line2D([], [], marker='o', color=colors[e], markeredgecolor='black', markeredgewidth=0.25,
                linestyle='', markersize=4, label='%d' % e) for e in ETAS]
bc_h  = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
        [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray', linestyle='',
                markersize=4, label='overturned')]
H_h   = [Line2D([], [], marker='o', color='gray', markeredgecolor='black', markeredgewidth=0.25, linestyle='',
                markersize=np.sqrt(MEAN_S[h]) * 0.75, label='%d km' % h) for h in HS]
hs = [header("viscosity ratio \u03b7':")] + eta_h + [header('plates:')] + bc_h + [header('thickness:')] + H_h
leg = ax.legend(handles=hs, ncol=1, loc='center left', bbox_to_anchor=(1.04, 0.5), fontsize=6.5,
                frameon=True, fancybox=False, edgecolor='black', facecolor='white', framealpha=1,
                borderpad=0.6, labelspacing=0.35, handletextpad=0.1, handlelength=1.0, borderaxespad=0)
leg.get_frame().set_linewidth(0.5)
for h_, t in zip(hs, leg.get_texts()):
    if h_.get_marker() in ('', 'None', None):
        t.set_fontweight('bold')
# the headers have an empty handle slot; shift them left so they start where
# the markers do
fig.canvas.draw()
r_ = fig.canvas.get_renderer()
hx0 = min(h.get_window_extent(r_).x0 for h in leg.legend_handles if h.get_window_extent(r_).width > 0)
for h_, t in zip(hs, leg.get_texts()):
    if h_.get_marker() in ('', 'None', None):
        dx = (t.get_window_extent(r_).x0 - hx0) * 72 / fig.dpi
        t.set_transform(mpl.transforms.offset_copy(t.get_transform(), fig=fig, x=-dx, y=0, units='points'))

out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
os.makedirs(out_dir, exist_ok=True)
plot_name = out_dir + 'H-scaling_supp' + suffix + '.tmin%d.L%d' % (tmin, L_lin)
plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/.pdf')
