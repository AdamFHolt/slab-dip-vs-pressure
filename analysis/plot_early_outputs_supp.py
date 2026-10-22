"""R2 comment 3, option 1: force balance and Fig. 5 metrics for the early outputs
(t = 8-10, tip ~470-650 km, not in the figures) against the figure outputs (t >= 11).
Same loading as plot_forces.varH.py (withT_varH, coverage flag, LAST_OUTPUT), n = 3, L_K = 200."""
import sys, numpy as np
from scipy.stats import pearsonr
sys.path.insert(0, '.')
from functions import nominal_plate_thickness

z = sys.argv[1] if len(sys.argv) > 1 else '250.0'
L_K = 200.
suffix = '.z%s.shear-dz10.0.ds10.0.prof-dz1.0km' % z
K_ind, ss_ind, sn_ind, anal_ind, DP_ind, vc_ind, cov_ind = 11, 17, 6, 4, 3, 19, 38
mant_visc, cmyr = 2.5e20, 0.01/(365.25*24*3600)
H100 = {'free': {50: 'lower-res_new_50plates', 250: 'lower-res_new_250plates', 375: 'lower-res_new_375plates',
                 500: 'lower-res_new2', 1000: 'lower-res_new_1000plates'},
        'FixedSP': {50: 'lower-res_new_FixedSP_50plates', 250: 'lower-res_new_FixedSP_250plates',
                    375: 'lower-res_new_FixedSP_375plates', 500: 'FixedSP_lower-res_new2',
                    1000: 'lower-res_new_FixedSP_1000plates2'},
        'FixedOP': {50: 'lower-res_new_FixedOP_50plates', 250: 'lower-res_new_FixedOP_250plates',
                    375: 'lower-res_new_FixedOP_375plates', 500: 'FixedOP_lower-res_new',
                    1000: 'lower-res_new_FixedOP_1000plates'}}
def name(bc, eta, H):
    if H == 100: return '2D_compositional_subd_' + H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130), ('FixedOP', 500, 100),
              ('FixedOP', 1000, 100), ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70), ('FixedOP', 375, 100)}
LAST_OUTPUT = {('FixedSP', 50, 100): 22}

tips = {}
for fn in ['plots/tmp/slab_tip_depths.txt', sys.argv[2]]:
    for line in open(fn):
        m, t, dd = line.split(); tips[(m, int(t))] = float(dd)

rows = []   # (key, t, tip, DP, Bslab, B, dqds, scal, full)
for bc in H100:
    for eta in [50, 250, 375, 500, 1000]:
        for H in [70, 100, 130]:
            if H != 100 and eta in (50, 375): continue
            m = name(bc, eta, H)
            try: d = np.loadtxt('text_files/withT_varH/' + m + suffix + '.txt')
            except OSError: print('missing', m); continue
            d = d[(d[:, cov_ind] == 1) & (d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9))]
            # extraction failures: midplane flipped (dip ~241 deg, negative thickness) -- fixed-SP H = 130
            # at t = 9-10, before the t >= 11 the figures use; dropped (2026-10-22)
            flipped = d[:, 5] > 180
            if flipped.any():
                print('dropped flipped-midplane outputs:', (bc, eta, H), d[flipped, 0].astype(int).tolist())
            d = d[~flipped]
            Hn = nominal_plate_thickness(m); B = 50*9.81*Hn*1e3*1e-6
            scal = (Hn/L_K)**3/3 * d[:, K_ind]*eta*mant_visc*d[:, vc_ind]*cmyr*1e-6
            full = (d[:, DP_ind] - d[:, anal_ind] + d[:, ss_ind] - d[:, sn_ind])/1e6
            for i in range(len(d)):
                t = int(d[i, 0])
                rows.append(((bc, eta, H), t, tips.get((m, t), np.nan), d[i, DP_ind]/1e6, d[i, anal_ind]/1e6,
                             B, -d[i, sn_ind]/1e6, scal[i], full[i]))

def summary(label, sel):
    r = [x for x in rows if sel(x)]
    if not r: return
    a = lambda j: np.array([x[j] for x in r])
    tip, DP, Bs, B, q, s, full = a(2), a(3), a(4), a(5), a(6), a(7), a(8)
    lam = np.abs(s)/B; mis = np.abs(Bs - DP)/B
    low = lam < 0.1
    print('%-34s N=%3d (%2d models) tip %3.0f-%3.0f | closure med|res|/B %.3f | DP~Bslab r %.2f, med|Bslab-DP|/B %.3f'
          ' | dQds~scal r %.2f slope %.2f | Lam<0.1: N=%3d, %3.0f%% within 10%%'
          % (label, len(r), len({x[0] for x in r}), np.nanmin(tip), np.nanmax(tip), np.median(np.abs(full)/B),
             pearsonr(Bs, DP)[0], np.median(mis), pearsonr(s, q)[0], np.sum(s*q)/np.sum(s*s),
             low.sum(), 100*np.mean(mis[low] < 0.1) if low.any() else np.nan))

print('z = %s km, n = 3, L_K = %.0f km' % (z, L_K))
for ot, otl in [(False, 'normal'), (True, 'overturned')]:
    summary('%s, t 8-10 (early)' % otl, lambda x: (x[0] in OVERTURNED) == ot and x[1] <= 10)
    summary('%s, t >= 11 (figures)' % otl, lambda x: (x[0] in OVERTURNED) == ot and x[1] >= 11)
    summary('%s, early, tip < 600 km' % otl, lambda x: (x[0] in OVERTURNED) == ot and x[1] <= 10 and x[2] < 600)
    summary('%s, early, tip < 550 km' % otl, lambda x: (x[0] in OVERTURNED) == ot and x[1] <= 10 and x[2] < 550)

print('\nnormal models, dQ/ds vs scaling spread and misfit (MPa):')
for lab, sel in [('t 8-10', lambda x: x[1] <= 10), ('t >= 11', lambda x: x[1] >= 11)]:
    r = [x for x in rows if x[0] not in OVERTURNED and sel(x)]
    s = np.array([x[7] for x in r]); q = np.array([x[6] for x in r])
    print('  %-8s scal range %.1f-%.1f (sd %.2f)  dQ/ds sd %.2f  rms(dQ/ds - scal) %.2f' % (lab, s.min(), s.max(), s.std(), q.std(), np.sqrt(np.mean((q-s)**2))))

# --- figure, in the style of Figs 3 and 5 (plot_forces.varH.py closure / fig5), every early output
# (t <= 10) coloured by slab tip depth, marker = plate configuration, overturned hollow; the t >= 11
# outputs of the manuscript figures pale grey behind.  a) closure DP + sigma_slab vs B_slab (Fig. 3),
# b) DP vs B_slab (Fig. 5a), c) dQ/ds vs the scaling (Fig. 5b).
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from cmcrameri import cm as cmc
mpl.rcParams.update({'font.family': 'Myriad Pro', 'font.size': 7, 'axes.labelsize': 9,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^'}
bc_label = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
cmap = cmc.navia_r
early_tips = [x[2] for x in rows if x[1] <= 10]
normT = mpl.colors.Normalize(vmin=50*np.floor(np.nanmin(early_tips)/50), vmax=700)
fig, (axC, axP, axQ) = plt.subplots(1, 3, figsize=(7.0, 2.6))
panels = [(axC, lambda x: x[4], lambda x: x[8] + x[4], (-2.5, 40, -2.5, 40),
           r'$B_{slab}$   [MPa]', r'$\Delta P + \sigma_{slab}$   [MPa]', 'a'),
          (axP, lambda x: x[4], lambda x: x[3], (-3, 40, -3, 40),
           r'$B_{slab}$   [MPa]', r'$\Delta P$   [MPa]', 'b'),
          (axQ, lambda x: x[7], lambda x: x[6], (-5, 40, -10, 30),
           r'$\eta H^3 K v_c / (3 L_K^3)$   [MPa]', r'$-dQ_{s}/ds$   [MPa]', 'c')]
for ax, fx, fy, (x0, x1, y0, y1), xlab, ylab, lab in panels:
    n_off = 0
    for (bc, eta, H) in sorted({x[0] for x in rows}):
        rk = [x for x in rows if x[0] == (bc, eta, H)]
        ot, m = (bc, eta, H) in OVERTURNED, bc_marker[bc]
        late = [x for x in rk if x[1] >= 11]; early = [x for x in rk if x[1] <= 10]
        xl, yl = [fx(x) for x in late], [fy(x) for x in late]
        xe, ye = np.array([fx(x) for x in early]), np.array([fy(x) for x in early])
        te = np.array([x[2] for x in early])
        if ot:
            ax.scatter(xl, yl, s=5, facecolors='none', edgecolor='0.8', marker=m, linewidth=0.3, zorder=1)
            ax.scatter(xe, ye, s=10, facecolors='none', edgecolor=cmap(normT(te)), marker=m, linewidth=0.5,
                       alpha=0.8, zorder=2)
        else:
            ax.scatter(xl, yl, s=4, color='0.8', marker=m, linewidth=0, zorder=1)
            sc = ax.scatter(xe, ye, s=10, c=te, cmap=cmap, norm=normT, edgecolor='black', marker=m,
                            linewidth=0.25, zorder=3)
            n_off += np.sum((xe < x0) | (xe > x1) | (ye < y0) | (ye > y1))
    print('%s) early normal outputs off axes: %d' % (lab, n_off))
    if ax is axP:
        ax.plot([x0, x1], [x0, x1], color='black', linewidth=1, zorder=2)
    else:
        ax.plot([x0, x1], [x0, x1], color='gray', linestyle='--', linewidth=1, zorder=1)
    ax.set_xticks([0, 10, 20, 30, 40]); ax.set_yticks([-10, 0, 10, 20, 30] if ax is axQ else [0, 10, 20, 30, 40])
    ax.axhline(0, color='lightgray', linewidth=1, zorder=0); ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
            bbox=dict(facecolor='white', edgecolor='none', pad=1))
plt.tight_layout()
# legend: one boxed row above the panels -- tip-depth colorbar, then plate configuration, overturned, t >= 11
handles = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                  linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
          [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray', linestyle='', markersize=4,
                  label='overturned'),
           Line2D([], [], marker='o', color='0.8', linestyle='', markersize=3, label='t ≥ 11 (Figs 3, 5)')]
leg = fig.legend(handles=handles, ncol=len(handles), loc='lower left', bbox_to_anchor=(0.47, 0.93), fontsize=6.5,
                 frameon=True, fancybox=False, edgecolor='black', facecolor='white', framealpha=1, borderpad=0.5,
                 handletextpad=0.1, handlelength=1.2, columnspacing=1.0)
leg.get_frame().set_linewidth(0.5)
cax = fig.add_axes([0.13, 0.985, 0.22, 0.025])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=normT, cmap=cmap), cax=cax, orientation='horizontal', extend='max')
cb.ax.tick_params(labelsize=6, length=2.5, pad=1)
cb.set_label('slab tip depth  [km]', fontsize=7, labelpad=2)
out = 'plots/DP-comparisons/compilations/manuscript_figs/early_outputs_supp.z%s' % z
plt.savefig(out + '.png', dpi=300, bbox_inches='tight'); plt.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
