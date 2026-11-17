"""Script 06 line-review sheet and population diagnostics."""
from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .bl_population_io import LINES
from .global_disk_plotting import COLORS, _orient_sky_map
from .rh3_refit import best_fraction
from .rh3_refit_plotting import save_pdf


# Guard the legacy checkpoint mapping against later edits to the fit driver.
LEGACY_COMPATIBLE_DRIVER_SHA256 = '04483604e7e30d66d01cb8897da0c0e15eef9972f37b9dc4ff613498d6b4bc3f'


def line_review(path, inputs, experiment, review):
    settings = inputs['contract']['settings']
    redshift = inputs['contract']['redshift']
    eligible = inputs['eligible']
    flux = inputs['flux'][eligible]
    good = experiment['native_good'][eligible]
    rest = inputs['rest_wave']
    counts = good.sum(axis=0)
    summed = np.where(counts > 0, np.where(good, flux, 0).sum(axis=0), np.nan)
    fit_lo, fit_hi = experiment['effective_range']
    in_band = (rest >= fit_lo) & (rest <= fit_hi)
    scale = np.nanmedian(np.abs(summed[in_band]))
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.
    y = summed/scale
    # Preserve the pre-review valid flux under the shaded exclusions to aid inspection.
    original_good = inputs['good'][eligible]
    original_sum = np.where(original_good, flux, 0).sum(axis=0)/scale
    original_sum[~original_good.any(axis=0)] = np.nan

    names = list(LINES)
    displayed = [name for name in names if name in review['lines']]
    optional = [name for name in names if name not in displayed]
    ncols = 4
    rows = (len(names)+ncols-1)//ncols
    fig = plt.figure(figsize=(15, 4.2+1.95*rows), layout='constrained')
    grid = fig.add_gridspec(rows+2, ncols, height_ratios=[2.25, .50] + [1.6]*rows)
    ax = fig.add_subplot(grid[0, :])
    ax.plot(rest, original_sum, color='.70', lw=.7, label='Original usable BL flux')
    ax.plot(rest, y, color='.15', lw=.7, label='Proposed usable BL flux')
    ax.fill_between(rest, 0, 1, where=(counts == 0), transform=ax.get_xaxis_transform(),
                    color='.35', alpha=.15, step='mid', label='No contributing bins')
    ax.axvspan(*settings['reference'], color='tab:blue', alpha=.12, label='Stellar fraction reference')
    ax.set(xlim=(fit_lo-20, fit_hi+20), ylabel='Summed flux density / display scale')
    span = np.nanpercentile(y[in_band], [1, 99])
    ax.set_ylim(span[0]-.1*np.ptp(span), span[1]+.37*np.ptp(span)+1e-8)
    ax.legend(loc='lower center', fontsize=8, ncol=4)
    coverage_ax = fig.add_subplot(grid[1, :], sharex=ax)
    coverage_ax.fill_between(rest, 0, counts, color='.45', step='mid')
    coverage_ax.set(ylabel='Bins', ylim=(0, eligible.sum()*1.05),
                    xlabel=f'Rest wavelength ({settings["template_medium"]}, Å)')
    for edge in (fit_lo, fit_hi):
        ax.axvline(edge, color='.35', lw=.8, ls=':')
    for feature, center in (('Dn4000', 4000), ('Mg b', 5175), ('Fe5270', 5270)):
        ax.axvline(center, lw=.7, color='tab:purple', alpha=.65)
        ax.annotate(feature, (center, .03), xycoords=('data', 'axes fraction'),
                    xytext=(2, 1), textcoords='offset points', fontsize=8, color='tab:purple')
    pairs = {'OII3726': ('OII3729', 'O II doublet'),
             'NI5198': ('NI5200', 'N I doublet')}
    for order, name in enumerate(displayed):
        i = names.index(name)
        center = np.mean(experiment['core'][i])
        if np.isfinite(center):
            ax.axvline(center, lw=.6, color='tab:orange')
            if not any(name == other and primary in displayed
                       for primary, (other, _) in pairs.items()):
                companion, group = pairs.get(name, (None, name))
                label = group if companion in displayed else name
                ax.annotate(label, (center, .81+.06*(order % 2)), xycoords=('data', 'axes fraction'),
                            rotation=45, fontsize=8, color='saddlebrown')
    for j, name in enumerate(displayed+optional):
        i = names.index(name)
        panel = fig.add_subplot(grid[2+j//ncols, j % ncols])
        nominal = experiment['line_wave'][i]
        core = experiment['core'][i]
        allowed = experiment['search'][i]
        center = np.mean(core) if np.isfinite(core).all() else nominal
        left = min(center-18, allowed[0]-4)
        right = max(center+18, allowed[1]+4)
        panel.axvspan(*allowed, color='tab:orange', alpha=.12)
        if np.isfinite(core).all():
            panel.axvspan(*core, color='tab:orange', alpha=.3)
        panel.axvline(center, color='saddlebrown', lw=.65, ls='--')
        panel.plot(rest, original_sum, color='.66', lw=.8)
        panel.plot(rest, y, color='.12', lw=.8)
        masked = (counts < eligible.sum()) & (rest >= left) & (rest <= right)
        panel.fill_between(rest, 0, 1, where=masked, color='.35', alpha=.16, step='mid',
                           transform=panel.get_xaxis_transform())
        selected = (rest >= left) & (rest <= right)
        values = original_sum[selected & np.isfinite(original_sum)]
        if values.size:
            lo, hi = np.percentile(values, [0, 100])
            margin = max(.025, (hi-lo)*.10)
            panel.set_ylim(lo-margin, hi+margin)
        frac = experiment['coverage'][eligible, i]
        n_use = experiment['usable'][:, i].sum()
        n_full = np.count_nonzero(frac >= 1-1e-8)
        status = (f'{n_use} usable / {n_full} full' if name in review['lines'] and settings['gas']['enabled']
                  else 'not proposed')
        panel.set(xlim=(left, right), title=f'{name}  |  {status}')
        panel.title.set_fontsize(9)
        panel.tick_params(labelsize=7)
    target = inputs['contract']['target']
    fig.suptitle(f'{target}  |  Script 06 BL line and mask review  |  z = {redshift:.6f}\n'
                 'Dark orange: nominal core. Pale orange: allowed gas line centers. '
                 'Gray shading: one or more eligible bins excluded at that wavelength.',
                 fontsize=12)
    save_pdf(fig, path)


def bin_diagnostic(path, bid, profile, wave, galaxy, noise, good, settings):
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 8), layout='constrained',
                             gridspec_kw={'height_ratios': [1, 1.5, .65]})
    f = profile['fraction']
    best = best_fraction(profile)
    valid = np.isfinite(profile['chi2'])
    score = profile['chi2_calibrated'] if profile['residual_scale'] > 0 else profile['chi2']
    if best is not None:
        axes[0].plot(f, score-score[best], '.-', color='.25', lw=.7, ms=3)
        axes[0].scatter(f[best], 0, color='tab:red', zorder=3, label='Best available')
    if not valid.all():
        axes[0].scatter(f[~valid], np.full((~valid).sum(), .07), marker='x',
                        transform=axes[0].get_xaxis_transform(), color='tab:red', label='Missing fit')
    axes[0].set(xlim=(0, 1), xlabel='Disk A stellar light fraction (4800–5300 Å convention)'
                if tuple(settings['reference']) == (4800., 5300.) else 'Disk A stellar light fraction',
                ylabel=r'$\Delta\chi^2$')
    axes[0].legend(loc='best', fontsize=8)
    axes[1].plot(wave, np.where(good, galaxy, np.nan), color='.25', lw=.7, label='BL spectrum')
    axes[1].set(ylabel='Normalized flux density')
    if best is not None:
        for key, label, color in zip(('model_total', 'model_A', 'model_B', 'model_gas'),
                                     ('Total', 'Disk A', 'Disk B', 'Gas'),
                                     ('tab:green', *COLORS, 'tab:purple')):
            axes[1].plot(wave, profile[key][best], label=label, color=color, lw=.75)
        residual = (galaxy[good]-profile['model_total'][best, good])/noise[good]
        q = float(profile['residual_scale'])
        if np.isfinite(q) and q > 0:
            residual /= q
        axes[2].scatter(wave[good], residual, s=3, color='.25', rasterized=True)
        axes[1].set_title(f'fA = {f[best]:.2f} | residual noise factor {q:.3g}', fontsize=9)
    else:
        axes[1].set_title('No scorable fraction', fontsize=9)
    axes[1].legend(ncol=5, fontsize=8)
    axes[2].axhline(0, lw=.7, color='.6')
    axes[2].set(xlabel=f'Rest wavelength ({settings["template_medium"]}, Å)',
                ylabel='Residual / error')
    fig.suptitle(f'BL bin {bid}  |  Fixed Script 05 stellar kinematics')
    save_pdf(fig, path)


def population_maps(path, table, inputs, band):
    fig, axes = plt.subplots(4, 3, figsize=(13.5, 14.5), layout='constrained')
    bid = inputs['bin_map']
    covered = bid >= 0
    quantities = (
        (('A_LIGHT_AGE_GYR', 'B_LIGHT_AGE_GYR', 'DELTA_LIGHT_AGE_GYR'), 'Light-weighted age (Gyr)'),
        (('A_LIGHT_METAL_DEX', 'B_LIGHT_METAL_DEX', 'DELTA_LIGHT_METAL_DEX'), 'Light-weighted metallicity (dex)'),
        (('A_MASS_TO_LIGHT', 'B_MASS_TO_LIGHT', 'DELTA_MASS_TO_LIGHT'), f'Mass to light ({band})'),
    )
    rows = []
    for names, label in quantities:
        values = tuple(np.asarray(table[key], dtype=float) for key in names)
        both = np.concatenate(values[:2])
        finite = both[np.isfinite(both)]
        limits = (float(finite.min()), float(finite.max())) if len(finite) else (0., 1.)
        delta = np.abs(values[2][np.isfinite(values[2])])
        limit = float(delta.max()) if delta.size and delta.max() > 0 else 1.
        rows.append((values, (label, label, label), (limits, limits, (-limit, limit))))

    f_a = np.where(table['USABLE'], np.asarray(table['FA_BL'], dtype=float), np.nan)
    rows.append(((100*f_a, 100*(1-f_a), 100*(2*f_a-1)),
                 ('BL stellar light (%)', 'BL stellar light (%)', 'Light difference (percentage points)'),
                 ((0., 100.), (0., 100.), (-100., 100.))))
    context = np.full(bid.shape, np.nan)
    context[covered] = np.where(np.asarray(table['RH3_ELIGIBLE'])[bid[covered]], .52, .18)
    for row_axes, (values, labels, limits) in zip(axes, rows):
        for col, (ax, value, label, (lo, hi)) in enumerate(zip(row_axes, values, labels, limits)):
            image = np.full(bid.shape, np.nan)
            image[covered] = value[bid[covered]]
            ax.pcolormesh(inputs['east'], inputs['north'], context, shading='nearest', cmap='Greys',
                          vmin=0, vmax=1, rasterized=True)
            plotted = ax.pcolormesh(inputs['east'], inputs['north'], image, shading='nearest',
                                   cmap='RdBu_r' if col == 2 else 'viridis', vmin=lo,
                                   vmax=hi if hi > lo else lo+1e-8, rasterized=True)
            fig.colorbar(plotted, ax=ax, label=label)
            ax.set(aspect='equal', xlabel='East offset (arcsec)', ylabel='North offset (arcsec)',
                   title=('Disk A', 'Disk B', 'A minus B')[col])
            _orient_sky_map(ax)
    reference = table.meta['fraction_reference_angstrom']
    fig.suptitle('Script 06 BL stellar populations at fixed RH3 kinematics\n'
                 f'Light fractions use the {reference[0]:g}-{reference[1]:g} Å rest-frame band.\n'
                 'Dark gray: eligible without a usable BL fit. Pale gray: unavailable RH3 anchor.')
    save_pdf(fig, path)
