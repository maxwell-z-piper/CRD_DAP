"""Script 05 diagnostics drawn from saved fraction profiles and bin membership."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .global_disk_plotting import COLORS, _orient_sky_map
from .rh3_refit import best_fraction


def save_pdf(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.pdf.tmp')
    try:
        fig.savefig(temp, format='pdf', dpi=180)
        temp.replace(path)
    finally:
        plt.close(fig)
        temp.unlink(missing_ok=True)


def bin_diagnostic(path, profile, spectrum, wavelength):
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), layout='constrained',
                             gridspec_kw=dict(height_ratios=[1, 1.5, .65]))
    fraction = profile['fraction']
    chi2 = profile['chi2']
    best = best_fraction(profile)
    finite = np.isfinite(chi2)
    ax = axes[0]
    if best is not None:
        ax.plot(fraction, chi2-chi2[best], '.-', ms=3, lw=.7, color='0.3')
        ax.scatter(fraction[best], 0, color='tab:red', zorder=3, label='Best available')
    if not finite.all():
        ax.scatter(fraction[~finite], np.full((~finite).sum(), .07), marker='x',
                   transform=ax.get_xaxis_transform(), color='tab:red', label='Missing fit')
    ax.set(xlabel='RH3 stellar light fraction of Disk A', ylabel=r'$\Delta\chi^2$',
           xlim=(0, 1))
    if best is not None or not finite.all():
        ax.legend(loc='best', fontsize=8)
    va, vb = profile['velocities']
    fig.suptitle(f'RH3 bin {int(profile["bin_id"])} | fixed LOS velocities {va:.2f}, {vb:.2f} km/s')
    good = spectrum['good']
    ax = axes[1]
    ax.plot(wavelength, np.where(good, spectrum['galaxy'], np.nan), color='0.25',
            lw=.8, label='Fitted pixels')
    ax.set(ylabel='Normalized flux')
    if best is not None:
        ax.plot(wavelength, profile['model_total'][best], color='tab:green', lw=.8, label='Total fit')
        for label, color in zip('AB', COLORS):
            ax.plot(wavelength, profile[f'model_{label}'][best], color=color, lw=.8,
                    label=f'Disk {label}')
        sa, sb = profile['sigma'][best]
        winner = profile['attempts'][best][profile['winning_attempt'][best]]
        flags = winner['flags'][:]
        if best in (0, len(fraction)-1):
            flags.append('FRACTION_EDGE')
        if not finite.all():
            flags.append('INCOMPLETE_PROFILE')
        ax.set_title(f'fA = {fraction[best]:.2f} | sigma = {sa:.1f}, {sb:.1f} km/s | '
                     f'total chi-square = {chi2[best]:.2f}\n'
                     + (', '.join(flags) if flags else 'No best-fit flags'), fontsize=9)
        residual = (spectrum['galaxy'][good] - profile['model_total'][best, good]) / spectrum['noise'][good]
        axes[2].scatter(wavelength[good], residual, color='0.25', s=4, rasterized=True)
    else:
        ax.set_title('No scorable fraction')
    ax.legend(ncol=2, fontsize=8)
    axes[2].axhline(0, color='0.6', lw=.8)
    axes[2].set(xlabel='Rest wavelength (Angstrom)', ylabel='Residual / error')
    save_pdf(fig, path)


def kinematic_maps(path, table, inputs):
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout='constrained')
    east, north, bin_map = (inputs[name] for name in ('east', 'north', 'bin_map'))
    vsys = inputs['fit']['theta'][4]
    covered = bin_map >= 0
    velocity_names = ('VA_GLOBAL_KMS', 'VB_GLOBAL_KMS', 'V_COMBINED_KMS')
    sigma_names = ('SIGMA_A_KMS', 'SIGMA_B_KMS', 'SIGMA_COMBINED_KMS')
    usable = np.asarray(table['USABLE'], dtype=bool)
    velocity = [np.where(usable, table[name]-vsys, np.nan) for name in velocity_names]
    sigma = [np.where(usable, table[name], np.nan) for name in sigma_names]
    vlim = max(1., np.max(np.abs(np.concatenate([v[np.isfinite(v)] for v in velocity])), initial=0))
    slim = max(1., np.max(np.concatenate([s[np.isfinite(s)] for s in sigma]), initial=0))
    masked_bin = bin_map[covered]
    for row, values in enumerate((velocity, sigma)):
        for column, (ax, value) in enumerate(zip(axes[row], values)):
            # Pale gray is outside the selected domain. Dark gray is a
            # selected bin without a usable RH3 solution.
            context = np.full(bin_map.shape, np.nan)
            context[covered] = np.where(np.asarray(table['SELECTED'])[masked_bin], .55, .17)
            ax.pcolormesh(east, north, context, cmap='Greys', vmin=0, vmax=1,
                          shading='nearest', rasterized=True)
            image = np.full(bin_map.shape, np.nan)
            image[covered] = value[masked_bin]
            artist = ax.pcolormesh(east, north, image, shading='nearest', rasterized=True,
                                   cmap='RdBu_r' if row == 0 else 'viridis',
                                   vmin=-vlim if row == 0 else 0, vmax=vlim if row == 0 else slim)
            title = ('Disk A', 'Disk B', 'Combined light')[column]
            ax.set(title=f'{title} | {"LOS velocity" if row == 0 else "dispersion"}',
                   xlabel='East offset (arcsec)', ylabel='North offset (arcsec)', aspect='equal')
            _orient_sky_map(ax)
            fig.colorbar(artist, ax=ax, label='km/s relative to systemic' if row == 0 else 'km/s')
    fig.suptitle(f'{inputs["fit"]["name"]} | best RH3 fraction in each bin\n'
                 'Combined light shows reconstructed LOSVD moments. '
                 'Pale gray = unmodeled, dark gray = no usable fit.')
    save_pdf(fig, path)


def major_axis_profile(path, inputs):
    fit = inputs['fit']
    nodes = np.asarray(fit['nodes'])
    x0, y0, pa, inc, vsys = fit['theta'][:5]
    theta = np.deg2rad(pa)
    major = ((inputs['east']-x0)*np.sin(theta) + (inputs['north']-y0)*np.cos(theta))
    member = np.isin(inputs['bin_map'], inputs['ids'])
    positions = major[member]
    members = inputs['bin_map'][member]
    count = np.bincount(members)
    centroids = np.bincount(members, weights=positions)[inputs['ids']] / count[inputs['ids']]
    extent = max(nodes[-1], np.max(np.abs(positions)))
    offset = np.linspace(-extent, extent, 501)
    fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
    ax.axvspan(positions.min(), positions.max(), color='0.93',
               label='Projected member-pixel span')
    ax.axhline(0, color='0.6', lw=.7)
    ax.axvline(0, color='0.6', lw=.7)
    projection = np.sin(np.deg2rad(inc))
    for disk, sense in enumerate((fit['sign_a'], -fit['sign_a'])):
        begin = 5 + disk*len(nodes)
        speed = fit['theta'][begin:begin+len(nodes)]
        los = sense*np.sign(offset)*projection*np.interp(np.abs(offset), np.r_[0, nodes], np.r_[0, speed])
        color = COLORS[disk]
        ax.plot(offset, los, color=color, label=f'Disk {"AB"[disk]} unblurred major-axis model')
        ax.scatter(centroids, inputs['velocities'][:, disk]-vsys, color=color, s=9,
                   alpha=.35, rasterized=True, label=f'Disk {"AB"[disk]} bin predictions')
    ax.set(xlabel='Signed major-axis offset (arcsec)', ylabel='LOS velocity relative to systemic (km/s)',
           title='Selected coherent model across the galaxy\n'
                 'Bin predictions plotted at projected centroids, including off-axis bins')
    ax.legend(fontsize=8, ncol=2)
    save_pdf(fig, path)
