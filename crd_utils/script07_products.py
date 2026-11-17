"""Compact draw summaries, aperture totals and final Script 07 science figures."""
from __future__ import annotations

import os
import warnings

from astropy.table import Table
import numpy as np

from . import global_disk, populations, templates
from .bl_population_products import write_npz
from .global_disk_io import save_json
from .script07_fits import BL_STATS, bl_population

RH3_COLUMNS = ('VA_GLOBAL_KMS', 'VB_GLOBAL_KMS', 'SIGMA_A_KMS', 'SIGMA_B_KMS',
               'FA_RH3', 'FA_RH3_ACHIEVED', 'CHI2_RH3')
DISK = {label: rf'$\mathrm{{Disk}}_{{{label}}}$' for label in 'AB'}

BL_COLUMNS = ('FA_BL', 'FA_BL_ACHIEVED', 'CHI2_BL', 'GAS_V_KMS', 'GAS_SIGMA_KMS',
              'A_STELLAR_MASS_FRACTION') + tuple(
    f'{d}_{k}' for d in 'AB' for k in BL_STATS) + tuple(f'DELTA_{k}' for k in BL_STATS)


def _contrasts(row):
    row.update({f'DELTA_{k}': row[f'A_{k}']-row[f'B_{k}'] for k in BL_STATS})


def aperture_statistics(weights, setup):
    """Aggregate flux-restored template coefficients over the common footprint."""
    usable = np.isfinite(weights).all(axis=(1, 2))
    shape = (2, len(setup['population_metadata']['ages']), len(setup['population_metadata']['metals']))
    if not np.any(usable):
        return {}, np.full(shape, np.nan), 0
    summed = weights[usable].sum(axis=0)
    disk_light = summed.sum(axis=1)
    if np.any(disk_light <= 0):
        raise ValueError('A usable aperture must contain positive stellar weight in each disk')
    stats, mass = {}, []
    for disk, coefficients, light in zip('AB', summed, disk_light):
        values, grids = populations.summarize_ssp_weights(
            coefficients, setup['stellar_normalization'], setup['population_metadata'])
        stats.update({f'{disk}_{k}': values[k] for k in BL_STATS})
        stats[f'{disk}_LIGHT_FRACTION'] = light/disk_light.sum()
        mass.append(grids['total_surviving_mass'])
    stats['A_MASS_FRACTION'] = mass[0]/sum(mass)
    _contrasts(stats)
    return stats, summed.reshape(shape), int(usable.sum())


def measurements(nominal, disk_fit=None, rh3=None, bl=None):
    """Return maps of numerical samples and absolute-weight BL SSP grids."""
    nbin = len(nominal.problem.bin_ids)
    result = {name: np.full(nbin, np.nan) for name in (*RH3_COLUMNS, *BL_COLUMNS, 'DELTA_V_KMS')}
    setup = nominal.setup['BL']
    weights = np.full((nbin, 2, setup['n_stellar']), np.nan)
    fit = disk_fit if disk_fit is not None else nominal.fit
    if 'theta' in fit:
        ids = np.asarray(fit['bin_ids'], dtype=int)
        result['VA_GLOBAL_KMS'][ids] = fit['va']
        result['VB_GLOBAL_KMS'][ids] = fit['vb']
        result['DELTA_V_KMS'][ids] = abs(result['VA_GLOBAL_KMS'][ids]-result['VB_GLOBAL_KMS'][ids])
    if disk_fit is None:
        for arm, columns in (('RH3', RH3_COLUMNS), ('BL', BL_COLUMNS)):
            table = nominal.tables[arm]
            ids = nominal.ids[arm]
            for name in columns:
                source = 'CHI2_TOTAL' if name == f'CHI2_{arm}' else name
                if source not in table.colnames:
                    continue
                result[name][ids] = np.ma.asarray(table[source][ids], dtype=float).filled(np.nan)
        for bid in nominal.ids['BL']:
            winning = nominal.spectra['BL']['model_weights'][int(bid)]
            weights[bid] = winning[:2*setup['n_stellar']].reshape(2, -1)*nominal.spectra['BL']['normalization_scale'][bid]
    else:
        for bid, profile in rh3.items():
            fit_row = profile['winner']
            if fit_row is not None:
                for name, value in zip(RH3_COLUMNS[2:], [*fit_row['sigma'], fit_row['fraction'],
                                                        fit_row['achieved'], fit_row['chi2']]):
                    result[name][bid] = value
        for bid, profile in bl.items():
            if profile['winner'] is not None:
                row, weights[bid] = bl_population(profile['winner'], setup,
                                                  nominal.spectra['BL']['normalization_scale'][bid])
                _contrasts(row)
                for key, value in row.items():
                    result[key][bid] = value
    return result, weights


def combined_rh3_moments(values):
    """Match Script 05's first two moments of the two-disk stellar LOSVD."""
    va, vb = (values[key] for key in ('VA_GLOBAL_KMS', 'VB_GLOBAL_KMS'))
    sa, sb = (values[key] for key in ('SIGMA_A_KMS', 'SIGMA_B_KMS'))
    f = values['FA_RH3_ACHIEVED']
    mean = f*va + (1-f)*vb
    variance = f*sa**2 + (1-f)*sb**2 + f*(1-f)*(va-vb)**2
    return dict(V_COMBINED_KMS=mean, SIGMA_COMBINED_KMS=np.sqrt(variance))


def summarize(values):
    """Sample standard deviation, percentiles and finite sample count by axis."""
    arr = np.asarray(values, dtype=float)
    finite = np.isfinite(arr)
    count = finite.sum(axis=0)
    arr = np.where(finite, arr, np.nan)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Degrees of freedom <= 0', category=RuntimeWarning)
        warnings.filterwarnings('ignore', message='All-NaN slice encountered', category=RuntimeWarning)
        std = np.nanstd(arr, axis=0, ddof=1)
        p16, p84 = np.nanpercentile(arr, [16, 84], axis=0)
    return {'STD': std, 'P16': p16, 'P84': p84, 'N': count}


def _table(path, table):
    tmp = path.with_suffix('.tmp')
    table.write(tmp, format='ascii.ecsv', overwrite=True)
    tmp.replace(path)


def _nominal_positions(nominal):
    p = nominal.problem
    x, y, pa, inc = nominal.fit['theta'][:4]
    radius, cosine = global_disk.disk_coordinates(p.east-x, p.north-y, pa, inc)
    px = (p.member_y, p.member_x)
    counts = np.bincount(p.member_bin, minlength=len(p.bin_ids))
    return tuple(np.bincount(p.member_bin, weights=values[px], minlength=len(p.bin_ids))/counts
                 for values in (radius, radius*cosine))


def export(directory, nominal, draws, attempts):
    output = directory/'products'
    figures = directory/'figures'
    output.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    points, nominal_weights = measurements(nominal)
    meta = nominal.setup['BL']['population_metadata']
    profiles = []
    theta = []
    aperture_samples = []
    ssp_samples = []
    for record in draws:
        with np.load(directory/'draws'/record['data_file'], allow_pickle=False) as data:
            profiles.append({key: data[key] for key in points})
            theta.append(data['theta'])
            ssp_samples.append(data['aperture_ssp'])
        aperture_samples.append(record['aperture'])
    points.update(combined_rh3_moments(points))
    for profile in profiles:
        profile.update(combined_rh3_moments(profile))
    radius, offset = _nominal_positions(nominal)
    bins = Table.read(nominal.paths['bin_table'], format='ascii.ecsv')
    table = Table(dict(BIN_ID=nominal.problem.bin_ids, RADIUS_ARCSEC=radius, SIGNED_MAJOR_ARCSEC=offset,
                       **{k: bins[k] for k in ('BL_SN', 'RH3_SN', 'RA_DEG', 'DEC_DEG') if k in bins.colnames}))
    rh3_table, bl_table = nominal.tables['RH3'], nominal.tables['BL']
    for arm, source in (('RH3', rh3_table), ('BL', bl_table)):
        table[f'{arm}_ELIGIBLE'] = np.isin(table['BIN_ID'], nominal.ids[arm])
        table[f'{arm}_NOMINAL_FLAGS'] = np.ma.asarray(source['FLAGS']).filled('')
        boundary = np.zeros(len(table), dtype=int)
        failed = np.zeros(len(table), dtype=int)
        possible = np.isin(table['BIN_ID'], nominal.ids[arm])
        key = 'FA_RH3' if arm == 'RH3' else 'FA_BL'
        for record, vals in zip(draws, profiles):
            failed += possible & ~np.isfinite(vals[key])
            for bid, flags in record.get('winning_flags', {}).get(arm, {}).items():
                boundary[int(bid)] += bool(any('BOUNDARY' in flag or flag == 'FRACTION_EDGE' for flag in flags))
        table[f'{arm}_FAILED_DRAWS'] = failed
        table[f'{arm}_BOUNDARY_DRAWS'] = boundary
    for name, value in points.items():
        table[name] = value
        for statistic, measures in summarize([prof[name] for prof in profiles]).items():
            table[f'{name}_{statistic}'] = measures
    table.meta.update(target=nominal.provenance['target'], model=nominal.fit['name'],
                      accepted_draws=len(draws), attempted_draws=len(attempts),
                      bl_reference_angstrom=list(nominal.bl_settings['reference']),
                      velocity_frame='Original Script 03 LOS frame', uncertainty='Sample standard deviation')
    _table(output/'stellar_kinematics_and_populations.ecsv', table)

    aperture, nominal_ssp, coverage = aperture_statistics(nominal_weights, nominal.setup['BL'])
    report = Table({'DISK': ['A', 'B', 'A minus B']})
    for metric in (*BL_STATS, 'LIGHT_FRACTION', 'MASS_FRACTION'):
        if metric == 'MASS_FRACTION':
            keys = ['A_MASS_FRACTION', 'B_MASS_FRACTION', 'DELTA_MASS_FRACTION']
        else:
            keys = [f'{prefix}_{metric}' for prefix in ('A', 'B', 'DELTA')]

        def vector(row):
            a, b, delta = (row.get(k) for k in keys)
            if metric == 'MASS_FRACTION' and a is not None:
                b = 1-a
            if delta is None and a is not None and b is not None:
                delta = a-b
            return np.asarray([a, b, delta], dtype=float)

        report[metric] = vector(aperture)
        for name, value in summarize([vector(row) for row in aperture_samples]).items():
            report[f'{metric}_{name}'] = value
    report.meta.update(table.meta, aperture='Nominal fitted KCWI BL footprint',
                       nominal_bins=coverage, eligible_bins=nominal.ids['BL'].tolist(),
                       mass_to_light_filter=nominal.bl_settings['ml_band'], imf=nominal.bl_settings['imf'])
    _table(output/'aperture_stellar_populations.ecsv', report)
    write_npz(output/'aperture_ssp_weights.npz', nominal=nominal_ssp, draws=np.stack(ssp_samples),
              ages=meta['ages'], metallicities=meta['metals'],
              draw_ids=np.asarray([r['draw'] for r in draws]))
    formed = np.stack([nominal_ssp, *ssp_samples])/nominal.setup['BL']['stellar_normalization'].reshape(
        1, 1, len(meta['ages']), len(meta['metals']))
    by_age = formed.sum(axis=-1)
    sfh = by_age/by_age.sum(axis=-1, keepdims=True)
    sfh_table = Table(dict(AGE_GYR=meta['ages']))
    for i, disk in enumerate('AB'):
        sfh_table[f'{disk}_FORMED_MASS_FRACTION'] = sfh[0, i]
        for statistic, values in summarize(sfh[1:, i]).items():
            sfh_table[f'{disk}_{statistic}'] = values
    _table(output/'aperture_star_formation_history.ecsv', sfh_table)
    export_kinematics(output, figures, nominal, theta, offset, profiles, table)
    export_map_figures(figures, nominal, table)
    export_radial_profiles(output, figures, nominal, radius, nominal_weights, draws, directory)
    export_spectrum(output, figures, nominal)
    save_json(output/'ensemble_summary.json', dict(accepted=len(draws), attempted=len(attempts),
              rejected=len(attempts)-len(draws), accepted_draw_ids=[r['draw'] for r in draws],
              stalled=sum('KINEMATIC_RESTART_STALLED' in r.get('disk_flags', []) for r in attempts),
              capped=sum('KINEMATIC_RESTART_LIMIT_REACHED' in r.get('disk_flags', []) for r in attempts)))
    export_index(directory, nominal)


def export_index(directory, nominal):
    references = {
        'Master spatial bins and S/N': nominal.paths['bin_table'],
        'Nominal RH3 residual-noise calibration': nominal.paths['spectra'].parent/'RH3_scalar_noise_calibration.npz',
        'Nominal coherent fit': nominal.paths['fit'],
        'Nominal RH3 stellar PSF QA': nominal.paths['run4']/'figures/stellar_PSF_QA.pdf',
        'Nominal BL stellar PSF QA': nominal.paths['run6']/'figures/BL_stellar_PSF_QA.png',
        'Nominal BL line review': nominal.paths['bl_review'],
        'Nominal RH3 kinematics': nominal.paths['nominal_rh3'],
        'Nominal BL populations': nominal.paths['nominal_bl'],
    }
    lines = ['# Script 07 science package', '', '## Publication figures', '']
    for path in sorted((directory/'figures').glob('*.pdf')):
        lines.append(f'- [{path.stem.replace("_", " ")}](figures/{path.name})')
    lines.extend(['', '## Machine-readable results', ''])
    for path in sorted((directory/'products').iterdir()):
        if path.is_file() and path.suffix in ('.ecsv', '.npz', '.json'):
            lines.append(f'- [{path.name}](products/{path.name})')
    lines.extend(['', '## Nominal analysis', '', 'The original files supply the nominal calibration, masks and full fit diagnostics.', ''])
    for title, path in references.items():
        if path.is_file():
            try:
                target = os.path.relpath(path, directory).replace(os.sep, '/')
            except ValueError:
                target = path.resolve().as_uri()
            lines.append(f'- [{title}](<{target}>)')
    path = directory/'README.md'
    tmp = path.with_suffix('.tmp')
    tmp.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    tmp.replace(path)


def export_kinematics(output, figures, nominal, theta, offset, profiles, table):
    import matplotlib.pyplot as plt

    original = np.asarray(nominal.fit['theta'])
    draws = np.asarray(theta)
    nodes = nominal.problem.nodes
    parameters = ['EAST0_ARCSEC', 'NORTH0_ARCSEC', 'PA_DEG', 'INCL_DEG', 'VSYS_KMS']
    geo = Table(dict(PARAMETER=parameters, VALUE=original[:5], **summarize(draws[:, :5])))
    _table(output/'disk_geometry.ecsv', geo)
    n = len(nodes)
    node_table = Table({'RADIUS_ARCSEC': nodes})
    for k, disk in enumerate('AB'):
        cut = slice(5+k*n, 5+(k+1)*n)
        node_table[f'SPEED_{disk}_KMS'] = original[cut]
        for key, value in summarize(draws[:, cut]).items():
            node_table[f'SPEED_{disk}_{key}'] = value
    _table(output/'rotation_nodes.ecsv', node_table)

    ids = np.asarray(nominal.fit['bin_ids'], dtype=int)
    offsets = offset[ids]
    lo, hi = float(np.min(offsets)), float(np.max(offsets))
    x = np.linspace(lo, hi, 501)
    fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
    ax.axvline(0, color='0.65', lw=.7)
    ax.axhline(0, color='0.65', lw=.7)
    rows = {'SIGNED_MAJOR_ARCSEC': x}
    for i, (disk, color) in enumerate(zip('AB', ('tab:blue', 'tab:red'))):
        cut = slice(5+i*n, 5+(i+1)*n)
        sense = nominal.problem.sign_a * (-1)**i

        def curve(state):
            speed = np.interp(abs(x), np.r_[0., nodes], np.r_[0., state[cut]])
            return sense*np.sign(x)*np.sin(np.deg2rad(state[3]))*speed

        projected = curve(original)
        spread = summarize([curve(state) for state in draws])
        rows[f'V{disk}_KMS'] = projected
        rows.update({f'V{disk}_{k}': v for k, v in spread.items()})
        ax.fill_between(x, spread['P16'], spread['P84'], color=color, alpha=.2)
        ax.plot(x, projected, color=color, label=f'{DISK[disk]} coherent curve')
        k = f'V{disk}_GLOBAL_KMS'
        sample = np.stack([p[k] for p in profiles])-draws[:, 4, None]
        errors = summarize(sample)['STD']
        vel = np.asarray(table[k])-original[4]
        use = ids[np.isfinite(errors[ids])]
        ax.errorbar(offset[use], vel[use], yerr=errors[use], color='0.35',
                    alpha=.35, lw=.75, fmt='none', zorder=2)
        ax.scatter(offset[ids], vel[ids], color=color, s=10, alpha=.65, zorder=3)
    if lo < hi:
        ax.set_xlim(lo, hi)
    ax.set(xlabel='Signed major-axis offset (arcsec)',
           ylabel='Line-of-sight velocity relative to systemic (km/s)')
    ax.legend(title='Bands show the 16th–84th percentiles', fontsize=8)
    _save(figures/'major_axis_rotation_curves.pdf', fig)
    _table(output/'major_axis_rotation_curves.ecsv', Table(rows))


def _save(path, fig):
    import matplotlib.pyplot as plt

    tmp = path.with_suffix('.tmp')
    try:
        fig.savefig(tmp, format='pdf')
        tmp.replace(path)
    finally:
        plt.close(fig)


def export_map_figures(directory, nominal, table):
    import matplotlib.pyplot as plt

    disk_a, disk_b = DISK['A'], DISK['B']
    difference = f'{disk_a} − {disk_b}'
    # RH3 shares a color scale across both disks and the combined panel.
    layouts = {
        'RH3': (
            ('km/s', (('VA_GLOBAL_KMS', f'{disk_a} stellar velocity'),
                      ('VB_GLOBAL_KMS', f'{disk_b} stellar velocity'),
                      ('V_COMBINED_KMS', 'Combined stellar velocity'))),
            ('km/s', (('SIGMA_A_KMS', f'{disk_a} velocity dispersion'),
                      ('SIGMA_B_KMS', f'{disk_b} velocity dispersion'),
                      ('SIGMA_COMBINED_KMS', 'Combined velocity dispersion'))),
        ),
        'BL': (
            ('Gyr', (('A_LIGHT_AGE_GYR', f'{disk_a} light-weighted age'),
                     ('B_LIGHT_AGE_GYR', f'{disk_b} light-weighted age'),
                     ('DELTA_LIGHT_AGE_GYR', f'Age difference\n({difference})'))),
            ('dex', (('A_LIGHT_METAL_DEX', f'{disk_a} light-weighted [M/H]'),
                     ('B_LIGHT_METAL_DEX', f'{disk_b} light-weighted [M/H]'),
                     ('DELTA_LIGHT_METAL_DEX', f'[M/H] difference\n({difference})'))),
            (r'$M_\odot/L_\odot$', (('A_MASS_TO_LIGHT', f'{disk_a} mass-to-light ratio'),
                                      ('B_MASS_TO_LIGHT', f'{disk_b} mass-to-light ratio'),
                                      ('DELTA_MASS_TO_LIGHT', f'Mass-to-light difference\n({difference})'))),
            ('Fraction', (('FA_BL', f'{disk_a} BL light fraction'),
                          ('FB_BL', f'{disk_b} BL light fraction'),
                          ('DELTA_F_BL', f'Light-fraction difference\n({difference})'))),
        ),
    }
    bl_science_colors = ('cool', 'winter', 'pink', 'bone')
    f_a = np.asarray(table['FA_BL'])
    f_std = np.asarray(table['FA_BL_STD'])
    light = {'FB_BL': 1-f_a, 'FB_BL_STD': f_std, 'DELTA_F_BL': 2*f_a-1, 'DELTA_F_BL_STD': 2*f_std}
    for arm, rows in layouts.items():
        east, north, ids = nominal.spatial[arm]
        valid = ids >= 0
        for suffix in ('', '_STD'):
            fig, axes = plt.subplots(len(rows), 3, figsize=(12, 3.6*len(rows)), layout='constrained')
            for row_index, (row, (unit, panels)) in enumerate(zip(axes, rows)):
                arrays = [np.asarray(table[k+suffix], dtype=float) if k+suffix not in light else light[k+suffix]
                          for k, _ in panels]
                shared = np.concatenate(arrays if arm == 'RH3' else arrays[:2])
                finite = shared[np.isfinite(shared)]
                low, high = (float(finite.min()), float(finite.max())) if finite.size else (0., 1.)
                # Use blue for negative velocities and red for positive ones.
                if panels[0][0] == 'VA_GLOBAL_KMS' and not suffix:
                    limit = max(1e-12, abs(low), abs(high))
                    low, high = -limit, limit
                for i, (ax, (_, title), values) in enumerate(zip(row, panels, arrays)):
                    img = np.full(ids.shape, np.nan)
                    img[valid] = values[ids[valid]]
                    cmap = 'viridis'
                    limits = {}
                    if arm == 'RH3':
                        cmap = 'viridis' if suffix else 'RdBu_r'
                        limits = dict(vmin=low, vmax=max(high, low+1e-10))
                    elif i == 2:
                        finite = values[np.isfinite(values)]
                        limit = max(1e-12, np.max(abs(finite))) if finite.size else 1.
                        cmap = 'cividis'
                        limits = dict(vmin=0 if suffix else -limit, vmax=limit)
                    else:
                        cmap = 'viridis' if suffix else bl_science_colors[row_index]
                        limits = dict(vmin=low, vmax=max(high, low+1e-10))
                    color = ax.pcolormesh(east, north, img, shading='nearest', cmap=cmap,
                                          rasterized=True, **limits)
                    fig.colorbar(color, ax=ax, shrink=.83, label=unit)
                    ax.set(title=title, aspect='equal', xlabel='East offset (arcsec)', ylabel='North offset (arcsec)')
                    ax.set_xlim(np.nanmax(east), np.nanmin(east))
            if suffix:
                fig.suptitle('Monte Carlo sample standard deviations')
            _save(directory/f'{arm}_{"uncertainty" if suffix else "science"}_maps.pdf', fig)

def export_radial_profiles(output, figures, nominal, radius, weights, draws, directory):
    import matplotlib.pyplot as plt

    nodes = nominal.problem.nodes
    edges = np.r_[0., (nodes[:-1]+nodes[1:])/2,
                  max(nodes[-1]+(nodes[-1]-nodes[-2])/2 if len(nodes) > 1 else 2*nodes[-1],
                      np.nextafter(radius[nominal.ids['BL']].max(), np.inf))]
    selectors = [(radius >= lo) & (radius < hi) for lo, hi in zip(edges[:-1], edges[1:])]
    keys = ('A_LIGHT_AGE_GYR', 'B_LIGHT_AGE_GYR', 'A_LIGHT_METAL_DEX', 'B_LIGHT_METAL_DEX')
    setup = nominal.setup['BL']

    def stats(matrix):
        rows = [aperture_statistics(matrix[selection], setup)[0] for selection in selectors]
        return np.asarray([[row.get(k, np.nan) for k in keys] for row in rows])

    points = stats(weights)
    samples = []
    for record in draws:
        with np.load(directory/'draws'/record['data_file'], allow_pickle=False) as stored:
            samples.append(stats(stored['bl_weights']))
    spread = summarize(samples)
    table = Table(dict(RADIUS_ARCSEC=nodes, INNER_ARCSEC=edges[:-1], OUTER_ARCSEC=edges[1:]))
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True, layout='constrained')
    for i, name in enumerate(keys):
        table[name] = points[:, i]
        for suffix, value in spread.items():
            table[f'{name}_{suffix}'] = value[:, i]
        ax = axes[i//2]
        color = ['tab:blue', 'tab:red'][i % 2]
        errors = spread['STD'][:, i]
        use = np.isfinite(points[:, i]) & np.isfinite(errors)
        ax.fill_between(nodes, spread['P16'][:, i], spread['P84'][:, i], color=color, alpha=.2)
        ax.errorbar(nodes[use], points[use, i], yerr=errors[use], color='0.35',
                    alpha=.45, lw=.75, fmt='none', zorder=2)
        ax.plot(nodes, points[:, i], 'o-', color=color, label=DISK[name[0]], zorder=3)
    axes[0].set(ylabel='Light-weighted age (Gyr)')
    axes[1].set(xlabel='Nominal deprojected radius (arcsec)', ylabel='Light-weighted [M/H] (dex)')
    axes[0].legend(title='Bands show the 16th–84th percentiles', fontsize=9)
    _table(output/'radial_stellar_populations.ecsv', table)
    _save(figures/'radial_stellar_populations.pdf', fig)


def export_spectrum(output, figures, nominal):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, MultipleLocator

    ids = np.intersect1d(nominal.ids['RH3'], nominal.ids['BL'])
    saved = {'bin_ids': ids}
    # Reference wavelengths in vacuum. Script 03 fits the full unmasked
    # interval, so these markers do not select or separately weight pixels.
    cat = ((8500.4, 'Ca II 8500'), (8544.4, 'Ca II 8544'), (8664.5, 'Ca II 8665'))
    lines = {
        'BL': ((3934.8, 'Ca K'), (4102.9, 'Hδ'), (4341.7, 'Hγ'),
               (4862.7, 'Hβ'), (5176.7, 'Mg b'), (5271.0, 'Fe')),
        'RH3': ((5891.6, 'Na D₂'), (5897.6, 'Na D₁'), (6564.6, 'Hα'),
                (8185.5, 'Na I 8186'), (8197.1, 'Na I 8197'), *cat, (8809.2, 'Mg I 8809')),
    }
    arms = {}
    for arm in ('BL', 'RH3'):
        spec = nominal.spectra[arm]
        medium = nominal.bl_settings['template_medium'] if arm == 'BL' else load_rh3_medium(nominal)
        wave = templates.convert_wavelength_medium(spec['wavelength'], medium, 'vacuum')
        good = spec['good'][ids]
        scale = spec['normalization_scale'][ids, None]
        covered = good.any(axis=0)
        channels = dict(observed=spec['galaxy'], model=spec['model'], **spec['components'])
        fluxes = {name: np.sum(np.where(good, arr[ids]*scale, 0), axis=0) for name, arr in channels.items()}
        for flux in fluxes.values():
            flux[~covered] = np.nan
        residuals = fluxes['observed']-fluxes['model']
        arms[arm] = (wave, covered, fluxes, residuals)
        saved.update({f'{arm}_{name}': flux for name, flux in fluxes.items()})
        saved[f'{arm}_wavelength_vacuum_angstrom'] = wave
        saved[f'{arm}_contributors'] = good.sum(axis=0)
        saved[f'{arm}_residual'] = residuals
        if arm == 'RH3':
            saved['RH3_additive'] = fluxes['model']-fluxes['A']-fluxes['B']

    colors = dict(observed='0.3', model='black', A='tab:blue', B='tab:red', gas='#da790a')
    labels = dict(observed='Observed', model='Total fit', A=DISK['A'], B=DISK['B'], gas='Ionized gas')
    fig = plt.figure(figsize=(13, 8.4), layout='constrained')
    grid = fig.add_gridspec(3, 1, height_ratios=(1, 1, .14))
    top = grid[0].subgridspec(1, 2)
    spans = [np.ptp(arm[0]) for arm in arms.values()]
    # Width ratios preserve the wavelength scale across the break.
    wide = grid[1].subgridspec(1, 2, width_ratios=spans, wspace=.04)
    blue = fig.add_subplot(top[0])
    red = fig.add_subplot(top[1], sharey=blue)
    long_blue = fig.add_subplot(wide[0], sharey=blue)
    long_red = fig.add_subplot(wide[1], sharey=blue)
    step = np.diff(MaxNLocator(nbins=4).tick_values(0, max(spans)))[0]
    for arm, ax, combined in (('BL', blue, long_blue), ('RH3', red, long_red)):
        wave, covered, fluxes, residuals = arms[arm]
        for panel in (ax, combined):
            panel.axhline(0, color='0.58', lw=.6, ls='--', zorder=0)
            for name, values in fluxes.items():
                panel.plot(wave, values, label=labels[name], lw=1.0 if name == 'gas' else .8,
                           color=colors[name], zorder=4 if name == 'gas' else 2)
            panel.plot(wave, residuals, '.', ls='none', ms=2.6, color='tab:green',
                       alpha=.5, label='Residuals', zorder=3)
            panel.set_xlim(wave[0], wave[-1])
            for center, label in lines[arm]:
                if wave[0] <= center <= wave[-1]:
                    masked = not covered[np.argmin(abs(wave-center))]
                    panel.axvline(center, color='0.65', lw=.65, ls=':' if masked else '-', zorder=0)
                    align = 'left' if label in ('Na D₁', 'Na I 8197') else 'right'
                    panel.text(center, .97, label, rotation=90, transform=panel.get_xaxis_transform(),
                               va='top', ha=align, fontsize=8, color='0.45' if masked else '0.2')
        if arm == 'RH3' and not any(wave[0] <= center <= wave[-1] for center, _ in cat):
            ax.text(.025, .26, 'Ca II triplet (8500, 8544, 8665 Å) outside this fit',
                    transform=ax.transAxes, fontsize=8, color='0.35', va='bottom')
        ax.set(title=f'{arm} arm')
        combined.xaxis.set_major_locator(MultipleLocator(step))

    long_blue.spines['right'].set_visible(False)
    long_red.spines['left'].set_visible(False)
    long_red.tick_params(left=False, labelleft=False)
    # A diagonal stroke at both edges signals the omitted wavelength span.
    for ax, side in ((long_blue, 1), (long_red, 0)):
        ax.plot([side]*2, [0, 1], marker=[(-1, -.6), (1, .6)], markersize=8,
                linestyle='none', color='0.2', transform=ax.transAxes, clip_on=False, zorder=5)
    for ax in (blue, long_blue):
        ax.set_ylabel('Summed aperture flux density')
    red.tick_params(labelleft=False)
    footer = fig.add_subplot(grid[2])
    footer.set_axis_off()
    footer.text(.5, .97, 'Rest-frame vacuum wavelength (Å)', ha='center', va='top')
    handles, series = blue.get_legend_handles_labels()
    footer.legend(handles, series, loc='lower center', ncol=len(series), fontsize=9, frameon=False)
    fig.suptitle('Integrated spectrum in the matched KCWI aperture')
    write_npz(output/'integrated_BL_RH3_spectrum.npz', **saved)
    _save(figures/'integrated_BL_RH3_spectrum.pdf', fig)


def load_rh3_medium(nominal):
    from .config import load_config

    cfg = load_config(nominal.paths['script3_config'], validate=False)
    return cfg.TEMPLATE_WAVELENGTH_MEDIUM
