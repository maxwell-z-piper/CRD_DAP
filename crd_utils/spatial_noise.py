"""Empirical scalar spatial-noise correction from KcwiKit stack replays.

Each mock starts with independent reduced-exposure UNCERT draws and is stacked
with the fixed science recipe. The resulting bin sums are measured with the
actual Script-1 masks and Script-2 bin membership. This calibrates the variance
of each aperture sum; it does not estimate spectral lag covariance.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import binning


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _check_prepared_variance(cube, reference_variance):
    """Verify the replay variance, allowing only Stage-1 zero-padding NaNs.

    Script 1 saves positive variance as a one-sigma UNCERT extension. Where
    KcwiKit writes zero variance for uncovered pixels, Script 1 writes NaN
    uncertainty instead. Squaring UNCERT in Script 2 cannot recover those
    zeros, but the pixels are excluded from the hard-good mask.
    """
    prepared = cube.variance
    if reference_variance.shape != prepared.shape:
        raise ValueError(
            f'{cube.arm}: prepared/reference variance shapes differ: '
            f'{prepared.shape} versus {reference_variance.shape}'
        )
    same = np.isclose(reference_variance, prepared,
                      equal_nan=True, rtol=2e-6, atol=1e-8)
    expected_padding = (reference_variance == 0) & np.isnan(prepared) & ~cube.good
    mismatch = ~(same | expected_padding)
    if np.any(mismatch):
        raise ValueError(
            f'{cube.arm}: prepared variance differs from noise-replay reference '
            f'at {int(np.count_nonzero(mismatch))} samples '
            f'({int(np.count_nonzero(mismatch & cube.good))} marked good)'
        )


def calibrate(cube, bin_map, formal, *, manifest_path, n_bins,
              pixel_area_arcsec2, min_member_fraction, rest_range, redshift,
              output_path, min_channels=20):
    """Return corrected spectra and one spatial multiplier per actual bin."""
    path = Path(manifest_path).expanduser().resolve()
    record = json.loads(path.read_text())
    if record.get('schema') != 'crd-noise-stack-ensemble-v1' or not record.get('complete'):
        raise ValueError(f'Incomplete or unsupported noise stack ensemble: {path}')
    if record.get('arm') != cube.arm or not record.get('science_replay_passed'):
        raise ValueError(f'Noise stack arm/replay mismatch: {path}')
    noise_paths = [Path(p) for p in record['noise_icubes']]
    if len(noise_paths) < 4 or len(set(noise_paths)) != len(noise_paths):
        raise ValueError('Need at least four distinct stack noise realizations')
    for letter in 'ivme':
        reference_product = Path(record['reference_products'][letter])
        if _sha256(reference_product) != record['reference_sha256'][letter]:
            raise ValueError(f'Science {letter} reference changed since noise replay')
    reference = Path(record['reference_products']['i'])
    # The prepared cube may add masks, but cannot change science flux, grid, or
    # formal variance; otherwise mock aperture sums would refer to another data set.
    with fits.open(reference, memmap=True) as hdus:
        ref = np.moveaxis(hdus[0].data, 0, -1)
        if ref.shape != cube.flux.shape or not np.allclose(ref, cube.flux, equal_nan=True, rtol=2e-6, atol=1e-8):
            raise ValueError('Prepared science cube differs from noise-replay reference')
    with fits.open(record['reference_products']['v'], memmap=True) as hdus:
        ref_var = np.moveaxis(hdus[0].data, 0, -1)
        _check_prepared_variance(cube, ref_var)
    selector = (cube.wavelength >= float(rest_range[0]) * (1 + redshift)) & (cube.wavelength <= float(rest_range[1]) * (1 + redshift))
    if not np.any(selector):
        raise ValueError('Spatial calibration rest range has no wavelength coverage')
    count = 0
    mean = np.zeros_like(formal.flux)
    m2 = np.zeros_like(formal.flux)
    for i, mock_path in enumerate(noise_paths):
        expected = record.get('noise_sha256', {}).get(str(mock_path))
        if not expected or _sha256(mock_path) != expected:
            raise ValueError(f'Noise stack modified after generation: {mock_path}')
        with fits.open(mock_path, memmap=True) as hdus:
            mock = np.moveaxis(np.asarray(hdus[0].data), 0, -1)
            if mock.shape != cube.flux.shape:
                raise ValueError(f'Noise stack grid mismatch: {mock_path}')
            replay_cube = replace(cube, flux=mock)
            sums = binning.coadd_bin_spectra(replay_cube, bin_map,
                n_bins=n_bins, pixel_area_arcsec2=pixel_area_arcsec2,
                min_member_fraction=min_member_fraction)
            if not np.array_equal(sums.good, formal.good):
                raise ValueError(f'Noise replay changed the good-sample mask: {mock_path}')
            if not np.array_equal(sums.contributing_spaxels, formal.contributing_spaxels):
                raise ValueError(f'Noise replay changed aperture membership: {mock_path}')
            sample = np.where(formal.good, sums.flux, 0.)
        count += 1
        delta = sample - mean
        mean += delta / count
        m2 += delta * (sample - mean)
    empirical_var = m2 / (count - 1)
    formal_var = formal.uncertainty ** 2
    ratio = np.divide(empirical_var, formal_var, out=np.full_like(empirical_var, np.nan),
                      where=formal.good & np.isfinite(formal_var) & (formal_var > 0))
    scales = np.full(n_bins, np.nan)
    used = np.zeros(n_bins, dtype=int)
    for bid in range(n_bins):
        good = formal.good[bid] & selector & np.isfinite(ratio[bid]) & (ratio[bid] > 0)
        used[bid] = np.count_nonzero(good)
        if used[bid] < min_channels:
            raise ValueError(f'Bin {bid}: only {used[bid]} valid spatial-noise calibration channels')
        # A median across channels is robust to isolated bad variance estimates.
        scales[bid] = np.sqrt(np.median(ratio[bid, good]))
    corrected = replace(formal, uncertainty=formal.uncertainty * scales[:, None])
    output_path = Path(output_path)
    np.savez_compressed(output_path, scale=scales, channels_used=used,
        empirical_variance=empirical_var, formal_variance=formal_var,
        wave=cube.wavelength, manifest_sha256=np.asarray(_sha256(path)),
        n_realizations=np.asarray(count))
    return corrected, scales
