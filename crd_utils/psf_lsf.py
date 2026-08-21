"""Empirical PSF and instrumental LSF characterization for CRD_DAP.

The production LSF is measured from the KCWI DRP master arc plus its geometry
maps rather than inferred only from a nominal resolving power.  The DRP uses the
master-arc root to generate wavelength, slice, and position maps; Script 1 uses
those same maps to collapse detector pixels into wavelength-calibrated arc
spectra and fit unresolved line widths.

The PSF is treated separately.  An extended galaxy cube cannot, by itself,
provide an unbiased seeing measurement.  CRD_DAP therefore prefers an explicit
configuration value or a point-source reference cube; a recognizable header
seeing value is an allowed fallback and its provenance is recorded.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from .io import discover_arc_sidecars


@dataclass(frozen=True)
class ArcLineMeasurement:
    wavelength: float
    fwhm_angstrom: float
    sigma_angstrom: float
    amplitude: float
    slice_id: int
    position_bin: int
    position_center: float
    reduced_chi2_proxy: float


@dataclass(frozen=True)
class ArcLSFResult:
    """Empirical line-width measurements and the adopted wavelength model."""

    wavelength: np.ndarray
    fwhm_angstrom: np.ndarray
    sigma_angstrom: np.ndarray
    slice_id: np.ndarray
    position_bin: np.ndarray
    position_center: np.ndarray
    reduced_chi2_proxy: np.ndarray
    polynomial_coefficients: np.ndarray
    polynomial_order: int
    # Instrument-good/model-requested wavelength interval, normally WAVGOOD0/1.
    # This is *not* necessarily the interval directly constrained by accepted
    # arc-line measurements.
    wavelength_min: float
    wavelength_max: float
    # Empirical support of the fitted polynomial: min/max accepted line centers
    # after sigma clipping.  Downstream science code must not silently evaluate
    # the LSF outside this interval.
    measurement_wavelength_min: float
    measurement_wavelength_max: float
    spectral_sampling_angstrom: float
    spatial_fractional_rms: float
    measurement_fractional_rms: float
    slice_fractional_rms: float
    position_fractional_rms: float
    n_lines_total: int
    n_lines_used: int

    def evaluate_fwhm(
        self,
        wavelength: np.ndarray | float,
        *,
        allow_extrapolation: bool = False,
    ) -> np.ndarray:
        """Evaluate the wavelength-only LSF model.

        By default values outside the empirically constrained arc-line interval
        are returned as NaN.  This is intentional: a smooth polynomial can look
        plausible far outside the available line measurements, but using that
        extrapolation in stellar-kinematic work would be scientifically unsafe.
        Set ``allow_extrapolation=True`` only for an explicitly labeled QC plot
        or sensitivity test.
        """
        wave = np.asarray(wavelength, dtype=float)
        values = np.polyval(self.polynomial_coefficients, wave)
        if allow_extrapolation:
            return values
        supported = (
            (wave >= self.measurement_wavelength_min)
            & (wave <= self.measurement_wavelength_max)
        )
        return np.where(supported, values, np.nan)

    @property
    def blue_edge_unconstrained_angstrom(self) -> float:
        return float(max(0.0, self.measurement_wavelength_min - self.wavelength_min))

    @property
    def red_edge_unconstrained_angstrom(self) -> float:
        return float(max(0.0, self.wavelength_max - self.measurement_wavelength_max))


@dataclass(frozen=True)
class PSFEstimate:
    fwhm_arcsec: float
    source: str
    detail: str


# -----------------------------------------------------------------------------
# Gaussian width conversion
# -----------------------------------------------------------------------------


def gaussian_fwhm_to_sigma(fwhm: np.ndarray | float) -> np.ndarray:
    """Convert Gaussian FWHM to sigma in the same units."""
    return np.asarray(fwhm, dtype=float) / np.sqrt(8.0 * np.log(2.0))


def gaussian_sigma_to_fwhm(sigma: np.ndarray | float) -> np.ndarray:
    """Convert Gaussian sigma to FWHM in the same units."""
    return np.asarray(sigma, dtype=float) * np.sqrt(8.0 * np.log(2.0))


def required_template_convolution_sigma(
    data_sigma: np.ndarray,
    template_sigma: np.ndarray,
) -> np.ndarray:
    r"""Return :math:`\sqrt{\sigma_\mathrm{data}^2-\sigma_\mathrm{temp}^2}`.

    Negative values indicate that the templates are broader than the data and
    cannot be matched by further smoothing the templates.  The caller must
    report that incompatibility rather than taking an absolute value.
    """
    data = np.asarray(data_sigma, dtype=float)
    template = np.asarray(template_sigma, dtype=float)
    diff2 = data**2 - template**2
    result = np.full(np.broadcast(data, template).shape, np.nan, dtype=float)
    good = diff2 >= 0
    result[good] = np.sqrt(diff2[good])
    return result


# -----------------------------------------------------------------------------
# Master-arc loading and 1-D spectrum construction
# -----------------------------------------------------------------------------


def read_primary_header(path: str | Path) -> fits.Header:
    """Read and copy the primary FITS header from a calibration/science file."""
    return fits.getheader(Path(path).expanduser().resolve(), 0).copy()


def _validate_geometry_map_header(
    map_header: fits.Header,
    arc_header: fits.Header,
    *,
    map_name: str,
) -> None:
    """Check that an explicitly supplied/discovered geometry map matches the arc.

    Filename roots alone are not reliable on every KCWI red-side reduction.  We
    therefore use ``OFNAME`` and instrumental setup metadata whenever present.
    """
    for key in ("OFNAME", "CAMERA", "IFUNAM", "BINNING", "BGRATNAM", "RGRATNAM"):
        a = str(arc_header.get(key, "")).strip().upper()
        b = str(map_header.get(key, "")).strip().upper()
        if a and b and a != b:
            raise ValueError(
                f"Master arc and {map_name} disagree in FITS {key}: {a!r} versus {b!r}. "
                "This geometry map does not appear to belong to the supplied master arc."
            )


def load_master_arc_maps(
    master_arc: str | Path,
    *,
    wavemap: str | Path | None = None,
    slicemap: str | Path | None = None,
    posmap: str | Path | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, fits.Header, dict[str, Path]]:
    """Load a master arc and its wavelength/slice/position geometry maps."""
    arc_path = Path(master_arc).expanduser().resolve()
    supplied = {"wavemap": wavemap, "slicemap": slicemap, "posmap": posmap}
    discovered = discover_arc_sidecars(arc_path) if any(v is None for v in supplied.values()) else {}
    paths = {
        key: (Path(value).expanduser().resolve() if value is not None else discovered[key])
        for key, value in supplied.items()
    }
    for key, value in paths.items():
        if not value.exists():
            raise FileNotFoundError(f"Configured/discovered master-arc {key} does not exist: {value}")

    with fits.open(arc_path, memmap=False) as hdul:
        arc = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()
    with fits.open(paths["wavemap"], memmap=False) as hdul:
        wave_header = hdul[0].header.copy()
        _validate_geometry_map_header(wave_header, header, map_name="wavemap")
        wave = np.asarray(hdul[0].data, dtype=float)
        # WAVGOOD values can live on the geometry product even if absent on the
        # master arc, so copy them into the returned header when useful.
        for key in ("WAVGOOD0", "WAVGOOD1", "WAVALL0", "WAVALL1", "WAVMID"):
            if key not in header and key in wave_header:
                header[key] = wave_header[key]
    with fits.open(paths["slicemap"], memmap=False) as hdul:
        slice_header = hdul[0].header.copy()
        _validate_geometry_map_header(slice_header, header, map_name="slicemap")
        slices = np.asarray(hdul[0].data)
    with fits.open(paths["posmap"], memmap=False) as hdul:
        pos_header = hdul[0].header.copy()
        _validate_geometry_map_header(pos_header, header, map_name="posmap")
        positions = np.asarray(hdul[0].data, dtype=float)

    if not (arc.shape == wave.shape == slices.shape == positions.shape):
        raise ValueError(
            "Master arc, wavemap, slicemap, and posmap must have identical detector shapes: "
            f"arc={arc.shape}, wave={wave.shape}, slice={slices.shape}, pos={positions.shape}"
        )
    return arc, wave, slices, positions, header, paths


def _estimate_detector_dispersion(wavemap: np.ndarray, valid: np.ndarray) -> float:
    """Estimate Angstrom per detector pixel from the axis with the larger wave gradient."""
    candidates: list[float] = []
    for axis in (0, 1):
        dw = np.diff(wavemap, axis=axis)
        pair_valid = np.diff(valid.astype(np.int8), axis=axis) == 0
        # ``pair_valid`` alone also includes adjacent invalid pairs.  Explicitly
        # build the two-side mask for clarity.
        if axis == 0:
            both = valid[:-1, :] & valid[1:, :]
        else:
            both = valid[:, :-1] & valid[:, 1:]
        vals = np.abs(dw[both & np.isfinite(dw)])
        vals = vals[vals > 0]
        if vals.size:
            candidates.append(float(np.nanmedian(vals)))
        else:
            candidates.append(np.nan)

    finite = [v for v in candidates if np.isfinite(v) and v > 0]
    if not finite:
        raise ValueError("Could not infer detector wavelength sampling from wavemap")
    # The spectral direction has a substantially larger wavelength increment;
    # the cross-dispersion direction is ideally close to zero.
    return float(max(finite))


def _collapse_arc_region(
    arc: np.ndarray,
    wavemap: np.ndarray,
    region: np.ndarray,
    *,
    wave_min: float,
    wave_max: float,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Median-bin detector pixels in one slice/position region onto wavelength."""
    use = (
        region
        & np.isfinite(arc)
        & np.isfinite(wavemap)
        & (wavemap >= wave_min)
        & (wavemap <= wave_max)
    )
    if np.sum(use) < 20:
        return np.array([]), np.array([])

    w = wavemap[use]
    f = arc[use]
    edges = np.arange(wave_min - 0.5 * step, wave_max + 1.5 * step, step)
    if edges.size < 5:
        return np.array([]), np.array([])
    idx = np.digitize(w, edges) - 1
    centers = 0.5 * (edges[:-1] + edges[1:])
    spec = np.full(centers.size, np.nan, dtype=float)

    for j in np.unique(idx[(idx >= 0) & (idx < centers.size)]):
        vals = f[idx == j]
        if vals.size:
            spec[j] = np.nanmedian(vals)

    valid = np.isfinite(spec)
    if np.sum(valid) < 10:
        return np.array([]), np.array([])

    # Interpolate only small missing gaps so Gaussian fitting sees a regular
    # wavelength grid.  Large uncovered edges remain removed by trimming to the
    # first/last valid sample.
    first, last = np.where(valid)[0][[0, -1]]
    x = centers[first : last + 1]
    y = spec[first : last + 1]
    good = np.isfinite(y)
    if np.sum(good) < 10:
        return np.array([]), np.array([])
    y = np.interp(x, x[good], y[good])
    return x, y


# -----------------------------------------------------------------------------
# Arc-line fitting
# -----------------------------------------------------------------------------


def _gaussian_linear(x: np.ndarray, amp: float, mu: float, sigma: float, c0: float, c1: float) -> np.ndarray:
    return c0 + c1 * (x - mu) + amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _fit_one_arc_line(
    wavelength: np.ndarray,
    spectrum: np.ndarray,
    peak_index: int,
    *,
    half_width_pix: int,
    sampling: float,
) -> tuple[float, float, float, float] | None:
    lo = max(0, int(peak_index) - int(half_width_pix))
    hi = min(spectrum.size, int(peak_index) + int(half_width_pix) + 1)
    if hi - lo < 7:
        return None

    x = wavelength[lo:hi]
    y = spectrum[lo:hi]
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        return None

    edge_background = float(np.median(np.r_[y[:2], y[-2:]]))
    amp0 = float(y[peak_index - lo] - edge_background)
    if not np.isfinite(amp0) or amp0 <= 0:
        return None

    mu0 = float(wavelength[peak_index])
    sigma0 = max(float(sampling), 0.5 * float(sampling))
    slope0 = 0.0

    lower = [0.0, mu0 - 1.5 * sampling, 0.15 * sampling, -np.inf, -np.inf]
    upper = [np.inf, mu0 + 1.5 * sampling, 8.0 * sampling, np.inf, np.inf]

    try:
        pars, _ = curve_fit(
            _gaussian_linear,
            x,
            y,
            p0=[amp0, mu0, sigma0, edge_background, slope0],
            bounds=(lower, upper),
            maxfev=10000,
        )
    except Exception:
        return None

    amp, mu, sigma, c0, c1 = [float(v) for v in pars]
    model = _gaussian_linear(x, *pars)
    resid = y - model
    robust_noise = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    if not np.isfinite(robust_noise) or robust_noise <= 0:
        robust_noise = np.std(resid)
    if not np.isfinite(robust_noise) or robust_noise <= 0:
        reduced_proxy = np.nan
    else:
        reduced_proxy = float(np.mean((resid / robust_noise) ** 2))

    if sigma <= 0 or not np.isfinite(sigma) or abs(mu - mu0) > 1.5 * sampling:
        return None
    return mu, sigma, amp, reduced_proxy


def _robust_polyfit(
    x: np.ndarray,
    y: np.ndarray,
    order: int,
    *,
    sigma_clip: float = 3.0,
    max_iter: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterative sigma-clipped polynomial fit; returns coefficients and keep mask."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if np.sum(keep) < order + 2:
        raise ValueError("Too few finite LSF measurements for polynomial model")

    for _ in range(max_iter):
        coeff = np.polyfit(x[keep], y[keep], order)
        resid = y - np.polyval(coeff, x)
        med = np.nanmedian(resid[keep])
        mad = np.nanmedian(np.abs(resid[keep] - med))
        scale = 1.4826 * mad
        if not np.isfinite(scale) or scale <= 0:
            break
        new_keep = keep & (np.abs(resid - med) <= sigma_clip * scale)
        if np.array_equal(new_keep, keep):
            keep = new_keep
            break
        if np.sum(new_keep) < order + 2:
            break
        keep = new_keep

    coeff = np.polyfit(x[keep], y[keep], order)
    return coeff, keep


def measure_arc_lsf(
    arc: np.ndarray,
    wavemap: np.ndarray,
    slicemap: np.ndarray,
    posmap: np.ndarray,
    *,
    wavegood0: float | None = None,
    wavegood1: float | None = None,
    polynomial_order: int = 2,
    measure_spatial_variation: bool = True,
    spatial_bins: int = 3,
    peak_prominence_fraction: float = 0.05,
    peak_height_percentile: float = 70.0,
    min_peak_distance_pix: int = 4,
    line_fit_half_width_pix: int = 5,
    min_good_lines: int = 6,
    sigma_clip: float = 3.0,
) -> ArcLSFResult:
    """Measure unresolved arc-line FWHM as a function of wavelength and slice.

    The detector data are grouped by KCWI slice and, optionally, by position
    within each slice.  Each group is rebinned using the DRP wavelength map,
    emission peaks are detected, and isolated peaks are fit with a Gaussian plus
    local linear background.  A sigma-clipped polynomial provides the baseline
    wavelength-only LSF model; retained per-slice/position measurements quantify
    spatial variation for QC and later bin-specific LSF work.
    """
    arc = np.asarray(arc, dtype=float)
    wavemap = np.asarray(wavemap, dtype=float)
    slicemap = np.asarray(slicemap)
    posmap = np.asarray(posmap, dtype=float)
    if not (arc.shape == wavemap.shape == slicemap.shape == posmap.shape) or arc.ndim != 2:
        raise ValueError("arc/wavemap/slicemap/posmap must be matching 2-D arrays")

    valid = (
        np.isfinite(arc)
        & np.isfinite(wavemap)
        & np.isfinite(posmap)
        & np.isfinite(slicemap)
        & (slicemap >= 0)
    )
    if not np.any(valid):
        raise ValueError("No valid master-arc geometry pixels")

    wmin = float(np.nanmin(wavemap[valid])) if wavegood0 is None else float(wavegood0)
    wmax = float(np.nanmax(wavemap[valid])) if wavegood1 is None else float(wavegood1)
    if wmax <= wmin:
        raise ValueError("Invalid master-arc wavelength range")

    sampling = _estimate_detector_dispersion(wavemap, valid)
    measurements: list[ArcLineMeasurement] = []
    slice_int = np.full(slicemap.shape, -999, dtype=int)
    slice_int[valid] = np.asarray(slicemap[valid], dtype=int)
    slice_values = np.unique(slice_int[valid])

    for slice_id in slice_values:
        slice_mask = valid & (slice_int == int(slice_id))
        if np.sum(slice_mask) < 20:
            continue

        if measure_spatial_variation and spatial_bins > 1:
            pos_values = posmap[slice_mask]
            edges = np.quantile(pos_values, np.linspace(0.0, 1.0, spatial_bins + 1))
            # Repeated quantile edges can occur for sparse/quantized position
            # maps.  Fall back to one full-slice region in that case.
            if np.unique(edges).size < spatial_bins + 1:
                regions = [(0, slice_mask, float(np.nanmedian(pos_values)))]
            else:
                regions = []
                for pos_bin in range(spatial_bins):
                    lo, hi = edges[pos_bin], edges[pos_bin + 1]
                    if pos_bin == spatial_bins - 1:
                        region = slice_mask & (posmap >= lo) & (posmap <= hi)
                    else:
                        region = slice_mask & (posmap >= lo) & (posmap < hi)
                    regions.append((pos_bin, region, float(np.nanmedian(posmap[region]))))
        else:
            regions = [(0, slice_mask, float(np.nanmedian(posmap[slice_mask])))]

        for pos_bin, region, pos_center in regions:
            x, y = _collapse_arc_region(
                arc,
                wavemap,
                region,
                wave_min=wmin,
                wave_max=wmax,
                step=sampling,
            )
            if x.size < 15:
                continue

            # Remove broad background structure before peak detection while
            # leaving the original y-values untouched for Gaussian fitting.
            filt_size = max(9, 2 * line_fit_half_width_pix + 5)
            if filt_size % 2 == 0:
                filt_size += 1
            baseline = median_filter(y, size=filt_size, mode="nearest")
            high = y - baseline
            positive = high[np.isfinite(high) & (high > 0)]
            if positive.size < 5:
                continue
            height = np.percentile(positive, peak_height_percentile)
            prominence = max(
                peak_prominence_fraction * float(np.nanmax(high)),
                0.1 * float(height),
            )
            peaks, _ = find_peaks(
                high,
                height=height,
                prominence=prominence,
                distance=max(1, int(min_peak_distance_pix)),
            )

            for peak in peaks:
                fit = _fit_one_arc_line(
                    x,
                    y,
                    int(peak),
                    half_width_pix=line_fit_half_width_pix,
                    sampling=sampling,
                )
                if fit is None:
                    continue
                mu, sigma, amp, chi_proxy = fit
                fwhm = float(gaussian_sigma_to_fwhm(sigma))
                # Reject clearly pathological widths.  This is deliberately
                # generous; the sigma-clipped wavelength model performs the
                # finer outlier rejection.
                if not (0.25 * sampling <= fwhm <= 20.0 * sampling):
                    continue
                measurements.append(
                    ArcLineMeasurement(
                        wavelength=mu,
                        fwhm_angstrom=fwhm,
                        sigma_angstrom=sigma,
                        amplitude=amp,
                        slice_id=int(slice_id),
                        position_bin=int(pos_bin),
                        position_center=float(pos_center),
                        reduced_chi2_proxy=float(chi_proxy),
                    )
                )

    if len(measurements) < int(min_good_lines):
        raise RuntimeError(
            f"Only {len(measurements)} usable master-arc line measurements were found; "
            f"at least {min_good_lines} are required. Inspect arc/geometry products "
            "and tune the configurable peak-detection thresholds if necessary."
        )

    wave = np.array([m.wavelength for m in measurements], dtype=float)
    fwhm = np.array([m.fwhm_angstrom for m in measurements], dtype=float)
    sigma = np.array([m.sigma_angstrom for m in measurements], dtype=float)
    sid = np.array([m.slice_id for m in measurements], dtype=int)
    pbin = np.array([m.position_bin for m in measurements], dtype=int)
    pcenter = np.array([m.position_center for m in measurements], dtype=float)
    chi = np.array([m.reduced_chi2_proxy for m in measurements], dtype=float)

    coeff, keep = _robust_polyfit(
        wave,
        fwhm,
        int(polynomial_order),
        sigma_clip=float(sigma_clip),
    )

    model = np.polyval(coeff, wave[keep])
    frac_resid = (fwhm[keep] - model) / model

    # Separate raw line-to-line scatter from coherent spatial structure.  The
    # previous prototype called the RMS of *all* residuals "spatial RMS", which
    # mixes measurement scatter, imperfect Gaussian line modeling, and genuine
    # slice/position dependence.  Real KCWI arcs show thousands of line fits, so
    # coherent group medians provide a much cleaner spatial diagnostic.
    measurement_fractional_rms = float(np.sqrt(np.nanmean(frac_resid**2)))

    sid_keep = sid[keep]
    pbin_keep = pbin[keep]
    group_medians = []
    slice_medians = []
    position_medians = []

    for slice_value in np.unique(sid_keep):
        m = sid_keep == slice_value
        if np.any(m):
            slice_medians.append(float(np.nanmedian(frac_resid[m])))
        for pos_value in np.unique(pbin_keep[m]):
            mp = m & (pbin_keep == pos_value)
            if np.any(mp):
                group_medians.append(float(np.nanmedian(frac_resid[mp])))

    for pos_value in np.unique(pbin_keep):
        m = pbin_keep == pos_value
        if np.any(m):
            position_medians.append(float(np.nanmedian(frac_resid[m])))

    def _rms(values: list[float]) -> float:
        arr = np.asarray(values, dtype=float)
        if arr.size == 0 or not np.any(np.isfinite(arr)):
            return np.nan
        return float(np.sqrt(np.nanmean(arr**2)))

    spatial_fractional_rms = _rms(group_medians)
    slice_fractional_rms = _rms(slice_medians)
    position_fractional_rms = _rms(position_medians)

    return ArcLSFResult(
        wavelength=wave[keep],
        fwhm_angstrom=fwhm[keep],
        sigma_angstrom=sigma[keep],
        slice_id=sid[keep],
        position_bin=pbin[keep],
        position_center=pcenter[keep],
        reduced_chi2_proxy=chi[keep],
        polynomial_coefficients=np.asarray(coeff, dtype=float),
        polynomial_order=int(polynomial_order),
        wavelength_min=wmin,
        wavelength_max=wmax,
        measurement_wavelength_min=float(np.nanmin(wave[keep])),
        measurement_wavelength_max=float(np.nanmax(wave[keep])),
        spectral_sampling_angstrom=float(sampling),
        spatial_fractional_rms=spatial_fractional_rms,
        measurement_fractional_rms=measurement_fractional_rms,
        slice_fractional_rms=slice_fractional_rms,
        position_fractional_rms=position_fractional_rms,
        n_lines_total=len(measurements),
        n_lines_used=int(np.sum(keep)),
    )


def measure_arc_lsf_from_files(
    master_arc: str | Path,
    **kwargs,
) -> tuple[ArcLSFResult, dict[str, Path]]:
    """File-oriented wrapper around :func:`measure_arc_lsf`."""
    wavemap = kwargs.pop("wavemap", None)
    slicemap = kwargs.pop("slicemap", None)
    posmap = kwargs.pop("posmap", None)
    arc, wave, slices, positions, header, paths = load_master_arc_maps(
        master_arc,
        wavemap=wavemap,
        slicemap=slicemap,
        posmap=posmap,
    )
    kwargs.setdefault("wavegood0", header.get("WAVGOOD0"))
    kwargs.setdefault("wavegood1", header.get("WAVGOOD1"))
    return measure_arc_lsf(arc, wave, slices, positions, **kwargs), paths


def save_arc_lsf_result(result: ArcLSFResult, npz_path: str | Path, csv_path: str | Path) -> tuple[Path, Path]:
    """Save both machine-efficient and human-readable empirical LSF products."""
    npz_path = Path(npz_path)
    csv_path = Path(csv_path)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        npz_path,
        wavelength=result.wavelength,
        fwhm_angstrom=result.fwhm_angstrom,
        sigma_angstrom=result.sigma_angstrom,
        slice_id=result.slice_id,
        position_bin=result.position_bin,
        position_center=result.position_center,
        reduced_chi2_proxy=result.reduced_chi2_proxy,
        polynomial_coefficients=result.polynomial_coefficients,
        polynomial_order=result.polynomial_order,
        wavelength_min=result.wavelength_min,
        wavelength_max=result.wavelength_max,
        measurement_wavelength_min=result.measurement_wavelength_min,
        measurement_wavelength_max=result.measurement_wavelength_max,
        blue_edge_unconstrained_angstrom=result.blue_edge_unconstrained_angstrom,
        red_edge_unconstrained_angstrom=result.red_edge_unconstrained_angstrom,
        spectral_sampling_angstrom=result.spectral_sampling_angstrom,
        spatial_fractional_rms=result.spatial_fractional_rms,
        measurement_fractional_rms=result.measurement_fractional_rms,
        slice_fractional_rms=result.slice_fractional_rms,
        position_fractional_rms=result.position_fractional_rms,
        n_lines_total=result.n_lines_total,
        n_lines_used=result.n_lines_used,
    )

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "wavelength_A",
                "fwhm_A",
                "sigma_A",
                "slice_id",
                "position_bin",
                "position_center",
                "reduced_chi2_proxy",
            ]
        )
        for row in zip(
            result.wavelength,
            result.fwhm_angstrom,
            result.sigma_angstrom,
            result.slice_id,
            result.position_bin,
            result.position_center,
            result.reduced_chi2_proxy,
        ):
            writer.writerow(row)
    return npz_path, csv_path


# -----------------------------------------------------------------------------
# PSF helpers
# -----------------------------------------------------------------------------


def _header_psf_candidate(header: fits.Header, keys: tuple[str, ...]) -> tuple[float, str] | None:
    for key in keys:
        if key not in header:
            continue
        try:
            value = float(header[key])
        except Exception:
            continue
        if np.isfinite(value) and value > 0:
            return value, key
    return None


def estimate_psf(
    *,
    configured_fwhm_arcsec: float | None,
    header: fits.Header,
    header_keys: tuple[str, ...] = ("SEEING", "FWHM", "GUIDFWHM"),
) -> PSFEstimate:
    """Return the best available PSF estimate without fitting the galaxy itself."""
    if configured_fwhm_arcsec is not None:
        value = float(configured_fwhm_arcsec)
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Configured PSF FWHM must be positive and finite")
        return PSFEstimate(value, "configuration", "Explicit target configuration value")

    candidate = _header_psf_candidate(header, tuple(header_keys))
    if candidate is not None:
        value, key = candidate
        return PSFEstimate(value, "header", f"Primary-header keyword {key}")

    return PSFEstimate(
        np.nan,
        "unavailable",
        "No explicit PSF was configured and no recognized seeing/FWHM header keyword was found",
    )
