"""XSL SSP template preparation for RH3 kinematics and BL populations.

Script 3 deliberately gives both stellar components the *same complete XSL SSP
basis*.  Before duplicating that basis, each SSP spectrum is independently
normalized to unit mean flux density over the fixed RH3 fitting interval.  This
normalization is what makes pPXF's two-component ``fraction`` constraint an
explicit passband light fraction for the RH3 likelihood cube rather than a
ratio contaminated by the arbitrary native normalization of different SSPs.

Important interpretation
------------------------
After this normalization, ``f_A`` in Script 3 is the fraction of the fitted
stellar-template light assigned to component A over ``RH3_FIT_REST_RANGE_ANGSTROM``.
It is *not* a mass fraction, and it is not expected to equal the BL light
fraction.  The additive polynomial is a nuisance continuum term and is not
included in the definition of component stellar light.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

C_KMS = 299792.458


@dataclass
class PreparedTemplates:
    templates: np.ndarray
    wavelength: np.ndarray
    velscale: float
    native_template_fwhm: np.ndarray
    target_fwhm: np.ndarray
    convolution_fwhm: np.ndarray
    age: np.ndarray | None
    metallicity: np.ndarray | None
    normalization: np.ndarray
    normalization_range: tuple[float, float]
    template_medium: str


def _morton2000_refractive_index(vacuum_angstrom: np.ndarray) -> np.ndarray:
    """Dry-air refractive index using the Morton (2000) optical expression."""
    wave = np.asarray(vacuum_angstrom, dtype=float)
    if np.any(wave < 2000.0):
        raise ValueError("Morton2000 air/vacuum conversion is restricted here to lambda >= 2000 A.")
    sigma2 = (1.0e4 / wave) ** 2  # inverse micrometre squared
    return (
        1.0
        + 8.34254e-5
        + 2.406147e-2 / (130.0 - sigma2)
        + 1.5998e-4 / (38.9 - sigma2)
    )


def vacuum_to_air(wavelength_angstrom: np.ndarray) -> np.ndarray:
    """Convert optical vacuum wavelengths to standard-air wavelengths."""
    wave = np.asarray(wavelength_angstrom, dtype=float)
    return wave / _morton2000_refractive_index(wave)


def air_to_vacuum(wavelength_angstrom: np.ndarray, *, maxiter: int = 12) -> np.ndarray:
    """Convert optical air wavelengths to vacuum by fixed-point iteration."""
    air = np.asarray(wavelength_angstrom, dtype=float)
    vac = air * 1.00028
    for _ in range(int(maxiter)):
        new = air * _morton2000_refractive_index(vac)
        if np.all(np.abs(new - vac) <= 1.0e-10 * np.maximum(vac, 1.0)):
            vac = new
            break
        vac = new
    return vac


def convert_wavelength_medium(
    wavelength_angstrom: np.ndarray,
    from_medium: str,
    to_medium: str,
) -> np.ndarray:
    """Convert an optical wavelength vector between air and vacuum."""
    src = str(from_medium).strip().lower()
    dst = str(to_medium).strip().lower()
    if src not in {"air", "vacuum"} or dst not in {"air", "vacuum"}:
        raise ValueError("Wavelength medium must be 'air' or 'vacuum'.")
    wave = np.asarray(wavelength_angstrom, dtype=float)
    if src == dst:
        return wave.copy()
    if src == "vacuum" and dst == "air":
        return vacuum_to_air(wave)
    return air_to_vacuum(wave)


def infer_science_wavelength_medium(header) -> str:
    """Infer air/vacuum convention from a prepared KCWI FITS header.

    Script 1 is the authoritative validator.  This function is intentionally
    conservative and raises when the prepared product does not expose enough
    metadata to reconstruct the convention.
    """
    # Explicit CRD_DAP/header conventions take priority when present.
    for key in ("WAVEMED", "WAVEMEDM", "SPECMED", "MEDIUM"):
        value = header.get(key)
        if value is not None:
            text = str(value).strip().lower()
            if "vac" in text:
                return "vacuum"
            if "air" in text:
                return "air"

    for idx in (1, 2, 3):
        ctype = str(header.get(f"CTYPE{idx}", "")).upper()
        comment = ""
        try:
            comment = str(header.comments[f"CTYPE{idx}"]).upper()
        except Exception:
            pass
        if "AWAV" in ctype or "AIR" in comment:
            return "air"
        if "VAC" in ctype or "VACUUM" in comment:
            return "vacuum"

    raise ValueError(
        "Could not infer the science wavelength medium from the prepared RH3 "
        "FITS header. Script 3 will not guess because an air/vacuum mismatch "
        "can mimic a velocity zero-point offset."
    )


def _pixel_edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 1 or centers.size < 2 or np.any(np.diff(centers) <= 0):
        raise ValueError("Wavelength centers must be a strictly increasing 1-D array.")
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def make_log_wavelength_grid(
    minimum: float,
    maximum: float,
    *,
    velscale: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return logarithmically spaced pixel centers and edges."""
    if not (minimum > 0 and maximum > minimum and velscale > 0):
        raise ValueError("Invalid log-wavelength grid bounds/velscale.")
    dlog = float(velscale) / C_KMS
    lo = np.log(float(minimum))
    hi = np.log(float(maximum))
    n = int(np.floor((hi - lo) / dlog))
    if n < 10:
        raise ValueError("Requested logarithmic wavelength grid is too short.")
    log_edges = lo + np.arange(n + 1, dtype=float) * dlog
    centers = np.exp(0.5 * (log_edges[:-1] + log_edges[1:]))
    return centers, np.exp(log_edges)


def _rebin_density_to_edges(
    wave: np.ndarray,
    values: np.ndarray,
    out_edges: np.ndarray,
) -> np.ndarray:
    """Flux-density preserving overlap rebin for one or many columns."""
    wave = np.asarray(wave, dtype=float)
    arr = np.asarray(values, dtype=float)
    if arr.shape[0] != wave.size:
        raise ValueError("values first axis must match wave.")
    in_edges = _pixel_edges(wave)
    nout = len(out_edges) - 1
    tail_shape = arr.shape[1:]
    out = np.full((nout,) + tail_shape, np.nan, dtype=float)

    i = 0
    for j in range(nout):
        left, right = out_edges[j], out_edges[j + 1]
        while i + 1 < in_edges.size and in_edges[i + 1] <= left:
            i += 1
        k = i
        pieces = []
        weights = []
        while k < wave.size and in_edges[k] < right:
            overlap = max(0.0, min(right, in_edges[k + 1]) - max(left, in_edges[k]))
            if overlap > 0:
                pieces.append(arr[k])
                weights.append(overlap)
            if in_edges[k + 1] >= right:
                break
            k += 1
        if weights:
            w = np.asarray(weights, dtype=float)
            p = np.asarray(pieces, dtype=float)
            out[j] = np.tensordot(w / np.sum(w), p, axes=(0, 0))
    return out


def normalize_template_light_fraction_basis(
    templates: np.ndarray,
    wavelength: np.ndarray,
    normalization_range: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize each SSP independently to unit mean flux over one fixed band.

    This is the defining convention for ``f_A,RH3`` in Script 3.  Independent
    SSP normalization is deliberate here: Script 3 needs a transparent
    passband-light fraction, not mass-normalized population weights.  Detailed
    population inference occurs later in the BL stage with its own normalization
    convention.
    """
    wave = np.asarray(wavelength, dtype=float)
    temp = np.asarray(templates, dtype=float)
    lo, hi = map(float, normalization_range)
    use = (wave >= lo) & (wave <= hi)
    if np.sum(use) < 10:
        raise ValueError("Template normalization range contains fewer than 10 log pixels.")
    means = np.nanmean(temp[use, :], axis=0)
    bad = (~np.isfinite(means)) | (means <= 0)
    if np.any(bad):
        ids = np.flatnonzero(bad)[:10].tolist()
        raise ValueError(
            "One or more XSL SSP templates have non-positive/non-finite mean "
            f"flux in the RH3 normalization band; example columns: {ids}."
        )
    return temp / means[None, :], means


def _load_xsl_arrays(path: str | Path):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as data:
        keys = set(data.files)
        if "lam" not in keys or "templates" not in keys:
            raise KeyError(
                f"XSL NPZ must contain 'lam' and 'templates'; found keys={sorted(keys)}"
            )
        lam = np.asarray(data["lam"], dtype=float)
        templates_nd = np.asarray(data["templates"], dtype=float)
        if templates_nd.shape[0] != lam.size:
            raise ValueError("XSL templates first axis must match data['lam'].")
        templates = templates_nd.reshape(lam.size, -1)

        if "fwhm" not in keys:
            raise KeyError(
                "XSL NPZ does not contain its native 'fwhm' resolution array. "
                "CRD_DAP will not substitute a nominal resolving power."
            )
        native_fwhm = np.asarray(data["fwhm"], dtype=float)
        if native_fwhm.ndim == 0:
            native_fwhm = np.full(lam.size, float(native_fwhm))
        elif native_fwhm.size != lam.size:
            native_fwhm = np.full(lam.size, float(np.nanmedian(native_fwhm)))

        age = None
        metallicity = None
        for key in ("age_grid", "ages", "age"):
            if key in keys:
                candidate = np.asarray(data[key], dtype=float)
                if candidate.size == templates.shape[1]:
                    age = candidate.reshape(-1)
                break
        for key in ("metal_grid", "metals", "metallicity", "metal"):
            if key in keys:
                candidate = np.asarray(data[key], dtype=float)
                if candidate.size == templates.shape[1]:
                    metallicity = candidate.reshape(-1)
                break
    return lam, templates, native_fwhm, age, metallicity


def _variable_gaussian_smooth(
    wavelength: np.ndarray,
    templates: np.ndarray,
    sigma_angstrom: np.ndarray,
) -> np.ndarray:
    """Apply a wavelength-dependent Gaussian kernel using pPXF ``varsmooth``."""
    try:
        from ppxf.ppxf_util import varsmooth
    except Exception as exc:  # pragma: no cover - optional science dependency
        raise ImportError("Template LSF matching requires pPXF's ppxf_util.varsmooth.") from exc

    wave = np.asarray(wavelength, dtype=float)
    temp = np.asarray(templates, dtype=float)
    sig = np.asarray(sigma_angstrom, dtype=float)
    if np.any(~np.isfinite(sig)) or np.any(sig < 0):
        raise ValueError("Variable smoothing sigma must be finite and non-negative.")

    smoothed = np.empty_like(temp, dtype=float)
    for j in range(temp.shape[1]):
        smoothed[:, j] = np.asarray(varsmooth(wave, temp[:, j], sig), dtype=float)
    return smoothed


def prepare_xsl_rh3_templates(
    *,
    xsl_path: str | Path,
    fit_range: tuple[float, float],
    velscale: float,
    target_fwhm_rest: Callable[[np.ndarray], np.ndarray],
    template_medium: str,
    velocity_padding_kms: float,
) -> PreparedTemplates:
    """Prepare the complete XSL SSP grid for Script-3 RH3 likelihood fitting.

    Parameters
    ----------
    target_fwhm_rest
        Callable evaluated on template-medium *rest-frame* wavelengths and
        returning the measured galaxy instrumental FWHM in the same frame and
        wavelength medium. Values outside Script-1 empirical LSF support must be
        NaN; this routine then refuses to extrapolate.
    """
    lo, hi = map(float, fit_range)
    if hi <= lo:
        raise ValueError("fit_range must be increasing.")
    pad_factor = np.exp(float(velocity_padding_kms) / C_KMS)
    padded_lo = lo / pad_factor
    padded_hi = hi * pad_factor

    lam, raw, native_fwhm, age, metal = _load_xsl_arrays(xsl_path)
    keep = (lam >= padded_lo - 5.0) & (lam <= padded_hi + 5.0)
    if np.sum(keep) < 30:
        raise ValueError(
            "XSL library does not provide enough wavelength coverage for the RH3 "
            "fit plus configured velocity/dispersion padding."
        )
    lam = lam[keep]
    raw = raw[keep, :]
    native_fwhm = native_fwhm[keep]

    target = np.asarray(target_fwhm_rest(lam), dtype=float)
    if target.shape != lam.shape:
        raise ValueError("target_fwhm_rest must return one FWHM value per wavelength.")
    finite = np.isfinite(target)
    if not np.all(finite):
        bad_lo = float(np.nanmin(lam[~finite])) if np.any(~finite) else np.nan
        bad_hi = float(np.nanmax(lam[~finite])) if np.any(~finite) else np.nan
        raise ValueError(
            "Requested template/fit wavelength support extends outside the empirically "
            f"measured Script-1 LSF (unsupported template rest wavelengths include {bad_lo:.1f}--{bad_hi:.1f} A)."
        )

    diff2 = target**2 - native_fwhm**2
    scale = np.maximum(target**2, native_fwhm**2)
    materially_negative = diff2 < -1.0e-6 * np.maximum(scale, 1.0e-12)
    if np.any(materially_negative):
        ii = np.flatnonzero(materially_negative)
        raise ValueError(
            "XSL templates are lower spectral resolution than the observed galaxy LSF "
            "at part of the requested range; a real-valued convolution kernel does not "
            f"exist. First offending rest wavelength={lam[ii[0]]:.2f} A, "
            f"template FWHM={native_fwhm[ii[0]]:.4f} A, target FWHM={target[ii[0]]:.4f} A."
        )
    convolution_fwhm = np.sqrt(np.clip(diff2, 0.0, None))
    sigma_add = convolution_fwhm / 2.354820045
    smoothed = _variable_gaussian_smooth(lam, raw, sigma_add)

    log_wave, log_edges = make_log_wavelength_grid(
        padded_lo, padded_hi, velscale=float(velscale)
    )
    log_templates = _rebin_density_to_edges(lam, smoothed, log_edges)
    if np.any(~np.isfinite(log_templates)):
        raise ValueError("Template log-rebin produced unsupported/non-finite pixels.")

    # The light-fraction convention is defined on the *science fitting interval*,
    # not the wider velocity-padding interval.
    log_templates, norms = normalize_template_light_fraction_basis(
        log_templates, log_wave, (lo, hi)
    )
    target_log = np.interp(log_wave, lam, target)
    native_log = np.interp(log_wave, lam, native_fwhm)
    conv_log = np.interp(log_wave, lam, convolution_fwhm)

    return PreparedTemplates(
        templates=log_templates,
        wavelength=log_wave,
        velscale=float(velscale),
        native_template_fwhm=native_log,
        target_fwhm=target_log,
        convolution_fwhm=conv_log,
        age=age,
        metallicity=metal,
        normalization=norms,
        normalization_range=(lo, hi),
        template_medium=str(template_medium).lower(),
    )


def prepare_xsl_bl_templates(*args, **kwargs):
    """Prepare BL population templates (implemented with Script 6)."""
    raise NotImplementedError("Implemented with Script 6.")


def rebin_spectrum_with_diagonal_noise(
    *,
    wavelength: np.ndarray,
    flux: np.ndarray,
    uncertainty: np.ndarray,
    good: np.ndarray,
    out_wavelength: np.ndarray,
    out_edges: np.ndarray,
    min_valid_fraction: float = 0.80,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Overlap-rebin one spectrum and propagate *formal diagonal* variance.

    The input flux is treated as a flux density.  For each log-wavelength output
    pixel, valid native samples are combined with wavelength-overlap weights.
    Formal variance is propagated with the squared normalized weights.  This is
    intentionally only the diagonal part of the rebinned noise model; the
    interpolation/resampling covariance is recorded by Script 3 as unresolved
    until residual-based calibration is performed.

    Returns
    -------
    flux_log, uncertainty_log, good_log, valid_fraction
    """
    wave = np.asarray(wavelength, dtype=float)
    flux = np.asarray(flux, dtype=float)
    unc = np.asarray(uncertainty, dtype=float)
    good = np.asarray(good, dtype=bool)
    out_wave = np.asarray(out_wavelength, dtype=float)
    out_edges = np.asarray(out_edges, dtype=float)
    if not (wave.ndim == flux.ndim == unc.ndim == good.ndim == 1):
        raise ValueError("wavelength, flux, uncertainty, and good must be 1-D.")
    if not (wave.size == flux.size == unc.size == good.size):
        raise ValueError("Input spectrum arrays must have identical lengths.")
    if out_edges.size != out_wave.size + 1:
        raise ValueError("out_edges must have exactly one more element than out_wavelength.")
    if not 0.0 < float(min_valid_fraction) <= 1.0:
        raise ValueError("min_valid_fraction must lie in (0, 1].")

    in_edges = _pixel_edges(wave)
    nout = out_wave.size
    fout = np.full(nout, np.nan, dtype=float)
    vout = np.full(nout, np.nan, dtype=float)
    coverage = np.zeros(nout, dtype=float)

    i = 0
    for j in range(nout):
        left, right = out_edges[j], out_edges[j + 1]
        width = right - left
        while i + 1 < in_edges.size and in_edges[i + 1] <= left:
            i += 1
        k = i
        weighted_flux = 0.0
        weighted_var = 0.0
        valid_width = 0.0
        while k < wave.size and in_edges[k] < right:
            overlap = max(0.0, min(right, in_edges[k + 1]) - max(left, in_edges[k]))
            is_valid = (
                overlap > 0
                and good[k]
                and np.isfinite(flux[k])
                and np.isfinite(unc[k])
                and unc[k] > 0
            )
            if is_valid:
                valid_width += overlap
                weighted_flux += overlap * flux[k]
                weighted_var += (overlap * unc[k]) ** 2
            if in_edges[k + 1] >= right:
                break
            k += 1
        if width > 0:
            coverage[j] = valid_width / width
        if valid_width > 0:
            fout[j] = weighted_flux / valid_width
            vout[j] = weighted_var / valid_width**2

    good_out = (
        (coverage >= float(min_valid_fraction))
        & np.isfinite(fout)
        & np.isfinite(vout)
        & (vout > 0)
    )
    return fout, np.sqrt(vout), good_out, coverage


@dataclass
class SavedLSFModel:
    """Numerical wavelength-only LSF model recovered from a Script-1 product."""

    wavelength_observed: np.ndarray
    fwhm_observed: np.ndarray
    empirical_min: float
    empirical_max: float
    source_path: Path

    def evaluate(self, wavelength_observed: np.ndarray) -> np.ndarray:
        wave = np.asarray(wavelength_observed, dtype=float)
        out = np.full(wave.shape, np.nan, dtype=float)
        use = (
            np.isfinite(wave)
            & (wave >= self.empirical_min)
            & (wave <= self.empirical_max)
        )
        if np.any(use):
            out[use] = np.interp(
                wave[use], self.wavelength_observed, self.fwhm_observed,
                left=np.nan, right=np.nan,
            )
        return out


def discover_script1_lsf_product(
    script1_run: str | Path,
    *,
    arm: str = "RH3",
    explicit: str | Path | None = None,
) -> Path:
    """Locate the saved numerical Script-1 LSF product without guessing science.

    An explicit path wins.  Otherwise conventional historical names are tried
    first and then the products directory is searched for exactly one NPZ whose
    filename contains both the arm name and ``lsf``.
    """
    if explicit is not None:
        path = Path(explicit).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    products = Path(script1_run).expanduser().resolve() / "products"
    arm = str(arm).upper()
    conventional = [
        products / f"lsf_{arm}.npz",
        products / f"{arm}_lsf.npz",
        products / f"{arm}_LSF.npz",
        products / f"LSF_{arm}.npz",
    ]
    for path in conventional:
        if path.is_file():
            return path

    candidates = sorted(
        p for p in products.glob("*.npz")
        if "lsf" in p.name.lower() and arm.lower() in p.name.lower()
    )
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(
            f"No numerical {arm} LSF .npz product found in {products}. "
            "Script 3 requires Script 1's empirically measured LSF and will not "
            "substitute a nominal resolving power."
        )
    raise RuntimeError(
        f"Multiple possible {arm} LSF products found in {products}: "
        + ", ".join(p.name for p in candidates)
        + ". Set RH3_LSF_PRODUCT explicitly in the target configuration."
    )


def _pick_key(keys: set[str], candidates: tuple[str, ...]) -> str | None:
    lower = {k.lower(): k for k in keys}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def load_saved_lsf_model(path: str | Path, *, polynomial_order: int = 2) -> SavedLSFModel:
    """Load the numerical Script-1 LSF without extrapolating its empirical support.

    The current CRD_DAP Script-1 product stores the accepted arc-line samples,
    the fitted polynomial coefficients, and ``measurement_wavelength_min/max``.
    Script 3 reconstructs that same polynomial only inside those measured
    boundaries.  A few older development key names are accepted as fallbacks,
    but a nominal resolving power is never substituted.
    """
    path = Path(path).expanduser().resolve()
    with np.load(path, allow_pickle=False) as data:
        keys = set(data.files)

        # Current CRD_DAP Script-1 format: this is the authoritative path.
        coeff_key = _pick_key(keys, ("polynomial_coefficients", "poly_coefficients", "coefficients"))
        empirical_min_key = _pick_key(keys, (
            "measurement_wavelength_min", "empirical_wavelength_min",
            "empirical_wave_min", "empirical_min", "empirical_lsf_min",
        ))
        empirical_max_key = _pick_key(keys, (
            "measurement_wavelength_max", "empirical_wavelength_max",
            "empirical_wave_max", "empirical_max", "empirical_lsf_max",
        ))

        if coeff_key is not None and empirical_min_key is not None and empirical_max_key is not None:
            coeff = np.asarray(data[coeff_key], dtype=float).reshape(-1)
            empirical_min = float(np.asarray(data[empirical_min_key]).reshape(-1)[0])
            empirical_max = float(np.asarray(data[empirical_max_key]).reshape(-1)[0])
            if not (np.isfinite(empirical_min) and np.isfinite(empirical_max) and empirical_max > empirical_min):
                raise ValueError(f"Invalid empirical LSF support in {path.name}.")
            if coeff.size < 1 or np.any(~np.isfinite(coeff)):
                raise ValueError(f"Invalid LSF polynomial coefficients in {path.name}.")
            model_wave = np.linspace(empirical_min, empirical_max, 2048)
            model_fwhm = np.polyval(coeff, model_wave)
        else:
            # Backward-compatible fallback: use a saved model grid if one exists.
            wave_model_key = _pick_key(keys, (
                "model_wavelength", "model_wave", "wavelength_model",
                "wavelength_grid", "wave_grid",
            ))
            fwhm_model_key = _pick_key(keys, (
                "model_fwhm", "fwhm_model", "fwhm_grid",
            ))
            model_wave = None
            model_fwhm = None
            if wave_model_key is not None and fwhm_model_key is not None:
                w = np.asarray(data[wave_model_key], dtype=float).reshape(-1)
                f = np.asarray(data[fwhm_model_key], dtype=float).reshape(-1)
                if w.size == f.size and w.size >= 4:
                    model_wave, model_fwhm = w, f

            # Final fallback: reconstruct a polynomial from accepted line samples.
            wave_line_key = _pick_key(keys, (
                "wavelength", "line_wavelength", "line_wavelengths",
                "accepted_wavelength", "accepted_wavelengths", "wavelengths", "wave",
            ))
            fwhm_line_key = _pick_key(keys, (
                "fwhm_angstrom", "line_fwhm", "line_fwhms",
                "accepted_fwhm", "accepted_fwhms", "fwhm_values", "fwhm",
            ))
            line_wave = None
            line_fwhm = None
            if wave_line_key is not None and fwhm_line_key is not None:
                w = np.asarray(data[wave_line_key], dtype=float).reshape(-1)
                f = np.asarray(data[fwhm_line_key], dtype=float).reshape(-1)
                if w.size == f.size:
                    valid = np.isfinite(w) & np.isfinite(f) & (f > 0)
                    if np.sum(valid) >= max(4, int(polynomial_order) + 2):
                        line_wave, line_fwhm = w[valid], f[valid]

            if empirical_min_key is not None and empirical_max_key is not None:
                empirical_min = float(np.asarray(data[empirical_min_key]).reshape(-1)[0])
                empirical_max = float(np.asarray(data[empirical_max_key]).reshape(-1)[0])
            elif line_wave is not None:
                empirical_min = float(np.min(line_wave))
                empirical_max = float(np.max(line_wave))
            elif model_wave is not None:
                finite = np.isfinite(model_wave) & np.isfinite(model_fwhm)
                empirical_min = float(np.min(model_wave[finite]))
                empirical_max = float(np.max(model_wave[finite]))
            else:
                raise KeyError(
                    f"Could not identify empirical LSF support in {path.name}; keys={sorted(keys)}"
                )

            if model_wave is None or model_fwhm is None:
                if line_wave is None:
                    raise KeyError(
                        f"Could not identify wavelength/FWHM arrays in {path.name}; keys={sorted(keys)}"
                    )
                order = min(int(polynomial_order), max(1, line_wave.size - 2))
                coeff = np.polyfit(line_wave, line_fwhm, order)
                model_wave = np.linspace(empirical_min, empirical_max, 2048)
                model_fwhm = np.polyval(coeff, model_wave)

    valid = (
        np.isfinite(model_wave)
        & np.isfinite(model_fwhm)
        & (model_fwhm > 0)
        & (model_wave >= empirical_min)
        & (model_wave <= empirical_max)
    )
    if np.sum(valid) < 4:
        raise ValueError(f"Saved LSF model {path} has fewer than four valid empirical samples.")
    order = np.argsort(model_wave[valid])
    return SavedLSFModel(
        wavelength_observed=np.asarray(model_wave[valid][order], dtype=float),
        fwhm_observed=np.asarray(model_fwhm[valid][order], dtype=float),
        empirical_min=float(empirical_min),
        empirical_max=float(empirical_max),
        source_path=path,
    )

def observed_lsf_to_template_rest_fwhm(
    rest_wavelength_template_medium: np.ndarray,
    *,
    lsf_model: SavedLSFModel,
    redshift: float,
    science_medium: str,
    template_medium: str,
) -> np.ndarray:
    """Evaluate observed-frame LSF and express its FWHM in template rest frame."""
    lam_t = np.asarray(rest_wavelength_template_medium, dtype=float)
    rest_science = convert_wavelength_medium(lam_t, template_medium, science_medium)
    obs_science = rest_science * (1.0 + float(redshift))
    fwhm_obs = lsf_model.evaluate(obs_science)

    result = np.full(lam_t.shape, np.nan, dtype=float)
    good = np.isfinite(fwhm_obs) & (fwhm_obs > 0)
    if np.any(good):
        lo_obs = obs_science[good] - 0.5 * fwhm_obs[good]
        hi_obs = obs_science[good] + 0.5 * fwhm_obs[good]
        lo_rest_sci = lo_obs / (1.0 + float(redshift))
        hi_rest_sci = hi_obs / (1.0 + float(redshift))
        lo_t = convert_wavelength_medium(lo_rest_sci, science_medium, template_medium)
        hi_t = convert_wavelength_medium(hi_rest_sci, science_medium, template_medium)
        result[good] = hi_t - lo_t
    return result
