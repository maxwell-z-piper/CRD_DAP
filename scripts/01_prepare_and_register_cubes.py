#!/usr/bin/env python3
"""CRD_DAP Script 1: prepare and register the BL and RH3 science cubes.

This stage performs the data-quality and calibration bookkeeping that every
later likelihood calculation relies on.  It intentionally does *not* fit the
counter-rotating stellar components.  Its job is to make the two reduced KCWI
cubes scientifically safe and mutually interpretable before PowerBin.

Primary responsibilities
------------------------
1. Load matched KcwiKit four-file BL/RH3 stacks (or legacy DRP cubes) and
   standardize their internal axis order.
2. Build hard sample/spaxel/wavelength masks from variance, stack-mask,
   effective-exposure, wavelength-validity information, and the validated
   CRD_DRP atmospheric wavelength masks.
3. Create robust collapsed-continuum images and center diagnostics.
4. Verify BL/RH3 WCS registration without resampling the science cubes.
5. Save common tangent-plane spatial coordinates for both arms.
6. Measure empirical instrumental LSFs from the required master arcs and their
   DRP wavelength/slice/position maps.
7. Record the best available PSF information without attempting to infer seeing
   from an extended galaxy image.
8. Run preliminary variance-scale and spectral-correlation diagnostics.
9. Save prepared FITS products, numerical metadata, plots, and full logs.

Important statistical note
--------------------------
The Script-1 spectral-covariance diagnostic is preliminary.  It uses high-pass
residuals because no pPXF stellar model exists yet.  The final noise/covariance
model must be revisited using stellar-fit residuals before profile-likelihood
widths are interpreted quantitatively.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import hashlib
import platform
import sys
import time

import numpy as np
from astropy.io import fits

import crd_utils as crd
from crd_utils import cube_utils, io, noise, plotting, psf_lsf, validation


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare/register stacked BL and RH3 cubes for CRD_DAP."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to a target-specific config derived from config/target_config_template.py",
    )
    parser.add_argument(
        "--reduction-manifest",
        default=None,
        help=(
            "CRD_DRP_reduction_manifest.json produced by the validated CRD_DRP "
            "handoff. When supplied, its BLUE/RED i/v/m/e stacks and atmospheric "
            "masks replace the science-cube paths in the target config."
        ),
    )
    parser.add_argument(
        "--skip-reduction-hash-check",
        action="store_true",
        help=(
            "Skip SHA256 verification of CRD_DRP files against the manifest. "
            "Production runs should normally leave verification enabled."
        ),
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional explicit run-directory name. Default is TARGET_YYYYMMDD_HHMMSS.",
    )
    return parser


def _config_path(cfg, name: str) -> Path | None:
    value = getattr(cfg, name, None)
    if value is None:
        return None
    return Path(value).expanduser().resolve()


def _quality_flag(flags: list[str], name: str, logger, message: str) -> None:
    if name not in flags:
        flags.append(name)
    logger.warning("%s | %s", name, message)


def _save_spatial_coordinates(
    cube: io.KCWICube,
    center_sky,
    destination: Path,
) -> tuple[np.ndarray, np.ndarray]:
    x_arcsec, y_arcsec = cube_utils.spatial_offset_grids_arcsec(
        cube.celestial_wcs,
        cube.good_spaxel.shape,
        center_sky,
    )
    np.savez_compressed(
        destination,
        x_arcsec=x_arcsec,
        y_arcsec=y_arcsec,
        good_spaxel=cube.good_spaxel,
        center_ra_deg=float(center_sky.ra.deg),
        center_dec_deg=float(center_sky.dec.deg),
    )
    return x_arcsec, y_arcsec



def _sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    """Return a streaming SHA256 digest without loading a cube into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(block_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_record_path(record, *, label: str) -> tuple[Path, str | None]:
    """Resolve one CRD_DRP manifest file record and return path + expected hash."""
    if isinstance(record, str):
        path = Path(record).expanduser().resolve()
        expected = None
    elif isinstance(record, dict):
        value = record.get("path")
        if value in (None, ""):
            raise ValueError(f"CRD_DRP manifest record {label!r} has no path.")
        path = Path(value).expanduser().resolve()
        expected = record.get("sha256")
    else:
        raise ValueError(f"CRD_DRP manifest record {label!r} has invalid structure.")

    if not path.is_file():
        raise FileNotFoundError(f"CRD_DRP manifest file does not exist: {label}={path}")
    return path, None if expected in (None, "") else str(expected).lower()


def _load_reduction_manifest(
    manifest_path: str | Path,
    *,
    verify_hashes: bool,
    logger,
) -> tuple[Path, dict, dict[str, dict[str, object]]]:
    """Validate the CRD_DRP handoff and resolve BL/RH3 science inputs.

    CRD_DRP calls the arms BLUE/RED while CRD_DAP's analysis streams are BL/RH3.
    The mapping is therefore explicit and recorded in Script-1 provenance.
    """
    path = Path(manifest_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"CRD_DRP reduction manifest does not exist: {path}")

    manifest = dict(io.read_json(path))
    if str(manifest.get("pipeline", "")).strip() != "CRD_DRP":
        raise ValueError(f"{path} is not identified as a CRD_DRP manifest.")
    schema = str(manifest.get("schema", "")).strip()
    if not schema.startswith("CRD_DRP_reduction_manifest"):
        raise ValueError(
            f"Unsupported/unknown CRD_DRP manifest schema {schema!r}; "
            "expected CRD_DRP_reduction_manifest_*."
        )
    if str(manifest.get("validation_status", "")).upper() != "PASS":
        raise RuntimeError(
            "CRD_DRP reduction package is not validated PASS. "
            f"validation_status={manifest.get('validation_status')!r}"
        )

    arm_map = {"BL": "blue", "RH3": "red"}
    resolved: dict[str, dict[str, object]] = {}

    for dap_arm, drp_arm in arm_map.items():
        arm_payload = manifest.get(drp_arm)
        if not isinstance(arm_payload, dict):
            raise ValueError(f"CRD_DRP manifest lacks arm block {drp_arm!r}.")
        if str(arm_payload.get("status", "")).upper() != "PASS":
            raise RuntimeError(
                f"CRD_DRP {drp_arm.upper()} arm is not PASS: "
                f"{arm_payload.get('status')!r}"
            )
        files = arm_payload.get("files")
        if not isinstance(files, dict):
            raise ValueError(f"CRD_DRP {drp_arm} block lacks a files dictionary.")

        role_map = {
            "icube": "icube",
            "vcube": "vcube",
            "mcube": "mcube",
            "ecube": "ecube",
            "atmospheric_mask": "atmospheric_mask_fits",
            "atmospheric_intervals": "atmospheric_intervals_ecsv",
        }
        arm_resolved: dict[str, object] = {"manifest_arm": drp_arm}
        for role, manifest_key in role_map.items():
            rec = files.get(manifest_key)
            if rec is None:
                if role == "atmospheric_intervals":
                    arm_resolved[role] = None
                    continue
                raise ValueError(
                    f"CRD_DRP {drp_arm} files block lacks required {manifest_key!r}."
                )
            file_path, expected_hash = _manifest_record_path(
                rec, label=f"{drp_arm}.{manifest_key}"
            )
            if verify_hashes and expected_hash is not None:
                actual = _sha256_file(file_path)
                if actual.lower() != expected_hash:
                    raise RuntimeError(
                        "CRD_DRP handoff hash mismatch: "
                        f"{drp_arm}.{manifest_key}={file_path}; "
                        f"manifest={expected_hash}, actual={actual}"
                    )
            arm_resolved[role] = file_path
            arm_resolved[f"{role}_sha256"] = expected_hash

        resolved[dap_arm] = arm_resolved

    logger.info("CRD_DRP manifest validated: %s", path)
    logger.info(
        "CRD_DRP arm mapping: BLUE -> BL; RED -> RH3 | SHA256 verification=%s",
        "enabled" if verify_hashes else "SKIPPED BY USER",
    )
    return path, manifest, resolved


def _load_crd_drp_atmospheric_mask(
    mask_path: Path,
    cube: io.KCWICube,
    *,
    expected_science_medium: str,
) -> tuple[np.ndarray, dict[str, object]]:
    """Load a native-grid CRD_DRP atmospheric mask and verify exact alignment."""
    with fits.open(mask_path, memmap=True) as hdul:
        if hdul[0].data is None:
            raise ValueError(f"CRD_DRP atmospheric mask contains no primary data: {mask_path}")
        raw = np.asarray(hdul[0].data)
        header = hdul[0].header.copy()

    mask = np.asarray(raw).squeeze()
    if mask.ndim != 1:
        raise ValueError(
            f"CRD_DRP atmospheric mask must be 1-D; got shape {raw.shape} in {mask_path}"
        )
    if mask.size != cube.nwave:
        raise ValueError(
            f"CRD_DRP atmospheric-mask length {mask.size} does not match "
            f"{cube.arm} science wavelength length {cube.nwave}."
        )
    finite = np.isfinite(mask)
    if not np.all(finite):
        raise ValueError(f"CRD_DRP atmospheric mask contains non-finite samples: {mask_path}")
    unique = set(np.unique(mask).tolist())
    if not unique.issubset({0, 1, False, True}):
        raise ValueError(
            f"CRD_DRP atmospheric mask is not binary in {mask_path}; values={sorted(unique)[:12]}"
        )
    mask = mask.astype(bool)

    mask_wave = io.wavelength_axis_from_header(header, mask.size, fits_axis=1)
    if not np.allclose(mask_wave, cube.wavelength, rtol=0.0, atol=1.0e-6):
        diff = float(np.nanmax(np.abs(mask_wave - cube.wavelength)))
        raise ValueError(
            f"CRD_DRP atmospheric-mask wavelength grid does not match {cube.arm} "
            f"science grid; max |delta lambda|={diff:.6g} A."
        )

    medium = str(header.get("WAVEMED", "")).strip().lower()
    if medium not in {"air", "vacuum"}:
        raise ValueError(
            f"CRD_DRP atmospheric mask {mask_path} does not declare WAVEMED=air/vacuum."
        )
    if medium != str(expected_science_medium).lower():
        raise ValueError(
            f"CRD_DRP atmospheric mask is observed {medium}, but {cube.arm} science "
            f"data are observed {expected_science_medium}. Regenerate the CRD_DRP "
            "package rather than silently converting the finalized mask."
        )

    return mask, {
        "source": str(mask_path),
        "wavelength_medium": medium,
        "n_native_pixels": int(mask.size),
        "n_masked_pixels": int(np.count_nonzero(mask)),
        "masked_fraction": float(np.mean(mask)),
    }


def _apply_crd_drp_atmospheric_mask(
    cube: io.KCWICube,
    atmospheric_mask: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    """Union the CRD_DRP wavelength mask with Script-1's native hard-good mask.

    The original sample-level/native quality state is intentionally retained in
    ``base_good``, ``good_spaxel``, and ``good_fraction_spaxel``.  Only the
    authoritative downstream wavelength/sample masks are changed:

        combined GOODWAVE = native GOODWAVE & ~ATMMASK
        combined GOODMASK = native GOODMASK & ~ATMMASK

    A global atmospheric exclusion therefore cannot by itself invalidate a
    spatial spaxel.
    """
    atm = np.asarray(atmospheric_mask, dtype=bool)
    if atm.ndim != 1 or atm.size != cube.nwave:
        raise ValueError("Atmospheric mask does not match cube wavelength dimension.")

    native_good_wavelength = np.asarray(cube.good_wavelength, dtype=bool).copy()
    native_good_count = int(np.count_nonzero(native_good_wavelength))
    atmospheric_on_native_good = atm & native_good_wavelength

    cube.good_wavelength = native_good_wavelength & ~atm
    cube.good = np.asarray(cube.good, dtype=bool) & (~atm)[None, None, :]

    info = {
        "n_native_good_wavelength_pixels_before_atmosphere": native_good_count,
        "n_atmospheric_mask_pixels": int(np.count_nonzero(atm)),
        "n_atmospheric_mask_pixels_overlapping_native_good": int(
            np.count_nonzero(atmospheric_on_native_good)
        ),
        "n_combined_good_wavelength_pixels": int(np.count_nonzero(cube.good_wavelength)),
        "combined_good_wavelength_fraction": float(np.mean(cube.good_wavelength)),
    }
    return native_good_wavelength, info


def _append_mask_provenance_to_prepared_cube(
    prepared_path: Path,
    *,
    native_good_wavelength: np.ndarray,
    atmospheric_mask: np.ndarray,
    reduction_manifest_path: Path,
) -> None:
    """Embed native-vs-atmospheric wavelength provenance in a prepared cube."""
    with fits.open(prepared_path, mode="update", memmap=False) as hdul:
        for name in ("NATIVEGW", "ATMMASK"):
            if name in hdul:
                raise RuntimeError(
                    f"Prepared cube already contains {name}; refusing ambiguous overwrite: {prepared_path}"
                )
        hdul.append(
            fits.ImageHDU(
                np.asarray(native_good_wavelength, dtype=np.uint8),
                name="NATIVEGW",
            )
        )
        atm_hdu = fits.ImageHDU(
            np.asarray(atmospheric_mask, dtype=np.uint8),
            name="ATMMASK",
        )
        atm_hdu.header["MASKTRUE"] = (1, "1 means exclude from stellar analysis")
        hdul.append(atm_hdu)
        hdul[0].header["CRDRPMF"] = (
            reduction_manifest_path.name[:68],
            "Source CRD_DRP reduction manifest basename",
        )
        hdul[0].header["CRDATMSK"] = (True, "CRD_DRP atmospheric mask applied")
        hdul.flush()

def _load_science_arm(arm: str, cfg, logger, reduction_input: dict[str, object] | None = None) -> io.KCWICube:
    """Load one science arm using the configured input layout.

    KcwiKit is the production path because it preserves the stacked flux,
    variance, mask, and effective-exposure products separately.  The legacy DRP
    path remains useful for single-cube validation and backwards compatibility.
    """
    prefix = arm.upper()

    if reduction_input is not None:
        cube = io.load_kcwikit_stack(
            reduction_input["icube"],
            reduction_input["vcube"],
            reduction_input["mcube"],
            reduction_input["ecube"],
            arm=arm,
            min_good_wavelength_fraction=float(cfg.MIN_GOOD_WAVELENGTH_FRACTION),
            bad_channel_fraction_threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
            float_dtype=str(cfg.STACK_FLOAT_DTYPE),
        )
        logger.info(
            "%s science input: validated CRD_DRP manifest arm=%s",
            arm,
            reduction_input.get("manifest_arm"),
        )
        for role, path in (cube.source_paths or {}).items():
            logger.info("%s %s: %s", arm, role, path)
        logger.info("%s KcwiKit stack-mask counts: %s", arm, io.summarize_binary_mask(cube.drp_mask))
        logger.info("%s effective-exposure summary: %s", arm, io.summarize_effective_exposure(cube.exposure))
        return cube

    input_format = str(cfg.SCIENCE_INPUT_FORMAT).strip().lower()

    if input_format == "kcwikit":
        cube = io.load_kcwikit_stack(
            getattr(cfg, f"{prefix}_ICUBE"),
            getattr(cfg, f"{prefix}_VCUBE"),
            getattr(cfg, f"{prefix}_MCUBE"),
            getattr(cfg, f"{prefix}_ECUBE"),
            arm=arm,
            min_good_wavelength_fraction=float(cfg.MIN_GOOD_WAVELENGTH_FRACTION),
            bad_channel_fraction_threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
            float_dtype=str(cfg.STACK_FLOAT_DTYPE),
        )
        logger.info("%s science input: KcwiKit four-file stack", arm)
        for role, path in (cube.source_paths or {}).items():
            logger.info("%s %s: %s", arm, role, path)
        logger.info("%s KcwiKit stack-mask counts: %s", arm, io.summarize_binary_mask(cube.drp_mask))
        logger.info("%s effective-exposure summary: %s", arm, io.summarize_effective_exposure(cube.exposure))
        logger.info(
            "%s note: original PyDRP FLAGS were consumed by KcwiKit during stacking; "
            "the final mcube is binary and no bit-level FLAGS are reconstructed.",
            arm,
        )
        return cube

    if input_format == "drp":
        cube = io.load_kcwi_cube(
            getattr(cfg, f"{prefix}_CUBE"),
            arm=arm,
            reject_any_nonzero_flag=bool(cfg.REJECT_ANY_NONZERO_DRP_FLAG),
            min_good_wavelength_fraction=float(cfg.MIN_GOOD_WAVELENGTH_FRACTION),
            bad_channel_fraction_threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
        )
        logger.info("%s science input: native/legacy DRP multi-extension cube", arm)
        logger.info("%s exact DRP FLAGS value counts: %s", arm, io.summarize_integer_flags(cube.flags))
        return cube

    raise ValueError(f"Unsupported SCIENCE_INPUT_FORMAT={cfg.SCIENCE_INPUT_FORMAT!r}")


def _save_exposure_diagnostic(arm: str, cube: io.KCWICube, run, logger) -> None:
    """Save a 2-D effective-exposure / coverage diagnostic for KcwiKit stacks."""
    if cube.exposure is None:
        logger.info("%s has no explicit exposure cube; skipping effective-exposure plot", arm)
        return

    # Keep zeros in the statistic so a spatial pixel covered at only a subset of
    # wavelengths is visibly distinguished from one with full wavelength
    # coverage.  Restrict the collapse to instrument/global-good wavelengths.
    exp = np.asarray(cube.exposure, dtype=float)
    use_wave = np.asarray(cube.good_wavelength, dtype=bool)
    if not np.any(use_wave):
        logger.warning("%s has no good wavelengths for exposure diagnostic", arm)
        return
    median_exp = np.nanmedian(exp[..., use_wave], axis=-1)
    coverage = (
        cube.coverage_fraction_spaxel
        if cube.coverage_fraction_spaxel is not None
        else np.mean(exp[..., use_wave] > 0, axis=-1)
    )
    plotting.plot_effective_exposure(
        median_exp,
        coverage,
        run.figures_dir / f"{arm.upper()}_effective_exposure.png",
        title=f"{arm.upper()} KcwiKit effective exposure / wavelength coverage",
    )


def _registration_images(
    bl: io.KCWICube,
    rh3: io.KCWICube,
    bl_default_image: np.ndarray,
    rh3_default_image: np.ndarray,
    cfg,
    logger,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Choose scientifically comparable 2-D images for registration QC.

    When the two arms share a sufficiently wide instrument-good wavelength
    interval, collapse *that same observed-frame interval* in both cubes.  This
    minimizes color/morphology differences between passbands.  If no useful
    overlap exists (the normal BL+RH3 production case may have none), retain the
    full-arm continuum images and rely on the explicit morphology-contrast test
    before trusting cross-correlation.
    """
    use_common = bool(
        getattr(cfg, "REGISTRATION_USE_COMMON_WAVELENGTH_IF_AVAILABLE", True)
    )
    min_width = float(getattr(cfg, "REGISTRATION_MIN_COMMON_RANGE_ANGSTROM", 50.0))
    min_channels = int(getattr(cfg, "REGISTRATION_MIN_COMMON_CHANNELS", 20))

    bl_good_wave = np.asarray(bl.good_wavelength, dtype=bool)
    rh3_good_wave = np.asarray(rh3.good_wavelength, dtype=bool)
    if not np.any(bl_good_wave) or not np.any(rh3_good_wave):
        return bl_default_image, rh3_default_image, {
            "mode": "full_arm_continuum",
            "reason": "one or both arms have no instrument-good wavelengths",
            "wavelength_min_A": None,
            "wavelength_max_A": None,
        }

    overlap_min = max(
        float(np.nanmin(bl.wavelength[bl_good_wave])),
        float(np.nanmin(rh3.wavelength[rh3_good_wave])),
    )
    overlap_max = min(
        float(np.nanmax(bl.wavelength[bl_good_wave])),
        float(np.nanmax(rh3.wavelength[rh3_good_wave])),
    )
    width = overlap_max - overlap_min
    bl_sel = bl_good_wave & (bl.wavelength >= overlap_min) & (bl.wavelength <= overlap_max)
    rh3_sel = rh3_good_wave & (rh3.wavelength >= overlap_min) & (rh3.wavelength <= overlap_max)

    enough = (
        use_common
        and np.isfinite(width)
        and width >= min_width
        and int(np.sum(bl_sel)) >= min_channels
        and int(np.sum(rh3_sel)) >= min_channels
    )
    if enough:
        bl_image = cube_utils.collapsed_continuum(
            bl.flux[..., bl_sel], bl.good[..., bl_sel], statistic="median"
        )
        rh3_image = cube_utils.collapsed_continuum(
            rh3.flux[..., rh3_sel], rh3.good[..., rh3_sel], statistic="median"
        )
        logger.info(
            "Registration will use common instrument-good wavelength interval "
            "%.2f--%.2f A (BL channels=%d, RH3 channels=%d)",
            overlap_min,
            overlap_max,
            int(np.sum(bl_sel)),
            int(np.sum(rh3_sel)),
        )
        return bl_image, rh3_image, {
            "mode": "common_wavelength_continuum",
            "reason": "sufficient shared instrument-good wavelength coverage",
            "wavelength_min_A": overlap_min,
            "wavelength_max_A": overlap_max,
            "BL_channels": int(np.sum(bl_sel)),
            "RH3_channels": int(np.sum(rh3_sel)),
        }

    reason = (
        f"common interval width={width:.2f} A, BL channels={int(np.sum(bl_sel))}, "
        f"RH3 channels={int(np.sum(rh3_sel))}; requirements are width>={min_width:.2f} A "
        f"and >= {min_channels} channels per arm"
    )
    logger.info("Registration will use full-arm continuum images: %s", reason)
    return bl_default_image, rh3_default_image, {
        "mode": "full_arm_continuum",
        "reason": reason,
        "wavelength_min_A": overlap_min if overlap_max > overlap_min else None,
        "wavelength_max_A": overlap_max if overlap_max > overlap_min else None,
        "BL_channels": int(np.sum(bl_sel)),
        "RH3_channels": int(np.sum(rh3_sel)),
    }


def _run_lsf(
    arm: str,
    cube: io.KCWICube,
    cfg,
    run,
    logger,
) -> tuple[psf_lsf.ArcLSFResult, dict[str, Path], validation.CalibrationMatch]:
    prefix = arm.upper()
    master_arc = Path(getattr(cfg, f"{prefix}_MASTER_ARC")).expanduser().resolve()

    # The pipeline arm labels remain BL/RH3 because the downstream science
    # architecture is organized around those two data streams.  The *actual*
    # grating expected in each stream is configurable so Script 1 can also be
    # integration-tested on other KCWI/KCRM setups (for example the existing
    # RL red-side test data) without disabling any calibration safety checks.
    canonical_grating = "BL" if prefix == "BL" else "RH3"
    expected_grating = str(
        getattr(cfg, f"{prefix}_EXPECTED_GRATING", canonical_grating)
    ).strip().upper()

    logger.info(
        "Measuring %s LSF from master arc: %s | configured expected grating=%s",
        arm,
        master_arc,
        expected_grating,
    )

    if expected_grating != canonical_grating:
        logger.warning(
            "NONSTANDARD_%s_GRATING | %s pipeline arm is configured for grating %s "
            "rather than canonical %s. Script 1 will validate and characterize "
            "this setup normally, but do not interpret later %s-specific science "
            "stages as a production %s analysis unless they are explicitly adapted.",
            prefix,
            arm,
            expected_grating,
            canonical_grating,
            canonical_grating,
            canonical_grating,
        )

    # Before measuring any widths, verify that the calibration actually belongs
    # to the same instrumental setup as the science cube. Configurability changes
    # only the requested grating name: science and arc still have to match that
    # grating, the same camera, IFU, binning, and central wavelength.
    arc_header = psf_lsf.read_primary_header(master_arc)
    calibration_match = validation.validate_arc_science_configuration(
        cube.header,
        arc_header,
        arm=arm,
        expected_grating=expected_grating,
        central_wavelength_tolerance_angstrom=float(
            getattr(cfg, "ARC_SCIENCE_CWAVE_TOLERANCE_ANGSTROM", 0.5)
        ),
    )
    logger.info(
        "%s science/master-arc configuration verified: camera=%s, expected grating=%s, "
        "science grating=%s, arc grating=%s, IFU=%s, binning=%s, central wavelength=%s A",
        arm,
        calibration_match.arc_camera,
        calibration_match.expected_grating,
        calibration_match.science_grating,
        calibration_match.arc_grating,
        calibration_match.arc_ifu,
        calibration_match.arc_binning,
        calibration_match.arc_cwave_angstrom,
    )

    result, sidecars = psf_lsf.measure_arc_lsf_from_files(
        master_arc,
        wavemap=_config_path(cfg, f"{prefix}_ARC_WAVEMAP"),
        slicemap=_config_path(cfg, f"{prefix}_ARC_SLICEMAP"),
        posmap=_config_path(cfg, f"{prefix}_ARC_POSMAP"),
        polynomial_order=int(cfg.LSF_MODEL_WAVELENGTH_ORDER),
        measure_spatial_variation=bool(cfg.LSF_MEASURE_SPATIAL_VARIATION),
        spatial_bins=int(cfg.LSF_SPATIAL_BINS),
        peak_prominence_fraction=float(cfg.ARC_PEAK_PROMINENCE_FRACTION),
        peak_height_percentile=float(cfg.ARC_PEAK_HEIGHT_PERCENTILE),
        min_peak_distance_pix=int(cfg.ARC_MIN_PEAK_DISTANCE_PIX),
        line_fit_half_width_pix=int(cfg.ARC_LINE_FIT_HALF_WIDTH_PIX),
        min_good_lines=int(cfg.ARC_MIN_GOOD_LINES),
        sigma_clip=float(cfg.ARC_LSF_SIGMA_CLIP),
    )

    npz = run.products_dir / f"{prefix}_lsf.npz"
    csv = run.products_dir / f"{prefix}_lsf_measurements.csv"
    psf_lsf.save_arc_lsf_result(result, npz, csv)

    plotting.plot_lsf(
        result,
        run.figures_dir / f"{prefix}_LSF.png",
        title=f"{prefix} empirical master-arc LSF",
    )
    plotting.plot_lsf_spatial_variation(
        result,
        run.figures_dir / f"{prefix}_LSF_spatial_variation.png",
        title=f"{prefix} LSF spatial/slice variation: individual line fits",
    )
    plotting.plot_lsf_spatial_summary(
        result,
        run.figures_dir / f"{prefix}_LSF_spatial_summary.png",
        title=f"{prefix} coherent LSF spatial summary",
    )

    logger.info(
        "%s LSF: %d accepted line measurements (%d initially fit); "
        "instrument-good/model domain %.1f--%.1f A; empirical accepted-line "
        "support %.1f--%.1f A; median FWHM %.4f A",
        arm,
        result.n_lines_used,
        result.n_lines_total,
        result.wavelength_min,
        result.wavelength_max,
        result.measurement_wavelength_min,
        result.measurement_wavelength_max,
        float(np.nanmedian(result.fwhm_angstrom)),
    )
    logger.info(
        "%s LSF unconstrained instrument-good edge widths: blue=%.1f A, red=%.1f A. "
        "The saved LSF evaluator returns NaN outside the empirical accepted-line interval by default.",
        arm,
        result.blue_edge_unconstrained_angstrom,
        result.red_edge_unconstrained_angstrom,
    )
    logger.info(
        "%s LSF fractional scatter: raw line-to-line RMS=%.4f; "
        "slice RMS=%.4f; position-bin RMS=%.4f; slice+position group RMS=%.4f",
        arm,
        result.measurement_fractional_rms,
        result.slice_fractional_rms,
        result.position_fractional_rms,
        result.spatial_fractional_rms,
    )
    for name, path in sidecars.items():
        logger.info("%s master-arc sidecar %s: %s", arm, name, path)
    return result, sidecars, calibration_match


def _run_noise_diagnostic(arm: str, cube: io.KCWICube, image: np.ndarray, cfg, run, logger):
    prefix = arm.upper()
    logger.info("Running preliminary %s noise/covariance diagnostic", arm)
    result = noise.characterize_preliminary_noise(
        cube.flux,
        cube.uncertainty,
        cube.good,
        image,
        cube.good_spaxel,
        max_spaxels=int(cfg.NOISE_DIAGNOSTIC_MAX_SPAXELS),
        low_flux_percentile=float(cfg.NOISE_LOW_FLUX_PERCENTILE),
        savgol_window=int(cfg.NOISE_SAVGOL_WINDOW),
        savgol_polyorder=int(cfg.NOISE_SAVGOL_POLYORDER),
        max_lag=int(cfg.NOISE_MAX_SPECTRAL_LAG),
    )
    noise.save_noise_diagnostic(result, str(run.products_dir / f"{prefix}_noise_diagnostic.npz"))
    plotting.plot_normalized_residuals(
        result,
        run.figures_dir / f"noise_normalized_residuals_{prefix}.png",
        title=f"{prefix} preliminary normalized residuals",
    )
    plotting.plot_spectral_covariance(
        result,
        run.figures_dir / f"spectral_covariance_{prefix}.png",
        title=f"{prefix} preliminary spectral correlation",
    )
    logger.info(
        "%s preliminary noise diagnostic: variance scale=%.4f, spaxels=%d, samples=%d",
        arm,
        result.variance_scale_factor,
        result.n_spaxels_used,
        result.n_samples_used,
    )
    return result


def main() -> int:
    args = _parser().parse_args()

    # When a CRD_DRP manifest is supplied, science-cube paths come from that
    # validated handoff instead of the target config.  Keep structural config
    # validation, but validate only the non-science external paths here.
    using_reduction_manifest = args.reduction_manifest is not None
    cfg = crd.load_config(
        args.config,
        validate=True,
        strict_paths=not using_reduction_manifest,
    )
    if using_reduction_manifest:
        crd.validate_input_paths(
            cfg,
            (
                "BL_MASTER_ARC",
                "RH3_MASTER_ARC",
                "PYMORPH_VAC",
                "XSL_TEMPLATE_LIBRARY",
            ),
        )

    run = crd.create_run_context(cfg, run_name=args.run_name)
    crd.setup_pipeline_logger(run)
    logger = crd.setup_step_logger(run, "01_prepare_and_register_cubes")
    quality_flags: list[str] = []
    start = time.perf_counter()

    reduction_manifest_path: Path | None = None
    reduction_manifest: dict | None = None
    reduction_inputs: dict[str, dict[str, object]] | None = None
    atmospheric_masks: dict[str, np.ndarray] = {}
    native_good_wavelengths: dict[str, np.ndarray] = {}
    atmospheric_provenance: dict[str, dict[str, object]] = {}

    try:
        crd.log_section(logger, "CRD_DAP SCRIPT 1: PREPARE AND REGISTER CUBES", "=")
        logger.info("Run directory: %s", run.run_dir)
        logger.info("Python: %s", sys.version.replace("\n", " "))
        logger.info("Platform: %s", platform.platform())
        logger.info("Target: %s", cfg.TARGET_NAME)
        logger.info("Redshift: %.8f", float(cfg.REDSHIFT))

        if using_reduction_manifest:
            (
                reduction_manifest_path,
                reduction_manifest,
                reduction_inputs,
            ) = _load_reduction_manifest(
                args.reduction_manifest,
                verify_hashes=not bool(args.skip_reduction_hash_check),
                logger=logger,
            )
        else:
            logger.warning(
                "LEGACY_SCIENCE_INPUT_PATH | No --reduction-manifest supplied. "
                "Script 1 will use science-cube paths from the target config and "
                "no CRD_DRP atmospheric masks will be applied."
            )

        # ------------------------------------------------------------------
        # 1. Read science cubes and establish hard data-quality masks.
        # ------------------------------------------------------------------
        crd.log_section(logger, "1. Load KCWI/KCRM science cubes")
        bl = _load_science_arm(
            "BL",
            cfg,
            logger,
            None if reduction_inputs is None else reduction_inputs["BL"],
        )
        rh3 = _load_science_arm(
            "RH3",
            cfg,
            logger,
            None if reduction_inputs is None else reduction_inputs["RH3"],
        )

        for cube in (bl, rh3):
            logger.info(
                "%s cube %s: standardized shape=%s, wavelength=%.2f--%.2f A, "
                "usable spaxels=%d/%d, usable channels=%d/%d",
                cube.arm,
                cube.path,
                cube.shape,
                float(np.nanmin(cube.wavelength)),
                float(np.nanmax(cube.wavelength)),
                int(np.sum(cube.good_spaxel)),
                int(cube.good_spaxel.size),
                int(np.sum(cube.good_wavelength)),
                int(cube.good_wavelength.size),
            )
            logger.info("%s primary science FITS structure: %s", cube.arm, io.inspect_fits_extensions(cube.path))
            logger.info(
                "%s spatial wavelength-coverage fraction: median=%.3f, min=%.3f, max=%.3f",
                cube.arm,
                float(np.nanmedian(cube.coverage_fraction_spaxel)),
                float(np.nanmin(cube.coverage_fraction_spaxel)),
                float(np.nanmax(cube.coverage_fraction_spaxel)),
            )

        # Wavelength medium/frame bookkeeping is explicit. The XSL template
        # library may currently be in a different medium; that is recorded as a
        # required conversion before pPXF rather than silently ignored.
        for cube in (bl, rh3):
            conventions = validation.validate_script1_conventions(
                cube.header,
                science_medium=cfg.SCIENCE_WAVELENGTH_MEDIUM,
                template_medium=cfg.TEMPLATE_WAVELENGTH_MEDIUM,
                science_velocity_frame=cfg.SCIENCE_VELOCITY_FRAME,
            )
            logger.info(
                "%s conventions: science=%s wavelengths (header=%s), velocity frame=%s "
                "(header=%s), template=%s wavelengths, template conversion required=%s",
                cube.arm,
                conventions.science_medium,
                conventions.header_wavelength_medium,
                conventions.science_velocity_frame,
                conventions.header_velocity_frame,
                conventions.template_medium,
                conventions.template_conversion_required,
            )


        # ------------------------------------------------------------------
        # 1b. Apply the validated CRD_DRP atmospheric wavelength masks.
        # ------------------------------------------------------------------
        if reduction_inputs is not None:
            crd.log_section(logger, "1b. Apply CRD_DRP atmospheric masks")
            for cube in (bl, rh3):
                prefix = cube.arm.upper()
                conventions = validation.validate_script1_conventions(
                    cube.header,
                    science_medium=cfg.SCIENCE_WAVELENGTH_MEDIUM,
                    template_medium=cfg.TEMPLATE_WAVELENGTH_MEDIUM,
                    science_velocity_frame=cfg.SCIENCE_VELOCITY_FRAME,
                )
                mask_path = Path(
                    reduction_inputs[prefix]["atmospheric_mask"]
                ).expanduser().resolve()
                atm_mask, mask_info = _load_crd_drp_atmospheric_mask(
                    mask_path,
                    cube,
                    expected_science_medium=conventions.science_medium,
                )
                native_good_wave, apply_info = _apply_crd_drp_atmospheric_mask(
                    cube,
                    atm_mask,
                )
                atmospheric_masks[prefix] = atm_mask
                native_good_wavelengths[prefix] = native_good_wave

                provenance = {
                    **mask_info,
                    **apply_info,
                    "manifest_arm": reduction_inputs[prefix]["manifest_arm"],
                    "mask_sha256": reduction_inputs[prefix].get(
                        "atmospheric_mask_sha256"
                    ),
                    "interval_table": (
                        None
                        if reduction_inputs[prefix].get("atmospheric_intervals") is None
                        else str(reduction_inputs[prefix]["atmospheric_intervals"])
                    ),
                    "semantics": "True/1 = exclude from stellar analysis",
                    "applied_in_script": "01_prepare_and_register_cubes.py",
                }
                atmospheric_provenance[prefix] = provenance

                logger.info(
                    "%s CRD_DRP atmospheric mask: %d/%d native channels masked "
                    "(%.2f%%); overlap with native Script-1 GOODWAVE=%d; "
                    "combined GOODWAVE=%d/%d",
                    prefix,
                    provenance["n_masked_pixels"],
                    provenance["n_native_pixels"],
                    100.0 * provenance["masked_fraction"],
                    provenance["n_atmospheric_mask_pixels_overlapping_native_good"],
                    provenance["n_combined_good_wavelength_pixels"],
                    provenance["n_native_pixels"],
                )

        # ------------------------------------------------------------------
        # 2. Continuum images and center diagnostics.
        # ------------------------------------------------------------------
        crd.log_section(logger, "2. Collapsed continuum and center diagnostics")
        bl_image = cube_utils.collapsed_continuum(bl.flux, bl.good, statistic="median")
        rh3_image = cube_utils.collapsed_continuum(rh3.flux, rh3.good, statistic="median")

        bl_center = cube_utils.estimate_continuum_center(
            bl_image,
            bl.celestial_wcs,
            smooth_sigma_pix=float(cfg.CENTER_SMOOTH_SIGMA_PIX),
            centroid_min_percentile=float(cfg.CENTER_CENTROID_MIN_PERCENTILE),
        )
        rh3_center = cube_utils.estimate_continuum_center(
            rh3_image,
            rh3.celestial_wcs,
            smooth_sigma_pix=float(cfg.CENTER_SMOOTH_SIGMA_PIX),
            centroid_min_percentile=float(cfg.CENTER_CENTROID_MIN_PERCENTILE),
        )

        logger.info("BL peak pixel (y,x)=%s; centroid=%s", bl_center.peak_yx, bl_center.centroid_yx)
        logger.info("RH3 peak pixel (y,x)=%s; centroid=%s", rh3_center.peak_yx, rh3_center.centroid_yx)
        peak_sep = validation.sky_separation_arcsec(bl_center.peak_sky, rh3_center.peak_sky)
        logger.info("BL/RH3 continuum-peak sky separation: %.4f arcsec", peak_sep)
        if peak_sep > float(cfg.CENTER_WARNING_ARCSEC):
            _quality_flag(
                quality_flags,
                "CENTER_DISAGREEMENT",
                logger,
                f"BL/RH3 continuum peaks differ by {peak_sep:.3f} arcsec.",
            )

        plotting.plot_collapsed_continuum(
            bl_image,
            run.figures_dir / "BL_collapsed_continuum.png",
            title="BL collapsed continuum",
            peak_yx=bl_center.peak_yx,
            centroid_yx=bl_center.centroid_yx,
        )
        plotting.plot_collapsed_continuum(
            rh3_image,
            run.figures_dir / "RH3_collapsed_continuum.png",
            title="RH3 collapsed continuum",
            peak_yx=rh3_center.peak_yx,
            centroid_yx=rh3_center.centroid_yx,
        )
        plotting.plot_center_comparison(
            bl_image,
            bl.celestial_wcs,
            bl_peak_sky=bl_center.peak_sky,
            bl_centroid_sky=bl_center.centroid_sky,
            rh3_peak_sky=rh3_center.peak_sky,
            rh3_centroid_sky=rh3_center.centroid_sky,
            path=run.figures_dir / "geometry_center_comparison.png",
        )

        # ------------------------------------------------------------------
        # 3. Registration diagnostic. The science cubes themselves are NOT
        #    resampled here; Script 2 will use the saved physical coordinates.
        # ------------------------------------------------------------------
        crd.log_section(logger, "3. BL/RH3 WCS registration")
        bl_registration_image, rh3_registration_image, registration_image_info = _registration_images(
            bl, rh3, bl_image, rh3_image, cfg, logger
        )
        bl_scale = bl.pixel_scales_arcsec()
        registration = cube_utils.register_cube_pair(
            bl_registration_image,
            bl.celestial_wcs,
            rh3_registration_image,
            rh3.celestial_wcs,
            reference_pixel_scale_arcsec=bl_scale,
            min_contrast_snr=float(getattr(cfg, "REGISTRATION_MIN_CONTRAST_SNR", 5.0)),
            contrast_smooth_sigma_pix=float(
                getattr(cfg, "REGISTRATION_CONTRAST_SMOOTH_SIGMA_PIX", 1.0)
            ),
            max_residual_shift_arcsec=float(
                getattr(cfg, "REGISTRATION_MAX_RESIDUAL_SHIFT_ARCSEC", 2.0)
            ),
        )
        logger.info(
            "Registration morphology contrast: BL=%.3f, RH3=%.3f; valid=%s | %s",
            registration.reference_contrast_snr,
            registration.moving_contrast_snr,
            registration.cross_correlation_valid,
            registration.status_reason,
        )
        if registration.cross_correlation_valid:
            logger.info(
                "Residual RH3->BL shift after WCS reprojection: dy=%.4f pix, dx=%.4f pix; "
                "dy=%.4f arcsec, dx=%.4f arcsec; radius=%.4f arcsec",
                registration.residual_shift_yx_pix[0],
                registration.residual_shift_yx_pix[1],
                registration.residual_shift_arcsec[0],
                registration.residual_shift_arcsec[1],
                registration.residual_shift_radius_arcsec,
            )
            if registration.residual_shift_radius_arcsec > float(cfg.REGISTRATION_WARNING_ARCSEC):
                _quality_flag(
                    quality_flags,
                    "REGISTRATION_OFFSET_WARNING",
                    logger,
                    "Residual morphology shift exceeds configured registration tolerance; "
                    "inspect BL_RH3_registration.png before Script 2.",
                )
        else:
            _quality_flag(
                quality_flags,
                "REGISTRATION_INCONCLUSIVE",
                logger,
                "Morphology cross-correlation was not trusted. This can occur because the "
                "registration images lack sufficient contrast or because the best correlation "
                "peak runs to the edge of the allowed local residual-search window. Use the "
                "independent center/WCS comparison and inspect BL_RH3_registration.png; no "
                "numerical residual shift is adopted.",
            )

        if registration_image_info.get("mode") == "common_wavelength_continuum":
            wave_label = (
                f"common observed wavelength "
                f"{registration_image_info['wavelength_min_A']:.1f}--"
                f"{registration_image_info['wavelength_max_A']:.1f} A"
            )
        else:
            wave_label = "full-arm continuum images"
        plotting.plot_registration(
            bl_registration_image,
            registration.moving_on_reference,
            registration.difference,
            run.figures_dir / "BL_RH3_registration.png",
            overlap=registration.overlap,
            residual_shift_arcsec=registration.residual_shift_arcsec,
            cross_correlation_valid=registration.cross_correlation_valid,
            status_reason=registration.status_reason,
            wavelength_label=wave_label,
            reference_contrast_snr=registration.reference_contrast_snr,
            moving_contrast_snr=registration.moving_contrast_snr,
        )

        # The common coordinate origin is bookkeeping only; Script 4 later fits
        # the kinematic center. BL_peak is the default because BL defines bins.
        center_source = str(cfg.COMMON_CENTER_SOURCE).strip().lower()
        if center_source == "bl_peak":
            common_center = bl_center.peak_sky
        elif center_source == "rh3_peak":
            common_center = rh3_center.peak_sky
        else:
            raise ValueError(
                "COMMON_CENTER_SOURCE must currently be 'BL_peak' or 'RH3_peak'"
            )

        _save_spatial_coordinates(
            bl,
            common_center,
            run.products_dir / "BL_spatial_coordinates.npz",
        )
        _save_spatial_coordinates(
            rh3,
            common_center,
            run.products_dir / "RH3_spatial_coordinates.npz",
        )

        io.write_json(
            {
                "common_center_source": cfg.COMMON_CENTER_SOURCE,
                "common_center_ra_deg": float(common_center.ra.deg),
                "common_center_dec_deg": float(common_center.dec.deg),
                "BL_peak_yx": bl_center.peak_yx,
                "BL_centroid_yx": bl_center.centroid_yx,
                "BL_peak_ra_deg": float(bl_center.peak_sky.ra.deg),
                "BL_peak_dec_deg": float(bl_center.peak_sky.dec.deg),
                "RH3_peak_yx": rh3_center.peak_yx,
                "RH3_centroid_yx": rh3_center.centroid_yx,
                "RH3_peak_ra_deg": float(rh3_center.peak_sky.ra.deg),
                "RH3_peak_dec_deg": float(rh3_center.peak_sky.dec.deg),
                "peak_separation_arcsec": peak_sep,
            },
            run.metadata_dir / "centers.json",
        )
        io.write_json(
            {
                "residual_shift_dy_pix": registration.residual_shift_yx_pix[0],
                "residual_shift_dx_pix": registration.residual_shift_yx_pix[1],
                "residual_shift_dy_arcsec": registration.residual_shift_arcsec[0],
                "residual_shift_dx_arcsec": registration.residual_shift_arcsec[1],
                "residual_shift_radius_arcsec": registration.residual_shift_radius_arcsec,
                "cross_correlation_valid": registration.cross_correlation_valid,
                "status_reason": registration.status_reason,
                "BL_registration_contrast_snr": registration.reference_contrast_snr,
                "RH3_registration_contrast_snr": registration.moving_contrast_snr,
                "registration_image": registration_image_info,
                "overlap_fraction_of_BL_grid": float(np.mean(registration.overlap)),
                "science_cubes_resampled": False,
            },
            run.metadata_dir / "registration.json",
        )

        # ------------------------------------------------------------------
        # 4. Data-quality diagnostics.
        # ------------------------------------------------------------------
        crd.log_section(logger, "4. Data-quality diagnostics")
        _save_exposure_diagnostic("BL", bl, run, logger)
        _save_exposure_diagnostic("RH3", rh3, run, logger)
        for cube in (bl, rh3):
            prefix = cube.arm.upper()
            plotting.plot_valid_spaxels(
                cube.good_spaxel,
                cube.good_fraction_spaxel,
                run.figures_dir / f"{prefix}_valid_spaxels.png",
                title=f"{prefix} valid spatial samples",
            )
            plotting.plot_bad_wavelength_fraction(
                cube.wavelength,
                cube.bad_fraction_wavelength,
                run.figures_dir / f"{prefix}_bad_wavelength_fraction.png",
                title=f"{prefix} bad spatial-sample fraction by wavelength",
                threshold=float(cfg.BAD_CHANNEL_FRACTION_THRESHOLD),
                wavegood0=cube.header.get("WAVGOOD0"),
                wavegood1=cube.header.get("WAVGOOD1"),
            )

        # ------------------------------------------------------------------
        # 5. PSF provenance. Never estimate seeing from the extended galaxy.
        # ------------------------------------------------------------------
        crd.log_section(logger, "5. PSF characterization/provenance")
        bl_psf = psf_lsf.estimate_psf(
            configured_fwhm_arcsec=cfg.BL_PSF_FWHM_ARCSEC,
            header=bl.header,
            header_keys=tuple(cfg.PSF_HEADER_KEYS),
        )
        rh3_psf = psf_lsf.estimate_psf(
            configured_fwhm_arcsec=cfg.RH3_PSF_FWHM_ARCSEC,
            header=rh3.header,
            header_keys=tuple(cfg.PSF_HEADER_KEYS),
        )
        for arm, estimate in (("BL", bl_psf), ("RH3", rh3_psf)):
            logger.info(
                "%s PSF FWHM: %s arcsec | source=%s | %s",
                arm,
                "nan" if not np.isfinite(estimate.fwhm_arcsec) else f"{estimate.fwhm_arcsec:.4f}",
                estimate.source,
                estimate.detail,
            )
            if not np.isfinite(estimate.fwhm_arcsec):
                _quality_flag(
                    quality_flags,
                    "PSF_NOT_MEASURED",
                    logger,
                    f"{arm} PSF is not known yet. Supply a measured/configured value before PSF-dependent modeling.",
                )
        plotting.plot_psf_summary(bl_psf, run.figures_dir / "BL_PSF.png", title="BL PSF provenance")
        plotting.plot_psf_summary(rh3_psf, run.figures_dir / "RH3_PSF.png", title="RH3 PSF provenance")
        plotting.plot_psf_comparison(bl_psf, rh3_psf, run.figures_dir / "BL_RH3_PSF_comparison.png")
        io.write_json(
            {
                "BL": bl_psf.__dict__,
                "RH3": rh3_psf.__dict__,
            },
            run.metadata_dir / "psf_summary.json",
        )

        # ------------------------------------------------------------------
        # 6. Empirical master-arc LSFs.
        # ------------------------------------------------------------------
        crd.log_section(logger, "6. Empirical master-arc LSF")
        bl_lsf, bl_sidecars, bl_cal_match = _run_lsf("BL", bl, cfg, run, logger)
        rh3_lsf, rh3_sidecars, rh3_cal_match = _run_lsf("RH3", rh3, cfg, run, logger)
        for arm, result in (("BL", bl_lsf), ("RH3", rh3_lsf)):
            if result.spatial_fractional_rms > float(cfg.LSF_SPATIAL_VARIATION_WARNING_FRACTION):
                _quality_flag(
                    quality_flags,
                    "LSF_SPATIAL_VARIATION",
                    logger,
                    f"{arm} LSF fractional spatial RMS={result.spatial_fractional_rms:.3f} exceeds "
                    f"{float(cfg.LSF_SPATIAL_VARIATION_WARNING_FRACTION):.3f}.",
                )
            edge_limit = float(
                getattr(cfg, "LSF_EDGE_EXTRAPOLATION_WARNING_ANGSTROM", 100.0)
            )
            if max(
                result.blue_edge_unconstrained_angstrom,
                result.red_edge_unconstrained_angstrom,
            ) > edge_limit:
                _quality_flag(
                    quality_flags,
                    "LSF_EMPIRICAL_COVERAGE_GAP",
                    logger,
                    f"{arm} has instrument-good wavelength coverage outside the accepted arc-line "
                    f"support by {result.blue_edge_unconstrained_angstrom:.1f} A on the blue edge "
                    f"and {result.red_edge_unconstrained_angstrom:.1f} A on the red edge. "
                    "The LSF evaluator returns NaN in unsupported regions by default; later fitting "
                    "must mask them or supply an independently validated LSF model.",
                )

        io.write_json(
            {
                "BL": {
                    "master_arc": str(Path(cfg.BL_MASTER_ARC).expanduser().resolve()),
                    "sidecars": {k: str(v) for k, v in bl_sidecars.items()},
                    "n_lines_used": bl_lsf.n_lines_used,
                    "instrument_good_wavelength_min_A": bl_lsf.wavelength_min,
                    "instrument_good_wavelength_max_A": bl_lsf.wavelength_max,
                    "empirical_measurement_wavelength_min_A": bl_lsf.measurement_wavelength_min,
                    "empirical_measurement_wavelength_max_A": bl_lsf.measurement_wavelength_max,
                    "blue_edge_unconstrained_A": bl_lsf.blue_edge_unconstrained_angstrom,
                    "red_edge_unconstrained_A": bl_lsf.red_edge_unconstrained_angstrom,
                    "median_fwhm_A": float(np.nanmedian(bl_lsf.fwhm_angstrom)),
                    "measurement_fractional_rms": bl_lsf.measurement_fractional_rms,
                    "slice_fractional_rms": bl_lsf.slice_fractional_rms,
                    "position_fractional_rms": bl_lsf.position_fractional_rms,
                    "spatial_fractional_rms": bl_lsf.spatial_fractional_rms,
                    "calibration_match": bl_cal_match.__dict__,
                },
                "RH3": {
                    "master_arc": str(Path(cfg.RH3_MASTER_ARC).expanduser().resolve()),
                    "sidecars": {k: str(v) for k, v in rh3_sidecars.items()},
                    "n_lines_used": rh3_lsf.n_lines_used,
                    "instrument_good_wavelength_min_A": rh3_lsf.wavelength_min,
                    "instrument_good_wavelength_max_A": rh3_lsf.wavelength_max,
                    "empirical_measurement_wavelength_min_A": rh3_lsf.measurement_wavelength_min,
                    "empirical_measurement_wavelength_max_A": rh3_lsf.measurement_wavelength_max,
                    "blue_edge_unconstrained_A": rh3_lsf.blue_edge_unconstrained_angstrom,
                    "red_edge_unconstrained_A": rh3_lsf.red_edge_unconstrained_angstrom,
                    "median_fwhm_A": float(np.nanmedian(rh3_lsf.fwhm_angstrom)),
                    "measurement_fractional_rms": rh3_lsf.measurement_fractional_rms,
                    "slice_fractional_rms": rh3_lsf.slice_fractional_rms,
                    "position_fractional_rms": rh3_lsf.position_fractional_rms,
                    "spatial_fractional_rms": rh3_lsf.spatial_fractional_rms,
                    "calibration_match": rh3_cal_match.__dict__,
                },
            },
            run.metadata_dir / "lsf_summary.json",
        )

        # ------------------------------------------------------------------
        # 7. Preliminary noise/covariance diagnostics.
        # ------------------------------------------------------------------
        crd.log_section(logger, "7. Preliminary noise and covariance diagnostics")
        noise_results = {}
        for arm, cube, image in (("BL", bl, bl_image), ("RH3", rh3, rh3_image)):
            try:
                result = _run_noise_diagnostic(arm, cube, image, cfg, run, logger)
                noise_results[arm] = result
                if abs(result.variance_scale_factor - 1.0) > float(cfg.NOISE_VARIANCE_SCALE_WARNING_FRACTION):
                    _quality_flag(
                        quality_flags,
                        "NOISE_SCALE_WARNING",
                        logger,
                        f"{arm} preliminary variance scale={result.variance_scale_factor:.3f}; "
                        "do not automatically rescale from this high-pass diagnostic alone.",
                    )
                if result.correlation.size > 1:
                    max_corr = float(np.nanmax(np.abs(result.correlation[1:])))
                    if max_corr > float(cfg.NOISE_CORRELATION_WARNING_ABS):
                        _quality_flag(
                            quality_flags,
                            "SPECTRAL_COVARIANCE_WARNING",
                            logger,
                            f"{arm} preliminary |rho(lag>0)| reaches {max_corr:.3f}.",
                        )
            except Exception as exc:
                logger.warning(
                    "%s preliminary noise diagnostic could not be completed: %s. "
                    "This does not modify the science cube, but noise calibration must be revisited after pPXF.",
                    arm,
                    exc,
                )

        if bool(cfg.APPLY_PRELIMINARY_VARIANCE_RESCALING):
            raise ValueError(
                "APPLY_PRELIMINARY_VARIANCE_RESCALING=True is intentionally unsupported in the first "
                "production Script-1 implementation. The high-pass estimate is QC-only; perform final "
                "rescaling from validated stellar-fit residuals."
            )

        # ------------------------------------------------------------------
        # 8. Prepared products and manifest.
        # ------------------------------------------------------------------
        crd.log_section(logger, "8. Save prepared products")
        overwrite = bool(cfg.OVERWRITE_EXISTING_PRODUCTS)
        bl_prepared = io.save_prepared_cube(
            bl,
            run.products_dir / "prepared_BL.fits",
            overwrite=overwrite,
        )
        rh3_prepared = io.save_prepared_cube(
            rh3,
            run.products_dir / "prepared_RH3.fits",
            overwrite=overwrite,
        )
        if reduction_manifest_path is not None:
            _append_mask_provenance_to_prepared_cube(
                Path(bl_prepared),
                native_good_wavelength=native_good_wavelengths["BL"],
                atmospheric_mask=atmospheric_masks["BL"],
                reduction_manifest_path=reduction_manifest_path,
            )
            _append_mask_provenance_to_prepared_cube(
                Path(rh3_prepared),
                native_good_wavelength=native_good_wavelengths["RH3"],
                atmospheric_mask=atmospheric_masks["RH3"],
                reduction_manifest_path=reduction_manifest_path,
            )
            logger.info(
                "Prepared cubes include NATIVEGW and ATMMASK extensions; "
                "GOODWAVE/GOODMASK are the authoritative combined masks."
            )

        logger.info("Prepared BL cube: %s", bl_prepared)
        logger.info("Prepared RH3 cube: %s", rh3_prepared)

        elapsed = time.perf_counter() - start
        manifest = {
            "script": "01_prepare_and_register_cubes.py",
            "target": cfg.TARGET_NAME,
            "elapsed_seconds": elapsed,
            "prepared_BL": str(bl_prepared),
            "prepared_RH3": str(rh3_prepared),
            "science_input_format": (
                "crd_drp_manifest"
                if reduction_manifest_path is not None
                else str(cfg.SCIENCE_INPUT_FORMAT)
            ),
            "source_reduction_manifest": (
                None if reduction_manifest_path is None else str(reduction_manifest_path)
            ),
            "source_reduction_manifest_sha256": (
                None
                if reduction_manifest_path is None
                else _sha256_file(reduction_manifest_path)
            ),
            "source_reduction_validation_status": (
                None
                if reduction_manifest is None
                else reduction_manifest.get("validation_status")
            ),
            "source_reduction_schema": (
                None if reduction_manifest is None else reduction_manifest.get("schema")
            ),
            "reduction_manifest_arm_map": (
                None if reduction_manifest_path is None else {"BL": "blue", "RH3": "red"}
            ),
            "BL_atmospheric_mask": atmospheric_provenance.get("BL"),
            "RH3_atmospheric_mask": atmospheric_provenance.get("RH3"),
            "prepared_mask_contract": (
                None
                if reduction_manifest_path is None
                else {
                    "GOODWAVE": "native Script-1 wavelength quality AND NOT CRD_DRP ATMMASK",
                    "GOODMASK": "native Script-1 sample quality AND combined GOODWAVE",
                    "NATIVEGW": "Script-1 native wavelength quality before CRD_DRP atmosphere",
                    "ATMMASK": "CRD_DRP atmospheric mask; 1 means excluded",
                }
            ),
            "BL_source_paths": {k: str(v) for k, v in (bl.source_paths or {}).items()},
            "RH3_source_paths": {k: str(v) for k, v in (rh3.source_paths or {}).items()},
            "BL_effective_exposure": io.summarize_effective_exposure(bl.exposure),
            "RH3_effective_exposure": io.summarize_effective_exposure(rh3.exposure),
            "quality_flags": quality_flags,
            "common_center_ra_deg": float(common_center.ra.deg),
            "common_center_dec_deg": float(common_center.dec.deg),
            "registration_cross_correlation_valid": registration.cross_correlation_valid,
            "registration_residual_arcsec": registration.residual_shift_radius_arcsec,
            "registration_image_mode": registration_image_info.get("mode"),
            "BL_PSF_FWHM_arcsec": bl_psf.fwhm_arcsec,
            "RH3_PSF_FWHM_arcsec": rh3_psf.fwhm_arcsec,
            "BL_LSF_median_FWHM_A": float(np.nanmedian(bl_lsf.fwhm_angstrom)),
            "RH3_LSF_median_FWHM_A": float(np.nanmedian(rh3_lsf.fwhm_angstrom)),
            "BL_LSF_empirical_range_A": [
                bl_lsf.measurement_wavelength_min, bl_lsf.measurement_wavelength_max
            ],
            "RH3_LSF_empirical_range_A": [
                rh3_lsf.measurement_wavelength_min, rh3_lsf.measurement_wavelength_max
            ],
        }
        io.write_json(manifest, run.metadata_dir / "script01_manifest.json")

        crd.log_section(logger, "SCRIPT 1 COMPLETE", "=")
        logger.info("Elapsed time: %.2f min", elapsed / 60.0)
        logger.info("Quality flags: %s", quality_flags if quality_flags else "none")
        logger.info("Inspect Script-1 diagnostics before proceeding to Script 2.")

        if quality_flags and bool(cfg.FAIL_ON_WARNING):
            logger.error("FAIL_ON_WARNING=True and Script 1 produced quality flags.")
            return 2
        return 0

    except Exception:
        logger.exception("Script 1 failed. Partial logs/products remain in %s", run.run_dir)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
