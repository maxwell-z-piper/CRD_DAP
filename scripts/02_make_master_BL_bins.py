#!/usr/bin/env python3
"""CRD_DAP Script 2: create the BL-defined master PowerBins.

This stage runs Cappellari's PowerBin exactly once on the BL cube and transfers
that same physical sky membership to the red/RH3 native grid.  It does not fit
stellar kinematics or populations.  Its job is to define the immutable spatial
apertures and extract one matched BL/red spectrum per aperture for all later
science stages.

Primary responsibilities
------------------------
1. Load the exact Script-1 prepared BL and red/RH3 cubes plus tangent-plane
   coordinate products.
2. Measure a robust BL continuum S/N proxy in a configured rest-frame window.
3. Define a useful stellar-body aperture without discarding low-S/N spaxels
   inside that aperture.
4. Run mandatory PowerBin on BL only at the configured target S/N.
5. Transfer BL membership to red/RH3 through the celestial WCS, without
   independently re-binning the red arm.
6. Geometrically coadd BL and red/RH3 spectra, formal diagonal variances, and
   continuum-light spatial weights for every bin.
7. Save bin maps, membership, per-bin metadata, spectra, inherited preliminary
   Script-1 spectral-correlation diagnostics, figures, and full provenance.

Statistical caveat
------------------
Script 1 demonstrated spectral correlation from resampling/stacking, but it did
not calibrate a spatial covariance law between neighboring KcwiKit spaxels.
Therefore Script 2 uses formal diagonal spatial variances by default and records
that limitation explicitly.  PowerBin's callable capacity is retained so a
validated non-additive spatial covariance law can be inserted later without
changing the tessellation architecture.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import platform
import sys
import time

import numpy as np
from astropy.io import fits
from astropy.table import Table

import crd_utils as crd
from crd_utils import binning, io, plotting


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create BL-defined master PowerBins and transfer them to RH3."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to a target-specific config derived from config/target_config_template.py",
    )
    parser.add_argument(
        "--script1-run",
        default=None,
        help=(
            "Explicit Script-1 run directory containing products/prepared_BL.fits. "
            "If omitted, the newest complete Script-1 run for the target is used."
        ),
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional explicit Script-2 output run-directory name.",
    )
    return parser


def _quality_flag(flags: list[str], name: str, logger, message: str) -> None:
    if name not in flags:
        flags.append(name)
    logger.warning("%s | %s", name, message)


def _log_sn_summary(logger, arm: str, diag: binning.BinSNDiagnostics) -> None:
    """Log only scientifically defined positive-continuum achieved S/N values."""
    valid = np.asarray(diag.sn, dtype=float)
    finite = np.isfinite(valid)
    if np.any(finite):
        logger.info(
            "Achieved %s bin S/N (robust ratio-of-medians): valid=%d/%d, median=%.2f, min=%.2f, max=%.2f",
            arm,
            int(np.sum(finite)),
            int(valid.size),
            float(np.nanmedian(valid)),
            float(np.nanmin(valid)),
            float(np.nanmax(valid)),
        )
    else:
        logger.warning(
            "Achieved %s bin S/N: no bins have a valid positive-continuum S/N in the configured window",
            arm,
        )


def _sn_qc_masks(diag: binning.BinSNDiagnostics, cfg) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return non-positive, extreme-value, and estimator-disagreement bin masks."""
    enough = np.asarray(diag.n_good_channels, dtype=int) >= int(cfg.BIN_SN_MIN_GOOD_CHANNELS)
    nonpositive = enough & np.isfinite(diag.signed_sn) & (~np.asarray(diag.positive_continuum, dtype=bool))

    extreme_limit = float(cfg.BIN_SN_EXTREME_ABS_WARNING)
    extreme = (
        (np.isfinite(diag.signed_sn) & (np.abs(diag.signed_sn) > extreme_limit))
        | (np.isfinite(diag.legacy_median_ratio) & (np.abs(diag.legacy_median_ratio) > extreme_limit))
    )

    signed_abs = np.abs(np.asarray(diag.signed_sn, dtype=float))
    legacy_abs = np.abs(np.asarray(diag.legacy_median_ratio, dtype=float))
    small = np.minimum(signed_abs, legacy_abs)
    large = np.maximum(signed_abs, legacy_abs)
    factor = float(cfg.BIN_SN_ESTIMATOR_DISAGREEMENT_FACTOR)
    disagreement = (
        np.isfinite(signed_abs)
        & np.isfinite(legacy_abs)
        & (large > 1.0)
        & (large > factor * np.maximum(small, 1.0e-12))
    )
    return nonpositive, extreme, disagreement


def _format_bin_examples(diag: binning.BinSNDiagnostics, mask: np.ndarray, max_bins: int = 8) -> str:
    """Compact diagnostic string for logs without flooding a large-bin run."""
    ids = np.flatnonzero(mask)[: int(max_bins)]
    if ids.size == 0:
        return ""
    parts = []
    for bid in ids:
        parts.append(
            f"{int(bid)}(signed={diag.signed_sn[bid]:.3g}, legacy={diag.legacy_median_ratio[bid]:.3g}, "
            f"medflux={diag.median_flux[bid]:.3g}, medunc={diag.median_uncertainty[bid]:.3g}, "
            f"nchan={int(diag.n_good_channels[bid])})"
        )
    return ", ".join(parts)


def _required_script1_products(run_dir: Path) -> dict[str, Path]:
    products = run_dir / "products"
    return {
        "prepared_bl": products / "prepared_BL.fits",
        "prepared_rh3": products / "prepared_RH3.fits",
        "bl_coords": products / "BL_spatial_coordinates.npz",
        "rh3_coords": products / "RH3_spatial_coordinates.npz",
        "bl_noise": products / "BL_noise_diagnostic.npz",
        "rh3_noise": products / "RH3_noise_diagnostic.npz",
    }


def _is_complete_script1_run(run_dir: Path) -> bool:
    required = _required_script1_products(run_dir)
    core = ("prepared_bl", "prepared_rh3", "bl_coords", "rh3_coords")
    return all(required[key].is_file() for key in core)


def _find_script1_run(cfg, explicit: str | None) -> Path:
    if explicit is not None:
        run = Path(explicit).expanduser().resolve()
        if not _is_complete_script1_run(run):
            missing = [
                str(path)
                for key, path in _required_script1_products(run).items()
                if key in {"prepared_bl", "prepared_rh3", "bl_coords", "rh3_coords"}
                and not path.is_file()
            ]
            raise FileNotFoundError(
                "Explicit Script-1 run is incomplete. Missing:\n  - " + "\n  - ".join(missing)
            )
        return run

    root = Path(cfg.RUNS_ROOT).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"RUNS_ROOT does not exist: {root}")
    safe_target = str(cfg.TARGET_NAME).replace(" ", "_").replace("/", "-")
    candidates = [
        p for p in root.iterdir()
        if p.is_dir() and p.name.startswith(safe_target) and _is_complete_script1_run(p)
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No complete Script-1 run found for target {cfg.TARGET_NAME!r} under {root}. "
            "Pass --script1-run explicitly."
        )
    # Timestamped run names sort chronologically, but mtime is a safer fallback
    # if a run was copied or renamed.
    candidates.sort(key=lambda p: (p.stat().st_mtime, p.name), reverse=True)
    return candidates[0]


def _create_script2_run(cfg, run_name: str | None):
    if run_name is None:
        safe_target = str(cfg.TARGET_NAME).replace(" ", "_").replace("/", "-")
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_name = f"{safe_target}_S02_{stamp}"
    return crd.create_run_context(cfg, run_name=run_name)


def _load_spatial_coordinates(path: Path, expected_shape: tuple[int, int]) -> dict[str, np.ndarray | float]:
    with np.load(path) as data:
        x = np.asarray(data["x_arcsec"], dtype=float)
        y = np.asarray(data["y_arcsec"], dtype=float)
        if x.shape != expected_shape or y.shape != expected_shape:
            raise ValueError(
                f"Spatial-coordinate product {path} has shape {x.shape}/{y.shape}, "
                f"expected {expected_shape}"
            )
        return {
            "x_arcsec": x,
            "y_arcsec": y,
            "center_ra_deg": float(data["center_ra_deg"]),
            "center_dec_deg": float(data["center_dec_deg"]),
        }


def _nearest_zero_pixel(x_arcsec: np.ndarray, y_arcsec: np.ndarray) -> tuple[float, float]:
    radius = np.hypot(np.asarray(x_arcsec, dtype=float), np.asarray(y_arcsec, dtype=float))
    if not np.any(np.isfinite(radius)):
        raise ValueError("No finite tangent-plane coordinates available to locate adopted center")
    idx = int(np.nanargmin(radius))
    y, x = np.unravel_index(idx, radius.shape)
    return float(y), float(x)


def _pixel_area_arcsec2(cube: io.KCWICube) -> tuple[float, float, float]:
    sx, sy = cube.pixel_scales_arcsec()
    return float(sx * sy), float(sx), float(sy)


def _check_spatial_sampling(bl: io.KCWICube, rh3: io.KCWICube, cfg) -> None:
    bl_area, bl_sx, bl_sy = _pixel_area_arcsec2(bl)
    rh_area, rh_sx, rh_sy = _pixel_area_arcsec2(rh3)
    tol = float(cfg.BIN_TRANSFER_MAX_PIXEL_SCALE_FRACTIONAL_DIFFERENCE)
    diffs = [
        abs(rh_sx - bl_sx) / bl_sx,
        abs(rh_sy - bl_sy) / bl_sy,
        abs(rh_area - bl_area) / bl_area,
    ]
    if max(diffs) > tol:
        raise ValueError(
            "BL and RH3 prepared cubes have materially different spatial pixel scales. "
            f"BL=({bl_sx:.4f}, {bl_sy:.4f}) arcsec/pix, "
            f"RH3=({rh_sx:.4f}, {rh_sy:.4f}) arcsec/pix, "
            f"allowed fractional difference={tol:.3f}. "
            "Nearest-pixel membership transfer is not valid for this geometry; implement "
            "area-overlap apertures rather than silently proceeding."
        )


def _load_preliminary_corr(path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    if not path.is_file():
        return None
    with np.load(path) as data:
        return np.asarray(data["lags"], dtype=int), np.asarray(data["correlation"], dtype=float)


def _load_script1_manifest(run_dir: Path) -> dict:
    """Read Script-1 machine-readable provenance when available.

    Script 2 does not reinterpret Script-1 quality flags as automatic failures.
    It records them so the exact upstream QC state travels with the master bins.
    In particular, an inconclusive optional morphology cross-correlation does not
    erase the celestial WCS; Script 2 still performs its own strict WCS membership
    transfer and saves the transfer-distance diagnostic.
    """
    path = run_dir / "metadata" / "script01_manifest.json"
    if not path.is_file():
        return {}
    return dict(io.read_json(path))


def _bin_pixel_centroids(bin_map: np.ndarray, weights: np.ndarray, n_bins: int):
    y_flux = np.full(n_bins, np.nan)
    x_flux = np.full(n_bins, np.nan)
    y_geom = np.full(n_bins, np.nan)
    x_geom = np.full(n_bins, np.nan)
    for bid in range(n_bins):
        yy, xx = np.where(bin_map == bid)
        if yy.size == 0:
            continue
        y_geom[bid] = float(np.mean(yy))
        x_geom[bid] = float(np.mean(xx))
        w = np.asarray(weights[yy, xx], dtype=float)
        sw = float(np.sum(w))
        if sw > 0:
            y_flux[bid] = float(np.sum(yy * w) / sw)
            x_flux[bid] = float(np.sum(xx * w) / sw)
    return y_geom, x_geom, y_flux, x_flux


def _save_master_maps(
    run,
    bl: io.KCWICube,
    result: binning.PowerBinResult,
    aperture: binning.ApertureResult,
    transfer: binning.TransferResult,
    bl_sn: np.ndarray,
    rh3_sn: np.ndarray,
    bin_area: np.ndarray,
    *,
    source_script1_run: Path,
) -> Path:
    header = bl.celestial_wcs.to_header()
    header["CRDDAP"] = (True, "Processed by CRD_DAP")
    header["CRDSTEP"] = (2, "CRD_DAP master spatial binning")
    header["CRDARM"] = ("BL", "Spatial tessellation defined from BL")
    header["NBIN"] = (int(result.n_bins), "Number of BL master PowerBins")
    header["SRCSTEP1"] = (source_script1_run.name[:68], "Source Script-1 run directory")
    primary = fits.PrimaryHDU(np.asarray(result.bin_map, dtype=np.int32), header=header)
    primary.header["EXTNAME"] = "BINMAP"
    hdus = [primary]
    hdus.append(fits.ImageHDU(aperture.mask.astype(np.uint8), name="APERTURE"))
    hdus.append(
        fits.ImageHDU(
            np.where(result.bin_map >= 0, bl_sn[result.bin_map.clip(min=0)], np.nan).astype(np.float32),
            name="BL_SN",
        )
    )
    hdus.append(
        fits.ImageHDU(
            np.where(result.bin_map >= 0, rh3_sn[result.bin_map.clip(min=0)], np.nan).astype(np.float32),
            name="RH3_SN",
        )
    )
    hdus.append(
        fits.ImageHDU(
            np.where(result.bin_map >= 0, bin_area[result.bin_map.clip(min=0)], np.nan).astype(np.float32),
            name="BINAREA",
        )
    )
    hdus[-1].header["BUNIT"] = "arcsec2"
    hdus.append(fits.ImageHDU(np.asarray(transfer.bin_map, dtype=np.int32), name="RH3BIN"))
    hdus.append(
        fits.ImageHDU(np.asarray(transfer.match_distance_arcsec, dtype=np.float32), name="RH3DIST")
    )
    hdus[-1].header["BUNIT"] = "arcsec"
    path = run.products_dir / "master_bins.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def _save_spectra(
    run,
    bl: io.KCWICube,
    rh3: io.KCWICube,
    bl_spec: binning.CoaddedBinSpectra,
    rh3_spec: binning.CoaddedBinSpectra,
    bl_corr,
    rh3_corr,
) -> Path:
    header = fits.Header()
    header["CRDDAP"] = True
    header["CRDSTEP"] = 2
    header["NBIN"] = int(bl_spec.flux.shape[0])
    header["BLSCALE"] = (float(bl_spec.spatial_scale_factor), "BL spatial sum scale per spaxel")
    header["RHSCALE"] = (float(rh3_spec.spatial_scale_factor), "RH3 spatial sum scale per spaxel")
    hdus = [fits.PrimaryHDU(header=header)]

    def add_arm(prefix, cube, spec):
        hdus.append(fits.ImageHDU(np.asarray(cube.wavelength, dtype=np.float64), name=f"{prefix}_WAVE"))
        hdus[-1].header["BUNIT"] = "Angstrom"
        hdus.append(fits.ImageHDU(np.asarray(spec.flux, dtype=np.float32), name=f"{prefix}_FLUX"))
        hdus.append(fits.ImageHDU(np.asarray(spec.uncertainty, dtype=np.float32), name=f"{prefix}_UNCERT"))
        hdus.append(fits.ImageHDU(np.asarray(spec.good, dtype=np.uint8), name=f"{prefix}_GOOD"))
        hdus.append(
            fits.ImageHDU(np.asarray(spec.contributing_spaxels, dtype=np.int16), name=f"{prefix}_NCONTRIB")
        )
        bunit = str(cube.header.get("BUNIT", ""))
        hdus[-4].header["SRCBUNIT"] = bunit[:68]
        hdus[-3].header["SRCBUNIT"] = bunit[:68]

    add_arm("BL", bl, bl_spec)
    add_arm("RH3", rh3, rh3_spec)
    for prefix, corr in (("BL", bl_corr), ("RH3", rh3_corr)):
        if corr is None:
            continue
        lags, values = corr
        hdus.append(fits.ImageHDU(np.asarray(lags, dtype=np.int16), name=f"{prefix}_CORRLAG"))
        hdus.append(fits.ImageHDU(np.asarray(values, dtype=np.float32), name=f"{prefix}_CORR"))
        hdus[-1].header["COMMENT"] = "Preliminary Script-1 high-pass spectral correlation; not a final covariance matrix"
    path = run.products_dir / "master_bin_spectra.fits"
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def main() -> int:
    args = _parser().parse_args()
    cfg = crd.load_config(args.config, validate=True, strict_paths=True)
    source_run = _find_script1_run(cfg, args.script1_run)
    source_products = _required_script1_products(source_run)
    source_manifest = _load_script1_manifest(source_run)
    run = _create_script2_run(cfg, args.run_name)
    crd.setup_pipeline_logger(run)
    logger = crd.setup_step_logger(run, "02_make_master_BL_bins")
    quality_flags: list[str] = []
    start = time.perf_counter()

    try:
        crd.log_section(logger, "CRD_DAP SCRIPT 2: MAKE MASTER BL POWERBINS", "=")
        logger.info("Run directory: %s", run.run_dir)
        logger.info("Source Script-1 run: %s", source_run)
        logger.info("Python: %s", sys.version.replace("\n", " "))
        logger.info("Platform: %s", platform.platform())
        logger.info("Target: %s", cfg.TARGET_NAME)
        logger.info("Redshift: %.8f", float(cfg.REDSHIFT))
        source_flags = list(source_manifest.get("quality_flags", []))
        logger.info("Source Script-1 quality flags: %s", source_flags if source_flags else "none/manifest unavailable")
        if "REGISTRATION_INCONCLUSIVE" in source_flags:
            logger.warning(
                "SOURCE_REGISTRATION_INCONCLUSIVE | Script 1 did not adopt an optional morphology-based residual shift. "
                "Script 2 will therefore use the saved celestial WCS directly and will QC every RH3-to-BL membership "
                "assignment through BL_RH3_bin_transfer.png and the configured transfer-distance tolerance."
            )

        crd.log_section(logger, "1. Load Script-1 prepared products")
        bl = io.load_prepared_cube(source_products["prepared_bl"], expected_arm="BL")
        rh3 = io.load_prepared_cube(source_products["prepared_rh3"], expected_arm="RH3")
        _check_spatial_sampling(bl, rh3, cfg)
        bl_coords = _load_spatial_coordinates(source_products["bl_coords"], bl.good_spaxel.shape)
        rh3_coords = _load_spatial_coordinates(source_products["rh3_coords"], rh3.good_spaxel.shape)
        bl_x = bl_coords["x_arcsec"]
        bl_y = bl_coords["y_arcsec"]
        rh3_x = rh3_coords["x_arcsec"]
        rh3_y = rh3_coords["y_arcsec"]
        center_yx = _nearest_zero_pixel(bl_x, bl_y)
        bl_area, bl_sx, bl_sy = _pixel_area_arcsec2(bl)
        rh3_area, rh3_sx, rh3_sy = _pixel_area_arcsec2(rh3)
        logger.info(
            "Prepared shapes: BL=%s, RH3=%s | pixel scales BL=(%.4f, %.4f), RH3=(%.4f, %.4f) arcsec/pix",
            bl.shape,
            rh3.shape,
            bl_sx,
            bl_sy,
            rh3_sx,
            rh3_sy,
        )
        logger.info("Adopted Script-1 tangent-plane center maps to BL pixel (y,x)=(%.2f, %.2f)", *center_yx)

        crd.log_section(logger, "2. Continuum S/N maps and stellar-body aperture")
        bl_window = binning.continuum_window_maps(
            bl,
            rest_range=tuple(cfg.BL_BINNING_REST_RANGE_ANGSTROM),
            redshift=float(cfg.REDSHIFT),
            min_valid_fraction=float(cfg.BINNING_MIN_VALID_WINDOW_FRACTION),
        )
        rh3_window = binning.continuum_window_maps(
            rh3,
            rest_range=tuple(cfg.RH3_SN_REST_RANGE_ANGSTROM),
            redshift=float(cfg.REDSHIFT),
            min_valid_fraction=float(cfg.BINNING_MIN_VALID_WINDOW_FRACTION),
        )
        logger.info(
            "BL binning S/N window: observed %.1f--%.1f A (%d good channels)",
            bl_window.observed_min,
            bl_window.observed_max,
            bl_window.n_channels,
        )
        logger.info(
            "RH3 S/N window: observed %.1f--%.1f A (%d good channels)",
            rh3_window.observed_min,
            rh3_window.observed_max,
            rh3_window.n_channels,
        )
        aperture = binning.make_analysis_aperture(
            bl_window.signal,
            bl_window.noise,
            bl.good_spaxel,
            center_yx=center_yx,
            mode=str(cfg.BINNING_APERTURE_MODE),
            smooth_sigma_pix=float(cfg.BINNING_APERTURE_SMOOTH_SIGMA_PIX),
            threshold=float(cfg.BINNING_APERTURE_SN_THRESHOLD),
            dilate_pix=int(cfg.BINNING_APERTURE_DILATE_PIX),
            center_max_distance_pix=float(cfg.BINNING_APERTURE_CENTER_MAX_DISTANCE_PIX),
            min_pixels=int(cfg.BINNING_APERTURE_MIN_PIXELS),
            max_radius_arcsec=getattr(cfg, "BINNING_APERTURE_MAX_RADIUS_ARCSEC", None),
            x_arcsec=bl_x,
            y_arcsec=bl_y,
        )
        logger.info(
            "BL binning aperture: %d pixels | mode=%s | threshold=%.3f | nearest detected component to center=%.2f pix",
            aperture.n_pixels,
            cfg.BINNING_APERTURE_MODE,
            aperture.threshold,
            aperture.nearest_component_distance_pix,
        )
        plotting.plot_binning_aperture(
            aperture.significance_proxy,
            aperture.mask,
            bl_x,
            bl_y,
            run.figures_dir / "binning_aperture.png",
            threshold=aperture.threshold,
        )

        crd.log_section(logger, "3. Run BL PowerBin tessellation")
        cov_mode = str(cfg.POWERBIN_SPATIAL_COVARIANCE_MODE).lower()
        cov_alpha = float(cfg.POWERBIN_SPATIAL_COVARIANCE_ALPHA)
        if cov_mode == "none":
            _quality_flag(
                quality_flags,
                "SPATIAL_COVARIANCE_UNCALIBRATED",
                logger,
                "PowerBin currently uses formal diagonal spatial variance. Script 1 measured spectral correlation, not a validated spatial covariance law; no empirical correction is invented here.",
            )
        pb = binning.run_powerbin(
            bl_x,
            bl_y,
            bl_window.signal,
            bl_window.noise,
            aperture.mask,
            target_sn=float(cfg.BL_TARGET_SN),
            pixel_size_arcsec=float(np.sqrt(bl_area)),
            covariance_mode=cov_mode,
            covariance_alpha=cov_alpha,
            regul=bool(cfg.POWERBIN_REGUL),
            maxiter=int(cfg.POWERBIN_MAXITER),
            verbose=int(cfg.POWERBIN_VERBOSE),
        )
        logger.info(
            "PowerBin complete: nbin=%d, input pixels=%d, version=%s, target S/N=%.2f, capacity RMS scatter=%.2f%%",
            pb.n_bins,
            aperture.n_pixels,
            pb.powerbin_version,
            float(cfg.BL_TARGET_SN),
            pb.rms_frac_percent,
        )
        plotting.plot_master_bins(
            pb.bin_map,
            bl_x,
            bl_y,
            run.figures_dir / "master_bins.png",
            pa_kin_deg=getattr(cfg, "PA_KIN_INITIAL_DEG", None),
        )

        crd.log_section(logger, "4. Transfer BL physical membership to RH3")
        transfer = binning.transfer_bin_map_by_wcs(
            pb.bin_map,
            bl.celestial_wcs,
            rh3.celestial_wcs,
            rh3.good_spaxel,
            bl_pixel_scale_arcsec=float(np.sqrt(bl_area)),
            max_distance_arcsec=float(cfg.BIN_TRANSFER_MAX_DISTANCE_ARCSEC),
        )
        logger.info(
            "RH3 membership transfer: assigned %d/%d candidate spaxels (%.2f%%) within %.3f arcsec",
            transfer.n_assigned_spaxels,
            transfer.n_candidate_spaxels,
            100.0 * transfer.assigned_fraction,
            float(cfg.BIN_TRANSFER_MAX_DISTANCE_ARCSEC),
        )
        if transfer.assigned_fraction < float(cfg.BIN_TRANSFER_MIN_ASSIGNED_FRACTION):
            _quality_flag(
                quality_flags,
                "BIN_TRANSFER_INCOMPLETE",
                logger,
                f"Only {100.0*transfer.assigned_fraction:.2f}% of RH3 candidate spaxels were assigned to BL bins; configured minimum is {100.0*float(cfg.BIN_TRANSFER_MIN_ASSIGNED_FRACTION):.2f}%.",
            )
        plotting.plot_bin_transfer(
            pb.bin_map,
            transfer.bin_map,
            transfer.match_distance_arcsec,
            run.figures_dir / "BL_RH3_bin_transfer.png",
            max_distance_arcsec=float(cfg.BIN_TRANSFER_MAX_DISTANCE_ARCSEC),
        )

        crd.log_section(logger, "5. Coadd matched BL and RH3 spectra")
        bl_spec = binning.coadd_bin_spectra(
            bl,
            pb.bin_map,
            n_bins=pb.n_bins,
            pixel_area_arcsec2=bl_area,
            min_member_fraction=float(cfg.BIN_SPECTRUM_MIN_MEMBER_FRACTION),
        )
        rh3_spec = binning.coadd_bin_spectra(
            rh3,
            transfer.bin_map,
            n_bins=pb.n_bins,
            pixel_area_arcsec2=rh3_area,
            min_member_fraction=float(cfg.BIN_SPECTRUM_MIN_MEMBER_FRACTION),
        )
        # Achieved continuum S/N is a QC/reporting quantity, not the quantity
        # used to define the PowerBin tessellation.  Use the robust
        # ratio-of-medians estimator and retain the older median(flux/uncertainty)
        # statistic only as an audit diagnostic.  This prevents a negative
        # continuum combined with tiny formal uncertainties from appearing as a
        # physically meaningful enormous negative "achieved S/N".
        bl_sn_diag = binning.achieved_sn_diagnostics_per_bin(
            bl_spec,
            bl.wavelength,
            rest_range=tuple(cfg.BL_BINNING_REST_RANGE_ANGSTROM),
            redshift=float(cfg.REDSHIFT),
            min_good_channels=int(cfg.BIN_SN_MIN_GOOD_CHANNELS),
            require_positive_continuum=bool(cfg.BIN_SN_REQUIRE_POSITIVE_CONTINUUM),
        )
        rh3_sn_diag = binning.achieved_sn_diagnostics_per_bin(
            rh3_spec,
            rh3.wavelength,
            rest_range=tuple(cfg.RH3_SN_REST_RANGE_ANGSTROM),
            redshift=float(cfg.REDSHIFT),
            min_good_channels=int(cfg.BIN_SN_MIN_GOOD_CHANNELS),
            require_positive_continuum=bool(cfg.BIN_SN_REQUIRE_POSITIVE_CONTINUUM),
        )
        bl_sn = bl_sn_diag.sn
        rh3_sn = rh3_sn_diag.sn

        _log_sn_summary(logger, "BL", bl_sn_diag)
        _log_sn_summary(logger, "RH3", rh3_sn_diag)

        # Preserve unusual signed values as explicit QC rather than silently
        # clipping or modifying the underlying spectra/variance.
        for arm, diag in (("BL", bl_sn_diag), ("RH3", rh3_sn_diag)):
            nonpositive, extreme, disagreement = _sn_qc_masks(diag, cfg)
            if np.any(nonpositive):
                _quality_flag(
                    quality_flags,
                    "NONPOSITIVE_BIN_CONTINUUM",
                    logger,
                    f"{arm}: {int(np.sum(nonpositive))}/{pb.n_bins} bins have non-positive median continuum in the configured S/N window; production-facing achieved S/N is NaN for those bins. "
                    f"Examples: {_format_bin_examples(diag, nonpositive)}",
                )
            if np.any(extreme):
                _quality_flag(
                    quality_flags,
                    "EXTREME_BIN_SN_DIAGNOSTIC",
                    logger,
                    f"{arm}: {int(np.sum(extreme))}/{pb.n_bins} bins exceed |S/N|>{float(cfg.BIN_SN_EXTREME_ABS_WARNING):g} in the signed robust or legacy estimator. "
                    f"No data are clipped; inspect formal uncertainties/window placement. Examples: {_format_bin_examples(diag, extreme)}",
                )
            if np.any(disagreement):
                _quality_flag(
                    quality_flags,
                    "BIN_SN_ESTIMATOR_DISAGREEMENT",
                    logger,
                    f"{arm}: {int(np.sum(disagreement))}/{pb.n_bins} bins show >{float(cfg.BIN_SN_ESTIMATOR_DISAGREEMENT_FACTOR):g}x disagreement between robust ratio-of-medians and legacy median(flux/uncertainty). "
                    f"Examples: {_format_bin_examples(diag, disagreement)}",
                )

        low_limit = float(cfg.BL_TARGET_SN) * float(cfg.BINNING_LOW_SN_WARNING_FRACTION)
        low_bins = np.flatnonzero(np.isfinite(bl_sn) & (bl_sn < low_limit))
        if low_bins.size:
            _quality_flag(
                quality_flags,
                "LOW_BL_BIN_SN",
                logger,
                f"{low_bins.size}/{pb.n_bins} BL bins have measured S/N below {low_limit:.2f} ({float(cfg.BINNING_LOW_SN_WARNING_FRACTION):.2f} x target).",
            )

        plotting.plot_bin_value_map(
            pb.bin_map,
            bl_sn,
            run.figures_dir / "BL_SN_per_bin.png",
            title="BL achieved S/N per master PowerBin",
            colorbar_label="Robust continuum S/N per spectral pixel",
            vmin=0.0,
        )
        plotting.plot_bin_value_map(
            pb.bin_map,
            rh3_sn,
            run.figures_dir / "RH3_SN_per_bin.png",
            title="RH3 achieved S/N in BL-defined PowerBins",
            colorbar_label="Robust continuum S/N per spectral pixel",
            vmin=0.0,
        )
        plotting.plot_bl_rh3_sn_comparison(
            pb.bin_map,
            bl_sn,
            rh3_sn,
            run.figures_dir / "BL_RH3_SN_comparison.png",
            upper_percentile=float(cfg.SN_PLOT_UPPER_PERCENTILE),
        )

        crd.log_section(logger, "6. Geometry, light weights, and per-bin table")
        bl_weights = binning.normalized_flux_weights(bl_window.signal, pb.bin_map, n_bins=pb.n_bins)
        rh3_weights = binning.normalized_flux_weights(
            rh3_window.signal,
            transfer.bin_map,
            n_bins=pb.n_bins,
        )
        cent = binning.bin_centroids(pb.bin_map, bl_x, bl_y, bl_weights, n_bins=pb.n_bins)
        y_geom_pix, x_geom_pix, y_flux_pix, x_flux_pix = _bin_pixel_centroids(
            pb.bin_map, bl_weights, pb.n_bins
        )
        sky = bl.celestial_wcs.pixel_to_world(x_flux_pix, y_flux_pix)
        n_bl = np.asarray(pb.npix_per_bin, dtype=int)
        n_rh = np.asarray([np.sum(transfer.bin_map == bid) for bid in range(pb.n_bins)], dtype=int)
        area = n_bl.astype(float) * bl_area
        rh3_coverage = np.divide(n_rh, n_bl, out=np.full(pb.n_bins, np.nan), where=n_bl > 0)
        table = Table()
        table["BIN_ID"] = np.arange(pb.n_bins, dtype=int)
        table["NPIX_BL"] = n_bl
        table["NPIX_RH3"] = n_rh
        table["AREA_ARCSEC2"] = area
        table["X_GEOM_ARCSEC"] = cent["x_geom"]
        table["Y_GEOM_ARCSEC"] = cent["y_geom"]
        table["X_FLUX_ARCSEC"] = cent["x_flux"]
        table["Y_FLUX_ARCSEC"] = cent["y_flux"]
        table["RA_DEG"] = np.asarray(sky.ra.deg, dtype=float)
        table["DEC_DEG"] = np.asarray(sky.dec.deg, dtype=float)
        table["BL_SN"] = bl_sn
        table["RH3_SN"] = rh3_sn
        # Signed/audit columns are intentionally verbose: they make a future
        # pathological S/N immediately traceable without reopening the cubes.
        table["BL_SN_SIGNED"] = bl_sn_diag.signed_sn
        table["RH3_SN_SIGNED"] = rh3_sn_diag.signed_sn
        table["BL_SN_LEGACY"] = bl_sn_diag.legacy_median_ratio
        table["RH3_SN_LEGACY"] = rh3_sn_diag.legacy_median_ratio
        table["BL_SN_MEDFLUX"] = bl_sn_diag.median_flux
        table["RH3_SN_MEDFLUX"] = rh3_sn_diag.median_flux
        table["BL_SN_MEDUNC"] = bl_sn_diag.median_uncertainty
        table["RH3_SN_MEDUNC"] = rh3_sn_diag.median_uncertainty
        table["BL_SN_MINUNC"] = bl_sn_diag.min_uncertainty
        table["RH3_SN_MINUNC"] = rh3_sn_diag.min_uncertainty
        table["BL_SN_NCHAN"] = bl_sn_diag.n_good_channels
        table["RH3_SN_NCHAN"] = rh3_sn_diag.n_good_channels
        table["BL_NEGFLUX_FRAC"] = bl_sn_diag.negative_flux_fraction
        table["RH3_NEGFLUX_FRAC"] = rh3_sn_diag.negative_flux_fraction
        table["POWERBIN_SN"] = np.sqrt(np.clip(pb.bin_capacity, 0.0, None))
        table["POWERBIN_CAPACITY"] = pb.bin_capacity
        table["RH3_PIXEL_COVERAGE"] = rh3_coverage
        table["SINGLE_PIXEL"] = pb.single
        table.write(run.products_dir / "master_bin_table.ecsv", format="ascii.ecsv", overwrite=True)
        plotting.plot_bin_value_map(
            pb.bin_map,
            area,
            run.figures_dir / "bin_area_map.png",
            title="Master PowerBin physical area",
            colorbar_label="Bin area (arcsec²)",
        )

        # Save complete per-spaxel membership/weights in a compact numerical file
        # rather than a variable-length FITS table.  This is the authoritative
        # input for later bin-integrated disk-model calculations.
        bly, blx = np.where(pb.bin_map >= 0)
        rhy, rhx = np.where(transfer.bin_map >= 0)
        np.savez_compressed(
            run.products_dir / "master_bin_membership.npz",
            bl_y_pix=bly.astype(np.int32),
            bl_x_pix=blx.astype(np.int32),
            bl_bin_id=pb.bin_map[bly, blx].astype(np.int32),
            bl_x_arcsec=np.asarray(bl_x[bly, blx], dtype=float),
            bl_y_arcsec=np.asarray(bl_y[bly, blx], dtype=float),
            bl_flux_weight=np.asarray(bl_weights[bly, blx], dtype=float),
            rh3_y_pix=rhy.astype(np.int32),
            rh3_x_pix=rhx.astype(np.int32),
            rh3_bin_id=transfer.bin_map[rhy, rhx].astype(np.int32),
            rh3_x_arcsec=np.asarray(rh3_x[rhy, rhx], dtype=float),
            rh3_y_arcsec=np.asarray(rh3_y[rhy, rhx], dtype=float),
            rh3_flux_weight=np.asarray(rh3_weights[rhy, rhx], dtype=float),
            rh3_match_distance_arcsec=np.asarray(transfer.match_distance_arcsec[rhy, rhx], dtype=float),
            generator_xy_arcsec=np.asarray(pb.generator_xy, dtype=float),
            generator_radius_arcsec=np.asarray(pb.generator_radius, dtype=float),
            powerbin_capacity=np.asarray(pb.bin_capacity, dtype=float),
        )

        crd.log_section(logger, "7. Save master-bin products and provenance")
        maps_path = _save_master_maps(
            run,
            bl,
            pb,
            aperture,
            transfer,
            bl_sn,
            rh3_sn,
            area,
            source_script1_run=source_run,
        )
        bl_corr = _load_preliminary_corr(source_products["bl_noise"])
        rh3_corr = _load_preliminary_corr(source_products["rh3_noise"])
        spectra_path = _save_spectra(run, bl, rh3, bl_spec, rh3_spec, bl_corr, rh3_corr)

        manifest = {
            "script": "02_make_master_BL_bins",
            "target": str(cfg.TARGET_NAME),
            "source_script1_run": str(source_run),
            "source_script1_quality_flags": list(source_manifest.get("quality_flags", [])),
            "prepared_bl": str(source_products["prepared_bl"]),
            "prepared_rh3": str(source_products["prepared_rh3"]),
            "n_bins": int(pb.n_bins),
            "powerbin_version": str(pb.powerbin_version),
            "powerbin_capacity_rms_frac_percent": float(pb.rms_frac_percent),
            "n_bl_aperture_pixels": int(aperture.n_pixels),
            "bl_target_sn": float(cfg.BL_TARGET_SN),
            "bl_window_observed_angstrom": [bl_window.observed_min, bl_window.observed_max],
            "rh3_window_observed_angstrom": [rh3_window.observed_min, rh3_window.observed_max],
            "bin_sn_estimator": "median(flux) / median(uncertainty)",
            "bin_sn_require_positive_continuum": bool(cfg.BIN_SN_REQUIRE_POSITIVE_CONTINUUM),
            "bin_sn_min_good_channels": int(cfg.BIN_SN_MIN_GOOD_CHANNELS),
            "n_bl_valid_positive_sn": int(np.sum(np.isfinite(bl_sn))),
            "n_rh3_valid_positive_sn": int(np.sum(np.isfinite(rh3_sn))),
            "n_bl_nonpositive_continuum": int(np.sum(~bl_sn_diag.positive_continuum & (bl_sn_diag.n_good_channels >= int(cfg.BIN_SN_MIN_GOOD_CHANNELS)))),
            "n_rh3_nonpositive_continuum": int(np.sum(~rh3_sn_diag.positive_continuum & (rh3_sn_diag.n_good_channels >= int(cfg.BIN_SN_MIN_GOOD_CHANNELS)))),
            "powerbin_spatial_covariance_mode": cov_mode,
            "powerbin_spatial_covariance_alpha": cov_alpha,
            "rh3_transfer_assigned_fraction": float(transfer.assigned_fraction),
            "bl_spatial_scale_factor": float(bl_spec.spatial_scale_factor),
            "bl_spatial_scale_reason": bl_spec.spatial_scale_reason,
            "rh3_spatial_scale_factor": float(rh3_spec.spatial_scale_factor),
            "rh3_spatial_scale_reason": rh3_spec.spatial_scale_reason,
            "formal_spatial_variance_only": True,
            "preliminary_spectral_correlation_carried_forward": bool(bl_corr is not None or rh3_corr is not None),
            "quality_flags": quality_flags,
            "products": {
                "master_bins_fits": str(maps_path),
                "master_bin_spectra_fits": str(spectra_path),
                "master_bin_table_ecsv": str(run.products_dir / "master_bin_table.ecsv"),
                "master_bin_membership_npz": str(run.products_dir / "master_bin_membership.npz"),
            },
        }
        io.write_json(manifest, run.metadata_dir / "script02_manifest.json")
        logger.info("Master bin map: %s", maps_path)
        logger.info("Master bin spectra: %s", spectra_path)
        logger.info("Master bin table: %s", run.products_dir / "master_bin_table.ecsv")
        logger.info("Master membership: %s", run.products_dir / "master_bin_membership.npz")

        crd.log_section(logger, "SCRIPT 2 COMPLETE", "=")
        logger.info("Elapsed time: %.2f min", (time.perf_counter() - start) / 60.0)
        logger.info("Quality flags: %s", quality_flags)
        logger.info("Inspect Script-2 aperture, tessellation, transfer, and S/N diagnostics before Script 3.")
        return 0

    except Exception:
        logger.exception("Script 2 failed. Partial logs/products remain in %s", run.run_dir)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
