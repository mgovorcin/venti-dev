# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Configuration handling for Venti workflows.

This module provides configuration management split into two files:
1. runconfig.yaml: Inputs, outputs, product version (changes per run)
2. algorithm_parameters.yaml: Algorithm settings with defaults (rarely changes)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal

import yaml  # type: ignore[import-untyped]
from pydantic import BaseModel, Field, field_validator, model_validator

# ============================================================================
# Algorithm Parameters (algorithm_parameters.yaml)
# ============================================================================


# Schema version of algorithm_parameters.yaml written by this Venti (plan T19).
ALGORITHM_SCHEMA_VERSION = 2
logger = logging.getLogger(__name__)


class ProcessingOptions(BaseModel):
    """Common processing options used by multiple workflows.

    Attributes
    ----------
    cal_downsample_factor : int
        Downsample factor for plane fitting
    downsample_method : str
        Method for aggregating pixels during downsampling
    downsample_weighted : bool
        Weight downsampling by each product's ``temporal_coherence``
    vlm_output_posting_meters : float
        Output grid posting in meters for decomposed ENU components

    """

    cal_downsample_factor: int = Field(
        1,
        ge=1,
        description=(
            "Downsample factor for performing plane fitting at lower resolution"
        ),
    )
    downsample_method: Literal["mean", "median"] = Field(
        "mean",
        description="Method for aggregating pixels during downsampling",
    )
    downsample_weighted: bool = Field(
        False,
        description="Weight downsampling by each product's temporal_coherence",
    )
    vlm_output_posting_meters: float = Field(
        120.0,
        gt=0,
        description="Output grid posting in meters for decomposed ENU components",
    )


class SavitzkyGolayOptions(BaseModel):
    """Savitzky-Golay filter options.

    Attributes
    ----------
    window_length : int
        Window length in pixels (must be odd)
    polyorder : int
        Polynomial order for fitting

    """

    window_length: int = Field(
        51, ge=3, description="Window length in pixels (must be odd)"
    )
    polyorder: int = Field(3, ge=0, description="Polynomial order for fitting")


class SurfaceOptions(BaseModel):
    """Calibration-surface estimator (PRD R-S1..R-S3).

    Defaults reproduce the gamma 0.3 algorithm: one windowed plane fit, no gap
    filling, single pass. The trade-study algorithm is ``method='loclin'``
    with ``cutoff_wavelength_meters=50000``, ``fill_gaps=True``,
    ``two_pass=True`` (plan T29-T33).
    """

    method: Literal["windowed_plane", "loclin"] = Field(
        "windowed_plane",
        description=(
            "'windowed_plane': the gamma moving-window plane fit sized by "
            "window_size_meters; 'loclin': local-linear kernel with a physical "
            "half-response cutoff (cutoff_wavelength_meters)"
        ),
    )
    cutoff_wavelength_meters: float = Field(
        50000.0,
        gt=0,
        description=(
            "Half-response wavelength of the 'loclin' kernel in meters: signal "
            "longer than this is attributed to GNSS, shorter stays InSAR"
        ),
    )
    fill_gaps: bool = Field(
        False,
        description=(
            "Fill masked/water cells before the fit so the surface is continuous "
            "everywhere and never 0 on masked cells (R-S3)"
        ),
    )
    two_pass: bool = Field(
        False,
        description=(
            "Robust frame-wide tie first, unwrap correction on DISP - CAL1, then "
            "the final surface (R-S1)"
        ),
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}


class WeightOptions(BaseModel):
    """Fit weights (PRD R-S4). Defaults: unweighted, as in gamma."""

    coherence_power: float = Field(
        0.0,
        ge=0,
        description=(
            "Weight pixels by temporal_coherence ** p; 0 disables. The trade "
            "studies recommend 8 (TS-B1 may revise)"
        ),
    )
    robust: bool = Field(
        False,
        description="Local robust (Huber/MAD) re-weighting inside the kernel fit",
    )
    filled_pixel_weight: float = Field(
        0.02,
        ge=0,
        le=1,
        description="Weight given to gap-filled pixels when fill_gaps is on",
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}


class TropoOptions(BaseModel):
    """Tropospheric correction mode (PRD R-T1).

    ``'legacy'`` keeps the gamma behaviour: `apply_tropo_correction` decides,
    the full ZTD is removed. ``'auto'`` chooses from the DEM relief
    (p5-p95): off below `relief_off_meters`, stratified at or above
    `relief_stratified_meters`, off in between until trade study TS-T1.
    """

    mode: Literal["legacy", "off", "stratified", "full", "auto"] = Field(
        "legacy", description="Tropo mode; 'legacy' defers to apply_tropo_correction"
    )
    relief_off_meters: float = Field(
        300.0, ge=0, description="'auto': below this relief the correction is off"
    )
    relief_stratified_meters: float = Field(
        1500.0,
        ge=0,
        description="'auto': at or above this relief the stratified model is used",
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}

    @model_validator(mode="after")
    def _bands_ordered(self) -> TropoOptions:
        if self.relief_off_meters > self.relief_stratified_meters:
            msg = (
                "relief_off_meters must not exceed relief_stratified_meters "
                f"({self.relief_off_meters} > {self.relief_stratified_meters})"
            )
            raise ValueError(msg)
        return self


class GnssOptions(BaseModel):
    """GNSS grid handling (PRD R-G1, R-G2, R-G4, R-G5). Defaults as in gamma."""

    buffer_meters: float = Field(
        0.0,
        ge=0,
        description=(
            "Use grid nodes this far outside the frame (R-G2; v0.5 uses 50000). "
            "sample_gnss_enu interpolates E/N/U and projects per pixel, so nodes "
            "outside the swath need no LOS look"
        ),
    )
    exclude_defo_nodes: bool = Field(
        False,
        description="Drop grid nodes inside the defo/event areas (R-G5)",
    )
    reinterpolate_excluded: bool = Field(
        False,
        description=(
            "Re-estimate the excluded nodes from the surrounding ones with GPS "
            "Imaging (geepers.gps_imaging.reinterpolate_nodes)"
        ),
    )
    reprocessing: bool = Field(
        False,
        description=(
            "Reprocessing mode: the only mode in which grid_type='variable' is "
            "allowed (R-G1)"
        ),
    )
    snapshot_id: str | None = Field(
        None,
        description=(
            "Identifier of the frozen UNR grid snapshot in use (R-G4); recorded "
            "in the product metadata"
        ),
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}


class UncertaintyOptions(BaseModel):
    """calibration_std model (PRD R-E1, R-E2)."""

    k_grid: float | Literal["frame_table"] = Field(
        1.0,
        description=(
            "Inflation of the GNSS grid sigma: a number, or 'frame_table' to "
            "take the per-frame value from the frame-parameter table (TS-G1)"
        ),
    )
    inflate_inside_areas: bool = Field(
        True,
        description="Grow sigma with distance inside interpolated defo/event areas",
    )
    sigma_nonsecular_meters: float = Field(
        0.0,
        ge=0,
        description=(
            "Optional per-pair term for non-secular station motion (seasonal "
            "loading etc.) that a constant-velocity grid cannot carry; added in "
            "quadrature to sigma_CAL. 0 (default): sigma_CAL is the uncertainty of "
            "the calibration surface only. TS-G1 measured 2.6-4.9 mm on the four "
            "validation frames. Units of the displacement (metres); loclin route"
        ),
    )
    sigma_disp_placeholder_mm: float = Field(
        10.0,
        ge=0,
        description=(
            "Documented DISP noise below the cutoff that users add themselves; "
            "replaced by the TS-S1 model"
        ),
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}


class FrameOptions(BaseModel):
    """Frame-level facts that are not algorithm tuning (PRD §4.3, D3).

    Filled from the frame-parameter table (`venti.frames`); recorded in the
    product metadata.
    """

    plate: str = Field(
        "NA",
        description=(
            "Tectonic plate of the plate_motion layer: a UNR two-letter code "
            "(NA, PA, CA, ...) or an ITRF PMM name (NOAM, PCFC, CARB, ...)"
        ),
    )
    name: str | None = Field(None, description="Human-readable frame name")
    benchmark_category: str | None = Field(
        None, description="CalVal benchmark category of the frame (PRD 2.7), if any"
    )

    @field_validator("plate")
    @classmethod
    def _known_plate(cls, value: str) -> str:
        from geepers.euler import PLATE_CODES, load_plate_motion_model

        code = value.upper()
        if code in PLATE_CODES or code in load_plate_motion_model():
            return code
        msg = (
            f"Unknown plate {value!r}; codes: {sorted(PLATE_CODES)}, "
            f"names: {sorted(load_plate_motion_model())}"
        )
        raise ValueError(msg)

    model_config = {"validate_assignment": True, "extra": "forbid"}


class UnwrapOptions(BaseModel):
    """Unwrap-error correction details (PRD R-U1).

    `CalibrationOptions.unwrap_error_correction` is the on/off switch.
    """

    region_source: Literal["mask", "water_mask"] = Field(
        "water_mask",
        description=(
            "How regions are segmented: 'water_mask' watershed (islands, cut-off "
            "peninsulas) or the gamma 'mask' islands"
        ),
    )
    whole_cycles_only: bool = Field(
        True, description="Only integer multiples of the cycle (lambda/2) are applied"
    )
    free_offsets: bool = Field(
        False,
        description=(
            "Also allow non-integer offsets >= 0.3 cycle (removes real island "
            "motion too; not for v0.5)"
        ),
    )
    gnss_veto: bool = Field(
        True, description="A shift must agree in sign with the GNSS residual direction"
    )
    min_region_area: int = Field(20, gt=0, description="Minimum region size in pixels")
    anchor_distance_meters: float = Field(
        12_000.0,
        gt=0,
        description=(
            "A region is measured against anchored coherent land within this "
            "distance of it (edge medians on both sides)"
        ),
    )
    cycle_tolerance: float = Field(
        0.15,
        gt=0,
        lt=0.5,
        description="|jump - round(jump)| must be below this (cycles) to shift",
    )
    min_coherent_area_km2: float = Field(
        0.9,
        gt=0,
        description="Regions with less coherent area than this are not measured",
    )
    min_edge_area_km2: float = Field(
        0.18,
        gt=0,
        description="Minimum coherent area on each side of the water for a jump",
    )
    residual_gate_cycles: float | None = Field(
        None,
        description=(
            "If set, only regions whose median timeseries_inversion_residuals "
            "exceed this (cycles) may be shifted; necessary, not sufficient"
        ),
    )

    model_config = {"validate_assignment": True, "extra": "forbid"}


class CalibrationOptions(BaseModel):
    """Calibration algorithm options.

    Attributes
    ----------
    grid_type : str
        UNR GNSS grid: ``'constant'`` (linear rates) or ``'variable'``
        (per-epoch positions).
    reference_frame : str
        GNSS reference frame (``'IGS14'`` or ``'IGS20'``).
    unwrap_error_correction : bool
        Correct unwrapping-error islands by whole-cycle offsets.
    apply_tropo_correction : bool
        Apply tropospheric correction when `tropo_files` is configured.
    apply_solid_earth_tide_correction : bool
        Remove each product's ``/corrections/solid_earth_tide`` layer before
        the fit. Products without it are calibrated without it.
    recompute_gnss : bool
        Rebuild GNSS LOS caches even if they exist.
    window_size_meters : float
        Fit window size in meters.
    posting_meters : float
        Input pixel spacing in meters.
    event_mask_buffer_pixels : int
        Pixels to dilate event regions by before filling; ``0`` disables.
    residual_outlier_mad_threshold : float or None
        Automatic event detection when an epoch has no event mask: pixels
        whose downsampled ``disp - gnss_los`` residual is more than this
        many scaled MADs from the median are filled before the fit.
        ``None`` disables it.
    residual_region_mad_threshold : float or None
        Automatic detection of wide, coherent deformation (e.g. a subsidence
        bowl) that the per-pixel test misses: pixels above this looser MAD
        threshold that form regions of at least `residual_region_min_pixels`
        are filled. Combined with `residual_outlier_mad_threshold` when both
        are set. ``None`` disables it.
    residual_region_min_pixels : int
        Minimum region size (downsampled pixels) for
        `residual_region_mad_threshold`.
    mask_fit_residual_outliers : bool
        Discard the 15% most extreme residuals at each end (~30% of pixels
        per window, real deformation included) before fitting.
    weight_fit_by_gnss_uncertainty : bool
        Weight the fit by the GNSS LOS uncertainty (Govorcin et al. 2025):
        the rate uncertainty for ``'constant'``, the per-epoch position
        uncertainty for ``'variable'``.
    calibration_surface_smoothing_method : str
        Low-pass filter for the assembled surface: ``'gaussian'``,
        ``'gaussian_fft'``, ``'hanning_fft'`` or ``'savitzky_golay'``.
    calibration_surface_smoothing_sigma : float or None
        Smoothing sigma in pixels (ignored for ``'savitzky_golay'``).
        ``None`` uses window size / 8; ``0`` disables smoothing.
    savitzky_golay : SavitzkyGolayOptions
        Savitzky-Golay filter parameters.

    """

    grid_type: Literal["constant", "variable"] = Field(
        "constant",
        description=(
            "UNR GNSS grid: 'constant' (linear rates) or 'variable' "
            "(per-epoch positions)"
        ),
    )
    reference_frame: str = Field(
        "IGS20", description="GNSS reference frame (IGS14 or IGS20)"
    )
    unwrap_error_correction: bool = Field(
        False,
        description=(
            "Correct unwrapping-error islands by whole-cycle offsets. Off until "
            "trade study TS-U1 passes (PRD R-U1); details in `unwrap`"
        ),
    )
    apply_tropo_correction: bool = Field(
        True,
        description="Apply tropospheric correction when tropo_files is configured",
    )
    apply_solid_earth_tide_correction: bool = Field(
        True,
        description=(
            "Remove each product's /corrections/solid_earth_tide layer before "
            "the fit; products without it are calibrated without it"
        ),
    )
    recompute_gnss: bool = Field(
        True, description="Rebuild GNSS LOS caches even if they exist"
    )
    window_size_meters: float = Field(
        30000.0, gt=0, description="Fit window size in meters"
    )
    posting_meters: float = Field(
        30.0, gt=0, description="Input pixel spacing in meters"
    )
    downsample_factor: int = Field(
        1, ge=1, description="Block-average the inputs by this factor before the fit"
    )
    downsample_method: Literal["mean", "median"] = Field(
        "mean", description="Reducer used when downsampling"
    )
    downsample_weighted: bool = Field(
        False, description="Weight the block average by the GNSS LOS uncertainty"
    )
    event_mask_buffer_pixels: int = Field(
        0,
        ge=0,
        description="Pixels to dilate event regions by before filling; 0 disables",
    )
    residual_outlier_mad_threshold: float | None = Field(
        None,
        gt=0,
        description=(
            "Automatic event detection when an epoch has no event mask: fill "
            "pixels whose downsampled (disp - gnss_los) residual is more than "
            "this many scaled MADs from the median. None disables it"
        ),
    )
    residual_region_mad_threshold: float | None = Field(
        None,
        gt=0,
        description=(
            "Automatic detection of wide, coherent deformation: fill regions of "
            "at least residual_region_min_pixels above this looser MAD "
            "threshold. None disables it"
        ),
    )
    residual_region_min_pixels: int = Field(
        20,
        gt=0,
        description=(
            "Minimum region size (downsampled pixels) for residual_region_mad_threshold"
        ),
    )
    mask_fit_residual_outliers: bool = Field(
        True,
        description=(
            "Discard the 15% most extreme residuals at each end (~30% of "
            "pixels per window, real deformation included) before fitting"
        ),
    )
    weight_fit_by_gnss_uncertainty: bool = Field(
        False,
        description=(
            "Weight the fit by the GNSS LOS uncertainty (Govorcin et al. 2025): "
            "rate uncertainty for 'constant', per-epoch position uncertainty "
            "for 'variable'"
        ),
    )
    calibration_surface_smoothing_method: Literal[
        "gaussian", "gaussian_fft", "hanning_fft", "savitzky_golay"
    ] = Field(
        "gaussian",
        description=(
            "Low-pass filter for the assembled surface: 'gaussian', "
            "'gaussian_fft', 'hanning_fft' or 'savitzky_golay'"
        ),
    )
    calibration_surface_smoothing_sigma: float | None = Field(
        None,
        ge=0,
        description=(
            "Smoothing sigma in pixels (ignored for 'savitzky_golay'); None uses "
            "window size / 8, 0 disables smoothing"
        ),
    )
    savitzky_golay: SavitzkyGolayOptions = Field(
        default_factory=SavitzkyGolayOptions,
        description="Savitzky-Golay filter parameters",
    )
    # Schema v2 (plan T19): the trade-study algorithm, all defaulting to the
    # gamma behaviour so a v1 file loads unchanged.
    surface: SurfaceOptions = Field(
        default_factory=SurfaceOptions, description="Surface estimator (R-S1..R-S3)"
    )
    weights: WeightOptions = Field(
        default_factory=WeightOptions, description="Fit weights (R-S4)"
    )
    tropo: TropoOptions = Field(
        default_factory=TropoOptions, description="Tropospheric correction mode (R-T1)"
    )
    gnss: GnssOptions = Field(
        default_factory=GnssOptions, description="GNSS grid handling (R-G1..R-G5)"
    )
    uncertainty: UncertaintyOptions = Field(
        default_factory=UncertaintyOptions, description="calibration_std model (R-E1)"
    )
    unwrap: UnwrapOptions = Field(
        default_factory=UnwrapOptions,
        description="Unwrap-error correction details (R-U1)",
    )
    frame: FrameOptions = Field(
        default_factory=FrameOptions,
        description=(
            "Frame-level facts from the frame-parameter table (plate, category)"
        ),
    )

    # A typo in a cal-disp algorithm_parameters.yaml must fail loudly, not be
    # ignored (plan T19; previously extra keys were dropped silently).
    model_config = {"validate_assignment": True, "extra": "forbid"}


class DecompositionOptions(BaseModel):
    """Decomposition algorithm options.

    Attributes
    ----------
    inversion_method : str
        Method for least-squares inversion
    uncertainty_propagation : bool
        Whether to compute and propagate uncertainties
    min_geometries : int
        Minimum number of viewing geometries required
    vertical_only : bool
        Only solve for vertical component (ignore horizontal)
    geometry_weighting : str
        How to weight different geometries
    quality_threshold : float
        Minimum quality metric for including pixels

    """

    inversion_method: Literal["weighted_least_squares", "ridge", "lasso"] = Field(
        "weighted_least_squares", description="Method for least-squares inversion"
    )
    uncertainty_propagation: bool = Field(
        True, description="Whether to compute and propagate uncertainties to output"
    )
    min_geometries: int = Field(
        2,
        ge=1,
        description="Minimum number of viewing geometries required for decomposition",
    )
    vertical_only: bool = Field(
        False,
        description=(
            "Only solve for vertical component (assume horizontal motion is zero)"
        ),
    )
    geometry_weighting: Literal["uniform", "temporal_coherence", "uncertainty"] = Field(
        "temporal_coherence",
        description="How to weight different viewing geometries in inversion",
    )
    quality_threshold: float = Field(
        0.5,
        ge=0.0,
        le=1.0,
        description="Minimum quality metric for including pixels in decomposition",
    )


class OutputOptions(BaseModel):
    """Output file options.

    Attributes
    ----------
    gtiff_creation_options : list
        GDAL creation options for GeoTIFF files
    add_overviews : bool
        Whether to add overviews to output GeoTIFFs
    compression : str
        Compression method for output files

    """

    gtiff_creation_options: list[str] = Field(
        default_factory=lambda: [
            "COMPRESS=lzw",
            "ZLEVEL=4",
            "BIGTIFF=yes",
            "TILED=yes",
            "BLOCKXSIZE=128",
            "BLOCKYSIZE=128",
        ],
        description="GDAL creation options for GeoTIFF output files",
    )
    add_overviews: bool = Field(
        False, description="Whether to add overviews (pyramids) to output GeoTIFFs"
    )
    compression: str = Field("lzw", description="Compression method for output files")


class AlgorithmParameters(BaseModel):
    """Algorithm parameters configuration.

    This contains algorithm-specific settings that typically don't change
    between runs. These are saved in algorithm_parameters.yaml.
    """

    schema_version: Literal[1, 2] = Field(
        ALGORITHM_SCHEMA_VERSION,
        description=(
            "Version of this file's schema. 1 = gamma 0.3 (no nested option "
            "groups); 2 adds surface/weights/tropo/gnss/uncertainty/unwrap. A "
            "file without the key is treated as version 1 and upgraded in memory"
        ),
    )
    processing_options: ProcessingOptions = Field(
        default_factory=ProcessingOptions,
        description="Common processing settings used by multiple workflows",
    )
    calibration_options: CalibrationOptions = Field(
        default_factory=CalibrationOptions,
        description="Settings for calibration workflow",
    )
    decomposition_options: DecompositionOptions = Field(
        default_factory=DecompositionOptions,
        description="Settings for decomposition workflow",
    )
    output_options: OutputOptions = Field(
        default_factory=OutputOptions, description="Output file format options"
    )

    model_config = {
        "validate_assignment": True,
        "extra": "forbid",
    }

    @classmethod
    def from_yaml(cls, yaml_path: str | Path) -> AlgorithmParameters:
        """Load algorithm parameters from YAML file."""
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            msg = f"Algorithm parameters file not found: {yaml_path}"
            raise FileNotFoundError(msg)

        with open(yaml_path) as f:
            data = yaml.safe_load(f) or {}

        return cls.from_dict(data)

    @classmethod
    def from_dict(cls, data: dict) -> AlgorithmParameters:
        """Build from a parsed YAML mapping, upgrading version-1 files in memory.

        A version-1 file (gamma 0.3) has no ``schema_version`` and none of the
        nested option groups; every new field takes the default that reproduces
        the gamma behaviour, so the upgrade changes no result.
        """
        data = dict(data)
        version = data.setdefault("schema_version", 1)
        if version not in (1, ALGORITHM_SCHEMA_VERSION):
            msg = (
                f"algorithm_parameters schema_version {version!r} is not supported; "
                f"this Venti reads versions 1 and {ALGORITHM_SCHEMA_VERSION}"
            )
            raise ValueError(msg)
        if version == 1:
            logger.info(
                "algorithm_parameters has schema_version 1 (gamma 0.3); "
                "upgrading to %d in memory with gamma-equivalent defaults",
                ALGORITHM_SCHEMA_VERSION,
            )
            data["schema_version"] = ALGORITHM_SCHEMA_VERSION
        return cls(**data)

    def to_yaml(self, yaml_path: str | Path) -> None:
        """Save algorithm parameters to YAML file."""
        yaml_path = Path(yaml_path)
        data = self.model_dump(mode="python")

        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        with open(yaml_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


# ============================================================================
# Run Configuration (runconfig.yaml)
# ============================================================================


class CalibrationInputGroup(BaseModel):
    """Input file group for calibration workflow.

    Attributes
    ----------
    input_files : Path
        Directory containing input NetCDF displacement files
    los_file : Path
        Path to LOS unit vector file (3-band GeoTIFF: east, north, up)
    water_mask : Path
        Path to water mask file (GeoTIFF)
    custom_mask : Path, optional
        Path to custom mask file (GeoTIFF, 1=valid, 0=invalid)
    frame_bounds : Path, optional
        Path to GeoJSON file defining frame boundaries
    tropo_files : Path, optional
        Directory with tropospheric correction files
    event_mask_dir : Path, optional
        Directory with per-epoch event mask GeoTIFFs (1=valid, 0=event region)
    reference_point : tuple[int, int], optional
        Reference point (row, col), None for auto-select

    """

    input_files: Path = Field(
        ...,
        description=(
            "Directory containing input NetCDF displacement files from OPERA DISP"
            " products"
        ),
    )
    los_file: Path = Field(
        ...,
        description=(
            "Path to LOS unit vector file (3-band GeoTIFF: east, north, up components)"
        ),
    )
    water_mask: Path = Field(
        ...,
        description="Path to water mask file (GeoTIFF, 1=valid land, 0=invalid/water)",
    )
    custom_mask: Path | None = Field(
        None,
        description=(
            "Path to custom mask file (GeoTIFF, 1=valid, 0=invalid). "
            "Combined with water_mask if provided"
        ),
    )
    frame_bounds: Path | None = Field(
        None,
        description="Path to GeoJSON file defining frame boundaries for processing",
    )
    tropo_files: Path | None = Field(
        None,
        description="Directory with tropospheric correction NetCDF files (optional)",
    )
    event_mask_dir: Path | None = Field(
        None,
        description=(
            "Directory containing per-epoch event mask GeoTIFFs (1=valid, 0=event "
            "region). Each mask file must be named with the displacement file stem as "
            "a prefix, e.g. as produced by generate_event_mask.py. When provided, "
            "event-region pixels are filled with nearest valid neighbors before "
            "calibration surface estimation, then the surface is removed from the "
            "full (unmasked) displacement."
        ),
    )
    wavelength_m: float = Field(
        0.05546,
        gt=0,
        description=(
            "Radar wavelength in meters used to convert displacement to phase. "
            "Default is ~0.05546 m (Sentinel-1 C-band)."
        ),
    )
    reference_point: tuple[int, int] | None = Field(
        None,
        description=(
            "Reference point as [row, col] in pixel coordinates, None for automatic"
            " selection"
        ),
    )

    @field_validator("input_files", "los_file", "water_mask", mode="before")
    @classmethod
    def convert_to_path(cls, v):
        """Convert string paths to Path objects."""
        if v is not None:
            return Path(v)
        return v

    @field_validator(
        "custom_mask", "frame_bounds", "tropo_files", "event_mask_dir", mode="before"
    )
    @classmethod
    def convert_optional_to_path(cls, v):
        """Convert optional string paths to Path objects."""
        if v is not None and v != "":
            return Path(v)
        return None if v == "" else v

    @field_validator("input_files", "los_file", "water_mask")
    @classmethod
    def validate_exists(cls, v, info):
        """Validate that required paths exist (skip for template placeholders)."""
        # Skip validation for template placeholder paths
        if str(v).startswith("path/to/"):
            return v
        if not v.exists():
            msg = f"{info.field_name} does not exist: {v}"
            raise ValueError(msg)
        return v

    @field_validator("custom_mask", "frame_bounds", "tropo_files", "event_mask_dir")
    @classmethod
    def validate_optional_exists(cls, v):
        """Validate that optional paths exist if provided.

        Skips validation for template placeholders.
        """
        if v is None:
            return v
        # Skip validation for template placeholder paths
        if str(v).startswith("path/to/"):
            return v
        if not v.exists():
            msg = f"Path does not exist: {v}"
            raise ValueError(msg)
        return v


class DecompositionInputGroup(BaseModel):
    """Input file group for decomposition workflow.

    Attributes
    ----------
    asc_displacement_files : Path
        Directory with ascending geometry displacement files (NetCDF format)
    desc_displacement_files : Path
        Directory with descending geometry displacement files (NetCDF format)
    asc_los_file : Path
        Path to ascending LOS unit vector file (3-band GeoTIFF: east, north, up)
    desc_los_file : Path
        Path to descending LOS unit vector file (3-band GeoTIFF: east, north, up)
    water_mask : Path
        Path to water mask file (GeoTIFF)
    custom_mask : Path, optional
        Path to custom mask file (GeoTIFF, 1=valid, 0=invalid)
    frame_bounds : Path, optional
        Path to GeoJSON file defining frame boundaries
    asc_static_layers : Path, optional
        Path to ascending static layers file (GeoTIFF with temporal coherence, etc.)
    desc_static_layers : Path, optional
        Path to descending static layers file (GeoTIFF with temporal coherence, etc.)
    reference_point : tuple[int, int], optional
        Reference point (row, col), None for auto-select

    """

    asc_displacement_files: Path = Field(
        ...,
        description=(
            "Directory with ascending geometry displacement files in NetCDF format"
        ),
    )
    desc_displacement_files: Path = Field(
        ...,
        description=(
            "Directory with descending geometry displacement files in NetCDF format"
        ),
    )
    asc_los_file: Path = Field(
        ...,
        description=(
            "Path to ascending LOS unit vector file (3-band GeoTIFF: east, north, up)"
        ),
    )
    desc_los_file: Path = Field(
        ...,
        description=(
            "Path to descending LOS unit vector file (3-band GeoTIFF: east, north, up)"
        ),
    )
    water_mask: Path = Field(
        ...,
        description="Path to water mask file (GeoTIFF, 1=valid land, 0=invalid/water)",
    )
    custom_mask: Path | None = Field(
        None,
        description=(
            "Path to custom mask file (GeoTIFF, 1=valid, 0=invalid). "
            "Combined with water_mask if provided"
        ),
    )
    frame_bounds: Path | None = Field(
        None,
        description="Path to GeoJSON file defining frame boundaries for processing",
    )
    asc_static_layers: Path | None = Field(
        None,
        description=(
            "Path to ascending static layers GeoTIFF (temporal coherence, amplitude"
            " dispersion, etc.)"
        ),
    )
    desc_static_layers: Path | None = Field(
        None,
        description=(
            "Path to descending static layers GeoTIFF (temporal coherence, amplitude"
            " dispersion, etc.)"
        ),
    )
    reference_point: tuple[int, int] | None = Field(
        None,
        description=(
            "Reference point as [row, col] in pixel coordinates, None for automatic"
            " selection"
        ),
    )

    @field_validator(
        "asc_displacement_files",
        "desc_displacement_files",
        "asc_los_file",
        "desc_los_file",
        "water_mask",
        mode="before",
    )
    @classmethod
    def convert_to_path(cls, v):
        """Convert string paths to Path objects."""
        if v is not None:
            return Path(v)
        return v

    @field_validator(
        "custom_mask",
        "frame_bounds",
        "asc_static_layers",
        "desc_static_layers",
        mode="before",
    )
    @classmethod
    def convert_optional_to_path(cls, v):
        """Convert optional string paths to Path objects."""
        if v is not None and v != "":
            return Path(v)
        return None if v == "" else v

    @field_validator(
        "asc_displacement_files",
        "desc_displacement_files",
        "asc_los_file",
        "desc_los_file",
        "water_mask",
    )
    @classmethod
    def validate_exists(cls, v, info):
        """Validate that required paths exist (skip for template placeholders)."""
        # Skip validation for template placeholder paths
        if str(v).startswith("path/to/"):
            return v
        if not v.exists():
            msg = f"{info.field_name} does not exist: {v}"
            raise ValueError(msg)
        return v

    @field_validator(
        "custom_mask", "frame_bounds", "asc_static_layers", "desc_static_layers"
    )
    @classmethod
    def validate_optional_exists(cls, v):
        """Validate that optional paths exist if provided.

        Skips validation for template placeholders.
        """
        if v is None:
            return v
        # Skip validation for template placeholder paths
        if str(v).startswith("path/to/"):
            return v
        if not v.exists():
            msg = f"Path does not exist: {v}"
            raise ValueError(msg)
        return v


class ProductPathGroup(BaseModel):
    """Product path group configuration.

    Attributes
    ----------
    product_path : Path
        Directory where products will be placed
    scratch_path : Path
        Scratch directory; during a run all temporary files (Python and GDAL)
        go to ``<scratch_path>/tmp``, so it must be writable
    sas_output_path : Path
        Path to SAS output directory
    product_version : str
        Version of the product in <major>.<minor> format

    """

    product_path: Path = Field(
        Path("output"), description="Directory where products will be placed"
    )
    scratch_path: Path = Field(
        Path("scratch"),
        description=(
            "Scratch directory; during a run all temporary files go to "
            "<scratch_path>/tmp, so it must be writable"
        ),
    )
    sas_output_path: Path = Field(
        Path("output"), description="Path to SAS output directory"
    )
    product_version: str = Field(
        "1.0", description="Version of the product in <major>.<minor> format"
    )

    @field_validator("product_path", "scratch_path", "sas_output_path", mode="before")
    @classmethod
    def convert_to_path(cls, v):
        """Convert string paths to Path objects."""
        if isinstance(v, str):
            return Path(v)
        return v


class WorkerSettings(BaseModel):
    """Worker configuration for processing.

    Attributes
    ----------
    gpu_enabled : bool
        Whether to use GPU for processing (if available)
    threads_per_worker : int
        Number of threads to use per worker (sets OMP_NUM_THREADS)
    block_shape : list[int]
        Size (rows, columns) of blocks of data to load at a time

    """

    gpu_enabled: bool = Field(
        False, description="Whether to use GPU for processing (if available)"
    )
    threads_per_worker: int = Field(
        1,
        ge=1,
        description="Number of threads to use per worker (sets OMP_NUM_THREADS)",
    )
    block_shape: list[int] = Field(
        [512, 512],
        min_length=2,
        max_length=2,
        description="Size (rows, columns) of blocks of data to load at a time",
    )

    @field_validator("block_shape")
    @classmethod
    def validate_block_shape(cls, v):
        """Validate block shape dimensions."""
        if len(v) != 2:
            msg = "block_shape must have exactly 2 elements"
            raise ValueError(msg)
        if any(dim <= 0 for dim in v):
            msg = "block_shape dimensions must be positive"
            raise ValueError(msg)
        return v


class PrimaryExecutable(BaseModel):
    """Primary executable configuration.

    Attributes
    ----------
    product_type : str
        Product type of the workflow
    workflow_name : str
        Name of the workflow to execute

    """

    product_type: Literal["CAL", "VLM"] = Field(
        "CAL", description="Product type of the workflow"
    )
    workflow_name: Literal["calibrate", "decompose"] = Field(
        "calibrate", description="Name of the workflow to execute"
    )


class RunConfig(BaseModel):
    """Run configuration.

    This contains run-specific settings that change between runs.
    These are saved in runconfig.yaml.

    Note: Provide either calibration_input_group OR decomposition_input_group
    based on the workflow type.
    """

    calibration_input_group: CalibrationInputGroup | None = Field(
        None,
        description=(
            "Input files for calibration workflow (required if workflow_name is"
            " 'calibrate')"
        ),
    )
    decomposition_input_group: DecompositionInputGroup | None = Field(
        None,
        description=(
            "Input files for decomposition workflow (required if workflow_name is"
            " 'decompose')"
        ),
    )
    product_path_group: ProductPathGroup = Field(
        default_factory=ProductPathGroup, description="Output product paths and version"
    )
    primary_executable: PrimaryExecutable = Field(
        default_factory=PrimaryExecutable,
        description="Primary executable configuration",
    )
    worker_settings: WorkerSettings = Field(
        default_factory=WorkerSettings,
        description="Worker configuration for processing",
    )
    log_file: str | None = Field(
        None, description="Path to output log file (in addition to logging to stderr)"
    )
    keep_paths_relative: bool = Field(
        False,
        description="Don't resolve filepaths that are given as relative to be absolute",
    )

    model_config = {
        "validate_assignment": True,
        "extra": "forbid",
    }

    def model_post_init(self, __context: Any, /) -> None:
        """Validate that the correct input group is provided for the workflow type."""
        super().model_post_init(__context)
        workflow_name = self.primary_executable.workflow_name

        if workflow_name == "calibrate" and self.calibration_input_group is None:
            msg = (
                "calibration_input_group is required when workflow_name is"
                " 'calibrate'. Fill in the calibration_input_group section with"
                " your data paths."
            )
            raise ValueError(msg)
        elif workflow_name == "decompose" and self.decomposition_input_group is None:
            msg = (
                "decomposition_input_group is required when workflow_name is"
                " 'decompose'. Fill in the decomposition_input_group section with"
                " your data paths."
            )
            raise ValueError(msg)

    @property
    def input_file_group(self):
        """Get the active input file group based on workflow type."""
        if self.primary_executable.workflow_name == "calibrate":
            return self.calibration_input_group
        else:
            return self.decomposition_input_group

    @classmethod
    def from_yaml(cls, yaml_path: str | Path) -> RunConfig:
        """Load run configuration from YAML file."""
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            msg = f"Run configuration file not found: {yaml_path}"
            raise FileNotFoundError(msg)

        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        return cls(**data)

    def to_yaml(self, yaml_path: str | Path) -> None:
        """Save run configuration to YAML file."""
        yaml_path = Path(yaml_path)
        data = self.model_dump(mode="python")

        # Convert Path objects to strings for YAML serialization
        def convert_paths(obj):
            if isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            elif isinstance(obj, Path):
                return str(obj)
            return obj

        data = convert_paths(data)

        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        with open(yaml_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)


# ============================================================================
# Combined Configuration
# ============================================================================


class VentiConfig(BaseModel):
    """Combined Venti configuration.

    This combines both run configuration and algorithm parameters.
    """

    run_config: RunConfig
    algorithm_parameters: AlgorithmParameters

    @property
    def input_options(self):
        """Shortcut to the active input group (calibration, else decomposition)."""
        return (
            self.run_config.calibration_input_group
            or self.run_config.decomposition_input_group
        )

    @property
    def worker_settings(self):
        """Shortcut to ``run_config.worker_settings``."""
        return self.run_config.worker_settings

    @property
    def grid_settings(self):
        """Flat view of calibration and processing options used by the workflow."""

        # Create a dynamic object that combines calibration_options
        # and processing_options
        class GridSettings:
            def __init__(self, config):
                self.config = config

            @property
            def grid_type(self):
                return self.config.algorithm_parameters.calibration_options.grid_type

            @property
            def reference_frame(self):
                return (
                    self.config.algorithm_parameters.calibration_options.reference_frame
                )

            @property
            def downsample_factor(self):
                proc_opts = self.config.algorithm_parameters.processing_options
                return proc_opts.cal_downsample_factor

            @property
            def downsample_method(self):
                proc_opts = self.config.algorithm_parameters.processing_options
                return proc_opts.downsample_method

            @property
            def downsample_weighted(self):
                proc_opts = self.config.algorithm_parameters.processing_options
                return proc_opts.downsample_weighted

            @property
            def window_size_meters(self):
                cal_opts = self.config.algorithm_parameters.calibration_options
                return cal_opts.window_size_meters

            @property
            def posting_meters(self):
                cal_opts = self.config.algorithm_parameters.calibration_options
                return cal_opts.posting_meters

            @property
            def recompute_gnss(self):
                cal_opts = self.config.algorithm_parameters.calibration_options
                return cal_opts.recompute_gnss

            @property
            def apply_tropo_correction(self):
                cal_opts = self.config.algorithm_parameters.calibration_options
                return cal_opts.apply_tropo_correction

            @property
            def output_posting_meters(self):
                proc_opts = self.config.algorithm_parameters.processing_options
                return proc_opts.vlm_output_posting_meters

        return GridSettings(self)

    @classmethod
    def from_yaml_files(
        cls,
        runconfig_path: str | Path,
        algorithm_params_path: str | Path,
    ) -> VentiConfig:
        """Load configuration from separate YAML files.

        Parameters
        ----------
        runconfig_path : str or Path
            Path to runconfig.yaml file
        algorithm_params_path : str or Path
            Path to algorithm_parameters.yaml file

        Returns
        -------
        VentiConfig
            Combined configuration

        """
        run_config = RunConfig.from_yaml(runconfig_path)
        algorithm_params = AlgorithmParameters.from_yaml(algorithm_params_path)

        return cls(
            run_config=run_config,
            algorithm_parameters=algorithm_params,
        )


# ============================================================================
# Helper Functions
# ============================================================================


def load_config(
    runconfig_path: str | Path,
    algorithm_params_path: str | Path | None = None,
) -> VentiConfig:
    """Load Venti configuration from YAML files.

    Parameters
    ----------
    runconfig_path : str or Path
        Path to runconfig.yaml file
    algorithm_params_path : str or Path, optional
        Path to algorithm_parameters.yaml. Defaults to the file of that name
        next to `runconfig_path`.

    Returns
    -------
    VentiConfig
        Combined configuration

    """
    runconfig_path = Path(runconfig_path)

    if algorithm_params_path is None:
        algorithm_params_path = runconfig_path.parent / "algorithm_parameters.yaml"
    return VentiConfig.from_yaml_files(runconfig_path, algorithm_params_path)


def _write_yaml_with_comments(
    data: dict, model: type[BaseModel], output_path: Path
) -> None:
    """Write YAML with inline comments from Field descriptions.

    Parameters
    ----------
    data : dict
        Data to write
    model : type[BaseModel]
        Pydantic model class to extract descriptions from
    output_path : Path
        Path to write YAML file

    """

    def write_dict_with_comments(
        obj: dict, model_class: type[BaseModel], indent: int = 0
    ) -> str:
        """Recursively write dict with comments."""
        lines = []
        indent_str = "  " * indent

        for key, value in obj.items():
            # Get field description if available
            if hasattr(model_class, "model_fields") and key in model_class.model_fields:
                field_info = model_class.model_fields[key]
                description = field_info.description
                if description:
                    # Write comment above field
                    lines.append(f"{indent_str}# {description}")

            # Handle nested dicts (nested models)
            if isinstance(value, dict):
                lines.append(f"{indent_str}{key}:")
                # Get nested model class if available
                if (
                    hasattr(model_class, "model_fields")
                    and key in model_class.model_fields
                ):
                    field_info = model_class.model_fields[key]
                    # Get the actual type (handle Optional, etc.)
                    field_type = field_info.annotation
                    if hasattr(field_type, "__origin__"):  # Handle Optional, Union
                        args = getattr(field_type, "__args__", ())
                        field_type = args[0] if args else field_type
                    if isinstance(field_type, type) and issubclass(
                        field_type, BaseModel
                    ):
                        lines.append(
                            write_dict_with_comments(value, field_type, indent + 1)
                        )
                    else:
                        lines.append(
                            write_dict_with_comments(value, model_class, indent + 1)
                        )
                else:
                    lines.append(
                        write_dict_with_comments(value, model_class, indent + 1)
                    )
            # Handle lists
            elif isinstance(value, list):
                lines.append(f"{indent_str}{key}:")
                for item in value:
                    if isinstance(item, dict):
                        lines.append(f"{indent_str}  -")
                        lines.append(
                            write_dict_with_comments(item, model_class, indent + 2)
                        )
                    else:
                        lines.append(f"{indent_str}  - {item}")
            # Handle None
            elif value is None:
                lines.append(f"{indent_str}{key}: null")
            # Handle strings with special characters
            elif isinstance(value, str):
                if any(c in value for c in [":", "#", "@", "`"]):
                    lines.append(f"{indent_str}{key}: '{value}'")
                else:
                    lines.append(f"{indent_str}{key}: {value}")
            # Handle Path objects
            elif isinstance(value, Path):
                lines.append(f"{indent_str}{key}: {value!s}")
            # Handle everything else
            else:
                lines.append(f"{indent_str}{key}: {value}")

        return "\n".join(lines)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    content = write_dict_with_comments(data, model)

    with open(output_path, "w") as f:
        f.write(content)


def create_config_templates(output_dir: str | Path = ".") -> tuple[Path, Path]:
    """Create template configuration files.

    Generates two template files:
    1. runconfig.yaml - Contains both calibration and decomposition input groups
    2. algorithm_parameters.yaml - Algorithm settings (shared)

    The runconfig.yaml includes both input groups - users only fill in the one
    they need based on their workflow type (calibrate or decompose).

    Parameters
    ----------
    output_dir : str or Path
        Directory to write template files

    Returns
    -------
    tuple of Path
        (runconfig_path, algorithm_params_path)

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    runconfig_path = output_dir / "runconfig.yaml"
    algorithm_params_path = output_dir / "algorithm_parameters.yaml"

    # Create config with BOTH input groups as examples
    # Users will fill in only the one they need based on workflow_name
    config = RunConfig(
        calibration_input_group=CalibrationInputGroup(
            input_files=Path("path/to/displacement/files"),
            los_file=Path("path/to/los_vectors.tif"),
            water_mask=Path("path/to/water_mask.tif"),
            custom_mask=None,
            frame_bounds=None,
            tropo_files=None,
            reference_point=None,
        ),
        decomposition_input_group=DecompositionInputGroup(
            asc_displacement_files=Path("path/to/asc_displacement/"),
            desc_displacement_files=Path("path/to/desc_displacement/"),
            asc_los_file=Path("path/to/asc_los_vectors.tif"),
            desc_los_file=Path("path/to/desc_los_vectors.tif"),
            water_mask=Path("path/to/water_mask.tif"),
            custom_mask=None,
            frame_bounds=None,
            asc_static_layers=None,
            desc_static_layers=None,
            reference_point=None,
        ),
        primary_executable=PrimaryExecutable(
            product_type="CAL",
            workflow_name="calibrate",
        ),
    )

    algorithm_params = AlgorithmParameters()

    # Convert configs to dicts
    config_data = config.model_dump(mode="python")
    algorithm_params_data = algorithm_params.model_dump(mode="python")

    # Convert Path objects to strings
    def convert_paths(obj):
        if isinstance(obj, dict):
            return {k: convert_paths(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_paths(item) for item in obj]
        elif isinstance(obj, Path):
            return str(obj)
        return obj

    config_data = convert_paths(config_data)
    algorithm_params_data = convert_paths(algorithm_params_data)

    # Write YAML files with comments
    _write_yaml_with_comments(config_data, RunConfig, runconfig_path)
    _write_yaml_with_comments(
        algorithm_params_data, AlgorithmParameters, algorithm_params_path
    )

    return runconfig_path, algorithm_params_path


# Shorter name used throughout the workflow API.
WorkflowConfig = VentiConfig
