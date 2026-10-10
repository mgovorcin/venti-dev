# Algorithm parameters

`algorithm_parameters.yaml` is the only place new science options go (the
cal-disp CLI and runconfig are frozen, PRD §2.2). Every key has a default that
reproduces the gamma 0.3 result, so a version-1 file loads unchanged.

| Key | Default | Unit | Requirement | Since |
|---|---|---|---|---|
| `schema_version` | 2 | — | — | v2 (absent = 1, upgraded in memory) |
| `calibration_options.grid_type` | `constant` | — | R-G1 | 1 |
| `calibration_options.reference_frame` | `IGS20` | — | A6 | 1 |
| `calibration_options.unwrap_error_correction` | `false` | — | R-U1 | 1 (default was `true`; off since v2) |
| `calibration_options.apply_tropo_correction` | `true` | — | R-T1 | 1 |
| `calibration_options.apply_solid_earth_tide_correction` | `true` | — | R-T2 | 1 |
| `calibration_options.recompute_gnss` | `true` | — | — | 1 |
| `calibration_options.window_size_meters` | 30000 | m | gamma plane window | 1 |
| `calibration_options.posting_meters` | 30 | m | — | 1 |
| `calibration_options.downsample_factor` | 1 | — | — | 2 (gamma file: 6) |
| `calibration_options.downsample_method` | `mean` | — | — | 2 |
| `calibration_options.downsample_weighted` | `false` | — | — | 2 |
| `calibration_options.event_mask_buffer_pixels` | 0 | px | A5 | 1 |
| `calibration_options.residual_outlier_mad_threshold` | `null` | MAD | — | 1 |
| `calibration_options.residual_region_mad_threshold` | `null` | MAD | — | 1 |
| `calibration_options.residual_region_min_pixels` | 20 | px | — | 1 |
| `calibration_options.mask_fit_residual_outliers` | `true` | — | R-S4 (the global quantile mask; `false` for v0.5) | 1 |
| `calibration_options.weight_fit_by_gnss_uncertainty` | `false` | — | — | 1 |
| `calibration_options.calibration_surface_smoothing_method` | `gaussian` | — | — | 1 |
| `calibration_options.calibration_surface_smoothing_sigma` | `null` | px | — | 1 |
| `calibration_options.savitzky_golay.window_length` | 51 | px | — | 1 |
| `calibration_options.savitzky_golay.polyorder` | 3 | — | — | 1 |
| `calibration_options.surface.method` | `windowed_plane` | — | R-S1/R-S2 (`loclin` for v0.5) | 2 |
| `calibration_options.surface.cutoff_wavelength_meters` | 50000 | m | R-S2 | 2 |
| `calibration_options.surface.fill_gaps` | `false` | — | R-S3 (`true` for v0.5) | 2 |
| `calibration_options.surface.two_pass` | `false` | — | R-S1 (`true` for v0.5) | 2 |
| `calibration_options.weights.coherence_power` | 0 | — | R-S4 (8 for v0.5) | 2 |
| `calibration_options.weights.robust` | `false` | — | R-S4 | 2 |
| `calibration_options.weights.filled_pixel_weight` | 0.02 | — | R-S4 | 2 |
| `calibration_options.tropo.mode` | `legacy` | — | R-T1 (`auto` for v0.5) | 2 |
| `calibration_options.tropo.relief_off_meters` | 300 | m | R-T1 | 2 |
| `calibration_options.tropo.relief_stratified_meters` | 1500 | m | R-T1 | 2 |
| `calibration_options.gnss.buffer_meters` | 0 | m | R-G2 (50000 for v0.5) | 2 |
| `calibration_options.gnss.exclude_defo_nodes` | `false` | — | R-G5 | 2 |
| `calibration_options.gnss.reinterpolate_excluded` | `false` | — | R-G5 | 2 |
| `calibration_options.gnss.reprocessing` | `false` | — | R-G1 (`variable` grid only here) | 2 |
| `calibration_options.gnss.snapshot_id` | `null` | — | R-G4 | 2 |
| `calibration_options.uncertainty.k_grid` | 1.0 | — | R-E1 (`frame_table` for v0.5) | 2 |
| `calibration_options.uncertainty.inflate_inside_areas` | `true` | — | R-E1 | 2 |
| `calibration_options.uncertainty.sigma_disp_placeholder_mm` | 10 | mm | R-E2 | 2 |
| `calibration_options.unwrap.region_source` | `water_mask` | — | R-U1 | 2 |
| `calibration_options.unwrap.whole_cycles_only` | `true` | — | R-U1 | 2 |
| `calibration_options.unwrap.free_offsets` | `false` | — | R-U1 (never for v0.5) | 2 |
| `calibration_options.unwrap.gnss_veto` | `true` | — | R-U1 | 2 |
| `calibration_options.unwrap.min_region_area` | 20 | px | R-U1 | 2 |

Groups `processing_options`, `decomposition_options` and `output_options` are
unchanged from version 1 (see `venti.workflow.config`).

## Rules

- **Additive only.** New keys get defaults equal to the previous behaviour;
  an existing key never changes meaning. Changing a default is a value change
  and follows the golden policy (`CONTRIBUTING.md`).
- **Unknown keys are errors** (`extra = "forbid"`), so a typo in a cal-disp
  file is caught at load time instead of silently falling back to a default.
- **Versioning.** `schema_version` is written by `to_yaml`; `from_yaml`
  accepts 1 (gamma 0.3, no nested groups, no `schema_version`) and 2, and
  upgrades 1 in memory. Reading a newer version than the library knows is an
  error.
- The "v0.5" column values are what the CalVal release (plan T49) will flip
  the shipped `algorithm_parameters.yaml` to, together with a deliberate golden
  regeneration; the schema defaults stay at gamma so the 0.3 line remains
  reproducible.
