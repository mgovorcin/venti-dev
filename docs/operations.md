# Operations: the frozen UNR grid snapshot

PRD R-G4 and R-O4 (plan T48). Operations never read the live UNR server: every
product is calibrated against a **frozen snapshot** of the UNR gridded GNSS
product, identified in the runconfig ancillary paths and in the product
metadata (`gnss_provenance.snapshot_id`, `lookup_sha256`, `digest`). The same
inputs, snapshot and configuration give the same product to 1e-6.

## Layout

```
unr_grid_<version>_<YYYYMMDD>/
    snapshot.json            id, UNR version, frame, grid types, node counts,
                             data span, lookup hash, source URLs, tool versions
    MANIFEST.sha256          one line per file; `sha256sum -c MANIFEST.sha256`
    grid_latlon_lookup.txt   the UNR lookup, byte for byte
    nodes/<id:06d>_<frame>_<grid_type>.tenv8
```

`nodes/` is the layout `venti.gnss.sampling.GnssGridConfig` reads;
`venti.gnss.snapshot.grid_config_from_snapshot(dir, utm_epsg, ...)` builds the
config and carries the snapshot id into the GNSS provenance of every product.

## Making a snapshot

```bash
# whole IGS20 constant grid (operations)
python scripts/snapshot_unr_grid.py /data/unr_snapshots --notes "v0.5 release candidate"

# constant + variable, CONUS only (reprocessing / trade studies)
python scripts/snapshot_unr_grid.py /data/unr_snapshots \
    --grid-types constant variable --bounds 24 50 -125 -66

# verify against the manifest
python scripts/snapshot_unr_grid.py --verify /data/unr_snapshots/unr_grid_0.3_20261006
```

Retrieval goes through geepers (`UnrGridSource`, plan D10). A node that fails
to download is retried alone and, if it still fails, listed in
`snapshot.json: failed_ids`; a snapshot with failures is not released.
Snapshots are immutable: the script refuses to write into an existing
directory, and `verify_snapshot` reports changed, missing and unlisted files.

## Storage

The release snapshot is mirrored to controlled storage (S3 bucket, PRD §7.5
Q6: **bucket and IAM to be decided**, plan T48.2). Until then the snapshot
lives on the processing host and its id and hashes are what the product
records. The bucket layout will mirror the directory above, one prefix per
snapshot id, read-only for the processing role.

## Roll-forward procedure (plan T48.4)

A new UNR release (or a refreshed download of the same version) is adopted
deliberately, about every six months, never automatically:

1. **Snapshot.** Run the script on the new UNR version/date; keep the old
   snapshot in place. Record `snapshot_id`, `lookup_sha256` and the node
   counts in the release notes.
2. **Gate.** Run the e2e validation (plan T46: the 8 benchmark frames and
   the edge cases) with the new snapshot and the frozen algorithm
   parameters. The gate is the same as for an algorithm change: the
   double-difference sill must not rise and the MIDAS velocity comparison
   must not degrade on any frame.
3. **Frame table.** If TS-G1 (plan T39) shows the grid-sigma inflation `k`
   or a plate assignment changed, update `frame_parameters.json` and bump its
   `version`.
4. **Release.** A snapshot change is a **minor** version of cal-disp
   (products are not bit-identical to the previous snapshot's); the runconfig
   ancillary path points at the new snapshot id, the Docker image pins it,
   and the golden dataset is regenerated with the golden-stability rule
   (`TODO.md`, section Algorithm).
5. **Retire.** The previous snapshot stays available for reprocessing and
   for reproducing earlier products; it is never deleted from controlled
   storage.

What is *not* a roll-forward: re-running the script on the same UNR version
and date gives the same files and hashes (the UNR constant product is static
per version); a differing hash means UNR changed the files in place, which is
treated as a new snapshot and goes through the gate.

## Where this moves

`snapshot_unr_grid.py` lives in Venti because Venti owns the sampler that
reads the layout. Plan T48.3 wires cal-disp's `download unr` to read from a
snapshot path (`unr_timeseries_dir`) and record `gnss_snapshot_id`; the
operational copy of this page then belongs to cal-disp's documentation.
