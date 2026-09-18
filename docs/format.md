# The `.tdfs` container

One directory per converted timsTOF run. Everything a reader needs is inside; the original `.d` is not consulted.

```
<name>.tdfs/
  meta.json      format_version, source, calibration, converter parameters, counts (see below)
  frames.arrow   one row per converted frame (MS1 and diaPASEF MS2 frames, in acquisition order)
  slices.arrow   one row per slice (= one Pioneer scan), in frame order, then window, then IM scan
  blocks.bin     one block per frame, in the same order as frames.arrow
```

## blocks.bin

```
block = u32 block_size (bytes, header included) | u32 n_words | zstd frame of 4·n_words bytes
```

The zstd payload decompresses to four **byte planes** of `n_words` bytes each (plane `p` holds byte `p` of every
word, little-endian), exactly as Bruker's `analysis.tdf_bin` codec 2. Reassembled, the `UInt32` word stream is:

```
word 1                 n_slices
words 2 .. n_slices+1  peak count of slice j
then, per slice, per peak:
    bin_delta          fixed-point bin minus the previous peak's fixed-point bin in the slice; the accumulator
                       starts at 0xFFFFFFFF, so the first peak's delta is bin + 1 (Bruker's convention)
    intensity          UInt32
```

`n_words = 1 + n_slices + 2 · n_peaks`. Peaks in a slice are sorted by fixed-point bin and unique (centroids
that round to the same bin are merged by summing intensities at encode time); a peak whose rounded intensity is 0
is dropped. A frame with no surviving centroids still has a block (`n_words = 1`).

- fixed-point bin = `round(bin_scale · bin_position)`; `bin_scale` is in `meta.json` (1 = integer TOF bins).
- m/z = `(a + b · bin_position)²` with `a = mz_cal_sqrt_intercept`, `b = mz_cal_sqrt_slope` from `meta.json`.
- intensity = `round(int_scale · centroid_intensity)`; `int_scale` in `meta.json` (default 1). With `sum_scale`
  (default) a centroid's intensity is the ion's summed intensity over the `stride` scans the slice stands for.

## frames.arrow

| column | type | meaning |
|---|---|---|
| frame_id | Int32 | `Frames.Id` of the source |
| ms_order | UInt8 | 1 = MS1, 2 = diaPASEF MS2 |
| cycle_idx | Int32 | increments at every MS1 frame (0 before the first) |
| window_group | UInt8 | `DiaFrameMsMsInfo.WindowGroup` (0 for MS1) |
| rt_s | Float64 | `Frames.Time` in seconds (start of the ramp) |
| n_scans | Int32 | `Frames.NumScans` |
| n_slices, n_peaks | Int32, Int64 | slices and peaks in the block |
| block_offset, block_size, n_words | Int64 | the block in `blocks.bin` |
| first_slice | Int64 | 1-based row of the frame's first slice in `slices.arrow` |

## slices.arrow

| column | type | meaning |
|---|---|---|
| frame_row | Int32 | row in `frames.arrow` |
| slice_in_frame | Int32 | 1-based index of the slice in its frame's block |
| frame_id, ms_order, cycle_idx, window_group | | copied from the frame row |
| im_scan | UInt16 | IM scan the slice is centred on (0-based) |
| window | Int32 | window index within the frame (1-based; 1 for MS1) |
| retention_time | Float32 | minutes, at the slice scan: `(rt_s + im_scan · RampTime / NumScans) / 60` |
| center_mz, isolation_width | Float32 | the quad window (NaN for MS1) |
| collision_energy_ev | Float32 | the CE ramp evaluated at the slice scan (0 for MS1) |
| window_ce | Float32 | the window table's CE (NaN for MS1) |
| tic | Float32 | sum of the slice's stored intensities |
| n_peaks | Int32 | peaks in the slice |
| peak_offset | Int64 | 1-based start of the slice's peaks in the frame's decoded peak arrays |

1/K0 of a slice: `im_scan0_1overK0 + im_slope_1overK0_per_scan · im_scan` (`meta.json`); Pioneer fits its own line.

## meta.json

`format_version` (1), `source`, `source_bin_bytes`, `instrument`, `mz_cal_sqrt_intercept`, `mz_cal_sqrt_slope`,
`im_scan0_1overK0`, `im_slope_1overK0_per_scan`, `ce_ev_intercept`, `ce_ev_slope_per_scan`, `NumScans`, `n_bins`,
`mz_lo`, `mz_hi`, `OneOverK0AcqRangeLower/Upper`, `params` (every `ConvertParams` field), `cull_thr_ms1`,
`cull_thr_ms2` (intensity thresholds actually applied), `bin_scale`, `int_scale`, `zstd_level`, `converter`,
`converted_at`, `n_frames`, `n_slices`, `n_peaks`, `blocks_bytes`.

## Variant A

With `centroid = none` the m/z kernel and centroiding are skipped: each slice holds the IM-accumulated per-bin
sums (integer bins, `bin_scale` still applies) and `params.centroid` in `meta.json` says so. A reader that wants
centroids applies the m/z kernel and the centroider itself.
