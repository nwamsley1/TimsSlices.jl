# TimsSlices.jl

Converts Bruker timsTOF diaPASEF `.d` bundles into ion-mobility-smoothed, m/z-centroided slices stored in
Bruker's own block layout (`.tdfs`, see `docs/format.md`), and expands them back into Pioneer's slice Arrow
schema. The smoothing is paid once at conversion; a search decodes frames on the fly.

```
julia -t 12 --project=. bin/tims_convert.jl <run.d> <out_dir> [options]
julia --project=. bin/tims_convert.jl expand <name.tdfs> <out.arrow>
julia --project=. bin/tims_convert.jl            # prints the options
```

The ion-mobility scale is set in 1/K0 and converted to scans per run with the run's own 1/K0-per-scan slope:
slices every `stride_k0` = 0.0065 1/K0 (rounded up to whole scans) with an IM Gaussian of `im_sigma_k0` = 0.004325
1/K0, i.e. 8 scans and sigma 5 scans on the timsTOF Ultra ramps (1/K0 0.64-1.45 over ~940 scans), 7 scans and
sigma 4 on a 0.60-1.60 timsTOF Pro ramp. `stride` / `im_sigma` (scans) override the conversion.

Pipeline per frame: decode the raw block (zstd, byte planes, per-scan prefix sums) → for every window and every
`stride`-th IM scan, accumulate a Gaussian (`im_sigma` scans) over the scans in reach per TOF bin → Gaussian
(`mz_sigma` bins) along the TOF axis over runs of nearby bins → local maxima with a footprint walk
(`max_half`), intensity-weighted mean position (`wmean`) or Gaussian apex (`gauss`), footprint sum as intensity
→ culls (`min_scans` persistence; `max_peaks` keeps the 1,500 most intense centroids per MS2 slice, MS1
uncapped) → fixed-point bins (`bin_scale`),
delta coding, byte planes, zstd.

Development:

```
julia --project=. -e 'using Pkg; Pkg.test()'                 # synthetic + HeLa real-data tests (~40 s)
TIMSSLICES_TEST_BIG=1 julia --project=. -e 'using Pkg; Pkg.test()'   # + E. coli / 250 pg equivalence
julia --project=bench bench/stages.jl <run.d> [n_frames] [--options]   # per-stage ns/peak
bench/wholefile.sh <run.d> <out_dir> [--options]                        # wall time, peak RSS
```

`test/reference/` holds the frozen prototype (`~/BrukerTims/proto/tdf_centroid_to_arrow.jl` of 2026-09-18); the
equivalence test requires identical bin positions and intensities within Float32 rounding.
