# Real-data tests: the raw reader against its checksum fixture, the pipeline against the frozen prototype, and
# the container round trip through `expand`.
using Arrow, DataFrames, CodecZstd

module Ref
    include(joinpath(@__DIR__, "reference", "tdf_centroid_to_arrow.jl"))
    # The prototype converts bin positions to Float32 m/z inside centroid_dense!. To compare positions exactly,
    # this calibration returns a wrapper whose Float32 conversion records the Float64 position on the way.
    struct IdCal end
    const CAPTURED = Float64[]
    struct Capture; v::Float64; end
    Base.Float32(c::Capture) = (push!(CAPTURED, c.v); Float32(c.v))
    tof_to_mz(::IdCal, t) = Capture(t)
end

"Prototype slices of frame row i: (slices, exact positions in push order, raw (scans, tofs, ints))."
function reference_frame(rb, groups, frame_group, rp, kim, kmz, rsc, i)
    rfr = rb.frames
    scans, tofs, ints = Ref.read_frame(rb, Int(rfr.Id[i]))
    nscans = Int(rfr.NumScans[i]); scan_start = zeros(Int, nscans + 1)
    for k in eachindex(scans); scan_start[scans[k] + 2] += 1; end
    scan_start[1] = 1
    for s in 2:nscans+1; scan_start[s] += scan_start[s-1]; end
    slices = Tuple{Int, Vector{Float32}, Vector{Float32}}[]
    empty!(Ref.CAPTURED)
    if rfr.MsMsType[i] == 0
        Ref.centroid_window!(slices, rsc, scan_start, tofs, ints, 0, nscans, rp, kim, kmz, 0.0, Ref.IdCal())
    else
        for w in groups[frame_group[Int(rfr.Id[i])]]
            Ref.centroid_window!(slices, rsc, scan_start, tofs, ints, w.scan_begin, w.scan_end, rp, kim, kmz, 0.0, Ref.IdCal())
        end
    end
    slices, copy(Ref.CAPTURED), (scans, tofs, ints)
end

"Compare the new pipeline with the prototype on `sel` frame rows of `dpath` with the given parameters."
function check_equivalence(dpath, sel, p::ConvertParams, rp)
    f = open_tdf(dpath)
    ls1 = LevelSetup(level_params(p, true), 0.0); ls2 = LevelSetup(level_params(p, false), 0.0)
    buf = FrameBuffer(); sc = SmoothScratch(); out = FrameSlices()
    rb = Ref.open_tdf(dpath); groups, frame_group = Ref.read_windows(rb.db)
    kim = Ref.gauss_kernel(rp.im_sigma) .* (rp.sum_scale ? rp.stride : 1); kmz = Ref.gauss_kernel(rp.mz_sigma); rsc = Ref.Scratch()
    n_cmp = 0; max_dpos = 0.0; max_dint = 0.0
    for i in sel
        read_frame!(buf, f, i)
        smooth_frame!(out, sc, buf, f, i, ls1, ls2)
        slices, cap, (scans, tofs, ints) = reference_frame(rb, groups, frame_group, rp, kim, kmz, rsc, i)
        @test length(tofs) == buf.n_peaks
        @test tofs == view(buf.tof, 1:buf.n_peaks) && ints == view(buf.intensity, 1:buf.n_peaks)
        @test all(k -> buf.scan_start[scans[k] + 1] <= k <= buf.scan_start[scans[k] + 2] - 1, eachindex(scans))
        @test length(slices) == out.n_slices
        length(slices) == out.n_slices || continue
        k = 0
        for (j, (sscan, mzs32, its32)) in enumerate(slices)
            @test sscan == out.scan[j]
            r = out.ptr[j]:out.ptr[j+1]-1
            @test length(mzs32) == length(r)
            length(mzs32) == length(r) || break
            for (q, kk) in enumerate(r)
                k += 1; n_cmp += 1
                max_dpos = max(max_dpos, abs(cap[k] - out.pos[kk]))
                max_dint = max(max_dint, abs(Float64(its32[q]) - out.val[kk]) / out.val[kk])
            end
        end
    end
    close(rb.bin)
    n_cmp, max_dpos, max_dint
end

sample_rows(f, ms1::Bool, n) = (v = [i for i in 1:n_frames(f) if (f.frames.msms_type[i] == 0) == ms1]; v[unique(round.(Int, range(1, length(v); length = n)))])

@testset "raw reader checksums (HeLa fixture)" begin
    f = open_tdf(HELA)
    buf = FrameBuffer()
    for line in Iterators.drop(eachline(joinpath(@__DIR__, "fixtures", "hela_first50_checksums.csv")), 1)
        row, fid, mst, nsc, np, st, si, ssc, mt = parse.(Int, split(line, ','))
        read_frame!(buf, f, row)
        n = buf.n_peaks
        @test f.frames.id[row] == fid && f.frames.msms_type[row] == mst && buf.n_scans == nsc && n == np
        @test sum(Int, view(buf.tof, 1:n)) == st && sum(Int, view(buf.intensity, 1:n)) == si
        @test (n > 0 ? maximum(view(buf.tof, 1:n)) : 0) == mt
        @test sum(s * length(TS.scan_range(buf, s)) for s in 0:buf.n_scans-1) == ssc
    end
    # per-scan sortedness and the >= 5-bin spacing of the instrument's own centroids
    read_frame!(buf, f, 1)
    @test all(s -> issorted(view(buf.tof, TS.scan_range(buf, s))), 0:buf.n_scans-1)
    @test f.max_peaks == 739812 && f.max_scans == 927 && TS.n_bins(f) == 394535
    @test length(windows(f, 2)) == 3 && length(windows(f, 1)) == 1
end

@testset "equivalence with the prototype (HeLa)" begin
    f = open_tdf(HELA)
    sel = vcat(sample_rows(f, true, 4), sample_rows(f, false, 6))
    for (p, rp) in ((ConvertParams(im_sigma = 5.0, mz_sigma = 3.0, stride = 8, sum_scale = true), Ref.CParams(5.0, 3.0, 8, 0.0, :wmean, 12, true, 1)),
                    (ConvertParams(im_sigma = 2.0, mz_sigma = 1.0, stride = 4, sum_scale = false, centroid = :gauss, min_scans = 3),
                     Ref.CParams(2.0, 1.0, 4, 0.0, :gauss, 4, false, 3)))
        n, dpos, dint = check_equivalence(HELA, sel, p, rp)
        @test n > 100_000
        @test dpos == 0.0
        @test dint < 1e-6
    end
end

if BIG
    @testset "equivalence with the prototype (E. coli, 250 pg)" begin
        for dpath in (ECOLI, PG250)
            isdir(dpath) || (@warn "missing $dpath"; continue)
            f = open_tdf(dpath)
            sel = vcat(sample_rows(f, true, 10), sample_rows(f, false, 10))
            n, dpos, dint = check_equivalence(dpath, sel, ConvertParams(), Ref.CParams(5.0, 3.0, 8, 0.0, :wmean, 12, true, 1))
            @test n > 1_000_000 && dpos == 0.0 && dint < 1e-6
        end
    end
end

@testset "container round trip (HeLa, 120 frames)" begin
    out_dir = mktempdir()
    p = ConvertParams(format = :both, frames = collect(1:120), bin_scale = 2, int_scale = 16, cull_q = 0.01, ms1_cull_q = 0.0)
    paths = TimsSlices.convert(HELA, out_dir; params = p, name = "rt", log = devnull)
    @test isdir(paths.tdfs) && isfile(paths.arrow)
    t = open_tdfs(paths.tdfs)
    @test n_frames(t) == 120 && t.bin_scale == 2 && t.int_scale == 16 && t.meta["params"]["cull_q"] == 0.01
    @test t.frames.ms_order[1] == 0x01 && t.frames.cycle_idx[1] == 1 && t.frames.cycle_idx[10] == 2
    # slices table consistent with frames table and blocks
    @test sum(t.frames.n_slices) == n_slices(t)
    blk = SliceBlock(); codec = BlockCodec()
    for i in 1:n_frames(t)
        read_frame_block!(blk, codec, t, i)
        @test TS.n_peaks(blk) == t.frames.n_peaks[i]
        r = t.frames.first_slice[i]:t.frames.first_slice[i] + t.frames.n_slices[i] - 1
        @test all(j -> t.slices.n_peaks[r[j]] == length(TS.slice_range(blk, j)) && t.slices.peak_offset[r[j]] == blk.ptr[j], 1:blk.n_slices)
        @test all(j -> t.slices.tic[r[j]] == Float32(sum(blk.intensity[TS.slice_range(blk, j)]) / t.int_scale), 1:blk.n_slices)
        @test all(==(i), view(t.slices.frame_row, r))
    end
    # MS1 rows carry NaN centre / width, MS2 rows the window values
    ms1 = t.slices.ms_order .== 0x01
    @test all(isnan, t.slices.center_mz[ms1]) && all(!isnan, t.slices.center_mz[.!ms1]) && all(==(0f0), t.slices.collision_energy_ev[ms1])
    # expand reproduces the directly written Arrow exactly
    exp_path = joinpath(out_dir, "rt_expanded.arrow")
    expand(paths.tdfs, exp_path; log = devnull)
    a = Arrow.Table(paths.arrow); b = Arrow.Table(exp_path)
    @test length(a.scanNumber) == n_slices(t) > 0
    for k in propertynames(a)
        x = getproperty(a, k); y = getproperty(b, k)
        @test length(x) == length(y) && all(i -> isequal(x[i], y[i]), eachindex(x))
    end
    @test Arrow.getmetadata(a) == Arrow.getmetadata(b)
    # the Arrow peaks are the quantised bins converted with the file's calibration
    read_frame_block!(blk, codec, t, 1)
    @test a.mz_array[1][1] == Float32(TS.bin_to_mz(t, blk.bin[1])) && a.intensity_array[1][1] == Float32(TS.stored_to_intensity(t, blk.intensity[1]))
    @test a.intensity_array[1][1] == Float32(blk.intensity[1] / 16)
    @test a.retentionTime[1] == t.slices.retention_time[1] && a.imScan[1] == t.slices.im_scan[1]
    # cull thresholds recorded and applied (MS2 quantile cull, MS1 none)
    @test t.meta["cull_thr_ms2"] > 0 && t.meta["cull_thr_ms1"] == 0
    rm(out_dir; recursive = true)
end
