# The converter driver: cull sampling, frame batches over worker tasks with per-task scratch, ordered writers.
using Statistics, Printf

"Everything one worker task owns."
mutable struct Worker
    buf::FrameBuffer
    sc::SmoothScratch
    out::FrameSlices
    blk::SliceBlock
    codec::BlockCodec
    rows::SliceRows
    t_decode::Float64; t_smooth::Float64; t_encode::Float64
end
Worker() = Worker(FrameBuffer(), SmoothScratch(), FrameSlices(), SliceBlock(), BlockCodec(), SliceRows(), 0.0, 0.0, 0.0)

"What one frame yields for the writers (owned by the main task after the batch)."
struct FrameResult
    seq::Int                            # position in the conversion order (for the ordered writer)
    row::Int
    fm::FrameMeta
    rows::SliceRows
    blk::Union{Nothing, SliceBlock}     # kept only when the Arrow output needs the peaks
    n_peaks::Int
    zbytes::Vector{UInt8}
    n_words::Int
end

function frame_meta(f::TdfFile, i::Integer)
    FrameMeta(f.frames.id[i], is_ms1(f, i) ? 0x01 : 0x02, UInt8(f.dia.frame_group[i]), f.frames.time[i], f.frames.ramp_time[i], f.frames.num_scans[i])
end

"Raw-intensity quantile thresholds (MS1, MS2) from sample frames; 0 when the quantile is 0."
function cull_thresholds(f::TdfFile, rows::Vector{Int}, p::ConvertParams)
    (p.cull_q == 0 && p.ms1_cull_q == 0) && return 0.0, 0.0
    buf = FrameBuffer()
    sample_of(idxs) = isempty(idxs) ? Int[] : idxs[unique(round.(Int, range(1, length(idxs); length = min(p.cull_sample_frames, length(idxs)))))]
    function raw_of(idxs)
        v = Float64[]
        for i in idxs
            read_frame!(buf, f, i); append!(v, Float64.(view(buf.intensity, 1:buf.n_peaks)))
        end
        v
    end
    if p.split_cull
        raw1 = raw_of(sample_of(filter(i -> is_ms1(f, i), rows))); raw2 = raw_of(sample_of(filter(i -> !is_ms1(f, i), rows)))
        thr1 = p.ms1_cull_q > 0 && !isempty(raw1) ? quantile(raw1, p.ms1_cull_q) : 0.0
        thr2 = p.cull_q > 0 && !isempty(raw2) ? quantile(raw2, p.cull_q) : 0.0
        return thr1, thr2
    else
        raw = raw_of(sample_of(rows))
        thr = p.cull_q > 0 && !isempty(raw) ? quantile(raw, p.cull_q) : 0.0
        return thr, thr
    end
end

function process_frame!(wk::Worker, f::TdfFile, i::Int, ls1::LevelSetup, ls2::LevelSetup, p::ConvertParams, keep_blk::Bool, seq::Int = 0)
    t0 = time_ns()
    read_frame!(wk.buf, f, i)
    t1 = time_ns()
    smooth_frame!(wk.out, wk.sc, wk.buf, f, i, ls1, ls2)
    t2 = time_ns()
    quantize!(wk.blk, wk.out, p.bin_scale, p.int_scale)
    n_words, nb = encode_block!(wk.codec, wk.blk, p.zstd_level)
    fm = frame_meta(f, i)
    slice_rows!(wk.rows, fm, wk.blk, wk.out.scan, wk.out.window, windows(f, i), f.ce_ramp, p.int_scale)
    t3 = time_ns()
    wk.t_decode += (t1 - t0) / 1e9; wk.t_smooth += (t2 - t1) / 1e9; wk.t_encode += (t3 - t2) / 1e9
    FrameResult(seq, i, fm, deepcopy(wk.rows), keep_blk ? deepcopy(wk.blk) : nothing, n_peaks(wk.blk), wk.codec.zbuf[1:nb], n_words)
end

"""
    convert(dir, out_dir; params = ConvertParams(), name = output_name(dir, params), log = stdout) -> (tdfs, arrow)

Convert a `.d` bundle. Returns the paths written (`nothing` for a format not requested).
"""
function convert(dir::AbstractString, out_dir::AbstractString; params::ConvertParams = ConvertParams(),
                 name::AbstractString = output_name(dir, params), log::IO = stdout, _fail_frames::Bool = false)
    p = validate(params)
    t_start = time()
    f = open_tdf(dir)
    rows = p.frames === nothing ? valid_frames(f) : sort(p.frames)
    all(i -> 1 <= i <= n_frames(f) && (is_ms1(f, i) || is_dia(f, i)), rows) || throw(ArgumentError("frames must be rows of MS1 / diaPASEF frames"))
    @printf(log, "source %s: %d frames (%d MS1, %d MS2), %d raw peaks; mz cal residual ppm %s\n", basename(rstrip(dir, '/')), length(rows),
            count(i -> is_ms1(f, i), rows), count(i -> !is_ms1(f, i), rows), sum(Int, f.frames.num_peaks[rows]), string(round.(f.mz_cal_resid_ppm, digits = 2)))
    thr1, thr2 = cull_thresholds(f, rows, p)
    ls1 = LevelSetup(level_params(p, true), thr1); ls2 = LevelSetup(level_params(p, false), thr2)
    @printf(log, "params %s\ncull thresholds: MS1 %.1f (q%g), MS2 %.1f (q%g)\n", string(Dict(p)), thr1, p.ms1_cull_q, thr2, p.cull_q)

    mz_lo = parse(Float64, f.meta["MzAcqRangeLower"]); mz_hi = parse(Float64, f.meta["MzAcqRangeUpper"])
    meta = Dict{String, Any}(
        "source" => basename(rstrip(dir, '/')), "source_bin_bytes" => f.bin_size, "instrument" => get(f.meta, "InstrumentName", ""),
        "mz_cal_sqrt_intercept" => f.mz_cal.intercept, "mz_cal_sqrt_slope" => f.mz_cal.slope,
        "im_scan0_1overK0" => f.im_cal.intercept, "im_slope_1overK0_per_scan" => f.im_cal.slope,
        "ce_ev_intercept" => f.ce_ramp.intercept, "ce_ev_slope_per_scan" => f.ce_ramp.slope,
        "NumScans" => f.max_scans, "n_bins" => n_bins(f), "mz_lo" => mz_lo, "mz_hi" => mz_hi,
        "OneOverK0AcqRangeLower" => f.meta["OneOverK0AcqRangeLower"], "OneOverK0AcqRangeUpper" => f.meta["OneOverK0AcqRangeUpper"],
        "params" => Dict(p), "cull_thr_ms1" => thr1, "cull_thr_ms2" => thr2,
        "bin_scale" => p.bin_scale, "int_scale" => p.int_scale, "zstd_level" => p.zstd_level,
        "converter" => "TimsSlices.jl $(pkgversion(TimsSlices))", "converted_at" => string(now_utc()))
    mkpath(out_dir)
    want_tdfs = p.format in (:tdfs, :both); want_arrow = p.format in (:arrow, :both)
    tdfs_path = want_tdfs ? joinpath(out_dir, name * ".tdfs") : nothing
    arrow_path = want_arrow ? joinpath(out_dir, name * ".arrow") : nothing
    tw = want_tdfs ? TdfsWriter(tdfs_path, meta) : nothing
    aw = want_arrow ? SliceArrowWriter(arrow_path, arrow_metadata(meta, p, thr1, thr2), f.mz_cal, p.bin_scale, p.int_scale, mz_lo, mz_hi) : nothing

    nt = Threads.nthreads()
    workers = [Worker() for _ in 1:nt]
    inflight = p.batch_frames > 0 ? p.batch_frames : 16nt
    n_rows = length(rows)
    # workers pull frame indices from a counter and push results into a bounded channel (memory bound = inflight
    # results); the main task reorders them and writes frames in order
    next = Threads.Atomic{Int}(1)
    results = Channel{FrameResult}(inflight)
    worker_tasks = [Threads.@spawn begin
        wk = workers[t]
        try
            while true
                q = Threads.atomic_add!(next, 1)
                q > n_rows && break
                _fail_frames && error("injected failure (test)")
                put!(results, process_frame!(wk, f, rows[q], ls1, ls2, p, want_arrow, q))
            end
        catch e
            # a failed worker must not leave the writer waiting forever: close the channel with the error
            isopen(results) && close(results, ErrorException("worker failed on frame row $(rows[min(Threads.atomic_add!(next, 0) - 1, n_rows)]): $(sprint(showerror, e))"))
            rethrow()
        end
    end for t in 1:nt]
    done = 0; t_loop = time(); n_slices_tot = 0; n_peaks_tot = Int[0, 0]
    pending = Dict{Int, FrameResult}()
    next_write = 1
    while next_write <= n_rows
        r = try
            take!(results)
        catch e
            # the channel was closed by a failing worker: surface that worker's exception
            for t in worker_tasks; istaskfailed(t) && wait(t); end
            rethrow(e)
        end
        pending[r.seq] = r
        while haskey(pending, next_write)
            r = pop!(pending, next_write)
            want_tdfs && write_frame!(tw, r.fm, r.rows, r.n_peaks, r.zbytes, r.n_words)
            want_arrow && write_frame!(aw, r.fm, r.rows, r.blk)
            n_slices_tot += length(r.rows); n_peaks_tot[r.fm.ms_order] += r.n_peaks
            next_write += 1; done += 1
            if done % 2000 == 0 || done == n_rows
                @printf(log, "  %d / %d frames, %.1f s\n", done, n_rows, time() - t_loop)
            end
        end
    end
    foreach(wait, worker_tasks)
    t_proc = time() - t_loop
    want_tdfs && close(tw); want_arrow && close(aw)
    close(f)
    td = sum(w.t_decode for w in workers); ts = sum(w.t_smooth for w in workers); te = sum(w.t_encode for w in workers)
    @printf(log, "frames %d -> slices %d, centroids MS1 %d + MS2 %d = %d\n", length(rows), n_slices_tot, n_peaks_tot[1], n_peaks_tot[2], sum(n_peaks_tot))
    @printf(log, "CPU seconds: decode %.1f, smooth %.1f, encode+rows %.1f; wall %.1f s on %d threads (write included); total %.1f s\n",
            td, ts, te, t_proc, nt, time() - t_start)
    want_tdfs && @printf(log, "wrote %s  %.3f GB\n", tdfs_path, dir_size(tdfs_path) / 1e9)
    want_arrow && @printf(log, "wrote %s  %.3f GB\n", arrow_path, filesize(arrow_path) / 1e9)
    (tdfs = tdfs_path, arrow = arrow_path)
end

function arrow_metadata(meta::Dict{String, Any}, p::ConvertParams, thr1, thr2)
    m = Dict{String, String}()
    for k in ("source", "instrument", "mz_cal_sqrt_intercept", "mz_cal_sqrt_slope", "im_scan0_1overK0", "im_slope_1overK0_per_scan",
              "ce_ev_intercept", "ce_ev_slope_per_scan", "NumScans", "OneOverK0AcqRangeLower", "OneOverK0AcqRangeUpper")
        m[k] = string(meta[k])
    end
    m["centroid_im_sigma"] = string(p.im_sigma); m["centroid_mz_sigma"] = string(p.mz_sigma); m["centroid_stride"] = string(p.stride)
    m["centroid_cull_q"] = string(p.cull_q); m["centroid_cull_thr"] = string(thr2); m["centroid_method"] = string(p.centroid)
    m["centroid_sum_scale"] = string(p.sum_scale); m["centroid_min_scans"] = string(p.min_scans)
    m["centroid_ms1_stride"] = string(p.ms1_stride); m["centroid_ms1_cull_q"] = string(p.ms1_cull_q)
    m["centroid_split_cull"] = string(p.split_cull); m["centroid_cull_thr_ms1"] = string(thr1); m["centroid_cull_thr_ms2"] = string(thr2)
    m["bin_scale"] = string(p.bin_scale); m["int_scale"] = string(p.int_scale)
    m
end

dir_size(d) = sum(filesize(joinpath(d, x)) for x in readdir(d); init = 0)
now_utc() = Libc.strftime("%Y-%m-%dT%H:%M:%SZ", time())
