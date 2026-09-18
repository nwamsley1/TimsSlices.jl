# The per-window pipeline: slices every `stride` scans of window scans s0:s1-1 of one decoded frame.

"Per-level constants derived from the parameters once per conversion."
struct LevelSetup
    lp::LevelParams
    kim::Vector{Float64}
    kmz::Vector{Float64}
    h_im::Int
    h_mz::Int
    thr::Float64          # cull threshold in intensity units (0 = off)
end
function LevelSetup(lp::LevelParams, thr::Float64)
    kim = im_kernel(lp); kmz = mz_kernel(lp)
    LevelSetup(lp, kim, kmz, length(kim) ÷ 2, length(kmz) ÷ 2, thr)
end

"""
    smooth_window!(out, sc, buf, s0, s1, widx, ls)

Append the slices of window `widx` (scans `s0:s1-1`) of the frame in `buf` to `out`. A slice with no raw
entries in reach, or no surviving centroids, is not recorded.
"""
function smooth_window!(out::FrameSlices, sc::SmoothScratch, buf::FrameBuffer, s0::Integer, s1::Integer, widx::Integer, ls::LevelSetup)
    lp = ls.lp
    for c in s0:lp.stride:s1-1
        im_accumulate!(sc, buf, s0, s1, c, ls.kim, ls.h_im)
        isempty(sc.sp_bin) && continue
        if lp.centroid == :none
            emit_sparse!(out, sc, ls.thr, lp.min_scans)
        else
            centroid_slice!(out, sc, lp, ls.kmz, ls.h_mz, ls.thr)
        end
        end_slice!(out, c, widx)
    end
    out
end

"""
    smooth_frame!(out, sc, buf, file, i, setup_ms1, setup_ms2)

All windows of frame row `i` (already decoded into `buf`) -> `out` (reset first).
"""
function smooth_frame!(out::FrameSlices, sc::SmoothScratch, buf::FrameBuffer, f::TdfFile, i::Integer, ls1::LevelSetup, ls2::LevelSetup)
    reset!(out)
    buf.n_peaks == 0 && return out
    ensure_bins!(sc, max(n_bins(f), Int(maximum(view(buf.tof, 1:buf.n_peaks))) + 1))
    ls = is_ms1(f, i) ? ls1 : ls2
    for (widx, w) in enumerate(windows(f, i))
        smooth_window!(out, sc, buf, w.scan_begin, w.scan_end, widx, ls)
    end
    out
end
