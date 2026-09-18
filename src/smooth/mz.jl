# Step 2: m/z kernel over runs of nearby bins of a slice's sparse list, by scattering each nonzero bin across
# the kernel taps into a dense run buffer padded by h_mz on each side, then centroiding the run.

"""
    centroid_slice!(out, sc, lp, kmz, h_mz, thr)

Consume `sc.sp_*` (one slice, sorted by bin) and append the slice's centroids to `out`.
"""
function centroid_slice!(out::FrameSlices, sc::SmoothScratch, lp::LevelParams, kmz::Vector{Float64}, h_mz::Int, thr::Float64)
    bins = sc.sp_bin; vals = sc.sp_val; cnts = sc.sp_cnt
    n = length(bins); gap_max = 2h_mz + 2
    ntaps = length(kmz)
    i = 1
    @inbounds while i <= n
        j = i
        while j < n && bins[j + 1] - bins[j] <= gap_max; j += 1; end
        t0 = Int(bins[i]) - h_mz; L = Int(bins[j]) - Int(bins[i]) + 1 + 2h_mz
        if length(sc.dense) < L
            resize!(sc.dense, L); resize!(sc.dcnt, L)
        end
        dense = sc.dense; dcnt = sc.dcnt
        @simd for a in 1:L
            dense[a] = 0.0; dcnt[a] = Int32(0)
        end
        for q in i:j
            p = Int(bins[q]) - t0 + 1; v = vals[q]
            dcnt[p] = cnts[q]
            base = p - h_mz - 1
            @simd for t in 1:ntaps
                dense[base + t] += kmz[t] * v
            end
        end
        centroid_dense!(out, dense, dcnt, L, t0, thr, lp, h_mz + 1, L - h_mz)
        i = j + 1
    end
    out
end

"Variant A (`centroid = :none`): emit the IM-accumulated per-bin sums as peaks, subject to the culls."
function emit_sparse!(out::FrameSlices, sc::SmoothScratch, thr::Float64, min_scans::Int)
    @inbounds for q in eachindex(sc.sp_bin)
        v = sc.sp_val[q]
        (v >= thr && sc.sp_cnt[q] >= min_scans) || continue
        push_peak!(out, Float64(sc.sp_bin[q]), v)
    end
    out
end
