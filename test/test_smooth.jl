# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Synthetic tests of the smoothing pipeline: kernels, IM accumulation, m/z kernel + centroiding.

"A FrameBuffer from per-scan (bin, intensity) lists."
function make_frame(scans::Vector{Vector{Tuple{Int, Int}}})
    buf = FrameBuffer()
    ns = length(scans); np = sum(length, scans)
    buf.n_scans = ns; buf.n_peaks = np
    resize!(buf.scan_start, ns + 1); resize!(buf.tof, np); resize!(buf.intensity, np)
    k = 0
    for s in 1:ns
        buf.scan_start[s] = Int32(k + 1)
        for (b, i) in scans[s]; k += 1; buf.tof[k] = UInt32(b); buf.intensity[k] = UInt32(i); end
    end
    buf.scan_start[ns + 1] = Int32(np + 1)
    buf
end

lp(; im_sigma = 5.0, extent = 3.0, stride = 8, sum_scale = true, mz_sigma = 3.0, centroid = :wmean, max_half = 12, cull_q = 0.0, min_scans = 1, max_peaks = 0) =
    LevelParams(im_sigma, extent, stride, sum_scale, mz_sigma, centroid, max_half, cull_q, min_scans, max_peaks)

@testset "kernels" begin
    @test gauss_kernel(0.0) == [1.0]
    for (σ, e) in ((1.0, 3.0), (5.0, 3.0), (3.0, 4.0), (2.5, 3.0))
        k = gauss_kernel(σ, e)
        @test length(k) == 2 * ceil(Int, e * σ) + 1
        @test sum(k) ≈ 1.0 && k == reverse(k) && argmax(k) == length(k) ÷ 2 + 1
    end
    @test sum(TS.im_kernel(lp(stride = 8, sum_scale = true))) ≈ 8.0
    @test sum(TS.im_kernel(lp(stride = 8, sum_scale = false))) ≈ 1.0
    @test length(TS.mz_kernel(lp(mz_sigma = 3.0))) == 19
end

@testset "IM accumulation" begin
    # three scans (0, 1, 2); bin 100 in all of them, bin 105 in scan 0 only, bin 64 (word boundary) in scan 1
    buf = make_frame([[(100, 10), (105, 4)], [(64, 1), (100, 20)], [(100, 30)]])
    sc = SmoothScratch(); TS.ensure_bins!(sc, 200)
    l = lp(im_sigma = 1.0, stride = 1)
    kim = TS.im_kernel(l); h = length(kim) ÷ 2
    TS.im_accumulate!(sc, buf, 0, 3, 1, kim, h)
    @test sc.sp_bin == Int32[64, 100, 105]
    w = kim[h+1]; wm = kim[h]; wm2 = kim[h-1]     # centre weight, distance-1 and distance-2 weights (symmetric)
    @test sc.sp_val ≈ [w * 1, wm * 10 + w * 20 + wm * 30, wm * 4]
    @test sc.sp_cnt == Int32[1, 3, 1]
    # accumulators left clean
    @test all(==(0.0), sc.acc) && all(==(0), sc.cnt) && all(==(0), sc.bitmap)
    # slice at scan 0: scans 0, 1, 2 at distances 0, 1, 2; a window s0:s1 restricts the scans
    TS.im_accumulate!(sc, buf, 0, 3, 0, kim, h)
    @test sc.sp_bin == Int32[64, 100, 105] && sc.sp_val ≈ [wm * 1, w * 10 + wm * 20 + wm2 * 30, w * 4]
    TS.im_accumulate!(sc, buf, 2, 3, 2, kim, h)
    @test sc.sp_bin == Int32[100] && sc.sp_val ≈ [w * 30] && sc.sp_cnt == Int32[1]
    # empty reach
    TS.im_accumulate!(sc, buf, 1, 1, 1, kim, h)
    @test isempty(sc.sp_bin)
    # ensure_bins! grows without disturbing state
    TS.ensure_bins!(sc, 1_000_000)
    @test sc.n_bins == 1_000_000 && length(sc.bitmap) == cld(1_000_000, 64)
end

"Run the m/z stage + centroiding on one sparse slice (bins, values, counts)."
function centroid_sparse(bins, vals, cnts, l::LevelParams; thr = 0.0)
    sc = SmoothScratch()
    append!(sc.sp_bin, Int32.(bins)); append!(sc.sp_val, Float64.(vals)); append!(sc.sp_cnt, Int32.(cnts))
    out = FrameSlices()
    kmz = TS.mz_kernel(l)
    TS.centroid_slice!(out, sc, l, kmz, length(kmz) ÷ 2, thr)
    collect(zip(out.pos, out.val))
end

@testset "centroiding" begin
    l = lp(mz_sigma = 1.0, max_half = 4)
    # a single input bin: the smoothed profile is one Gaussian -> one centroid at the bin, footprint sum = kernel mass
    # within max_half (kernel half-width 3 <= max_half 4 -> the full mass)
    c = centroid_sparse([1000], [10.0], [1], l)
    @test length(c) == 1 && c[1][1] ≈ 1000.0 && c[1][2] ≈ 10.0
    # two ions 10 bins apart -> two centroids, symmetric footprints
    c = centroid_sparse([1000, 1010], [10.0, 20.0], [1, 1], l)
    @test length(c) == 2 && c[1][1] ≈ 1000.0 && c[2][1] ≈ 1010.0 && c[1][2] ≈ 10.0 && c[2][2] ≈ 20.0
    # two adjacent bins of equal intensity: one apex (the plateau rule sm[j] >= sm[j+1] takes the left one) whose
    # footprint does not cross the plateau (the walk needs a strict decrease), so the centroid sits left of 1000.5
    c = centroid_sparse([1000, 1001], [10.0, 10.0], [1, 1], l)
    @test length(c) == 1 && 999.0 < c[1][1] < 1000.5
    # a shoulder: 1000 (10) and 1001 (4) merge into one apex; the weighted mean sits between
    c = centroid_sparse([1000, 1001], [10.0, 4.0], [1, 1], l)
    @test length(c) == 1 && 1000.0 < c[1][1] < 1001.0
    # 3 bins apart at sigma 1 is resolved: two apices
    @test length(centroid_sparse([1000, 1003], [10.0, 4.0], [1, 1], l)) == 2
    # max_half truncates the footprint: with max_half 1 only apex +/- 1 is summed
    l1 = lp(mz_sigma = 1.0, max_half = 1)
    c = centroid_sparse([1000], [10.0], [1], l1)
    k = gauss_kernel(1.0); h = length(k) ÷ 2
    @test length(c) == 1 && c[1][2] ≈ 10.0 * (k[h] + k[h+1] + k[h+2])
    # threshold cull and min_scans cull
    @test isempty(centroid_sparse([1000], [10.0], [1], l; thr = 10.1))
    @test length(centroid_sparse([1000], [10.0], [1], l; thr = 9.999)) == 1
    @test isempty(centroid_sparse([1000], [10.0], [2], lp(mz_sigma = 1.0, max_half = 4, min_scans = 3)))
    @test length(centroid_sparse([1000, 1001], [10.0, 1.0], [2, 1], lp(mz_sigma = 1.0, max_half = 4, min_scans = 3))) == 1
    # gauss apex on a sampled Gaussian recovers the offset
    lg = lp(mz_sigma = 0.0, max_half = 4, centroid = :gauss)      # no extra smoothing: the input is the profile
    off = 0.3; bins = 995:1005; vals = [100 * exp(-((b - (1000 + off))^2) / 2) for b in bins]
    c = centroid_sparse(collect(bins), vals, ones(Int, length(bins)), lg)
    @test length(c) == 1 && abs(c[1][1] - (1000 + off)) < 0.01
    # variant A emits the sparse bins themselves
    out = FrameSlices(); sc = SmoothScratch()
    append!(sc.sp_bin, Int32[5, 9]); append!(sc.sp_val, [1.0, 0.4]); append!(sc.sp_cnt, Int32[2, 1])
    TS.emit_sparse!(out, sc, 0.5, 1)
    @test out.pos == [5.0] && out.val == [1.0]
end

@testset "window pipeline" begin
    # one ion persisting over scans 0..45 at bin 5000 (+/-1 jitter), one at bin 9000 in scans 40..50 only
    scans = Vector{Vector{Tuple{Int, Int}}}()
    for s in 0:63
        v = Tuple{Int, Int}[]
        0 <= s <= 45 && push!(v, (5000 + (s % 3) - 1, 100))
        40 <= s <= 50 && push!(v, (9000, 50))
        push!(scans, v)
    end
    buf = make_frame(scans)
    l = lp(im_sigma = 5.0, stride = 8, mz_sigma = 3.0)
    ls = LevelSetup(l, 0.0)
    out = FrameSlices(); sc = SmoothScratch(); TS.ensure_bins!(sc, 10_000)
    smooth_window!(out, sc, buf, 0, 64, 1, ls)
    # slices at 0, 8, 16, 24, 32, 40, 48, 56; 56 is out of reach of both ions (kernel half-width 15 -> reaches 41.. no: 56-15 = 41 -> in reach of 9000)
    @test out.scan == Int32[0, 8, 16, 24, 32, 40, 48, 56]
    @test all(==(1), out.window)
    for j in 1:out.n_slices
        r = out.ptr[j]:out.ptr[j+1]-1
        pos = out.pos[r]
        s = out.scan[j]
        has5000 = any(p -> abs(p - 5000) < 1.5, pos); has9000 = any(p -> abs(p - 9000) < 0.01, pos)
        @test has5000 == (s - 15 <= 45)          # reach of the first ion
        @test has9000 == (s + 15 >= 40)          # reach of the second
        @test length(pos) == has5000 + has9000    # the jittering ion is merged into ONE centroid
    end
    # the sum-scaled centroid intensity of a fully covered slice (scan 16 or 24: kernel reach inside 0..45) is stride * 100
    for sref in (16, 24)
        j = findfirst(==(sref), out.scan); r = out.ptr[j]:out.ptr[j+1]-1
        @test out.val[r][argmin(abs.(out.pos[r] .- 5000))] ≈ 8 * 100
    end
    # a window restricted to scans 40:64: slices start at the window start and only its scans contribute
    TS.reset!(out); smooth_window!(out, sc, buf, 40, 64, 2, ls)
    @test out.scan == Int32[40, 48, 56] && all(==(2), out.window)
    @test all(p -> abs(p - 9000) < 0.01 || abs(p - 5000) < 1.5, out.pos)
    # at slice 40 the first ion gets scans 40..45 only here, against 25..45 in the full window: a smaller sum
    val_at(o, sref, bin) = (j = findfirst(==(sref), o.scan); r = o.ptr[j]:o.ptr[j+1]-1; o.val[r][argmin(abs.(o.pos[r] .- bin))])
    v_restricted = val_at(out, 40, 5000)
    TS.reset!(out); smooth_window!(out, sc, buf, 0, 64, 1, ls)
    @test v_restricted < 0.6 * val_at(out, 40, 5000)
    @test val_at(out, 40, 9000) ≈ val_at(out, 40, 9000)

    # stride 1 without sum scaling: per-scan average intensity of the constant ion is 100
    l2 = lp(im_sigma = 5.0, stride = 1, sum_scale = false, mz_sigma = 3.0)
    TS.reset!(out); smooth_window!(out, sc, buf, 0, 64, 1, LevelSetup(l2, 0.0))
    j = findfirst(==(15), out.scan); r = out.ptr[j]:out.ptr[j+1]-1
    @test out.val[r][argmin(abs.(out.pos[r] .- 5000))] ≈ 100 atol = 1e-6
end

@testset "top-N cap per slice" begin
    fs = FrameSlices()
    for (pos, val) in ((10.0, 1.0), (20.0, 9.0), (30.0, 3.0), (40.0, 7.0), (50.0, 5.0)); TS.push_peak!(fs, pos, val); end
    TS.cap_slice!(fs, 3, Float64[])
    @test fs.pos == [20.0, 40.0, 50.0] && fs.val == [9.0, 7.0, 5.0]     # the 3 most intense, still in position order
    TS.end_slice!(fs, 0, 1)
    for (pos, val) in ((1.0, 2.0), (2.0, 2.0), (3.0, 1.0)); TS.push_peak!(fs, pos, val); end
    TS.cap_slice!(fs, 1, Float64[])                                        # ties at the threshold are kept
    @test fs.pos[fs.ptr[end]:end] == [1.0, 2.0]
    TS.cap_slice!(fs, 5, Float64[])                                        # no-op below the cap
    @test length(fs.pos) == 5 && fs.ptr == Int32[1, 4]
    # through the pipeline: the capped level keeps the n most intense centroids of every slice
    scans = Vector{Vector{Tuple{Int, Int}}}()
    for s in 0:31
        push!(scans, [(1000 + 20i, 10 + i) for i in 1:30])      # 30 ions, intensities 11..40
    end
    buf = make_frame(scans); sc = SmoothScratch(); TS.ensure_bins!(sc, 10_000); out = FrameSlices()
    smooth_window!(out, sc, buf, 0, 32, 1, LevelSetup(lp(mz_sigma = 1.0, max_peaks = 5), 0.0))
    @test all(j -> out.ptr[j+1] - out.ptr[j] == 5, 1:out.n_slices)
    @test all(j -> all(p -> p > 1000 + 20 * 25 - 1, out.pos[out.ptr[j]:out.ptr[j+1]-1]), 1:out.n_slices)   # the 5 brightest ions
    @test output_name("r.d", ConvertParams(max_peaks = 1000)) == "r_cen_s5_m3_k8_q0_wmean_sum_top1000"
    # the cap must not disturb the m/z-stage buffers: a wide slice (many peaks) capped, then a long run
    # (regression: the cap once grew sc.dense past sc.dcnt and a later run wrote out of bounds)
    @test length(sc.dense) == length(sc.dcnt)
    wide = [[(500 + 3i, 100 + i) for i in 1:4000] for _ in 1:8]          # 4000 peaks per scan, 3 bins apart (one long run), distinct intensities
    buf2 = make_frame(wide); TS.reset!(out); TS.ensure_bins!(sc, 20_000)
    smooth_window!(out, sc, buf2, 0, 8, 1, LevelSetup(lp(mz_sigma = 1.0, stride = 8, max_peaks = 100), 0.0))
    @test length(sc.dense) == length(sc.dcnt) && out.n_slices == 1 && length(out.pos) == 100
end

@testset "params" begin
    p = ConvertParams()
    @test p.max_half == 12 && p.ms1_stride == 8 && p.ms1_cull_q == 0.0 && !p.split_cull
    @test ConvertParams(cull_q = 0.05, ms1_cull_q = 0.0).split_cull
    @test ConvertParams(mz_sigma = 1.0).max_half == 4
    @test_throws ArgumentError TS.validate(ConvertParams(centroid = :apex))
    @test_throws ArgumentError TS.validate(ConvertParams(bin_scale = 0))
    @test_throws ArgumentError TS.validate(ConvertParams(format = :csv))
    @test output_name("x/run.d", ConvertParams()) == "run_cen_s5_m3_k8_q0_wmean_sum"
    @test output_name("run.d", ConvertParams(cull_q = 0.05, ms1_cull_q = 0.0, bin_scale = 256, min_scans = 3)) == "run_cen_s5_m3_k8_q0.05_wmean_sum_n3_ms1q0_b256"
    p2, _ = TS.parse_cli(["a.d", "out", "--mz-sigma", "2", "--centroid", "none", "--no-sum-scale", "--frames", "1:10", "--format", "both"])
    @test p2.mz_sigma == 2.0 && p2.centroid == :none && !p2.sum_scale && p2.frames == collect(1:10) && p2.format == :both
    @test_throws ErrorException TS.parse_cli(["a.d", "out", "--bogus", "1"])
end
