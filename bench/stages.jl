# Per-stage cost on sampled frames, single thread.
# Usage: julia --project=bench bench/stages.jl <run.d> [n_frames=40] [--mz-sigma 3] ...
using TimsSlices, Statistics, Printf
const TS = TimsSlices
dpath = ARGS[1]; nf = length(ARGS) >= 2 && !startswith(ARGS[2], "--") ? parse(Int, ARGS[2]) : 40
p, _ = TS.parse_cli(vcat([dpath, "x"], ARGS[findfirst(a -> startswith(a, "--"), ARGS) === nothing ? (1:0) : (findfirst(a -> startswith(a, "--"), ARGS):end)]))
f = open_tdf(dpath)
ls1 = LevelSetup(level_params(p, true), 0.0); ls2 = LevelSetup(level_params(p, false), 0.0)
buf = FrameBuffer(); sc = SmoothScratch(); out = FrameSlices(); blk = SliceBlock(); codec = BlockCodec(); rblk = SliceBlock()
sample(ms1) = (v = [i for i in 1:n_frames(f) if (f.frames.msms_type[i] == 0) == ms1]; v[unique(round.(Int, range(1, length(v); length = min(nf, length(v)))))])
function stage_times(i)
    t0 = time_ns(); read_frame!(buf, f, i); t1 = time_ns()
    smooth_frame!(out, sc, buf, f, i, ls1, ls2); t2 = time_ns()
    quantize!(blk, out, p.bin_scale, p.int_scale); t3 = time_ns()
    nw, nb = encode_block!(codec, blk, p.zstd_level); t4 = time_ns()
    decode_block!(rblk, codec, codec.zbuf[1:nb], nw); t5 = time_ns()
    (dec = (t1 - t0) / 1e6, smooth = (t2 - t1) / 1e6, quant = (t3 - t2) / 1e6, enc = (t4 - t3) / 1e6, bdec = (t5 - t4) / 1e6,
     raw = buf.n_peaks, cen = TS.n_peaks(blk), bytes = nb, slices = blk.n_slices)
end
stage_times(sample(true)[1]); stage_times(sample(false)[1])   # JIT
println("file ", basename(dpath), "  params ", p)
tot = Dict{Symbol, Float64}()
for (name, ms1) in (("MS1", true), ("MS2", false))
    rows = sample(ms1); r = [stage_times(i) for i in rows]
    med(k) = median(getindex.(r, k))
    n_level = count(i -> (f.frames.msms_type[i] == 0) == ms1, 1:n_frames(f))
    @printf("%s (%d frames sampled, %d in file): raw %d peaks -> %d slices, %d centroids, %.2f B/centroid\n", name, length(rows), n_level, med(:raw), med(:slices), med(:cen), med(:bytes) / med(:cen))
    @printf("  decode .d   %6.1f ms  %5.1f ns/raw peak\n", med(:dec), med(:dec) * 1e6 / med(:raw))
    @printf("  smooth      %6.1f ms  %5.1f ns/raw peak  %5.1f ns/centroid\n", med(:smooth), med(:smooth) * 1e6 / med(:raw), med(:smooth) * 1e6 / med(:cen))
    @printf("  quantize    %6.1f ms  %5.1f ns/centroid\n", med(:quant), med(:quant) * 1e6 / med(:cen))
    @printf("  encode      %6.1f ms  %5.1f ns/centroid\n", med(:enc), med(:enc) * 1e6 / med(:cen))
    @printf("  decode tdfs %6.1f ms  %5.1f ns/centroid\n", med(:bdec), med(:bdec) * 1e6 / med(:cen))
    for k in (:dec, :smooth, :quant, :enc, :bdec); tot[k] = get(tot, k, 0.0) + n_level * med(k) / 1e3; end
end
@printf("whole file, single thread (medians x counts): decode %.0f s, smooth %.0f s, quantize %.0f s, encode %.0f s; tdfs decode per pass %.0f s\n",
        tot[:dec], tot[:smooth], tot[:quant], tot[:enc], tot[:bdec])
