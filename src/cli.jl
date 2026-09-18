# Command line: julia -t N bin/tims_convert.jl <run.d> <out_dir> [--key value ...]
#               julia bin/tims_convert.jl expand <name.tdfs> <out.arrow>

const CLI_HELP = """
tims_convert.jl <run.d> <out_dir> [options]      convert a .d bundle
tims_convert.jl expand <name.tdfs> <out.arrow>   expand a tdfs into Pioneer's slice Arrow

Options (defaults in brackets):
  --im-sigma S [5]  --ms1-im-sigma S  --kernel-extent E [3]  --stride K [8]  --ms1-stride K  --no-sum-scale
  --mz-sigma S [3]  --centroid wmean|gauss|none [wmean]  --max-half H [max(4, 4*mz-sigma)]
  --cull-q Q [0]  --ms1-cull-q Q  --split-cull  --min-scans N [1]  --cull-sample-frames N [40]
  --bin-scale K [1]  --int-scale X [1]  --zstd-level L [3]  --format tdfs|arrow|both [tdfs]
  --batch-frames N [16*threads, results in flight]  --frames a:b (frame rows)  --name NAME
"""

function parse_cli(args::Vector{String})
    length(args) >= 2 || (println(CLI_HELP); error("not enough arguments"))
    kw = Dict{Symbol, Any}(); name = nothing
    i = 3
    flags = Dict("--im-sigma" => (:im_sigma, Float64), "--ms1-im-sigma" => (:ms1_im_sigma, Float64), "--kernel-extent" => (:kernel_extent, Float64),
                 "--stride" => (:stride, Int), "--ms1-stride" => (:ms1_stride, Int), "--mz-sigma" => (:mz_sigma, Float64),
                 "--centroid" => (:centroid, Symbol), "--max-half" => (:max_half, Int), "--cull-q" => (:cull_q, Float64),
                 "--ms1-cull-q" => (:ms1_cull_q, Float64), "--min-scans" => (:min_scans, Int), "--cull-sample-frames" => (:cull_sample_frames, Int),
                 "--bin-scale" => (:bin_scale, Int), "--int-scale" => (:int_scale, Float64), "--zstd-level" => (:zstd_level, Int),
                 "--format" => (:format, Symbol), "--batch-frames" => (:batch_frames, Int))
    while i <= length(args)
        a = args[i]
        if a == "--no-sum-scale"; kw[:sum_scale] = false; i += 1
        elseif a == "--split-cull"; kw[:split_cull] = true; i += 1
        elseif a == "--frames"
            lo, hi = split(args[i + 1], ':'); kw[:frames] = collect(parse(Int, lo):parse(Int, hi)); i += 2
        elseif a == "--name"; name = args[i + 1]; i += 2
        elseif haskey(flags, a)
            k, T = flags[a]; v = args[i + 1]
            kw[k] = T === Symbol ? Symbol(v) : T === Int ? parse(Int, v) : parse(Float64, v); i += 2
        else
            println(CLI_HELP); error("unknown option $a")
        end
    end
    ConvertParams(; kw...), name
end

function main(args::Vector{String} = ARGS)
    if !isempty(args) && args[1] == "expand"
        length(args) == 3 || (println(CLI_HELP); error("expand takes <name.tdfs> <out.arrow>"))
        expand(args[2], args[3])
        return
    end
    p, name = parse_cli(args)
    # `--frames` filter needs the file to validate; convert does that
    if name === nothing
        convert(args[1], args[2]; params = p)
    else
        convert(args[1], args[2]; params = p, name = name)
    end
    nothing
end
