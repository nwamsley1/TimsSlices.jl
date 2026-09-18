# The .tdfs container: blocks.bin (u32 block_size, u32 n_words, zstd payload per frame), frames.arrow,
# slices.arrow, meta.json. See docs/format.md.
using Arrow, JSON3

const TDFS_FORMAT_VERSION = 1

"Per-frame metadata carried from the .d (or the tdfs) to the writers."
struct FrameMeta
    frame_id::Int32
    ms_order::UInt8
    window_group::UInt8
    rt_s::Float64
    ramp_ms::Float64
    n_scans::Int32
end

"Per-slice metadata of one frame (struct of vectors, appended by the converter)."
mutable struct SliceRows
    im_scan::Vector{UInt16}
    window::Vector{Int32}              # window index within the frame (1-based)
    retention_time::Vector{Float32}    # minutes
    center_mz::Vector{Float32}         # NaN for MS1
    isolation_width::Vector{Float32}   # NaN for MS1
    collision_energy_ev::Vector{Float32}   # CE ramp at the slice scan (0 for MS1)
    window_ce::Vector{Float32}             # the window table's CE (NaN for MS1)
    tic::Vector{Float32}
    n_peaks::Vector{Int32}
    peak_offset::Vector{Int64}         # 1-based start of the slice's peaks in the frame's peak arrays
end
SliceRows() = SliceRows(UInt16[], Int32[], Float32[], Float32[], Float32[], Float32[], Float32[], Float32[], Int32[], Int64[])
Base.length(r::SliceRows) = length(r.im_scan)
function Base.empty!(r::SliceRows)
    for n in fieldnames(SliceRows); empty!(getfield(r, n)); end
    r
end
function Base.append!(a::SliceRows, b::SliceRows)
    for n in fieldnames(SliceRows); append!(getfield(a, n), getfield(b, n)); end
    a
end

"""
    slice_rows!(rows, fm, blk, scan, window, wins, ce_ramp)

Fill `rows` for one frame from its quantised block and per-slice scan / window indices.
"""
function slice_rows!(rows::SliceRows, fm::FrameMeta, blk::SliceBlock, scan::Vector{Int32}, window::Vector{Int32}, wins, ce::CeRamp)
    empty!(rows)
    dt = fm.ramp_ms / 1000 / fm.n_scans        # seconds per scan
    ms1 = fm.ms_order == 0x01
    @inbounds for j in 1:blk.n_slices
        s = scan[j]; w = wins[window[j]]
        r = slice_range(blk, j)
        tic = 0.0
        for k in r; tic += blk.intensity[k]; end
        push!(rows.im_scan, UInt16(s)); push!(rows.window, window[j])
        push!(rows.retention_time, Float32((fm.rt_s + s * dt) / 60))
        push!(rows.center_mz, ms1 ? NaN32 : w.center); push!(rows.isolation_width, ms1 ? NaN32 : w.width)
        push!(rows.collision_energy_ev, ms1 ? 0f0 : Float32(ce_at(ce, s))); push!(rows.window_ce, ms1 ? NaN32 : w.ce)
        push!(rows.tic, Float32(tic)); push!(rows.n_peaks, Int32(length(r))); push!(rows.peak_offset, Int64(first(r)))
    end
    rows
end

# --- writer -----------------------------------------------------------------------------------------------
mutable struct TdfsWriter
    dir::String
    blocks::IOStream
    offset::Int64
    cycle::Int32
    # frames.arrow columns
    f_frame_id::Vector{Int32}; f_ms_order::Vector{UInt8}; f_cycle_idx::Vector{Int32}; f_window_group::Vector{UInt8}
    f_rt_s::Vector{Float64}; f_n_scans::Vector{Int32}; f_n_slices::Vector{Int32}; f_n_peaks::Vector{Int64}
    f_block_offset::Vector{Int64}; f_block_size::Vector{Int64}; f_n_words::Vector{Int64}; f_first_slice::Vector{Int64}
    # slices.arrow columns
    s_frame_row::Vector{Int32}; s_slice_in_frame::Vector{Int32}; s_ms_order::Vector{UInt8}; s_cycle_idx::Vector{Int32}
    s_window_group::Vector{UInt8}; s_frame_id::Vector{Int32}
    s_rows::SliceRows
    meta::Dict{String, Any}
end

function TdfsWriter(dir::AbstractString, meta::Dict{String, Any})
    mkpath(dir)
    io = open(joinpath(dir, "blocks.bin"), "w")
    TdfsWriter(String(dir), io, 0, 0,
               Int32[], UInt8[], Int32[], UInt8[], Float64[], Int32[], Int32[], Int64[], Int64[], Int64[], Int64[], Int64[],
               Int32[], Int32[], UInt8[], Int32[], UInt8[], Int32[], SliceRows(), meta)
end

"Append one frame: its compressed block bytes, its slice rows, its frame row. Frames must arrive in order."
function write_frame!(w::TdfsWriter, fm::FrameMeta, rows::SliceRows, n_peaks::Integer, zbytes::AbstractVector{UInt8}, n_words::Integer)
    fm.ms_order == 0x01 && (w.cycle += Int32(1))
    nb = length(zbytes)
    block_size = 8 + nb
    write(w.blocks, UInt32(block_size)); write(w.blocks, UInt32(n_words)); write(w.blocks, zbytes)
    frow = Int32(length(w.f_frame_id) + 1)
    push!(w.f_frame_id, fm.frame_id); push!(w.f_ms_order, fm.ms_order); push!(w.f_cycle_idx, w.cycle); push!(w.f_window_group, fm.window_group)
    push!(w.f_rt_s, fm.rt_s); push!(w.f_n_scans, fm.n_scans); push!(w.f_n_slices, Int32(length(rows))); push!(w.f_n_peaks, Int64(n_peaks))
    push!(w.f_block_offset, w.offset); push!(w.f_block_size, block_size); push!(w.f_n_words, Int64(n_words)); push!(w.f_first_slice, Int64(length(w.s_frame_row) + 1))
    for j in 1:length(rows)
        push!(w.s_frame_row, frow); push!(w.s_slice_in_frame, Int32(j)); push!(w.s_ms_order, fm.ms_order)
        push!(w.s_cycle_idx, w.cycle); push!(w.s_window_group, fm.window_group); push!(w.s_frame_id, fm.frame_id)
    end
    append!(w.s_rows, rows)
    w.offset += block_size
    w
end

function Base.close(w::TdfsWriter)
    close(w.blocks)
    Arrow.write(joinpath(w.dir, "frames.arrow"),
        (frame_id = w.f_frame_id, ms_order = w.f_ms_order, cycle_idx = w.f_cycle_idx, window_group = w.f_window_group,
         rt_s = w.f_rt_s, n_scans = w.f_n_scans, n_slices = w.f_n_slices, n_peaks = w.f_n_peaks,
         block_offset = w.f_block_offset, block_size = w.f_block_size, n_words = w.f_n_words, first_slice = w.f_first_slice))
    r = w.s_rows
    Arrow.write(joinpath(w.dir, "slices.arrow"),
        (frame_row = w.s_frame_row, slice_in_frame = w.s_slice_in_frame, frame_id = w.s_frame_id, im_scan = r.im_scan,
         window = r.window, ms_order = w.s_ms_order, cycle_idx = w.s_cycle_idx, window_group = w.s_window_group,
         retention_time = r.retention_time, center_mz = r.center_mz, isolation_width = r.isolation_width,
         collision_energy_ev = r.collision_energy_ev, window_ce = r.window_ce, tic = r.tic, n_peaks = r.n_peaks, peak_offset = r.peak_offset))
    m = copy(w.meta)
    m["format_version"] = TDFS_FORMAT_VERSION
    m["n_frames"] = length(w.f_frame_id); m["n_slices"] = length(w.s_frame_row); m["n_peaks"] = sum(w.f_n_peaks; init = 0)
    m["blocks_bytes"] = w.offset
    open(joinpath(w.dir, "meta.json"), "w") do io
        JSON3.pretty(io, m)
    end
    nothing
end

# --- reader -----------------------------------------------------------------------------------------------
struct TdfsFile
    dir::String
    meta::Dict{String, Any}
    frames::NamedTuple
    slices::NamedTuple
    blocks::Vector{UInt8}
    mz_cal::LinearMzCal
    bin_scale::Int
    int_scale::Float64
end

function open_tdfs(dir::AbstractString)
    meta = Dict{String, Any}(JSON3.read(read(joinpath(dir, "meta.json"), String), Dict{String, Any}))
    Int(meta["format_version"]) == TDFS_FORMAT_VERSION || error("tdfs format version $(meta["format_version"]) (reader is $TDFS_FORMAT_VERSION)")
    load(name) = (t = Arrow.Table(joinpath(dir, name)); NamedTuple(k => collect(getproperty(t, k)) for k in propertynames(t)))
    frames = load("frames.arrow"); slices = load("slices.arrow")
    blocks = open(joinpath(dir, "blocks.bin"), "r") do io
        Mmap.mmap(io, Vector{UInt8}, filesize(io))
    end
    cal = LinearMzCal(Float64(meta["mz_cal_sqrt_intercept"]), Float64(meta["mz_cal_sqrt_slope"]))
    TdfsFile(String(dir), meta, frames, slices, blocks, cal, Int(meta["bin_scale"]), Float64(meta["int_scale"]))
end
n_frames(t::TdfsFile) = length(t.frames.frame_id)
n_slices(t::TdfsFile) = length(t.slices.frame_row)

"Decode frame row i of a tdfs into `blk`."
function read_frame_block!(blk::SliceBlock, codec::BlockCodec, t::TdfsFile, i::Integer)
    off = t.frames.block_offset[i]; bs = t.frames.block_size[i]
    p = pointer(t.blocks) + off
    hdr_size = GC.@preserve t unsafe_load(Ptr{UInt32}(p)); n_words = GC.@preserve t unsafe_load(Ptr{UInt32}(p + 4))
    Int(hdr_size) == bs || error("tdfs frame $i: header block_size $hdr_size != frames table $bs")
    Int(n_words) == t.frames.n_words[i] || error("tdfs frame $i: header n_words $n_words != frames table $(t.frames.n_words[i])")
    decode_block!(blk, codec, view(t.blocks, off + 9:off + bs), Int(n_words))
    blk.n_slices == t.frames.n_slices[i] || error("tdfs frame $i: decoded $(blk.n_slices) slices, table says $(t.frames.n_slices[i])")
    blk
end

"Fixed-point bin -> m/z."
@inline bin_to_mz(t::TdfsFile, b::UInt32) = bin_to_mz(t.mz_cal, Float64(b) / t.bin_scale)
