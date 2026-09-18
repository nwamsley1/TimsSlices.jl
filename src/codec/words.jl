# The slice word stream of one frame block (our format, see docs/format.md):
#   word 1            n_slices
#   words 2..n+1      peak count of slice j
#   then per slice, per peak: (bin_delta, intensity); bin_delta is the fixed-point bin minus the previous
#   peak's bin in the slice, the accumulator starting at 0xFFFFFFFF (so the first delta is bin + 1).

"Float centroids of one frame before quantisation: per slice, (bin position, intensity) pairs."
mutable struct FrameSlices
    n_slices::Int
    scan::Vector{Int32}      # IM scan of each slice
    window::Vector{Int32}    # window index within the frame (1-based; 1 for MS1)
    ptr::Vector{Int32}       # n_slices + 1 entries; peaks of slice j are ptr[j]:ptr[j+1]-1
    pos::Vector{Float64}     # bin position (fractional)
    val::Vector{Float64}     # intensity
end
FrameSlices() = FrameSlices(0, Int32[], Int32[], Int32[1], Float64[], Float64[])

function reset!(fs::FrameSlices)
    fs.n_slices = 0
    empty!(fs.scan); empty!(fs.window); empty!(fs.pos); empty!(fs.val)
    resize!(fs.ptr, 1); fs.ptr[1] = 1
    fs
end
n_peaks(fs::FrameSlices) = length(fs.pos)
@inline function push_peak!(fs::FrameSlices, pos::Float64, val::Float64)
    push!(fs.pos, pos); push!(fs.val, val); nothing
end
"Close the current slice if it received any peaks since the last `end_slice!`; returns true if kept."
function end_slice!(fs::FrameSlices, scan::Integer, window::Integer)
    length(fs.pos) == fs.ptr[end] - 1 && return false
    fs.n_slices += 1
    push!(fs.scan, Int32(scan)); push!(fs.window, Int32(window)); push!(fs.ptr, Int32(length(fs.pos) + 1))
    true
end
"Drop the peaks of the slice currently being filled (those after the last `end_slice!`)."
function discard_slice!(fs::FrameSlices)
    resize!(fs.pos, fs.ptr[end] - 1); resize!(fs.val, fs.ptr[end] - 1); fs
end

"Quantised slices of one frame: exactly what a block holds."
mutable struct SliceBlock
    n_slices::Int
    ptr::Vector{Int32}
    bin::Vector{UInt32}
    intensity::Vector{UInt32}
end
SliceBlock() = SliceBlock(0, Int32[1], UInt32[], UInt32[])
n_peaks(b::SliceBlock) = length(b.bin)
function reset!(b::SliceBlock)
    b.n_slices = 0; resize!(b.ptr, 1); b.ptr[1] = 1; empty!(b.bin); empty!(b.intensity); b
end
slice_range(b::SliceBlock, j::Integer) = (b.ptr[j]):(b.ptr[j + 1] - 1)

"""
    quantize!(blk, fs, bin_scale, int_scale) -> blk

Round bin positions to `round(bin_scale * pos)` and intensities to `round(int_scale * val)`, merging peaks that
land on the same fixed-point bin (intensities summed) and dropping peaks whose rounded intensity is 0.
Slices keep their identity (a slice may become empty). Peaks within a slice must be sorted by position.
"""
function quantize!(blk::SliceBlock, fs::FrameSlices, bin_scale::Integer, int_scale::Real)
    reset!(blk)
    blk.n_slices = fs.n_slices
    resize!(blk.ptr, fs.n_slices + 1)
    sizehint!(blk.bin, n_peaks(fs)); sizehint!(blk.intensity, n_peaks(fs))
    @inbounds for j in 1:fs.n_slices
        blk.ptr[j] = Int32(length(blk.bin) + 1)
        last_bin = typemax(UInt32); acc = UInt64(0)
        for k in fs.ptr[j]:fs.ptr[j + 1] - 1
            b = round(UInt32, bin_scale * fs.pos[k])
            v = UInt64(round(Int, int_scale * fs.val[k]))
            b < last_bin && last_bin != typemax(UInt32) && throw(ArgumentError("slice $j peaks not sorted by position"))
            if b == last_bin
                acc += v
                blk.intensity[end] = UInt32(min(acc, UInt64(typemax(UInt32))))
            else
                if v > 0
                    push!(blk.bin, b); push!(blk.intensity, UInt32(min(v, UInt64(typemax(UInt32)))))
                    last_bin = b; acc = v
                end
            end
        end
    end
    blk.ptr[end] = Int32(length(blk.bin) + 1)
    blk
end

n_words(blk::SliceBlock) = 1 + blk.n_slices + 2 * n_peaks(blk)

"Serialise `blk` into `words` (resized); returns the word count."
function encode_words!(words::Vector{UInt32}, blk::SliceBlock)
    n = n_words(blk)
    length(words) < n && resize!(words, n)
    @inbounds begin
        words[1] = UInt32(blk.n_slices)
        for j in 1:blk.n_slices
            words[1 + j] = UInt32(blk.ptr[j + 1] - blk.ptr[j])
        end
        pos = 2 + blk.n_slices
        for j in 1:blk.n_slices
            acc = typemax(UInt32)
            for k in blk.ptr[j]:blk.ptr[j + 1] - 1
                b = blk.bin[k]
                words[pos] = b - acc; words[pos + 1] = blk.intensity[k]   # wrapping subtraction
                acc = b; pos += 2
            end
        end
    end
    n
end

"Inverse of `encode_words!`: fills `blk` from `words[1:n]`."
function decode_words!(blk::SliceBlock, words::AbstractVector{UInt32}, n::Integer)
    reset!(blk)
    n >= 1 || throw(ArgumentError("empty word stream"))
    ns = Int(words[1])
    n >= 1 + ns || throw(ArgumentError("word stream too short for $ns slice counts"))
    blk.n_slices = ns
    resize!(blk.ptr, ns + 1)
    total = 0
    @inbounds for j in 1:ns
        blk.ptr[j] = Int32(total + 1); total += Int(words[1 + j])
    end
    blk.ptr[end] = Int32(total + 1)
    n == 1 + ns + 2total || throw(ArgumentError("word stream has $n words, expected $(1 + ns + 2total)"))
    resize!(blk.bin, total); resize!(blk.intensity, total)
    pos = 2 + ns
    @inbounds for j in 1:ns
        acc = typemax(UInt32)
        for k in blk.ptr[j]:blk.ptr[j + 1] - 1
            acc += words[pos]; blk.bin[k] = acc; blk.intensity[k] = words[pos + 1]; pos += 2
        end
    end
    blk
end

"Scratch for encoding / decoding blocks: word buffer, plane buffer, zstd buffer and contexts."
mutable struct BlockCodec
    words::Vector{UInt32}
    planes::Vector{UInt8}
    zbuf::Vector{UInt8}
    zstd::ZstdCtx
end
BlockCodec() = BlockCodec(UInt32[], UInt8[], UInt8[], ZstdCtx())

"""
    encode_block!(codec, blk, level) -> (n_words, nbytes)

Word stream -> byte planes -> zstd into `codec.zbuf[1:nbytes]`.
"""
function encode_block!(c::BlockCodec, blk::SliceBlock, level::Integer)
    n = encode_words!(c.words, blk)
    transpose!(c.planes, c.words, n)
    nb = zstd_compress!(c.zbuf, c.zstd, c.planes, 4n, level)
    n, nb
end

"""
    decode_block!(blk, codec, payload, n_words) -> blk

zstd payload -> byte planes -> word stream -> `blk`. `n_words` is the word count stored in the block header.
"""
function decode_block!(blk::SliceBlock, c::BlockCodec, payload::AbstractVector{UInt8}, n_words::Integer)
    zstd_decompress!(c.planes, c.zstd, payload, 4n_words)
    untranspose!(c.words, c.planes, n_words)
    decode_words!(blk, c.words, n_words)
end
