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

# Step 1: IM accumulation. For one slice centred on scan c, accumulate kernel-weighted intensities of the
# window's scans within the kernel reach, per TOF bin, into a dense per-thread accumulator; a bitmap of
# touched bins is then swept to emit the slice's sparse (bin, value, count) list in ascending bin order.

mutable struct SmoothScratch
    n_bins::Int
    acc::Vector{Float64}
    cnt::Vector{UInt16}
    bitmap::Vector{UInt64}
    sp_bin::Vector{Int32}
    sp_val::Vector{Float64}
    sp_cnt::Vector{Int32}
    dense::Vector{Float64}
    dcnt::Vector{Int32}
    tmp::Vector{Float64}       # cap_slice! scratch (dense / dcnt are resized together and must stay paired)
end
SmoothScratch() = SmoothScratch(0, Float64[], UInt16[], UInt64[], Int32[], Float64[], Int32[], Float64[], Int32[], Float64[])

"Grow the dense accumulators to cover bins 0 .. n_bins-1."
function ensure_bins!(sc::SmoothScratch, n_bins::Integer)
    n_bins <= sc.n_bins && return sc
    resize!(sc.acc, n_bins); fill!(view(sc.acc, sc.n_bins + 1:n_bins), 0.0)
    resize!(sc.cnt, n_bins); fill!(view(sc.cnt, sc.n_bins + 1:n_bins), UInt16(0))
    nw = cld(n_bins, 64)
    ow = length(sc.bitmap)
    resize!(sc.bitmap, nw); fill!(view(sc.bitmap, ow + 1:nw), UInt64(0))
    sc.n_bins = n_bins
    sc
end

"""
    im_accumulate!(sc, buf, s0, s1, c, kim, h_im)

Slice at scan `c` of window scans `s0:s1-1`: `acc[bin] += kim[s - c + h_im + 1] * intensity` over scans
`s` in `max(s0, c-h_im) : min(s1-1, c+h_im)`. Emits `sc.sp_*` sorted by bin and leaves the accumulators clean.
"""
function im_accumulate!(sc::SmoothScratch, buf::FrameBuffer, s0::Integer, s1::Integer, c::Integer, kim::Vector{Float64}, h_im::Integer)
    acc = sc.acc; cnt = sc.cnt; bm = sc.bitmap
    tof = buf.tof; it = buf.intensity; ss = buf.scan_start
    lo_bin = typemax(Int); hi_bin = -1
    @inbounds for s in max(s0, c - h_im):min(s1 - 1, c + h_im)
        w = kim[s - c + h_im + 1]
        r = ss[s + 1]:ss[s + 2] - 1
        isempty(r) && continue
        for k in r
            b = Int(tof[k])
            acc[b + 1] += w * it[k]
            cnt[b + 1] += UInt16(1)
            bm[(b >> 6) + 1] |= UInt64(1) << (b & 63)
        end
        lo_bin = min(lo_bin, Int(tof[first(r)])); hi_bin = max(hi_bin, Int(tof[last(r)]))
    end
    empty!(sc.sp_bin); empty!(sc.sp_val); empty!(sc.sp_cnt)
    hi_bin < 0 && return sc
    @inbounds for wi in (lo_bin >> 6) + 1:(hi_bin >> 6) + 1
        word = bm[wi]
        word == 0 && continue
        bm[wi] = 0
        base = (wi - 1) << 6
        while word != 0
            t = trailing_zeros(word)
            b = base + t
            push!(sc.sp_bin, Int32(b)); push!(sc.sp_val, acc[b + 1]); push!(sc.sp_cnt, Int32(cnt[b + 1]))
            acc[b + 1] = 0.0; cnt[b + 1] = UInt16(0)
            word &= word - 1
        end
    end
    sc
end
