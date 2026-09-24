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

# All converter knobs. Per-level (MS1 / MS2) values are resolved by `level_params`.

"Parameters of the smoothing pipeline for one MS level."
struct LevelParams
    im_sigma::Float64       # scans
    kernel_extent::Float64  # kernel half-width in sigmas
    stride::Int             # slice spacing in scans
    sum_scale::Bool         # IM kernel sums to `stride` (else to 1)
    mz_sigma::Float64       # bins
    centroid::Symbol        # :wmean, :gauss, :none
    max_half::Int           # footprint walk limit (bins)
    min_scans::Int          # persistence cull (1 = off)
    max_peaks::Int          # keep only the N most intense centroids of a slice (0 = off)
end

Base.@kwdef struct ConvertParams
    # IM kernel
    im_sigma::Float64 = 5.0
    ms1_im_sigma::Float64 = im_sigma
    kernel_extent::Float64 = 3.0
    stride::Int = 8
    ms1_stride::Int = stride
    sum_scale::Bool = true
    # m/z kernel + centroid
    mz_sigma::Float64 = 3.0
    centroid::Symbol = :wmean
    max_half::Int = max(4, ceil(Int, 4 * mz_sigma))
    # culls
    min_scans::Int = 1
    max_peaks::Int = 1500            # per slice, MS2: keep the N most intense centroids (0 = off)
    ms1_max_peaks::Int = 0           # per slice, MS1 (0 = off)
    # format
    bin_scale::Int = 1
    int_scale::Float64 = 1.0
    zstd_level::Int = 3
    format::Symbol = :tdfs          # :tdfs, :arrow, :both
    # run
    batch_frames::Int = 0            # results in flight between workers and the writer; 0 = 16 * threads
    frames::Union{Nothing, Vector{Int}} = nothing   # subset of frame rows (testing / benchmarks)
end

function validate(p::ConvertParams)
    p.im_sigma >= 0 && p.ms1_im_sigma >= 0 || throw(ArgumentError("im_sigma must be >= 0"))
    p.kernel_extent > 0 || throw(ArgumentError("kernel_extent must be > 0"))
    p.stride >= 1 && p.ms1_stride >= 1 || throw(ArgumentError("stride must be >= 1"))
    p.mz_sigma >= 0 || throw(ArgumentError("mz_sigma must be >= 0"))
    p.centroid in (:wmean, :gauss, :none) || throw(ArgumentError("centroid must be :wmean, :gauss or :none"))
    p.max_half >= 1 || throw(ArgumentError("max_half must be >= 1"))
    p.min_scans >= 1 || throw(ArgumentError("min_scans must be >= 1"))
    p.max_peaks >= 0 && p.ms1_max_peaks >= 0 || throw(ArgumentError("max_peaks must be >= 0"))
    p.bin_scale >= 1 || throw(ArgumentError("bin_scale must be >= 1"))
    p.int_scale > 0 || throw(ArgumentError("int_scale must be > 0"))
    1 <= p.zstd_level <= 22 || throw(ArgumentError("zstd_level must be in 1:22"))
    p.format in (:tdfs, :arrow, :both) || throw(ArgumentError("format must be :tdfs, :arrow or :both"))
    p
end

level_params(p::ConvertParams, ms1::Bool) = ms1 ?
    LevelParams(p.ms1_im_sigma, p.kernel_extent, p.ms1_stride, p.sum_scale, p.mz_sigma, p.centroid, p.max_half, p.min_scans, p.ms1_max_peaks) :
    LevelParams(p.im_sigma, p.kernel_extent, p.stride, p.sum_scale, p.mz_sigma, p.centroid, p.max_half, p.min_scans, p.max_peaks)

"Output base name from the source name and the parameters (same convention as the prototype sweeps)."
function output_name(source::AbstractString, p::ConvertParams)
    fmt(x) = isinteger(x) ? string(Int(x)) : string(x)
    name = replace(basename(rstrip(source, '/')), r"\.d$" => "")
    name *= "_cen_s$(fmt(p.im_sigma))_m$(fmt(p.mz_sigma))_k$(p.stride)_$(p.centroid)"
    p.sum_scale && (name *= "_sum")
    p.min_scans > 1 && (name *= "_n$(p.min_scans)")
    p.ms1_im_sigma != p.im_sigma && (name *= "_ms1s$(fmt(p.ms1_im_sigma))")
    p.ms1_stride != p.stride && (name *= "_ms1k$(p.ms1_stride)")
    p.max_peaks > 0 && (name *= "_top$(p.max_peaks)")
    p.ms1_max_peaks > 0 && (name *= "_ms1top$(p.ms1_max_peaks)")
    p.bin_scale != 1 && (name *= "_b$(p.bin_scale)")
    p.int_scale != 1 && (name *= "_i$(fmt(p.int_scale))")
    name
end

function Base.Dict(p::ConvertParams)
    Dict{String, Any}(string(k) => (v = getfield(p, k); v isa Symbol ? string(v) : v) for k in fieldnames(ConvertParams))
end
