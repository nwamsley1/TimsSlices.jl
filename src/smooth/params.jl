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
#
# The ion-mobility scale (slice spacing and IM kernel width) is set in 1/K0 units, because the width of a mobility
# peak is a physical width in 1/K0 while the number of scans it spans depends on the method's TIMS ramp (1/K0 range
# / number of scans). Each run's scan counts are derived once, from its own 1/K0-per-scan slope, by
# `resolve_im_scale`; an explicit `stride` / `im_sigma` (in scans) overrides the derivation.

"""
Default slice spacing in 1/K0. Converted per run to `ceil(STRIDE_K0 / (1/K0 per scan))` scans, rounded UP so slices
are never closer than this. 0.0065 gives 8 scans on both ramps the defaults were validated on (timsTOF Ultra / Ultra 2,
1/K0 0.64-1.45 over 936 / 953 scans: 7.51 and 7.65 scans -> 8); the obvious 8 x 0.000865 = 0.0069 would give 9 on
the 953-scan ramp.
"""
const STRIDE_K0 = 0.0065

"""
Default IM Gaussian sigma in 1/K0: 5 scans at 0.000865 1/K0 per scan (the validated setting on the Ultra 2 ramp).
Converted per run to `IM_SIGMA_K0 / (1/K0 per scan)` scans and rounded to 0.01 scan (the kernel takes a fractional
sigma; the rounding only makes the validated ramps reproduce exactly: 0.81/936 per scan gives 4.998 -> 5.00).
"""
const IM_SIGMA_K0 = 0.004325

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
    # IM scale: targets in 1/K0, converted to scans per run (resolve_im_scale). `stride` / `im_sigma` (scans), when
    # given, override the conversion; `nothing` means derive. MS1 follows MS2 unless set.
    stride_k0::Float64 = STRIDE_K0
    im_sigma_k0::Float64 = IM_SIGMA_K0
    stride::Union{Nothing, Int} = nothing
    ms1_stride::Union{Nothing, Int} = stride
    im_sigma::Union{Nothing, Float64} = nothing
    ms1_im_sigma::Union{Nothing, Float64} = im_sigma
    kernel_extent::Float64 = 3.0
    sum_scale::Bool = true
    # m/z kernel + centroid
    mz_sigma::Float64 = 3.0
    mz_sigma_ppm::Float64 = 0.0      # > 0: m/z sigma constant in ppm (bins vary with m/z), replaces mz_sigma
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

"""
    stride_scans(stride_k0, k0_per_scan) -> Int

Slice spacing in scans for a target spacing in 1/K0: `ceil(stride_k0 / k0_per_scan)`, at least 1. A ratio within
1e-9 of an integer is taken as that integer (so float noise cannot push an exact 8 to 9).
"""
stride_scans(stride_k0::Real, k0_per_scan::Real) = max(1, ceil(Int, stride_k0 / abs(k0_per_scan) - 1e-9))

"""
    resolve_im_scale(p, k0_per_scan) -> ConvertParams

`p` with every IM-scale field in scans: `stride` / `ms1_stride` from `stride_k0` (rounded up, `stride_scans`) and
`im_sigma` / `ms1_im_sigma` from `im_sigma_k0` (to 0.01 scan), using the run's 1/K0 per scan (`k0_per_scan`, sign
ignored). Fields already set are kept (explicit overrides); unset MS1 values follow the MS2 ones.

Examples (defaults 0.0065 / 0.004325): 0.000865 1/K0 per scan -> stride 8, sigma 5.00; 0.000850 -> 8, 5.09;
0.001079 (a timsTOF Pro ramp, 0.60-1.60 over 927 scans) -> 7, 4.01.
"""
function resolve_im_scale(p::ConvertParams, k0_per_scan::Real)
    s = abs(Float64(k0_per_scan))
    needs = p.stride === nothing || p.ms1_stride === nothing || p.im_sigma === nothing || p.ms1_im_sigma === nothing
    needs && !(isfinite(s) && s > 0) &&
        throw(ArgumentError("cannot derive the IM scale: 1/K0 per scan is $k0_per_scan; set stride and im_sigma explicitly"))
    stride = something(p.stride, needs ? stride_scans(p.stride_k0, s) : 0)
    im_sigma = something(p.im_sigma, needs ? round(p.im_sigma_k0 / s; digits = 2) : 0.0)
    fields = NamedTuple{fieldnames(ConvertParams)}(Tuple(getfield(p, k) for k in fieldnames(ConvertParams)))
    ConvertParams(; merge(fields, (stride = stride, ms1_stride = something(p.ms1_stride, stride),
                                   im_sigma = im_sigma, ms1_im_sigma = something(p.ms1_im_sigma, im_sigma)))...)
end

"Check `p`. The IM scale may still be unresolved (`nothing`); `resolve_im_scale` fills it in per run."
function validate(p::ConvertParams)
    p.stride_k0 > 0 || throw(ArgumentError("stride_k0 must be > 0"))
    p.im_sigma_k0 >= 0 || throw(ArgumentError("im_sigma_k0 must be >= 0"))
    all(x -> x === nothing || x >= 0, (p.im_sigma, p.ms1_im_sigma)) || throw(ArgumentError("im_sigma must be >= 0"))
    p.kernel_extent > 0 || throw(ArgumentError("kernel_extent must be > 0"))
    all(x -> x === nothing || x >= 1, (p.stride, p.ms1_stride)) || throw(ArgumentError("stride must be >= 1"))
    p.mz_sigma >= 0 || throw(ArgumentError("mz_sigma must be >= 0"))
    p.mz_sigma_ppm >= 0 || throw(ArgumentError("mz_sigma_ppm must be >= 0"))
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

# `p` must be resolved (resolve_im_scale): the IM fields are converted to concrete Float64 / Int here.
level_params(p::ConvertParams, ms1::Bool) = ms1 ?
    LevelParams(p.ms1_im_sigma, p.kernel_extent, p.ms1_stride, p.sum_scale, p.mz_sigma, p.centroid, p.max_half, p.min_scans, p.ms1_max_peaks) :
    LevelParams(p.im_sigma, p.kernel_extent, p.stride, p.sum_scale, p.mz_sigma, p.centroid, p.max_half, p.min_scans, p.max_peaks)

"""
Output base name from the source name and the (resolved) parameters (same convention as the prototype sweeps).
IM sigma is written to two decimals (a sigma derived from 1/K0 is fractional).
"""
function output_name(source::AbstractString, p::ConvertParams)
    fmt(x) = (x = round(x; digits = 2); isinteger(x) ? string(Int(x)) : string(x))
    name = replace(basename(rstrip(source, '/')), r"\.d$" => "")
    mz = p.mz_sigma_ppm > 0 ? "$(fmt(p.mz_sigma_ppm))ppm" : fmt(p.mz_sigma)
    name *= "_cen_s$(fmt(p.im_sigma))_m$(mz)_k$(p.stride)_$(p.centroid)"
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
