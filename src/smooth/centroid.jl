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

# Step 3: centroids of one dense run `sm[1:L]` (bins t0 .. t0+L-1): local maxima with a footprint walked down
# both sides (<= max_half bins), position = intensity-weighted mean bin (:wmean) or 3-point Gaussian apex
# (:gauss), intensity = footprint sum. Culls: footprint sum >= thr, raw entries in the footprint >= min_scans.
#
# Logic identical to the prototype's centroid_dense!, with two mechanical changes: the left footprint edge is
# tracked during the forward scan (`rise`: where the current strictly increasing positive run began) instead of
# walked back from each apex, and `jlo:jhi` restricts where apices are looked for. The m/z stage passes the
# span of the run's input bins: the smoothed profile is monotone in the padding (a sum of Gaussian tails) and
# cannot hold an apex there, while footprints still walk into the padding.

function centroid_dense!(out::FrameSlices, sm::Vector{Float64}, cnt::Vector{Int32}, L::Int, t0::Int, thr::Float64, lp::LevelParams,
                         jlo::Int = 2, jhi::Int = L - 1)
    max_half = lp.max_half; gauss = lp.centroid == :gauss; min_scans = lp.min_scans
    jlo = max(2, jlo); jhi = min(L - 1, jhi)
    jlo > jhi && return out
    # the increasing run containing jlo: walk back at most max_half steps (once per run, not per apex)
    rise = jlo
    @inbounds while rise > 1 && rise > jlo - max_half && sm[rise-1] < sm[rise] && sm[rise-1] > 0; rise -= 1; end
    @inbounds for j in jlo:jhi
        v = sm[j]; vm = sm[j-1]
        j > jlo && !((vm < v) & (vm > 0.0)) && (rise = j)
        ((v > 0.0) & (v > vm) & (v >= sm[j+1])) || continue
        lo = max(rise, j - max_half)
        hi = j; while hi < L && hi < j + max_half && sm[hi+1] < sm[hi] && sm[hi+1] > 0; hi += 1; end
        s = 0.0; ws = 0.0; n_entries = Int32(0)
        for k in lo:hi; s += sm[k]; ws += sm[k] * (t0 + k - 1); n_entries += cnt[k]; end
        s >= thr || continue
        n_entries >= min_scans || continue
        if gauss
            l0 = log(v); lm = log(max(vm, 1e-9)); lp_ = log(max(sm[j+1], 1e-9))
            den = lm - 2l0 + lp_
            off = den < 0 ? clamp(0.5 * (lm - lp_) / den, -1.0, 1.0) : 0.0
            tpos = Float64(t0 + j - 1) + off
        else
            tpos = ws / s
        end
        push_peak!(out, tpos, s)
    end
    out
end
