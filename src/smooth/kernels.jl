"Normalised Gaussian over -h:h with h = ceil(extent * sigma); sigma <= 0 gives [1.0]."
function gauss_kernel(sigma::Real, extent::Real = 3.0)
    sigma <= 0 && return [1.0]
    h = ceil(Int, extent * sigma)
    k = [exp(-(x^2) / (2sigma^2)) for x in -h:h]
    k ./ sum(k)
end
"IM kernel of a level: Gaussian scaled to sum to the stride when `sum_scale`."
im_kernel(lp::LevelParams) = gauss_kernel(lp.im_sigma, lp.kernel_extent) .* (lp.sum_scale ? lp.stride : 1)
mz_kernel(lp::LevelParams) = gauss_kernel(lp.mz_sigma, lp.kernel_extent)
