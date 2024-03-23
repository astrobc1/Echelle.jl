export get_inds1d, polyfit1d, interp1d, repair_bad_pix1d, quantile_filter, robust_stddev, robust_mean, weighted_stddev, weighted_quantile

function get_inds1d(coords::Vector{<:CartesianIndex}; dim::Int)
    good = Int[coord.I[dim] for coord ∈ coords]
    return good
end

quantile_filter(x::AbstractArray; window, q::Real=0.5) = mapwindow((xx) -> nanquantile(xx, q), x, window)

function polyfit1d(
        x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, w::Union{AbstractVector{<:Real}, Nothing}=nothing;
        deg::Int, max_iterations::Int=5, nσ::Real=5
    )

    # Make sure we do 1 iteration
    @assert max_iterations > 0 "max_iterations=$max_iterations must be > 0"

    # Initial weights
    w = isnothing(w) ? ones(size(x)) : copy(w)

    # Julia scope
    pfit = nothing

    for _=1:max_iterations

        # Fit
        good = findall(@. isfinite(x) && isfinite(y) && isfinite(w) && (w > 0))
        pfit = Polynomials.fit(ArnoldiFit, x[good], y[good], deg, weights=w[good])

        # Flag
        res = pfit.(x) .- y
        σ = 1.4826 * nanmad(res)
        bad = findall(abs.(res) .> nσ * σ)
        if length(bad) == 0
            break
        end
        w[bad] .= 0

    end

    # good vals used in fit
    good = findall((w .> 0) .&& isfinite.(w))

    # Return
    return pfit, good

end

function weighted_stddev(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}; μ::Union{Real, Nothing}=nothing)
    good = findall(@. isfinite(x) && (w > 0) && isfinite(w))
    if length(good) == 0
        return NaN
    end
    xx = @view x[good]
    ww = w[good]
    ww ./= sum(ww)
    if isnothing(μ)
        μ = nansum(xx .* ww) / nansum(ww)
    end
    dev = xx .- μ
    bias_estimator = 1.0 - sum(ww.^2)
    σ = sqrt(sum(dev .^2 .* ww) / bias_estimator)
    return σ
end

function weighted_quantile(x::AbstractArray{<:Real}, w::AbstractArray{<:Real}; q::Real=0.5)
    good = findall(@. isfinite(x) && (w > 0) && isfinite(w))
    if length(good) > 0
        xx = @view x[good]
        ww = @view w[good]
        return quantile(xx, Weights(ww), q)
    else
        return NaN
    end
end


function robust_stddev(x::AbstractArray{<:Real}, w::Union{AbstractArray{<:Real}, Nothing}=nothing; nσ::Real=4)
    if isnothing(w)
        w = ones(size(x))
    end
    med = weighted_quantile(x, w)
    adevs = abs.(med .- x)
    mad = weighted_quantile(adevs, w)
    good = findall(adevs .< 1.4826 * mad * nσ)
    if length(good) > 1
        return @views weighted_stddev(x[good], w[good])
    else
        return NaN
    end
end


function robust_mean(x::AbstractArray{<:Real}, w::Union{AbstractArray{<:Real}, Nothing}=nothing; nσ::Real=4)
    if isnothing(w)
        w = ones(size(x))
    end
    med = weighted_quantile(x, w)
    adevs = abs.(med .- x)
    mad = weighted_quantile(adevs, w)
    good = findall(adevs .< 1.4826 * mad * nσ)
    if length(good) > 1
        return @views nansum(x[good] .* w[good]) / nansum(w[good])
    else
        return NaN
    end
end

function interp1d(
        x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xnew::AbstractVector{<:Real};
        check_finite::Bool=true, extrapolate::Bool=true
    )

    # Good/bounds
    if check_finite || extrapolate
        good = findall(@. isfinite(x) && isfinite(y))
        xx, yy = @views x[good], y[good]
        xi, xf = xx[1], xx[end]
        yi, yf = yy[1], yy[end]
    else
        xi, xf = x[1], x[end]
        xx, yy = x, y
    end

    # Interpolate
    itp = CubicSpline(yy, xx, extrapolate=true)
    ynew = itp(xnew)

    # Bounds are either constant or NaN
    if extrapolate
        ynew[xnew .< xi] .= yi
        ynew[xnew .> xf] .= yf
    else
        ynew[xnew .< xi] .= NaN
        ynew[xnew .> xf] .= NaN
    end

    # Return
    return ynew

end

function interp1d(
        itp::DataInterpolations.AbstractInterpolation, xnew::AbstractVector{<:Real};
        extrapolate::Bool=true
    )

    # Interpolate
    ynew = itp(xnew)

    # Bounds are either constant or NaN
    xi, xf = itp.t[1], itp.t[end]
    if extrapolate
        ynew[xnew .< xi] .= itp.u[1]
        ynew[xnew .> xf] .= itp.u[end]
    else
        ynew[xnew .< xi] .= NaN
        ynew[xnew .> xf] .= NaN
    end

    # Return
    return ynew

end


function repair_badpix1d(x, y)
    good = findall(isfinite.(y))
    bad = findall(.~isfinite.(y))
    yrep = copy(y)
    yrep[bad] .= @views interp1d(x[good], y[good], x[bad], check_finite, check_bounds)
    return yrep
end