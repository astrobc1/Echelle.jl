export interp1d, repair_bad_pix1d

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