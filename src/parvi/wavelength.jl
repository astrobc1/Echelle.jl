using Glob
using JLD2, FileIO, FITSIO
using Distributed
using AstroTime
using Distributed
using NaNStatistics
using Infiltrator
using PyCall

using EchelleBase
using EchelleSpectralModeling

########################
#### BASIC WLS INFO ####
########################

function SpectralData.get_λsolution_estimate(data::SpecData{:parvi}, sregion::SpecRegion1d)
    return data.data.λ
end

function get_λsolution_estimate_order(order, fiber; pixel_offset=0)
    if fiber == 1
        λ = λpolys_fiber1[order].([1:2048;] .+ pixel_offset)
        return λ
    else
        λ = λpolys_fiber3[order].([1:2048;] .+ pixel_offset)
        return λ
    end
end

function get_λsolution_cheb2d(fname::String; do_orders::AbstractVector{<:Real}, fiber::Int, xrange=[1, 2048], pixel_offset::Int=0, deg_inter_order::Int, deg_intra_order::Int, max_vel_cut=200, debug=false)

    # Load in the fits file and data
    f = FITSIO.FITS(fname)
    if fiber == 1
        data = FITSIO.read(f[1])
    else
        data = FITSIO.read(f[2])
    end

    # Storage arrays
    pixel_centers = Vector{Float64}[]
    amplitudes = Vector{Float64}[]
    peak_integers = Vector{Int}[]
    λ_centers = Vector{Float64}[]
    weights = Vector{Float64}[]
    orders = Vector{Float64}[]

    # Data grid
    nx = 2048

    # Loop over orders and get peaks
    for (i, order) ∈ enumerate(do_orders)

        println("Getting peaks for order $order, fiber $fiber")
        
        # Initial wavelength solution
        λ_estimate = parvi.get_λsolution_estimate_order(order, fiber; pixel_offset=pixel_offset)
        
        # Lfc flux for this order
        oi = order - parvi.echelle_orders[1] + 1
        lfc_flux = data[oi, :, 1]
        
        # lfc_centers_pix, lfc_centers_λ, peak_integers, amplitudes, σs, rms, offsets, slopes
        result = Wavelength.get_peaks(λ_estimate, lfc_flux, parvi.lfc_ν0, parvi.lfc_Δν, σ_guess=[0.9, 1.4, 2.5], μ_bounds=[-1, 1])

        # Store results
        push!(pixel_centers, result[1])
        push!(orders, fill(order, length(result[1])))
        push!(λ_centers, result[2])
        push!(peak_integers, result[3])
        push!(amplitudes, result[4])
        push!(weights, result[4])
    end

    # Fit peaks
    println("Fitting peaks, fiber $(fiber)")
    pixel_centers_flat = collect(Iterators.flatten(pixel_centers))
    orders_flat = collect(Iterators.flatten(orders))
    λ_centers_flat = collect(Iterators.flatten(λ_centers))
    weights_flat = collect(Iterators.flatten(weights))
    peak_integers_flat = collect(Iterators.flatten(peak_integers))
    amplitudes_flat = collect(Iterators.flatten(amplitudes))

    # Set bad peaks to zero
    bad = findall(pixel_centers_flat .< xrange[1] .|| pixel_centers_flat .> xrange[2])
    weights_flat[bad] .= 0

    max_pixel, max_order = nx, parvi.echelle_orders[2]
    coeffs_best, good_peaks = Wavelength.fit_peaks_cheb2d(pixel_centers_flat, orders_flat, λ_centers_flat, weights_flat, max_pixel, max_order, nx, deg_inter_order, deg_intra_order, 3, max_vel_cut)

    chebs_pixels, chebs_orders = Wavelength.get_chebvals(pixel_centers_flat, orders_flat, max_pixel, max_order, deg_intra_order, deg_inter_order)
    model_best = Wavelength.build_λsolution_chebyval2d_flat(chebs_pixels, chebs_orders, coeffs_best, orders_flat)
    residuals = maths.δλ2δv.(λ_centers_flat .- model_best, λ_centers_flat)

    # Debugging
    if debug
        @infiltrate
        #begin
        #using PyPlot
        #scatter(λ_centers_flat, residuals)
        #scatter(λ_centers_flat[good_peaks], residuals[good_peaks])
        #end
        #sqrt(nansum(residuals.^2 / length(residuals))) / sqrt(length(residuals))
        #sqrt(nansum(residuals[good_peaks].^2 / length(residuals[good_peaks]))) / sqrt(length(residuals[good_peaks]))
    end

    # Return
    return coeffs_best, pixel_centers_flat, orders_flat, λ_centers_flat, peak_integers_flat, amplitudes_flat, residuals, good_peaks
end

# LFC info
# 192.1852839999997
lfc_ν0 = 299792458.0 / (1559.91370 * 1E-9) # freq of pump line in Hz.
lfc_Δν = 10E9 # spacing of peaks [Hz]

λpolys_fiber1 = Dict{Int, Polynomial}(
    130 => Polynomial([158372.3437110153, -169.09554851093307, 0.03448497134928927]),
    129 => Polynomial([111242.70157509336, -119.32832012444352, 0.024775122904140665]),
    128 => Polynomial([77319.9485162296, -83.22960131335182, 0.01759877293043629]),
    127 => Polynomial([53176.044625458686, -57.32712204486438, 0.01234990258160693]),
    126 => Polynomial([36198.45284515843, -38.95552005852378, 0.008553584296151443]),
    125 => Polynomial([24414.767950758207, -26.087195202924764, 0.005840792843948004]),
    124 => Polynomial([16350.827250671682, -17.19456462663133, 0.00392741189300218]),
    123 => Polynomial([10916.766118516141, -11.138684261659257, 0.002596851341288826]),
    122 => Polynomial([7316.246769651352, -7.079874312700209, 0.0016857551778042014]),
    121 => Polynomial([4974.770416157323, -4.4065881330266015, 0.0010723390572417594]),
    120 => Polynomial([3483.5875937346723, -2.6793001869383706, 0.000666951321562449]),
    119 => Polynomial([2556.255279872216, -1.5866648103845395, 0.00040450110579618267]),
    118 => Polynomial([1995.358361333678, -0.9116180118170986, 0.0002384426501338046]),
    117 => Polynomial([1667.3227252653387, -0.505464215095493, 0.00013604622988823477]),
    116 => Polynomial([1483.603126383985, -0.26831303171111487, 7.472343190498676e-5]),
    115 => Polynomial([1386.836128297745, -0.1345120436923588, 3.920807157179006e-5]),
    114 => Polynomial([1340.8116583854583, -0.0619641518327604, 1.9424078029271408e-5]),
    113 => Polynomial([1323.340605818848, -0.02442605563944602, 8.898393875872561e-6]),
    112 => Polynomial([1321.284710985964, -0.006061430901160734, 3.600554778794373e-6]),
    111 => Polynomial([1327.172750518455, 0.0023283013025766656, 1.1113467341941746e-6]),
    110 => Polynomial([1336.957455692854, 0.005847930340710206, 4.099438501828187e-8]),
    109 => Polynomial([1348.5741880009512, 0.007176272270673457, -3.670800333159288e-7]),
    108 => Polynomial([1361.048346707343, 0.007622686357497526, -4.975670150444879e-7]),
    107 => Polynomial([1373.9667531027812, 0.007769196126295652, -5.290341419337081e-7]),
    106 => Polynomial([1387.1815441794995, 0.007841218711605938, -5.341433728564627e-7]),
    105 => Polynomial([1400.6558638184824, 0.007907911189757597, -5.366119267427536e-7]),
    104 => Polynomial([1414.3910666027655, 0.007981115958560313, -5.409902647522738e-7]),
    103 => Polynomial([1428.3972132530284, 0.008058081222591617, -5.464236866509022e-7]),
    102 => Polynomial([1442.6840700249277, 0.008135984545167473, -5.517926283003902e-7]),
    101 => Polynomial([1457.2601345852208, 0.008214415160809264, -5.568233403390091e-7]),
    100 => Polynomial([1472.1336883592094, 0.008294167766834698, -5.617656091137616e-7]),
    99 => Polynomial([1487.3136030824905, 0.008375930076163928, -5.668967827743897e-7]),
    98 => Polynomial([1502.8094934285753, 0.008459835027505087, -5.722907736218696e-7]),
    97 => Polynomial([1518.6315095103566, 0.00854567859674007, -5.778572911791951e-7]),
    96 => Polynomial([1534.790133142191, 0.0086333054953405, -5.834932695079237e-7]),
    95 => Polynomial([1551.29615496747, 0.008722805471192741, -5.891924078054319e-7]),
    94 => Polynomial([1568.1608032273941, 0.008814437309612835, -5.950430464316972e-7]),
    93 => Polynomial([1585.3958890001036, 0.008908415519800159, -6.011403931625662e-7]),
    92 => Polynomial([1603.0138528205102, 0.009004766195672786, -6.074979341452176e-7]),
    91 => Polynomial([1621.0277014719372, 0.009103378178057852, -6.140393321165461e-7]),
    90 => Polynomial([1639.450931546454, 0.0092042038598898, -6.206870145623612e-7]),
    89 => Polynomial([1658.2975652204982, 0.009307405952932403, -6.274674601123445e-7]),
    88 => Polynomial([1677.5823274587067, 0.009413227888232805, -6.344901140835397e-7]),
    87 => Polynomial([1697.3208088663648, 0.009521604354596036, -6.417295072369249e-7]),
    86 => Polynomial([1717.529356045263, 0.009632103207386491, -6.488948676537776e-7]),
    85 => Polynomial([1738.2247786677613, 0.009745702250847744, -6.565045816625439e-7])
)


λpolys_fiber3 = Dict{Int, Polynomial}(
    130 => Polynomial([116262.05814509898, -107.79777193015569, 0.0007792282626132959]),
    129 => Polynomial([82217.83369257161, -75.9759397691202, 0.0005996328269413608]),
    128 => Polynomial([57603.28227487466, -52.94915928634331, 0.0004574075287472836]),
    127 => Polynomial([39994.78213030599, -36.46088802840173, 0.00034565050187743807]),
    126 => Polynomial([27540.758259356197, -24.786316019189055, 0.00025856452053535996]),
    125 => Polynomial([18839.37204523135, -16.618847646404898, 0.00019130895145041363]),
    124 => Polynomial([12839.479859006824, -10.978149857206553, 0.00013986746139217187]),
    123 => Polynomial([8761.017202325827, -7.13620896045032, 0.00010093026745157466]),
    122 => Polynomial([6031.496878883733, -4.558330155297405, 7.178977232371458e-5]),
    121 => Polynomial([4235.783296610552, -2.8564518569716624, 5.024848107140877e-5]),
    120 => Polynomial([3076.724863927296, -1.7525351944707044, 3.45381495079559e-5]),
    119 => Polynomial([2344.5969322223254, -1.050131750274724, 2.324916739913134e-5]),
    118 => Polynomial([1893.6330281622575, -0.6125335514139577, 1.5269232127856652e-5]),
    117 => Polynomial([1624.206188785324, -0.3461721729499844, 9.73042027147273e-6]),
    116 => Polynomial([1469.4688354644686, -0.18816206366721716, 5.963815685296689e-6]),
    115 => Polynomial([1385.4723776050844, -0.09708015162702786, 3.4609031404210107e-6]),
    114 => Polynomial([1343.9700043196005, -0.046242554808243574, 1.8409862962568198e-6]),
    113 => Polynomial([1327.2610877735444, -0.018882749750701255, 8.23937763054757e-7]),
    112 => Polynomial([1324.56627695789, -0.004756600564068472, 2.0763718555428352e-7]),
    111 => Polynomial([1329.531504625838, 0.0021991936980603373, -1.5049939038483644e-7]),
    110 => Polynomial([1338.5493720644101, 0.005445004448326299, -3.4844917345978377e-7]),
    109 => Polynomial([1349.66013714818, 0.006873881260852106, -4.5147613408202823e-7]),
    108 => Polynomial([1361.8540459278552, 0.0074709745337994855, -5.014487393634689e-7]),
    107 => Polynomial([1374.6440690239297, 0.007718938162762526, -5.239989375901737e-7]),
    106 => Polynomial([1387.8151035550852, 0.007837089082143407, -5.338972409543047e-7]),
    105 => Polynomial([1401.284074928649, 0.007915704850266152, -5.389778583583993e-7]),
    104 => Polynomial([1415.0266435700319, 0.007987029571897446, -5.42908314607489e-7]),
    103 => Polynomial([1429.0417442355078, 0.008060079922214726, -5.470599763711291e-7]),
    102 => Polynomial([1443.336151145158, 0.008136095178519799, -5.516995210072799e-7]),
    101 => Polynomial([1457.918704054565, 0.00821450199911188, -5.566867804590695e-7]),
    100 => Polynomial([1472.7986301127385, 0.00829475019854931, -5.618317342839317e-7]),
    99 => Polynomial([1487.985290832551, 0.008376637431927198, -5.670329000401198e-7]),
    98 => Polynomial([1503.4882729668795, 0.008460222140724984, -5.722911835418997e-7]),
    97 => Polynomial([1519.3174991705807, 0.008545686777871934, -5.776676397430713e-7]),
    96 => Polynomial([1535.4833143482326, 0.008633227040489974, -5.832308188421313e-7]),
    95 => Polynomial([1551.996556227745, 0.008722975043874797, -5.890197252513548e-7]),
    94 => Polynomial([1568.8686013335346, 0.008814967916174818, -5.950322290457952e-7]),
    93 => Polynomial([1586.1113703912145, 0.008909172482820301, -6.012364104851879e-7]),
    92 => Polynomial([1603.7373006833486, 0.009005551393134689, -6.075942042253552e-7]),
    91 => Polynomial([1621.7593273602904, 0.009104123432951982, -6.140833090342702e-7]),
    90 => Polynomial([1640.19092417084, 0.009204964634064047, -6.207051676952014e-7]),
    89 => Polynomial([1659.046208060018, 0.009308143608295977, -6.274744942949382e-7]),
    88 => Polynomial([1678.3400214898013, 0.00941367610697494, -6.34400001424998e-7]),
    87 => Polynomial([1698.0878536767952, 0.009521645813106869, -6.414874139198141e-7]),
    86 => Polynomial([1718.3056417297012, 0.009632494143791028, -6.488254029001884e-7]),
    85 => Polynomial([1739.0102568955, 0.009746810598849067, -6.568537004840001e-7])
)