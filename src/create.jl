struct TmPiDefault{FT<:Real} <: TmPiDataset

    ts :: Array{Float32,2}
    td :: Array{Float32,2}
    sp :: Array{Float32,2}
    ta :: Array{Float32,3}
    sh :: Array{Float32,3}

    tmp2D :: Array{Int16,2}
    tmp3D :: Array{Int16,3}

    tm :: Array{Float32,3}
    Pi :: Array{Float32,3}

    lsd :: LandSea{FT}
    p   :: Vector{Int16}

end

struct TmPiPrecise{FT<:Real} <: TmPiDataset

    ts :: Array{Float32,2}
    td :: Array{Float32,2}
    sp :: Array{Float32,2}
    ta :: Array{Float32,3}
    sh :: Array{Float32,3}

    tmp2D :: Array{Int16,2}

    tm :: Array{Float32,3}
    Pi :: Array{Float32,3}

    lsd :: LandSea{FT}
    p   :: Vector{Int16}

end

function create(
    tmpi  :: TmPiDataset,
    dt    :: Date;
    eroot :: AbstractString = homedir(),
    verbose :: Bool = false
)

    e5ds = ERA5Hourly(dtbeg=dt,dtend=dt,eroot=eroot)

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    downloadERA5(e5ds,[psfc,tsfc,tdew])
    downloadERA5(e5ds,[tair,shum],tmpi)

    calculate(e5ds,tmpi,verbose)

end

function TmPiDataset(;
    isprecise :: Bool,
    FT = Float32
)

    p = era5Pressures(); p = p[p.>=50]

    @info "$(modulelog()) - Loading Global LandSea dataset (0.25º resolution)"
    lsd  = getLandSea(e5ds,ERA5Region(GeoRegion("GLB"),gres=0.25))
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    @info "$(modulelog()) - Preallocating arrays for downloaded ERA5 datasets"
    ts = zeros(Float32,nlon,nlat)
    td = zeros(Float32,nlon,nlat)
    sp = zeros(Float32,nlon,nlat)
    ta = zeros(Float32,nlon,nlat,np)
    sh = zeros(Float32,nlon,nlat,np)

    @info "$(modulelog()) - Preallocating temporary arrays to load raw data from NetCDF"
    tmp2D = zeros(Int16,nlon,nlat)
    tmp3D = zeros(Int16,nlon,nlat,np)

    @info "$(modulelog()) - Preallocating arrays for final Tm and Pi data"
    tm = zeros(Float32,nlon,nlat,744) # 31*24 = 744
    Pi = zeros(Float32,nlon,nlat,744) # 31*24 = 744

    if isprecise
        return TmPiPrecise{FT}(
            ts, td, sp, ta, sh,
            tmp2D,
            tm, Pi, lsd, p
        )
    else
        return TmPiDefault{FT}(
            ts, td, sp, ta, sh,
            tmp2D, tmp3D,
            tm, Pi, lsd, p
        )
    end

end