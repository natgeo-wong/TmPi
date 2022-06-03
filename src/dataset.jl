struct TmPiDefault{ST<:AbstractString} <: TmPiDataset

    ts :: Array{Float32,2}
    td :: Array{Float32,2}
    sp :: Array{Float32,2}
    ta :: Array{Float32,3}
    sh :: Array{Float32,3}

    tm :: Array{Float32,3}
    Pi :: Array{Float32,3}

    p  :: Vector{Float32}

    tmp2D :: Array{Int16,2}
    tmp3D :: Array{Int16,3}

    lname :: ST
    ptype :: ST
	sldoi :: ST
    eroot :: ST
    emask :: ST

end

struct TmPiPrecise{ST<:AbstractString} <: TmPiDataset

    ts :: Array{Float32,2}
    td :: Array{Float32,2}
    sp :: Array{Float32,2}
    ta :: Array{Float32,3}
    sh :: Array{Float32,3}

    tm :: Array{Float32,3}
    Pi :: Array{Float32,3}

    p  :: Vector{Float32}

    tmp2D :: Array{Int16,2}

    lname :: ST
    ptype :: ST
	sldoi :: ST
    eroot :: ST
    emask :: ST

end

function TmPiDataset(;
    eroot :: AbstractString = homedir(),
    isprecise :: Bool = false,
    FT = Float32,
    ST = String
)

    p = era5Pressures(); p = p[p.>=50]; np = length(p)
    p = Float32.(p*100)

    @info "$(modulelog()) - Loading Global LandSea dataset (0.25º resolution)"
    lsd  = getLandSea(ERA5Region(GeoRegion("GLB"),gres=0.25),eroot=eroot)
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
        return TmPiPrecise{FT,ST,DT}(
            ts, td, sp, ta, sh, tm, Pi, p, tmp2D,
            "ERA5 Hourly", "reanalysis", "10.24381/cds.adbb2d47",
            eroot, eroot
        )
    else
        return TmPiDefault{FT,ST,DT}(
            ts, td, sp, ta, sh, tm, Pi, p, tmp2D, tmp3D,
            "ERA5 Hourly", "reanalysis", "10.24381/cds.adbb2d47",
            eroot, eroot
        )
    end

end