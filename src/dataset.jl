struct TmPiDefault{ST<:AbstractString} <: TmPiDataset

    ts :: Array{Float64,2}
    td :: Array{Float64,2}
    sp :: Array{Float64,2}
    ta :: Array{Float64,3}
    sh :: Array{Float64,3}

    tm :: Array{Float64,3}
    Pi :: Array{Float64,3}

    lsd :: LandSea
    p   :: Vector{Float64}

    tmp2D :: Array{Int16,2}
    tmp3D :: Array{Int16,3}

    name  :: ST
    ptype :: ST
	sldoi :: ST
    path  :: ST
    emask :: ST

end

struct TmPiPrecise{ST<:AbstractString} <: TmPiDataset

    ts :: Array{Float64,2}
    td :: Array{Float64,2}
    sp :: Array{Float64,2}
    ta :: Array{Float64,3}
    sh :: Array{Float64,3}

    tm :: Array{Float64,3}
    Pi :: Array{Float64,3}

    lsd :: LandSea
    p  :: Vector{Float64}

    tmp2D :: Array{Int16,2}

    lname :: ST
    ptype :: ST
	sldoi :: ST
    path  :: ST
    emask :: ST

end

function TmPiDataset(;
    path  :: AbstractString = homedir(),
    isprecise :: Bool = false,
    ST = String
)

    p = era5Pressures(); np = length(p)
    p = Float64.(p*100)

    @info "$(modulelog()) - Loading Global LandSea dataset (0.25º resolution)"
    lsd  = getLandSea(ERA5Region(GeoRegion("GLB"),resolution=0.25),path=path)
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    @info "$(modulelog()) - Preallocating arrays for downloaded ERA5 datasets"
    ts = zeros(Float64,nlon,nlat)
    td = zeros(Float64,nlon,nlat)
    sp = zeros(Float64,nlon,nlat)
    ta = zeros(Float64,nlon,nlat,np)
    sh = zeros(Float64,nlon,nlat,np)

    @info "$(modulelog()) - Preallocating temporary arrays to load raw data from NetCDF"
    tmp2D = zeros(Int16,nlon,nlat)
    tmp3D = zeros(Int16,nlon,nlat,np)

    @info "$(modulelog()) - Preallocating arrays for final Tm and Pi data"
    tm = zeros(Float64,nlon,nlat,744) # 31*24 = 744
    Pi = zeros(Float64,nlon,nlat,744) # 31*24 = 744

    if isprecise
        return TmPiPrecise{ST}(
            ts, td, sp, ta, sh, tm, Pi, lsd, p, tmp2D,
            "ERA5 Hourly", "reanalysis", "10.24381/cds.adbb2d47",
            path, path
        )
    else
        return TmPiDefault{ST}(
            ts, td, sp, ta, sh, tm, Pi, lsd, p, tmp2D, tmp3D,
            "ERA5 Hourly", "reanalysis", "10.24381/cds.adbb2d47",
            path, path
        )
    end

end