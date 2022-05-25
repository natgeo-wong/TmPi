function create(
    dt :: Date,
    eroot :: AbstractString = homedir(),
)

    e5ds = ERA5Hourly(dtbeg=dt,dtend=dt,eroot=eroot)

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    download(e5ds,[psfc,tsfc,tdew])
    download(e5ds,[tair,shum])

    calculate(e5ds)

end

function download(
    e5ds :: ERA5Dataset,
    evar :: Vector{SingleVariable},
)

    ckeys = cdskey()
    dtii  = e5ds.dtbeg

    @info "$(modulelog()) - Using CDSAPI in Julia to download $(uppercase(e5ds.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(dtii)."

    fnc = joinpath(e5ds.eroot,"tmpnc-single.nc")
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

    e5dkey = Dict(
        "product_type" => e5ds.ptype,
        "year"         => year(dtii),
        "month"        => month(dtii),
        "day"          => collect(1:31),
        "variable"     => [evarii.lname for evarii in evar],
        "area"         => [90, 0, -90, 360],
        "grid"         => [0.25, 0.25],
        "time"         => [
            "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00", "21:00", "22:00", "23:00",
        ],
        "format"       => "netcdf",
    )
    retrieve("reanalysis-era5-single-levels",e5dkey,fnc,ckeys)

end

function download(
    e5ds :: ERA5Dataset,
    evar :: Vector{PressureVariable},
)

    ckeys = cdskey()
    dtii  = e5ds.dtbeg

    @info "$(modulelog()) - Using CDSAPI in Julia to download $(uppercase(e5ds.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(dtii)."

    for ip in era5Pressures()

        fnc = joinpath(e5ds.eroot,"tmpnc-pressure-$ip.nc")
        fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

        e5dkey = Dict(
            "product_type"   => e5ds.ptype,
            "year"           => year(dtii),
            "month"          => month(dtii),
            "day"            => collect(1:31),
            "variable"       => [evarii.lname for evarii in evar],
            "pressure_level" => ip,
            "area"           => [90, 0, -90, 360],
            "grid"           => [0.25, 0.25],
            "time"           => [
                "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00", "21:00", "22:00", "23:00",
            ],
            "format"         => "netcdf",
        )
        retrieve("reanalysis-era5-pressure-levels",e5dkey,fnc,ckeys)

    end

end