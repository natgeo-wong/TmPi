function create(
    dt :: Date,
    eroot :: AbstractString = homedir();
    isprecise :: Bool = false
)

    e5ds = ERA5Hourly(dtbeg=dt,dtend=dt,eroot=eroot)

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    downloadERA5(e5ds,[psfc,tsfc,tdew])
    downloadERA5(e5ds,[tair,shum],isprecise)

    calculate(e5ds,isprecise)

end

function downloadERA5(
    e5ds :: ERA5Dataset,
    evar :: Vector{SingleVariable{String}},
)

    ckeys = cdskey()
    dtii  = e5ds.dtbeg

    @info "$(modulelog()) - Using CDSAPI in Julia to download SINGLE-LEVEL $(uppercase(e5ds.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(dtii)."

    fnc = joinpath(e5ds.eroot,"tmpnc-single-$dtii.nc")
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
    
    if !isfile(fnc)
        retrieve("reanalysis-era5-single-levels",e5dkey,fnc,ckeys)
    end

end

function downloadERA5(
    e5ds :: ERA5Dataset,
    evar :: Vector{PressureVariable{String}},
    isprecise :: Bool
)

    ckeys = cdskey()
    dtii  = e5ds.dtbeg

    @info "$(modulelog()) - Using CDSAPI in Julia to download PRESSURE-LEVEL $(uppercase(e5ds.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(dtii)."

    plist = era5Pressures(); plist = plist[plist.>=50]

    if isprecise

        for ip in plist

            fnc = joinpath(e5ds.eroot,"tmpnc-pressure-$ip-$dtii.nc")
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

            if !isfile(fnc)
                retrieve("reanalysis-era5-pressure-levels",e5dkey,fnc,ckeys)
            end

        end

    else

        for evarii in evar

            fnc = joinpath(e5ds.eroot,"tmpnc-pressure-$(evar.varID)-$dtii.nc")
            fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

            e5dkey = Dict(
                "product_type"   => e5ds.ptype,
                "year"           => year(dtii),
                "month"          => month(dtii),
                "day"            => collect(1:31),
                "variable"       => evarii.lname,
                "pressure_level" => plist,
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
            
            if !isfile(fnc)
                retrieve("reanalysis-era5-pressure-levels",e5dkey,fnc,ckeys)
            end

        end

    end

end