function downloadERA5(
    tmpi :: TmPiDataset,
    evar :: Vector{SingleVariable{String}},
)

    ckeys = cdskey()

    @info "$(modulelog()) - Using CDSAPI in Julia to download SINGLE-LEVEL $(uppercase(tmpi.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(tmpi.dates)."

    fnc = joinpath(tmpi.eroot,"tmpnc-single-$(tmpi.dates).nc")
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

    e5dkey = Dict(
        "product_type" => tmpi.ptype,
        "year"         => year(tmpi.dates),
        "month"        => month(tmpi.dates),
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
    tmpi :: TmPiDefault,
    evar :: Vector{PressureVariable{String}}
)

    ckeys = cdskey()

    @info "$(modulelog()) - Using CDSAPI in Julia to download PRESSURE-LEVEL $(uppercase(tmpi.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(tmpi.dates)."

    for evarii in evar

        fnc = joinpath(tmpi.eroot,"tmpnc-pressure-$(evarii.varID)-$(tmpi.dates).nc")
        fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

        e5dkey = Dict(
            "product_type"   => tmpi.ptype,
            "year"           => year(tmpi.dates),
            "month"          => month(tmpi.dates),
            "day"            => collect(1:31),
            "variable"       => evarii.lname,
            "pressure_level" => tmpi.p,
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

function downloadERA5(
    tmpi :: TmPiPrecise,
    evar :: Vector{PressureVariable{String}}
)

    ckeys = cdskey()

    @info "$(modulelog()) - Using CDSAPI in Julia to download PRESSURE-LEVEL $(uppercase(tmpi.lname)) data in the Global Region (Horizontal Resolution: 0.25) for $(tmpi.dates)."

    for ip in tmpi.p

        fnc = joinpath(tmpi.eroot,"tmpnc-pressure-$ip-$(tmpi.dates).nc")
        fol = dirname(fnc); if !isdir(fol); mkpath(fol) end

        e5dkey = Dict(
            "product_type"   => tmpi.ptype,
            "year"           => year(tmpi.dates),
            "month"          => month(tmpi.dates),
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

end