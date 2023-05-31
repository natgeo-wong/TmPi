function readTm(
    dt   :: TimeType;
    path :: AbstractString,
    analysis :: Bool = false
)

    evar = SingleVariable("Tm")
    enc = e5dfnc(path,evar,dt)
    raw = true
    if analysis
        enc = e5danc(path,evar,dt)
        raw = false
    end

    if raw
        if !isfile(enc)
            error("$(modulelog()) - The ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt does not exist at $(enc).  Check if files exist at $(path) or download the files here")
        end
        @info "$(modulelog()) - Opening the ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt"
    end
    if analysis
        if !isfile(enc)
            error("$(modulelog()) - The annually analyzed ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt does not exist at $(enc).  Check if files exist at $(path) or download the files here")
        end
        @info "$(modulelog()) - Opening the annually analyzed ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion GeoRegion during Date $dt"
    end

    return NCDataset(enc)

end

function readPi(
    dt   :: TimeType;
    path :: AbstractString,
    analysis :: Bool = false
)

    evar = SingleVariable("Pi")
    enc = e5dfnc(path,evar,dt)
    raw = true
    if analysis
        enc = e5danc(path,evar,dt)
        raw = false
    end

    if raw
        if !isfile(enc)
            error("$(modulelog()) - The ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt does not exist at $(enc).  Check if files exist at $(path) or download the files here")
        end
        @info "$(modulelog()) - Opening the ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt"
    end
    if analysis
        if !isfile(enc)
            error("$(modulelog()) - The annually analyzed ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion during Date $dt does not exist at $(enc).  Check if files exist at $(path) or download the files here")
        end
        @info "$(modulelog()) - Opening the annually analyzed ERA5 Hourly Dataset for $(evar.name) in the Global (0.25º Resolution) GeoRegion GeoRegion during Date $dt"
    end

    return NCDataset(enc)

end