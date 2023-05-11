function e5dfnc(
    tmpi :: TmPiDataset,
	evar :: SingleLevel,
    date :: TimeType
)

    dts = yr2str(date)
    fol = joinpath(tmpi.path,evar.varID,dts)
    fnc = evar.varID * "-" * yrmo2str(date) * ".nc"
    return joinpath(fol,fnc)

end

function e5dfnc(
    path :: AbstractString,
	evar :: SingleLevel,
    date :: TimeType
)

    dts = yr2str(date)
    fol = joinpath(path,evar.varID,dts)
    fnc = evar.varID * "-" * yrmo2str(date) * ".nc"
    return joinpath(fol,fnc)

end

function e5danc(
    tmpi :: TmPiDataset,
	evar :: SingleLevel,
    date :: TimeType
)

    fnc = evar.varID * "-" * yr2str(date) * ".nc"
    return joinpath(tmpi.path,evar.varID,fnc)

end

function e5danc(
    path :: AbstractString,
	evar :: SingleLevel,
    date :: TimeType
)

    fnc = evar.varID * "-" * yr2str(date) * ".nc"
    return joinpath(path,evar.varID,fnc)

end

function e5dcnc(
    tmpi :: TmPiDataset,
	evar :: SingleLevel,
    dtbeg :: TimeType,
    dtend :: TimeType
)

    fnc = evar.varID * "-$(yr2str(dtbeg))_$(yr2str(dtend)).nc"
    return joinpath(tmpi.path,evar.varID,fnc)

end

function save(
    data :: AbstractArray{<:Real,3},
    date :: Date,
    tmpi :: TmPiDataset,
    evar :: ERA5Variable,
    ereg :: ERA5Region,
    lsd  :: LandSea
)

    @info "$(modulelog()) - Saving raw $(tmpi.lname) $(evar.vname) data in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) for $(year(date)) $(Dates.monthname(date)) ..."

    fnc = e5dfnc(tmpi,evar,date)
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end
    if isfile(fnc)
        @info "$(modulelog()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) with TmPi.jl",
        "comments"    => "TmPi.jl creates NetCDF files in the same format that data is saved on the Climate Data Store"
    ))
    ds.attrib["doi"] = tmpi.sldoi

    nhr = 24 * daysinmonth(date)

    ds.dim["longitude"] = length(lsd.lon);
    ds.dim["latitude"]  = length(lsd.lat);
    ds.dim["time"] = nhr

    nclon = defVar(ds,"longitude",Float32,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float32,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nctime = defVar(ds,"time",Int32,("time",),attrib = Dict(
        "units"     => "hours since $(date) 00:00:00.0",
        "long_name" => "time",
        "calendar"  => "gregorian",
    ))

    ncvar = defVar(ds,evar.varID,Float64,("longitude","latitude","time"),attrib = Dict(
        "long_name"     => evar.lname,
        "full_name"     => evar.vname,
        "units"         => evar.units,
        # "scale_factor"  => scale,
        # "add_offset"    => offset,
        # "_FillValue"    => Int16(-32767),
        # "missing_value" => Int16(-32767),
    ))

    nclon[:]  = lsd.lon
    nclat[:]  = lsd.lat
    nctime[:] = collect(1:nhr) .- 1
    ncvar[:]  = data

    # if iszero(sum(isnan.(data)))
    #       ncvar[:] = data
    # else; ncvar.var[:] = real2int16(data,scale,offset)
    # end

    close(ds)

    @info "$(modulelog()) - Raw $(uppercase(tmpi.lname)) $(evar.vname) in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) for $(year(date)) $(Dates.monthname(date)) has been saved into $(fnc)."

end