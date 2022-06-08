function compile(
    tmpi  :: TmPiDataset,
	evar  :: SingleLevel;
    dtbeg :: TimeType,
    dtend :: TimeType,
    verbose :: Bool = false
)

    yrbeg = year(dtbeg)
    yrend = year(dtend)
    nt    = yrend - yrbeg + 1

    @info "$(modulelog()) - Loading the Global (0.25º Resolution) LandSea Dataset"
    lsd = tmpi.lsd
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    @info "$(modulelog()) - Preallocating arrays ..."

    eavg = zeros(nlon,nlat)
    emax = zeros(nlon,nlat)
    emin = zeros(nlon,nlat)
    erng = zeros(nlon,nlat)
    edhr = zeros(nlon,nlat)
    eitr = zeros(nlon,nlat)
    esea = zeros(nlon,nlat)

    for yr in yrbeg : yrend

        eds = NCDataset(e5danc(tmpi,evar,date))

        eavg += eds["domain_yearly_mean_climatology"][:]
        erng += (eds["domain_yearly_maximum_climatology"][:] .- 
                 eds["domain_yearly_minimum_climatology"][:])
        esea += dropdims(maximum(eds["domain_monthly_mean_climatology"][:],dims=3) .- 
                         minimum(eds["domain_monthly_mean_climatology"][:],dims=3),dims=3)
        eitr += dropdims(mean(eds["domain_monthly_maximum_climatology"][:] .- 
                              eds["domain_monthly_minimum_climatology"][:],dims=3),dims=3)
        edhr += dropdims(maximum(eds["domain_yearly_mean_hourly"][:],dims=3) .- 
                         minimum(eds["domain_yearly_mean_hourly"][:],dims=3),dims=3)

        if yr == yrbeg
            emax += eds["domain_yearly_mean_climatology"][:]
            emin += eds["domain_yearly_mean_climatology"][:]
        else
            emax .= max.(eds["domain_yearly_mean_climatology"][:],emax)
            emin .= min.(eds["domain_yearly_mean_climatology"][:],emax)
        end

        close(eds)

    end

    @info "$(modulelog()) - Calculating yearly mean, and diurnal, seasonal and interannual variability ..."
    eavg = eavg / nt
    edhr = edhr / nt
    eitr = eitr / nt
    esea = esea / nt
    erng = erng / nt .- (esea .+ eitr)
    eian = emax .- emin

    save(
        eavg, edhr, eitr, esea,  erng, eian,
        tmpi, evar, ERA5Region(GeoRegion("GLB"),gres=0.25), lsd
    )

end

function save(
    eavg :: Array{<:Real,2},
    edhr :: Array{<:Real,2},
    eitr :: Array{<:Real,2},
    esea :: Array{<:Real,2},
    erng :: Array{<:Real,2},
    eian :: Array{<:Real,2},
    tmpi :: TmPiDataset,
    evar :: ERA5Variable,
    ereg :: ERA5Region,
    lsd  :: LandSea
)

    @info "$(modulelog()) - Saving analyzed $(tmpi.lname) $(evar.vname) data in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) for $(year(date)) ..."
    fnc = e5danc(tmpi,evar,date)
    fol = dirname(fnc); if !isdir(fol); mkpath(fol) end
    if isfile(fnc)
        @info "$(modulelog()) - Stale NetCDF file $(fnc) detected.  Overwriting ..."
        rm(fnc);
    end
    ds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(modulelog()) with ERA5Reanalysis.jl",
        "comments"    => "ERA5Reanalysis.jl creates NetCDF files in the same format that data is saved on the Climate Data Store"
    ))
    ds.attrib["doi"] = tmpi.sldoi

    ds.dim["longitude"] = length(lsd.lon)
    ds.dim["latitude"]  = length(lsd.lat)
    ds.dim["hour"]  = 24
    ds.dim["month"] = 12

    nclon = defVar(ds,"longitude",Float64,("longitude",),attrib = Dict(
        "units"     => "degrees_east",
        "long_name" => "longitude",
    ))

    nclat = defVar(ds,"latitude",Float64,("latitude",),attrib = Dict(
        "units"     => "degrees_north",
        "long_name" => "latitude",
    ))

    nclon[:] = lsd.lon
    nclat[:] = lsd.lat
    
    attr_var = Dict(
        "long_name"     => evar.lname,
        "full_name"     => evar.vname,
        "units"         => evar.units,
        "_FillValue"    => Int16(-32767),
        "missing_value" => Int16(-32767),
    )

    ## DOMAIN YEARLY CLIMATOLOGY

    scale,offset = ncoffsetscale(view(davg,:,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_mean_climatology",Int16,
        ("longitude","latitude"),attrib=attr_var)
    ncvar.var[:] = real2int16(view(davg,:,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(dstd,:,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_std_climatology",Int16,
        ("longitude","latitude"),attrib=attr_var)
    ncvar.var[:] = real2int16(view(dstd,:,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(dmax,:,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_maximum_climatology",Int16,
        ("longitude","latitude"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmax,:,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(dmin,:,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_minimum_climatology",Int16,
        ("longitude","latitude"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmin,:,:,25,13),scale,offset)

    ## DOMAIN YEARLY DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(davg,:,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_mean_hourly",Int16,
        ("longitude","latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(davg,:,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(dstd,:,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_std_hourly",Int16,
    ("longitude","latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dstd,:,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(dmax,:,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_maximum_hourly",Int16,
        ("longitude","latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmax,:,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(dmin,:,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_yearly_minimum_hourly",Int16,
        ("longitude","latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmin,:,:,1:24,13),scale,offset)

    ## DOMAIN MONTHLY CLIMATOLOGY

    scale,offset = ncoffsetscale(view(davg,:,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_mean_climatology",Int16,
        ("longitude","latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(davg,:,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dstd,:,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_std_climatology",Int16,
        ("longitude","latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dstd,:,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dmax,:,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_maximum_climatology",Int16,
        ("longitude","latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmax,:,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dmin,:,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_minimum_climatology",Int16,
        ("longitude","latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmin,:,:,25,1:12),scale,offset)

    ## DOMAIN MONTHLY DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(davg,:,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_mean_hourly",Int16,
        ("longitude","latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(davg,:,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dstd,:,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_std_hourly",Int16,
        ("longitude","latitude","hour","month"),attrib=attr_var);
        ncvar.var[:] = real2int16(view(dstd,:,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dmax,:,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_maximum_hourly",Int16,
        ("longitude","latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmax,:,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(dmin,:,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"domain_monthly_minimum_hourly",Int16,
        ("longitude","latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(dmin,:,:,1:24,1:12),scale,offset)

    ## DOMAIN YEARLY ZONAL-MEAN CLIMATOLOGY

    scale,offset = ncoffsetscale(view(zavg,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_mean_climatology",Int16,
        ("latitude",),attrib=attr_var)
    ncvar.var[:] = real2int16(view(zavg,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(zstd,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_std_climatology",Int16,
        ("latitude",),attrib=attr_var)
    ncvar.var[:] = real2int16(view(zstd,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(zmax,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_maximum_climatology",Int16,
        ("latitude",),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmax,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(zmin,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_minimum_climatology",Int16,
        ("latitude",),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmin,:,25,13),scale,offset)

    ## DOMAIN YEARLY ZONAL-MEAN DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(zavg,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_mean_hourly",Int16,
        ("latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zavg,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(zstd,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_std_hourly",Int16,
    ("latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zstd,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(zmax,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_maximum_hourly",Int16,
        ("latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmax,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(zmin,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_yearly_minimum_hourly",Int16,
        ("latitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmin,:,1:24,13),scale,offset)

    ## DOMAIN MONTHLY ZONAL-MEAN CLIMATOLOGY

    scale,offset = ncoffsetscale(view(zavg,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_mean_climatology",Int16,
        ("latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zavg,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zstd,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_std_climatology",Int16,
        ("latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zstd,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zmax,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_maximum_climatology",Int16,
        ("latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmax,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zmin,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_minimum_climatology",Int16,
        ("latitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmin,:,25,1:12),scale,offset)

    ## DOMAIN MONTHLY ZONAL-MEAN DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(zavg,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_mean_hourly",Int16,
        ("latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zavg,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zstd,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_std_hourly",Int16,
        ("latitude","hour","month"),attrib=attr_var);
        ncvar.var[:] = real2int16(view(zstd,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zmax,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_maximum_hourly",Int16,
        ("latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmax,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(zmin,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"zonalavg_monthly_minimum_hourly",Int16,
        ("latitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(zmin,:,1:24,1:12),scale,offset)

    ## DOMAIN YEARLY MERIDIONAL-MEAN CLIMATOLOGY

    scale,offset = ncoffsetscale(view(mavg,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_mean_climatology",Int16,
        ("longitude",),attrib=attr_var)
    ncvar.var[:] = real2int16(view(mavg,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(mstd,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_std_climatology",Int16,
        ("longitude",),attrib=attr_var)
    ncvar.var[:] = real2int16(view(mstd,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(mmax,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_maximum_climatology",Int16,
        ("longitude",),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmax,:,25,13),scale,offset)

    scale,offset = ncoffsetscale(view(mmin,:,25,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_minimum_climatology",Int16,
        ("longitude",),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmin,:,25,13),scale,offset)

    ## DOMAIN YEARLY MERIDIONAL-MEAN DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(mavg,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_mean_hourly",Int16,
        ("longitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mavg,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(mstd,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_std_hourly",Int16,
        ("longitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mstd,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(mmax,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_maximum_hourly",Int16,
        ("longitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmax,:,1:24,13),scale,offset)

    scale,offset = ncoffsetscale(view(mmin,:,1:24,13))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_yearly_minimum_hourly",Int16,
        ("longitude","hour"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmin,:,1:24,13),scale,offset)

    ## DOMAIN MONTHLY MERIDIONAL-MEAN CLIMATOLOGY

    scale,offset = ncoffsetscale(view(mavg,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_mean_climatology",Int16,
        ("longitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mavg,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mstd,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_std_climatology",Int16,
        ("longitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mstd,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mmax,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_maximum_climatology",Int16,
        ("longitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmax,:,25,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mmin,:,25,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_minimum_climatology",Int16,
        ("longitude","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmin,:,25,1:12),scale,offset)

    ## DOMAIN MONTHLY MERIDIONAL-MEAN DIURNAL STATISTICS

    scale,offset = ncoffsetscale(view(mavg,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_mean_hourly",Int16,
        ("longitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mavg,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mstd,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_std_hourly",Int16,
        ("longitude","hour","month"),attrib=attr_var);
        ncvar.var[:] = real2int16(view(mstd,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mmax,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_maximum_hourly",Int16,
        ("longitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmax,:,1:24,1:12),scale,offset)

    scale,offset = ncoffsetscale(view(mmin,:,1:24,1:12))
    attr_var["scale_factor"] = scale
    attr_var["add_offset"]   = offset
    ncvar = defVar(ds,"meridionalavg_monthly_minimum_hourly",Int16,
        ("longitude","hour","month"),attrib=attr_var);
    ncvar.var[:] = real2int16(view(mmin,:,1:24,1:12),scale,offset)

    close(ds)

    @info "$(modulelog()) - Analyzed $(uppercase(tmpi.lname)) $(evar.vname) in $(ereg.geo.name) (Horizontal Resolution: $(ereg.gres)) for $(year(date)) has been saved into $(fnc)."

end