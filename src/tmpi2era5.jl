"""
These set of scripts are meant to convert the raw TmPi.jl output into the ERA5Reanalysis.jl filesystem so that the full range of ERA5Reanalysis.jl functionality can be used without needing to duplicate all the functionality of ERA5Reanalysis.jl

Note: you cannot directly move, because TmPi.jl saves data in Float64. There is a need to convert to Int16 format.  Must copy instead.
"""

function tmpisave2era5(
    e5ds :: ERA5Hourly;
    path :: AbstractString,
)

    evar_Tm = SingleVariable("Tm")
    evar_Pi = SingleVariable("Pi")
    egeo    = ERA5Region("GLB",gres=0.25)
    lsd     = getLandSea(e5ds,egeo)

    dtbeg = e5ds.start
    dtend = e5ds.stop

    data = zeros(length(lsd.lon),length(lsd.lat),744)

    for idt in dtbeg : Month(1) : dtend

        nhr = daysinmonth(idt) * 24
        dataii = @view data[:,:,1:nhr]
        ds = readTm(idt,path=path)
        NCDatasets.load!(ds["Tm"].var,dataii,:,:,:)
        close(ds)
        ERA5Reanalysis.save(dataii,idt,e5ds,evar_Tm,egeo,lsd)
        ds = readPi(idt,path=path)
        NCDatasets.load!(ds["Pi"].var,dataii,:,:,:)
        close(ds)
        ERA5Reanalysis.save(dataii,idt,e5ds,evar_Pi,egeo,lsd)

    end

end

function tmpimove2era5(
    e5ds :: ERA5Hourly;
    path :: AbstractString,
    force :: Bool = false
)

    evar_Tm = SingleVariable("Tm")
    evar_Pi = SingleVariable("Pi")
    egeo    = ERA5Region("GLB",gres=0.25)

    dtbeg = e5ds.start
    dtend = e5ds.stop

    for idt in dtbeg : Month(1) : dtend

        Tmold = e5dfnc(path,evar_Tm,idt)
        Tmnew = ERA5Reanalysis.e5dfnc(e5ds,evar_Tm,egeo,idt)
        newfol = dirname(Tmnew)
        if !isdir(newfol); mkpath(newfol) end
        if isfile(Tmold); mv(Tmold,Tmnew,force=force) end

        Piold = e5dfnc(path,evar_Pi,idt)
        Pinew = ERA5Reanalysis.e5dfnc(e5ds,evar_Pi,egeo,idt)
        newfol = dirname(Pinew)
        if !isdir(newfol); mkpath(newfol) end
        if isfile(Piold); mv(Piold,Pinew,force=force) end

    end

end

function era5move2tmpi(
    e5ds :: ERA5Hourly;
    path :: AbstractString,
    force :: Bool = false
)

    evar_Tm = SingleVariable("Tm")
    evar_Pi = SingleVariable("Pi")
    egeo    = ERA5Region("GLB",gres=0.25)

    dtbeg = e5ds.start
    dtend = e5ds.stop

    for idt in dtbeg : Month(1) : dtend

        Tmnew = e5dfnc(path,evar_Tm,idt)
        Tmold = ERA5Reanalysis.e5dfnc(e5ds,evar_Tm,egeo,idt)
        newfol = dirname(Tmnew)
        if !isdir(newfol); mkpath(newfol) end
        if isfile(Tmold); mv(Tmold,Tmnew,force=force) end

        Pinew = e5dfnc(path,evar_Pi,idt)
        Piold = ERA5Reanalysis.e5dfnc(e5ds,evar_Pi,egeo,idt)
        newfol = dirname(Pinew)
        if !isdir(newfol); mkpath(newfol) end
        if isfile(Piold); mv(Piold,Pinew,force=force) end

    end

end