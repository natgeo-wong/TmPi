"""
These set of scripts are meant to convert the raw TmPi.jl output into the ERA5Reanalysis.jl filesystem so that the full range of ERA5Reanalysis.jl functionality can be used without needing to duplicate all the functionality of ERA5Reanalysis.jl
"""

function tmpisave2era5(
    e5ds :: ERA5Hourly;
    path :: AbstractString
)

    evar_Tm = SingleVariable("Tm")
    evar_Pi = SingleVariable("Pi")
    egeo    = ERA5Region("GLB",gres=0.25)
    lsd     = getLandSea(e5ds,egeo)

    dtbeg = e5ds.start
    dtend = e5ds.stop

    for idt in dtbeg : Month(1) : dtend

        dsTm = readTm(idt,path=path); Tm = dsTm["Tm"][:]; close(dsTm)
        dsPi = readPi(idt,path=path); Pi = dsPi["Pi"][:]; close(dsPi)

        ERA5Reanalysis.save(Tm,idt,e5ds,evar_Tm,egeo,lsd)
        ERA5Reanalysis.save(Pi,idt,e5ds,evar_Pi,egeo,lsd)

    end

end

function tmpimove2era5(
    e5ds :: ERA5Hourly;
    path :: AbstractString
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
        if isfile(Tmold); mv(Tmold,Tmnew) end

        Piold = e5dfnc(path,evar_Pi,idt)
        Pinew = ERA5Reanalysis.e5dfnc(e5ds,evar_Pi,egeo,idt)
        newfol = dirname(Pinew)
        if !isdir(newfol); mkpath(newfol) end
        if isfile(Piold); mv(Piold,Pinew) end

    end

end