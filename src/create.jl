function create(
    tmpi :: ERA5Dataset;
    verbose :: Bool = false,
    keepraw :: Bool = false
)

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    downloadERA5(tmpi,[psfc,tsfc,tdew])
    downloadERA5(tmpi,[tair,shum])

    calculate(tmpi,verbose,keepraw)

end