function create(
    tmpi :: ERA5Dataset;
    date :: Date,
    verbose :: Bool = false,
    keepraw :: Bool = false
)

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    download(tmpi,date,[psfc,tsfc,tdew])
    download(tmpi,date,[tair,shum])

    calculate(tmpi,date,verbose,keepraw)

end