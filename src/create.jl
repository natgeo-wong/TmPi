function create(
    tmpi :: TmPiDataset,
    e5ds :: ERA5Dataset;
    dt   :: Date,
    verbose :: Bool = false,
    keepraw :: Bool = false
)

    e5ds = ERA5Hourly(dtbeg=dt,dtend=dt,eroot=joinpath(e5ds.eroot,".."))

    psfc = SingleVariable("sp")
    tsfc = SingleVariable("t2m")
    tdew = SingleVariable("d2m")
    tair = PressureVariable("t",hPa=1)
    shum = PressureVariable("q",hPa=1)

    downloadERA5(e5ds,[psfc,tsfc,tdew])
    downloadERA5(e5ds,[tair,shum],tmpi)

    calculate(e5ds,tmpi,verbose,keepraw)

end