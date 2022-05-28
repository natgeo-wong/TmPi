calcTd2e(Td::Float32) = 6.1078 * exp((2.5e6/461.5181) * (1/273.16 - 1/Td))
calce2q(e::Float32,p::Float32) = e * 0.6219838793551742 / (p - e * 0.3780161206448258)
calcTm2Pi(Tm::Real) = 10^6 / ((3.739e3 / Tm + 0.221) * 461.5181) / 1000

function calculate(
    e5ds :: ERA5Dataset,
    tmpi :: TmPiDefault{FT},
    verbose :: Bool
) where FT <: Real
    
    dt = e5ds.dtbeg; ndt = daysinmonth(dt) * 24
    p = tmpi.p; np = length(p)
    
    @info "$(modulelog()) - Preallocating arrays for numerical integration to calculate Tm"
    ind = zeros(Bool,np+2)
    bot = zeros(Float32,np+2)
    ita = zeros(Float32,np+2)
    ish = zeros(Float32,np+2)
    ipv = Float32.(vcat(0,p,0))

    @info "$(modulelog()) - Opening NetCDF files for Single Level datasets"
    sds = NCDataset(joinpath(e5ds.eroot,"tmpnc-single-$dt.nc"))
    
    @info "$(modulelog()) - Opening NetCDF files for Pressure Level datasets"
    tds = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-t-$dt.nc"))
    qds = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-q-$dt.nc"))
    
    p = Float32.(p*100)

    for it in 1 : ndt

        if verbose
            @info "$(modulelog()) - Calculating Tm and Pi for step $it out of $ndt"
        end

        sc = sds["t2m"].attrib["scale_factor"]
        of = sds["t2m"].attrib["add_offset"]
        mv = sds["t2m"].attrib["missing_value"]
        fv = sds["t2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["t2m"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.ts,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["d2m"].attrib["scale_factor"]
        of = sds["d2m"].attrib["add_offset"]
        mv = sds["d2m"].attrib["missing_value"]
        fv = sds["d2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["d2m"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.td,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["sp"].attrib["scale_factor"]
        of = sds["sp"].attrib["add_offset"]
        mv = sds["sp"].attrib["missing_value"]
        fv = sds["sp"].attrib["_FillValue"]
        NCDatasets.load!(sds["sp"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.sp,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = tds["t"].attrib["scale_factor"]
        of = tds["t"].attrib["add_offset"]
        mv = tds["t"].attrib["missing_value"]
        fv = tds["t"].attrib["_FillValue"]
        NCDatasets.load!(tds["t"].var,tmpi.tmp3D,:,:,:,it)
        int2real!(tmpi.ta,tmpi.tmp3D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = qds["q"].attrib["scale_factor"]
        of = qds["q"].attrib["add_offset"]
        mv = qds["q"].attrib["missing_value"]
        fv = qds["q"].attrib["_FillValue"]
        NCDatasets.load!(qds["q"].var,tmpi.tmp3D,:,:,:,it)
        int2real!(tmpi.sh,tmpi.tmp3D,scale=sc,offset=of,mvalue=mv,fvalue=fv)


        for ilat = 1 : nlat, ilon = 1 : nlon

            its = tmpi.ts[ilon,ilat]
            itd = tmpi.td[ilon,ilat]
            isp = tmpi.sp[ilon,ilat]

            for ip = 2 : (np+1)
                ita[ip] = tmpi.ta[ilon,ilat,ip-1]
                ish[ip] = tmpi.sh[ilon,ilat,ip-1]
            end

            ita[end] = its
            ish[end] = Float32(calce2q(Float32(calcTd2e(itd)),isp))
            ipv[end] = isp

            for ip = 2 : (np+2)
                bot[ip] = ish[ip] / ita[ip]
            end

            for ip = 1 : (np+2)
                ind[ip] = ipv[ip] < isp
            end
            ind[end] = true

            top = @view ish[ind]
            btm = @view bot[ind]
            ipp = @view ipv[ind]

            tmpi.tm[ilon,ilat,it] = integrate(ipp,top) / integrate(ipp,btm)
            tmpi.Pi[ilon,ilat,it] = calcTm2Pi(tm[ilon,ilat,it])

        end

    end
    
    close(sds)
    close(tds)
    close(qds)

    @info "$(modulelog()) - Saving Tm and Pi data for $dt"

    save(
        view(tmpi.tm,:,:,1:ndt),dt,e5ds,
        SingleVariable("t_qwm"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    save(
        view(tmpi.Pi,:,:,1:ndt),dt,e5ds,
        SingleVariable("Pi"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

end

function calculate(
    e5ds :: ERA5Dataset,
    tmpi :: TmPiPrecise{FT}
) where FT <: Real
    
    dt = e5ds.dtbeg; ndt = daysinmonth(dt) * 24
    p = tmpi.p; np = length(p)
    
    @info "$(modulelog()) - Preallocating arrays for numerical integration to calculate Tm"
    ind = zeros(Bool,np+2)
    bot = zeros(Float32,np+2)
    ita = zeros(Float32,np+2)
    ish = zeros(Float32,np+2)
    ipv = Float32.(vcat(0,p,0))

    @info "$(modulelog()) - Opening NetCDF files for Single Level datasets"
    sds = NCDataset(joinpath(e5ds.eroot,"tmpnc-single-$dt.nc"))
    
    @info "$(modulelog()) - Opening NetCDF files for Pressure Level datasets"
    pds = Vector{NCDataset}(undef,np)
    for ip in 1 : np
        pds[ip] = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-$(p[ip])-$dt.nc"))
    end
    
    p = Float32.(p*100)

    for it in 1 : ndt

        if verbose
            @info "$(modulelog()) - Calculating Tm and Pi for step $it out of $ndt"
        end

        sc = sds["t2m"].attrib["scale_factor"]
        of = sds["t2m"].attrib["add_offset"]
        mv = sds["t2m"].attrib["missing_value"]
        fv = sds["t2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["t2m"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.ts,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["d2m"].attrib["scale_factor"]
        of = sds["d2m"].attrib["add_offset"]
        mv = sds["d2m"].attrib["missing_value"]
        fv = sds["d2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["d2m"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.td,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["sp"].attrib["scale_factor"]
        of = sds["sp"].attrib["add_offset"]
        mv = sds["sp"].attrib["missing_value"]
        fv = sds["sp"].attrib["_FillValue"]
        NCDatasets.load!(sds["sp"].var,tmpi.tmp2D,:,:,it)
        int2real!(tmpi.sp,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        for ip = 1 : np

            taip = @view tmpi.ta[:,:,ip]
            ship = @view tmpi.sh[:,:,ip]

            sc = pds[ip]["t"].attrib["scale_factor"]
            of = pds[ip]["t"].attrib["add_offset"]
            mv = pds[ip]["t"].attrib["missing_value"]
            fv = pds[ip]["t"].attrib["_FillValue"]
            NCDatasets.load!(pds[ip]["t"].var,tmpi.tmp2D,:,:,it)
            int2real!(taip,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            sc = pds[ip]["q"].attrib["scale_factor"]
            of = pds[ip]["q"].attrib["add_offset"]
            mv = pds[ip]["q"].attrib["missing_value"]
            fv = pds[ip]["q"].attrib["_FillValue"]
            NCDatasets.load!(pds[ip]["q"].var,tmpi.tmp2D,:,:,it)
            int2real!(ship,tmpi.tmp2D,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        end

        for ilat = 1 : nlat, ilon = 1 : nlon

            its = tmpi.ts[ilon,ilat]
            itd = tmpi.td[ilon,ilat]
            isp = tmpi.sp[ilon,ilat]

            for ip = 2 : (np+1)
                ita[ip] = tmpi.ta[ilon,ilat,ip-1]
                ish[ip] = tmpi.sh[ilon,ilat,ip-1]
            end

            ita[end] = its
            ish[end] = Float32(calce2q(Float32(calcTd2e(itd)),isp))
            ipv[end] = isp

            for ip = 2 : (np+2)
                bot[ip] = ish[ip] / ita[ip]
            end

            for ip = 1 : (np+2)
                ind[ip] = ipv[ip] < isp
            end
            ind[end] = true

            top = @view ish[ind]
            btm = @view bot[ind]
            ipp = @view ipv[ind]

            tmpi.tm[ilon,ilat,it] = integrate(ipp,top) / integrate(ipp,btm)
            tmpi.Pi[ilon,ilat,it] = calcTm2Pi(tm[ilon,ilat,it])

        end

    end
    
    close(sds)
    for pdsii in pds
        close(pdsii)
    end

    @info "$(modulelog()) - Saving Tm and Pi data for $dt"

    save(
        view(tmpi.tm,:,:,1:ndt),dt,e5ds,
        SingleVariable("t_qwm"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    save(
        view(tmpi.Pi,:,:,1:ndt),dt,e5ds,
        SingleVariable("Pi"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

end