calcTd2e(Td::Float32) = 611.657 * exp((2.5e6/461.5181) * (1/273.16 - 1/Td))

function calcTd2ef(Td::Float32)

    if Td >= 273.16
    	return Float32(611.21 * exp(17.502 * (Td-273.16) / (Td-32.19)))
    elseif Td <= 250.16
    	return Float32(611.21 * exp(22.587 * (Td-273.16) / (Td+0.7)))
    else
        α  = ((Td - 250.16) / (273.16 - 250.16))^2
        ei = 611.21 * exp(22.587 * (Td-273.16) / (Td+0.7))
        ew = 611.21 * exp(17.502 * (Td-273.16) / (Td-32.19))
        return Float32(α * ew + (1-α) * ei)
    end

end

calce2q(e::Float32,p::Float32) = e * 0.6219838793551742 / (p - e * 0.3780161206448258)
calcTm2Pi(Tm::Real) = 10^6 / ((3.739e3 / Tm + 0.221) * 461.5181) / 1000

function calculate(
    tmpi :: TmPiDefault,
    date :: Date,
    verbose :: Bool,
    keepraw :: Bool
)
    
    ndt = daysinmonth(date) * 24
    p = tmpi.p; np = length(p)
    
    @info "$(modulelog()) - Preallocating arrays for numerical integration to calculate Tm"
    ind = zeros(Bool,np+2)
    bot = zeros(Float32,np+2)
    ita = zeros(Float32,np+2)
    ish = zeros(Float32,np+2)
    ipv = Float32.(vcat(0,p,0))

    @info "$(modulelog()) - Opening NetCDF files for Single Level datasets in $(year(date)) $(monthname(date))"
    sds = NCDataset(joinpath(tmpi.path,"tmpnc-single-$date.nc"))
    
    @info "$(modulelog()) - Opening NetCDF files for Pressure Level datasets in $(year(date)) $(monthname(date))"
    tds = NCDataset(joinpath(tmpi.path,"tmpnc-pressure-t-$date.nc"))
    qds = NCDataset(joinpath(tmpi.path,"tmpnc-pressure-q-$date.nc"))
    
    @info "$(modulelog()) - Loading the Global (0.25º Resolution) LandSea Dataset"
    lsd = tmpi.lsd
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    @info "$(modulelog()) - Calculating Tm and Pi for $(year(date)) $(monthname(date))"

    for it in 1 : ndt

        if verbose
            @info "$(modulelog()) - Calculating Tm and Pi for step $it out of $ndt in $(year(date)) $(monthname(date))"
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
            ish[end] = Float32(calce2q(Float32(calcTd2ef(itd)),isp))
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
            tmpi.Pi[ilon,ilat,it] = calcTm2Pi(tmpi.tm[ilon,ilat,it])

        end

    end
    
    close(sds)
    close(tds)
    close(qds)

    @info "$(modulelog()) - Saving Tm and Pi data for $(year(date)) $(monthname(date))"

    save(
        view(tmpi.tm,:,:,1:ndt),date,tmpi,
        SingleVariable("t_qwm"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    save(
        view(tmpi.Pi,:,:,1:ndt),date,tmpi,
        SingleVariable("Pi"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    if !keepraw
        @info "$(modulelog()) - Removing temporary data downloaded from ERA5 to save space"
        rm(joinpath(tmpi.path,"tmpnc-single-$date.nc"),force=true)
        rm(joinpath(tmpi.path,"tmpnc-pressure-t-$date.nc"),force=true)
        rm(joinpath(tmpi.path,"tmpnc-pressure-q-$date.nc"),force=true)
    end
end

function calculate(
    tmpi :: TmPiPrecise,
    date :: Date,
    verbose :: Bool,
    keepraw :: Bool
)
    
    ndt = daysinmonth(date) * 24
    p = tmpi.p; np = length(p)
    
    @info "$(modulelog()) - Preallocating arrays for numerical integration to calculate Tm"
    ind = zeros(Bool,np+2)
    bot = zeros(Float32,np+2)
    ita = zeros(Float32,np+2)
    ish = zeros(Float32,np+2)
    ipv = Float32.(vcat(0,p,0))

    @info "$(modulelog()) - Opening NetCDF files for Single Level datasets in $(year(date)) $(monthname(date))"
    sds = NCDataset(joinpath(tmpi.path,"tmpnc-single-$date.nc"))
    
    @info "$(modulelog()) - Opening NetCDF files for Pressure Level datasets in $(year(date)) $(monthname(date))"
    pds = Vector{NCDataset}(undef,np)
    for ip in 1 : np
        pds[ip] = NCDataset(joinpath(tmpi.path,"tmpnc-pressure-$(p[ip])-$date.nc"))
    end
    
    @info "$(modulelog()) - Loading the Global (0.25º Resolution) LandSea Dataset"
    lsd = tmpi.lsd
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    @info "$(modulelog()) - Calculating Tm and Pi for $(year(date)) $(monthname(date))"

    for it in 1 : ndt

        if verbose
            @info "$(modulelog()) - Calculating Tm and Pi for step $it out of $ndt in $(year(date)) $(monthname(date))"
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
            tmpi.Pi[ilon,ilat,it] = calcTm2Pi(tmpi.tm[ilon,ilat,it])

        end

    end
    
    close(sds)
    for pdsii in pds
        close(pdsii)
    end

    @info "$(modulelog()) - Saving Tm and Pi data for $(year(date)) $(monthname(date))"

    save(
        view(tmpi.tm,:,:,1:ndt),date,tmpi,
        SingleVariable("t_qwm"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    save(
        view(tmpi.Pi,:,:,1:ndt),date,tmpi,
        SingleVariable("Pi"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    if !keepraw
        @info "$(modulelog()) - Removing temporary data downloaded from ERA5 to save space"
        rm(joinpath(tmpi.path,"tmpnc-single-$date.nc"),force=true)
        for ip in 1 : np
           rm(joinpath(tmpi.path,"tmpnc-pressure-$(p[ip])-$date.nc"),force=true)
        end
    end

end