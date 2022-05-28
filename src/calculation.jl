calcTd2e(Td::Float32) = 6.1078 * exp((2.5e6/461.5181) * (1/273.16 - 1/Td))
calce2q(e::Float32,p::Float32) = e * 0.6219838793551742 / (p - e * 0.3780161206448258)
calcTm2Pi(Tm::Real) = 10^6 / ((3.739e3 / Tm + 0.221) * 461.5181) / 1000

function calculate(e5ds::ERA5Dataset,isprecise)
    
    dt = e5ds.dtbeg; ndt = daysinmonth(dt) * 24
    p = era5Pressures(); p = p[p.>=50]; np = length(p)
    
    ind = zeros(Bool,np+2)
    bot = zeros(Float32,np+2)
    ita = zeros(Float32,np+2)
    ish = zeros(Float32,np+2)
    ipv = vcat(0,p,0)

    sds = NCDataset(joinpath(e5ds.eroot,"tmpnc-single-$dt.nc"))
    
    if isprecise
        pds = Vector{NCDataset}(undef,np)
        for ip in 1 : np
            pds[ip] = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-$(p[ip])-$dt.nc"))
        end
    else
        tds = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-t-$dt.nc"))
        qds = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-q-$dt.nc"))
    end

    lsd  = getLandSea(e5ds,ERA5Region(GeoRegion("GLB"),gres=0.25))
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    ts = zeros(Float32,nlon,nlat)
    td = zeros(Float32,nlon,nlat)
    sp = zeros(Float32,nlon,nlat)
    ta = zeros(Float32,nlon,nlat,np)
    sh = zeros(Float32,nlon,nlat,np)
    
    tmp = zeros(Int16,nlon,nlat)

    tm = zeros(Float32,nlon,nlat,ndt*24)
    Pi = zeros(Float32,nlon,nlat,ndt*24)
    
    p = Float32.(p*100)

    for it in 1 : ndt

        sc = sds["t2m"].attrib["scale_factor"]
        of = sds["t2m"].attrib["add_offset"]
        mv = sds["t2m"].attrib["missing_value"]
        fv = sds["t2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["t2m"].var,tmp,:,:,it)
        int2real!(ts,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["d2m"].attrib["scale_factor"]
        of = sds["d2m"].attrib["add_offset"]
        mv = sds["d2m"].attrib["missing_value"]
        fv = sds["d2m"].attrib["_FillValue"]
        NCDatasets.load!(sds["d2m"].var,tmp,:,:,it)
        int2real!(td,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        sc = sds["sp"].attrib["scale_factor"]
        of = sds["sp"].attrib["add_offset"]
        mv = sds["sp"].attrib["missing_value"]
        fv = sds["sp"].attrib["_FillValue"]
        NCDatasets.load!(sds["sp"].var,tmp,:,:,it)
        int2real!(sp,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)
        
        if isprecise

            for ip = 1 : np

                taip = @view ta[:,:,ip]
                ship = @view sh[:,:,ip]

                sc = pds[ip]["t"].attrib["scale_factor"]
                of = pds[ip]["t"].attrib["add_offset"]
                mv = pds[ip]["t"].attrib["missing_value"]
                fv = pds[ip]["t"].attrib["_FillValue"]
                NCDatasets.load!(pds[ip]["t"].var,tmp,:,:,it)
                int2real!(taip,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

                sc = pds[ip]["q"].attrib["scale_factor"]
                of = pds[ip]["q"].attrib["add_offset"]
                mv = pds[ip]["q"].attrib["missing_value"]
                fv = pds[ip]["q"].attrib["_FillValue"]
                NCDatasets.load!(pds[ip]["q"].var,tmp,:,:,it)
                int2real!(ship,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            end

        else

            sc = tds["t"].attrib["scale_factor"]
            of = tds["t"].attrib["add_offset"]
            mv = tds["t"].attrib["missing_value"]
            fv = tds["t"].attrib["_FillValue"]
            NCDatasets.load!(tds["t"].var,tmp,:,:,:,it)
            int2real!(ta,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            sc = qds["q"].attrib["scale_factor"]
            of = qds["q"].attrib["add_offset"]
            mv = qds["q"].attrib["missing_value"]
            fv = qds["q"].attrib["_FillValue"]
            NCDatasets.load!(qds["q"].var,tmp,:,:,:,it)
            int2real!(sh,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

        end


        for ilat = 1 : nlat, ilon = 1 : nlon

            its = ts[ilon,ilat]
            itd = td[ilon,ilat]
            isp = sp[ilon,ilat]

            for ip = 2 : (np+1)
                ita[ip] = ta[ilon,ilat,ip-1]
                ish[ip] = sh[ilon,ilat,ip-1]
            end

            ita[end] = its
            ish[end] = calce2q(calcTd2e(itd),isp)
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

            tm[ilon,ilat,it] = integrate(ipp,top) / integrate(ipp,btm)
            Pi[ilon,ilat,it] = calcTm2Pi(tm[ilon,ilat,it])

        end

    end
    
    close(sds)
    close(pds)

    ERA5Reanalysis.save(
        tm,dt,e5ds,
        SingleVariable("t_qwm"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

    ERA5Reanalysis.save(
        Pi,dt,e5ds,
        SingleVariable("Pi"),
        ERA5Region(GeoRegion("GLB"),gres=0.25),
        lsd
    )

end