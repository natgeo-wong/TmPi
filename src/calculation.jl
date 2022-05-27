calcTd2e(Td::Float32) = 6.1078 * exp((2.5e6/461.5181) * (1/273.16 - 1/Td))
calce2q(e::Float32,p::Float32) = e * 0.6219838793551742 / (p - e * 0.3780161206448258)
calcTm2Pi(Tm::Real) = 10^6 / ((3.739e3 / Tm + 0.221) * 461.5181) / 1000

function calculate(e5ds::ERA5Dataset,isprecise)
    
    dt = e5ds.dtbeg; ndt = daysinmonth(dt) * 24
    p = era5Pressures(); p = p .> 50; np = length(p)
    
    ind = Vector{Bool}(0,np+2)
    bot = Vector{Float32}(0,np+2)
    ita = Vector{Float32}(0,np+2)
    ish = Vector{Float32}(0,np+2)
    ipv = vcat(0,p,0)

    sds = NCDataset(joinpath(e5ds.eroot,"tmpnc-single.nc"))
    
    if isprecise
        pds = Vector{NCDataset}(undef,np)
        for ip in 1 : np
            pds[ip] = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure-$(p[ip]).nc"))
        end
    else
        pds = NCDataset(joinpath(e5ds.eroot,"tmpnc-pressure.nc"))
    end

    lsd  = getLandSea(e5ds,ERA5Region(GeoRegion("GLB"),gres=0.25))
    nlon = length(lsd.lon)
    nlat = length(lsd.lat)

    ts = Array{Float32,2}(undef,nlon,nlat)
    td = Array{Float32,2}(undef,nlon,nlat)
    sp = Array{Float32,2}(undef,nlon,nlat)
    ta = Array{Float32,3}(0    ,nlon,nlat,np)
    sh = Array{Float32,3}(0    ,nlon,nlat,np)
    
    tmp = Array{Int16,2}(undef,nlon,nlat)

    tm = Array{Float32,3}(undef,nlon,nlat,ndt*24)
    Pi = Array{Float32,3}(undef,nlon,nlat,ndt*24)
    
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

            sc = pds["t"].attrib["scale_factor"]
            of = pds["t"].attrib["add_offset"]
            mv = pds["t"].attrib["missing_value"]
            fv = pds["t"].attrib["_FillValue"]
            NCDatasets.load!(pds["t"].var,tmp,:,:,:,it)
            int2real!(ta,tmp,scale=sc,offset=of,mvalue=mv,fvalue=fv)

            sc = pds["q"].attrib["scale_factor"]
            of = pds["q"].attrib["add_offset"]
            mv = pds["q"].attrib["missing_value"]
            fv = pds["q"].attrib["_FillValue"]
            NCDatasets.load!(pds["q"].var,tmp,:,:,:,it)
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