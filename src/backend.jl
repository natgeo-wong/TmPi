yr2str(date::TimeType) = Dates.format(date,dateformat"yyyy")
yrmo2str(date::TimeType) = Dates.format(date,dateformat"yyyymm")

function int2real!(
    oarray :: AbstractArray{FT},
    iarray :: AbstractArray{Int16};
    scale  :: Real,
    offset :: Real,
    fvalue :: Int16,
    mvalue :: Int16
) where FT <: Real

    for ii = 1 : length(iarray)

        if (iarray[ii] == fvalue) || (iarray[ii] == mvalue)
              oarray[ii] = FT(NaN)
        else; oarray[ii] = iarray[ii] * scale + offset
        end

    end

    return

end

function ncoffsetscale(data::AbstractArray{<:Real})

    dmax = data[findfirst(!isnan,data)]
    dmin = data[findfirst(!isnan,data)]
    for ii = 1 : length(data)
        dataii = data[ii]
        if !isnan(dataii)
            if dataii > dmax; dmax = dataii end
            if dataii < dmin; dmin = dataii end
        end
    end

    scale = (dmax-dmin) / 65531;
    offset = (dmax+dmin-scale) / 2;

    return scale,offset

end

function nanmean(
    data :: AbstractArray,
    dNaN :: AbstractArray
)
    nNaN = length(dNaN)
    for iNaN in 1 : nNaN
        dNaN[iNaN] = !isnan(data[iNaN])
    end
    dataii = @view data[dNaN]
    if isempty(dataii); return mean(dataii); else; return NaN; end
end