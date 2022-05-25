module CreateTmPiDataset

## Modules Used
using Logging
using NCDatasets
using Printf
using Statistics

## Reexporting exported functions within these modules
using Reexport
@reexport using Dates
@reexport using GeoRegions
@reexport using ERA5Reanalysis: ERA5Dataset, ERA5Variable, ERA5Region, ERA5Hourly
@reexport using ERA5Reanalysis: SingleVariable, PressureVariable, isSingle
@reexport using ERA5Reanalysis: era5Pressures

import Base: show, read, download

## Exporting the following functions:
export
        create

## CreateTmPiDataset.jl logging preface

modulelog() = "$(now()) - CreateTmPiDataset.jl"

## Creating Tm and Pi Single-Level variables

function __init__()
    
    if !isSingle("t_qwm",throw=false)
        SingleVariable(
            varID = "t_qwm",
            lname = "water_vapour_weighted_mean_temperature",
            vname = "Water Vapour Weighted Mean Temperature",
            units = "K",
            inCDS = false
        )
    end

    if !isSingle("Pi",throw=false)
        SingleVariable(
            varID = "Pi",
            lname = "pi_conversion_constant",
            vname = "Pi Conversion Constant [Askne and Nordius 1987]",
            units = "N/A",
            inCDS = false
        )
    end

end

## Including other files in the module

include("backend.jl")
include("calculation.jl")
include("cdsapi.jl")
include("create.jl")

end
