module TmPi

## Base Modules Used
using Base64
using Logging
using Printf

## Modules Used
using HTTP
using JSON
using NCDatasets
using Trapz
using Statistics

## Import relevant functions from ERA5Reanalysis.jl
using ERA5Reanalysis: ERA5Dataset, ERA5Variable, ERA5Region, ERA5Hourly
using ERA5Reanalysis: SingleLevel, PressureLevel, SingleVariable, PressureVariable
using ERA5Reanalysis: LandSea, getLandSea
using ERA5Reanalysis: isSingle, era5Pressures

## Reexporting exported functions within these modules
using Reexport
@reexport using Dates
@reexport using GeoRegions

import Downloads: download
import ERA5Reanalysis

## Exporting the following functions:
export
        create, analysis, readTm, readPi, TmPiDataset

## Abstract SuperTypes
"""
    TmPiDataset

Abstract supertype for temporary arrays used to calculate the Tm and Pi datasets.
"""
abstract type TmPiDataset <: ERA5Dataset end

## TmPi.jl logging preface

modulelog() = "$(now()) - TmPi.jl"

## Creating Tm and Pi Single-Level variables

function __init__()
    
    if !isSingle("Tm",throw=false)
        SingleVariable(
            varID = "Tm",
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

include("dataset.jl")
include("create.jl")
include("cdsapi.jl")
include("download.jl")
include("calculation.jl")
include("analysis.jl")
include("filesystem.jl")
include("read.jl")
include("tmpi2era5.jl")
include("backend.jl")

end
