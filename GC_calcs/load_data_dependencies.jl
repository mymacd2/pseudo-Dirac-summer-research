include("read_nu_dist_sphere.jl")
include("read_eff_area.jl")
include("read_energy_res.jl")
include("functions_utils.jl")
include("MC_smear.jl")

# Packages

using DelimitedFiles
using StatsBase
using Interpolations
using Plots
using LaTeXStrings
using FFTW
using Distributions
using Profile
using PyCall
using SpecialFunctions
using Optim
using StaticArrays
using Healpix
using BenchmarkTools
using CurveFit
using ApproxFun
using LinearAlgebra

# Importing the necessary python libraries
@pyimport matplotlib.pyplot as plt
@pyimport numpy as np
@pyimport healpy as hp



# Loading the models and necessary data

model = "galnu_GALcase2"

# Getting spatial data
rs0, ls0, bs0, rws0, us0, logCRes0 = read_models(model)
oneweights0 = Weights(ones(length(rws0)))

# Getting the (previously unneeded) energy gc_data
minu, maxu = minimum(us0), maximum(us0)
logνebins = range(minu, maxu, 100)
;