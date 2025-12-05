# Import all required packages for the double pendulum simulation
# Include this file in your main script with: include("imports.jl")

using DifferentialEquations  # For solving differential equations
using Plots                  # For visualization
using CSV                    # For reading .trk or .csv files
using DataFrames            # For data manipulation
using Optim                 # For optimization (finding masses and damping)
using LaTeXStrings          # For LaTeX formatting in plots
using LinearAlgebra         # For matrix operations (standard library)
using Statistics            # For statistical analysis (standard library)
using Markdown              # For markdown documentation blocks (standard library)

println("All packages loaded successfully!")
