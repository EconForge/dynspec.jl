#!/usr/bin/env julia
using Pkg
Pkg.activate(@__DIR__)

using Pluto

port = parse(Int, get(ENV, "PLUTO_PORT", "1234"))
notebook_file = joinpath(@__DIR__, "interactive_model.jl")

println("="^60)
println(" Dynspec.jl — Pluto Notebook Environment")
println(" Fast Plotting Backend: PlotlyLight (< 0.2s load time)")
println(" Directory: ", @__DIR__)
println(" Default Notebook: ", basename(notebook_file))
println(" Launching on port: ", port)
println("="^60)

cd(@__DIR__)
Pluto.run(; port=port, notebook=notebook_file, notebook_path_suggestion=@__DIR__)
