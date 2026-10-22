using LiveServer

println("Starting LiveServer doc watcher...")
println("Open your browser at: http://localhost:8000")
println("Any edits in docs/src/ or src/ will trigger a rebuild and live browser reload.")

servedocs(
    doc_env = true,
    makejl = joinpath(@__DIR__, "make.jl")
)
