module Dynspec

include("types.jl")
include("parser.jl")
include("evaluator.jl")
include("analysis.jl")
include("dyno.jl")

using .Types
using .Parser
using .Evaluator
using .Analysis
using .Dyno
using YAML

# Re-export key types
export Normal, N, TimedVar, Complementarity, Equation, Context, AbstractModelBlock, ModelBlock, ModelSpec

# Re-export key functions
export parse_model, read_model, list_variables, list_constants, extract_variables, eval_residuals, residuals, filter_equations, equations
export jacobian, derivatives

# Submodule for DSGE / Dyno-specific helpers
export Dyno

"""
    parse_model(text::AbstractString; filename::AbstractString="<string>") -> ModelBlock

Parse a DynSpec model from a string, interpreting assignments and equations.
"""
function parse_model(text::AbstractString; filename::AbstractString="<string>")
    tree = Parser.parse_tree(text; start="free_block")
    return Evaluator.interpret_model(tree, String(filename), String(text))
end

"""
    read_model(filename::AbstractString) -> ModelBlock

Read and interpret a DynSpec model file (`.dyno`).
"""
function read_model(filename::AbstractString)
    content = read(filename, String)
    return parse_model(content; filename=filename)
end

end # module Dynspec
