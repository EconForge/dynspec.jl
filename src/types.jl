module Types

using LinearAlgebra

export Normal, TimedVar, Complementarity, Equation, Context, ContextBuilder, build_context
export AbstractModelBlock, ModelBlock, ModelSpec, N, equations

"""
    Normal

Gaussian distribution with mean vector `Μ` and covariance matrix `Σ`.
"""
struct Normal
    Σ::Matrix{Float64}
    Μ::Vector{Float64}

    function Normal(Σ::AbstractMatrix{<:Real}, Μ::Union{Nothing, AbstractVector{<:Real}}=nothing)
        S = Matrix{Float64}(Σ)
        d = size(S, 1)
        size(S, 2) == d || error("Covariance matrix must be square")
        issymmetric(S) || isapprox(S, S'; atol=1e-10) || error("Covariance matrix must be symmetric")
        m = Μ === nothing ? zeros(Float64, d) : Vector{Float64}(Μ)
        length(m) == d || error("Mean vector dimension does not match covariance matrix")
        new(S, m)
    end
end

Normal(σ::Real) = Normal(fill(Float64(σ)^2, 1, 1), [0.0])
Normal(μ::Real, σ::Real) = Normal(fill(Float64(σ)^2, 1, 1), [Float64(μ)])

function N(args...)
    if length(args) == 1
        return Normal(args[1])
    elseif length(args) == 2
        return Normal(args[1], args[2])
    else
        throw(ArgumentError("N() expects 1 or 2 arguments, got $(length(args))"))
    end
end

"""
    TimedVar(name, shift)

A variable with a time shift, e.g. `x[t]` (shift 0), `x[t-1]` (shift -1), or `x[~]` (shift :~).
"""
struct TimedVar
    name::Symbol
    shift::Union{Int, Symbol}
end

Base.show(io::IO, v::TimedVar) = print(io, v.name, "[", v.shift == :~ ? "~" : (v.shift > 0 ? "t+$(v.shift)" : (v.shift < 0 ? "t$(v.shift)" : "t")), "]")

"""
    Complementarity(variable, lb, ub)

Complementarity condition `residual ⟂ lb <= variable <= ub`.
"""
struct Complementarity
    variable::TimedVar
    lb::Any
    ub::Any
end

"""
    Equation(residual, raw, metadata, complementarity, tree=nothing)

A model equation with its residual expression, source string, metadata, optional complementarity,
and optional parsed AST node.
"""
struct Equation
    residual::Any
    raw::String
    metadata::Dict{Symbol, Any}
    complementarity::Union{Nothing, Complementarity}
    tree::Any
end

Equation(residual, raw::String, metadata::Dict{Symbol, Any}, comp::Union{Nothing, Complementarity}) =
    Equation(residual, raw, metadata, comp, nothing)

"""
    Context

Model context holding parameters, steady states, processes, time-dated values, and metadata
as dictionaries.
"""
struct Context
    constants::Dict{Symbol, Any}
    steady_states::Dict{Symbol, Any}
    processes::Dict{Symbol, Any}
    values::Dict{Symbol, Dict{Int, Any}}
    metadata::Dict{Symbol, Any}
end

"""
    ContextBuilder

Mutable container used during sequential parsing and evaluation.
"""
mutable struct ContextBuilder
    constants::Dict{Symbol, Any}
    steady_states::Dict{Symbol, Any}
    processes::Dict{Symbol, Any}
    values::Dict{Symbol, Dict{Int, Any}}
    metadata::Dict{Symbol, Any}

    ContextBuilder() = new(
        Dict{Symbol, Any}(),
        Dict{Symbol, Any}(),
        Dict{Symbol, Any}(),
        Dict{Symbol, Dict{Int, Any}}(),
        Dict{Symbol, Any}()
    )
end

"""
    build_context(builder::ContextBuilder) -> Context

Convert a mutable `ContextBuilder` into an immutable `Context`.
"""
function build_context(builder::ContextBuilder)
    return Context(
        builder.constants,
        builder.steady_states,
        builder.processes,
        builder.values,
        builder.metadata
    )
end

"""
    AbstractModelBlock

Abstract interface for a model block. Any `AbstractModelBlock` has an associated
`context` and can evaluate residuals and derivatives.
"""
abstract type AbstractModelBlock end

"""
    ModelBlock(context, equations; metadata=Dict{Symbol, Any}(), filename="", source="", tree=nothing)

A block of equations associated with a model evaluation context and preserving the parsed AST `tree`.
"""
struct ModelBlock{CTX<:Context} <: AbstractModelBlock
    context::CTX
    equations::Vector{Equation}
    metadata::Dict{Symbol, Any}
    filename::String
    source::String
    tree::Any
end

ModelBlock(context::CTX, equations::Vector{Equation};
           metadata=Dict{Symbol, Any}(), filename="", source="", tree=nothing) where {CTX<:Context} =
    ModelBlock(context, equations, Dict{Symbol, Any}(metadata), String(filename), String(source), tree)

const ModelSpec = ModelBlock

"""
    equations(b::AbstractModelBlock) -> Vector{Equation}

Return the list of equations in `b`.
"""
equations(b::ModelBlock) = b.equations

# Collection interface on ModelBlock
Base.length(b::ModelBlock) = length(b.equations)
Base.iterate(b::ModelBlock, state...) = iterate(b.equations, state...)
Base.getindex(b::ModelBlock, i::Int) = b.equations[i]
Base.getindex(b::ModelBlock, inds) = b.equations[inds]
Base.eltype(::Type{<:ModelBlock}) = Equation

Base.show(io::IO, b::ModelBlock) = print(io, "ModelBlock($(length(b)) equations, context=$(length(keys(b.context.constants))) constants, $(length(keys(b.context.steady_states))) steady states)")

end # module Types
