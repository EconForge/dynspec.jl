module Analysis

using ..Types
using ..Evaluator
using ..Parser
using ForwardDiff
using LinearAlgebra

export list_variables, list_constants, extract_variables, eval_residuals, residuals, filter_equations
export jacobian, derivatives

"""
    list_variables(block::AbstractModelBlock; shifts::Bool=false) -> Vector{Symbol} or Dict{Symbol, Set{Int}}
    list_variables(eq::Equation; shifts::Bool=false)

List all variables occurring in equations.
- `shifts=false` (default): returns a sorted `Vector{Symbol}` of variable names.
- `shifts=true`: returns `Dict{Symbol, Set{Int}}` mapping variables to time shifts.
"""
function list_variables(eq::Equation; shifts::Bool=false)
    dict = Dict{Symbol, Set{Int}}()
    Parser.collect_variables!(dict, eq.residual)
    return shifts ? dict : sort(collect(keys(dict)))
end

function list_variables(block::AbstractModelBlock; shifts::Bool=false)
    dict = Dict{Symbol, Set{Int}}()
    for eq in block.equations
        Parser.collect_variables!(dict, eq.residual)
    end
    return shifts ? dict : sort(collect(keys(dict)))
end

"""
    list_constants(block::AbstractModelBlock) -> Vector{Symbol}

List all constant/parameter names defined in the block's context, sorted alphabetically.
"""
list_constants(block::AbstractModelBlock) = sort(collect(keys(block.context.constants)))

function extract_variables(target)
    if target isa AbstractModelBlock || target isa Equation
        return list_variables(target; shifts=true)
    else
        dict = Dict{Symbol, Set{Int}}()
        Parser.collect_variables!(dict, target)
        return dict
    end
end

function _filter_equations(equations::Vector{Equation}, filter_by)
    if filter_by === nothing
        return equations
    elseif filter_by isa Symbol
        return filter(eq -> filter_by in get(eq.metadata, :tags, Symbol[]), equations)
    elseif filter_by isa AbstractString
        sym = Symbol(filter_by)
        return filter(eq -> sym in get(eq.metadata, :tags, Symbol[]), equations)
    elseif filter_by isa AbstractVector
        syms = Set(Symbol.(filter_by))
        return filter(eq -> any(t in syms for t in get(eq.metadata, :tags, Symbol[])), equations)
    elseif filter_by isa Function
        return filter(filter_by, equations)
    else
        error("Unsupported filter_by specification: $filter_by")
    end
end

"""
    residuals(block::AbstractModelBlock; filter_by=nothing, at=nothing) -> Vector{Float64}

Compute equation residuals directly on the parsed AST equations.
- `filter_by`: filter equations by tag (`:transition`), tags (`[:transition, :arbitrage]`), or predicate.
- `at`: overrides for variables (e.g. `(:z, 0) => 0.1` or `:z => 0.1`) or constants (`:α => 0.36`).
"""
function residuals(block::AbstractModelBlock; filter_by=nothing, at=nothing)::Vector{Float64}
    eqs = _filter_equations(block.equations, filter_by)
    env = Dict{Tuple{Symbol, Int}, Float64}()
    constants_override = Dict{Symbol, Float64}()

    if at isa AbstractDict
        for (k, v) in at
            if k isa Tuple{Symbol, Int}
                env[k] = Float64(v)
            elseif k isa Symbol
                if haskey(block.context.constants, k)
                    constants_override[k] = Float64(v)
                else
                    env[(k, 0)] = Float64(v)
                end
            end
        end
    elseif at isa NamedTuple
        for (k, v) in pairs(at)
            if haskey(block.context.constants, k)
                constants_override[k] = Float64(v)
            else
                env[(k, 0)] = Float64(v)
            end
        end
    end

    res = Vector{Float64}(undef, length(eqs))
    for (i, eq) in enumerate(eqs)
        res[i] = Float64(Evaluator.eval_formula(eq.residual, block.context; env=env, constants=constants_override))
    end
    return res
end

eval_residuals(block::AbstractModelBlock; filter_by=nothing, at=nothing) = residuals(block; filter_by=filter_by, at=at)

"""
    filter_equations(block::AbstractModelBlock, tag::Union{Symbol, String}) -> Vector{Equation}

Filter equations belonging to a specific tag (e.g. `:transition`, `:arbitrage`).
"""
function filter_equations(block::AbstractModelBlock, tag::Union{Symbol, String})
    t_sym = Symbol(tag)
    return filter(eq -> t_sym in get(eq.metadata, :tags, Symbol[]), block.equations)
end

# Target normalization:
# Bare symbol -> constant/parameter
# Tuple{Symbol, Int} or TimedVar -> dynamic variable
_normalize_target(s::Symbol) = s
_normalize_target(t::Tuple{Symbol, Int}) = t
_normalize_target(v::TimedVar) = (v.name, v.shift isa Int ? v.shift : 0)
_normalize_target(p::Pair{Symbol, Int}) = (p.first, p.second)
_normalize_target(x) = error("Unsupported target in wrt: $x. Expected Symbol (constant), Tuple{Symbol, Int} (variable), or TimedVar.")

function _derivatives_flat(block::AbstractModelBlock, wrt_norm, filter_by, at)
    eqs = _filter_equations(block.equations, filter_by)
    n_targets = length(wrt_norm)
    n_eqs = length(eqs)

    x0 = zeros(Float64, n_targets)
    if at isa AbstractVector
        length(at) == n_targets || error("Length of `at` ($(length(at))) does not match length of `wrt` ($n_targets)")
        x0 .= Float64.(at)
    else
        for (i, target) in enumerate(wrt_norm)
            val = if target isa Symbol
                # Constant
                if at isa AbstractDict && haskey(at, target)
                    Float64(at[target])
                elseif at isa NamedTuple && haskey(at, target)
                    Float64(getfield(at, target))
                elseif haskey(block.context.constants, target)
                    Float64(block.context.constants[target])
                else
                    0.0
                end
            else
                # Variable (varname, shift)
                var, shift = target
                if at isa AbstractDict && haskey(at, (var, shift))
                    Float64(at[(var, shift)])
                elseif at isa AbstractDict && shift == 0 && haskey(at, var)
                    Float64(at[var])
                elseif at isa NamedTuple && shift == 0 && haskey(at, var)
                    Float64(getfield(at, var))
                elseif haskey(block.context.steady_states, var)
                    Float64(block.context.steady_states[var])
                else
                    0.0
                end
            end
            x0[i] = val
        end
    end

    base_env = Dict{Tuple{Symbol, Int}, Float64}()
    base_constants = Dict{Symbol, Float64}()
    if at isa AbstractDict
        for (k, v) in at
            if k isa Tuple{Symbol, Int}
                base_env[k] = Float64(v)
            elseif k isa Symbol
                if haskey(block.context.constants, k)
                    base_constants[k] = Float64(v)
                else
                    base_env[(k, 0)] = Float64(v)
                end
            end
        end
    elseif at isa NamedTuple
        for (k, v) in pairs(at)
            if haskey(block.context.constants, k)
                base_constants[k] = Float64(v)
            else
                base_env[(k, 0)] = Float64(v)
            end
        end
    end

    function f_residual(x::AbstractVector{T}) where {T}
        env = Dict{Tuple{Symbol, Int}, T}()
        constants_override = Dict{Symbol, T}()

        for (k, v) in base_env
            env[k] = T(v)
        end
        for (k, v) in base_constants
            constants_override[k] = T(v)
        end

        for (i, target) in enumerate(wrt_norm)
            if target isa Symbol
                constants_override[target] = x[i]
            else
                env[target] = x[i]
            end
        end

        out = Vector{T}(undef, n_eqs)
        for (i, eq) in enumerate(eqs)
            out[i] = Evaluator.eval_formula(eq.residual, block.context; env=env, constants=constants_override)
        end
        return out
    end

    res = Vector{Float64}(f_residual(x0))
    J = Matrix{Float64}(ForwardDiff.jacobian(f_residual, x0))
    return (res, J)
end

_is_target(x::Symbol) = true
_is_target(x::Tuple{Symbol, Int}) = true
_is_target(x::TimedVar) = true
_is_target(x::Pair{Symbol, Int}) = true
_is_target(x) = false

"""
    derivatives(block::AbstractModelBlock; wrt, filter_by=nothing, at=nothing)

Compute equation residuals and analytical Jacobians with respect to `wrt`.
- Bare symbols (`:α`, `:β`) represent constants/parameters.
- Tuples (`(:z, 0)`, `(:k, -1)`) or `TimedVar`s represent variables with time shifts.
- If `wrt` is grouped (e.g. `[ [(:z,0), (:x,0)], [(:z,1), (:x,1)] ]`), returns `(residuals, (J1, J2, ...))`.
"""
function derivatives(block::AbstractModelBlock; wrt, filter_by=nothing, at=nothing)
    is_grouped = (wrt isa AbstractVector || wrt isa Tuple) &&
                 !isempty(wrt) &&
                 !_is_target(first(wrt)) &&
                 (first(wrt) isa AbstractVector || first(wrt) isa Tuple)
    if is_grouped
        groups = [[_normalize_target(t) for t in g] for g in wrt]
        flat_wrt = reduce(vcat, groups)
        res, J_full = _derivatives_flat(block, flat_wrt, filter_by, at)
        col_idx = 1
        jac_groups = Matrix{Float64}[]
        for g in groups
            n_cols = length(g)
            push!(jac_groups, J_full[:, col_idx:(col_idx + n_cols - 1)])
            col_idx += n_cols
        end
        return (res, Tuple(jac_groups))
    else
        wrt_norm = [_normalize_target(t) for t in wrt]
        return _derivatives_flat(block, wrt_norm, filter_by, at)
    end
end

"""
    jacobian(block::AbstractModelBlock; wrt, filter_by=nothing, at=nothing)

Compute analytical Jacobian(s) with respect to `wrt`.
If `wrt` is a flat vector of targets, returns a single `Matrix{Float64}`.
If `wrt` is grouped (vector of vectors), returns a `Tuple` of `Matrix{Float64}`.
"""
function jacobian(block::AbstractModelBlock; wrt, filter_by=nothing, at=nothing)
    return derivatives(block; wrt=wrt, filter_by=filter_by, at=at)[2]
end

end # module Analysis
