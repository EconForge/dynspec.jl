"""
    Dyno

Submodule containing DSGE perturbation matrix extractors and model helpers.
"""
module Dyno

using ..Types
using ..Evaluator
using ..Parser
using ..Analysis
using LinearAlgebra

export compute_dsge_jacobians, compute_jacobians, compute_jacobians!
export DynamicModel, compile_dynamic, first_order_model
export endogenous_exogenous, extract_endogenous_exogenous, list_endogenous_exogenous
export solve_dsge, simulate_linear

"""
    endogenous_exogenous(block::AbstractModelBlock; endogenous=nothing, exogenous=nothing)
        -> (endogenous = Vector{Symbol}, exogenous = Vector{Symbol})

Extract the ordered lists of endogenous and exogenous variables from `block`.
If not explicitly provided, exogenous variables are inferred from declared processes
and endogenous variables from steady states and equation variable occurrences.
"""
function endogenous_exogenous(block::AbstractModelBlock; endogenous=nothing, exogenous=nothing)
    exo_syms = if exogenous !== nothing
        Symbol.(exogenous)
    else
        collect(keys(block.context.processes))
    end
    exo_set = Set{Symbol}(exo_syms)

    all_vars = list_variables(block; shifts=true)
    endo_syms = if endogenous !== nothing
        Symbol.(endogenous)
    else
        ordered = Symbol[]
        for s in keys(block.context.steady_states)
            if s ∉ exo_set && haskey(all_vars, s) && s ∉ ordered
                push!(ordered, s)
            end
        end
        for (v, _) in all_vars
            if v ∉ exo_set && v ∉ ordered
                push!(ordered, v)
            end
        end
        ordered
    end
    return (endogenous = endo_syms, exogenous = exo_syms)
end

const extract_endogenous_exogenous = endogenous_exogenous
const list_endogenous_exogenous = endogenous_exogenous

"""
    first_order_model(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing)

Construct a first-order DSGE model representation as a NamedTuple:
- `symbols`: NamedTuple of variable and parameter groups:
  - `endogenous`: ordered list of endogenous variable symbols
  - `exogenous`: ordered list of exogenous shock symbols
  - `constants`: ordered list of constant / parameter symbols
- `residuals`: residual vector evaluated at steady state (or `at`)
- `jacobians`: NamedTuple of first-order Jacobians:
  - `A`: ∂f / ∂y_{t+1} (lead)
  - `B`: ∂f / ∂y_t     (current)
  - `C`: ∂f / ∂y_{t-1} (lag)
  - `D`: ∂f / ∂e_t     (shocks)
"""
function first_order_model(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing)
    endo_syms, exo_syms = endogenous_exogenous(block; endogenous=endogenous, exogenous=exogenous)
    const_syms = sort(collect(keys(block.context.constants)))

    res, A, B, C, D = compute_dsge_jacobians(block; filter_by=filter_by, endogenous=endo_syms, exogenous=exo_syms, at=at)

    return (
        symbols = (
            endogenous = endo_syms,
            exogenous  = exo_syms,
            constants  = const_syms,
        ),
        residuals = res,
        jacobians = (
            A = A,
            B = B,
            C = C,
            D = D,
        )
    )
end

"""
    compute_dsge_jacobians(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing)
        -> (residuals, A, B, C, D)

[TEMPORARY / CONVENIENCE HELPER]
Evaluate residuals and compute standard first-order perturbation matrices for DSGE solvers:
- `A` = ∂f / ∂y_{t+1}
- `B` = ∂f / ∂y_t
- `C` = ∂f / ∂y_{t-1}
- `D` = ∂f / ∂e_t
"""
function compute_dsge_jacobians(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing)::Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    endo_syms, exo_syms = endogenous_exogenous(block; endogenous=endogenous, exogenous=exogenous)

    wrt_grouped = [
        [(v, 1) for v in endo_syms],
        [(v, 0) for v in endo_syms],
        [(v, -1) for v in endo_syms],
        [(v, 0) for v in exo_syms]
    ]

    res, (A, B, C, D) = residuals_and_jacobian(block; wrt=wrt_grouped, filter_by=filter_by, at=at)
    return (res, A, B, C, D)
end

# Backwards compatibility alias
const compute_jacobians = compute_dsge_jacobians

# DynamicModel wrapper
struct DynamicModel{M<:ModelBlock}
    model::M
    endogenous::Vector{Symbol}
    exogenous::Vector{Symbol}
end

"""
    compile_dynamic(block::ModelBlock; endogenous=nothing, exogenous=nothing) -> DynamicModel

Compile a model block into a `DynamicModel` with resolved endogenous and exogenous variable orderings.
"""
function compile_dynamic(block::ModelBlock; endogenous=nothing, exogenous=nothing)
    endo_syms, exo_syms = endogenous_exogenous(block; endogenous=endogenous, exogenous=exogenous)
    return DynamicModel(block, endo_syms, exo_syms)
end

Analysis.residuals(dm::DynamicModel) = residuals(dm.model)
Analysis.residuals(dm::DynamicModel, y_next, y_curr, y_prev, e) = compute_dsge_jacobians(dm, y_next, y_curr, y_prev, e)[1]

function compute_dsge_jacobians(dm::DynamicModel)
    return compute_dsge_jacobians(dm.model; endogenous=dm.endogenous, exogenous=dm.exogenous)
end

function compute_dsge_jacobians(dm::DynamicModel, y_next, y_curr, y_prev, e)
    at = Dict{Tuple{Symbol, Int}, Float64}()
    for (i, v) in enumerate(dm.endogenous)
        at[(v, 1)] = Float64(y_next[i])
        at[(v, 0)] = Float64(y_curr[i])
        at[(v, -1)] = Float64(y_prev[i])
    end
    for (i, v) in enumerate(dm.exogenous)
        at[(v, 0)] = Float64(e[i])
    end
    return compute_dsge_jacobians(dm.model; endogenous=dm.endogenous, exogenous=dm.exogenous, at=at)
end

"""
    solve_dsge(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix; max_iter::Int=50, tol::Float64=1e-11) -> (P, Q)

Solve the first-order linear DSGE system:
    A * E_t[y_{t+1}] + B * y_t + C * y_{t-1} + D * e_t = 0
using cyclic reduction to find the stable recursive law of motion matrices P and Q:
    y_t = P * y_{t-1} + Q * e_t
"""
function solve_dsge(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, D::AbstractMatrix; max_iter::Int=50, tol::Float64=1e-11)
    A0, B0, C0, B_hat = copy(A), copy(B), copy(C), copy(B)
    for _ in 1:max_iter
        invB = B0 \ I
        A1 = A0 * invB * A0
        C1 = C0 * invB * C0
        B1 = B0 - A0 * invB * C0 - C0 * invB * A0
        B_hat = B_hat - A0 * invB * C0
        A0, C0, B0 = A1, C1, B1
        (norm(A0, Inf) < tol || norm(C0, Inf) < tol) && break
    end
    P = - (B_hat \ C)
    Q = - (A*P + B) \ D
    return (P=P, Q=Q)
end

"""
    solve_dsge(model::NamedTuple; max_iter::Int=50, tol::Float64=1e-11)

Solve a first-order DSGE model (as returned by `first_order_model(block)`) for the recursive decision rules:
    y_t = P * y_{t-1} + Q * e_t
"""
function solve_dsge(m::NamedTuple; max_iter::Int=50, tol::Float64=1e-11)
    if haskey(m, :jacobians) && haskey(m, :symbols)
        j = m.jacobians
        sol = solve_dsge(j.A, j.B, j.C, j.D; max_iter=max_iter, tol=tol)
        return (P=sol.P, Q=sol.Q, symbols=m.symbols, residuals=m.residuals, jacobians=m.jacobians,
                endogenous=m.symbols.endogenous, exogenous=m.symbols.exogenous,
                A=j.A, B=j.B, C=j.C, D=j.D)
    else
        error("Expected a first_order_model NamedTuple with :jacobians and :symbols")
    end
end

"""
    solve_dsge(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing, max_iter::Int=50, tol::Float64=1e-11)

Compute first-order perturbation matrices and solve for the recursive decision rules `P` and `Q` (`y_t = P * y_{t-1} + Q * e_t`).
"""
function solve_dsge(block::AbstractModelBlock; filter_by=nothing, endogenous=nothing, exogenous=nothing, at=nothing, max_iter::Int=50, tol::Float64=1e-11)
    m = first_order_model(block; filter_by=filter_by, endogenous=endogenous, exogenous=exogenous, at=at)
    return solve_dsge(m; max_iter=max_iter, tol=tol)
end

"""
    simulate_linear(P, Q; T=40, shocks=nothing, shock_idx=1, shock_size=1.0) -> Matrix{Float64}

Simulate the linear system `y_t = P * y_{t-1} + Q * e_t` for `T` periods.
If `shocks` is omitted, computes the Impulse Response Function (IRF) to a one-time shock of `shock_size` on shock `shock_idx` at `t=1`.
If `shocks` is an `n_shocks × T` matrix, simulates under that shock sequence.
"""
function simulate_linear(P::AbstractMatrix, Q::AbstractMatrix; T::Int=40, shocks::Union{Nothing, AbstractMatrix}=nothing, shock_idx::Int=1, shock_size::Float64=1.0)
    n_vars = size(P, 1)
    n_shocks = size(Q, 2)
    Y = zeros(n_vars, T)
    if shocks === nothing
        e1 = zeros(n_shocks)
        e1[shock_idx] = shock_size
        Y[:, 1] = Q * e1
        for t in 2:T
            Y[:, t] = P * Y[:, t-1]
        end
    else
        for t in 1:T
            e_t = shocks[:, t]
            y_prev = (t == 1) ? zeros(n_vars) : Y[:, t-1]
            Y[:, t] = P * y_prev + Q * e_t
        end
    end
    return Y
end

const compute_jacobians! = compute_dsge_jacobians

end # module Dyno
