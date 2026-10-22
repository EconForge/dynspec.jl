# Function & Option Guide

This guide provides a comprehensive walkthrough of all main functions, types, and options provided by **Dynspec.jl**.

---

## 1. Model Ingestion

Dynspec models can be parsed from inline strings or read directly from `.dyno` specification files.

### `parse_model`
Parse a model directly from a string:
```julia
using Dynspec

model_str = """
@name: Neoclassical

# Parameters
α <- 0.36
β <- 0.96
δ <- 0.08
ρ <- 0.90

# Calibrated steady states
z[~] <- 0.0
k[~] <- ((1/β - (1-δ))/α)^(1/(α-1))
y[~] <- k[~]^α
i[~] <- δ*k[~]
c[~] <- y[~] - i[~]

# Stochastic processes
e_z[t] <- N(0, 0.01^2)

# Tagged equation blocks
[transition] :: {
    z[t] = ρ*z[t-1] + e_z[t]           :: "Productivity Process"
    k[t] = (1-δ)*k[t-1] + i[t]         :: "Capital Accumulation"
}

[definitions] :: {
    y[t] = exp(z[t])*k[t-1]^α          :: "Production"
    c[t] = y[t] - i[t]                 :: "Resource Constraint"
}

[arbitrage] :: {
    1 - β*(c[t+1]/c[t])^(-1)*(1-δ + α*y[t+1]/k[t]) ⟂ 0 <= c[t] <= inf :: "Euler Equation"
}
"""

block = parse_model(model_str; filename="<string>")
```

### `read_model`
Read and interpret a `.dyno` model file from disk:
```julia
block = read_model("test/models/neo.dyno")
```

Both functions return an instance of `ModelBlock <: AbstractModelBlock`.

---

## 2. Inspecting the Context & Parameters

The `context` holds all parameters, calibrated steady states, processes, time-dated values, and metadata as standard, typed dictionaries.

### Structured Access via Dictionaries
```julia
block.context.constants       # Dict{Symbol, Any} of parameters: Dict(:α => 0.36, :β => 0.96, ...)
block.context.steady_states   # Dict{Symbol, Any} of steady states: Dict(:k => 5.447, :y => 1.841, ...)
block.context.processes       # Dict{Symbol, Any} of shock distributions: Dict(:e_z => Normal(...))
block.context.values          # Dict{Symbol, Dict{Int, Any}} of time-dated values: values[:k][0]
block.context.metadata        # Dict{Symbol, Any} of metadata: Dict(:name => "Neoclassical", ...)
```

Direct dictionary indexing provides unambiguous, explicit access:
```julia
block.context.constants[:α]       # 0.36
block.context.steady_states[:k]   # 5.4468
```

> [!NOTE]
> Values are stored in standard dictionaries without shorthand property forwarding (e.g. `block.context.α`), avoiding namespace collisions between context fields and model identifiers.

### `model_properties`
To inspect all key model properties in a single call:
```julia
props = model_properties(block)
# Returns a Dict with :name, :parameters, :steady_states, :processes,
# :variables, :dynamic_variables, :equations_count, and :tags
```

---

## 3. Preserved AST Parse Trees

Dynspec preserves the complete syntax parse tree on `ModelBlock` and individual AST statement nodes on each equation:

```julia
# The full Lark parse tree (Lerche.Tree)
block.tree

# Individual equation AST node
block[1].tree
```

This enables re-evaluating model variants (`@run: variants: {β: [0.95, 0.96]}`) or targeting downstream compilers (Symbolics.jl, CasADi, C) without re-parsing raw strings.

---

## 4. Variable Listing (`list_variables`)

### `list_variables(block; [shifts=false])`
Inspect variable names across the entire block or within a single equation:

```julia
# 1. Unique variable names:
list_variables(block)
# Vector{Symbol}: [:c, :e_z, :i, :k, :y, :z]

# 2. Variables with their dynamic time shifts:
list_variables(block; shifts=true)
# Dict{Symbol, Set{Int}} with:
#   :c   => Set([0, 1])
#   :k   => Set([-1, 0])
#   :y   => Set([0, 1])
#   :z   => Set([-1, 0])
#   :e_z => Set([0])

# On a single equation:
list_variables(block[1])
```

*(Note: `extract_variables(block)` is provided as a backward-compatible alias for `list_variables(block; shifts=true)`).*

### `endogenous_exogenous(block; [endogenous], [exogenous])`
Extract the ordered lists of endogenous and exogenous variables from declared processes and equation occurrences:

```julia
# Named tuple access:
vars = endogenous_exogenous(block)
vars.endogenous   # [:z, :k, :y, :i, :c]
vars.exogenous    # [:e_z]

# Or with tuple destructuring:
endo, exo = endogenous_exogenous(block)
```

*(Note: `extract_endogenous_exogenous` and `list_endogenous_exogenous` are provided as aliases).*

---

## 5. Inspecting Equations & Complementarities

`ModelBlock` implements the standard Julia collection interface (`length(block)`, `block[i]`, `for eq in block`), so you do not need to access internal fields:

```julia
# Access equation by index:
eq = block[5]
println(eq.raw)               # "1 - β*(c[t+1]/c[t])^(-1)..."
println(eq.complementarity)   # Complementarity(c[t], 0.0, Inf)

# Or retrieve all equations:
eqs = equations(block)
```

Each `Equation` object contains:
- `residual`: The parsed AST tree node for the residual ($rhs - lhs$).
- `raw::String`: Original textual formula.
- `metadata::Dict{Symbol, Any}`: Tags and annotations.
- `complementarity::Union{Nothing, Complementarity}`: Normalized complementarity conditions (`lb <= x[t] <= ub`).
- `tree::Any`: Preserved syntax AST node.

---

## 6. Residual Evaluation (`residuals`)

### Direct AST Evaluation
Dynspec evaluates equation residuals directly on the parsed Lark syntax tree without intermediate conversion to Julia expressions.

```julia
# 1. Steady-state residuals for all equations:
residuals(block)

# 2. Filter by a single block tag:
residuals(block; filter_by=:transition)

# 3. Filter by multiple block tags:
residuals(block; filter_by=[:transition, :arbitrage])

# 4. Filter by a custom predicate function:
residuals(block; filter_by=eq -> :transition in get(eq.metadata, :tags, Symbol[]))

# 5. Evaluate with custom variable or parameter overrides:
residuals(block; filter_by=:transition, at=Dict(:z => 0.1))
residuals(block; at=(z=0.1, k=5.0))
```

---

## 7. Analytical Jacobians via ForwardDiff (`jacobian` & `derivatives`)

Dynspec computes analytical Jacobians using `ForwardDiff.Dual` numbers evaluated directly on the equation AST trees.

### Semantics of Targets in `wrt`
- **Bare Symbols** (`:α`, `:β`): represent **constants/parameters** (e.g. for parameter sensitivity $\partial f / \partial \alpha$).
- **Time-dated tuples** (`(:z, 0)`, `(:k, -1)`, `(:c, 1)`): represent **dynamic variables** at specific time shifts ($z_t, k_{t-1}, c_{t+1}$).
- **`TimedVar`s**: `[TimedVar(:k, -1), TimedVar(:c, 1)]`.

### Grouped Arguments
You can supply nested groups in `wrt`. Dynspec computes the full Jacobian in a single forward pass and returns a `Tuple` of matrices:

```julia
wrt_grouped = [
    [(:z, 0), (:k, 0)],    # Group 1: current variables (y_t)
    [(:c, 1), (:y, 1)]     # Group 2: lead variables (y_{t+1})
]
J_curr, J_lead = jacobian(block; wrt=wrt_grouped)
# J_curr: 5×2 Matrix{Float64}
# J_lead: 5×2 Matrix{Float64}
```

### Standard Differentiation
```julia
# 1. Differentiate all equations with respect to current endogenous variables:
endo = [(:z, 0), (:k, 0), (:y, 0), (:i, 0), (:c, 0)]
J = jacobian(block; wrt=endo) # 5×5 Matrix{Float64}

# 2. Parameter sensitivity w.r.t bare symbols:
J_params = jacobian(block; wrt=[:α, :β, :δ]) # 5×3 Matrix{Float64}

# 3. Differentiate a filtered block:
J_trans = jacobian(block; wrt=[(:z, 0), (:k, -1), (:i, 0)], filter_by=:transition) # 2×3 Matrix{Float64}
```

### Simultaneous Residuals & Jacobians (`derivatives`)
```julia
r, J = derivatives(block; wrt=[(:c, 0), (:c, 1)], filter_by=:arbitrage)
# r: Vector{Float64} of length 1
# J: 1×2 Matrix{Float64}

# With grouped wrt:
r, (J1, J2) = derivatives(block; wrt=wrt_grouped)
```

---

## 8. DSGE Perturbation Matrices & Helpers (`Dyno` Submodule)

All DSGE perturbation and model inspection helpers are strictly organized under the `Dyno` submodule:

```julia
using Dynspec
using Dynspec.Dyno
# or access as Dyno.<function>
```

### `Dyno.endogenous_exogenous(block; [endogenous], [exogenous])`
Extract the ordered lists of endogenous and exogenous variables:
```julia
vars = Dyno.endogenous_exogenous(block)
vars.endogenous   # [:z, :k, :y, :i, :c]
vars.exogenous    # [:e_z]

# Or with tuple destructuring:
endo, exo = Dyno.endogenous_exogenous(block)
```

### `Dyno.compute_dsge_jacobians(block; [filter_by], [endogenous], [exogenous], [at])`

> [!NOTE]
> `Dyno.compute_dsge_jacobians` is a convenience helper for first-order DSGE perturbation solvers ($A \cdot y_{t+1} + B \cdot y_t + C \cdot y_{t-1} + D \cdot e_t = 0$).

```julia
r, A, B, C, D = Dyno.compute_dsge_jacobians(block;
    endogenous = [:z, :k, :y, :i, :c],
    exogenous  = [:e_z]
)

# Returns:
# r: Residuals vector (length 5)
# A: ∂f / ∂y_{t+1} (5×5 Matrix{Float64})
# B: ∂f / ∂y_t     (5×5 Matrix{Float64})
# C: ∂f / ∂y_{t-1} (5×5 Matrix{Float64})
# D: ∂f / ∂e_t     (5×1 Matrix{Float64})
```

If `endogenous` and `exogenous` are omitted, `Dyno` automatically infers them from declared processes and equations. Under the hood, this helper calls `derivatives(block; wrt=[...])` using grouped arguments.

### `Dyno.model_properties(block)`
Inspect all model components:
```julia
props = Dyno.model_properties(block)
# Dict with :name, :parameters, :steady_states, :processes,
# :variables, :dynamic_variables, :endogenous, :exogenous, :equations_count, and :tags
```

---

---

## 9. Interactive Pluto Notebooks (`notebooks/`)

Dynspec provides an optional sub-environment in `notebooks/` with **Pluto.jl**, **PlutoUI**, and **PlotlyLight** (a lightweight plotting package with sub-second TTFX and responsive browser-side interactive charts).

### Starting the Pluto Server
```bash
# Using the executable runner script:
./notebooks/run.jl

# Or via Julia CLI:
julia --project=notebooks notebooks/run.jl
```

### Included Interactive Notebook
The default notebook `notebooks/interactive_model.jl` showcases:
- Real-time parameter tweaking via PlutoUI sliders (`α`, `β`, `δ`).
- Dynamic steady-state recalculation and model parsing.
- Interactive response and sensitivity charts using `PlotlyLight.Plot`.

---

## Summary Table of Functions & Options

| Function | Primary Arguments | Key Options / Keywords | Return Type | Description |
| :--- | :--- | :--- | :--- | :--- |
| `parse_model` | `text::AbstractString` | `filename="<string>"` | `ModelBlock` | Parse DSL string into a model block |
| `read_model` | `filepath::AbstractString` | — | `ModelBlock` | Read `.dyno` model file |
| `model_properties` | `block::AbstractModelBlock` | — | `Dict{Symbol, Any}` | Overview of parameters, variables, shocks, tags |
| `list_variables` | `block` or `equation` | `shifts=false` | `Vector{Symbol}` or `Dict` | List variables (and optional time shifts) |
| `endogenous_exogenous` | `block` | `endogenous=nothing`, `exogenous=nothing` | `NamedTuple` | Extract lists of endogenous & exogenous variables |
| `residuals` | `block::AbstractModelBlock` | `filter_by=nothing`, `at=nothing` | `Vector{Float64}` | Evaluate equation residuals directly on AST |
| `jacobian` | `block::AbstractModelBlock` | `wrt`, `filter_by=nothing`, `at=nothing` | `Matrix` or `Tuple{Matrix...}` | Compute analytical Jacobian matrix / matrices |
| `derivatives` | `block::AbstractModelBlock` | `wrt`, `filter_by=nothing`, `at=nothing` | `(Vector, Matrix)` or `(Vector, Tuple)` | Simultaneous residuals and Jacobian(s) |
| `compute_dsge_jacobians` | `block::AbstractModelBlock` | `filter_by=nothing`, `endogenous=nothing`, `exogenous=nothing`, `at=nothing` | `(r, A, B, C, D)` | First-order perturbation matrices |
| `equations` | `block::AbstractModelBlock` | — | `Vector{Equation}` | Access equation list |
| `filter_equations` | `block, tag` | — | `Vector{Equation}` | Filter equations by metadata tag |

