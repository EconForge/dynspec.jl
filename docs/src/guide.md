# Function & Option Guide

This guide walks through the user-facing functions, types, and workflows in **Dynspec.jl**.

---

## Part I: User-Facing API

### 1. Model Ingestion

Dynspec models can be parsed from strings or read from `.dyno` specification files.

#### `parse_model`
Parse a model directly from a specification string:
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

#### `read_model`
Read and interpret a `.dyno` model file from disk:
```julia
block = read_model("test/models/neo.dyno")
```

Both return an instance of `ModelBlock <: AbstractModelBlock`.

---

### 2. Context & Parameters

The model `context` holds all parameters, calibrated steady states, processes, time-dated values, and metadata as standard typed dictionaries.

#### Structured Dictionary Access
```julia
block.context.constants       # Dict{Symbol, Any} of parameters: Dict(:α => 0.36, :β => 0.96, ...)
block.context.steady_states   # Dict{Symbol, Any} of steady states: Dict(:k => 5.447, :y => 1.841, ...)
block.context.processes       # Dict{Symbol, Any} of shock distributions: Dict(:e_z => Normal(...))
block.context.values          # Dict{Symbol, Dict{Int, Any}} of time-dated values: values[:k][0]
block.context.metadata        # Dict{Symbol, Any} of metadata: Dict(:name => "Neoclassical", ...)
```

Direct indexing provides explicit, unambiguous access:
```julia
block.context.constants[:α]       # 0.36
block.context.steady_states[:k]   # 5.4468
```

> [!NOTE]
> Values are stored in standard dictionaries without shorthand property forwarding (e.g. `block.context.α`), avoiding namespace collisions between dictionary fields and model identifiers.

---

### 3. Symbols: Variables & Constants

Dynspec provides clean, complementary functions to inspect the symbols in a model:

#### `list_constants(block)`
List all constant/parameter symbols defined in the block context:
```julia
list_constants(block)
# Vector{Symbol}: [:α, :β, :γ, :δ, :ρ]
```

#### `list_variables(block; [shifts=false])`
List dynamic variable symbols occurring across equations:
```julia
# 1. Unique variable names:
list_variables(block)
# Vector{Symbol}: [:c, :e_z, :i, :k, :y, :z]

# 2. Variables with their dynamic time shifts:
list_variables(block; shifts=true)
# Dict{Symbol, Set{Int}} mapping each symbol to its observed time lags/leads:
#   :c   => Set([0, 1])
#   :k   => Set([-1, 0])
#   :y   => Set([0, 1])
#   :z   => Set([-1, 0])
#   :e_z => Set([0])
```

*(Note: Domain-specific partitioning of variables into endogenous and exogenous sets is handled by the [`Dyno` submodule](dsge.md)).*

---

### 4. Inspecting Equations

To inspect the equations of a model block, use the standard `equations` accessor:

```julia
eqs = equations(block)
```

Each `Equation` object exposes user-facing properties:
- `eq.raw::String`: Original textual formula (e.g. `"1 - β*(c[t+1]/c[t])^(-1)..."`).
- `eq.metadata::Dict{Symbol, Any}`: Tags (e.g. `[:arbitrage]`) and annotations.
- `eq.complementarity::Union{Nothing, Complementarity}`: Normalized complementarity conditions (`lb <= x[t] <= ub`).

```julia
eq = equations(block)[5]
println(eq.raw)               # "1 - β*(c[t+1]/c[t])^(-1)*(1-δ + α*y[t+1]/k[t]) ⟂ 0 <= c[t] <= inf"
println(eq.metadata[:tags])   # [:arbitrage]
println(eq.complementarity)   # Complementarity(c[t], 0.0, Inf)
```

> [!WARNING]
> #### ⚠️ Design Note: `equations(block)[i]` vs `block[i]`
> A `ModelBlock` encapsulates both **assigned objects** (parameters, steady-state calibrations, processes in `context`) and **equations**.
> In principle, indexing `block[1]` could refer to the first assignment, the first instruction, or the first equation.
> While `block[i]` currently delegates to `block.equations[i]` for convenience, this behavior is open to discussion.
> **We strongly recommend using `equations(block)[i]` or `block.equations[i]`** to make your intent completely unambiguous.

#### Filtering Equations
Filter equations by block tag:
```julia
filter_equations(block, :transition)
filter_equations(block, "arbitrage")
```

---

### 5. Residual Evaluation (`residuals`)

Equation residuals are evaluated directly on parsed equations without code generation or slow compilation:

```julia
# 1. Steady-state residuals for all equations:
residuals(block)

# 2. Filter by tag:
residuals(block; filter_by=:transition)
residuals(block; filter_by=[:transition, :arbitrage])

# 3. Filter by custom predicate:
residuals(block; filter_by=eq -> :transition in get(eq.metadata, :tags, Symbol[]))
```

#### The `at` Option (Open Design / Undecided)
You can evaluate residuals at custom points using the `at` keyword:
```julia
residuals(block; at=Dict(:z => 0.1, :k => 5.2))
residuals(block; at=(z=0.1, k=5.2))
```

> [!NOTE]
> #### ❓ Semantics of `at` (Undecided / Subject to Change)
> By default (when `at=nothing`), equations are evaluated at the calibrated steady states (`block.context.steady_states`) and parameters (`block.context.constants`).
> When `at` is supplied, it overrides variable or parameter values.
> **Open Question**: Should `at` represent steady-state levels ($\bar{y}$), dynamic time-dated state vectors ($y_t, y_{t-1}$), or general overrides?
> This signature is currently experimental and may be refined into distinct keywords (e.g. `steady_state=...`, `point=...`) in future versions.

---

### 6. Analytical Derivatives & Jacobians (`jacobian` & `derivatives`)

Exact analytical derivatives are computed via `ForwardDiff.Dual` numbers:

#### `jacobian(block; wrt=...)`
Compute the Jacobian matrix (or matrices):

- **Bare Symbols** (`:α`, `:β`): parameters/constants.
- **Time-dated tuples** (`(:z, 0)`, `(:k, -1)`, `(:c, 1)`): dynamic variables at specific time shifts.

```julia
# 1. Jacobian w.r.t current variables:
endo = [(:z, 0), (:k, 0), (:y, 0), (:i, 0), (:c, 0)]
J = jacobian(block; wrt=endo) # 5×5 Matrix{Float64}

# 2. Parameter sensitivity Jacobian:
J_params = jacobian(block; wrt=[:α, :β, :δ]) # 5×3 Matrix{Float64}

# 3. Filtered equations:
J_trans = jacobian(block; wrt=[(:z, 0), (:k, -1), (:i, 0)], filter_by=:transition)
```

#### Grouped `wrt` Arguments
Pass nested groups in `wrt` to compute multiple Jacobian blocks simultaneously in a single forward pass:
```julia
wrt_grouped = [
    [(:z, 0), (:k, 0)],    # Group 1: current variables (y_t)
    [(:c, 1), (:y, 1)]     # Group 2: lead variables (y_{t+1})
]
J_curr, J_lead = jacobian(block; wrt=wrt_grouped)
```

#### `derivatives(block; wrt=...)` & The Taylor Expansion Hierarchy
`derivatives` returns both residuals and Jacobian(s) simultaneously:
```julia
r, J = derivatives(block; wrt=[(:c, 0), (:c, 1)], filter_by=:arbitrage)
# With grouped wrt:
r, (J_curr, J_lead) = derivatives(block; wrt=wrt_grouped)
```

> [!TIP]
> #### 💡 Why does `derivatives` return residuals?
> In Dynspec, `derivatives` is structured around **orders of expansion** around the reference point $x_0$:
> - **Order 0**: Residuals $f(x_0)$
> - **Order 1**: Jacobian $\nabla f(x_0)$
> - **Higher Orders (future)**: Hessians $\nabla^2 f(x_0)$, third-order tensors, etc.
> 
> Evaluating order 0 and order 1 together in `derivatives` ensures that both the baseline evaluation point and the sensitivity matrix are computed consistently. If you only need the first-order sensitivity matrix, use `jacobian(block; wrt=...)`.

---

## Part II: Internal Architecture & Developer API

### Preserved AST Parse Trees

For compiler developers, symbolic backends, and advanced analysis tools, Dynspec preserves the full grammar AST:

```julia
# Full Lark parse tree of the model block (Lerche.Tree)
block.tree

# Individual equation parse tree node
eq = equations(block)[1]
eq.tree
eq.residual  # AST node of the residual expression (rhs - lhs)
```

Preserving the AST avoids re-lexing and re-parsing when calibrating parameters, running variant sweeps, or exporting to other modeling engines.

---

## Interactive Pluto Notebooks (`notebooks/`)

Dynspec provides an optional sub-environment in `notebooks/` with **Pluto.jl**, **PlutoUI**, and **PlotlyLight**:

```bash
# Using the executable runner script:
./notebooks/run.jl

# Or via Julia CLI:
julia --project=notebooks notebooks/run.jl
```

Open `notebooks/interactive_model.jl` to explore parameter sliders, real-time steady-state recalculation, and interactive sensitivity plots.

---

## Summary Table of Functions & Options

| Function | Primary Arguments | Key Options / Keywords | Return Type | Description |
| :--- | :--- | :--- | :--- | :--- |
| `parse_model` | `text::AbstractString` | `filename="<string>"` | `ModelBlock` | Parse DSL string into a model block |
| `read_model` | `filepath::AbstractString` | — | `ModelBlock` | Read `.dyno` model file |
| `list_constants` | `block::AbstractModelBlock` | — | `Vector{Symbol}` | List all constant/parameter names |
| `list_variables` | `block` or `equation` | `shifts=false` | `Vector{Symbol}` or `Dict` | List variables (and optional time shifts) |
| `equations` | `block::AbstractModelBlock` | — | `Vector{Equation}` | Access equation collection |
| `filter_equations` | `block, tag` | — | `Vector{Equation}` | Filter equations by metadata tag |
| `residuals` | `block::AbstractModelBlock` | `filter_by=nothing`, `at=nothing` | `Vector{Float64}` | Evaluate residuals directly on AST |
| `jacobian` | `block::AbstractModelBlock` | `wrt`, `filter_by=nothing`, `at=nothing` | `Matrix` or `Tuple{Matrix...}` | Compute analytical Jacobian(s) |
| `derivatives` | `block::AbstractModelBlock` | `wrt`, `filter_by=nothing`, `at=nothing` | `(r, J)` or `(r, (J1, ...))` | Simultaneous 0th-order residuals & 1st-order Jacobians |
