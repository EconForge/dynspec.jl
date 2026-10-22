# Dynspec.jl

[![Status](https://img.shields.io/badge/status-work--in--progress-orange.svg)](https://github.com/EconForge/dynspec.jl)
[![CI](https://github.com/EconForge/dynspec.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/EconForge/dynspec.jl/actions)
[![Documentation](https://github.com/EconForge/dynspec.jl/actions/workflows/docs.yml/badge.svg)](https://econforge.github.io/dynspec.jl/dev/)

> [!WARNING]
> ## ⚠️ Work in Progress
> **Dynspec.jl is under active development.**
> Syntax, data structures, and APIs are experimental and subject to breaking changes.

A concise and elegant Julia package for reading, interpreting, and analyzing dynamic economic model specifications, replicating the `dynspec` engine of [dyno](https://github.com/EconForge/dyno.py).

## Features

- **Model File Interpretation**: Parses `.dyno` specification files and DSL strings.
- **Context Extraction**:
  - `constants`: parameters and sequentially evaluated constants (`Dict{Symbol, Any}`).
  - `steady_states`: calibrated steady states (`x[~] <- ...`) and process steady-state means.
  - `processes`: stochastic shocks and processes (`e[t] <- N(0, sig^2)`).
  - `values`: explicit time-dated values (`k[0] <- ...`) and quantified ranges (`∀ t, 0 <= t < T : ...`).
  - `metadata`: top-level `@key: value` declarations (e.g., `@name`, `@run`, `@variants`).
- **Equation & Block Extraction**:
  - Residual equations evaluated directly on the parsed syntax tree.
  - Block tags (e.g., `[definitions]`, `[transition]`, `[arbitrage]`) and inline statement metadata / labels.
  - Complementarity conditions (`res ⟂ lb <= x[t] <= ub` with automatic normalization for one-sided bounds).
  - Collection methods (`length(block)`, `block[i]`, `equations(block)`).
- **Unified ModelBlock & AST Retention**:
  - `ModelBlock <: AbstractModelBlock` encapsulates `context`, `equations`, `metadata`, and preserves the complete parsed AST (`tree`).
  - Equations retain their individual parse tree nodes (`block[i].tree`).
- **Model Analysis & Differentiation**:
  - `list_variables(block)`: lists variable names (and optional time shifts with `shifts=true`).
  - `list_constants(block)`: lists parameter/constant symbols defined in context.
  - `residuals(block; filter_by=..., at=...)`: evaluates equation residuals directly on AST trees.
  - `jacobian(block; wrt=..., filter_by=..., at=...)`: computes exact analytical Jacobians with `ForwardDiff` dual numbers; supports bare symbols as constants/parameters, time-dated variables, and grouped `wrt`.
  - `derivatives(block; wrt=..., filter_by=..., at=...)`: computes 0th-order residuals and 1st-order Jacobians simultaneously.
  - `Dyno.first_order_model(block)`: unified first-order DSGE model NamedTuple with symbols, residuals, and Jacobians `(A, B, C, D)`.

## Quick Start

```julia
using Dynspec

# Read and interpret a model file (returns a ModelBlock)
block = read_model("test/models/neo.dyno")

# Inspect context (typed dictionaries)
block.context.constants[:α]       # 0.387
block.context.steady_states[:k]   # 5.152
block.context.processes[:e_z]     # Normal(...)
block.context.metadata[:name]     # "Neoclassical"

# Access preserved parse tree
block.tree                        # Lark Tree AST

# Compute residuals of all equations or a filtered subset (direct AST evaluation):
residuals(block)                          # [0.0, 0.0, 0.0, 0.0, 0.0]
residuals(block; filter_by=:transition)   # only equations tagged :transition
residuals(block; filter_by=[:transition, :arbitrage])
residuals(block; filter_by=eq -> :transition in get(eq.metadata, :tags, Symbol[]))

# Compute Jacobians on all or filtered equations:
# wrt can be:
# - Tuples for dynamic variables: (:z, 0), (:k, -1)
# - Bare Symbols for parameter sensitivity: :α, :β
# - Grouped vectors: [ [(:z, 0), (:k, 0)], [(:c, 1), (:y, 1)] ] -> returns (J1, J2)
endo_vars = [(:z, 0), (:k, 0), (:y, 0), (:i, 0), (:c, 0)]
J_endo = jacobian(block; wrt=endo_vars)
J_params = jacobian(block; wrt=[:α, :β, :δ])
J_trans = jacobian(block; wrt=[(:z, 0), (:k, -1), (:i, 0)], filter_by=:transition)

# Simultaneous residuals and Jacobian evaluation
r, J = derivatives(block; wrt=[(:c, 0), (:c, 1)], filter_by=:arbitrage)

# List variables and constants
list_variables(block)             # [:c, :e_y, :e_z, :i, :k, :y, :z]
list_constants(block)             # [:α, :β, :γ, :δ, :ρ]
list_variables(block; shifts=true)# Dict(:k => Set([-1, 0]), :c => Set([0, 1]), ...)

# DSGE first-order model representation (via Dyno submodule)
using Dynspec.Dyno
fo = Dyno.first_order_model(block)
fo.symbols.endogenous             # [:z, :k, :y, :i, :c]
fo.jacobians.A                    # 5×5 Matrix{Float64}
```

## Running Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation

To build the static HTML documentation once:
```bash
julia --project=docs docs/make.jl
```

To live-serve the documentation with automatic real-time rebuild and browser refresh on edits:
```bash
julia --project=docs docs/serve.jl
```

## Interactive Pluto Notebooks

An optional sub-environment is provided under `notebooks/` equipped with **Pluto.jl**, **PlutoUI**, and **PlotlyLight** (an ultra-fast interactive plotting engine with sub-second load times):

```bash
# Launch Pluto with the preconfigured environment
./notebooks/run.jl

# or using Julia CLI:
julia --project=notebooks notebooks/run.jl
```

Open the default interactive notebook `interactive_model.jl` to explore parameter sliders, real-time steady-state calculations, and interactive sensitivity curves.
