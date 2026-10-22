# Dynspec.jl

!!! warning "⚠️ Work in Progress / Active Development"
    **Dynspec.jl is under active initial development.**
    APIs, syntax representations, and internal interfaces are experimental and subject to breaking changes.
    Follow progress, report issues, and contribute on [GitHub (EconForge/dynspec.jl)](https://github.com/EconForge/dynspec.jl).

A concise, type-stable Julia package for reading, interpreting, and analyzing dynamic economic model specifications, replicating and extending the `dynspec` engine of [dyno](https://github.com/EconForge/dyno.py).

## Core Principles

1. **Unified `AbstractModelBlock` Architecture**:
   Models are represented by `ModelBlock <: AbstractModelBlock`, holding a calibrated type-stable `Context`, an iterable collection of `Equation`s, metadata, and the complete preserved grammar parse tree (`Tree`).
2. **Behavioral Interface**:
   A model block is defined by its ability to evaluate `residuals` and analytical `jacobian` / `derivatives`, rather than merely being a vector of equations.
3. **Preserved Parse Tree**:
   The full AST is preserved on both the block and individual equations. This avoids re-lexing or re-parsing when calibrating parameters, running variant sweeps, or targeting alternate symbolic backends.
4. **Decoupled Equation Filtering (`filter_by`)**:
   Residual and derivative evaluations accept `filter_by` to filter equations by block tags (`:transition`, `:arbitrage`, etc.) or custom predicates.
5. **Exact Analytical Differentiation via ForwardDiff**:
   Computes exact machine-precision derivatives using `ForwardDiff.Dual` numbers without slow symbolic engines or custom derivative dictionaries.
6. **Clean, Predictable Context**:
   Model contexts hold typed dictionaries (`constants`, `steady_states`, `processes`, `values`), providing unambiguous, explicit access (`context.constants[:α]`) and avoiding namespace collisions.

## Quick Installation

```julia
using Pkg
Pkg.add(url="https://github.com/EconForge/dynspec.jl")
```

See the [Function & Option Guide](guide.md) for a comprehensive demonstration of all package functions and options.
