# API Reference

```@meta
CurrentModule = Dynspec
```

## Model Construction & Types

```@docs
AbstractModelBlock
ModelBlock
Equation
Context
Normal
TimedVar
Complementarity
equations
```

## Model Ingestion

```@docs
parse_model
read_model
```

## Analysis & Evaluation

```@docs
list_variables
list_constants
residuals
jacobian
derivatives
filter_equations
```

## Dyno Submodule (DSGE Perturbation)

```@docs
Dyno
Dyno.first_order_model
Dyno.endogenous_exogenous
Dyno.compute_dsge_jacobians
Dyno.solve_dsge
Dyno.simulate_linear
Dyno.compile_dynamic
```
