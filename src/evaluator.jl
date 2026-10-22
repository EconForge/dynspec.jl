module Evaluator

using Lerche
using YAML
using ..Types
using ..Parser

export interpret_model, eval_formula

const MATH_FUNCTIONS = Dict{Symbol, Function}(
    :exp => exp,
    :log => log,
    :sqrt => sqrt,
    :abs => abs,
    :sin => sin,
    :cos => cos,
    :tan => tan,
    :asin => asin,
    :acos => acos,
    :atan => atan,
    :sinh => sinh,
    :cosh => cosh,
    :tanh => tanh,
    :min => min,
    :max => max,
    :N => Types.N,
    :Normal => (args...) -> Types.Normal(args...)
)

"""
    eval_formula(node, context::Union{Context, ContextBuilder};
                 steady_state::Bool=false,
                 time::Union{Nothing, Int}=nothing,
                 env::Union{Nothing, AbstractDict}=nothing,
                 constants::Union{Nothing, AbstractDict}=nothing)

Evaluate a formula AST node directly within the given context.
"""
function eval_formula(node, context::Union{Context, ContextBuilder};
                      steady_state::Bool=false,
                      time::Union{Nothing, Int}=nothing,
                      env::Union{Nothing, AbstractDict}=nothing,
                      constants::Union{Nothing, AbstractDict}=nothing)
    if node isa Token
        val = node.value
        n = tryparse(Int, val)
        n !== nothing && return n
        f = tryparse(Float64, val)
        f !== nothing && return f
        return val
    elseif node isa Number
        return node
    elseif !(node isa Tree)
        return node
    end

    data = node.data
    children = node.children

    if data == "number"
        val = children[1].value
        n = tryparse(Int, val)
        n !== nothing && return n
        return parse(Float64, val)

    elseif data == "constant"
        name = Symbol(children[1].children[1].value)
        if constants !== nothing && haskey(constants, name)
            return constants[name]
        elseif name == :inf || name == :Inf
            return Inf
        elseif haskey(context.constants, name)
            return context.constants[name]
        else
            return NaN
        end

    elseif data == "variable"
        v = extract_var_info(node)
        if steady_state || v.shift == :~
            return get(context.steady_states, v.name, NaN)
        end

        shift = if time !== nothing
            time + (v.shift isa Int ? v.shift : 0)
        else
            v.shift isa Int ? v.shift : 0
        end

        if env !== nothing && haskey(env, (v.name, shift))
            return env[(v.name, shift)]
        elseif time !== nothing && haskey(context.values, v.name) && haskey(context.values[v.name], shift)
            return context.values[v.name][shift]
        elseif haskey(context.steady_states, v.name)
            return context.steady_states[v.name]
        else
            return NaN
        end

    elseif data == "value"
        name = Symbol(children[1].children[1].value)
        t_val = parse(Int, children[2].children[1].value)
        if env !== nothing && haskey(env, (name, t_val))
            return env[(name, t_val)]
        elseif haskey(context.values, name) && haskey(context.values[name], t_val)
            return context.values[name][t_val]
        else
            return NaN
        end

    elseif data == "add"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants) +
               eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "sub"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants) -
               eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "mul"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants) *
               eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "div"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants) /
               eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "pow"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants) ^
               eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "neg"
        return -eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "call"
        fname = Symbol(children[1].children[1].value)
        if fname == :steady_state
            return eval_formula(children[2], context; steady_state=true, time=nothing, env=env, constants=constants)
        elseif haskey(MATH_FUNCTIONS, fname)
            args = [eval_formula(c, context; steady_state=steady_state, time=time, env=env, constants=constants) for c in children[2:end]]
            return MATH_FUNCTIONS[fname](args...)
        else
            error("Undefined function: $(fname)")
        end

    elseif data == "equality"
        # Equation residual: rhs - lhs
        r = eval_formula(children[2], context; steady_state=steady_state, time=time, env=env, constants=constants)
        l = eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants)
        return r - l

    elseif data == "bare_formula"
        return eval_formula(children[1], context; steady_state=steady_state, time=time, env=env, constants=constants)

    elseif data == "complementarity"
        res_tree, _ = normalize_complementarity(node)
        return eval_formula(res_tree, context; steady_state=steady_state, time=time, env=env, constants=constants)

    else
        error("Unknown node in formula evaluation: $(data)")
    end
end

function _parse_yaml_scalar(text::AbstractString)
    clean = strip(text)
    try
        return YAML.load(clean)
    catch
        return clean
    end
end

function _parse_annotation(node::Tree)
    meta = Dict{Symbol, Any}(:tags => Symbol[])
    for entry in node.children
        if entry.data == "kv"
            k = Symbol(entry.children[1].children[1].value)
            v_node = entry.children[2]
            v = if v_node isa Tree
                _parse_yaml_scalar(v_node.children[1].value)
            elseif v_node isa Token
                _parse_yaml_scalar(v_node.value)
            else
                _parse_yaml_scalar(string(v_node))
            end
            meta[k] = v
        elseif entry.data == "baretag"
            tag = Symbol(entry.children[1].children[1].value)
            push!(meta[:tags], tag)
        elseif entry.data == "strtag"
            raw_s = entry.children[1].value
            s = strip(raw_s, ['"', '\''])
            push!(meta[:tags], Symbol(s))
        end
    end
    return meta
end

function _parse_statement_metadata(node::Tree)
    child = node.children[1]
    if child isa Tree
        return _parse_annotation(child)
    else
        # META_TEXT token
        raw = strip(string(child))
        if startswith(raw, "\"") && endswith(raw, "\"") || startswith(raw, "'") && endswith(raw, "'")
            return Dict{Symbol, Any}(:label => strip(raw, ['"', '\'']), :tags => Symbol[])
        else
            # Comma-separated tags
            tags = [Symbol(strip(part)) for part in split(raw, ",") if !isempty(strip(part))]
            return Dict{Symbol, Any}(:tags => tags)
        end
    end
end

function _merge_metadata(parent::Dict{Symbol, Any}, child::Dict{Symbol, Any})
    res = copy(parent)
    parent_tags = copy(get(parent, :tags, Symbol[]))
    child_tags = copy(get(child, :tags, Symbol[]))
    merged_tags = unique(vcat(parent_tags, child_tags))
    res[:tags] = merged_tags

    for (k, v) in child
        k == :tags && continue
        res[k] = v
    end
    return res
end

function _process_assignment!(tree::Tree, context::ContextBuilder)
    sym_node = tree.children[1]
    val_node = tree.children[2]
    val = eval_formula(val_node, context)

    data = sym_node.data
    if data == "constant"
        name = Symbol(sym_node.children[1].children[1].value)
        context.constants[name] = val
    elseif data == "variable"
        var = extract_var_info(sym_node)
        if var.shift == :~
            context.steady_states[var.name] = val
        else
            # Process / shock
            context.processes[var.name] = val
            if val isa Normal
                context.steady_states[var.name] = val.Μ[1]
            end
        end
    elseif data == "value"
        name = Symbol(sym_node.children[1].children[1].value)
        t = parse(Int, sym_node.children[2].children[1].value)
        if !haskey(context.values, name)
            context.values[name] = Dict{Int, Any}()
        end
        context.values[name][t] = val
    end
    return val
end

function _process_quantified_assignment!(tree::Tree, context::ContextBuilder)
    bounds_node = tree.children[1]
    @assert bounds_node.data == "t_double_bound"
    lower = Int(eval_formula(bounds_node.children[1], context))
    upper = Int(eval_formula(bounds_node.children[2], context))

    sym_node = tree.children[2]
    rhs_node = tree.children[3]
    var = extract_var_info(sym_node)

    if !haskey(context.values, var.name)
        context.values[var.name] = Dict{Int, Any}()
    end

    for t in lower:(upper - 1)
        context.constants[:t] = t
        val = eval_formula(rhs_node, context; time=t)
        context.values[var.name][t] = val
    end
    delete!(context.constants, :t)
end

function _process_block_item!(item::Tree, context::ContextBuilder, equations::Vector{Equation}, meta_stack::Vector{Dict{Symbol, Any}})
    if item.data == "model_metadata"
        k = Symbol(item.children[1].value)
        raw_val = string(item.children[2].value)
        mute = false
        if endswith(strip(raw_val), ";")
            mute = true
            raw_val = rstrip(strip(raw_val), ';')
        end
        val = _parse_yaml_scalar(raw_val)

        if k == :run
            cmd = if val isa String
                Dict{Symbol, Any}(:command => Symbol(val), :options => Dict{Symbol, Any}(), :mute => mute)
            elseif val isa Dict
                if haskey(val, "command")
                    Dict{Symbol, Any}(:command => Symbol(val["command"]), :options => get(val, "options", Dict{Symbol, Any}()), :mute => mute)
                else
                    k_cmd = first(keys(val))
                    Dict{Symbol, Any}(:command => Symbol(k_cmd), :options => val[k_cmd] === nothing ? Dict{Symbol, Any}() : val[k_cmd], :mute => mute)
                end
            else
                Dict{Symbol, Any}(:command => Symbol(val), :options => Dict{Symbol, Any}(), :mute => mute)
            end

            if haskey(context.metadata, :run)
                curr = context.metadata[:run]
                if curr isa Vector
                    push!(curr, cmd)
                else
                    context.metadata[:run] = [curr, cmd]
                end
            else
                context.metadata[:run] = [cmd]
            end
        else
            context.metadata[k] = val
            if mute
                context.metadata[Symbol("_muted_$(k)")] = true
            end
        end

    elseif item.data == "annotated_statement"
        stmt = item.children[1]
        stmt_meta = length(item.children) > 1 ? _parse_statement_metadata(item.children[2]) : Dict{Symbol, Any}()
        merged = _merge_metadata(meta_stack[end], stmt_meta)

        if stmt.data in ("equality", "bare_formula", "complementarity")
            res_tree, comp = if stmt.data == "complementarity"
                normalize_complementarity(stmt)
            else
                (stmt, nothing)
            end
            raw = to_string(stmt)
            push!(equations, Equation(res_tree, raw, merged, comp, stmt))
        elseif stmt.data == "assignment"
            _process_assignment!(stmt, context)
        elseif stmt.data == "quantified_assignment"
            _process_quantified_assignment!(stmt, context)
        end

    elseif item.data == "annotated_block"
        block_meta = _parse_annotation(item.children[1].children[1])
        merged = _merge_metadata(meta_stack[end], block_meta)
        push!(meta_stack, merged)
        block_node = item.children[2]
        for child in block_node.children
            if child isa Tree
                _process_block_item!(child, context, equations, meta_stack)
            end
        end
        pop!(meta_stack)
    end
end

"""
    interpret_model(tree::Tree, filename::String="<string>", source::String="") -> ModelBlock

Interpret an AST parsed from a DynSpec model file, extracting context and equations.
"""
function interpret_model(tree::Tree, filename::String="<string>", source::String="")
    builder = ContextBuilder()
    equations = Equation[]
    meta_stack = [Dict{Symbol, Any}(:tags => Symbol[])]

    for item in tree.children
        if item isa Tree
            _process_block_item!(item, builder, equations, meta_stack)
        end
    end

    context = build_context(builder)
    return ModelBlock(context, equations; metadata=meta_stack[1], filename=filename, source=source, tree=tree)
end

end # module Evaluator
