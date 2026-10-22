module Parser

using Lerche
using ..Types

export get_parser, parse_tree, is_variable_node, extract_var_info, normalize_complementarity, to_string, to_expr, collect_variables!

const GRAMMAR_PATH = joinpath(dirname(@__DIR__), "src", "grammar.lark")

function _load_grammar()
    txt = read(GRAMMAR_PATH, String)
    return Lark(txt; start=["free_block", "equation_block", "formula"], parser="lalr")
end

const _PARSER_CACHE = Ref{Union{Nothing, Lark}}(nothing)

function get_parser()
    if _PARSER_CACHE[] === nothing
        _PARSER_CACHE[] = _load_grammar()
    end
    return _PARSER_CACHE[]
end

function parse_tree(text::AbstractString; start::String="free_block")
    p = get_parser()
    # ensure trailing newline for Lark line grammar
    txt = endswith(text, "\n") ? text : text * "\n"
    return Lerche.parse(p, txt; start=start)
end

function is_variable_node(tree)
    return tree isa Tree && tree.data == "variable"
end

function extract_var_info(tree::Tree)
    name = Symbol(tree.children[1].children[1].value)
    idx_token = tree.children[2].children[1].value
    if idx_token == "~"
        return TimedVar(name, :~)
    end
    shift = 0
    if length(tree.children) >= 3 && tree.children[3] !== nothing
        c3 = tree.children[3]
        str_shift = c3 isa Token ? c3.value : (c3 isa Tree ? c3.children[1].value : string(c3))
        shift = parse(Int, str_shift)
    end
    return TimedVar(name, shift)
end

function _is_current(tree::Tree)
    v = extract_var_info(tree)
    return v.shift == 0
end

function collect_variables!(res::Dict{Symbol, Set{Int}}, node)
    if node isa Tree
        if node.data == "variable"
            v = extract_var_info(node)
            if v.shift isa Int
                if !haskey(res, v.name)
                    res[v.name] = Set{Int}()
                end
                push!(res[v.name], v.shift)
            end
        else
            for child in node.children
                collect_variables!(res, child)
            end
        end
    end
    return res
end

"""
    normalize_complementarity(tree::Tree)

Normalize `res ⟂ lb <= x <= ub` (or one-sided bounds) into `(residual_tree, Complementarity(var, lb, ub))`.
"""
function normalize_complementarity(tree::Tree)
    @assert tree.data == "complementarity"
    children = tree.children
    res = children[1]
    if length(children) == 4 && children[4] !== nothing
        lb_node = children[2]
        x_node = children[3]
        ub_node = children[4]
        is_variable_node(x_node) || error("The middle term of lb <= x <= ub must be a variable")
        x_var = extract_var_info(x_node)
        return res, Complementarity(x_var, lb_node, ub_node)
    else
        a = children[2]
        b = children[3]
        is_a_var = is_variable_node(a)
        is_b_var = is_variable_node(b)

        bounded_is_rhs = false
        if is_a_var && is_b_var
            _is_current(a) == _is_current(b) && error("Ambiguous one-sided bound: both sides of <= are variables; write lb <= x <= ub explicitly")
            bounded_is_rhs = _is_current(b)
        elseif is_a_var || is_b_var
            bounded_is_rhs = is_b_var
        else
            error("The constraint must bound a variable (lb <= x[t] or x[t] <= ub)")
        end

        if bounded_is_rhs
            x_var = extract_var_info(b)
            return res, Complementarity(x_var, a, Inf)
        else
            x_var = extract_var_info(a)
            return res, Complementarity(x_var, -Inf, b)
        end
    end
end

"""
    to_string(tree_or_token)

Convert an AST node into a clean formula string.
"""
function to_string(node)
    if node isa Token
        return node.value
    elseif node isa Number
        return string(node)
    elseif !(node isa Tree)
        return string(node)
    end

    data = node.data
    children = node.children

    if data == "constant"
        return children[1].children[1].value
    elseif data == "variable"
        v = extract_var_info(node)
        return string(v)
    elseif data == "value"
        name = children[1].children[1].value
        time = children[2].children[1].value
        return "$(name)[$(time)]"
    elseif data == "number"
        return children[1].value
    elseif data == "add"
        return "$(to_string(children[1])) + $(to_string(children[2]))"
    elseif data == "sub"
        return "$(to_string(children[1])) - ($(to_string(children[2])))"
    elseif data == "mul"
        return "($(to_string(children[1]))) * ($(to_string(children[2])))"
    elseif data == "div"
        return "($(to_string(children[1]))) / ($(to_string(children[2])))"
    elseif data == "pow"
        return "($(to_string(children[1])))^($(to_string(children[2])))"
    elseif data == "neg"
        return "-($(to_string(children[1])))"
    elseif data == "call"
        fname = children[1].children[1].value
        args = join([to_string(c) for c in children[2:end]], ", ")
        return "$(fname)($(args))"
    elseif data == "equality"
        return "$(to_string(children[1])) = $(to_string(children[2]))"
    elseif data == "bare_formula"
        return to_string(children[1])
    elseif data == "complementarity"
        res, comp = normalize_complementarity(node)
        return "$(to_string(res)) ⟂ $(to_string(comp.lb)) <= $(comp.variable) <= $(to_string(comp.ub))"
    elseif data == "assignment"
        return "$(to_string(children[1])) = $(to_string(children[2]))"
    else
        return join([to_string(c) for c in children], " ")
    end
end

function Base.show(io::IO, c::Types.Complementarity)
    lb_str = c.lb isa Number ? (c.lb == -Inf ? "-inf" : (c.lb == Inf ? "inf" : string(c.lb))) : to_string(c.lb)
    ub_str = c.ub isa Number ? (c.ub == -Inf ? "-inf" : (c.ub == Inf ? "inf" : string(c.ub))) : to_string(c.ub)
    print(io, lb_str, " <= ", c.variable, " <= ", ub_str)
end

function Base.show(io::IO, eq::Types.Equation)
    print(io, eq.raw)
end

"""
    to_expr(tree_or_token)

Convert an AST node into a Julia `Expr`.
"""
function to_expr(node)
    if node isa Token
        val = node.value
        n = tryparse(Int, val)
        n !== nothing && return n
        f = tryparse(Float64, val)
        f !== nothing && return f
        return Symbol(val)
    elseif node isa Number
        return node
    elseif !(node isa Tree)
        return node
    end

    data = node.data
    children = node.children

    if data == "constant"
        sym = Symbol(children[1].children[1].value)
        return sym == :inf ? Inf : sym
    elseif data == "variable"
        v = extract_var_info(node)
        shift_expr = v.shift == :~ ? :~ : (v.shift == 0 ? :t : Expr(:call, v.shift > 0 ? :+ : :-, :t, abs(v.shift)))
        return Expr(:ref, v.name, shift_expr)
    elseif data == "value"
        name = Symbol(children[1].children[1].value)
        time = parse(Int, children[2].children[1].value)
        return Expr(:ref, name, time)
    elseif data == "number"
        val = children[1].value
        n = tryparse(Int, val)
        n !== nothing && return n
        return parse(Float64, val)
    elseif data == "add"
        return Expr(:call, :+, to_expr(children[1]), to_expr(children[2]))
    elseif data == "sub"
        return Expr(:call, :-, to_expr(children[1]), to_expr(children[2]))
    elseif data == "mul"
        return Expr(:call, :*, to_expr(children[1]), to_expr(children[2]))
    elseif data == "div"
        return Expr(:call, :/, to_expr(children[1]), to_expr(children[2]))
    elseif data == "pow"
        return Expr(:call, :^, to_expr(children[1]), to_expr(children[2]))
    elseif data == "neg"
        return Expr(:call, :-, to_expr(children[1]))
    elseif data == "call"
        fname = Symbol(children[1].children[1].value)
        args = [to_expr(c) for c in children[2:end]]
        return Expr(:call, fname, args...)
    elseif data == "equality"
        # Residual form: rhs - lhs
        return Expr(:call, :-, to_expr(children[2]), to_expr(children[1]))
    elseif data == "bare_formula"
        return to_expr(children[1])
    elseif data == "complementarity"
        res, _ = normalize_complementarity(node)
        return to_expr(res)
    else
        error("Cannot convert AST node $(data) to Expr")
    end
end

to_expr(eq::Equation) = to_expr(eq.residual)

end # module Parser
