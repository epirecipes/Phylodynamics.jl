# ==============================================================================
# types_and_trees.jl
# PhyloTree representation, Newick / ReMASTER extended Newick
# parser & writer, singleton collapsing, cherry assortativity, and Sackin index.
# ==============================================================================

"""
    PhyloTree

Rooted phylogenetic tree with integer tip states indexing arbitrary demes,
optional internal node deme states, and node times. The `true_jumps_12`,
`true_jumps_21` and corresponding time fields are pairwise diagnostics for
states 1 and 2, not a restriction on the number or identity of demes.

Conventions match APE `phylo`:
- Tips are numbered `1:ntips`.
- Root is numbered `ntips + 1`.
- Other internal (and optional singleton) nodes are `ntips + 2 : ntips + nnode`.
- `edge` is an `nedges × 2` matrix where `edge[e, 1]` is parent and `edge[e, 2]` is child.
"""
mutable struct PhyloTree
    edge::Matrix{Int}
    edge_length::Vector{Float64}
    tip_label::Vector{String}
    tip_state::Vector{Int}
    node_state::Vector{Int}
    node_time::Vector{Float64}
    nnode::Int
    true_jumps_12::Int
    true_jumps_21::Int
    true_jump_times_12::Vector{Float64}
    true_jump_times_21::Vector{Float64}
end

function PhyloTree(
    edge::Matrix{Int},
    edge_length::Vector{Float64},
    tip_label::Vector{String},
    tip_state::Vector{Int};
    node_state::Vector{Int} = Int[],
    node_time::Vector{Float64} = Float64[],
    nnode::Int = maximum(edge) - length(tip_label),
    true_jumps_12::Int = 0,
    true_jumps_21::Int = 0,
    true_jump_times_12::Vector{Float64} = Float64[],
    true_jump_times_21::Vector{Float64} = Float64[]
)
    ntips = length(tip_label)
    ntotal = ntips + nnode
    if isempty(node_state)
        node_state = zeros(Int, ntotal)
        node_state[1:ntips] .= tip_state
    end
    if isempty(node_time)
        node_time = compute_node_times_from_root(edge, edge_length, ntips, nnode)
    end
    return PhyloTree(
        edge, edge_length, tip_label, tip_state, node_state, node_time,
        nnode, true_jumps_12, true_jumps_21, true_jump_times_12, true_jump_times_21
    )
end

"""Return the number of sampled tips in a `PhyloTree`."""
ntips(tree::PhyloTree) = length(tree.tip_label)
"""Return the number of tips plus internal nodes in a `PhyloTree`."""
ntotal_nodes(tree::PhyloTree) = length(tree.tip_label) + tree.nnode
"""Return the root index (`ntips(tree) + 1`) under the package's APE numbering."""
root_node(tree::PhyloTree) = length(tree.tip_label) + 1

function compute_node_times_from_root(edge::Matrix{Int}, edge_length::Vector{Float64}, ntips::Int, nnode::Int)
    ntotal = ntips + nnode
    children_map = [Tuple{Int,Float64}[] for _ in 1:ntotal]
    has_parent = falses(ntotal)
    for e in 1:size(edge, 1)
        p = edge[e, 1]
        c = edge[e, 2]
        push!(children_map[p], (c, edge_length[e]))
        has_parent[c] = true
    end
    r = findfirst(!, has_parent[(ntips+1):ntotal])
    root = r === nothing ? (ntips + 1) : (ntips + r)
    times = zeros(Float64, ntotal)
    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        for (v, bl) in children_map[u]
            times[v] = times[u] + bl
            push!(stack, v)
        end
    end
    return times
end

"""
    infer_tip_state(label::AbstractString, meta_type::AbstractString="") -> Int

Infer host deme state (1 = animal reservoir, 2 = human) from a tip label or
ReMASTER `[&type="..."]` metadata annotation.
"""
function infer_tip_state(label::AbstractString, meta_type::AbstractString="")
    if !isempty(meta_type)
        state = tryparse(Int, meta_type)
        state !== nothing && state > 0 && return state
        if occursin("Y2", meta_type) || occursin("sampleY", meta_type) || meta_type == "Y" || meta_type == "2"
            return 2
        elseif occursin("Y1", meta_type) || occursin("sampleI", meta_type) || meta_type == "I" || meta_type == "1"
            return 1
        end
    end
    if endswith(label, "_I") || occursin("_Y1", label) || occursin("sampleI", label) || startswith(label, "I")
        return 1
    elseif endswith(label, "_Y") || occursin("_Y_", label) || occursin("sampleY", label) || occursin("Y2", label) || startswith(label, "Y")
        return 2
    else
        return 1
    end
end

# Internal AST node for parsing Newick / extended Newick strings
mutable struct _NewickAST
    label::String
    meta_type::String
    meta_time::Float64
    branch_length::Float64
    children::Vector{_NewickAST}
end

"""
    read_newick(newick_str::AbstractString; T_sample::Float64=12.85) -> PhyloTree

Parse a standard Newick or ReMASTER extended Newick (`[&type="Y1",time=...]`) string
into a `PhyloTree`. If internal/singleton nodes carry `[&type=...]` annotations,
computes exact ground-truth directional jumps (`true_jumps_12`, `true_jumps_21`).
"""
function read_newick(newick_str::AbstractString; T_sample::Float64=12.85)
    s = strip(String(newick_str))
    isempty(s) && error("Empty Newick string")
    if startswith(lowercase(s), "tree ")
        eq_idx = findfirst('=', s)
        if eq_idx !== nothing
            s = strip(s[(eq_idx+1):end])
        end
    end
    if endswith(s, ";")
        s = chop(s)
    end
    chars = collect(s)
    pos = Ref(1)
    nchars = length(chars)

    function skip_ws()
        while pos[] <= nchars && isspace(chars[pos[]])
            pos[] += 1
        end
    end

    function parse_metadata()
        mtype = ""
        mtime = NaN
        skip_ws()
        if pos[] <= nchars && chars[pos[]] == '['
            start_idx = pos[]
            while pos[] <= nchars && chars[pos[]] != ']'
                pos[] += 1
            end
            meta_str = String(chars[start_idx:min(pos[], nchars)])
            if pos[] <= nchars && chars[pos[]] == ']'
                pos[] += 1
            end
            m1 = match(r"type=\"?([^\",\]\s]+)\"?", meta_str)
            if m1 !== nothing
                mtype = m1.captures[1]
            end
            m2 = match(r"time=([0-9.eE+-]+)", meta_str)
            if m2 !== nothing
                mtime = tryparse(Float64, m2.captures[1]) === nothing ? NaN : parse(Float64, m2.captures[1])
            end
        end
        return mtype, mtime
    end

    function parse_label()
        skip_ws()
        buf = Char[]
        while pos[] <= nchars && !(chars[pos[]] in ('(', ')', ',', ':', '[', ';'))
            push!(buf, chars[pos[]])
            pos[] += 1
        end
        return strip(String(buf))
    end

    function parse_Subtree()
        skip_ws()
        children = _NewickAST[]
        if pos[] <= nchars && chars[pos[]] == '('
            pos[] += 1 # consume '('
            while true
                push!(children, parse_Subtree())
                skip_ws()
                if pos[] <= nchars && chars[pos[]] == ','
                    pos[] += 1
                elseif pos[] <= nchars && chars[pos[]] == ')'
                    pos[] += 1
                    break
                else
                    break
                end
            end
        end
        lbl = parse_label()
        mtype1, mtime1 = parse_metadata()
        bl = 0.0
        skip_ws()
        if pos[] <= nchars && chars[pos[]] == ':'
            pos[] += 1
            mtype2, mtime2 = parse_metadata()
            if isempty(mtype1) && !isempty(mtype2)
                mtype1 = mtype2
            end
            if isnan(mtime1) && !isnan(mtime2)
                mtime1 = mtime2
            end
            num_str = parse_label()
            bl = isempty(num_str) ? 0.0 : parse(Float64, num_str)
            mtype3, mtime3 = parse_metadata()
            if isempty(mtype1) && !isempty(mtype3)
                mtype1 = mtype3
            end
            if isnan(mtime1) && !isnan(mtime3)
                mtime1 = mtime3
            end
        end
        return _NewickAST(lbl, mtype1, mtime1, bl, children)
    end

    ast_root = parse_Subtree()

    # Collect leaves and internal nodes
    leaves = _NewickAST[]
    internals = _NewickAST[]
    function collect_nodes(node::_NewickAST)
        if isempty(node.children)
            push!(leaves, node)
        else
            push!(internals, node)
            for c in node.children
                collect_nodes(c)
            end
        end
    end
    collect_nodes(ast_root)

    n_tips = length(leaves)
    n_int = length(internals)
    n_total = n_tips + n_int
    node_id = IdDict{_NewickAST, Int}()
    tip_label = Vector{String}(undef, n_tips)
    tip_state = Vector{Int}(undef, n_tips)
    node_state = zeros(Int, n_total)
    node_time = fill(NaN, n_total)

    for (i, lf) in enumerate(leaves)
        node_id[lf] = i
        st = infer_tip_state(lf.label, lf.meta_type)
        lbl = isempty(lf.label) ? "tip_$i" : lf.label
        tip_label[i] = lbl
        tip_state[i] = st
        node_state[i] = st
        node_time[i] = lf.meta_time
    end

    # Root gets n_tips + 1, remaining internals get n_tips + 2 ... n_total
    for (j, nd) in enumerate(internals)
        id = n_tips + j
        node_id[nd] = id
        if !isempty(nd.meta_type)
            node_state[id] = infer_tip_state("", nd.meta_type)
        end
        node_time[id] = nd.meta_time
    end

    edges_p = Int[]
    edges_c = Int[]
    edge_len = Float64[]
    function collect_edges(node::_NewickAST)
        p = node_id[node]
        for c in node.children
            cid = node_id[c]
            push!(edges_p, p)
            push!(edges_c, cid)
            push!(edge_len, c.branch_length)
            collect_edges(c)
        end
    end
    collect_edges(ast_root)

    edge = hcat(edges_p, edges_c)
    rel_times = compute_node_times_from_root(edge, edge_len, n_tips, n_int)
    max_rel = maximum(rel_times[1:n_tips])
    root_t = any(!isnan, node_time) ? minimum(filter(!isnan, node_time)) : (T_sample - max_rel)
    for u in 1:n_total
        if isnan(node_time[u])
            node_time[u] = root_t + rel_times[u]
        end
    end

    # Count true jumps if internal node states are annotated
    jumps_12 = 0
    jumps_21 = 0
    jtimes_12 = Float64[]
    jtimes_21 = Float64[]
    for e in 1:size(edge, 1)
        p = edge[e, 1]
        c = edge[e, 2]
        sp = node_state[p]
        sc = node_state[c]
        if sp > 0 && sc > 0 && sp != sc
            if sp == 1 && sc == 2
                jumps_12 += 1
                push!(jtimes_12, node_time[p])
            elseif sp == 2 && sc == 1
                jumps_21 += 1
                push!(jtimes_21, node_time[p])
            end
        end
    end

    return PhyloTree(
        edge, edge_len, tip_label, tip_state, node_state, node_time,
        n_int, jumps_12, jumps_21, sort!(jtimes_12), sort!(jtimes_21)
    )
end

"""
    read_newick_file(path::AbstractString; T_sample::Float64=12.85) -> Vector{PhyloTree}

Read all Newick or Nexus trees from `path`.
"""
function read_newick_file(path::AbstractString; T_sample::Float64=12.85)
    lines = readlines(path)
    trees = PhyloTree[]
    for line in lines
        s = strip(line)
        isempty(s) && continue
        if startswith(lowercase(s), "tree ")
            eq_idx = findfirst('=', s)
            eq_idx === nothing && continue
            s = strip(s[(eq_idx+1):end])
        end
        if startswith(s, "(")
            push!(trees, read_newick(s; T_sample=T_sample))
        end
    end
    return trees
end

"""
    collapse_singletons(tree::PhyloTree) -> PhyloTree

Collapse degree-2 singleton nodes (such as ReMASTER transmission/migration reactions along a branch)
into a strictly bifurcating tree while preserving total branch lengths, true node states at
coalescent events, and true directional jump counts (`true_jumps_12`, `true_jumps_21`).
"""
function collapse_singletons(tree::PhyloTree)
    n_tips = ntips(tree)
    n_tot = ntotal_nodes(tree)
    children_map = [Tuple{Int,Float64}[] for _ in 1:n_tot]
    has_parent = falses(n_tot)
    for e in 1:size(tree.edge, 1)
        p = tree.edge[e, 1]
        c = tree.edge[e, 2]
        push!(children_map[p], (c, tree.edge_length[e]))
        has_parent[c] = true
    end
    r_idx = findfirst(!, has_parent[(n_tips+1):n_tot])
    orig_root = r_idx === nothing ? (n_tips + 1) : (n_tips + r_idx)

    # Descend past any degree-1 root chain
    while length(children_map[orig_root]) == 1 && orig_root > n_tips
        orig_root = children_map[orig_root][1][1]
    end

    # Collect kept internal nodes (degree >= 2)
    kept_internals = Int[]
    function find_kept(u::Int)
        if u > n_tips
            if length(children_map[u]) >= 2
                push!(kept_internals, u)
            end
            for (v, _) in children_map[u]
                find_kept(v)
            end
        end
    end
    find_kept(orig_root)

    n_int_new = length(kept_internals)
    n_tot_new = n_tips + n_int_new
    old_to_new = Dict{Int,Int}()
    new_node_state = zeros(Int, n_tot_new)
    new_node_time = zeros(Float64, n_tot_new)

    for i in 1:n_tips
        old_to_new[i] = i
        new_node_state[i] = tree.node_state[i]
        new_node_time[i] = tree.node_time[i]
    end
    for (j, old_u) in enumerate(kept_internals)
        new_u = n_tips + j
        old_to_new[old_u] = new_u
        new_node_state[new_u] = tree.node_state[old_u]
        new_node_time[new_u] = tree.node_time[old_u]
    end

    edges_p = Int[]
    edges_c = Int[]
    edge_len = Float64[]

    for old_u in kept_internals
        p_new = old_to_new[old_u]
        for (child, bl0) in children_map[old_u]
            curr = child
            acc_bl = bl0
            while curr > n_tips && length(children_map[curr]) == 1
                curr, bl_next = children_map[curr][1]
                acc_bl += bl_next
            end
            c_new = old_to_new[curr]
            push!(edges_p, p_new)
            push!(edges_c, c_new)
            push!(edge_len, acc_bl)
        end
    end

    return PhyloTree(
        hcat(edges_p, edges_c),
        edge_len,
        copy(tree.tip_label),
        copy(tree.tip_state),
        new_node_state,
        new_node_time,
        n_int_new,
        tree.true_jumps_12,
        tree.true_jumps_21,
        copy(tree.true_jump_times_12),
        copy(tree.true_jump_times_21)
    )
end

"""
    write_newick(tree::PhyloTree; include_metadata::Bool=false) -> String

Serialize a `PhyloTree` to a Newick string (optionally including `[&type=...,time=...]` metadata).
"""
function write_newick(tree::PhyloTree; include_metadata::Bool=false)
    n_tips = ntips(tree)
    n_tot = ntotal_nodes(tree)
    children_map = [Tuple{Int,Float64}[] for _ in 1:n_tot]
    has_parent = falses(n_tot)
    for e in 1:size(tree.edge, 1)
        p = tree.edge[e, 1]
        c = tree.edge[e, 2]
        push!(children_map[p], (c, tree.edge_length[e]))
        has_parent[c] = true
    end
    r_idx = findfirst(!, has_parent[(n_tips+1):n_tot])
    root = r_idx === nothing ? (n_tips + 1) : (n_tips + r_idx)

    function fmt_node(u::Int, bl::Union{Nothing,Float64})
        meta = ""
        if include_metadata && tree.node_state[u] > 0
            st_str = string(tree.node_state[u])
            meta = @sprintf("[&type=\"%s\",time=%.6f]", st_str, tree.node_time[u])
        end
        if u <= n_tips
            lbl = tree.tip_label[u] * meta
            return bl === nothing ? lbl : @sprintf("%s:%.6f", lbl, bl)
        else
            sub = join([fmt_node(c, c_bl) for (c, c_bl) in children_map[u]], ",")
            s = "($sub)$meta"
            return bl === nothing ? s : @sprintf("%s:%.6f", s, bl)
        end
    end

    return fmt_node(root, nothing) * ";"
end

"""
    write_newick_file(tree::PhyloTree, path; include_metadata=false)

Write a Newick tree to `path` and return that path. Set `include_metadata=true`
to include node states and times, as in [`write_newick`](@ref).
"""
function write_newick_file(tree::PhyloTree, path::AbstractString; include_metadata::Bool=false)
    open(path, "w") do io
        println(io, write_newick(tree; include_metadata=include_metadata))
    end
    return path
end

"""
    postorder_edges(tree::PhyloTree) -> Vector{Int}

Return edge row indices in postorder (children visited before parents).
"""
function postorder_edges(tree::PhyloTree)
    n_tips = ntips(tree)
    n_tot = ntotal_nodes(tree)
    child_edges = [Int[] for _ in 1:n_tot]
    has_parent = falses(n_tot)
    for e in 1:size(tree.edge, 1)
        p = tree.edge[e, 1]
        c = tree.edge[e, 2]
        push!(child_edges[p], e)
        has_parent[c] = true
    end
    r_idx = findfirst(!, has_parent[(n_tips+1):n_tot])
    root = r_idx === nothing ? (n_tips + 1) : (n_tips + r_idx)

    order = Int[]
    sizehint!(order, size(tree.edge, 1))
    function visit(u::Int)
        for e in child_edges[u]
            c = tree.edge[e, 2]
            visit(c)
            push!(order, e)
        end
    end
    visit(root)
    return order
end

"""
    node_heights(tree::PhyloTree) -> Vector{Float64}

Compute node distances from the root (`height[root] == 0.0`).
"""
function node_heights(tree::PhyloTree)
    return compute_node_times_from_root(tree.edge, tree.edge_length, ntips(tree), tree.nnode)
end

"""
    node_depths(tree::PhyloTree) -> Vector{Int}

Compute topological depth (number of edges from the root) for each node.
"""
function node_depths(tree::PhyloTree)
    n_tips = ntips(tree)
    n_tot = ntotal_nodes(tree)
    children_map = [Int[] for _ in 1:n_tot]
    has_parent = falses(n_tot)
    for e in 1:size(tree.edge, 1)
        p = tree.edge[e, 1]
        c = tree.edge[e, 2]
        push!(children_map[p], c)
        has_parent[c] = true
    end
    r_idx = findfirst(!, has_parent[(n_tips+1):n_tot])
    root = r_idx === nothing ? (n_tips + 1) : (n_tips + r_idx)
    depths = zeros(Int, n_tot)
    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        for v in children_map[u]
            depths[v] = depths[u] + 1
            push!(stack, v)
        end
    end
    return depths
end

"""
    cherry_assortativity(tree::PhyloTree) -> NamedTuple{(:C11, :C12, :C22, :r)}

Compute cherry counts `(C11, C12, C22)` and Newman's assortativity coefficient `r`
for a binary `PhyloTree`.
"""
function cherry_assortativity(tree::PhyloTree)
    bt = collapse_singletons(tree)
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    children_map = [Int[] for _ in 1:n_tot]
    for e in 1:size(bt.edge, 1)
        push!(children_map[bt.edge[e, 1]], bt.edge[e, 2])
    end
    C11 = 0
    C12 = 0
    C22 = 0
    for u in (n_tips+1):n_tot
        ch = children_map[u]
        if length(ch) == 2 && ch[1] <= n_tips && ch[2] <= n_tips
            s1 = bt.tip_state[ch[1]]
            s2 = bt.tip_state[ch[2]]
            if s1 == 1 && s2 == 1
                C11 += 1
            elseif s1 == 2 && s2 == 2
                C22 += 1
            else
                C12 += 1
            end
        end
    end
    denom = (C11 + 0.5 * C12) * (C22 + 0.5 * C12)
    r = denom > 0 ? (C11 * C22 - 0.25 * C12^2) / denom : 0.0
    return (C11=C11, C12=C12, C22=C22, r=r)
end

"""
    sackin_index(tree::PhyloTree) -> NamedTuple{(:IS, :IS_norm, :E_yule)}

Compute raw Sackin's index `IS` and Yule-normalized Sackin's index `IS_norm`.
"""
function sackin_index(tree::PhyloTree)
    bt = collapse_singletons(tree)
    n = ntips(bt)
    depths = node_depths(bt)
    IS = Float64(sum(depths[1:n]))
    harm = sum(1.0 / j for j in 2:n)
    E_yule = 2.0 * n * harm
    IS_norm = (IS - E_yule) / E_yule
    return (IS=IS, IS_norm=IS_norm, E_yule=E_yule)
end
