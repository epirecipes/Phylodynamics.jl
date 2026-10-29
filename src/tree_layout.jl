# Helper to compute 2D coordinates (x = time from root, y = vertical tip layout) for plotting a PhyloTree in CairoMakie
function _tree_layout_coords(tree_in::PhyloTree)
    bt = collapse_singletons(tree_in)
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    root = n_tips + 1
    children_map = [Int[] for _ in 1:n_tot]
    for e in 1:size(bt.edge, 1)
        push!(children_map[bt.edge[e, 1]], bt.edge[e, 2])
    end

    x = node_heights(bt)
    y = zeros(Float64, n_tot)
    leaf_counter = Ref(0)
    function assign_y(u::Int)
        if u <= n_tips || isempty(children_map[u])
            leaf_counter[] += 1
            y[u] = Float64(leaf_counter[])
        else
            for c in children_map[u]
                assign_y(c)
            end
            y[u] = mean(y[c] for c in children_map[u])
        end
    end
    assign_y(root)
    return bt, x, y
end
