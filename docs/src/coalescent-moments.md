# Selectable multi-deme coalescent equations

`solve_coalescent_odes(tfgy, samples; states=...)` integrates analytical
coalescent moments for any number of demes. Unlike the two-host
`ZoonosisParams` method, it allocates only the selected groups and their
dependencies. These are **mean-field closure equations**, not exact
expectations of the finite-population stochastic process.

## Usage and state sizes

```@example moments
using Phylodynamics

F = [3.0 1.0 0.5; 0.2 2.0 0.7; 0.4 0.1 1.0]
G = zeros(3, 3)
Y = [50.0, 40.0, 30.0]
demography = TFGY(
    [2.0, 0.0], [F, F], [G, G], [Y, Y],
    [zeros(3), zeros(3)], ["wildlife", "livestock", "human"])

# Three ODE variables: no parsimony closure or topology states are constructed.
minimal = solve_coalescent_odes(demography, [6.0, 4.0, 3.0])

# Directed, weighted, three-character Sankoff parsimony.
result = solve_coalescent_odes(demography, [6.0, 4.0, 3.0];
    states=(:jumps, :coalescences, :branch_lengths, :cherries, :sackin, :sankoff),
    costs=[0 1 2; 3 0 1; 2 3 0],
    max_states=10_000, max_parsimony_classes=256)

result.summary.sankoff_score
result.summary.reached_mrca
result.layout.groups
```

```@example moments
@assert minimal.layout.nstates == 3
@assert sum(result.summary.sankoff; dims=2)[:] ≈ result.summary.lineages
(lineages=result.summary.lineages, weighted_score=result.summary.sankoff_score,
    ode_states=result.layout.nstates)
```

For a Catalyst/ModelingToolkit reaction system, use
`solve_coalescent_odes(rs, infected, u0, tspan, p, samples; states=..., res=1000)`.
The ordering matches `TreeProblem`; `samples` can be a count vector in infected
species order or a map such as `[Y1 => 20, Y2 => 10, Y3 => 5]`.
The reaction system is converted to a `TFGY` trajectory first.

Let ``m`` be the number of demographic demes and ``k`` the number of observed
character states. `character_states` maps each deme to a label in `1:k`;
the default is `1:m`. For example, `[1, 2, 2]` retains three demographic demes
but groups the last two into one observed character state. This changes the
character being reconstructed, not the demographic model.

| Selected group | Allocated ODE variables |
|:---------------|:-----------------------|
| `:lineages` (always included) | ``m`` |
| `:jumps` | ``m^2`` (diagonal remains zero) |
| `:coalescences` | ``m^2`` |
| `:branch_lengths` | ``m`` |
| `:cherries` | ``mk + k(k+1)/2`` (includes external lineages) |
| `:sackin` | ``m+1`` (includes descendant masses) |
| `:fitch` | ``m(2^k-1)+1`` |
| `:camin_sokal` | Two binary directional systems, ``2(2m+1)``; requires ``k=2`` |
| `:sankoff` | ``mq+1``, where ``q`` is the reachable cost-message count |
| `:parsimony_lengths` | An additional ``q`` for each selected parsimony program |

Parsimony class lengths are summed over demographic demes; this includes the
two-host pure-parsimony branch-length statistic as a special case. Request at
least one parsimony program alongside `:parsimony_lengths`.

Parsimony is opt-in. The total ODE budget and the per-automaton class budget are
checked before integration; exceeding either throws an error rather than
truncating, merging, or approximating classes. A ``q``-class merge table has
``q^2`` entries and its RHS costs ``O(m^2q^2)`` per evaluation. Grouping observed
states is often more effective than raising the budgets.

The result includes `sol` (the ODE solution), `traj` (only allocated variables),
`layout.groups` (flat state ranges), and `layout.programs` (parsimony messages,
leaf classes, merge tables, and score increments). Count matrices use deme rows
and class columns, with Julia's column-major order in the flat state ranges.
Cherry matrices duplicate each unordered off-diagonal count symmetrically;
do not sum both triangles to obtain total cherries.

## Demographic and sampling assumptions

Time ``s=T_{\mathrm{sample}}-t`` runs backwards. ``F_{ij}`` is forward
transmission from a parent in deme ``i`` to a child in deme ``j``.
``G_{ij}`` is forward physical migration from ``i`` to ``j``.
``Y_i>0`` is the infected population size and ``A_i`` the ancestral lineage
count. All samples are homochronous at `T_sample`; sample counts must not
exceed the corresponding populations, and their total must be greater than one.

`TFGY.times` must be strictly decreasing. Rates and populations are linearly
interpolated between its grid points. The selected interval must have strictly
positive populations, even if an excluded earlier part of the trajectory starts
at zero; use `t0` to restrict the interval when necessary. Death rates affect the
forward populations but do not produce separate events on surviving ancestral
lineages.

The probabilistic event interpretation additionally requires ``0\leq A_i/Y_i\leq1``.
The numerical solution is checked against this domain at saved times, allowing
solver-scale roundoff; a violation throws `DomainError` rather than returning
invalid occupancy probabilities. This grid check is not a continuous-time proof
that the domain is preserved.
Factorization of lineage occupancy, and of descendant mass within a deme, are
closure assumptions. Algebraic identities below do not prove those assumptions
or establish exact finite-population expectations.

Integration stops at ``\sum_i A_i=1`` if reached. Otherwise it ends at `t0` and
`summary.reached_mrca` is false: the statistics then describe a partial ancestral
forest, not a complete tree. The older `ZoonosisParams` method still uses its
original nearest-saved-time summary convention.

## Event-derived equations

Define ``p_i=A_i/Y_i`` and ``B_{ij}=F_{ij}(1-p_i)+G_{ij}``.
A transmission with both endpoints ancestral has rate ``F_{ij}p_ip_j``:
two lineages merge into one in parent deme ``i``. If only the child endpoint
is ancestral, its lineage moves backwards from ``j`` to ``i``.
Physical migration also moves it backwards, without coalescence.

### Lineages, jumps, coalescences and lengths

```math
\dot A_i =
-\sum_j F_{ji}p_jp_i
+\sum_{j\ne i}B_{ij}p_j
-\sum_{j\ne i}B_{ji}p_i,
\qquad
\sum_i\dot A_i=-\sum_{ij}F_{ij}p_ip_j.
```

Directional true jumps include both cross-deme transmission and physical
migration; ``J_{ii}=0``. Coalescences include within-deme events.

```math
\dot J_{ij}=(F_{ij}+G_{ij})p_j\quad(i\ne j),\qquad
\dot Q_{ij}=F_{ij}p_ip_j,\qquad
\dot{\mathcal L}_i=A_i.
```

### External lineages and cherries

Let ``L_{ia}`` be external lineages currently in deme ``i`` whose sampled
character label was ``a``, and ``\ell_{ia}=L_{ia}/Y_i``.
Coalescence removes both external endpoints; their merged ancestor is internal.

```math
\dot L_{ia}
=-\ell_{ia}\left(\sum_jF_{ij}p_j+\sum_jF_{ji}p_j\right)
+\sum_{j\ne i}B_{ij}\ell_{ja}
-\sum_{j\ne i}B_{ji}\ell_{ia}.
```

For unordered cherry labels ``a\le b``,

```math
\dot C_{ab}=\sum_{ij}F_{ij}
\left(\ell_{ia}\ell_{jb}+
\mathbf 1_{a\ne b}\ell_{ib}\ell_{ja}\right).
```

### Descendant masses and Sackin's index

``X_i`` is the total descendant-tip mass carried by lineages in deme ``i``,
and ``x_i=X_i/A_i`` is the mean mass per lineage (zero when both vanish).
At coalescence the two descendant masses add in the parent deme.

```math
\dot X_i=
\sum_{j\ne i}(F_{ij}+G_{ij})p_jx_j
-\sum_{j\ne i}(F_{ji}+G_{ji})p_ix_i,\qquad
\sum_i\dot X_i=0,
```

```math
\dot K=\sum_{ij}F_{ij}p_ip_j(x_i+x_j).
```

## General directed Sankoff parsimony

The `costs` matrix specifies **edge costs**, not shortest-path distances:
``C_{ab}`` is the cost of parent character ``a`` and child character ``b``.
Costs must be nonnegative machine-representable integers or `Inf`, with zero
diagonal. Asymmetry and forbidden transitions are supported. Integer additions
are checked for overflow. Internal message tables encode `Inf` as `typemax(Int)`.
The root is free;
the accumulated score is a weighted parsimony cost, not necessarily a number
of changes. Integer arithmetic makes class identification exact.

For a subtree with root-conditioned costs ``v``, retain the normalized
**stem message**

```math
a_i=\min_j(C_{ij}+v_j)-\min_jv_j.
```

A leaf observed in state ``b`` has message ``C_{\cdot b}``.
To merge messages ``a,b``, set

```math
t_j=a_j+b_j,\qquad
d(a,b)=\min_jt_j,\qquad
M(a,b)_i=\min_j(C_{ij}+t_j)-d(a,b).
```

The subtree score increases by ``d(a,b)`` and its new class is ``M(a,b)``.
This is the ordinary Sankoff recursion with the subtree's common cost
subtracted; keeping only an argmin label would lose information for asymmetric
costs. Zero diagonal and nonnegative costs imply
``\min_i\min_j(C_{ij}+t_j)=\min_jt_j``, so no root/stem charge is lost.
Induction on a binary tree therefore gives the ordinary free-root Sankoff score.

Starting from leaf messages, the implementation enumerates their reachable merge
closure. For finite costs bounded by ``c``, every normalized message entry is an
integer in ``[0,c]``, so the closure is finite (a coarse bound is
``(c+1)^k-c^k``). With forbidden transitions this finite bound no longer applies;
the compiler rejects incompatible subtrees or stops at the class budget. No
finite-closure claim is made for every matrix containing `Inf`.

For each retained message class ``a``, let ``H_{ia}`` count its lineages in deme
``i`` and ``h_{ia}=H_{ia}/Y_i``. Enumerate ordered parent/child classes:

```math
(i,a),(j,b)\longrightarrow(i,M(a,b))
\quad\hbox{at rate }F_{ij}h_{ia}h_{jb}.
```

Subtract both input classes, add the merged class, and increase the parsimony
score at that rate times ``d(a,b)``. Child-only ancestry and migration transfer
class ``b`` from ``j`` to ``i`` at rate ``B_{ij}h_{jb}`` for ``i\ne j``.
Summing class derivatives over ``a`` recovers ``\dot A_i`` exactly on
``\sum_aH_{ia}=A_i``.

Unit off-diagonal costs recover Fitch's intersection/union rule and the
``2^k-1`` nonempty optimal-state subsets. Costs `[0 1; Inf 0]` and
`[0 Inf; 1 0]` recover the two binary directional Camin-Sokal priors.

## Verification and the corrected legacy equations

The Symbolics tests invoke the **actual general RHS** with symbolic state
variables, specialize to two demes, and project onto all 37 legacy derivatives.
They compare with the separately written two-host RHS, under positive-population
and class-partition assumptions. Additional three-deme identities verify class
row sums, coalescent loss, and descendant-mass conservation.

This comparison exposed a pre-existing error in **both directional** legacy
systems: their class row sums exceeded ``\dot A_1`` and ``\dot A_2`` by
``F_{12}(p_{11}p_{22}+p_{12}p_{21})`` and
``F_{21}(p_{11}p_{22}+p_{12}p_{21})``, respectively (with analogous reverse
labels). The erroneous coalescence terms are corrected. Fitch, true jumps,
cherries and Sackin equations are unchanged.

Consequently, old directional-parsimony scores and pure-parsimony branch-length
summaries must be regenerated; preserving those erroneous numerical results is
not a compatibility goal. The historical source is pinned by commit and SHA-256
in `verification/source.json`, and the counterexample remains an executable check.

Run the linked checks from the package root:

```sh
julia --project=. test/coalescent_symbolics.jl --export
python3 -m venv verification/.venv
verification/.venv/bin/python -m pip install -r verification/requirements.txt
verification/.venv/bin/python verification/verify_coalescent_equations.py \
  --generated verification/generated-equations.csv
julia --project=. test/runtests.jl
```

The local-venv SymPy verifier checks independent event identities, the exported
Julia equations, and normalized Sankoff scores against exhaustive ancestral-state
assignments on small trees. Julia tests also check the compiled automata against
ordinary tree dynamic programming, analytic one-deme solutions, migration-only
transport, relabeling, state selection and budget failures.
These checks establish algebraic consistency and the tested implementation
properties; they are not a Lean proof or a proof of the stochastic closure.
