#!/usr/bin/env python3
"""Fail-closed symbolic event audit and exact finite-tree Sankoff checks.

Run: verification/.venv/bin/python verification/verify_coalescent_equations.py
Add --generated verification/generated-equations.csv to check Julia exports.
Only the pinned historical source is read for the baseline discrepancy audit.
"""

import argparse
import csv
import hashlib
import io
import itertools
import json
import re
import subprocess
import sys
import tokenize
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CHECKS = 0


def equal(actual, expected=0):
    global CHECKS
    if actual == expected:
        CHECKS += 1
        return
    residual = sp.factor(actual - expected)
    assert residual == 0, f"Nonzero residual: {residual}"
    CHECKS += 1


def symbols(prefix, n, m=None):
    if m is None:
        return [sp.Symbol(f"{prefix}{i+1}", positive=True) for i in range(n)]
    return [[sp.Symbol(f"{prefix}{i+1}{j+1}", nonnegative=True)
             for j in range(m)] for i in range(n)]


def event_equations(n, m, merge):
    """Enumerate ordered parent/child events; never use target A/H formulas."""
    Y = symbols("Y", n)
    H = symbols("H", n, m)
    F, G = symbols("F", n, n), symbols("G", n, n)
    q = [[h / Y[i] for h in H[i]] for i in range(n)]
    p = [sum(row) for row in q]
    dH = [[sp.S.Zero for _ in range(m)] for _ in range(n)]
    migration = [[sp.S.Zero for _ in range(m)] for _ in range(n)]
    score = sp.S.Zero
    for i, j, a, b in itertools.product(range(n), range(n), range(m), range(m)):
        rate = F[i][j] * q[i][a] * q[j][b]
        c, increment = merge(a, b)
        assert 0 <= c < m and increment >= 0
        dH[i][a] -= rate
        dH[j][b] -= rate
        dH[i][c] += rate
        score += rate * increment
    for i, j, b in itertools.product(range(n), range(n), range(m)):
        if i == j:
            continue
        transmission = F[i][j] * (1 - p[i]) * q[j][b]
        movement = G[i][j] * q[j][b]
        dH[i][b] += transmission + movement
        dH[j][b] -= transmission + movement
        migration[i][b] += movement
        migration[j][b] -= movement
    return Y, H, F, G, q, p, dH, migration, score


def moments(n, Y, F, G, p):
    X = symbols("X", n)
    L = symbols("L", n, 2)
    dX, dL = [sp.S.Zero] * n, [[sp.S.Zero] * 2 for _ in range(n)]
    cherries = {(a, b): sp.S.Zero for a in range(2) for b in range(a, 2)}
    for i, j in itertools.product(range(n), repeat=2):
        # Both endpoints disappear from external-lineage counts. Descendant
        # masses add in the parent and therefore only child mass moves.
        dX[i] -= F[i][j] * X[i] / Y[i] * p[j]
        dX[j] -= F[i][j] * p[i] * X[j] / Y[j]
        dX[i] += F[i][j] * (X[i] / Y[i] * p[j] + p[i] * X[j] / Y[j])
        for a in range(2):
            dL[i][a] -= F[i][j] * L[i][a] / Y[i] * p[j]
            dL[j][a] -= F[i][j] * p[i] * L[j][a] / Y[j]
        for a, b in itertools.product(range(2), repeat=2):
            cherries[tuple(sorted((a, b)))] += F[i][j] * L[i][a] / Y[i] * L[j][b] / Y[j]
        if i == j:
            continue
        child_only = F[i][j] * (1 - p[i]) + G[i][j]
        mass_rate = child_only * X[j] / Y[j]
        dX[i] += mass_rate
        dX[j] -= mass_rate
        for a in range(2):
            rate = child_only * L[j][a] / Y[j]
            dL[i][a] += rate
            dL[j][a] -= rate
    equal(sum(dX))
    for i in range(n):
        equal(dX[i], sum((F[i][j] + G[i][j]) * X[j] / Y[j]
                         - (F[j][i] + G[j][i]) * X[i] / Y[i]
                         for j in range(n) if j != i))
        loss = sum((F[i][j] + F[j][i]) * p[j] for j in range(n))
        for a in range(2):
            transfer = sum(
                (F[i][j] * (1-p[i]) + G[i][j]) * L[j][a] / Y[j]
                - (F[j][i] * (1-p[j]) + G[j][i]) * L[i][a] / Y[i]
                for j in range(n) if j != i)
            equal(dL[i][a], -L[i][a] / Y[i] * loss + transfer)
    external_total = sum(sum(row) for row in dL)
    equal(external_total,
          -sum(F[i][j] * (sum(L[i])/Y[i] * p[j] + p[i] * sum(L[j])/Y[j])
               for i, j in itertools.product(range(n), repeat=2)))
    # Cherry creation consumes exactly two external lineages when all extant
    # lineages are external; other coalescences can consume only one or none.
    external_sub = {L[i][a]: sp.S.Zero for i in range(n) for a in range(2)}
    external_sub.update({L[i][0]: p[i] * Y[i] for i in range(n)})
    equal((external_total + 2 * sum(cherries.values())).subs(external_sub))
    return X, L, dX, dL, cherries


def check_events():
    registry = {}
    for n, m in ((1, 2), (2, 2), (3, 3)):
        # A deterministic total merge map is sufficient for row-sum identities;
        # the two-state specialization uses the Camin-Sokal class rule.
        merge = lambda a, b: (min(a, b), int(a != b))
        Y, H, F, G, q, p, dH, migration, score = event_equations(n, m, merge)
        dA = []
        for i in range(n):
            expected = -sum(F[j][i] * p[j] * p[i] for j in range(n))
            expected += sum(
                (F[i][j] * (1-p[i]) + G[i][j]) * p[j]
                - (F[j][i] * (1-p[j]) + G[j][i]) * p[i]
                for j in range(n) if j != i)
            equal(sum(dH[i]), expected)
            dA.append(expected)
        equal(sum(dA), -sum(F[i][j] * p[i] * p[j]
                            for i, j in itertools.product(range(n), repeat=2)))
        equal(sum(sum(row) for row in migration))
        X, L, dX, dL, cherries = moments(n, Y, F, G, p)
        if n == m == 2:
            registry = {f"H{i+1}{a+1}": dH[i][a] for i in range(n) for a in range(m)}
            registry.update({f"A{i+1}": dA[i] for i in range(n)})
            registry.update({f"X{i+1}": dX[i] for i in range(n)})
            registry.update({f"L{i+1}{a+1}": dL[i][a] for i in range(n) for a in range(2)})
            registry.update({f"C{a+1}{b+1}": rate for (a, b), rate in cherries.items()})
            registry["J"] = score
            equal(score, 2*F[0][0]*q[0][0]*q[0][1]
                  + 2*F[1][1]*q[1][0]*q[1][1]
                  + (F[0][1]+F[1][0])*(q[0][0]*q[1][1]+q[0][1]*q[1][0]))
    print("PASS event enumeration: 1/2/3 demes, class row sums, total coalescent loss,")
    print("     migration conservation, descendant mass, external endpoints and cherries.")
    return registry


def baseline_assignments(text, names):
    # Only one-line assignments in the pinned bwd_rhs! block are accepted.
    result = {}
    for name in names:
        matches = re.findall(rf"^\s*{re.escape(name)} = ([^\n]+)$", text, re.M)
        assert len(matches) == 1, (name, matches)
        result[name] = sp.sympify(matches[0].replace("^", "**"), rational=True)
    return result


def check_baseline():
    meta = json.loads((HERE / "source.json").read_text())
    raw = subprocess.check_output(
        ["git", "show", f"{meta['baseline_commit']}:{meta['baseline_path']}"], cwd=ROOT)
    assert hashlib.sha256(raw).hexdigest() == meta["baseline_sha256"]
    text = raw.decode().split("function bwd_rhs!", 1)[1]
    names = ("dA1", "dA2", "dA11", "dA12", "dA21", "dA22",
             "dB11", "dB12", "dB21", "dB22",
             "dF11", "dF12", "dF1s", "dF21", "dF22", "dF2s")
    e = baseline_assignments(text, names)
    f11, f12, f21, f22, p11, p12, p21, p22 = sp.symbols("f11 f12 f21 f22 p11 p12 p21 p22")
    p1, p2 = sp.symbols("p1 p2")
    substitutions = {p1: p11+p12, p2: p21+p22}
    expected = (f12*(p11*p22+p12*p21), f21*(p11*p22+p12*p21))
    for i, pair in enumerate((("dA11", "dA12"), ("dA21", "dA22"))):
        residual = sp.factor((e[pair[0]]+e[pair[1]]-e[f"dA{i+1}"]).subs(substitutions))
        equal(residual, expected[i])
        print(f"BASELINE KNOWN MISMATCH Prior 1 row {i+1}: sum(dH)-dA = {residual}")
    q11, q12, q21, q22 = sp.symbols("q11 q12 q21 q22")
    reverse = {p11: q12, p12: q11, p21: q22, p22: q21}
    for i, pair in enumerate((("dB11", "dB12"), ("dB21", "dB22"))):
        residual = sp.factor((e[pair[0]]+e[pair[1]]-e[f"dA{i+1}"]).subs(
            {p1: q11+q12, p2: q21+q22}))
        equal(residual, expected[i].subs(reverse, simultaneous=True))
        print(f"BASELINE KNOWN MISMATCH Prior 2 row {i+1}: sum(dH)-dA = {residual}")
    r11, r12, r1s, r21, r22, r2s = sp.symbols("r11 r12 r1s r21 r22 r2s")
    for i, row in enumerate((("dF11", "dF12", "dF1s"), ("dF21", "dF22", "dF2s"))):
        equal(sum(e[k] for k in row).subs(
            {sp.Symbol("r1"): r11+r12+r1s, sp.Symbol("r2"): r21+r22+r2s}),
            e[f"dA{i+1}"].subs(
            {p1: r11+r12+r1s, p2: r21+r22+r2s}))
    # Admissible pure-state example: Y1=Y2=2, H11=H22=1,
    # F12=1, all other transmission rates zero.
    example = {p11: sp.Rational(1, 2), p12: 0, p21: 0, p22: sp.Rational(1, 2),
               f11: 0, f12: 1, f21: 0, f22: 0}
    true_row = e["dA1"].subs(substitutions).subs(example)
    legacy_row = (e["dA11"]+e["dA12"]).subs(substitutions).subs(example)
    equal(true_row, sp.Rational(1, 4))
    equal(legacy_row, sp.Rational(1, 2))
    assert true_row != legacy_row
    print("BASELINE COUNTEREXAMPLE Y=(2,2), H11=H22=1, F12=1:")
    print("     Prior 1 deme-1 class sum derivative=1/2, lineage derivative=1/4.")
    print("PASS baseline Fitch row sums; baseline Camin-Sokal mismatches are asserted, not waived.")


def stem(v, C):
    minimum = min(v)
    assert minimum != sp.oo
    return tuple(min(C[i][j]+v[j] for j in range(len(C)))-minimum for i in range(len(C)))


def score_message(tree, C):
    if isinstance(tree, int):
        v = tuple(0 if i == tree else sp.oo for i in range(len(C)))
        return 0, stem(v, C)
    left_score, a = score_message(tree[0], C)
    right_score, b = score_message(tree[1], C)
    if left_score == sp.oo or right_score == sp.oo:
        return sp.oo, None
    t = tuple(x+y for x, y in zip(a, b))
    increment = min(t)
    if increment == sp.oo:
        return sp.oo, None
    return left_score+right_score+increment, stem(t, C)


def brute_score(tree, C):
    """Enumerate all internal-node assignments, independently of message recursion."""
    edges, labels = [], {}
    counter = itertools.count()
    internal = []

    def visit(t):
        node = next(counter)
        if isinstance(t, int):
            labels[node] = t
        else:
            internal.append(node)
            for child in t:
                edges.append((node, visit(child)))
        return node

    visit(tree)
    best = sp.oo
    for assignment in itertools.product(range(len(C)), repeat=len(internal)):
        states = labels | dict(zip(internal, assignment))
        score = sum(C[states[i]][states[j]] for i, j in edges)
        best = min(best, score)
    return best


def shapes(n):
    if n == 1:
        return [None]
    return [(a, b) for k in range(1, n)
            for a in shapes(k) for b in shapes(n-k)]


def label_shape(shape, labels):
    if shape is None:
        return next(labels)
    return tuple(label_shape(t, labels) for t in shape)


def wrong_raw_merge(tree, C):
    """Incorrect: combines unpropagated node vectors and loses edge costs."""
    if isinstance(tree, int):
        return 0, tuple(0 if i == tree else sp.oo for i in range(len(C)))
    sa, a = wrong_raw_merge(tree[0], C)
    sb, b = wrong_raw_merge(tree[1], C)
    t = tuple(x+y for x, y in zip(a, b))
    if min(t) == sp.oo:
        return sp.oo, t
    delta = min(t)
    return sa+sb+delta, tuple(x-delta for x in t)


def check_parsimony():
    matrices = [
        ((0, 1), (1, 0)),
        ((0, 1), (sp.oo, 0)),
        ((0, sp.oo), (1, 0)),
        ((0, 1, 4), (3, 0, 2), (1, 5, 0)),
        ((0, 1, 3), (sp.oo, 0, 1), (sp.oo, sp.oo, 0)),
        ((0, 0, 2), (3, 0, 1), (1, 4, 0)),
        ((0, sp.oo), (sp.oo, 0)),
        ((0, 1, sp.oo), (1, 0, sp.oo), (sp.oo, sp.oo, 0)),
    ]
    cases = 0
    for C in matrices:
        for n in range(1, 5):
            for shape in shapes(n):
                for labels in itertools.product(range(len(C)), repeat=n):
                    tree = label_shape(shape, iter(labels))
                    actual, message = score_message(tree, C)
                    expected = brute_score(tree, C)
                    equal(actual, expected)
                    if actual == sp.oo:
                        assert message is None
                    else:
                        assert min(message) == 0
                    cases += 1
    C = matrices[1]
    equal(score_message((1, 1), C)[0], 0)  # Free human root has no stem charge.
    assert C[0][1]+C[0][1] == 2  # Fixing a reservoir root is a different objective.
    equal(score_message((0, 1), C)[0], 1)
    assert wrong_raw_merge((0, 1), C)[0] == sp.oo
    # Retaining only argmin states loses directed cost information.
    asymmetric = matrices[3]
    a = stem((0, sp.oo, sp.oo), asymmetric)
    b = stem((sp.oo, 0, sp.oo), asymmetric)
    equal(min(x+y for x, y in zip(a, b)), 1)
    assert min(a)+min(b) == 0
    assert score_message((0, 1), matrices[6]) == (sp.oo, None)
    print(f"PASS Sankoff: {cases} labeled plane binary trees with 1–4 tips, eight cost matrices;")
    print("     exact ancestral assignments, free root, directed/irreversible/disconnected costs.")
    print("PASS counterexamples: no edge propagation and independent-minimum merging are wrong.")


def check_generated(path, registry):
    assert path.is_file(), f"Required generated CSV missing: {path}"
    local = {str(s): s for expr in registry.values() for s in expr.free_symbols}
    for row in symbols("G", 2, 2):
        for symbol in row:
            local.setdefault(str(symbol), symbol)
    rows = 0
    claims = set()
    independent = set()
    corrected = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == ["claim", "expression", "expected"], reader.fieldnames
        for row in reader:
            assert row["claim"] and row["expression"] and row["expected"]
            assert row["claim"] not in claims, f"Duplicate generated claim: {row['claim']}"
            claims.add(row["claim"])
            # Julia exports parameter names such as beta, which otherwise
            # resolve to SymPy functions instead of indeterminates.
            for value in (row["expression"], row["expected"]):
                for token in tokenize.generate_tokens(io.StringIO(value).readline):
                    if token.type == tokenize.NAME:
                        local.setdefault(token.string, sp.Symbol(token.string))
            expr = sp.sympify(row["expression"], locals=local, rational=True)
            expected = sp.sympify(row["expected"], locals=local, rational=True)
            equal(expr, expected)
            if row["claim"] in registry:
                equal(expr, registry[row["claim"]])
                independent.add(row["claim"])
            if row["claim"] in {f"legacy_{k}" for k in range(15, 23)}:
                index = int(row["claim"].split("_")[1])
                reverse = index >= 19
                offset = 19 if reverse else 15
                slot = index - offset
                deme, state = divmod(slot, 2)
                target_state = 1-state if reverse else state
                substitutions = {}
                for i, j in itertools.product(range(2), repeat=2):
                    actual_state = 1-j if reverse else j
                    substitutions[local[f"H{i+1}{j+1}"]] = local[f"z{offset+2*i+actual_state}"]
                    substitutions[local[f"G{i+1}{j+1}"]] = 0
                beta, epsilon, phi, psi = (local[name] for name in ("beta", "epsilon", "phi", "psi"))
                X1, X2, Y1, Y2 = (local[name] for name in ("X1", "X2", "Y1", "Y2"))
                substitutions.update({
                    local["F11"]: X1*Y1*beta,
                    local["F12"]: X2*Y1*beta*epsilon,
                    local["F21"]: X1*Y2*beta*psi,
                    local["F22"]: X2*Y2*beta*phi,
                })
                derived = registry[f"H{deme+1}{target_state+1}"].subs(substitutions, simultaneous=True)
                constraints = {
                    local["z1"]: local[f"z{offset}"]+local[f"z{offset+1}"],
                    local["z2"]: local[f"z{offset+2}"]+local[f"z{offset+3}"],
                }
                equal(expr.subs(constraints, simultaneous=True), derived)
                corrected.add(index)
            rows += 1
    assert rows > 0, "Generated CSV must not be empty"
    assert independent == set(registry), f"Missing independent claims: {set(registry)-independent}"
    assert {f"legacy_{k}" for k in range(1, 38)} <= claims, "Incomplete 37-state legacy export"
    assert corrected == set(range(15, 23)), f"Incomplete directional derivative export: {corrected}"
    print(f"PASS generated CSV: {rows} exact expression/expected identities;")
    print(f"     {len(independent)} distinct independently derived claims: {', '.join(sorted(independent)) or 'none'}.")
    print("PASS 37 legacy cross-generator identities; corrected directional CS derivatives 15–22")
    print("     additionally match independent events under class-row constraints and reverse relabeling.")
    return independent


def main():
    if not __debug__:
        raise RuntimeError("Run without -O: verification assertions must remain enabled")
    assert sys.prefix != sys.base_prefix, "Use verification/.venv/bin/python"
    assert Path(sys.prefix).resolve() == (HERE / ".venv").resolve(), "Wrong virtualenv"
    assert sp.__version__ == "1.14.0", "Install pinned verification/requirements.txt"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated", type=Path)
    args = parser.parse_args()
    print("SCOPE: factorized mean-field algebra and finite-tree dynamic programming only.")
    print("NOT a proof of stochastic expectation closure, numerical integration, or current Julia runtime behavior.")
    registry = check_events()
    check_baseline()
    check_parsimony()
    if args.generated:
        check_generated(args.generated.resolve(), registry)
    else:
        print("NOT CHECKED: Julia Symbolics export (supply --generated PATH to require it).")
    print(f"PASS: {CHECKS} exact symbolic identities; all required checks completed.")


if __name__ == "__main__":
    main()
