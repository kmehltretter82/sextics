"""Verify the cycle h^2 obstruction without using symmetry or matrix generation.

The adjoint functional is tested on every possible numerator coefficient
with normal degree below 12. Only Python's standard library is required.
"""
import argparse
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path

VARIABLES = ('X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2')
ZERO = (0,) * 6

def plus(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, 0) + coefficient
    return {e: c for e, c in result.items() if c}

def scale(poly, coefficient):
    return {e: c * coefficient for e, c in poly.items() if c * coefficient}

def times(left, right):
    result = {}
    for a, c in left.items():
        for b, d in right.items():
            exponent = tuple(x + y for x, y in zip(a, b))
            result[exponent] = result.get(exponent, 0) + c * d
    return {e: c for e, c in result.items() if c}

def power(poly, n):
    result = {ZERO: Fraction(1)}
    for _ in range(n):
        result = times(result, poly)
    return result

def derivative(poly, index):
    result = {}
    for exponent, coefficient in poly.items():
        if exponent[index]:
            e = list(exponent)
            e[index] -= 1
            result[tuple(e)] = coefficient * exponent[index]
    return result

def variable_times(poly, index):
    result = {}
    for exponent, coefficient in poly.items():
        e = list(exponent)
        e[index] += 1
        result[tuple(e)] = coefficient
    return result

def parse(node):
    if isinstance(node, ast.Expression):
        return parse(node.body)
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return {} if node.value == 0 else {ZERO: Fraction(node.value)}
    if isinstance(node, ast.Name) and node.id in VARIABLES:
        return {tuple(int(i == VARIABLES.index(node.id)) for i in range(6)): Fraction(1)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return scale(parse(node.operand), -1)
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Pow):
            assert isinstance(node.right, ast.Constant) and type(node.right.value) is int and node.right.value >= 0
            return power(parse(node.left), node.right.value)
        left, right = parse(node.left), parse(node.right)
        if isinstance(node.op, ast.Add):
            return plus(left, right)
        if isinstance(node.op, ast.Sub):
            return plus(left, scale(right, -1))
        if isinstance(node.op, ast.Mult):
            return times(left, right)
        if isinstance(node.op, ast.Div):
            assert set(right) == {ZERO} and right[ZERO] != 0
            return scale(left, 1 / right[ZERO])
    raise ValueError('Unsupported protected expression')

parser = argparse.ArgumentParser()
parser.add_argument('--instance', type=Path, required=True)
parser.add_argument('--certificate', type=Path, required=True)
parser.add_argument('--output', type=Path)
args = parser.parse_args()
record = args.instance.read_bytes()
fields = dict(line.split(':', 1) for line in record.decode().splitlines() if ':' in line and not line.startswith('#'))
assert fields['id'].strip() == 'cycle-q2' and fields['q'].strip() == '2'
F = parse(ast.parse(fields['F'].strip(), mode='eval'))
A = parse(ast.parse(fields['A'].strip(), mode='eval'))
assert all(sum(e) == 6 for e in F) and all(sum(e) == 12 for e in A)
assert F[(0, 0, 0, 0, 0, 6)] == 1
h = plus(variable_times(derivative(F, 0), 0), scale(variable_times(derivative(F, 3), 3), -1))
assert h == {(6, 0, 0, 0, 0, 0): 6, (0, 0, 0, 6, 0, 0): -6, (3, 1, 0, 1, 1, 0): 1}
assert all(e[5] == 0 for e in h)
certificate_bytes = args.certificate.read_bytes()
saved = json.loads(certificate_bytes)
assert saved['m'] == 2
dual = saved['exact_dual_certificate']
functional = {tuple(e): int(c) for e, c in dual['functional']}
assert len(functional) == len(dual['functional'])
assert all(len(e) == 6 and min(e) >= 0 and sum(e) == 30 for e in functional)

# 2 times the strict equation, after cancelling h:
# sum_i [ (4 h F_i + 4 F h_i) N_i - 2 F h partial_i N_i ] = 2 A h^3.
# We compute its full coefficient adjoint directly. No list of numerator
# monomials, diagonal character, group order, or claimed matrix rank is used.
derivative_multiplier = scale(times(F, h), 2)
candidate_counts = []
for index in range(6):
    multiplication = plus(scale(times(h, derivative(F, index)), 4), scale(times(F, derivative(h, index)), 4))
    adjoint = {}
    candidates = set()
    for row, weight in functional.items():
        for exponent, coefficient in multiplication.items():
            k = tuple(r - e for r, e in zip(row, exponent))
            if min(k) < 0 or k[5] >= 12:
                continue
            assert sum(k) == 19
            candidates.add(k)
            adjoint[k] = adjoint.get(k, 0) + weight * coefficient
        for exponent, coefficient in derivative_multiplier.items():
            k = tuple(r - e + int(j == index) for j, (r, e) in enumerate(zip(row, exponent)))
            if min(k) < 0 or k[5] >= 12 or k[index] == 0:
                continue
            assert sum(k) == 19
            candidates.add(k)
            adjoint[k] = adjoint.get(k, 0) - k[index] * weight * coefficient
    assert all(coefficient == 0 for coefficient in adjoint.values())
    candidate_counts.append(len(candidates))
rhs = scale(times(A, power(h, 3)), 2)
pairing = sum(functional.get(e, 0) * coefficient for e, coefficient in rhs.items())
assert pairing.denominator == 1 and pairing != 0
assert int(pairing) == int(dual['rhs_pairing'])
result = {
    'status': 'PASS',
    'instance_sha256': hashlib.sha256(record).hexdigest(),
    'certificate_sha256': hashlib.sha256(certificate_bytes).hexdigest(),
    'functional_terms': len(functional),
    'adjoint_nontrivial_candidates_per_component': candidate_counts,
    'normal_degree_bound': 12,
    'numerator_degree': 19,
    'rhs_pairing': str(int(pairing)),
    'interpretation': 'The functional vanishes on every degree-19 numerator monomial with Y2 degree below 12, without symmetry assumptions.'
}
if args.output:
    args.output.write_text(json.dumps(result, indent=2) + '\n')
print('PASS: all six full coefficient adjoints vanish; source pairing is nonzero.')
print('No symmetry projection or matrix-rank assertion was used.')
print('Nontrivial candidate counts:', candidate_counts)
