"""
Algebraic structure FunctionStructure example

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

from algebraic_structure import FunctionStructure
from operation_complex import operation_digraph

print('Create a structure for the functions on a set of order 3.')
a = FunctionStructure(3)
print(a)
print('')

print('Use `f` as a shorthand for function composition.')
f = a.operations['*']
print('Multiply two elements.')
s = a.canonical_order[4]
t = a.canonical_order[15]
u = f(s,t)
u.set_immutable()
for i in range(3): print(s[i],t[i],u[i]) # This is just so we can view the matrices side-by-side.
print(tuple(a.element_names[p] for p in (s,t,u)))
print('')

print('Examine the left action matrix for an element.')
m = a.action_matrix(lambda x: s*x)
for r in m:
    print(r)
print('')

print('View the operation digraph for left multiplication by `s`.')
operation_digraph(a, lambda x: s*x, render=True)