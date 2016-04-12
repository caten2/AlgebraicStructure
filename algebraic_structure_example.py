"""
Algebraic structure example

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from algebraic_structure import dict_to_lambda, AlgebraicStructure, FiniteFieldStructure, GroupStructure
from sage.groups.perm_gps.permgroup_named import AlternatingGroup

# Define a structure manually by calling AlgebraicStructure directly.
print('Create an algebraic structure `a` on the set {0,1} with one binary operation, called \'mult\'.')
dic = {(0,0): 0, (0,1): 0, (1,0): 0, (1,1): 1}
a = AlgebraicStructure({0,1},{'mult': dict_to_lambda(dic)})
print('Print the elements in `a`.')
print(list(a.elements))
print('')

# Use `f` as shorthand for the `mult` operation.
f = a.operations['mult']
print('Perform \'mult\' on the pair (0,1) and give the result.')
print(f(0,1))
print('')

# Define a structure using FiniteFieldStructure.
print('Create an algebraic structure `b` for the finite field with 2 elements.') 
b = FiniteFieldStructure(2)
print(b)
print('')

print('Print a list of the elements in `b`.')
print(list(b.elements))
print('')

# Use `g` as shorthand for the multiplication operation.
g = b.operations['*']
print('Perform multiplication on the pair (0,1) and give the result.')
print(g(0,1))
print('')

print('Check the canonical ordering on `b`.')
print(b.canonical_order)
print('')

print('This ordering is used to create the action matrix for addition by 2.')
# Write "add 2" as a lambda function.
f = lambda x: b.operations['+'](x,b.canonical_order[1])
print(b.action_matrix(f))
print('')

print('An example using the ``GroupStructure`` device with the alternating group on 5 letters follows.')
c = GroupStructure(AlternatingGroup(5))
print(c)
print('')

# Use `h` as shorthand for the group operation.
h = c.operations['*']
x = list(c.elements)[5]
y = list(c.elements)[26]
print('We multiply {} and {}.'.format(x,y))
print(h(x,y))