"""
Example use of the operation_complex tool

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from algebraic_structure import FiniteFieldStructure, GroupStructure, FiniteRingStructure
from operation_complex import operation_digraph, OperationComplex
from sage.groups.perm_gps.permgroup_named import DihedralGroup, SymmetricGroup

print('We use FiniteFieldStructure to create an AlgebraicStructure object for F_5.')
a = FiniteFieldStructure(5)
print(a)
print('')

print('We draw an operation digraph for adding 2 in F_5.')
operation_digraph(a,lambda x: x+2, render=True)

print('We the create a complex for the addition operation.')
elems = list(a.elements)
f = OperationComplex(a,'+')
print(f)
print('')
 
print('We print a flower starting with 1+1 with six sides up to its fourth layer.')
print(f.grow_flower((elems[1],elems[1]),6,4))
print('')
 
print('We print a spiral starting with 1+1 with six sides up to its fourth layer.')
print(f.grow_spiral((elems[1],elems[1]),6,4))
print('')
 
print('Draw the flower with 6 sides out to layer 4.')
f.render_flower((elems[1],elems[1]),6,4,'flower')
print('')
 
print('The default colors aren\'t the best, so we override them.')
f.render_flower((elems[1],elems[1]),6,4,'flower_better_colors',initial_colors=('blue','red','green','purple','orange'))
print('')
 
print('Let\'s examine a batch of \'flower\' covering space images for the dihedral group on 4 vertices.')
b = GroupStructure(DihedralGroup(4))
print(b)
print('')
g = OperationComplex(b,'*')
elems = list(b.elements)
x = elems[2]
y = elems[5]
print('We will start our complexes with the product of {} and {}.'.format(x,y))
print('We let the number of \'petals\' on our complex vary from 4 to 29.')
print('We also adjust the number of \'layers\' we draw in order to avoid overlapping triangles in the plane.')
cols = ('blue','red','green','purple','orange')
for i in range(4,30):
    g.render_flower((x,y),i,i,'dihedral{}'.format(i),initial_colors=cols)
print('')
  
print('Another way to avoid triangles overlapping is to separate them in 3-space.')
print('This is illustrated in the following example, which can be imported into a 3d graphics program such as Blender.')
g.render_flower((x,y),11,23,'3d_example',dimension=3,use_height=True,initial_colors=cols)
 
print('We can also display helical complexes.')
print('Combining these with the \'flower\' constructions can produce interesting 3d objects.')
f.render_flower((elems[1],elems[1]),6,4,'1+1_flower',dimension=3,initial_colors=cols)
f.render_spiral((elems[1],elems[1]),6,4,'1+1_spiral',initial_colors=cols)
 
print('The actual simplicial complex in question can be drawn, as well.')
print('In this case we examine F_3 under addition.')
c = FiniteFieldStructure(3)
h = OperationComplex(c,'+')
h.render_complex('operation_complex',cycles=1)
print('')
 
print('We can remove all of the degenerate parts of the complex for a clearner picture.')
h.render_complex('operation_complex_no_degens',cycles=1,degenerate_faces=False)