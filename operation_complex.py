"""
Tool for examining simplicial complexes corresponding to structures built with algebraic_structure

All of these methods are applicable only to unary or binary operations.

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

from __future__ import division
import sage.all
from sage.plot.polygon import polygon
from math import cos, sin, pi
from operator import add
from itertools import chain
from sage.plot.colors import colors
from numpy import cross, subtract
from sage.graphs.digraph import DiGraph
from algebraic_structure import FiniteRingStructure

def polar_location(radius, numerator, denominator, height=0, dimension=2):
    """
    Plot a point at distance `radius` from the origin at the appropriate angle for the "flower" construction.

    Args:
        radius (float): The distance from the origin.
        numerator (int): The number of steps around a cycle.
        denominator (int): The number of steps in one cycle.
        height (int): The height off the plane in 3-space.
        dimension (int): The number of entries in the resulting tuple.

    Returns:
        tuple: A tuple specifying a point in space.
    """

    tau = 2*pi
    frac = numerator/denominator
    x = radius * cos(frac * tau + (pi/denominator) * (radius-1))
    y = radius * sin(frac * tau + (pi/denominator) * (radius-1))
    z = height
    if dimension == 2:
        return (x,y)
    if dimension == 3:
        return (x,y,z)

def barycentric_subdivision(p1, p2, p3):
    """
    Find the vertices of the triangles in the barycentric subdivision of a given triangle.

    Args:
        p1, p2, p3 (tuple): Vertices of a triangle given as tuples.

    Returns:
        tuple of tuple of tuple: The vertices of each triangle in the subdivision.
    """

    avg = lambda p,q: tuple(map(lambda x: x/2,map(add,p,q)))
    return ((p1,avg(p1,p2),avg(p1,p3)),(p2,avg(p1,p2),avg(p2,p3)),(p3,avg(p1,p3),avg(p2,p3)),(avg(p1,p2),avg(p1,p3),avg(p2,p3)))

def operation_digraph(structure, operation, render=False, file_name='operation_digraph.asy', labels=True, size=100):
    """
    A digraph corresponding to a unary operation.

    Note that one might need to play with the `order` setting to get a decent image.

    Args:
        structure (AlgebraicStructure): The algebraic structure used in the construction of the complex.
        operation (function): The operation used in the construction of the complex.
        render (bool): Whether the digraph should be given as asymptote vector graphics source code.
        file_name (str): The name of the file generated is `render` is True.
        labels (bool): Whether to label elements in the asymptote source.
        order (int): The order of the resulting image.

    Returns:
        Digraph: The corresponding operation digraph.
    """

    m = structure.action_matrix(operation)
    if render == True:
        order = structure.order
        file = open(file_name, 'w')
        file.write('size({});'.format(size) + '\n' +
                  'dotfactor=10;')
        file.write('for(int i=0; i<{}; ++i){{'.format(order) + '\n' +
                  'dot(dir(90-360/{}*i));}}'.format(order) + '\n')
        elem_lis = structure.canonical_order
        for i in range(order):
                if labels:
                    file.write('label(\"${}$\",1.2*dir(90-360/{}*{}));'.format(structure.element_names[elem_lis[i]], order, i) + '\n')
                if list(m[i]).index(1) == i:
                        file.write('draw(dir(90-360/{0}*{1})..1.35*dir(90-360/{0}*{1})..cycle,MidArcArrow);'.format(order, i) + '\n')
                else:
                    file.write('draw(dir(90-360/{0}*{1})..dir(90-360/{0}*{2}),MidArrow);'.format(order, i, list(m[i]).index(1)) + '\n')
    return DiGraph(m)

class OperationComplex:
    """
    A simplicial complex corresponding to a binary operation.

    Attributes:
        structure (AlgebraicStructure): The algebraic structure used in the construction of the complex.
        operation (function): The operation used in the construction of the complex.
        operation_name (str): The name of the operation used in the construction of the complex.
        initial_colors (tuple): The first colors to use when displaying the complex.
        color_dic (dict): Dictionary of colors to be used when displaying the complex.
    """

    def __init__(self, structure, operation_name, initial_colors=(), color_dic=None):
        """
        Args:
            structure (AlgebraicStructure): The algebraic structure used in the construction of the complex.
            operation_name (str): The name of the operation used in the construction of the complex.
            initial_colors (tuple): The first colors to use when displaying the complex if.
            color_dic (dict): Dictionary of colors to be used when displaying the complex if.
        """

        self.structure = structure
        self.operation = structure.operations[operation_name]
        self.operation_name = operation_name
        self.initial_colors = initial_colors
        self.color_dic = color_dic

    def __repr__(self):
        return 'simplicial complex for {} under {}'.format(self.structure,self.operation_name)

    def grow_flower(self, initial_elements, sides, max_layer):
        """
        Create a flower complex. Currently this is only for 2-to-1 operations.

        Args:
            initial_elements (tuple): The first elements to be operated on.
            sides (int): The number of sides of the original polygon.
            max_layer (int): The number of steps out from the origin to construct.

        Returns:
            tuple of tuple: Each entry contains the vertex elements of `structure` listed in order. See examples.
        """

        f = self.operation
        lis = [[initial_elements[0]], [initial_elements[1]]]
        for petal in range(sides-1):
            lis[1].append(f(lis[0][0], lis[1][petal]))
        for layer in range(1,max_layer):
            lis.append([])
            for i in range(len(lis[layer])):
                lis[layer+1].append(f(lis[layer][i], lis[layer][(i+1)%sides]))
        return tuple(map(tuple,lis))

    def grow_spiral(self, initial_elements, sides, max_height):
        """
        Create a (left) spiral complex. Currently this is only for 2-to-1 operations.

        Args:
            initial_elements (tuple): The first elements to be operated on.
            sides (int): The number of sides of the original polygon.
            max_height (int): The number of steps up from the origin to construct.

        Returns:
            tuple of tuple: Each entry contains the vertex elements of `structure` listed in order. See examples.
        """

        f = self.operation
        lis = [[initial_elements[0]], [initial_elements[1]]]
        for height in range(max_height):
            for petal in range(sides-1): lis[height+1].append(f(lis[0][0], lis[height+1][petal]))
            lis.append([f(lis[0][0], lis[height+1][sides-1])])
        lis.pop(-1)
        return tuple(map(tuple,lis))

    def color_dictionary(self, values, initial_colors=()):
        """
        Create a dictionary mapping given values to colors.

        Args:
            values (tuple of tuple): A tuple such as that output by the `grow_flower` method.
            initial_colors (tuple): Optional tuple of colors to be used first.
            # TODO: Warn when a duplicate color is used.

        Returns:
            color_dic (dict): Dictionary mapping each value to a color.
        """

        val_tup = tuple(set(chain(*values)))
        col_tup = initial_colors + self.initial_colors + tuple(set(colors)-set(('white','black','grey')+tuple('light{}'.format(col) for col in colors)))
        color_dic = {}
        for val in val_tup:
            color_dic[val] = col_tup[val_tup.index(val)]
        return color_dic

    def render_flower(self, initial_elements, sides, max_layer, file_name, dimension=2, use_height=False, initial_colors=(), color_dic=None):
        """
        Create an image file depicting a flower complex.

        Note that the identification between colors and elements is not preserved across images with different
        `sides` and `max_layer` values unless a color dictionary is provided. Don't be misled by this.

        Also note that `max_layer` cannot be set much higher than `sides` without triangles overlapping. This can
        be dealt with in 3-space by using the `use_height` setting.

        Args:
            initial_elements (tuple): The first elements to be operated on.
            sides (int): The number of sides of the original polygon.
            max_layer (int): The number of steps out from the origin to construct.
            file_name (str): The name of the file generated.
            dim (int): The number of coordinate axes in the plot space.
            use_height (bool): If `dimension` is 3 then make use of the extra room in space.
            intial_colors (tuple): Colors used to display the complex. These are used before the colors in `self.initial_colors`.
            color_dic (dict): Dictionary of colors to be used when displaying the complex. This is used before the dictionary in `self.color_dic`.
        """

        values = self.grow_flower(initial_elements,sides,max_layer)
        while color_dic is None:
            color_dic = self.color_dic
            color_dic = self.color_dictionary(values,initial_colors)
        g = 0
        loc = lambda a,b,c=0: polar_location(a, b, sides, height=c, dimension=dimension)
        # initial points
        for n in range(sides-1):
            p1 = (0,)*dimension
            p2 = loc(1,n)
            p3 = loc(1,n+1)
            if dimension == 3 and use_height: p2 = loc(1,n,1); p3 = loc(1,n+1,1)
            c1 = color_dic[values[0][0]] # left input
            c2 = color_dic[values[1][n]] # right input
            c3 = color_dic[values[1][n+1]] # output
            sub_div = barycentric_subdivision(p1,p2,p3)
            cols = (c1,c2,c3,c3)
            for i in range(4):
                g+=polygon(sub_div[i],color=cols[i])
        # subsequent points
        for r in range(1,max_layer-1):
            for n in range(sides):
                p1 = loc(r,n)
                p2 = loc(r,n+1)
                p3 = loc(r+1,n)
                if dimension == 3 and use_height: p1 = loc(r,n,r); p2 = loc(r,n+1,r); p3 = loc(r+1,n,r+1)
                c1 = color_dic[values[r][n]] # left input
                c2 = color_dic[values[r][(n+1)%sides]] # right input
                c3 = color_dic[values[r+1][n]] # output
                sub_div = barycentric_subdivision(p1,p2,p3)
                cols = (c1,c2,c3,c3)
                for i in range(4):
                    g+=polygon(sub_div[i],color=cols[i])
        if dimension == 2: g.save(file_name+'.svg',axes=False)
        if dimension == 3: g.save(file_name+'.x3d')

    def render_spiral(self, initial_elements, sides, max_height, file_name, initial_colors=(), color_dic=None):
        """
        Create an image file depicting a flower complex.

        Note that the identification between colors and elements is not preserved across images with different
        `sides` and `max_layer` values unless a color dictionary is provided. Don't be misled by this.

        Args:
            initial_elements (tuple): The first elements to be operated on.
            sides (int): The number of sides of the original polygon.
            max_beight (int): The number of steps up from the origin to construct.
            file_name (str): The name of the file generated.
            intial_colors (tuple): Colors used to display the complex. These are used before the colors in `self.initial_colors`.
            color_dic (dict): Dictionary of colors to be used when displaying the complex. This is used before the dictionary in `self.color_dic`.
        """

        values = self.grow_spiral(initial_elements,sides,max_height)
        while color_dic is None:
            color_dic = self.color_dic
            color_dic = self.color_dictionary(values,initial_colors)
        g=0
        loc = lambda b,c,a=1: polar_location(a, b, sides, height=c, dimension=3)
        for h in range(max_height):
            for n in range(sides-1):
                p1 = loc(0,h,0)
                p2 = loc(n,h)
                p3 = loc(n+1,h)
                c1 = color_dic[values[0][0]] # left input
                c2 = color_dic[values[h+1][n]] # right input
                c3 = color_dic[values[h+1][n+1]] # output
                sub_div = barycentric_subdivision(p1,p2,p3)
                cols = (c1,c2,c3,c3)
                for i in range(4):
                    g+=polygon(sub_div[i],color=cols[i])
            if h != max_height-1:
                p1 = loc(0,h,0)
                p2 = loc(sides-1,h)
                p3 = loc(0,h+1)
                c1 = color_dic[values[0][0]] # left input
                c2 = color_dic[values[h+1][sides-1]] # right input
                c3 = color_dic[values[h+2][0]] # output
                sub_div = barycentric_subdivision(p1,p2,p3)
                cols = (c1,c2,c3,c3)
                for i in range(4):
                    g+=polygon(sub_div[i],color=cols[i])
        g.save(file_name+'.x3d',axes=False)

    def render_complex(self, file_name, cycles=1, degenerate_faces=True, initial_colors=(), color_dic=None):
        """
        Create 3d graphics depicting a binary operation complex.

        Vertices are placed on a helix which winds around the surface of the sphere with center (0,0,1/2) and radius 1/2.
        That is, the first vertex is placed at (0,0,0) and the last is placed at (0,0,1).

        Args:
            file_name (str): The name of the file generated.
            cycles (int): The number of times the determining helix winds around the line between the origin and (0,0,1).
            intial_colors (tuple): Colors used to display the complex. These are used before the colors in `self.initial_colors`.
            color_dic (dict): Dictionary of colors to be used when displaying the complex. This is used before the dictionary in `self.color_dic`.
        """

        elems = tuple(self.structure.elements)
        f = self.operation
        while color_dic is None:
            color_dic = self.color_dic
            color_dic = self.color_dictionary((elems,),initial_colors)
        g=0
        loc = lambda a,b,c: (sin(c*pi)*cos(a*2*pi/b),sin(c*pi)*sin(a*2*pi/b),c)
        locations = {}
        for h in range(len(elems)):
            locations[elems[h]] = loc(cycles*h,len(elems),h/len(elems))
        for x in elems:
            for y in elems:
                z = self.operation(x,y)
                if degenerate_faces:
                    if x == y:
                        # x=y=z case
                        if x == z:
                            if locations[x] == (0,0,0):
                                p1, p2, p3, p4 = (0,0,0), loc(0,3,-1/2), loc(1,3,-1/2), loc(2,3,-1/2)
                                for entry in ((p1,p2,p3),(p1,p2,p4),(p1,p3,p4),(p2,p3,p4)):
                                    g+=polygon(entry,colors=color_dic[x])
                            if locations[x] == (0,0,1):
                                p1, p2, p3, p4 = (0,0,1), loc(0,3,3/2), loc(1,3,3/2), loc(2,3,3/2)
                                for entry in ((p1,p2,p3),(p1,p2,p4),(p1,p3,p4),(p2,p3,p4)):
                                    g+=polygon(entry,colors=color_dic[x])
                            if locations[x] != (0,0,0) and locations[x] != (0,0,1):
                                p1 = locations[x]
                                p2 = tuple(map(add,map(lambda u: 14*u/10,locations[x]),map(lambda u: 6*u/10,locations[elems[elems.index(x)+1]])))
                                p3 = tuple(map(add,map(lambda u: 14*u/10,locations[x]),map(lambda u: 6*u/10,locations[elems[elems.index(x)-1]])))
                                p4 = tuple(map(add,map(lambda u: 14*u/10,locations[x]),map(lambda u: 6*u/10,map(add,locations[x],(0,0,1)))))
                                for entry in ((p1,p2,p3),(p1,p2,p4),(p1,p3,p4),(p2,p3,p4)):
                                    g+=polygon(entry,colors=color_dic[x])
                    # x=y or x=z or y=z case
                    if (x==y and x!=z) or (x==z and x!=y) or (y==z and x!=y):
                        if x==y: s,t = x,z
                        if x==z: s,t = x,y
                        if y==z: s,t = x,z
                        p1 = locations[s]
                        p2 = locations[t]
                        a = map(lambda u: u/2,map(add,p1,p2))
                        b = map(subtract,p1,p2)
                        p3 = map(add,a,cross(b,(0,0,1)))
                        p4 = map(add,p3,(0,0,1/len(elems)))
                        polys = ((p1,p2,p3),(p1,p2,p4),(p1,p3,p4),(p2,p3,p4))
                        cols = (color_dic[s],color_dic[t],color_dic[s],color_dic[t])
                        for i in range(4):
                            g+=polygon(polys[i],colors=cols[i])
                # x,y,z are distinct case
                if x!=y and x!=z and y!=z:
                    if f(x,y) == f(y,x):
                        p1, p2, p3 = locations[x], locations[y], locations[z]
                        p4 = cross(map(subtract,p1,p2),map(subtract,p1,p3))
                        c1, c2, c3 = color_dic[x], color_dic[y], color_dic[z]
                        cols = (c1,c2,c3)
                        faces = ((p1,p2,p4),(p1,p3,p4),(p2,p3,p4))
                        for i in range(3):
                                    g+=polygon(faces[i],color=cols[i])
                    else:
                        p1, p2, p3 = locations[x], locations[y], locations[z]
                        c1, c2, c3 = color_dic[x], color_dic[y], color_dic[z]
                        sub_div = barycentric_subdivision(p1,p2,p3)
                        cols = (c1,c2,c3,c3)
                        for i in range(4):
                            g+=polygon(sub_div[i],color=cols[i])
        g.save(file_name+'.x3d',axes=False)

def examine_structure(structure, operation_name, min_size=0, max_size=0):
    """
    Generate a collection of flower complexes for a specified structure.

    Args:
        structure (AlgebraicStructure): The algebraic structure used in the construction of the complex.
        operation_name (str): The name of the operation used in the construction of the complex.
        min_size (int): The minimum number of petals to place on a generated flower.
        max_size (int): The maximum number of petals to place on a generated flower.
    """

    n = structure.order
    if min_size == 0:
        min_size = n
    if max_size == 0:
        max_size = n
    elems = structure.canonical_order
    f = OperationComplex(structure, operation_name)
    for size in range(min_size,max_size+1):
        for i in range(n):
            for j in range(n):
                x = elems[i]
                y = elems[j]
                if j>=i or structure.operations[operation_name](x,y) != structure.operations[operation_name](y,x):
                    f.render_flower((x,y),size,size,'{}({},{},{}) {} petals'.format(structure,operation_name,str(x),str(y),size))