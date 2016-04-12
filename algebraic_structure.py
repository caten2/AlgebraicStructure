"""
Algebraic structure

Charlotte Aten (caten2@u.rochester.edu) 2016
"""

import sage.all
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.matrix.constructor import matrix
from itertools import product
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sympy.matrices.immutable import ImmutableSparseMatrix

def dict_to_lambda(dic):
    """
    Turn a dictionary into a lambda function.

    Args:
        dic (dict of tuple: object): An operation formatted as a dictionary. See examples.

    Returns:
        function: A function which takes components of a key tuple and returns the corresponding entry in `dic`.
    """

    return lambda *keys: dic[keys]

class AlgebraicStructure:
    """
    An algebraic structure with elements and operations.

    Attributes:
        elements (frozenset): The elements in the structure.
        canonical_order (tuple): The standard ordering used on `elements`. #TODO: Allow user to give order specification.
        operations (dict of str: func): The operations in the structure.
        order (int): The number of elements in the structure.
    """

    def __init__(self, elements, operations, element_names={}):
        """
        Args:
            elements (set): The elements in the structure.
            operations (dict of str: func): The operations in the structure.
            element_names (dic of Object: str): The human-readable names of the elements in the structure.
        """

        self.elements = frozenset(elements)
        self.canonical_order = tuple(self.elements)
        self.operations = operations
        self.order = len(self.elements)
        if element_names == {}:
            self.element_names = {elem: str(elem) for elem in self.elements}
        else:
            self.element_names = element_names

    def __repr__(self):
        return 'an algebraic structure ({})'.format(id(self))

    def action_matrix(self, operation):
        """
        Create a matrix encoding the action of `operation` on the elements of `self`.

        Args:
            operation (function): A unary operation defined on `self.elements`.

        Returns:
            sage.matrix.matrix_integer_sparse.Matrix_integer_sparse: The matrix obtained by applying `elem` to `self.elements`.
        """

        m = matrix(self.order,sparse=True)
        tup = self.canonical_order
        f = operation
        for i in range(self.order):
            m[i,tup.index(f(tup[i]))] = 1
        return m

class FiniteFieldStructure(AlgebraicStructure):
    """
    An AlgebraicStructure for the finite field of order `n`.
    """

    def __init__(self, n):
        """
        Args:
            n (int): The number of elements in the finite field.
        """

        f = FiniteField(n)
        operations = {'+': lambda x,y: x+y, '*': lambda x,y: x*y}
        AlgebraicStructure.__init__(self,set(f),operations)

    def __repr__(self):
        return 'finite field with {} elements'.format(self.order)

class FiniteRingStructure(AlgebraicStructure):
    """
    An AlgebraicStructure for the ring of integers modulo `n`.
    """

    def __init__(self, n):
        """
        Args:
            n (int): The number of elements in the finite ring.
        """

        f = Integers(n)
        operations = {'+': lambda x,y: x+y, '*': lambda x,y: x*y}
        AlgebraicStructure.__init__(self,set(f),operations)

    def __repr__(self):
        return 'ring of integers modulo {}'.format(self.order)

class GroupStructure(AlgebraicStructure):
    """
    An AlgebraicStructure for a given group.

    Attributes:
        group (SageMath permutation group): The group in question.
    """

    def __init__(self, group):
        """
        Args:
            group (SageMath permutation group): The group in question.
        """

        operation = {'*': lambda x,y: x*y}
        self.group = group
        AlgebraicStructure.__init__(self,set(group),operation)

    def __repr__(self):
        return self.group.__repr__()

class FunctionStructure(AlgebraicStructure):
    """
    An AlgebraicStructure for a monoid of functions under composition.
    The elements of this structure are actually matrices representing the actions of the functions in question.
    The operation is matrix multiplication, which naturally corresponds to function composition.

    Attributes:
        n (int): The size of the underlying set.
    """

    def __init__(self, n):
        """
        Args:
            n (int): The size of the underlying set.
        """

        self.n = n
        coeffs = product(range(n),repeat=n)
        polys = frozenset(PolynomialRing(Integers(n),'x')(r) for r in coeffs)
        mats = set()
        element_names = {}
        for poly in polys:
            m = matrix(n,sparse=True)
            for i in range(n):
                m[i,poly(i)%n] = 1
            m.set_immutable()
            mats.add(m)
            element_names[m] = poly
        operation = {'*': lambda x,y: x*y}
        AlgebraicStructure.__init__(self, mats, operation, element_names)

    def __repr__(self):
        return 'the monoid of functions on a set of order {}'.format(self.n)