#!/usr/bin/env python

"""Classes and modules for cylindrical BOUT++ Clebsch system"""

from collections import deque
from IPython.display import display
from sympy import symbols, simplify
from sympy import Eq, Derivative, Matrix
from sympy import sqrt
from sympy import Rational as R

from common import cartCoord, cartMap, cylCoord, cylMap, rho, theta, z

# TODO: Consider making a factory pattern as creating a new vector at
#       the moment takes some time

#{{{CartVec
class CartVec(object):
    """
    A vector class for cartesian coordinates.
    Inherits from 'object' class. This makes it a 'new-style class'
    """

    def __init__(self, x=0, y=0, z=0):
        """Constructor"""
        self.x = x
        self.y = y
        self.z = z

    def __mul__(self, other):
        """Muliplication overload"""
        out = 0
        for i in cartCoord:
            out += self.__getattribute__(str(i))*other.__getattribute__(str(i))
        return out

    def __truediv__(self, other):
        """Division overload"""
        if CartVec in type(other).__bases__:
            raise RuntimeError("Division by vector is udefined")
        else:
            # Create a new instance of the input object
            out = self.__new__(type(self))
            out.__init__()
            for i in cartCoord:
                out.__setattr__(str(i), self.__getattribute__(str(i))/ other)
            return out

    def len(self):
        """Lenght of a vector"""
        return simplify(sqrt(self*self))
#}}}

#{{{CylVec
class CylVec(object):
    """
    A base class for cylindrical vectors using the Clebsch system B =
    e^3 sub.
    Inherits from 'object' class. This makes it a 'new-style class'
    """

    #{{{__init__
    def __init__(self, rho=0, z=0, theta=0, covariant=True):
        """
        Constructor

        Covariant refers to the vector ELEMENTS
        """
        self.rho       = rho
        self.z         = z
        self.theta     = theta
        self.covariant = covariant
        # Call method from subclass
        self._createBasisVec()
        # Calculate the metrics
        self._calcMetrics()
    #}}}

    #{{{__mul__
    def __mul__(self, other):
        """Muliplication overload"""
        if "Vec" in type(other).__name__:
            #{{{ If other is a vector
            # Security check
            if type(self).__name__ != type(other).__name__:
                string = "Operation between {0} and {1} forbidden".\
                        format(type(self).__name__, type(other).__name__)
                raise RuntimeError(string)
            # If co and contravariant
            if (self.covariant and not(other.covariant)) or\
               (not(self.covariant) and other.covariant):
                out = 0
                for coord in cylCoord:
                    out += self. __getattribute__(str(coord))*\
                           other.__getattribute__(str(coord))
            # If covariant
            elif self.covariant and other.covariant:
                # Equation (2.5.33d)
                # Components are covariant
                # => Basis vectors are contravariant
                # => Contravariant metric tensor used
                out = 0
                for i in cylCoord:
                    for j in cylCoord:
                        out += self. __getattribute__(str(i))*\
                               other.__getattribute__(str(j))*\
                               self.g[i][j]
            # If contravariant
            elif not(self.covariant) and not(other.covariant):
                # Equation (2.5.33c)
                # Components are contravariant
                # => Basis vectors are covariant
                # => Covariant metric tensor used
                out = 0
                for i in cylCoord:
                    for j in cylCoord:
                        out += self. __getattribute__(str(i))*\
                               other.__getattribute__(str(j))*\
                               self.g_[i][j]
            #}}}
        elif type(type(other).__name__) == str:
            #{{{If other is a sympy object
            out = self.__new__(type(self))
            out.__init__()
            for i in cylCoord:
                out.__setattr__(str(i), self.__getattribute__(str(i))*other)

            out.covariant = self.covariant
            #}}}
        else:
            raise RuntimeError('Operation not implemented')

        return out
    #}}}

    #{{{__xor___
    def __xor__(self, other):
        """Cross product"""
        # Security check
        # Same type of vector
        if type(self).__name__ != type(other).__name__:
            string = "Operation between {0} and {1} forbidden".\
                     format(type(self).__name__, type(other).__name__)
            raise RuntimeError(string)
        # Conversion so that both are co/contravariant
        if self.covariant and not(other.covariant):
            other = other.convertToCovariant()
        elif not(self.covariant) and other.covariant:
            self = self.convertToCovariant()
        # Create a new instance of the input object
        out = self.__new__(type(self))
        out.__init__()
        # Equation (2.5.39) in D'Haeseleer
        # NOTE: Formula has covariant COMPONENTS as input
        #       => Projected onto a covariant BASIS VECTOR
        #       => Resulting vector COMPONENTS will be contravariant
        #       And vice versa
        if self.covariant:
            out.covariant = False
            factor = (1/self.J)
        else:
            out.covariant = True
            factor = self.J
        # Deque used for cyclic permutations
        index = deque(cylCoord)
        for _ in range(len(index)):
            out.__setattr__(str(index[2]),
                factor*(
                      self.__getattribute__(str(index[0]))*\
                      other.__getattribute__(str(index[1]))\
                      -\
                      self.__getattribute__(str(index[1]))*\
                      other.__getattribute__(str(index[0]))
                ))
            # Rotate the indices
            index.rotate(1)
        return out
    #}}}

    #{{{__truediv___
    def __truediv__(self, other):
        """Division overload"""
        if CylVec in type(other).__bases__:
            raise RuntimeError("Division by vector is udefined")
        else:
            # Create a new instance of the input object
            out = self.__new__(type(self))
            out.__init__()
            for i in cylCoord:
                out.__setattr__(str(i), self.__getattribute__(str(i))/ other)

            out.covariant = self.covariant
            return out
    #}}}

    #{{{__sub__
    def __sub__(self, other):
        """Subtraction overload"""
        if "Vec" in type(other).__name__:
            #{{{ If other is a vector
            # Security check
            if type(self).__name__ != type(other).__name__:
                string = "Operation between {0} and {1} forbidden".\
                        format(type(self).__name__, type(other).__name__)
                raise RuntimeError(string)
            # If co and contravariant
            if (self.covariant and not(other.covariant)) or\
               (not(self.covariant) and other.covariant):
                if self.covariant:
                    selfStr = 'covariant'
                else:
                    selfStr = 'contravariant'
                if other.covariant:
                    otherStr = 'covariant'
                else:
                    otherStr = 'contravariant'

                string = "Undefined operation for minus between {0} and {1}".\
                        format(selfStr, otherStr)
                raise RuntimeError(string)

            else:
                # Create a new instance of the input object
                out = self.__new__(type(self))
                out.__init__()
                out.covariant = self.covariant

                for i in cylCoord:
                    out.__setattr__(str(i),\
                                    self.__getattribute__(str(i)) -
                                    other.__getattribute__(str(i))\
                                    )
            #}}}
        else:
            raise RuntimeError('Operation not implemented')

        return out
    #}}}

    #{{{__neg__
    def __neg__(self):
        """Overload negation"""
        # Create a new instance of the input object
        out = self.__new__(type(self))
        out.__init__()
        for i in cylCoord:
            out.__setattr__(str(i), -self.__getattribute__(str(i)))

        out.covariant = self.covariant
        return out
    #}}}

    #{{{ _calcMetrics
    def _calcMetrics(self):
        """Calculates the co and contravariant metrics"""
        #{{{Covariant g
        # Initialize covariant g
        self.g_ = {}
        # Initialize list which will be the input in the matrix
        matrix = []
        for i in cylCoord:
            # Initialize second coord
            self.g_[i] = {}
            # Initialize a row in the matrix
            row = []
            for j in cylCoord:
                self.g_[i][j] = simplify(\
                                    self.__getattribute__(str(i) + 'CovBasis')*\
                                    self.__getattribute__(str(j) + 'CovBasis')\
                                    )
                # Append the element to the row
                row.append(self.g_[i][j])
            # Append the row to the matrix
            matrix.append(row)
        # Define the matrix (makes it easy to display)
        self.covMetrics = Matrix(matrix)
        #}}}

        #{{{Contravariant g
        # Initialize contravariant g
        self.g = {}
        # Initialize list which will be the input in the matrix
        matrix = []
        for i in cylCoord:
            # Initialize second coord
            self.g[i] = {}
            # Initialize a row in the matrix
            row = []
            for j in cylCoord:
                self.g[i][j] = simplify(\
                                self.__getattribute__(str(i) + 'ConBasis')*\
                                self.__getattribute__(str(j) + 'ConBasis')\
                               )
                # Append the element to the row
                row.append(self.g[i][j])
            # Append the row to the matrix
            matrix.append(row)
        # Define the matrix (makes it easy to display)
        self.conMetrics = Matrix(matrix)
        #}}}

        #{{{ Calculate J
        # Define the jacobian
        self.J = simplify(sqrt(self.covMetrics.det()))
        #}}}

        #{{{Calculate Christoffel
        #The Christoffel symbols reads
        #$$
        #\Gamma^i_{kl}=\frac{1}{2}g^{im}
        #\left(
        #       \frac{\partial g_{mk}}{\partial x^l} + \frac{\partial
        #       g_{ml}}{\partial x^k} - \frac{\partial g_{kl}}{\partial x^m}
        # \right)
        #$$
        self.C = dict.fromkeys(cylCoord)

        for i in cylCoord:
            self.C[i] = dict.fromkeys(cylCoord)
            # Initialize list which will be the input in the matrix
            matrix = []
            for j in cylCoord:
                # Initialize a row in the matrix
                row = []
                self.C[i][j] = dict.fromkeys(cylCoord)
                for k in cylCoord:
                    self.C[i][j][k] = 0
                    for l in cylCoord:
                        self.C[i][j][k] +=\
                            R(1/2)*self.g[i][l]\
                            *(\
                                self.g_[l][j].diff(k)\
                              + self.g_[l][k].diff(j)\
                              - self.g_[j][k].diff(l)\
                             )
                    self.C[i][j][k] = simplify(self.C[i][j][k])
                    # Append the element to the row
                    row.append(self.C[i][j][k])
                # Append the row to the matrix
                matrix.append(row)
            # Define the matrix (makes it easy to display)
            self.__setattr__('C'+str(i), Matrix(matrix))
        #}}}
    #}}}

    #{{{convertToCovariant
    def convertToCovariant(self):
        """Convert a vector to a contravariant vector"""
        if not(self.covariant):
            # Create a output vector
            out = self.__new__(type(self))
            out.__init__()
            # Equation (2.5.9) in D'Haeseleer
            for i in cylCoord:
                out.__setattr__(str(i), 0)
                for j in cylCoord:
                    out.__setattr__(str(i),\
                            out.__getattribute__(str(i)) +\
                            self.__getattribute__(str(j))*\
                            self.g_[i][j]
                            )
            # Make sure vector is set to covariant
            out.covariant = True

            # self = out will only be visible within this scope, thus
            # the output object must be returned
            return out
    #}}}

    #{{{convertToContravariant
    def convertToContravariant(self):
        """Convert a vector to a contravariant vector"""
        if self.covariant:
            # Create a output vector
            out = self.__new__(type(self))
            out.__init__()
            # Equation (2.5.10) in D'Haeseleer
            for i in cylCoord:
                out.__setattr__(str(i), 0)
                for j in cylCoord:
                    out.__setattr__(str(i),\
                            out.__getattribute__(str(i)) +\
                            self.__getattribute__(str(j))*\
                            self.g[i][j]
                            )
            # Make sure vector is set covariant to False
            out.covariant = False

            # self = out will only be visible within this scope, thus
            # the output object must be returned
            return out
    #}}}

    #{{{len
    def len(self):
        """Lenght of a vector"""
        return simplify(sqrt(self*self))
    #}}}

    #{{{doitVec
    def doitVec(self):
        """Sympy doit for a vector"""
        for coord in cylCoord:
            self.__setattr__(str(coord),\
                    self.__getattribute__(str(coord)).doit())
        return self
    #}}}

    #{{{simplifyVec
    def simplifyVec(self):
        """Sympy simplify for a vector"""
        for coord in cylCoord:
            self.__setattr__(str(coord),\
                    simplify(self.__getattribute__(str(coord))))
        return self
    #}}}
#}}}

#{{{ClebschVec
class ClebschVec(CylVec):
    """
    Child-class for non-normalized cylindrical Clebsch vectors.
    """

    def _createBasisVec(self):
        """
        Method for creating the basis vectors.
        That is the basis vector of a cylindrical basis vector written
        in terms of cartesian vectors.
        These will be used to calculate the metrics

        NOTE: Here covariant refers to the BASIS VECTORS, it usually
              refers to ELEMENTS OF A VECTOR
        """
        # Define the covariant basis vectors
        for cylC in cylCoord:
            # Initiate a cartesian vector
            self.__setattr__(str(cylC) + 'CovBasis', CartVec())
            for cartC in cartCoord:
                self.__getattribute__(str(cylC) + 'CovBasis').\
                                      __setattr__(str(cartC),\
                                      simplify(cartMap[cartC].diff(cylC)))

        # Define the contravariant basis vectors
        for cylC in cylCoord:
            # Initiate a cartesian vector
            self.__setattr__(str(cylC) + 'ConBasis', CartVec())
            for cartC in cartCoord:
                self.__getattribute__(str(cylC) + 'ConBasis').\
                                      __setattr__(str(cartC),\
                                      simplify(cylMap[cylC].diff(cartC)))
                # Rewrite to cylindrical coordinates
                for coord in cartMap.keys():
                    self.__getattribute__(str(cylC) + 'ConBasis').\
                         __setattr__(str(cartC),\
                            simplify(\
                                self.__getattribute__(str(cylC) + 'ConBasis').\
                                __getattribute__(str(cartC)).\
                                    subs(coord, cartMap[coord])\
                                ))
#}}}

#{{{div
def div(v):
    """
    Divergence of a vector.
    Equation (2.6.39) in D'Haeseleer
    """
    s = 0
    if v.covariant:
        v = v.convertToContravariant()
    for coord in cylCoord:
        s += Derivative(v.J*v.__getattribute__(str(coord)), coord)
    s *= (1/v.J)
    return s
#}}}

#{{{grad
def grad(s):
    """
    Gradient of a scalar.
    Equation (2.6.34) in D'Haeseleer
    """

    v = ClebschVec()
    # Basis vectors of gradient is contravariant, so elements are
    # covariant
    v.covariant = True
    for coord in cylCoord:
        v.__setattr__(str(coord), Derivative(s, coord))
    return v
#}}}

#{{{gradPerp
def gradPerp(s):
    """
    Perpendicular gradient of a scalar.
    Equation (2.6.34) in D'Haeseleer
    """

    v = ClebschVec()
    # Basis vectors of gradient is contravariant, so elements are
    # covariant
    v.covariant = True
    v.rho   = Derivative(s, rho)
    v.z     = 0
    v.theta = Derivative(s, theta)
    return v
#}}}

#{{{advVec
def advVec(v,a):
    """
    Advection of a vector

    v*grad(a)
    """

    if type(v).__name__ != "ClebschVec" or type(a).__name__ != "ClebschVec":
        message = "v and a must be ClebschVec, not {0} and {1}".\
                    format(type(v).__name__, type(a).__name__)
        raise RuntimeError(message)

    out = ClebschVec()
    # As we are dotting v and the gradient, 'a' will decide whether or not
    # out is covariant (see equation (2.6.31) and (2.6.32) in
    # D'Haeseleer
    if a.covariant:
        out.covariant = True
    else:
        out.covariant = False
    #{{{v covariant
    if v.covariant:
        # We are dealing with
        #   v_i e^i * e^j partial_j a@k e@k
        # = v_i g^{ij} partial_j a a@k e@k
        for k in cylCoord:
            v_i_gIJ_pj_a = 0
            for i in cylCoord:
                for j in cylCoord:
                    v_i_gIJ = v.__getattribute__(str(i))*v.g[i][j]
                    # Take the derivative of the vector element
                    pj_a = Derivative(a.__getattribute__(str(k)), j)
                    # Add the Christoffel symbol to pj_a
                    # NOTE:
                    # The indices here does not correspond to
                    # (2.6.31) or (2.6.32)
                    # We have:
                    # j -> k (looking at the k'th element in the vector)
                    # k -> j (derivative w.r.t j)
                    # i -> l (Christoffel sum)
                    if a.covariant:
                        # We are dealing with
                        # v_i g^{ij} partial_j a_k e^k
                        # See equation (2.6.32)
                        for l in cylCoord:
                            pj_a -= a.C[l][k][j]*\
                                     a.__getattribute__(str(l))
                    else:
                        # We are dealing with
                        # v_i g^{ij} partial_j a^k e_k
                        # See equation (2.6.31)
                        for l in cylCoord:
                            pj_a += a.C[k][l][j]*\
                                     a.__getattribute__(str(l))
                    # Add v_i_gIJ*pj_a to the sum
                    v_i_gIJ_pj_a += v_i_gIJ*pj_a
            # Insert in the vector
            out.__setattr__(str(k), v_i_gIJ_pj_a)
    #}}}
    #{{{ v contravariant
    else:
        # We are dealing with
        #   v^i e_i * e^j partial_j a@k e@k
        # = v^i partial_i a@j e@j
        for j in cylCoord:
            vI_pi_a = 0
            for i in cylCoord:
                vI = v.__getattribute__(str(i))
                # Take the derivative of the vector element
                pi_a = Derivative(a.__getattribute__(str(j)), i)
                # Add the Christoffel symbol to pi_a
                # NOTE:
                # The indices here does not correspond to
                # (2.6.31) or (2.6.32)
                # We have:
                # j -> j (looking at the j'th element in the vector)
                # k -> i (derivative w.r.t i)
                # i -> k (Christoffel sum)
                if a.covariant:
                    # We are dealing with
                    # v^i partial_i a_j e^j
                    # See equation (2.6.32)
                    for k in cylCoord:
                        pi_a -= a.C[k][j][i]*\
                                 a.__getattribute__(str(k))
                else:
                    # We are dealing with
                    # v^i partial_i a^j e_j
                    # See equation (2.6.31)
                    for k in cylCoord:
                        pi_a += a.C[j][k][i]*\
                                 a.__getattribute__(str(k))
                # Add v_i_gIJ*pj_a to the sum
                vI_pi_a += vI*pi_a
            # Insert in the vector
            out.__setattr__(str(j), vI_pi_a)
    #}}}

    return out
#}}}
