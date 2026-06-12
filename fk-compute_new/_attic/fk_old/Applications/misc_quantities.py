from ..braidstates_links_old import BraidStates
from .resolutions import *

# cf. arXiv:1002.0898v1 : (2010, Dasbach and Lowrance) : Turaev Genus, Knot Signature, and the Knot Homology Concordance Invariants
def DL(
    bs : BraidStates, 
    L : int
):
    '''
    Checks the Dasbach-Lowrance (DL) inequalities for DL concordance 
    (note: it has been called into question whether or not their DL
    conditions imply concordance invariance) invariants against the 
    quantity L given a BratidStates object bs.
    '''
    lower = b_resolution_components(bs)[0] - sum([x < 0 for x in bs.braid]) - 1
    upper = 1 + sum([x > 0 for x in bs.braid]) - a_resolution_components(bs)[0]
    return L <= upper and L >= lower

# cf. arXiv:0908.2745v2 : (2009, Lobb) : Computable Bounds for Rasmussen’s Concordance Invariant
def Lobb(
    bs : BraidStates, 
    L : int
):
    '''
    Checks Lobb's inequalities for the Rasmussen invariant 
    against the quantity L given a BratidStates object bs.
    '''
    lower = sum([-x not in bs.braid for x in range(1, bs.n_strands)])
    upper = sum([x not in bs.braid for x in range(1, bs.n_strands)])
    U = bs.writhe + bs.n_strands + 1 - 2 * lower
    Delta = bs.n_strands + 1 - upper - lower
    return L <= U and L >= U - 2 * Delta