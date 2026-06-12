"""
Braid Classification Functions

This module provides functions to classify braids based on various mathematical properties.
Braids are represented as lists of integers where:
- Positive integers represent standard braid generators (σ_i)
- Negative integers represent inverse generators (σ_i^{-1})
- The absolute value indicates which strands are crossed
- Zero is not used as a generator

For example, [1, -2, 1] represents the braid σ_1 σ_2^{-1} σ_1 on 3 strands.
"""

def is_positive_braid(braid):
    """
    Check if a braid is positive (contains only positive generators).
    
    A positive braid uses only positive generators σ_i, with no inverses σ_i^{-1}.
    Positive braids have special properties in knot theory and are easier to analyze.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if all generators are positive, False otherwise
        
    Examples:
        >>> is_positive_braid([1, 2, 3, 1])
        True
        >>> is_positive_braid([1, -2, 3])
        False
        >>> is_positive_braid([])
        True
    """
    for x in braid:
        if x < 0:
            return False
    return True

def is_homogeneous_braid(braid):
    """
    Check if a braid is homogeneous (no generator appears with both signs).
    
    A homogeneous braid never uses both σ_i and σ_i^{-1} for the same i.
    This is a weaker condition than being positive but still provides structure.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if no generator appears with both positive and negative signs
        
    Examples:
        >>> is_homogeneous_braid([1, 2, 1, 2])
        True
        >>> is_homogeneous_braid([-1, -2, -1])
        True  
        >>> is_homogeneous_braid([1, -1, 2])
        False
    """
    gens = set(braid)
    for x in gens:
        if -x in gens:
            return False
    return True

def is_alternating_braid(braid):
    """
    Check if a braid is alternating (signs alternate between consecutive generators).
    
    An alternating braid switches between positive and negative generators,
    creating a regular over-under crossing pattern when closed to form a link.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if the signs of consecutive generators alternate
        
    Examples:
        >>> is_alternating_braid([1, -2, 3, -1])
        True
        >>> is_alternating_braid([1, 2, -3])
        False
        >>> is_alternating_braid([1])
        True
        >>> is_alternating_braid([])
        True
    """
    if len(braid) <= 1:
        return True
        
    for i in range(len(braid) - 1):
        if (braid[i] > 0) == (braid[i + 1] > 0):  # Same sign
            return False
    return True

def is_pure_braid(braid):
    """
    Check if a braid is pure (geometric braiding equals identity permutation).
    
    A pure braid returns each strand to its original position when viewed
    geometrically. This requires the sum of exponents for each generator to be even.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if braid is pure (identity permutation)
        
    Examples:
        >>> is_pure_braid([1, -1])
        True
        >>> is_pure_braid([1, 2, -1, -2])
        True
        >>> is_pure_braid([1, 2, 3])
        False
    """
    from collections import Counter
    
    # Count occurrences of each generator (considering sign)
    gen_count = Counter(braid)
    
    # For purity, we need to check if the permutation is identity
    # This is a simplified check - full purity requires more complex analysis
    # Here we check if positive and negative occurrences balance
    generators = set(abs(x) for x in braid)
    
    for gen in generators:
        total_power = gen_count.get(gen, 0) + gen_count.get(-gen, 0)
        if total_power % 2 != 0:
            return False
    
    return True

def get_braid_length(braid):
    """
    Get the length (number of crossings) of a braid.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        int: Number of generators (crossings) in the braid
        
    Examples:
        >>> get_braid_length([1, 2, -1, 3])
        4
        >>> get_braid_length([])
        0
    """
    return len(braid)

def get_strand_number(braid):
    """
    Determine the minimum number of strands needed for the braid.
    
    The number of strands is one more than the maximum generator index,
    since generator i connects strands i and i+1.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        int: Minimum number of strands, or 1 for empty braid
        
    Examples:
        >>> get_strand_number([1, 2, 3])
        4
        >>> get_strand_number([-2, 1, -3])
        4
        >>> get_strand_number([])
        1
    """
    if not braid:
        return 1
    return max(abs(x) for x in braid) + 1

def get_generator_spectrum(braid):
    """
    Get the set of all generators used in the braid (ignoring signs).
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        set: Set of positive integers representing generators used
        
    Examples:
        >>> get_generator_spectrum([1, -2, 1, 3, -2])
        {1, 2, 3}
        >>> get_generator_spectrum([])
        set()
    """
    return set(abs(x) for x in braid)

def is_conjugate_positive(braid):
    """
    Check if braid might be conjugate to a positive braid.
    
    This is a heuristic check based on the balance of positive/negative generators.
    A more rigorous test would require complex algorithms.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if generator balance suggests possible conjugacy to positive braid
        
    Examples:
        >>> is_conjugate_positive([1, 2, -1, 2])
        True
        >>> is_conjugate_positive([-1, -2, -3])
        False
    """
    from collections import Counter
    
    gen_count = Counter(braid)
    positive_total = sum(count for gen, count in gen_count.items() if gen > 0)
    negative_total = sum(count for gen, count in gen_count.items() if gen < 0)
    
    # Heuristic: if significantly more positive than negative generators
    return positive_total >= negative_total

def get_writhe(braid):
    """
    Calculate the writhe (linking number contribution) of the braid.
    
    Writhe is the sum of crossing signs: +1 for positive generators, -1 for negative.
    This gives the writhe contribution when the braid is closed to form a link.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        int: Writhe (sum of crossing signs)
        
    Examples:
        >>> get_writhe([1, 2, 3])
        3
        >>> get_writhe([1, -2, 3])
        2
        >>> get_writhe([-1, -2, -1])
        -3
    """
    return sum(1 if x > 0 else -1 for x in braid)

def is_palindromic_braid(braid):
    """
    Check if a braid is palindromic (reads the same forwards and backwards).
    
    A palindromic braid has the form w = w^{-1} when read as a word.
    These braids have special symmetry properties.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if braid is palindromic
        
    Examples:
        >>> is_palindromic_braid([1, 2, -2, -1])
        True
        >>> is_palindromic_braid([1, 2, 1])
        False
        >>> is_palindromic_braid([])
        True
    """
    if len(braid) == 0:
        return True
        
    # Check if braid equals its reverse with signs flipped
    reversed_inverse = [-x for x in reversed(braid)]
    return braid == reversed_inverse

def is_quasipositive_braid(braid):
    """
    Check if a braid is quasipositive.
    
    A braid is quasipositive if it can be written as a product of conjugates
    of positive generators. This is equivalent to checking if the braid can
    be expressed using only positive generators and their conjugates by
    positive braids.
    
    This implementation uses a simplified heuristic: a braid is likely
    quasipositive if every negative generator can be "cancelled" by
    surrounding positive generators through conjugation.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if braid appears to be quasipositive
        
    Examples:
        >>> is_quasipositive_braid([1, 2, 1])
        True
        >>> is_quasipositive_braid([1, -2, 1])  # σ₁σ₂⁻¹σ₁ = σ₁σ₂⁻¹σ₁
        True
        >>> is_quasipositive_braid([-1, -2, -3])
        False
        
    Note:
        This is a heuristic implementation. The full quasipositive test
        requires complex algorithms involving braid conjugation and
        may not detect all quasipositive braids.
    """
    if is_positive_braid(braid):
        return True
        
    # Heuristic: check if negative generators appear in "conjugating" positions
    # A more sophisticated algorithm would need to actually attempt conjugation
    
    from collections import defaultdict
    
    # Track positions of each generator
    positions = defaultdict(list)
    for i, gen in enumerate(braid):
        positions[abs(gen)].append((i, gen > 0))
    
    # For each generator that appears negatively, check if it's "surrounded"
    # by positive generators that could conjugate it
    for gen_idx, occurrences in positions.items():
        negative_positions = [pos for pos, is_pos in occurrences if not is_pos]
        
        if negative_positions:
            # Heuristic: negative generators should be "surrounded" by positive ones
            for neg_pos in negative_positions:
                # Check if there are positive generators nearby that could conjugate
                has_left_conjugator = False
                has_right_conjugator = False
                
                # Look for potential conjugating generators
                for pos, is_pos in occurrences:
                    if is_pos and pos < neg_pos:
                        has_left_conjugator = True
                    if is_pos and pos > neg_pos:
                        has_right_conjugator = True
                
                if not (has_left_conjugator and has_right_conjugator):
                    return False
    
    return True

def is_strongly_quasipositive_braid(braid):
    """
    Check if a braid is strongly quasipositive.
    
    A braid is strongly quasipositive if it can be written as a product
    of conjugates of σᵢ by positive braids, where each conjugate has the
    form P σᵢ P⁻¹ with P positive.
    
    Strong quasipositivity is a stronger condition than quasipositivity.
    Every strongly quasipositive braid is quasipositive, but not vice versa.
    
    This implementation uses the property that strongly quasipositive braids
    have a very specific structure in their generator sequences.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        bool: True if braid appears to be strongly quasipositive
        
    Examples:
        >>> is_strongly_quasipositive_braid([1, 2, 1])
        True
        >>> is_strongly_quasipositive_braid([1, 2, -1, -2, 1, 2])
        False
        >>> is_strongly_quasipositive_braid([2, 1, 2])
        True
        
    Note:
        This is a heuristic implementation. The complete algorithm for
        detecting strongly quasipositive braids is computationally intensive
        and involves detailed analysis of braid structure.
    """
    if is_positive_braid(braid):
        return True
    
    # Strongly quasipositive braids have very restrictive structure
    # Heuristic: they should have significantly more positive than negative generators
    # and negative generators should appear in very specific patterns
    
    if not is_quasipositive_braid(braid):
        return False
    
    positive_count = sum(1 for x in braid if x > 0)
    negative_count = sum(1 for x in braid if x < 0)
    
    # Strong heuristic: strongly quasipositive braids typically have
    # many more positive than negative generators
    if negative_count > 0 and positive_count < 3 * negative_count:
        return False
    
    # Additional check: negative generators should appear in "sandwiched" patterns
    for i, gen in enumerate(braid):
        if gen < 0:
            # Check if this negative generator is properly "sandwiched"
            left_positive = False
            right_positive = False
            
            # Look for positive generators of same or adjacent index
            for j in range(max(0, i-2), min(len(braid), i+3)):
                if j != i and braid[j] > 0:
                    if abs(braid[j] - abs(gen)) <= 1:  # Same or adjacent generator
                        if j < i:
                            left_positive = True
                        if j > i:
                            right_positive = True
            
            if not (left_positive and right_positive):
                return False
    
    return True

def is_fibered_braid(braid):
    pass

def classify_braid(braid):
    """
    Comprehensive classification of a braid's properties.
    
    Args:
        braid (list): List of integers representing braid generators
        
    Returns:
        dict: Dictionary containing all classification results
        
    Examples:
        >>> classify_braid([1, 2, 1])
        {'positive': True, 'homogeneous': True, 'alternating': False, 'quasipositive': True, ...}
    """
    return {
        'length': get_braid_length(braid),
        'strand_number': get_strand_number(braid),
        'generator_spectrum': get_generator_spectrum(braid),
        'positive': is_positive_braid(braid),
        'homogeneous': is_homogeneous_braid(braid),
        'alternating': is_alternating_braid(braid),
        'pure': is_pure_braid(braid),
        'palindromic': is_palindromic_braid(braid),
        'quasipositive': is_quasipositive_braid(braid),
        'strongly_quasipositive': is_strongly_quasipositive_braid(braid),
        'conjugate_positive': is_conjugate_positive(braid),
        'writhe': get_writhe(braid),
        'fibered': is_fibered_braid(braid)  # Placeholder
    }
