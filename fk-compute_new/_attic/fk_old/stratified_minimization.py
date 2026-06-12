from relations import *
from typing import List, Tuple


def decompose_braid_to_band_generators(braid: List[int]) -> List[Tuple[int, int, int, str]]:
    """
    Decomposes a braid into its constituent band generators.

    A band generator is a positive crossing conjugated by incrementing or
    decrementing positive crossings: (σᵢ₊₁...σᵢ₊ₖ) σᵢ (σᵢ₊₁...σᵢ₊ₖ)⁻¹

    Args:
        braid: List of integers representing a braid as product of band generators

    Returns:
        List of tuples (center, min_crossing, max_crossing, sign) for each band generator,
        where center is the conjugated crossing and sign indicates '+' for ascending or '-' for descending conjugations
    """
    band_generators = []
    i = 0
    n = len(braid)

    while i < n:
        # Try to match a band generator pattern starting at position i
        # First try w⁻¹ σᵢ w pattern
        band_info = _extract_band_generator_at(braid, i)
        
        # If not found, try w σᵢ w⁻¹ pattern
        if band_info is None:
            band_info = _extract_band_generator_w_center_winv_at(braid, i)

        if band_info is not None:
            center, min_crossing, max_crossing, sign, consumed_length = band_info
            band_generators.append((center, min_crossing, max_crossing, sign))
            i += consumed_length
        else:
            # Single crossing - treat as degenerate band generator
            center = abs(braid[i])
            band_generators.append((center, center, center, "0"))
            i += 1

    return band_generators


def _extract_band_generator_at(braid: List[int], start: int) -> Tuple[int, int, int, str, int]:
    """
    Try to extract a band generator starting at the given position.

    A band generator has the form: w⁻¹ σᵢ w where w is either:
    - Ascending: σᵢ₊₁σᵢ₊₂...σᵢ₊ₖ (positive conjugation)
    - Descending: σᵢ₋₁σᵢ₋₂...σᵢ₋ₖ (negative conjugation)

    Returns:
        Tuple (center, min_crossing, max_crossing, sign, consumed_length) if a band
        generator is found, None otherwise. Sign is '+' for ascending, '-' for descending.
    """
    n = len(braid)
    if start >= n:
        return None

    # Try different band generator lengths
    for band_length in range(1, (n - start) // 2 + 1):
        if start + 2 * band_length >= n:
            break

        # Extract potential inverse conjugator (first part)
        inverse_conjugator = braid[start:start + band_length]

        # Check if there's a center crossing
        if start + band_length >= n:
            continue

        center = braid[start + band_length]
        if center <= 0:  # Center must be positive
            continue

        # Extract the forward conjugator (last part)
        forward_conjugator = braid[start + band_length + 1:start + band_length + 1 + band_length]

        if len(forward_conjugator) != band_length:
            continue

        # Try ascending conjugation (incrementing from center)
        if _is_valid_incremental_conjugator(forward_conjugator, center):
            expected_inverse = [-x for x in reversed(forward_conjugator)]
            if inverse_conjugator == expected_inverse:
                # Check if all crossings form a contiguous block and have immediate adjacency
                all_crossings = [center] + forward_conjugator
                if _forms_contiguous_block(all_crossings) and _has_immediate_adjacency_in_pattern(forward_conjugator, center, "w_inv_center_w"):
                    min_crossing = min(abs(x) for x in all_crossings)
                    max_crossing = max(abs(x) for x in all_crossings)
                    consumed_length = 2 * band_length + 1
                    return (abs(center), min_crossing, max_crossing, "+", consumed_length)

        # Try descending conjugation (decrementing from center)
        if _is_valid_decremental_conjugator(forward_conjugator, center):
            expected_inverse = [-x for x in reversed(forward_conjugator)]
            if inverse_conjugator == expected_inverse:
                # Check if all crossings form a contiguous block and have immediate adjacency
                all_crossings = [center] + forward_conjugator
                if _forms_contiguous_block(all_crossings) and _has_immediate_adjacency_in_pattern(forward_conjugator, center, "w_inv_center_w"):
                    min_crossing = min(abs(x) for x in all_crossings)
                    max_crossing = max(abs(x) for x in all_crossings)
                    consumed_length = 2 * band_length + 1
                    return (abs(center), min_crossing, max_crossing, "-", consumed_length)

    return None


def _is_valid_incremental_conjugator(sequence: List[int], center: int) -> bool:
    """
    Check if a sequence forms a valid incremental conjugator for a band generator.

    The conjugator must be σᵢ₊₁, σᵢ₊₂, ..., σᵢ₊ₖ where i is the center crossing.
    This ensures that all crossings (center + conjugator) form a contiguous block.
    """
    if not sequence:
        return True

    # All must be positive
    if any(x <= 0 for x in sequence):
        return False

    # Must start at center + 1 and increment sequentially to form contiguous block
    expected_sequence = list(range(center + 1, center + 1 + len(sequence)))

    return sequence == expected_sequence


def _is_valid_decremental_conjugator(sequence: List[int], center: int) -> bool:
    """
    Check if a sequence forms a valid decremental conjugator for a band generator.

    The conjugator must be σᵢ₊ₖ, σᵢ₊ₖ₋₁, ..., σᵢ₊₁ where i is the center crossing.
    This means each element must be positive and decrement sequentially from center+k down to center+1.
    """
    if not sequence:
        return True

    # All must be positive
    if any(x <= 0 for x in sequence):
        return False

    # Must start at center + len(sequence) and decrement sequentially down to center + 1
    expected_sequence = list(range(center + len(sequence), center, -1))

    return sequence == expected_sequence


def _forms_contiguous_block(crossings: List[int]) -> bool:
    """
    Check if a list of crossings forms a contiguous block (all adjacent to each other).
    
    Args:
        crossings: List of positive crossing indices
        
    Returns:
        True if the crossings form a contiguous block
    """
    if len(crossings) <= 1:
        return True
    
    # Sort the absolute values to check contiguity
    sorted_crossings = sorted(abs(x) for x in crossings)
    
    # Check if they form a consecutive sequence
    for i in range(1, len(sorted_crossings)):
        if sorted_crossings[i] != sorted_crossings[i-1] + 1:
            return False
    
    return True


def _has_immediate_adjacency_in_pattern(forward_conjugator: List[int], center: int, pattern_type: str) -> bool:
    """
    Check if the crossing immediately surrounding the center in the braid pattern is adjacent to the center.
    
    For a valid band generator:
    - In w⁻¹ σᵢ w pattern: the last element of w must be adjacent to center
    - In w σᵢ w⁻¹ pattern: the last element of w must be adjacent to center
    
    Args:
        forward_conjugator: The forward conjugator sequence  
        center: The center crossing
        pattern_type: Either "w_inv_center_w" or "w_center_w_inv"
        
    Returns:
        True if the immediate neighbor crossing is adjacent to center
    """
    if not forward_conjugator:
        return True
    
    # In both patterns, the crossing immediately before/after center comes from 
    # the last element of the forward conjugator
    immediate_neighbor = abs(forward_conjugator[-1])
    return abs(immediate_neighbor - center) == 1


def _extract_band_generator_w_center_winv_at(braid: List[int], start: int) -> Tuple[int, int, int, str, int]:
    """
    Try to extract a band generator with pattern w σᵢ w⁻¹ starting at the given position.

    This handles cases like σ₃ σ₂ σ₁ σ₂⁻¹ σ₃⁻¹ where w = σ₃ σ₂ and center = σ₁.

    Returns:
        Tuple (center, min_crossing, max_crossing, sign, consumed_length) if a band
        generator is found, None otherwise. Sign is '+' for ascending, '-' for descending.
    """
    n = len(braid)
    if start >= n:
        return None

    # Try different band generator lengths
    for band_length in range(1, (n - start) // 2 + 1):
        if start + 2 * band_length >= n:
            break

        # Extract potential forward conjugator (first part)
        forward_conjugator = braid[start:start + band_length]

        # Check if there's a center crossing
        if start + band_length >= n:
            continue

        center = braid[start + band_length]
        if center <= 0:  # Center must be positive
            continue

        # Extract the inverse conjugator (last part)
        inverse_conjugator = braid[start + band_length + 1:start + band_length + 1 + band_length]

        if len(inverse_conjugator) != band_length:
            continue

        # Try ascending conjugation (incrementing from center)
        if _is_valid_incremental_conjugator(forward_conjugator, center):
            expected_inverse = [-x for x in reversed(forward_conjugator)]
            if inverse_conjugator == expected_inverse:
                # Check if all crossings form a contiguous block and have immediate adjacency
                all_crossings = [center] + forward_conjugator
                if _forms_contiguous_block(all_crossings) and _has_immediate_adjacency_in_pattern(forward_conjugator, center, "w_center_w_inv"):
                    min_crossing = min(abs(x) for x in all_crossings)
                    max_crossing = max(abs(x) for x in all_crossings)
                    consumed_length = 2 * band_length + 1
                    return (abs(center), min_crossing, max_crossing, "+", consumed_length)

        # Try descending conjugation (decrementing from center)
        if _is_valid_decremental_conjugator(forward_conjugator, center):
            expected_inverse = [-x for x in reversed(forward_conjugator)]
            if inverse_conjugator == expected_inverse:
                # Check if all crossings form a contiguous block and have immediate adjacency
                all_crossings = [center] + forward_conjugator
                if _forms_contiguous_block(all_crossings) and _has_immediate_adjacency_in_pattern(forward_conjugator, center, "w_center_w_inv"):
                    min_crossing = min(abs(x) for x in all_crossings)
                    max_crossing = max(abs(x) for x in all_crossings)
                    consumed_length = 2 * band_length + 1
                    return (abs(center), min_crossing, max_crossing, "-", consumed_length)

    return None


if __name__ == "__main__":
    # Example usage (using w⁻¹gw conjugation convention)

    # Example 1: Simple band generator σ₂⁻¹ σ₁ σ₂ (conjugation of σ₁ by σ₂)
    braid1 = [-2, 1, 2]
    print("Braid 1:", braid1, "→ σ₂⁻¹ σ₁ σ₂")
    print("Band Generators:", decompose_braid_to_band_generators(braid1))
    print()

    # Example 2: NOT a band generator - σ₂ σ₃ σ₁ σ₃⁻¹ σ₂⁻¹ (wrong order for w⁻¹gw)
    braid2 = [2, 3, 1, -3, -2]
    print("Braid 2:", braid2, "→ σ₂ σ₃ σ₁ σ₃⁻¹ σ₂⁻¹ (NOT a band generator)")
    print("Band Generators:", decompose_braid_to_band_generators(braid2))
    print("DEBUG: _has_immediate_adjacency([2, 3], 1) =", _has_immediate_adjacency([2, 3], 1))
    print()

    # Example 3: Valid band generator σ₃⁻¹ σ₂⁻¹ σ₁ σ₂ σ₃ (conjugation of σ₁ by σ₂σ₃)
    braid3 = [-3, -2, 1, 2, 3]
    print("Braid 3:", braid3, "→ σ₃⁻¹ σ₂⁻¹ σ₁ σ₂ σ₃")
    print("Band Generators:", decompose_braid_to_band_generators(braid3))
    print()

    # Example 4: NOT a band generator - σ₃ σ₁ σ₃⁻¹ (conjugator doesn't increment from center)
    braid4 = [3, 1, -3]
    print("Braid 4:", braid4, "→ σ₃ σ₁ σ₃⁻¹ (NOT a band generator)")
    print("Band Generators:", decompose_braid_to_band_generators(braid4))
    print()

    # Example 5: Valid band generator σ₄⁻¹ σ₃⁻¹ σ₂ σ₃ σ₄ (conjugation of σ₂ by σ₃σ₄)
    braid5 = [-4, -3, 2, 3, 4]
    print("Braid 5:", braid5, "→ σ₄⁻¹ σ₃⁻¹ σ₂ σ₃ σ₄")
    print("Band Generators:", decompose_braid_to_band_generators(braid5))

    # Example 6: Valid band generator with descending conjugation
    braid6 = [3, 2, 1, -2, -3]
    print("Braid 6:", braid6, "→ σ₃ σ₂ σ₁ σ₂⁻¹ σ₃⁻¹")
    print("Band Generators:", decompose_braid_to_band_generators(braid6))
