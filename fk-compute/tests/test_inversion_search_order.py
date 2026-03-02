from fkcompute.inversion.search import _free_signs_from_rank


def _fmt(signs: list[int]) -> str:
    return "".join("+" if s == 1 else "-" for s in signs)


def test_free_signs_from_rank_n4_prefix() -> None:
    got = [_fmt(_free_signs_from_rank(i, 4)) for i in range(6)]
    assert got == ["++++", "+++-", "++-+", "+-++", "-+++", "++--"]


def test_free_signs_from_rank_orders_by_negative_count_then_numeric() -> None:
    n = 6
    total = 1 << n
    prev = None
    for r in range(total):
        signs = _free_signs_from_rank(r, n)
        bits = [0 if s == 1 else 1 for s in signs]
        weight = sum(bits)
        num = int("".join(str(b) for b in bits), 2)
        key = (weight, num)
        if prev is not None:
            assert prev <= key
        prev = key
