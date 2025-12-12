import yaml
import fkcompute
import multiprocessing as mp

# Load one knot configuration
with open("all_knots.yaml", "r") as f:
    all_knots = yaml.safe_load(f)

# Get the first knot
knot_name = list(all_knots.keys())[0]
knot_config = all_knots[knot_name]

print(f"Testing knot: {knot_name}")
print(f"Config: {knot_config}")

try:
    braid = knot_config.get("braid")
    inversion = {"inversion_data": knot_config.get("inversion")}
    threads = mp.cpu_count()
    degree = 4

    print(f"\nCalling fkcompute.fk with:")
    print(f"  braid: {braid}")
    print(f"  degree: {degree}")
    print(f"  name: {knot_name}")
    print(f"  threads: {threads}")
    print(f"  inversion: {inversion}")

    result = fkcompute.fk(
        braid,
        degree=degree,
        name=knot_name,
        threads=threads,
        save_data=False,
        inversion=inversion,
    )
    print(f"\nSuccess! Result: {result}")
except Exception as e:
    print(f"\nError: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
