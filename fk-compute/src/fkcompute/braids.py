### this file contains functions for working with braids (braid closures)
### a braid is a list of integers, where each integer represents a crossing
### CONVENTION: strand 0 is the top strand, strand 1 is the second strand, etc.
### CONVENTION: generators are indexed starting from 1, to number of strands - 1
### positive integers n represents strand n-1 crossing over strand n, so 

from IPython.display import HTML, SVG, display

def is_positive_braid(braid):
    for x in braid:
        if x < 0:
            return False
    return True

def is_homogeneous_braid(braid):
    gens = set(braid)
    for x in gens:
        if -x in gens:
            return False
    return True

def strand_path(strand, braid):
    path = []
    paths = []
    x = 15
    y = 15 + 30 * strand
    path.append((x,y))
    for i, crossing in enumerate(braid):
        if crossing == strand + 1:
            strand = strand + 1
            x = x + 30
            y = y + 30
        elif crossing == strand:
            strand = strand - 1
            path.append((x+10,y-10))
            paths.append(path)
            path = [(x+20,y-20)]
            x = x + 30
            y = y - 30
        elif crossing == -strand - 1:
            strand = strand + 1
            path.append((x+10,y+10))
            paths.append(path)
            path = [(x+20,y+20)]
            x = x + 30
            y = y + 30
        elif crossing == -strand:
            strand = strand - 1
            x = x + 30
            y = y - 30
        else:
            x = x + 30

        path.append((x,y))
    paths.append(path)
    return paths

                  
def braid_svg(braid, state_annotations=[], crossing_annotations=[]):
    nstrands = max(abs(x) for x in braid) + 1
    ncrossings = len(braid)
    paths = []
    for strand in range(0,nstrands):
        for path in strand_path(strand, braid):
            paths.append(path)

    def path_ds(path):
        return "M " + " L ".join(f"{x} {y}" for x,y in path)

    def path_svg(path):
        return f"""<path d="{path_ds(path)}" stroke="black" fill="none" stroke-width="2"/>"""


    width = 30 * (ncrossings + 1)
    height = 30 * (nstrands)

    annotation_svgs = []
    for (i,j,text) in state_annotations: 
        annotation_svgs.append(f"""<text x="{12+30*j}" y="{10+30*i}" font-family="monospace" font-size="8" fill="black">{text}</text>""")
    for (i,j,text) in crossing_annotations: 
        annotation_svgs.append(f"""<text x="{25+30*j}" y="{20+30*i}" font-family="monospace" font-size="8" fill="gray">{text}</text>""")

    
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{1.5*width}px">""" +  "\n".join(path_svg(p) for p in paths) +  "\n".join(annotation_svgs) + """</svg>"""
    return svg

def draw_braid(braid, state_annotations=[],crossing_annotations=[]):
    display(HTML(braid_svg(braid, state_annotations=state_annotations, crossing_annotations=crossing_annotations)))