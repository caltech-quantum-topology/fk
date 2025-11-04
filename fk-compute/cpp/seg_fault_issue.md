There's a segmentation fault issue when I run

./fk_main data/link_ilp out 

This is using rewritten code. Legacy code still works and gives the correct result:

./fk_segments_links data/link_ilp out

I believe the segmentation fault is raised in performOffsetAddition. Compare the two algorithms to find where the error lays.
