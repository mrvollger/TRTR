This tool removes tandem repeats from ends of unaligned sequencing reads (leaving one copy). This prevents reads that don't span the repeated region from overlapping and leading to innaccurate SNPs calls.

The maximimum repeat length is adjustable (use 1 to trim only homopolymers).

The "aggressive" option should not be touched in general. Setting to 0 will prevent the program from trimming to exactly 1 copy of the repeat, instead leaving between 1 and 2 copies.

This could also be a useful first step before assembly. More testing needs to be done.
