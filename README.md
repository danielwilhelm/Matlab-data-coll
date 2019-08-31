# Optimal Data Collection for RCTs

This replication package accompanies the paper "[Optimal Data Collection for Randomized Control Trials](https://www.ucl.ac.uk/~uctpdwi/papers/cwp211919.pdf)" by P. Carneiro, S. Lee, and D. Wilhelm.

To replicate Tables 5 - 8 in the paper, simply run the corresponding matlab files:

    Table5.m
    Table6.m
    Table7a.m
    Table7b.m
    Table7c.m
    Table7d.m
    Table8a.m
    Table8b.m
    Table8c.m
    Table8d.m

Each of these produces two latex tables in the folder "results-childcare" (Tables 5 and 6) or "results-schoolgrants" (Tables 7 and 8): one that reproduces the table in the paper and another that lists the names of the selected covariates (i.e. Table S4 in the supplement is the same as Table7a_varlist.tex).

Notice, however, that the code might run for a long time (especially for Tables 7 and 8). To speed up computation time you could decrease the value "nN" in each of these files, which governs the number sample sizes that the code searches over. Of course, the results might then slightly deviate from those in the paper.

Similarly, you can replicate the figure and the tables in the supplement by running the corresponding matlab files:

    TableTS1.m
    TableTS2.m
    TableTS3.m
    TableTS6.m
    TableTS7.m
    TableTS8a.m
    TableTS8b.m
    TableTS8c.m
    TableTS8d.m