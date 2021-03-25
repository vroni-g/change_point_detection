Notes
-----

Here we can keep track of what we are doing in the repo.

Currently I am using R version 4.0.2

To apply the methods we use the R package PerMuTe as the base code. Some
adaptations need to be done for improving the algorithm.

Note: Package implements an apply fn to a 3d array that is somehow
slower than an apply fn to a 2D array (what I curently do for the global
greening paper).

TODO: Get apply fn to ignore pixels with all NA.

TODO: Make the apply fn run in parallel - although for my own work I do

``` r
library(clustermq)
  results<- Q(tmp_fn, 
              i=1:NPERM, 
              const = list(perm_matrix = perm_matrix, 
                           my_data = my_data),
              export = list(ALPHA = ALPHA), 
              n_jobs = 100, 
              template = list(job_name = "perm", 
                              partition = "all", 
                              log_file = "perm.txt",
                              memory = 10000,
                              n_cpus = 1),
              fail_on_error = FALSE, 
              verbose = TRUE)
```

Because the above can spread the job to all servers, as opposed to
running the process in one server.

Introduction
------------

Methods
-------

Recent comparison of change point algorithms:

<a href="https://arxiv.org/pdf/2003.06222.pdf" class="uri">https://arxiv.org/pdf/2003.06222.pdf</a>

The above shows binary segmentation is not statistically different from
the others (see critical difference, CD, in Figure 4). Below we see a
paper that improves binary segementation and has an R package
`changepoint` for implementation:

<a href="https://www.researchgate.net/publication/48180788_Optimal_Detection_of_Changepoints_With_a_Linear_Computational_Cost" class="uri">https://www.researchgate.net/publication/48180788_Optimal_Detection_of_Changepoints_With_a_Linear_Computational_Cost</a>

\#\#\#Other sections…
