# 3D_contact_distribution

This R script calculates the Hi-C contact distribution for a subset of genes of interest (e.g., RXLRs or some other category of genes).  The code also calculates Monte Carlo statistics to determine whether the 3D contact distribution of the subset genes differs statistically from the distribution of an equivalent set of background genes randomly chosen from the entire gene set.


DEPENDENCIES

*3D_contact_distribution* requires command line access to [grep](https://www.gnu.org/software/grep/manual/grep.html), which will already be installed by default on most UNIX and macOS systems.  *3D_contact_distribution* was developed with grep v2.6.0-FreeBSD, but uses a basic function that should be compatible with most modern versions.

The script also requires the installation of one R package: data.table (v1.14.2).  All R packages can easily be installed from [CRAN](https://cran.r-project.org) within R.


INPUTS

*3D_contact_distribution* requires three input files: a GFF3 file of all genes in the genome ('all_gene_file' – the background gene set), a GFF3 file of the subset of genes being studied ('subset_gene_file' – the genes-of-interest gene set), and a sparse matrix of 3D contacts between regularly sized windows along the genome as calculated from Hi-C data ('example_matrix.dat'). Sparse matrices can be created from HIC files using [strawr](https://cran.r-project.org/web/packages/strawr/index.html).

Example files are included in this distribution.  Note that due to GitHub file size upload limitations, the sparse matrix example file provided only contains contacts for chromosome 1.  Contacts for genes outside chromosome 1 will thus return as zero.  For proper usage, a full sparse matrix for the entire genome should be used.

The script also requires four global variables: the window ('bin') size used to create the sparse matrix of 3D contacts; the maximum number of pairwise contacts to consider in the gene-of-interest subset and Monte Carlo simulated gene set; the number of iterations to be used for the Monte Carlo simulations (typically 10<sup>4</sup> to 10<sup>5</sup>); and the summary statistic for comparison (the default is 'mean' with alternative 'median').  The 'max.number' setting is needed because pairwise comparisons increase exponentially as each new gene is added; for large gene subsets, it may be necessary to limit the number of pairwise comparisons being considered.

Exemplars for all these inputs are shown in the worked example below.


USAGE

Program usage, from the UNIX command line, is as follows:

```
3D_contact_distribution.R
```

The script can also be run directly in an R console.


EXAMPLE

*3D_contact_distribution* produces one output: the probabilities that the mean (or median) number of contacts between the subset genes is greater than or less than that expected for the same number of randomly chosen genes (calculated via a bootstrap with replacement).  Here is an example with the 'median' setting:

```
Greater than 193 : p = 0
Less than 193 : p = 1
```

Note that in this example the subset genes are all adjacent along a chromosome, and therefore have a large number of contacts with each other (a median of 193).  In comparisons, randomly chosen genes across the genome typically have 0-1 contacts.  Hence the observed distribution for this subset of genes is notably different from the random distribution.  This may differ, likely quite markedly, for other subsets of genes.  Because most gene pairs have zero contacts, the 'median' setting is dominated by these zero cases and often washes out other signals.  The 'mean' setting is often much more informative.


RUNTIME SPEED

Note that the Monte Carlo simulation is *slow*. 

Profiling shows that ~99% of runtime is used to identify required rows in the sparse matrix.  In the 'data.table' package, this function is already heavily optimized, so cannot easily be improved upon.  It seems that the only feasible option to improve runtime is to reduce the size of the sparse matrix, by using larger window bins.  This effect is exponential; doubling the window size reduces the overall matrix size fourfold. However, because the data form is a sparse matrix (where zero contact windows are inferred but not explicitly written as a record), this improvement will diminish as the window size increases and more windows include at least one contact.  

In some scenarios, changing the 'max.number' setting may also decrease runtime. However, setting this variable too low will negatively affect statistical accuracy.
