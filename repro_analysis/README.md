
## Global analysis of iceplant flowering

This sub-folder holds the code, models, and manuscript outline/notes for
the global analysis of iceplant reproductive effort. Pretty much nothing
will work “off the shelf” to start, as most data isn’t ready to be made
public.

You can get the gist of what’s going on in the `R/` folder. The data are
not public yet, but that will change eventually. `R/00_main.R` runs the
project, and `R/01_depends.R` creates the correct environment for
subsequent scripts.

### Scripts

All contained in `R/` sub-folder:

1.  `01_depends.R`: loads dependencies and functions.

2.  `02_gbif_cleaning.R`: loads and cleans a Darwin Core Archive of all
    \_Carpobrotus spp.\*\_ occurrence data from GBIF. This archive is
    available with DOI: <https://doi.org/10.15468/dl.q8kmsy>. Due to the
    size, the bulk of the operations are commented out, and intermediate
    results are cached as csv files. Subsequent scripts work with the
    csv outputs, though un-commenting would allow one to re-run the
    complete pipeline using only the Darwin Core Archive.

3.  `03_clim_envelopes.R`: Creates Figure 1 A-D. These are plots use
    Chelsa data to generate a long-run climate envelope, and then using
    local weather data to look at actual experienced climate in the
    years immediately prior to sampling (when flowering decisions were
    likely to have been made). The idea is to see what “normal” years
    look like vs what our sampling actually captured. The GBIF data are
    used estimate the range of climate the species experience globally,
    and show the coverage of our sampling.

4.  `04_null_models.R`: This generates models for `repro ~ log_size`
    (Bernoulli) and `flower_n ~ log_size` (0-truncated Poisson). These
    are used as base models to compare against the climate models.

5.  `05_flower_models.R`: Contains code that can fit the
    `repro/flower_n ~ log_size + clim_vars` models. This code is
    technically runnable locally, but will take forever. See
    `R/slidingwin/r_scripts/shrinkage.R` and
    `R/slidingwin/submit_scripts/submit_shrinkage.sh` for cluster
    runnable code. Additionally, generates model summaries, and MCMC
    diagnostic plots + posterior predictive checks for all, which are
    stored in `model_summaries/` and `Manuscript/Figures/` folders,
    respectively.

6.  `06_model_comparison.R`: Contains code to run WAIC + LOO-CV model
    comparisons. Again, this code is technically runnable locally, but
    will take forever. See `R/slidingwin/r_scripts/shrink_compare.R` and
    `R/slidingwin/submit_scripts/submit_cv.sh` for cluster runnable
    code.

7.  `07_figures.R`: Creates figures for the publication’s main body and
    appendix.
