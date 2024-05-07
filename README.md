### Sensitivity Analysis for Causal Effects Under Missing Binary Data: Application to the Obesity Paradox

In this repository you can find the code necessary to reproduce the results from the paper ``Sensitivity Analysis for Causal Effects Under Missing Binary Data: Application to the Obesity Paradox''.


#### Reproducing the Results

The code used for reproducing the results of the paper is contained in the `scripts/` folder. In the below table, we point to the files used to generate the respective figures.

| Figure   | Code   | 
|:-----------|:------------:|
| [Fig. 1](#): Causal diagrams for proximal inference. | |
| [Fig. 2](nonpar-agnostic.png): Non-parametric inference of agnostic $\phi$-values. | `scripts/vignettes/nonpar.Rmd` <br> (chunk `r nonpar-agnostic`)
| [Fig. 3](nonpar-x.png/nonpar-y.png): Non-parametric inference of $x$- and $y$-specific $\phi$-values. | `scripts/vignettes/nonpar.Rmd` <br> (chunks `r x-grid`, `r y-grid`)
| [Fig. 4](nonpar-x-axis-ci.png): Non-parametric $x$-specific $\phi$-values with uncertainty quantification. | `scripts/vignettes/nonpar.Rmd` <br> (chunk `r axis-search-ci`)
| [Fig. 5](#): Inference of $\phi$-values for the parametric exponential family model. | `scripts/vignettes/par.Rmd` <br> (chunk `r fi-search`)
| [Fig. 6](#): Parametric $x$-specific $\phi$-values with uncertainty quantification. | `scripts/vignettes/par.Rmd` <br> (chunk `r x-axis-search`)
| [Fig. 7](#): Investigating First-Order & Gradient Stability on MIMIC-IV data. | `scripts/fos-gs-cond.R`
| [Fig. 8](#): Sensitivity of EM to initial parameter values. | `scripts/init-sense.R`

