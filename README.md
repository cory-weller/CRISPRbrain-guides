# CRISPRbrain-path

## Data retrieval
Retrievd the following tables and placed them in the `data/` directory:

[Glutamatergic Neuron-Day-14-Survival-CRISPRi](data/Tian_et_al_2019_2.csv) from [here](https://www.crisprbrain.org/simple-screen/?screen=Glutamatergic%20Neuron-Day-14-Survival-CRISPRi)

[Glutamatergic Neuron-Survival-CRISPRi](data/Tian_et_al_2020_15.csv) from [here](https://www.crisprbrain.org/simple-screen/?screen=Glutamatergic%20Neuron-Survival-CRISPRi)

Supplementary files [1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc1.xlsx)
[2](https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc2.xlsx)
and [3](https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc3.xlsx)
from [here](https://doi.org/10.1016/j.cell.2022.05.013)



## Pathway analysis
Ran as a single script; see [`pathwaysAnalysis.R`](src/pathwaysAnalysis.R) for details

```bash
module load R/4.2
Rscript src/pathwaysAnalysis.R
```

Producing the subsequent two outputs:
* [Glutamatergic Neuron-Day-14-Survival-CRISPRi](outputs/pathway_enrichment_2_all.tsv)
* [Glutamatergic Neuron-Survival-CRISPRi](outputs/pathway_enrichment_15_all.tsv)

## Guide RNA scoring for clusters
[`guideScoring.R`](src/guideScoring.R) imports the supplemental excel
file tables and combines them into a single summary file after doing
some minor filtering. See [`guideScoring.R`](src/guideScoring.R)

```bash
module load R/4.2
Rscript src/guideScoring.R
```



## Pathway and iNDI guide selection
See [`select_Pathway_iNDI_guides.R`](src/select_Pathway_iNDI_guides.R) 

```bash
module load R/4.2
Rscript src/select_Pathway_iNDI_guides.R
```
