
# Statistical approaches to study exposome-health associations in the context of repeated exposure data: a simulation study

by
Warembourg Charline,
Anguita-Ruiz Augusto,
Valérie Siroux,
Rémy Slama,
Vrijheid Martine,
Richiardi Lorenzo,
Basagaña Xavier

> This is an online repository gathering all simulated data and codes employed for the analyses presented in the paper (Warembourg et al. (2023)).
> The purpose of this repository is to allow researchers to reproduce all findings presented in the paper, as well as to adapt provided codes for the
> analyses of their own repeated exposome datasets.

> This paper has been submitted for publication in *Environmental Sciences & Technologies*.

> In the current paper, we conduct a simulation analysis to compare the performance of different exposome analysis techniques in the scenario 
> of repeated exposure data. The different tested approaches followed during our simulation analysis are illustrated in the next figure:

![](manuscript/figures/overview_methods.png)

*Figure. Overview of the statistical methods tested to estimate exposome-health associations in the context of repeated exposure.*


## Abstract

> The exposome concept has emerged quite recently and aims at considering all the environmental stressors simultaneously. The high-dimension of the data and the potential correlations that may exist 
> between exposures lead to various statistical challenges. Some methodological papers have reported recommendations on how to analyze exposome data. However, few studies tried to characterize the
> longitudinal relationship between repeated measures of the exposome and a health outcome. Here, we conduct a simulation study to compare the performance of different statistical approaches to 
> assess exposome-health associations in the context of multiple and repeated exposure variables. Two scenarios were tested: 1) assuming that all time points are associated with Y, and 2) assuming
> that only a single time point is associated with Y. An application study using real data collected within the INMA mother-child cohort (Spain) is also presented. Some methods such as sPLS and DSA 
> showed better performance than other methods in both scenarios. However, none of the tested methods provided good enough performances when the number of true predictors increased. The results of
> the simulation study call for the development of new statistical methods or approaches that are able to address both the issue of multiple (and correlated) variables and of repeated exposome data.


## Software implementation

> The software written to produce the results presented in this paper is organized in three main blocks (dataset generation, data analysis and results interpretation).

All source code used to generate the results and figures in the paper are in
the `src` folder.
The data used in this study is provided in `data` and the sources for the
manuscript text and figures are in `manuscript`.
See the `README.md` files in each directory for a full description.


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone ...

or [download a zip archive](https://github.com/.../archive/master.zip).


## Dependencies

You'll need a working R environment to run the code.
The required libraries for each analysis are specified at the beginning of each script.


## License

All source code is made available under a BSD 3-clause license. You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors. See `LICENSE.md` for the full license text.

The manuscript text is not open source. The authors reserve the rights to the
article content, which is currently submitted for publication in the
*Environmental Sciences & Technologies*.




