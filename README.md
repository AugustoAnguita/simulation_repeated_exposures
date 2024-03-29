
# Statistical approaches to study exposome-health associations in the context of repeated exposure data: a simulation study

by
Charline Warembourg,
Augusto Anguita-Ruiz,
Valérie Siroux,
Rémy Slama,
Martine Vrijheid,
Richiardi Lorenzo,
Basagaña Xavier

> This is an online repository gathering all simulated data and codes employed for the analyses presented in the paper (Warembourg et al. (2023)).
> The purpose of this repository is to allow researchers to reproduce all findings presented in the paper, as well as to adapt provided codes for the
> analyses of their own repeated exposome datasets.

> This paper has been submitted for publication in *Environmental Sciences & Technologies*.

> In the current paper, we conduct a simulation analysis to compare the performance of different exposome analysis techniques in the scenario 
> of repeated exposure data. The different tested approaches followed during our simulation analysis are illustrated in the next figure:


![](https://github.com/AugustoAnguita/simulation_repeated_exposures/blob/main/src/image/overview_methods.png)
*Figure. Overview of the statistical methods tested to estimate exposome-health associations in the context of repeated exposure.*


## Abstract

> The exposome concept aims at considering all environmental stressors simultaneously. The dimension of the data and the correlation that may exist between exposures lead to various statistical challenges. Some methodological studies have provided insight regarding the efficiency of specific modeling approaches in the context of exposome data assessed once in each subject. However, few studies considered the situation in which environmental exposures are assessed repeatedly. Here, we conduct a simulation study to compare the performance of statistical approaches to assess exposome-health associations in the context of multiple and repeated exposure variables. Different scenarios were tested assuming different types and numbers of exposure-outcome causal relationships. An application study using real data collected within the INMA mother-child cohort (Spain) is also presented. In the tested scenarios, methods such as sparse partial least squares (sPLS) and Deletion-Substitution-Addition algorithm (DSA) showed better performance than other methods. Performance of all methods decreased when the number of true predictors increased. Our results highlight that having repeated exposome data poses extra challenges to detecting true associations, and that care must be taken when choosing the statistical methods to use. Besides, results should be interpreted with caution, especially in contexts with limited sample size, given the elevated chance of reporting false positive or negative associations.


## Simulation results

![](https://github.com/AugustoAnguita/simulation_repeated_exposures/blob/main/src/image/figure_exp3_s1ands2.png)
*Figure 1. Scatter plot representing the sensitivity (x-axis) and the false discovery rate (y-axis) obtained from each method under Scenario 1 (all time points are associated with Y) and Scenario 2 (a single time point is associated with Y) for **K=3 number of true exposures**. Labels are in black if the performances correspond to a method performed using raw data and in red if the performances correspond to a method performed using averaged data (2-step). The green square represents the area under which the sensitivity is above 75% and the false discovery rate is below 25%. Abbreviations: Bon; Bonferroni multiple-testing correction; BY; Benjamini-Yekutieli multiple-testing correction (false discovery rate); DSA, Deletion-Substitution-Addition algorithm; ExWAS, Exposome-wide association study; ENET, Elastic net; FDR, False discovery rate; none, no multiple-testing correction is applied; sPLS, sparse partial least squares.*

![](https://github.com/AugustoAnguita/simulation_repeated_exposures/blob/main/src/image/figure_exp5_s1ands2.png)
*Figure 2. Scatter plot representing the sensitivity (x-axis) and the false discovery rate (y-axis) obtained from each method under Scenario 1 (all time points are associated with Y) and Scenario 2 (a single time point is associated with Y) for **K=5 number of true exposures**. Labels are in black if the performances correspond to a method performed using raw data and in red if the performances correspond to a method performed using averaged data (2-step). The green square represents the area under which the sensitivity is above 75% and the false discovery rate is below 25%. Abbreviations: Bon; Bonferroni multiple-testing correction; BY; Benjamini-Yekutieli multiple-testing correction (false discovery rate); DSA, Deletion-Substitution-Addition algorithm; ExWAS, Exposome-wide association study; ENET, Elastic net; FDR, False discovery rate; none, no multiple-testing correction is applied; sPLS, sparse partial least squares.*

![](https://github.com/AugustoAnguita/simulation_repeated_exposures/blob/main/src/image/figure_exp10_s1ands2.png)
*Figure 3. Scatter plot representing the sensitivity (x-axis) and the false discovery rate (y-axis) obtained from each method under Scenario 1 (all time points are associated with Y) and Scenario 2 (a single time point is associated with Y) for **K=10 number of true exposures**. Labels are in black if the performances correspond to a method performed using raw data and in red if the performances correspond to a method performed using averaged data (2-step). The green square represents the area under which the sensitivity is above 75% and the false discovery rate is below 25%. Abbreviations: Bon; Bonferroni multiple-testing correction; BY; Benjamini-Yekutieli multiple-testing correction (false discovery rate); DSA, Deletion-Substitution-Addition algorithm; ExWAS, Exposome-wide association study; ENET, Elastic net; FDR, False discovery rate; none, no multiple-testing correction is applied; sPLS, sparse partial least squares.*


## Software implementation

> The software written to produce the results presented in this paper is organized in three main blocks (dataset generation, data analysis and results interpretation).

All source code used to generate the data and the results in the paper are in
the `src` folder.
The simulated data used in this study is provided in `data`.
The selection of variables obtained after the application of each method are provided in `results` (sub-divided by scenario and number of true predictors)
See the `README.md` files in each directory for a full description.


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/AugustoAnguita/simulation_repeated_exposures

or [download a zip archive](https://github.com/AugustoAnguita/simulation_repeated_exposures/archive/refs/heads/main.zip).


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




