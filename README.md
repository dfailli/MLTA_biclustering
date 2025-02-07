## A novel approach for biclustering bipartite networks: an extension of finite mixtures of latent trait analyzers

we extend the mixture of latent trait analyzers (MLTA) model developed by Gollini and Murphy (2014) and Gollini (2020) to perform a joint clustering of the two disjoint sets of nodes of a bipartite network,
as in the biclustering framework. In detail, sending nodes are partitioned into clusters (called components) via a finite mixture of latent trait models. In each component, receiving nodes are partitioned into clusters (called segments) by adopting a flexible and parsimonious specification of the linear predictor. Residual dependence between receiving nodes is modeled via a multidimensional latent trait, as in the original MLTA specification. Furthermore, by incorporating nodal attributes into the model’s latent layer, we gain insight into how these attributes impact the formation of components. To estimate model parameters, an EM-type algorithm based on a Gauss-Hermite approximation of intractable integrals is proposed.

#### Files
* `Code.R` contains the code for:
  * estimating model parameters;
  * performing the simulation study;
  * analyzing data on pediatric patients affected by appendicitis with both the proposed biclustering approach and the original MLTA specification.

* `appl_pediatric2.RData` contains data and the R objects for the application.

* `logistic.cpp` and `logistic2.cpp` contain the C++ code to speed up the estimation algoithm.

#### Main References
* Gollini, I., & Murphy, T. B. (2014). Mixture of latent trait analyzers for model-based clustering of categorical data. Statistics and Computing, 24 (4), 569–588.

* Gollini, I. (2020). A mixture model approach for clustering bipartite networks. Challenges in Social Network Research Volume in the Lecture Notes in Social Networks (LNSN - Series of Springer).

* Failli, D., Marino, M. F., & Martella, F. (2024). Finite Mixtures of Latent Trait Analyzers With Concomitant Variables for Bipartite Networks: An Analysis of COVID-19 Data. Multivariate Behavioral Research, 59 (4), 801–817. https://doi.org/10.1080/00273171.2024.2335391

* Failli, D., Arpino, B., & Marino, M.F. (2024). A finite mixture approach for the analysis of digital skills in Bulgaria, Finland and Italy: the role of socioeconomic factors. Statistical Methods & Applications, 33, 1483–1511. https://doi.org/10.1007/s10260-024-00766-w

#### Data
* https://doi.org/10.5281/zenodo.7669442

<!---
* Failli, D., Marino, M.F., Martella, F. (2022) Extending finite mixtures of latent trait analyzers for bipartite networks. In Balzanella A., Bini M., Cavicchia C. and Verde R. (Eds.) Book of short Paper SIS 2022 (pp. 540-550), Pearson.
-->
