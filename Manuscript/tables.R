# table contents
concepts_ml <-
  data.frame(
    Concept = c(
      "Sample units",
      "Cluster",
      "Hierarchical data",
      "Level-1",
      "Level-2",
      "Hierarchical model",
      "Fixed effect",
      "Random effect",
      "Mixed effect",
      "ICC",
      "Stratified intercept"
    ),
    Details = c(
      "Units of the population from which measurements are taken in a sample, e.g., students.",
      "Variable that specify the cluster or agruppation, e.g., Classroom",
      "Data are grouped into clusters at different levels, observations belonging to the same cluster are expected to share certain characteristics.",
      "Variable that varies within a cluster, eg. Test score",
      "Variable that does not vary within a cluster but between, e.g. teacher experience.",
      "Model accounting for dependant observations relying on certain parameters ( within cluster) which in turn depend on other parameters (between cluster)" ,
      "Effects that are constant across all sample units, e.g. something that researchers control for and can repeat, such as  a teaching strategy (tutoring after class)",
      "Effects that are a source of random variation in the data, and whose levels are not fully sampled. e.g. test score tendency during academic year between students due to no controlled factors such as  genetic,family history",
      "Includes fixed and random effects, e.g. the fixed effect would be the treatment effect of a drug and the random effect would be the ID of the hospital where the patient is treated. Multilevel models typically accommodate for variability by including a separate group mean for each cluster e.g random intercept on hospitals. In addition to random intercepts, multilevel models can also include random coefficients and heterogeneous residual error variances across clusters (see e.g. @gelm06, @hox17 and @jong21).",
      "The variability due to clustering is often measured by means of the intraclass coefficient (ICC). The ICC can be seen as the percentage
of variance that can be attributed to the cluster-level, where a high ICC would indicate that a lot of variability is due to the cluster structure.",
""
    )
  )