# easybgm works for bgms

    Code
      summary(res_bgms)
    Output
      
       BAYESIAN ANALYSIS OF NETWORKS 
       Model type: ordinal 
       Number of nodes: 5 
       Fitting Package: bgms 
      --- 
       EDGE SPECIFIC OVERVIEW 
                Relation Estimate Posterior Incl. Prob. Inclusion BF     Category
        intrusion-dreams    0.635                 1.000          Inf     included
         intrusion-flash    0.712                 1.000          Inf     included
            dreams-flash    0.405                 1.000          Inf     included
         intrusion-upset    0.000                 0.029        0.030     excluded
            dreams-upset    0.064                 0.285        0.399 inconclusive
             flash-upset    0.354                 0.883        7.547 inconclusive
       intrusion-physior    0.337                 0.968       30.250     included
          dreams-physior    0.003                 0.028        0.029     excluded
           flash-physior    0.031                 0.152        0.179 inconclusive
           upset-physior    0.977                 1.000          Inf     included
      
       Bayes Factors larger than 10 were considered sufficient evidence for the classification 
       Bayes factors were obtained using Bayesian model-averaging. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of edges with sufficient evidence for inclusion: 5 
       Number of edges with insufficient evidence: 3 
       Number of edges with sufficient evidence for exclusion: 2 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 21 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.508 
      ---

# easybgm works for bdgraph

    Code
      summary(res_bdgraph)
    Output
      
       BAYESIAN ANALYSIS OF NETWORKS 
       Model type: ggm 
       Number of nodes: 5 
       Fitting Package: bdgraph 
      --- 
       EDGE SPECIFIC OVERVIEW 
                Relation Estimate Posterior Incl. Prob. Inclusion BF     Category
        intrusion-dreams    0.379                  1.00          Inf     included
         intrusion-flash    0.439                  1.00          Inf     included
            dreams-flash    0.262                  1.00          Inf     included
         intrusion-upset    0.277                  1.00          Inf     included
            dreams-upset    0.182                  1.00          Inf     included
             flash-upset    0.020                  0.21        0.266 inconclusive
       intrusion-physior    0.000                  0.00        0.000     excluded
          dreams-physior    0.013                  0.27        0.370 inconclusive
           flash-physior    0.285                  1.00          Inf     included
           upset-physior    0.595                  1.00          Inf     included
      
       Bayes Factors larger than 10 were considered sufficient evidence for the classification 
       Bayes factors were obtained using Bayesian model-averaging. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of edges with sufficient evidence for inclusion: 7 
       Number of edges with insufficient evidence: 2 
       Number of edges with sufficient evidence for exclusion: 1 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 4 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.592 
      ---

# easybgm works for sparse vs dense

    Code
      sparse_dense
    Output
      $log.BF
      [1] -1.443074
      
      $BF
      [1] 0.2362006
      
      $relative.complexity.sparse
      [1] 0.6344
      
      $relative.complexity.dense
      [1] 0.7019
      
      $relative.complexity.uniform
      [1] 0.633
      
      $no.hypotheses
      [1] 3
      

