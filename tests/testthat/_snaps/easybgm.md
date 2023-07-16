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
        intrusion-dreams    0.503                 1.000          Inf     included
         intrusion-flash    0.575                 1.000          Inf     included
            dreams-flash    0.237                 0.797        3.916 inconclusive
         intrusion-upset    0.037                 0.208        0.263 inconclusive
            dreams-upset    0.141                 0.572        1.335 inconclusive
             flash-upset    0.133                 0.504        1.015 inconclusive
       intrusion-physior    0.092                 0.394        0.649 inconclusive
          dreams-physior    0.010                 0.116        0.132 inconclusive
           flash-physior    0.138                 0.559        1.270 inconclusive
           upset-physior    0.775                 1.000          Inf     included
      
       Bayes Factors larger 10 were considered sufficient evidence for the categorization. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of included edges: 3 
       Number of inconclusive edges: 7 
       Number of excluded edges: 0 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 115 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.099 
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
                Relation Estimate Posterior Incl. Prob. Inclusion BF Category
        intrusion-dreams    0.346                     1          Inf included
         intrusion-flash    0.274                     1          Inf included
            dreams-flash    0.408                     1          Inf included
         intrusion-upset    0.140                     1          Inf included
            dreams-upset    0.283                     1          Inf included
             flash-upset    0.179                     1          Inf included
       intrusion-physior    0.242                     1          Inf included
          dreams-physior   -0.183                     1          Inf included
           flash-physior    0.222                     1          Inf included
           upset-physior    0.542                     1          Inf included
      
       Bayes Factors larger 10 were considered sufficient evidence for the categorization. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of included edges: 10 
       Number of inconclusive edges: 0 
       Number of excluded edges: 0 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 0 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 1 
      ---

