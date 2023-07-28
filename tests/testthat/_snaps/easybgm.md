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
        intrusion-dreams    0.527                 1.000          Inf     included
         intrusion-flash    0.563                 0.987       77.125     included
            dreams-flash    0.314                 0.960       23.876     included
         intrusion-upset    0.010                 0.075        0.081     excluded
            dreams-upset    0.051                 0.232        0.302 inconclusive
             flash-upset    0.177                 0.607        1.546 inconclusive
       intrusion-physior    0.115                 0.448        0.811 inconclusive
          dreams-physior    0.002                 0.030        0.031     excluded
           flash-physior    0.049                 0.215        0.273 inconclusive
           upset-physior    0.844                 1.000          Inf     included
      
       Bayes Factors larger than 10 were considered sufficient evidence for the categorization. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of included edges: 4 
       Number of inconclusive edges: 4 
       Number of excluded edges: 2 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 82 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.2425 
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
      
       Bayes Factors larger than 10 were considered sufficient evidence for the categorization. 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of included edges: 10 
       Number of inconclusive edges: 0 
       Number of excluded edges: 0 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 72 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 1 
      ---

