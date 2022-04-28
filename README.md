SwedishMusselAdaptivePotential
================
Mathias Wegner, AWI
31/03/2022

# SwedishMusselAdaptivePotential

This markdown file covers the analyses used in Ventura A., Wegner K.M.,
Dupont S. “Assessing adaptation potential to ocean acidification in blue
mussels, Mytilus edulis, from the Swedish west coast”

## Read Data

## Larval size

This fit larval size as a function of seawater acidification levels
(i.e. pH = 8.1 or pH = 7.5) as fixed effects and “animal” (relatedness
between individual larvae), “dam” (female parent), “block” (experimental
block) and “culture” (replicate larval culture within each treatment and
family) as random intercept effects. To test for family specific
reaction norms between the different seawater acidification levels we
also fit “animal” as a random slope and did the same for “dam”. We chose
the best fitting model based on the respective Deviance Information
Criterion (DIC). Only undeformed animals were used.

#test for normality

``` 
![](README_files/figure-gfm/LengthModels%20normality-1.png)<!-- -->

``` r
#within pH 7.5
norm.test_7.5 = ks.test(goodSize$Length[goodSize$pH == 7.5], 'pnorm', 
                        mean = mean(goodSize$Length[goodSize$pH == 7.5]),
                                    sd= sd(goodSize$Length[goodSize$pH == 7.5]))

```

    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  goodSize$Length[goodSize$pH == 7.5]
    ## D = 0.056987, p-value = 0.5851
    ## alternative hypothesis: two-sided

``` r
#within ph 8.1
norm.test_8.1 = ks.test(goodSize$Length[goodSize$pH == 8.1], 'pnorm', 
                        mean = mean(goodSize$Length[goodSize$pH == 8.1]),
                                    sd= sd(goodSize$Length[goodSize$pH == 8.1]))


```

    ## 
    ##  One-sample Kolmogorov-Smirnov test
    ## 
    ## data:  goodSize$Length[goodSize$pH == 8.1]
    ## D = 0.045697, p-value = 0.05793
    ## alternative hypothesis: two-sided

``` r
#remove deformed animals
goodSize=subset(size,size$Malformed == 'no')

LengthModels = makeMCMC_gaussian_with_pH(goodSize)

#give Model Fit
for(i in 1:length(LengthModels))
{
    print('######################################################################################################
')
    print(paste(c('Model:',i)))
    print(summary(LengthModels[[i]]))
}
```

    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "1"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1528.65 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.1946  0.04602   0.3613    2.585
    ## 
    ##                ~Dam
    ## 
    ##     post.mean  l-95% CI u-95% CI eff.samp
    ## Dam   0.03631 0.0003421   0.1503    130.1
    ## 
    ##                ~Block
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## Block   0.03776    3e-04   0.2055    153.1
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.04448 0.0006207   0.2636      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1978   0.1078   0.2743    2.353
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.896   16.547   17.337      180 <0.006 **
    ## pH8.1           2.134    2.039    2.205      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "2"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1579.544 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal   0.28855  0.10892   0.4874    7.224
    ## pH8.1:pH7.5.animal   0.04185 -0.01663   0.2278    4.657
    ## pH7.5:pH8.1.animal   0.04185 -0.01663   0.2278    4.657
    ## pH8.1:pH8.1.animal   0.06183  0.00368   0.2472    2.823
    ## 
    ##                ~Dam
    ## 
    ##     post.mean  l-95% CI u-95% CI eff.samp
    ## Dam    0.0584 0.0005334   0.1543    95.13
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.03179 0.0003042   0.1045      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture     0.159 0.0005176   0.1014      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.2349   0.1354   0.2797    3.983
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.807   16.310   17.235    140.8 <0.006 **
    ## pH8.1           2.218    1.925    2.520    180.0 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "3"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1564.081 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean  l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal    0.3332 0.1451559   0.4804   13.182
    ## pH8.1:pH7.5.animal    0.0778 0.0007085   0.2076   15.437
    ## pH7.5:pH8.1.animal    0.0778 0.0007085   0.2076   15.437
    ## pH8.1:pH8.1.animal    0.1094 0.0298107   0.2226    4.229
    ## 
    ##                ~us(pH):Dam
    ## 
    ##                 post.mean   l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.Dam   0.09099  0.0002886   0.3279    180.0
    ## pH8.1:pH7.5.Dam   0.02390 -0.0353395   0.1747    109.4
    ## pH7.5:pH8.1.Dam   0.02390 -0.0353395   0.1747    109.4
    ## pH8.1:pH8.1.Dam   0.05039  0.0014133   0.2049    180.0
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.03311 0.0001972   0.1149      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.04006 0.0004397   0.1795    148.1
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.2089   0.1401    0.263    4.627
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.784   16.432   17.198      180 <0.006 **
    ## pH8.1           2.220    1.916    2.518      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "4"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1186.936 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.327  0.09176   0.5173    3.432
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.04635 0.0003313   0.2124    62.35
    ## 
    ##                ~Culture
    ## 
    ##         post.mean l-95% CI u-95% CI eff.samp
    ## Culture   0.02227  0.00052   0.1194    26.57
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1335  0.03902   0.2512    3.787
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.858   16.456   17.262    75.81 <0.006 **
    ## pH8.1           2.136    2.042    2.215   180.00 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "5"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1529.421 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal    0.3573  0.22706   0.5178   14.063
    ## pH8.1:pH7.5.animal    0.1290 -0.02805   0.2861    4.555
    ## pH7.5:pH8.1.animal    0.1290 -0.02805   0.2861    4.555
    ## pH8.1:pH8.1.animal    0.1589  0.07973   0.2762    7.715
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.0279 0.0004165   0.1099      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.02941 0.0008771  0.08714      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1868   0.1306   0.2361    6.854
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.818   16.272   17.149      141 <0.006 **
    ## pH8.1           2.210    1.938    2.527      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "6"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1448.591 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.1939  0.04798   0.4587    2.207
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1148 0.0003271   0.1754      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.2036  0.07037   0.2904    2.054
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.865   16.641   17.148      180 <0.006 **
    ## pH8.1           2.129    2.041    2.227      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "7"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1578.779 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.1823  0.06357   0.3692     4.15
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.2076   0.1077   0.2692    5.275
    ## 
    ##  Location effects: Length ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    16.894   16.684   17.135    142.4 <0.006 **
    ## pH8.1           2.131    2.050    2.201    180.0 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#calculate models for Length as a response within each pH environment

lowpH = subset(goodSize, goodSize$pH == 7.5)

LengthModels_7.5 = makeMCMC_gaussian(lowpH)

#give Model Fit
for(i in 1:4)
{
    print('######################################################################################################')
    print(paste(c('Model:',i)))
    print(summary(LengthModels_7.5[[i]]))
}
```

    ## [1] "######################################################################################################"
    ## [1] "Model:" "1"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 420.4362 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.338  0.05486    1.182    2.825
    ## 
    ##                ~Dam
    ## 
    ##     post.mean l-95% CI u-95% CI eff.samp
    ## Dam    0.1556  0.00113    0.488    122.8
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1849 0.0006853   0.4995      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.05399 0.0003213   0.2898    106.7
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.5624  0.02658   0.7763    3.852
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     24.10    23.69    24.78      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "2"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 376.5695 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.7544   0.2552    1.305    11.67
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1964 0.0002669   0.7403      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.05764 0.0005065   0.3772      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.3485   0.0436   0.6074    11.21
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     24.12    23.52    24.85    77.62 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "3"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: -11.17553 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     1.213    0.119    1.613    5.707
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1419 0.0003607   0.7099     41.6
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.1157 0.007464   0.6059    6.831
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     24.04    23.38    24.63      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "4"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 377.342 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.7299   0.2153    1.405    7.424
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.3746  0.04518   0.6893    6.982
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     24.10    23.74    24.60      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#now for 8.1
highpH = subset(goodSize, goodSize$pH == 8.1)

LengthModels_8.1 = makeMCMC_gaussian(highpH)

#give Model Fit
for(i in 1:4)
{
    print('######################################################################################################')
    print(paste(c('Model:',i)))
    print(summary(LengthModels_8.1[[i]]))
}
```

    ## [1] "######################################################################################################"
    ## [1] "Model:" "1"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 2210.963 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.3521  0.08007   0.7826     3.35
    ## 
    ##                ~Dam
    ## 
    ##     post.mean  l-95% CI u-95% CI eff.samp
    ## Dam   0.09314 0.0003272   0.3514    13.69
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.08092 0.0004492   0.3403    74.86
    ## 
    ##                ~Culture
    ## 
    ##         post.mean l-95% CI u-95% CI eff.samp
    ## Culture    0.1139 0.002173   0.4025      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.6579    0.422   0.8445    2.287
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     33.77    33.22    34.31    142.9 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "2"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 2206.141 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.4465   0.1471   0.8458    16.17
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.06335 0.0003477   0.3532      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean l-95% CI u-95% CI eff.samp
    ## Culture   0.09975 0.001233   0.4305      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.6122   0.4078   0.7986    13.41
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     33.79    33.38    34.38      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "3"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 2195.46 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal      0.53    0.256   0.9017    7.263
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.09288 0.0002356   0.1471      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.5855   0.3943   0.7576    9.659
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     33.72    33.29    34.11      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################"
    ## [1] "Model:" "4"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 2109.667 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal    0.6541   0.1905    1.133    3.618
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units    0.5218   0.2619   0.7624    3.642
    ## 
    ##  Location effects: Length ~ 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     33.73    33.40    34.21      180 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
best_RS_Model = LengthModels[[5]]
best_75_Model = LengthModels_7.5[[2]]
best_81_Model = LengthModels_8.1[[2]]



autocorr.diag(best_RS_Model$Sol)
```

    ##          (Intercept)       pH8.1
    ## Lag 0    1.000000000  1.00000000
    ## Lag 5    0.118651215  0.05072907
    ## Lag 25  -0.122255864 -0.06539395
    ## Lag 50  -0.003409542  0.04565720
    ## Lag 250 -0.045264022  0.02853719

``` r
autocorr.diag(best_RS_Model$VCV)
```

    ##         pH7.5:pH7.5.animal pH8.1:pH7.5.animal pH7.5:pH8.1.animal
    ## Lag 0           1.00000000          1.0000000          1.0000000
    ## Lag 5           0.85430808          0.9385564          0.9385564
    ## Lag 25          0.54772898          0.7460908          0.7460908
    ## Lag 50          0.33899211          0.4809673          0.4809673
    ## Lag 250        -0.02749104         -0.1095597         -0.1095597
    ##         pH8.1:pH8.1.animal       Block      Culture      units
    ## Lag 0            1.0000000  1.00000000  1.000000000  1.0000000
    ## Lag 5            0.9173545  0.08664581 -0.029850513  0.8376214
    ## Lag 25           0.6614807 -0.05136609 -0.005332093  0.6122965
    ## Lag 50           0.4776442 -0.08742722 -0.031833655  0.4639455
    ## Lag 250         -0.1688314 -0.00833030  0.006719733 -0.1352499

``` r
#test convergence
heidel.diag(best_RS_Model$VCV)
```

    ##                                                  
    ##                    Stationarity start     p-value
    ##                    test         iteration        
    ## pH7.5:pH7.5.animal passed       1         0.3848 
    ## pH8.1:pH7.5.animal passed       1         0.0529 
    ## pH7.5:pH8.1.animal passed       1         0.0529 
    ## pH8.1:pH8.1.animal passed       1         0.3050 
    ## Block              passed       1         0.5132 
    ## Culture            passed       1         0.9512 
    ## units              passed       1         0.7042 
    ##                                              
    ##                    Halfwidth Mean   Halfwidth
    ##                    test                      
    ## pH7.5:pH7.5.animal failed    0.3573 0.04378  
    ## pH8.1:pH7.5.animal failed    0.1290 0.06645  
    ## pH7.5:pH8.1.animal failed    0.1290 0.06645  
    ## pH8.1:pH8.1.animal failed    0.1589 0.04197  
    ## Block              failed    0.0279 0.00912  
    ## Culture            failed    0.0294 0.01466  
    ## units              failed    0.1868 0.02297

``` r
h2_RS_75 = best_RS_Model$VCV[,'pH7.5:pH7.5.animal']/(best_RS_Model$VCV[,'pH7.5:pH7.5.animal']+best_RS_Model$VCV[,'pH8.1:pH7.5.animal']+best_RS_Model$VCV[,'pH8.1:pH8.1.animal']+best_RS_Model$VCV[,'Block']+best_RS_Model$VCV[,'Culture']+best_RS_Model$VCV[,'units'])
h2_RS_81 = best_RS_Model$VCV[,'pH8.1:pH8.1.animal']/(best_RS_Model$VCV[,'pH7.5:pH7.5.animal']+best_RS_Model$VCV[,'pH8.1:pH7.5.animal']+best_RS_Model$VCV[,'pH8.1:pH8.1.animal']+best_RS_Model$VCV[,'Block']+best_RS_Model$VCV[,'Culture']+best_RS_Model$VCV[,'units'])
h2_75 = best_75_Model$VCV[,'animal']/(best_75_Model$VCV[,'animal']+best_75_Model$VCV[,'Block']+best_75_Model$VCV[,'Culture']+best_75_Model$VCV[,'units'])
h2_81 = best_81_Model$VCV[,'animal']/(best_81_Model$VCV[,'animal']+best_81_Model$VCV[,'Block']+best_81_Model$VCV[,'Culture']+best_81_Model$VCV[,'units'])

posterior.mode(h2_RS_75)
```

    ##      var1 
    ## 0.4032025

``` r
posterior.mode(h2_RS_81)
```

    ##     var1 
    ## 0.133973

``` r
posterior.mode(h2_75)
```

    ##      var1 
    ## 0.6142119

``` r
posterior.mode(h2_81)
```

    ##      var1 
    ## 0.3278893

``` r
HPDinterval(h2_RS_75)
```

    ##          lower     upper
    ## var1 0.2784111 0.5554507
    ## attr(,"Probability")
    ## [1] 0.95

``` r
HPDinterval(h2_RS_81)
```

    ##           lower     upper
    ## var1 0.09301704 0.2918709
    ## attr(,"Probability")
    ## [1] 0.95

``` r
HPDinterval(h2_75)
```

    ##          lower     upper
    ## var1 0.2365988 0.9618591
    ## attr(,"Probability")
    ## [1] 0.95

``` r
HPDinterval(h2_81)
```

    ##          lower     upper
    ## var1 0.1527041 0.6261698
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'pH7.5:pH7.5.animal'])
```

    ##     var1 
    ## 0.299529

``` r
HPDinterval(best_RS_Model$VCV[,'pH7.5:pH7.5.animal'])
```

    ##          lower     upper
    ## var1 0.2270633 0.5177564
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'pH8.1:pH8.1.animal'])
```

    ##      var1 
    ## 0.1014164

``` r
HPDinterval(best_RS_Model$VCV[,'pH8.1:pH8.1.animal'])
```

    ##           lower     upper
    ## var1 0.07973273 0.2762089
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'pH7.5:pH8.1.animal'])
```

    ##      var1 
    ## 0.1469689

``` r
HPDinterval(best_RS_Model$VCV[,'pH7.5:pH8.1.animal'])
```

    ##            lower     upper
    ## var1 -0.02804693 0.2861256
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'Culture'])
```

    ##        var1 
    ## 0.004216391

``` r
HPDinterval(best_RS_Model$VCV[,'Culture'])
```

    ##             lower      upper
    ## var1 0.0008770832 0.08714077
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'Block'])
```

    ##        var1 
    ## 0.002075086

``` r
HPDinterval(best_RS_Model$VCV[,'Block'])
```

    ##             lower     upper
    ## var1 0.0004164715 0.1098782
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_RS_Model$VCV[,'units'])
```

    ##      var1 
    ## 0.2096505

``` r
HPDinterval(best_RS_Model$VCV[,'units'])
```

    ##         lower     upper
    ## var1 0.130577 0.2360593
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_75_Model$VCV[,'animal'])
```

    ##      var1 
    ## 0.4344261

``` r
HPDinterval(best_75_Model$VCV[,'animal'])
```

    ##         lower    upper
    ## var1 0.255193 1.304904
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_75_Model$VCV[,'Culture'])
```

    ##        var1 
    ## 0.003022448

``` r
HPDinterval(best_75_Model$VCV[,'Culture'])
```

    ##             lower     upper
    ## var1 0.0005064541 0.3771989
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_75_Model$VCV[,'Block'])
```

    ##       var1 
    ## 0.01508136

``` r
HPDinterval(best_75_Model$VCV[,'Block'])
```

    ##             lower     upper
    ## var1 0.0002669293 0.7403393
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_75_Model$VCV[,'units'])
```

    ##      var1 
    ## 0.3043381

``` r
HPDinterval(best_75_Model$VCV[,'units'])
```

    ##           lower     upper
    ## var1 0.04360211 0.6074497
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_81_Model$VCV[,'animal'])
```

    ##      var1 
    ## 0.3740626

``` r
HPDinterval(best_81_Model$VCV[,'animal'])
```

    ##          lower     upper
    ## var1 0.1470633 0.8457979
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_81_Model$VCV[,'Culture'])
```

    ##        var1 
    ## 0.007039994

``` r
HPDinterval(best_81_Model$VCV[,'Culture'])
```

    ##            lower     upper
    ## var1 0.001233452 0.4305328
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_81_Model$VCV[,'Block'])
```

    ##        var1 
    ## 0.001777947

``` r
HPDinterval(best_81_Model$VCV[,'Block'])
```

    ##             lower     upper
    ## var1 0.0003477375 0.3531913
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(best_81_Model$VCV[,'units'])
```

    ##      var1 
    ## 0.6314963

``` r
HPDinterval(best_81_Model$VCV[,'units'])
```

    ##          lower     upper
    ## var1 0.4078148 0.7986012
    ## attr(,"Probability")
    ## [1] 0.95

``` r
plot(best_RS_Model$VCV[,'pH7.5:pH7.5.animal'], main = 'animal RS', ylim = c(0,5))
```

![](/README_files/figure-gfm/LengthModels-1.png)<!-- -->

``` r
plot(best_81_Model$VCV[,'animal'], main = 'animal')
```

![](/README_files/figure-gfm/LengthModels-2.png)<!-- -->

## Larval Deformation

This fits a bionmial model analysing deformation as a function of pH,
and “animal” (relatedness between individual larvae), “dam” (female
parent), “block” (experimental block) and “culture” (replicate larval
culture within each treatment and family) as random intercept effects.

``` r
DeformModels = makeMCMC_binomial(size)
#give Model Fit
for(i in 1:length(DeformModels))
{
    print('######################################################################################################
')
    print(paste(c('Model:',i)))
    print(summary(DeformModels[[i]]))
}
```

    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "1"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1516.681 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.255   0.1468   0.3972    5.956
    ## 
    ##                ~Block
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## Block    0.5304  0.00047    1.454      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.04924 0.0002396     0.16      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     1.880    1.195    2.588   255.68 <0.006 **
    ## pH8.1          -4.401   -4.607   -4.157    27.52 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "2"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1488.434 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal    0.6304   0.3522   1.0723    7.753
    ## pH8.1:pH7.5.animal    0.2217  -0.3604   0.5103    6.180
    ## pH7.5:pH8.1.animal    0.2217  -0.3604   0.5103    6.180
    ## pH8.1:pH8.1.animal    1.1132   0.5163   2.0212    5.779
    ## 
    ##                ~Block
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## Block    0.9403 0.000327    3.238      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean l-95% CI u-95% CI eff.samp
    ## Culture   0.04171 0.000283   0.1957    110.9
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    1.8853   0.9716   2.9048   180.00 0.0111 * 
    ## pH8.1         -4.6093  -5.3233  -3.8994    41.64 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "3"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1504.794 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.585   0.3404   0.8279    7.772
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.7433 0.0004412    4.906    102.2
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     1.888    1.164    2.857   180.00 0.0111 * 
    ## pH8.1          -4.502   -4.794   -4.064    14.54 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "4"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1502.462 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.645   0.2642   0.9228    4.698
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     1.830    1.434    2.295   143.78 <0.006 **
    ## pH8.1          -4.504   -4.819   -4.224    22.07 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "5"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1507.369 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal   0.09956  0.02317    0.191    3.361
    ## 
    ##                ~Dam
    ## 
    ##     post.mean l-95% CI u-95% CI eff.samp
    ## Dam    0.4164  0.02015    1.084    144.1
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.2964 0.0004733    1.416    107.5
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture    0.2727 0.0002963   0.3997    116.5
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)     1.837    1.148    2.564   180.00 <0.006 **
    ## pH8.1          -4.460   -4.752   -4.150    20.48 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "6"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1488.705 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal    1.5609  0.50424   4.5768    2.705
    ## pH8.1:pH7.5.animal    0.5205  0.03454   1.6646    2.418
    ## pH7.5:pH8.1.animal    0.5205  0.03454   1.6646    2.418
    ## pH8.1:pH8.1.animal    0.2434  0.08055   0.6111    2.651
    ## 
    ##                ~Dam
    ## 
    ##     post.mean  l-95% CI u-95% CI eff.samp
    ## Dam    0.3348 0.0005971    1.157    49.78
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.4342 0.0003123    2.248      180
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture    0.1976 0.0002402   0.4704      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    1.9740   0.9568   3.2036   139.96 0.0222 * 
    ## pH8.1         -4.6075  -5.1964  -4.1109    16.68 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "7"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1493.341 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal    0.4553  0.23722   0.6806   13.082
    ## pH8.1:pH7.5.animal    0.1308  0.01455   0.3076    1.999
    ## pH7.5:pH8.1.animal    0.1308  0.01455   0.3076    1.999
    ## pH8.1:pH8.1.animal    0.1998  0.08865   0.3601    4.485
    ## 
    ##                ~us(pH):Dam
    ## 
    ##                 post.mean  l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.Dam   0.55022  0.001564   1.8135   115.72
    ## pH8.1:pH7.5.Dam  -0.02236 -0.487270   0.8765    18.69
    ## pH7.5:pH8.1.Dam  -0.02236 -0.487270   0.8765    18.69
    ## pH8.1:pH8.1.Dam   0.35776  0.008599   1.0338   122.11
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.3973 0.0002758    1.868    60.21
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture   0.03311 0.0003308   0.1128    122.2
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    1.8685   0.9097   2.9956   131.93 0.0111 * 
    ## pH8.1         -4.4694  -5.2347  -3.7317    79.17 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "8"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 1495.12 
    ## 
    ##  G-structure:  ~us(pH):Dam
    ## 
    ##                 post.mean  l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.Dam    1.8375  0.263719   4.3642    31.28
    ## pH8.1:pH7.5.Dam    0.5111 -0.163550   1.5496    22.46
    ## pH7.5:pH8.1.Dam    0.5111 -0.163550   1.5496    22.46
    ## pH8.1:pH8.1.Dam    0.1965  0.001148   0.6479    53.75
    ## 
    ##                ~animal
    ## 
    ##        post.mean l-95% CI u-95% CI eff.samp
    ## animal     0.119  0.04024   0.3333    8.709
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.4152 0.0003093    1.818    105.5
    ## 
    ##                ~Culture
    ## 
    ##         post.mean  l-95% CI u-95% CI eff.samp
    ## Culture    0.0982 0.0004929   0.2645      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: Malformed ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC   
    ## (Intercept)    1.9039   0.1516   3.0858    180.0 0.0111 * 
    ## pH8.1         -4.4125  -5.5866  -3.6517    116.3 <0.006 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
bestDeformedModel = DeformModels[[7]]


h2_def_75 = bestDeformedModel$VCV[,'pH7.5:pH7.5.animal']/(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH8.1.animal']+bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam']+bestDeformedModel$VCV[,'Block']+bestDeformedModel$VCV[,'Culture']+bestDeformedModel$VCV[,'units']+1)
h2_def_81 = bestDeformedModel$VCV[,'pH8.1:pH8.1.animal']/(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH8.1.animal']+bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam']+bestDeformedModel$VCV[,'Block']+bestDeformedModel$VCV[,'Culture']+bestDeformedModel$VCV[,'units']+1)


d2_def_75 = bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam']/(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH8.1.animal']+bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam']+bestDeformedModel$VCV[,'Block']+bestDeformedModel$VCV[,'Culture']+bestDeformedModel$VCV[,'units']+1)
d2_def_81 = bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam']/(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH7.5.animal']+bestDeformedModel$VCV[,'pH8.1:pH8.1.animal']+bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH7.5.Dam']+bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam']+bestDeformedModel$VCV[,'Block']+bestDeformedModel$VCV[,'Culture']+bestDeformedModel$VCV[,'units']+1)


posterior.mode(h2_def_75)
```

    ##     var1 
    ## 0.117951

``` r
posterior.mode(h2_def_81)
```

    ##       var1 
    ## 0.03643999

``` r
HPDinterval(h2_def_75)
```

    ##           lower   upper
    ## var1 0.03674129 0.18445
    ## attr(,"Probability")
    ## [1] 0.95

``` r
HPDinterval(h2_def_81)
```

    ##           lower      upper
    ## var1 0.01864063 0.09315467
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(d2_def_75)
```

    ##       var1 
    ## 0.01778438

``` r
posterior.mode(d2_def_81)
```

    ##       var1 
    ## 0.01955617

``` r
HPDinterval(d2_def_75)
```

    ##             lower     upper
    ## var1 0.0005490607 0.3636999
    ## attr(,"Probability")
    ## [1] 0.95

``` r
HPDinterval(d2_def_81)
```

    ##            lower     upper
    ## var1 0.001996691 0.2511025
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal'])
```

    ##      var1 
    ## 0.4506529

``` r
HPDinterval(bestDeformedModel$VCV[,'pH7.5:pH7.5.animal'])
```

    ##          lower     upper
    ## var1 0.2372228 0.6805943
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH8.1:pH8.1.animal'])
```

    ##      var1 
    ## 0.1968155

``` r
HPDinterval(bestDeformedModel$VCV[,'pH8.1:pH8.1.animal'])
```

    ##           lower     upper
    ## var1 0.08864516 0.3601258
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH7.5:pH8.1.animal'])
```

    ##       var1 
    ## 0.06649425

``` r
HPDinterval(bestDeformedModel$VCV[,'pH7.5:pH8.1.animal'])
```

    ##           lower     upper
    ## var1 0.01454756 0.3075689
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam'])
```

    ##       var1 
    ## 0.05090321

``` r
HPDinterval(bestDeformedModel$VCV[,'pH7.5:pH7.5.Dam'])
```

    ##            lower   upper
    ## var1 0.001563884 1.81347
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam'])
```

    ##       var1 
    ## 0.09117575

``` r
HPDinterval(bestDeformedModel$VCV[,'pH8.1:pH8.1.Dam'])
```

    ##            lower    upper
    ## var1 0.008598537 1.033834
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'pH7.5:pH8.1.Dam'])
```

    ##       var1 
    ## -0.2854262

``` r
HPDinterval(bestDeformedModel$VCV[,'pH7.5:pH8.1.Dam'])
```

    ##           lower     upper
    ## var1 -0.4872704 0.8765235
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'Culture'])
```

    ##        var1 
    ## 0.001821042

``` r
HPDinterval(bestDeformedModel$VCV[,'Culture'])
```

    ##             lower     upper
    ## var1 0.0003307814 0.1127737
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'Block'])
```

    ##       var1 
    ## 0.01022916

``` r
HPDinterval(bestDeformedModel$VCV[,'Block'])
```

    ##             lower   upper
    ## var1 0.0002757604 1.86759
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestDeformedModel$VCV[,'units'])
```

    ##     var1 
    ## 0.999813

``` r
HPDinterval(bestDeformedModel$VCV[,'units'])
```

    ##      lower upper
    ## var1     1     1
    ## attr(,"Probability")
    ## [1] 0.95

## Surivival

This fits a bionmial model analysing deformation as a function of pH,
and “animal” (relatedness between individual larvae), “dam” (female
parent), “block” (experimental block) and “culture” (replicate larval
culture within each treatment and family) as random intercept effects.

``` r
survMCMC = makeMCMC_binomial_surv(surv)
for(i in 1:length(survMCMC))
{
    print('######################################################################################################
')
    print(paste(c('Model:',i)))
    print(summary(survMCMC[[i]]))
}
```

    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "1"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30929.86 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean  l-95% CI u-95% CI eff.samp
    ## animal   0.01439 0.0003536  0.06097    5.427
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1375 0.0004972   0.5834      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC  
    ## (Intercept)   0.54684 -0.03515  1.11668      180 0.0444 *
    ## pH8.1         0.03546 -0.43857  0.43173      180 0.8444  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "2"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30930.79 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean   l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal 0.0411398  0.0004492  0.12289    22.01
    ## pH8.1:pH7.5.animal 0.0005458 -0.0610986  0.05823    11.35
    ## pH7.5:pH8.1.animal 0.0005458 -0.0610986  0.05823    11.35
    ## pH8.1:pH8.1.animal 0.0248974  0.0007067  0.09457    20.45
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1491 0.0006002    0.451      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC  
    ## (Intercept)    0.5171   0.1211   0.9323      180 0.0444 *
    ## pH8.1          0.0294  -0.4192   0.4562      180 0.8667  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "3"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30929.31 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean  l-95% CI u-95% CI eff.samp
    ## animal   0.02048 0.0008053  0.07799    4.901
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean  l-95% CI  u-95% CI eff.samp  pMCMC   
    ## (Intercept)  0.542508  0.277400  0.780991      180 <0.006 **
    ## pH8.1       -0.001928 -0.448470  0.427681      180  0.956   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "4"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30928.02 
    ## 
    ##  G-structure:  ~animal
    ## 
    ##        post.mean  l-95% CI u-95% CI eff.samp
    ## animal   0.02145 0.0002032  0.06422       17
    ## 
    ##                ~Dam
    ## 
    ##     post.mean l-95% CI u-95% CI eff.samp
    ## Dam    0.2794 0.001054   0.8563    82.27
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block   0.08977 0.0004361   0.4035      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC  
    ## (Intercept)   0.54079  0.05349  1.02827      180 0.0444 *
    ## pH8.1         0.01942 -0.39161  0.41212      180 0.9778  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "5"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30933.86 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean  l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal   0.02218  0.001695  0.07833   25.908
    ## pH8.1:pH7.5.animal   0.01689 -0.038574  0.06732    5.176
    ## pH7.5:pH8.1.animal   0.01689 -0.038574  0.06732    5.176
    ## pH8.1:pH8.1.animal   0.05140  0.001078  0.15095    4.638
    ## 
    ##                ~Dam
    ## 
    ##     post.mean l-95% CI u-95% CI eff.samp
    ## Dam    0.3181 0.007243   0.8721      180
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block    0.1703 0.0003692   0.2831      180
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC  
    ## (Intercept)   0.54445  0.02958  1.11573      180 0.0444 *
    ## pH8.1         0.04771 -0.37614  0.38652      180 0.8333  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "######################################################################################################\n"
    ## [1] "Model:" "6"     
    ## 
    ##  Iterations = 101:996
    ##  Thinning interval  = 5
    ##  Sample size  = 180 
    ## 
    ##  DIC: 30927 
    ## 
    ##  G-structure:  ~us(pH):animal
    ## 
    ##                    post.mean   l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.animal  0.033167  0.0006928  0.13859    7.879
    ## pH8.1:pH7.5.animal -0.004494 -0.0580952  0.02527    6.193
    ## pH7.5:pH8.1.animal -0.004494 -0.0580952  0.02527    6.193
    ## pH8.1:pH8.1.animal  0.019063  0.0008423  0.06526   12.916
    ## 
    ##                ~us(pH):Dam
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp
    ## pH7.5:pH7.5.Dam    0.6076 0.050575    1.752   180.00
    ## pH8.1:pH7.5.Dam    0.4011 0.008378    1.043   180.00
    ## pH7.5:pH8.1.Dam    0.4011 0.008378    1.043   180.00
    ## pH8.1:pH8.1.Dam    0.3377 0.006456    1.049    96.79
    ## 
    ##                ~Block
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## Block     0.113 0.0003576   0.5722    118.3
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean l-95% CI u-95% CI eff.samp
    ## units         1        1        1        0
    ## 
    ##  Location effects: cbind(Live, Dead) ~ pH 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC  
    ## (Intercept)   0.54978 -0.06394  1.28511      180 0.0778 .
    ## pH8.1         0.04534 -0.35954  0.58301      180 0.8000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
bestSurvMod = survMCMC[[4]]
plot(bestSurvMod$VCV)
```

![](/README_files/figure-gfm/survival-1.png)<!-- -->

``` r
autocorr.diag(bestSurvMod$Sol)
```

    ##         (Intercept)       pH8.1
    ## Lag 0    1.00000000  1.00000000
    ## Lag 5    0.02064267 -0.05242241
    ## Lag 25   0.15781063  0.10372951
    ## Lag 50   0.02384271 -0.04558470
    ## Lag 250  0.02863564  0.06861999

``` r
autocorr.diag(bestSurvMod$VCV)
```

    ##             animal         Dam        Block units
    ## Lag 0    1.0000000  1.00000000  1.000000000   NaN
    ## Lag 5    0.8265101  0.21363920  0.101310525   NaN
    ## Lag 25   0.4044246 -0.08967832 -0.035933445   NaN
    ## Lag 50   0.2773647 -0.06041102 -0.005232913   NaN
    ## Lag 250 -0.0678969  0.12083559 -0.012460618   NaN

``` r
#test convergence
heidel.diag(bestSurvMod$VCV)
```

    ##                                      
    ##        Stationarity start     p-value
    ##        test         iteration        
    ## animal passed        1        0.500  
    ## Dam    passed        1        0.517  
    ## Block  passed        1        0.249  
    ## units  failed       NA           NA  
    ##                                  
    ##        Halfwidth Mean   Halfwidth
    ##        test                      
    ## animal failed    0.0214 0.0116   
    ## Dam    failed    0.2794 0.0562   
    ## Block  failed    0.0898 0.0376   
    ## units  <NA>          NA     NA

``` r
h2_surv = bestSurvMod$VCV[,'animal']/(bestSurvMod$VCV[,'animal']+bestSurvMod$VCV[,'Dam']+bestSurvMod$VCV[,'Block']+bestSurvMod$VCV[,'units']+1)
posterior.mode(h2_surv)
```

    ##         var1 
    ## 0.0004467648

``` r
HPDinterval(h2_surv)
```

    ##             lower      upper
    ## var1 9.719049e-05 0.02930371
    ## attr(,"Probability")
    ## [1] 0.95

``` r
d2_surv = bestSurvMod$VCV[,'Dam']/(bestSurvMod$VCV[,'animal']+bestSurvMod$VCV[,'Dam']+bestSurvMod$VCV[,'Block']+bestSurvMod$VCV[,'units']+1)
posterior.mode(d2_surv)
```

    ##       var1 
    ## 0.03912818

``` r
HPDinterval(d2_surv)
```

    ##             lower     upper
    ## var1 0.0005224318 0.2954442
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestSurvMod$VCV[,'animal'])
```

    ##         var1 
    ## 0.0008026866

``` r
HPDinterval(bestSurvMod$VCV[,'animal'])
```

    ##             lower      upper
    ## var1 0.0002031989 0.06422089
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestSurvMod$VCV[,'Dam'])
```

    ##       var1 
    ## 0.08187572

``` r
HPDinterval(bestSurvMod$VCV[,'Dam'])
```

    ##            lower     upper
    ## var1 0.001053669 0.8562731
    ## attr(,"Probability")
    ## [1] 0.95

``` r
posterior.mode(bestSurvMod$VCV[,'Block'])
```

    ##        var1 
    ## 0.002452804

``` r
HPDinterval(bestSurvMod$VCV[,'Block'])
```

    ##             lower     upper
    ## var1 0.0004361478 0.4035403
    ## attr(,"Probability")
    ## [1] 0.95


