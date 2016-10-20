-   [[My directory of R packages for data analysis and visualization, Bayesian statistics, mapping, GIS, climate, Species Distribution Modelling, ecology, biogeography, evolution, reproducible science...](#my-directory-of-r-packages-for-data-analysis-and-visualization-bayesian-statistics-mapping-gis-climate-species-distribution-modelling-ecology-biogeography-evolution-reproducible-science...)](#[my-directory-of-r-packages-for-data-analysis-and-visualization,-bayesian-statistics,-mapping,-gis,-climate,-species-distribution-modelling,-ecology,-biogeography,-evolution,-reproducible-science...](#my-directory-of-r-packages-for-data-analysis-and-visualization-bayesian-statistics-mapping-gis-climate-species-distribution-modelling-ecology-biogeography-evolution-reproducible-science...))
-   [[Data sources](#data-sources)](#[data-sources](#data-sources))
-   [[Data management](#data-management)](#[data-management](#data-management))
-   [[General statistics](#general-statistics)](#[general-statistics](#general-statistics))
    -   [[Ordination & Multivariate Analysis](#ordination-multivariate-analysis)](#[ordination-&-multivariate-analysis](#ordination-multivariate-analysis))
    -   [[Survival analysis](#survival-analysis)](#[survival-analysis](#survival-analysis))
    -   [[Regression tools](#regression-tools)](#[regression-tools](#regression-tools))
    -   [[Bayesian/MCMC](#bayesianmcmc)](#[bayesian/mcmc](#bayesianmcmc))
-   [[Plotting & Visualisation](#plotting-visualisation)](#[plotting-&-visualisation](#plotting-visualisation))
    -   [[GGPLOT2](#ggplot2)](#[ggplot2](#ggplot2))
    -   [[Colour](#colour)](#[colour](#colour))
-   [[Reproducible Research - Report generation](#reproducible-research---report-generation)](#[reproducible-research---report-generation](#reproducible-research---report-generation))
-   [[Parallelisation & Big Data](#parallelisation-big-data)](#[parallelisation-&-big-data](#parallelisation-big-data))
-   [[Niche & Species Distribution Modelling](#niche-species-distribution-modelling)](#[niche-&-species-distribution-modelling](#niche-species-distribution-modelling))
    -   [[Occupancy modelling](#occupancy-modelling)](#[occupancy-modelling](#occupancy-modelling))
-   [[Climate](#climate)](#[climate](#climate))
-   [[GIS/spatial functionality](#gisspatial-functionality)](#[gis/spatial-functionality](#gisspatial-functionality))
-   [[Spatial Analysis](#spatial-analysis)](#[spatial-analysis](#spatial-analysis))
-   [[Networks](#networks)](#[networks](#networks))
-   [[Phylogenetics, phylogeography & comparative analysis](#phylogenetics-phylogeography-comparative-analysis)](#[phylogenetics,-phylogeography-&-comparative-analysis](#phylogenetics-phylogeography-comparative-analysis))
-   [[Palaeoecology](#palaeoecology)](#[palaeoecology](#palaeoecology))
-   [[Dendrochronology](#dendrochronology)](#[dendrochronology](#dendrochronology))
-   [[Ecological analyses (miscellaneous)](#ecological-analyses-miscellaneous)](#[ecological-analyses-(miscellaneous)](#ecological-analyses-miscellaneous))
-   [[Miscellaneous](#miscellaneous)](#[miscellaneous](#miscellaneous))
-   [[R programming](#r-programming)](#[r-programming](#r-programming))

My directory of R packages for data analysis and visualization, Bayesian statistics, mapping, GIS, climate, Species Distribution Modelling, ecology, biogeography, evolution, reproducible science...
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

[F. Rodriguez-Sanchez](http://bit.ly/frod_san)

Updated 2016-10-20

These are packages that I often use or, alternatively, I need only rarely but don't want to forget about. Of course, there are many other useful packages out there (e.g. at [CRAN](http://cran.r-project.org/web/packages/available_packages_by_name.html) or [GitHub](https://github.com/search?q=package&type=Repositories&ref=advsearch&l=R). Check also [CRAN task views](http://cran.r-project.org/web/views/).

Data sources
------------

Check [CRAN Task View on Web technologies](http://cran.r-project.org/web/views/WebTechnologies.html).

Data management
---------------

Check out this great cheatsheet: [Data wrangling with dplyr and tidyr](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)

-   [rio]() A Swiss-Army Knife for Data I/O
-   [tidyr]() Reshape data
-   [reshape2]() Reshape data from long to wide format and vice versa
-   [dplyr]() Data wrangling
-   [data.table]() Tools for managing data frames
-   [gdata]() Data manipulation
-   [Hmisc]() Frank Harrell's miscellaneous tools
-   [Kmisc]() Data reshaping, table and plot generation from RMarkdown
-   [lubridate]() Dates and times
-   [multitable]() Manipulate multiple arrays
-   [datacheck]() Tools for checking data consistency
-   [tabplot]() Large datasets viz.
-   [tabplotd3]() Interactive
-   [tableplot]()
-   [taRifx](http://cran.r-project.org/web/packages/taRifx/index.html) Useful functions
-   [DescTools](http://cran.r-project.org/web/packages/DescTools/index.html) Many useful functions
-   [summarytools]() Quickly summarize dataframes (inc. markdown output)

General statistics
------------------

### Ordination & Multivariate Analysis

Check also [Environmetrics](http://cran.r-project.org/web/views/Environmetrics.html) and [Multivariate](http://cran.r-project.org/web/views/Multivariate.html) CRAN Task Views.

-   [vegan]()
-   [ade4](http://cran.r-project.org/web/packages/ade4/index.html) Multivariate data analysis and display
-   [dave]() Data Analysis in Vegetation Ecology
-   [ecodist]() Dissimilarity-based analysis (ordination, Mantel tests...)
-   [labdsv]() Includes plotting functions
-   [mvabund]() Analysing multivariate data (upscaling from individual species models)
-   [boral]() Bayesian ordination and regression analysis
-   [eigenprcomp](http://cran.r-project.org/web/packages/eigenprcomp/) Computes confidence intervals for the proportion explained by the first 1,2,k principal components

### Survival analysis

-   [survival]()
-   [survMisc]()
-   [coxme]() Mixed-effects cox regression
-   [rawr]() Plotting functions: kmplot & ggsurv

### Regression tools

-   [aod]() Analysis of Overdispersed data
-   [AICcmodavg](http://cran.r-project.org/web/packages/AICcmodavg/index.html) Model selection and multimodel inference
-   [arm]() Gelman's package: includes bayesglm, sim, coefplot...
-   [bbmle]() Tools for MLE
-   [binomTools]() Diagnostics for binomial regression
-   [car]() Regression tools
-   [coefplot]() Plot model coefficients
-   [confReg]() Estimating confidence of individual regression predictions
-   [COUNT]() Regression models for count data (Poisson, Negative Binomial)
-   [DAAG]() Some regression tools from the book 'Data Analysis and Graphics using R'
-   [dhglm]() Hierarchical GLMs with random effects in both the mean and dispersion components
-   [dynlm]() Dynamic linear models and time series regression
-   [earth]() Multivariate Adaptive Regression Splines
-   [effects]() Displays effects estimated from regression models
-   [emdbook]() Tools from 'Ecological Models and Data' (inc. several plotting functions)
-   [fit.models]() Compare results from different models
-   [FME]() Modelling tools
-   [gam]() Generalised Additive Models (GAMs)
-   [gamm4]() GAMMs
-   [gbm]() Generalised Boosted regression models
-   [gee]() Generalized Estimating Equations
-   [glm2]() To fit GLMs with convergence problems
-   [MuMIn]() Model selection and averaging
-   [glmulti]() Model selection and multimodel inference
-   [gnm]() Generalized nonlinear models
-   [heatmapFit]() Checking logistic regression goodness of fit
-   [hier.part]() Variance partitioning to assess 'importance' of predictors. See also `relaimpo`.
-   [relaimpo]() Relative importance of predictors
-   [nlme]() Mixed models
-   [lme4]() Mixed models
-   [pamm]() Power analysis for mixed models
-   [odprism]() Power analysis for mixed models
-   [likelihood]() Maximum Likelihood Estimation
-   [lmtest]() Diagnostic checks for linear regression
-   [MARSS]() Multivariate Autoregressive State-Space Modeling
-   [MCMCglmm]() Mixed models fitted by MCMC
-   [mgcv]() GAM fitting
-   [miscF]() Multivariate Normal regression, spatial Bayesian mixed models, piecewise regression...
-   [mlogit]() Multinomial regression
-   [mvinfluence]() Influence Measures and Diagnostic Plots for Multivariate Linear Models
-   [msm]() Multi-state Markov and hidden Markov models in continuous time
-   [PCovR]() Principal Covariates Regression (Reducing collinear predictor variables to a few components and regressing on them)
-   [HLMdiag]() Diagnostic tools for hierarchical models (fitted with lme4)
-   [stremo]() Learning Structural Equation Models
-   [rockchalk]() Regression Estimation and Presentation
-   [R2STATS]() GUI for fitting and comparing GLM and GLMM in R
-   [LMERConvenienceFunctions]()
-   [R2admb]()
-   [rms]() Regression modeling strategies
-   [Zelig]()

### Bayesian/MCMC

Check also [Bayesian CRAN Task View](http://cran.r-project.org/web/views/Bayesian.html).

-   [blme]() Bayesian multilevel models
-   [BMA]() Model averaging
-   [dclone]() Data cloning and MCMC tools (inc. JAGS functions)
-   [dcmle]() Hierarchical models through data cloning
-   [dlm]() Dynamic models and MCMC tools
-   [dlmodeler]() State-space modelling
-   [sspir]() State-space models
-   [MCMCglmm]() Mixed models fitted by MCMC
-   [RSGHB]() Hierarchical Bayes models
-   [bayespref]() Analysis of count data
-   [predcomps](https://github.com/dchudz/predcomps) Average Predictive Comparisons

#### Generic MCMC samplers

-   [adaptMCMC](http://cran.r-project.org/web/packages/adaptMCMC/index.html) Generic MCMC sampler
-   [BRugs]() Interface to OpenBUGS
-   [BayesX]() Structured Additive Regression
-   [R2BayesX]()
-   [rjags]()
-   [runjags]()
-   [R2jags]()
-   [jagsUI]()
-   [R2OpenBUGS]()
-   [R2WinBUGS]()
-   [rstan]()
-   [rstanarm]()
-   [filzbach]()
-   [iBUGS]() Interface to BUGS/JAGS
-   [rube](http://stat.cmu.edu/~hseltman/rube/)
-   [INLA]() Integrated Nested Laplace Approximation
-   [LaplacesDemon]()
-   [LearnBayes]()
-   [MCMCpack]() MCMC samplers
-   [MHAdative]() MCMC sampler
-   [glmmBUGS]()
-   [mcmc]()
-   [datalist](http://cran.r-project.org/web/packages/datalist/index.html)

#### MCMC diagnostics

-   [coda]() MCMC output analysis
-   [boa]() Analyses of MCMC output (like `coda`)
-   [ggmcmc]() Graphic analysis of MCMC output
-   [mcmcmplots](https://cran.r-project.org/package=mcmcplots) Plot MCMC output
-   [plotMCMC]() Diagnostic plots
-   [bmk]() MCMC diagnostics
-   [superdiag]() testing MCMC noncovergence
-   [shinyStan]()
-   [MCMCvis]()
-   [rwty](https://cran.r-project.org/package=rwty)
-   [dMCMC](https://github.com/MarcoDVisser/dMCMC)

#### Spatial Bayes

-   [CARBayes]() Spatial models with CAR
-   [spBayes]() spatial Bayes
-   [geoCount](http://cran.r-project.org/web/packages/geoCount/index.html) generalized linear spatial models

#### ABC

-   [abc](http://cran.r-project.org/web/packages/abc/index.html) Approximate Bayesian Computation
-   [easyABC]() Approximate Bayesian Computation

Plotting & Visualisation
------------------------

Check also [CRAN Task View on Graphics](http://cran.r-project.org/web/views/Graphics.html).

-   [denstrip]() Plotting distributions (w uncertainty)
-   [visualize]() Graph Probability Distributions
-   [fanplot]() Visualise sequential probability distributions
-   [diagram]() Networks, flow diagrams, etc
-   [beanplot]() Bean plots
-   [vioplot]() Violin plots
-   [viopoints]() 1-D Scatter Plots with Jitter Using Kernel Density Estimates
-   [dagR]() Directed Acyclic graphs (DAGs)
-   [effects]() Displays effects estimated from regression models
-   [coefplot]() Plot model coefficients
-   [plotmo]() Plot model responses
-   [corrgram]() Plot correlation matrix
-   [gclus]() Clustering graphics
-   [vcd]() Viz categorical data
-   [extracat]() Viz categorical data
-   [plot3d](http://cran.r-project.org/web/packages/plot3D/index.html) Plotting multi-dimensional data
-   [VizOR](http://cran.r-project.org/web/packages/VizOR/index.html) visualization tools for complex observational data
-   [GrapheR]() GUI for base plots

### GGPLOT2

-   [ggplot2]()
-   [GGally]() extension to ggplot
-   [ggsubplot]() Embedding subplots within plots
-   [ggthemes]() Extra themes, scales and geoms for ggplot
-   [gplot2bdc](https://github.com/briandconnelly/ggplot2bdc)
-   [ggtern](http://www.ggtern.com/) Ternary diagrams
-   [ggthemr](https://github.com/cttobin/ggthemr) Extra themes for ggplot2

-   [plotly](https://github.com/ropensci/plotly) R interface to plotly
-   [googleVis]() Using google chart tools
-   [rCharts]() Interactive charts
-   [lattice]() Multivariate plots

-   [manipulate]() Interactive plots
-   [misc3d]() 3D plots
-   [mvtsplot]() Multivariate time series plot
-   [pca3d]() Three dimensional PCA plots
-   [longCatEDA]() Plot Categorical Longitudinal and Time-Series Data
-   [squash]() Color-based plots for multivariate visualization
-   [tourr]()
-   [latticist]() GUI for graphic exploratory data analysis using lattice
-   [tfplot]() Plot time series
-   [pathdiagram](http://cran.r-project.org/web/packages/pathdiagram/)
-   [RIGHT]() Interactive graphics via HTML
-   [clickme](https://github.com/nachocab/clickme) Interactive plots

-   [animation]() Create animations

-   [calibrate]() Axes calibration
-   [labeling]() Axis labeling
-   [directlabels]() Labels for multicolor plots. See also `Hmisc:::labcurve`.
-   [gplots]() Useful plotting tools
-   [magicaxis]() Pretty scientific plots (e.g. log scales)
-   [PlotRegionHighlighter]() Creates an envelope that surrounds a set of points
-   [compactr](http://cran.us.r-project.org/web/packages/compactr/index.html) Plots with compact axis notation
-   [plotflow](https://github.com/trinker/plotflow) Useful plotting functions
-   [sendplot]() Interactive plots with tool-tip content
-   [zoom]() Allow to zoom/navigate in any plot.
-   [Hmisc]()
-   [epade]() Easy plots
-   [prettyGraphs]() Publication-quality graphics.
-   [scagnostics]() Scatterplot diagnostics

### Colour

-   [colorRamps]() Colour palettes
-   [RColorBrewer]()
-   [colorspace]() HCL (perceptually-based) palettes
-   [colortools]()
-   [plotKML]() Colour palettes for mapping
-   [gplots]() Rich.colors palettes

-   [dataview]()
-   [tabplot]() Large datasets viz.
-   [tabplotd3]() Interactive
-   [R2SWF]() Convert R graphics to Flash
-   [sendplot]() Send interactive plots with tooltip content
-   [sjPlot]()
-   [squash]() Color-based plots for multivariate visualization

Reproducible Research - Report generation
-----------------------------------------

Check also CRAN Task View on [Reproducible Research](http://cran.r-project.org/web/views/ReproducibleResearch.html).

-   [brew]()
-   [knitr]() Dynamic report generation
-   [knitrBootstrap](https://github.com/jimhester/knitrBootstrap) create bootstrap styled HTML reports from Rmarkdown
-   [reports]() writing reports and presentations
-   [repmis]()
-   [rapport](http://cran.r-project.org/web/packages/rapport/index.html) report templating system
-   [pander]() Exploits pandoc to convert among multiple formats
-   [stargazer]() Easily create tables with regression outputs (directly from model objects)
-   [rtf](http://cran.r-project.org/web/packages/rtf/index.html) Produce Rich Text Format documents from R
-   [R2HTML](http://cran.r-project.org/web/packages/R2HTML/index.html) HTML reports from R
-   [xtable](http://cran.r-project.org/web/packages/xtable/index.html) Export R objects to HTML tables
-   [R2wd]() Write MS-Word documents from R
-   [rmarkdown]()
-   [sjPlot]()
-   [Rgitbook](https://github.com/jbryer/Rgitbook)
-   [table1xls](http://cran.r-project.org/web/packages/table1xls/index.html) Summary tables in Excel format.
-   [DescTools](http://cran.r-project.org/web/packages/DescTools/index.html)
-   [rctrack](http://www.stjuderesearch.org/site/depts/biostats/rctrack)

Parallelisation & Big Data
--------------------------

Check also [CRAN Task View on High Performance Computing](http://cran.r-project.org/web/views/HighPerformanceComputing.html).

-   [ff]() Big data management. See also package `ffbase`
-   [batch]() Batching routines in parallel
-   [bigdata]()

Niche & Species Distribution Modelling
--------------------------------------

-   [rgbif]() Access to GBIF data
-   [spocc]() Species occurrence data retrieval and mapping
-   [dismo]() Distribution modelling
-   [SDMTools](http://www.rdocumentation.org/packages/SDMTools)
-   [adehabitatHS](http://cran.r-project.org/web/packages/adehabitatHS/index.html)
-   [biomod2]() SDM ensembles
-   [hSDM]() Bayesian SDM
-   [maxlike]() SDM for presence-only data
-   [MigClim]() Implements dispersal in SDMs
-   [GRaF]() SDM using latent Gaussian random fields
-   [SightabilityModel]()
-   [EnvNicheR]()
-   [RinSp](http://cran.r-project.org/web/packages/RInSp/) Ecological niche metrics to measure individual specialization
-   [bdvis](https://github.com/vijaybarve/bdvis) Biodiversity data visualization
-   [phyloclim]() Includes functions for calculating niche overlap
-   [phylospacer]() Phyloclimates and phylomorphospaces
-   [PresenceAbsence]()
-   [mboost]() Decomposing environmental, spatial, and spatiotemporal components of species distributions (Hothorn et al. 2011 Ecol Monogr)
-   [usdm](http://www.rdocumentation.org/packages/usdm) Uncertainty analysis for species distribution models, esp. focused on positional uncertainty
-   [ModelMap](http://cran.r-project.org/web/packages/ModelMap/) Random Forest and Stochastic Gradient Boosting Models
-   [modEVA](http://modtools.wordpress.com)
-   [rangemapper]()
-   [comclim](http://cran.r-project.org/web/packages/comclim/index.html) Community climate statistics
-   [hypervolume](http://cran.r-project.org/web/packages/hypervolume/index.html) Modeling hypervolumes (species' niches)
-   [ppmlasso]() Point process models
-   [sdmvspecies]() Create virtual species for SDM
-   [virtualspecies](http://borisleroy.com/virtualspecies/) Create virtual species
-   [coenocline]() Simulate species presence and abundance along environmental gradients
-   [ecospat](http://cran.r-project.org/web/packages/ecospat/index.html) Many useful functions for niche & SDM (by A. Guisan's group)
-   [comclim]()
-   [Metadata]()
-   [ENiRG](http://cran.r-project.org/web/packages/ENiRG/index.html)

### Occupancy modelling

-   [detect]() Occupancy modelling
-   [marked]()
-   [unmarked]()
-   [stocc]() Occupancy modeling
-   [hierarchicalDS]() hierarchical analysis of distance sampling data

Climate
-------

-   [BerkeleyEarth](http://cran.r-project.org/web/packages/BerkeleyEarth/index.html) Climate data from Berkeley Earth database
-   [climates]() Tools for climate data (bioclim, downscaling, interpolation...)
-   [climatol]() Homogenisation of climate time series
-   [climstats]() Tools for climate data
-   [RMAWGEN]() Generate daily time series from monthly mean values
-   [climdex.pcic]() Computation of climate indices
-   [ClimClass](http://cran.r-project.org/web/packages/ClimClass/index.html) Climate Classification according to various indices
-   [RClimMAWGEN]() generate time series of climate indices
-   [chillR]() Climate and phenology analysis
-   [rWBclimate]() Lots of historical data and future projections
-   [rnoaa]()
-   [rghcnV3]() GHCN
-   [RNCEP]() weather data
-   [sirad]() Evapotranspiration etc
-   [SPEI]() Calculates Standardised Precipitation-Evaporation Index
-   [SPI]()
-   [weathermetrics]()
-   [raincpc]() Rainfall data
-   [Evapotranspiration](http://cran.r-project.org/web/packages/Evapotranspiration/index.html)
-   [GhcnDaily](http://cran.r-project.org/web/packages/GhcnDaily/index.html)
-   [crn](http://cran.r-project.org/web/packages/crn/index.html) Get data from Climate Reference Network
-   [iki.dataclim](http://cran.r-project.org/web/packages/iki.dataclim/index.html)
-   [FedData](http://cran.r-project.org/web/packages/FedData/index.html)
-   [stationaRy](https://github.com/rich-iannone/stationaRy) Hourly data worldwide
-   [tempcyclesdata](https://cran.r-project.org/web/packages/tempcyclesdata/index.html) Temperature data for ca. 8000 stations 1975-2013

GIS/spatial functionality
-------------------------

Check also CRAN Task View on [Spatial](http://cran.r-project.org/web/views/Spatial.html) and [Spatiotemporal](http://cran.r-project.org/web/views/SpatioTemporal.html) data. Check also CRAN Task View on [Web technologies](http://cran.r-project.org/web/views/WebTechnologies.html) for access to many GIS data.

-   [sp]()
-   [spacetime]()
-   [spatial]()
-   [rgeos]()
-   [rgdal]()
-   [raster]()
-   [rasterVis]()
-   [rts]() Raster time series
-   [spatial.tools]() Functions for raster processing, including parallel processing
-   [maptools]() Reading and manipulating geographic data
-   [PBSmapping]()
-   [taRifx.geo]() Useful spatial functions

-   [climstats]() Tools for climate and raster analyses

-   [countrycode]() Converts country codes
-   [cshapes]() Dataset of country boundaries
-   [gdistance]() Calculate distances and routes on grids
-   [geometry]() Mesh generation and surface tesselation
-   [geosphere]() Several GIS tools, esp. aimed at large scales (global)
-   [geonames]()
-   [mapproj]() Converts lat-lon data into projected coordinates
-   [geospacom]() Generates distance matrices from shape files and represents spatially weighted multilevel analysis results

-   [dismo]()
-   [fossil]()

-   [ggmap]() Plotting on Google Maps and OpenStreetMap
-   [GISTools]() Some further GIS tools (e.g. cloropleth maps)
-   [mapplots]() Data visualisation on maps
-   [maps]() Mapping tools
-   [OpenStreetMap]() Access to OpenStreetMap and Bing images
-   [osmar]()
-   [plotGoogleMaps]()
-   [RgoogleMaps]()
-   [R2G2]() Plot anything in Google Earth
-   [rCarto]() useful functions for mapping
-   [marmap]() working with bathymetric and topographic data
-   [plotKML]() Visualization of spatial and spatio-temporal objects in Google Earth
-   [MODISTools](http://cran.r-project.org/web/packages/MODISTools/)
-   [spGoogle](http://cran.r-project.org/web/packages/spGoogle) Interacting R with Google maps
-   [rMaps]() Interactive maps
-   [leafletR]()
-   [GEOmap]()
-   [rworldmap]()
-   [choroplethr]()
-   [mapmisc](http://cran.r-project.org/web/packages/mapmisc/index.html) Utilities for producing maps
-   [sinkr](https://github.com/menugget/sinkr) Some GIS functions and colour palettes

-   [rlandscape]() Generates spatial landscapes

-   [RSAGA]()
-   [spgrass6]()
-   [spsextante]()
-   [RPyGeo]()

Spatial Analysis
----------------

Check also CRAN Task View on [Spatial](http://cran.r-project.org/web/views/Spatial.html) and [Spatiotemporal](http://cran.r-project.org/web/views/SpatioTemporal.html) data.

-   [akima]() Interpolation
-   [automap]() Automatic interpolation (kriging)
-   [ecespa]() Point pattern analysis
-   [fields]() Splines, kriging, etc
-   [geoR]() Geostatistical analysis
-   [geoRglm]() Spatial GLMs
-   [geospt]() Geostatistics
-   [GeoXP]() Interactive exploratory analysis of spatial data
-   [gstat]() Geostatistics (variograms, kriging)
-   [intamap]() Automatic spatial interpolation
-   [geostatsp]() Geostatistics using SpatialPoints and rasters
-   [spatstat]()
-   [spBayes]()
-   [spdep]() Modeling spatial dependence
-   [stpp]() Point patterns
-   [rtop]() Spatial interpolation of data

Networks
--------

-   [bipartite]()
-   [igraph]() Network analysis and visualization
-   [enaR]()
-   [qgraph]()

Phylogenetics, phylogeography & comparative analysis
----------------------------------------------------

Check also [Phylogenetics CRAN task view](http://cran.r-project.org/web/views/Phylogenetics.html).

-   [ape]()
-   [adephylo](http://cran.r-project.org/web/packages/adephylo/index.html) Exploratory analysis
-   [Treesim]() Tree simulation. See also [treesimGM]()
-   [cladoRcpp]() Phylogenetic analysis of geographic ranges
-   [caper]() Comparative analysis
-   [PVR]() Phylogenetic Eigenvector Regression
-   [BEDASSLE](http://cran.r-project.org/web/packages/BEDASSLE/) Disentangles the effects of geographic and ecological isolation on genetic differentiation
-   [geiger]()
-   [geoscale]() Geological timescale
-   [phangorn]()
-   [phytools]()
-   [pegas]()
-   [picante]()
-   [stratigraph]()
-   [treebase]()
-   [MCMCglmm]()

Palaeoecology
-------------

-   [Bchron](https://cran.r-project.org/web/packages/Bchron/index.html) Bayesian chronologies
-   [Bclim](https://cran.r-project.org/web/packages/Bclim/index.html) Palaeoclimate reconstruction

Dendrochronology
----------------

-   [dplR](https://github.com/cran/dplR) Dendrochronology
-   [measuRing](https://github.com/cran/measuRing) Measuring ring width from scanned images.
-   [TRADER](https://github.com/pavel-fibich/TRADER) Tree Ring Analysis of Disturbance Events
-   [dendrobox](https://github.com/cszang/dendrobox) Interactive Tree-Ring Data Exploration Tool
-   [treeclim](https://github.com/cszang/treeclim) Modeling tree-climate relationships.
-   [bootRes](https://github.com/cszang/bootRes) Bootstrapped response functions.
-   [rwtocore](https://github.com/ltrr-arizona-edu/rwtocore) convert tree-ring measurements to drawings of cores (Ruby).
-   [triforce](https://github.com/ltrr-arizona-edu/triforce) read the TRiDaS dendrochronological standard and communicate with Tellervo servers.
-   [pointRes](https://cran.r-project.org/web/packages/pointRes/index.html) Calculate and plot pointer years
-   [BIOdry](https://cran.r-project.org/web/packages/BIOdry/index.html) Multilevel Modeling of Dendroclimatical Fluctuations
-   [sclero](https://cran.r-project.org/web/packages/sclero/) Measure Growth Patterns and Align Sampling Spots in Photographs
-   [dendrometeR](https://cran.r-project.org/web/packages/dendrometeR/) Analyzing Dendrometer Data.

Ecological analyses (miscellaneous)
-----------------------------------

Check also [Environmetrics CRAN Task View](http://cran.r-project.org/web/views/Environmetrics.html).

-   [betapart]() Beta diversity
-   [BiodiversityR]() GUI for biodiversity, suitability and community ecology analyses
-   [cheddar]() Analysis and visualisation of community data
-   [coexist]() Modelling species coexistence
-   [DSpat]() Distance Sampling
-   [fossil]() Species richness, species-area curves, beta-diversity
-   [indicspecies]() Assess associations between different species and sites (e.g. indicator species)
-   [IPMpack]() Integral Projection Models
-   [marked]() Mark-recapture analysis
-   [MBI]() Multiple Beta Diversity indices
-   [rmeta]() Meta-analysis
-   [neighlikeli]() Neighborhood models
-   [earlywarnings]()
-   [MAR1]() Multivariate Autoregressive Modeling for Analysis of Community Time-Series Data
-   [ncbit]() retrieve NBCI taxonomic data
-   [spacodiR]() Spatial and Phylogenetic Analysis of Community Diversity
-   [Reol](http://cran.r-project.org/web/packages/Reol/) Interface to Encyclopedia of Life
-   [SPECIES]() Species richness and diversity analysis
-   [BayesComm]() Bayesian multivariate binary regression models for analysis of ecological communities
-   [pom]() Patch occupancy models
-   [capwire]() Population size estimation from non-invasive sampling
-   [vegclust]() Fuzzy clustering of vegetation data
-   [popbio]()
-   [demoniche](http://demoniche.r-forge.r-project.org/) Spatially-explicit demographic modelling
-   [BEDASSLE]() Disentangles the effects of geographic and ecological isolation on genetic differentiation.
-   [rphylopic](https://github.com/sckott/rphylopic) Organisms silhouettes
-   [GSIF]() Global Soil Information
-   [Rarity]() Rarity indices
-   [taxize]() Taxonomy etc.
-   [Taxonstand]()
-   [sExtinct]() Extinction analysis based on sightings
-   [BEDASSLE](http://cran.r-project.org/web/packages/BEDASSLE/) Disentangles the effects of geographic and ecological isolation on genetic differentiation
-   [ecomodtools]() Simulation models (inc. dispersal)
-   [comclim](http://cran.r-project.org/web/packages/comclim/index.html) Community climate statistics
-   [kernelPop](http://cran.r-project.org/web/packages/kernelPop/index.html) Spatially explicit population genetic simulations
-   [siplab](http://cran.r-project.org/web/packages/siplab/index.html) spatially explicit individual-based vegetation models.

Miscellaneous
-------------

-   [digitize]() Extract data from plots
-   [downloader]() Download files from internet
-   [hwriter]() Outputs R objects in HTML format
-   [installr]() Automatically update R
-   [mail]() Send email notifications from R
-   [sendmailR]()
-   [mosaic]()
-   [audiolyzR]() Creates audio representations of common plots
-   [WebDevelopR]() Web development
-   [ProjectTemplate]()
-   [gpk]() Datasets for educational uses
-   [NCmisc](http://cran.r-project.org/web/packages/NCmisc/)
-   [PubMedWordcloud]() Create wordcloud with abstracts.
-   [alm]() Altmetrics
-   [rAltmetric]() Retrieves data from altmetrics.com
-   [RMendeley]()
-   [rcrossref]()
-   [RefManageR]() Bibliography Management.
-   [refnet]() Reads and works with data from ISI Web of Knowledge
-   [scholar]()
-   [slidify]()
-   [sos]() Search R packages
-   [source.gist]()
-   [twitteR]()
-   [DescTools](http://cran.r-project.org/web/packages/DescTools/index.html)

R programming
-------------

-   [codetools]() Code analysis tools
-   [compiler]() Compile R code (e.g. a function) to speed it up
-   [devtools]()
-   [formatR]() Format and tidy R code
-   [gtools]() Useful functions (e.g. mixedsort)
-   [inline]() Define functions in C code within R
-   [iterators]() Tools for iteration on many different types of objects (see also `itertools`)
-   [Kmisc]() Data reshaping, table and plot generation from RMarkdown
-   [operators]() Additional operators
-   [simFrame]() A framework to work with simulations
-   [simSummary]()
-   [staticdocs]() Create webpage with package documentation
-   [testit]() Testing R packages
-   [testthat]() Testing R code
-   [tester]() Test characteristics of R objects
-   [assertthat]().
-   [httr]() Working with URLs and HTTP
-   [RCurl]() HTTP interface
-   [stringr]()
-   [Rd2Roxygen]()
-   [Rxoygen2]()
-   [scrapeR]()
-   [aprof](http://cran.r-project.org/web/packages/aprof/index.html) Profiling code and visualizing results to evaluate where is best to focus code optimization.
-   [regex](https://github.com/richierocks/regex) Regular expressions made easier.