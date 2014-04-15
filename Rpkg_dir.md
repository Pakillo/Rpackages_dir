-   Data sources
-   Data management
-   General statistics
-   Regression tools
-   Ordination & Multivariate Analysis
-   Bayesian/MCMC
-   Plotting & Visualisation
-   Reproducible Research - Report generation
-   Parallelisation & Big Data
-   Niche & Species Distribution Modelling
-   Climate
-   GIS/spatial functionality
-   Spatial Analysis
-   Networks
-   Phylogenetics, phylogeography & comparative analysis
-   Palaeoecology
-   Ecological analyses (miscellaneous)
-   Miscellaneous
-   R programming

Updated 2014-04-15.

These are packages that I often use or, alternatively, I need only rarely but don't want to forget about. Of course, there are many other useful packages out there (e.g. at [CRAN](http://cran.r-project.org/web/packages/available_packages_by_name.html) or [GitHub](https://github.com/search?q=package&type=Repositories&ref=advsearch&l=R). Check also [CRAN task views](http://cran.r-project.org/web/views/).

[F. Rodriguez-Sanchez](http://tinyurl.com/frod-san)

Data sources
------------

Check CRAN Task View on [Web technologies](http://cran.r-project.org/web/views/WebTechnologies.html).

Data management
---------------

-   [data.table]() Tools for managing data frames
-   [datamerge]() Cleaning factors
-   [foreign]() Read data in many different formats
-   [gdata]() Data manipulation
-   [Hmisc]() Frank Harrell's miscellaneous tools
-   [Kmisc]() Data reshaping, table and plot generation from RMarkdown
-   [lubridate]() Dates and times
-   [multitable]() Manipulate multiple arrays
-   [datacheck]() Tools for checking data consistency
-   [dplyr]()
-   [reshape2]()
-   [tabplot]() Large datasets viz.
-   [tabplotd3]() Interactive

General statistics
------------------

Regression tools
----------------

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
-   [hier.part]() Variance partitioning
-   [relaimpo]() Relative importance of predictors
-   [nlme]() Mixed models
-   [lme4]() Mixed models
-   [likelihood]() Maximum Likelihood Estimation
-   [lmtest]() Diagnostic checks for linear regression
-   [MARSS]() Multivariate Autoregressive State-Space Modeling
-   [MCMCglmm]() Mixed models fitted by MCMC
-   [mgcv]() GAM fitting
-   [miscF]() Multivariate Normal regression, spatial Bayesian mixed models, piecewise regression...
-   [mlogit]() Multinomial regression
-   [mvinfluence]() Influence Measures and Diagnostic Plots for Multivariate Linear Models
-   [msm]() Multi-state Markov and hidden Markov models in continuous time
-   [odprism]() Power analysis for mixed models
-   [PCovR]() Principal Covariates Regression (Reducing collinear predictor variables to a few components and regressing on them)
-   [HLMdiag]() Diagnostic tools for hierarchical models (fitted with lme4)
-   [stremo]() Learning Structural Equation Models
-   [rockchalk]() Regression Estimation and Presentation
-   [R2STATS]() GUI for fitting and comparing GLM and GLMM in R
-   [LMERConvenienceFunctions]()
-   [R2admb]()
-   [rms]() Regression modeling strategies
-   [Zelig]()

Ordination & Multivariate Analysis
----------------------------------

Check also [Environmetrics CRAN Task View](http://cran.r-project.org/web/views/Environmetrics.html) and [Multivariate](http://cran.r-project.org/web/views/Multivariate.html)

-   [ade4](http://cran.r-project.org/web/packages/ade4/index.html) Multivariate data analysis and display
-   [dave]() Data Analysis in Vegetation Ecology
-   [ecodist]() Dissimilarity-based analysis (ordination, Mantel tests...)
-   [labdsv]() Includes plotting functions
-   [mvabund]() Analysing multivariate data (upscaling from individual species models)
-   [eigenprcomp](http://cran.r-project.org/web/packages/eigenprcomp/) Computes confidence intervals for the proportion explained by the first 1,2,k principal components
-   [vegan]()

Bayesian/MCMC
-------------

Check also [Bayesian CRAN Task View](http://cran.r-project.org/web/views/Bayesian.html).

-   [bayess]() Funcions used in the book 'Bayesian essentials with R'
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

-   [coda]() MCMC output analysis
-   [boa]() Analyses of MCMC output (like `coda`)
-   [ggmcmc]() Graphic analysis of MCMC output
-   [mcmc]()
-   [mcmcmplots]() Plot MCMC output
-   [bmk]() MCMC diagnostics
-   [scapeMCMC]() MCMC diagnostic plots
-   [superdiag]() testing MCMC noncovergence

-   [adaptMCMC](http://cran.r-project.org/web/packages/adaptMCMC/index.html) Generic MCMC sampler
-   [BRugs]() Interface to OpenBUGS
-   [BayesX]() Structured Additive Regression
-   [R2BayesX]()
-   [rjags]()
-   [runjags]()
-   [R2jags]()
-   [R2OpenBUGS]()
-   [R2WinBUGS]()
-   [rstan]()
-   [filzbach]()
-   [iBUGS]() Interface to BUGS/JAGS
-   [rube]()
-   [INLA]() Integrated Nested Laplace Approximation
-   [LaplacesDemon]()
-   [LearnBayes]()
-   [MCMCpack]() MCMC samplers
-   [MHAdative]() MCMC sampler
-   [glmmBUGS]()

-   [CARBayes]() Spatial models with CAR
-   [spBayes]() spatial Bayes

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

-   [ggplot2]()
-   [GGally]() extension to ggplot
-   [ggsubplot]() Embedding subplots within plots
-   [ggthemes]() Extra themes, scales and geoms for ggplot

-   [plotly]()
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

-   [animation]() Create animations

-   [calibrate]() Axes calibration
-   [labeling]() Axis labeling
-   [directlabels]() Labels for multicolor plots
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

-   [colorRamps]() Colour palettes
-   [RColorBrewer]()
-   [colorspace]()
-   [colortools]()

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
-   [detect]() Occupancy modelling
-   [hSDM]() Bayesian SDM
-   [maxlike]() SDM for presence-only data
-   [MigClim]() Implements dispersal in SDMs
-   [GRaF]() SDM using latent Gaussian random fields
-   [SightabilityModel]()
-   [marked]()
-   [unmarked]()
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
-   [stocc]() Occupancy modeling

Climate
-------

-   [BerkeleyEarth]() Climate data from Berkeley Earth database
-   [climates]() Tools for climate data (bioclim, downscaling, interpolation...)
-   [climatol]() Homogenisation of climate time series
-   [climstats]() Tools for climate data ()
-   [RMAWGEN]() Generate daily time series from monthly mean values
-   [climdex.pcic]() Computation of climate indices
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

-   [Bchron]() Bayesian chronologies
-   [Bclim]() Palaeoclimate reconstruction

Ecological analyses (miscellaneous)
-----------------------------------

Check also [Environmetrics CRAN Task View](http://cran.r-project.org/web/views/Environmetrics.html).

-   [betapart]() Beta diversity
-   [BiodiversityR]() GUI for biodiversity, suitability and community ecology analyses
-   [cheddar]() Analysis and visualisation of community data
-   [coexist]() Modelling species coexistence
-   [dplR]() Dendrochronology
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
-   [fylopic]() Organisms silhouettes
-   [GSIF]() Global Soil Information
-   [Rarity]() Rarity indices
-   [taxize]() Taxonomy etc.
-   [Taxonstand]()
-   [sExtinct]() Extinction analysis based on sightings
-   [BEDASSLE](http://cran.r-project.org/web/packages/BEDASSLE/) Disentangles the effects of geographic and ecological isolation on genetic differentiation
-   [ecomodtools]() Simulation models (inc. dispersal)

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
