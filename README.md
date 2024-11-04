# Hormonal contraceptive use and anaemia: A nation-wide pharmacoepidemiological study from Northern Europe 
This repository contains R code to reproduce the analyses for the manuscript Hormonal contraceptive use and anaemia: A nation-wide pharmacoepidemiological study from Northern Europe (DOI) by Ekroos S., Toffol E., Heikinheimo O., Haukka J. and Arvas M. 

The code is included in five separate files. The first part, "hormonal_contraception_DAG.Rmd", creates the Directed Acyclic Graph used for model selection. The second part, "mkdataAnemiaV20240515_MA_V2.R", preprocessed the cohort data and sets up the Nested Case Control study. The third part, "data_exploration_V2.Rmd", allows the user to describe the cohort and original study population. The fourth part, "mkAnemiaModels_V2.R", runs the conditional logistic regression models and sensitivity analysis, which can then be visualized in the fifth part, "anemia_and_hc_V2.Rmd".

Information on packages needed to run the code is included in the .Rmd:s and .R:s. 

For questions regarding the code or manuscript, please contact S Ekroos (sofie (dot) ekroos (at) helsinki (dot) fi).
