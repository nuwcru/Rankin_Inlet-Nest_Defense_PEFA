

# ![Picture2 (2)](https://user-images.githubusercontent.com/56907107/168872856-4bbd7c7d-eb2d-42dc-ab8c-79e7b97adba5.jpg)

# NuWCRU-- Does fluctuating selection maintain variation in nest defense behavior in Arctic peregrine falcons (_Falco peregrinus tundrius_)?

## Project Description
Within populations, individuals often exhibit repeatable among-individual differences in behavior, and in some cases, these are linked to different fitness outcomes. Several mechanisms have been proposed to contribute to the maintenance of repeatable amongindividual variation in behavior. Here we study nest defense behavior in Arctic peregrine falcons (_Falco peregrinus tundrius_) over two successive breeding seasons (2018 and 2019) to evaluate the importance of three potential mechanisms that could underlie the maintenance of among-individual variation in this trait; state-dependent behavior, assortative mating, and fluctuating selection. Nest defense was measured as the response made by peregrines towards human observers during standard nest visits; high nest defense was characterized by close approaches to the observer, and low nest defense was characterized by maintaining greater distance from the observer. Nest defense scores ranged from 0m (i.e., contact with observer) to 600m.

This work was conducted near the community of Rankin Inlet, Nunavut, Canada (62°49′N, 92°05′W), which is located on the western coast of Hudson Bay. The rugged coastline provides numerous potential nest sites for cliff nesting species (e.g., Peregrine falcons, Rough-legged hawks (_Buteo lagopus_) , Canada goose (_Branta canadensis_), Common eider (_Somateria mollissima_), and Common ravens (_Corvus corax_). The inland portion of the study area is dominated by rolling upland hills and eskers that are also suitable for nesting. 


## General Project File Structure

```
├── README.md                                                      <- The top-level README including general project descriptions
|
├── data
│   ├── data_2018.csv                                              <- Data ready for modeling
│   ├── data_2019.csv                                              <- Data ready for modeling
│   ├── nest_defense_test_final.csv                                <- Data ready for modeling
│   └── nest_defense_test_final_20220318.csv                       <- Data ready for modeling
│
└── scripts
    ├── State_dependence_repeatability_NestDefense_Rmd_files       <- folder that holds figures displayed in markdown file
    ├── State_dependence_repeatability_NestDefense_PEFA.R          <- all models combined used in main text and ESM
    └── State_dependence_repeatability_NestDefense_PEFA.md         <- markdown file with model output and figures

```

## Data 
  * **data_2018**- data used for assortative mating and selection gradients for 2018
  * **data_2019**- data used for assortative mating and selection gradients for 2018 
  * **nest_defense_test_final**- data used for main effects (Raw data that is used to produce data_2018 and data_2019 data frames)
  * **nest_defense_test_final_20220318**- data used for correlations across breeding contexts

## Markdown file
* Click the link below to access markdown file with full workflow and output from each model. 
* https://github.com/nuwcru/Rankin_Inlet-Nest_Defense_PEFA/blob/main/scripts/State_dependence_repeatability_NestDefense_Rmd.md

## Authors
* Nick Gulotta, University of Georgia (Nickolas.Gulotta@uga.edu)
* Kimberley Mathot, University of Alberta (mathot@ualberta.ca)


