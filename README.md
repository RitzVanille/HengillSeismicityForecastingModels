# Hengill Seismicity Forecasting Models
This project is in the frame of the COSEISMIQ project (http://www.coseismiq.ethz.ch/en/home/) and comprises the models and data prepared for deliverables 6 and 7.

## Folder structure

![alt text](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/FolderStructure.png "Folder structure")

## Run down of streams
The figure below shows a simplified visualisation of the interactions of the different folders and files of the project
![alt text](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/StructureHengillModels.png "representation of the interactions between folders and files")
with model blocks in yellow, static data in green, analysis scripts in red and model outputs in purple.

## Data
### Seismic data
The seismic catalogue can be acquired openly through FDSN servers: http://www.coseismiq.ethz.ch/en/dissemination/catalogs/.
We use the medium quality catalogue (SED_auto_MQ, event score >-5, within the Hengill regional limits).
A statistical analysis of the whole dataset shows that the overall magnitude of completeness (Mc) is 0.3, with a b-value of 0.93. For all models, we assume the Mc=0.3 and filter the catalogue above this magnitude. 
![alt text](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/Data/Catalogue_Statistics.png "Statistical analysis of the SED_auto_MQ catalogue")
Supporting information for the catalogues, seismic data acquisition and processing can be found in Grigoli et al., 2021 (submitted).

### Hydraulic data
Injection and production rates obtained from Reykjavik Energy are distributed in a symmetrical bivariate Gaussian distribution around the well-head to account for fluid migration in the subsurface with daily granularity.
The rates, volumes and cumulative volumes (and their associated times) are provided as is. In the future a automated script to create the distribution will be made available as well as the raw hydraulic data once they are public.

## Models

Two families of model are implemented: a Seismogenic Index type and a class of ETAS (Epidemic Type Aftershock Sequence) models. 
Details can be found for each model in their respective folder's README

### Seismogenic Index Model
Based on Mignan et al., 2016; Shapiro, 2018; Broccardo et al., 2020.
This model is implemented without a shut-in phase, with cumulative volumes spatially distributed and productivity `a_{fbi}` and the Gutenberg-Richter b-value estimated by minimizing the residuals with a maximum likelihood optimisation function (modified from Broccardo et al., 2017).

### ETAS type models
Three version of ETAS are implemented:
- the null hypothesis `ETAS_{0}` 
- an Advanced ETAS model with spatially varying background rate
- a prototype of ETAS with hydraulic coupling based on injection and production rate

## Learning and forecasting periods
We take as a base Learning Period (LP) all the data recorded from Dec. 1st 2018 up to Feb. 1st, 2020. Using this data, we calibrate all the competing models. Then we perform a daily forecast and update the LP with the additional information. We do this process of training and forecasting repeatedly until the end of the dataset on Feb. 1st 2021.

## Model comparison and evaluation scheme
### Visualisation
From each model's output, maps of the forecasted seismicity rates and log likelihood, temporal evolution of the forecasted number of events and temporal evolution of the log likelihood of the forecast can be computed with the VisualisationModels.m script.
### Models comparison and evaluation
For each model, the N-test (as defined by the CSEP: https://cseptesting.org/) for each forecast in computed.
Using ETAS0 as the null hypothesis, ComparisonModels.m computes the cumulative information gain of the models relative to this null hypothesis.

## References

- Broccardo, M., Mignan, A., Wiemer, S., Stojadinovic, B., & Giardini, D. (2017). Hierarchical Bayesian modeling of fluid‐induced seismicity. Geophysical Research Letters, 44(22), 11-357.
- Broccardo, M., Mignan, A., Grigoli, F., Karvounis, D., Rinaldi, A. P., Danciu, L., ... & Wiemer, S. (2020). Induced seismicity risk analysis of the hydraulic stimulation of a geothermal well on Geldinganes, Iceland. distributed in a symmetrical bivariate Gaussian distribution around the well-head to account for fluid migration in the subsurface.
- Grigoli, F., Clinton, J., Diehl, T., Kaestli, P., Scarabello, L., ..., Rinaldi, A., Ritz, V., Sanchez-Pastor, P., & Wiemer, S. (submitted). Monitoring microseismicity in the Hengill Geothermal Field, Iceland. Scientific Data.
- Mignan, A., Broccardo, M., Wiemer, S., & Giardini, D. (2017). Induced seismicity closed-form traffic light system for actuarial decision-making during deep fluid injections. Scientific reports, 7(1), 1-10.
- Nandan, S. (2017). Towards a physics based epidemic type aftershock sequence model (Doctoral dissertation, ETH Zurich).
- Shapiro, S. A. (2018). Seismogenic index of underground fluid injections and productions. Journal of Geophysical Research: Solid Earth, 123(9), 7983-7997.

## License:
This project is free software: You can redistribute it and/or modify it under the terms of the European Union Public Licence (EUPL v.1.2) or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence"); You may not use this work except in compliance with the Licence.

The project is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the Licence for the specific language governing permissions and limitations under the Licence.
