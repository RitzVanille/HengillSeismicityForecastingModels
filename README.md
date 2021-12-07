# Hengill Seismicity Forecasting Models

## Folder structure

![alt text](https://raw.githubusercontent.com/RitzVanille/HengillSeismicityForecastingModels/main/FolderStructure.png "folder structure od project")

## Run down of streams

![alt text](https://raw.githubusercontent.com/RitzVanille/HengillSeismicityForecastingModels/main/StructureHengillModels.png "folders interactions")
with model blocks in yellow, static data in green, analysis scripts in red and model outputs in purple.

## Data
### Seismic data

### Hydraulic data

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

## Model comparison and evaluation scheme


## References

- Broccardo, M., Mignan, A., Wiemer, S., Stojadinovic, B., & Giardini, D. (2017). Hierarchical Bayesian modeling of fluid‐induced seismicity. Geophysical Research Letters, 44(22), 11-357.
- Broccardo, M., Mignan, A., Grigoli, F., Karvounis, D., Rinaldi, A. P., Danciu, L., ... & Wiemer, S. (2019). Induced seismicity risk analysis of the hydraulic stimulation of a geothermal well on Geldinganes, Iceland.
- Mignan, A., Broccardo, M., Wiemer, S., & Giardini, D. (2017). Induced seismicity closed-form traffic light system for actuarial decision-making during deep fluid injections. Scientific reports, 7(1), 1-10.
- Nandan, S. (2017). Towards a physics based epidemic type aftershock sequence model (Doctoral dissertation, ETH Zurich).
- Shapiro, S. A. (2018). Seismogenic index of underground fluid injections and productions. Journal of Geophysical Research: Solid Earth, 123(9), 7983-7997.

## License:
This project is free software: You can redistribute it and/or modify it under the terms of the European Union Public Licence (EUPL v.1.2) or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence"); You may not use this work except in compliance with the Licence.

The project is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the Licence for the specific language governing permissions and limitations under the Licence.
