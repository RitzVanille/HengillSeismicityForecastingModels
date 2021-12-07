## Seismogenic Index model
This model is based on a purely empirical relationship between the injected volume and the 
seismicity (e.g. Shapiro, 2018).  Other forms accounts, for example for the injection rate rather than the injected volume, as well as for an exponential decay of the seismicity after shut-in (Mignan et al., 2017; Broccardo et al., 2017).
The proposed model provides a deterministic forecast with no flow model, and the expected number of events follow the relationship:

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/SeismogenicIndex/Nexpected_SeismogenicIndex.png)

The minimization function is based on the Maximum Likelihood Estimate of Broccardo et al., 2017, modified for continuous operation (no shut-in): 

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/SeismogenicIndex/MLE_SeismogenicIndex.png)
