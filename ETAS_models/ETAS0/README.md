## ETAS0 the null hypothesis

ETAS0 is the simplest version of the Epidemic Type Aftershock Sequence (ETAS) model. This model views the seismicity as composed of background earthquakes and triggered earthquakes. The background earthquakes are driven by the forces of plate tectonics and anthropogenic factors such fluid injection and extraction. These earthquakes trigger a cohort of aftershocks which then trigger their own aftershocks and so on. In its simplest from the ETAS model describes the conditional seismicity rate of magnitude m events, $\lambda$(t,x,y,m | Ht) at any location (x,y) and time t as:
![alt text](https://raw.githubusercontent.com/RitzVanille/HengillSeismicityForecastingModels/main/ETAS_models/seismicityrate_ETAS0.png "Seismicity Rate computation by ETAS_0")
where $\mu$ is the background intensity function, assumed to be independent of time, and Ht stands for the history of the process up to time t. 
The triggering kernel g quantifies the temporal and spatial influence of pas events onto future events:
![alt text](https://raw.githubusercontent.com/RitzVanille/HengillSeismicityForecastingModels/main/ETAS_models/kernel_ETAS0.png "ETAS_0 triggering kernel")
