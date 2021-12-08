## ETAS0 the null hypothesis

ETAS0 is the simplest version of the Epidemic Type Aftershock Sequence (ETAS) model. This model views the seismicity as composed of background earthquakes and triggered earthquakes. The background earthquakes are driven by the forces of plate tectonics and anthropogenic factors such fluid injection and extraction. These earthquakes trigger a cohort of aftershocks which then trigger their own aftershocks and so on. In its simplest from the ETAS model describes the conditional seismicity rate of magnitude m events, λ(t,x,y,m | H<sub>t</sub>) at any location (x,y) and time t as:

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/ETAS_models/seismicityrate_ETAS0.png)

where µ is the background intensity function, assumed to be independent of time, and H<sub>t</sub> stands for the history of the process up to time t. 
The triggering kernel g quantifies the temporal and spatial influence of pas events onto future events:

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/ETAS_models/kernel_ETAS0.png)
