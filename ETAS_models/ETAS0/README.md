## ETAS0 the null hypothesis

ETAS0 is the simplest version of the Epidemic Type Aftershock Sequence (ETAS) model. This model views the seismicity as composed of background earthquakes and triggered earthquakes. The background earthquakes are driven by the forces of plate tectonics and anthropogenic factors such fluid injection and extraction. These earthquakes trigger a cohort of aftershocks which then trigger their own aftershocks and so on. In its simplest from the ETAS model describes the conditional seismicity rate of magnitude m events, λ(t,x,y,m | H<sub>t</sub>) at any location (x,y) and time t as:

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/ETAS_models/seismicityrate_ETAS0.png)

where µ is the background intensity function, assumed to be independent of time, and H<sub>t</sub> stands for the history of the process up to time t. 
The triggering kernel g quantifies the temporal and spatial influence of pas events onto future events:

![equation](https://github.com/RitzVanille/HengillSeismicityForecastingModels/raw/main/ETAS_models/kernel_ETAS0.png)

### Requirements

This model involves inversions of parameters. The current version of the codes provided runs these inversions in batch on the ETH cluster EULER (https://scicomp.ethz.ch/wiki/Euler).

### Computational steps

The ETAS model is decomposed in several steps governed by the ```Wrapper_Experiment.m``` script which invokes the relevant wrappers and scripts for each step. The steps are as follows:
1. Data loading
2. Inversion of parameters with ```Wrapper_ETAS_Invert_verCoSeismiq.m``` (on the cluster, can be done in parallel)
3. Simulation of the sythetic catalogues with ```Wrapper_Sim_ETAS_verCoSeismiq.m``` (multiple realisations, on the cluster)
4. Model evaluations to perform listing with ```Wrapper_FailedEval_ver3.m```
5. Evaluation of the models with ```Wrapper_Eval_ETAS_verCoSeismiq.m```
After these steps are done, the results can be analysed with the visualisation and model comparison scripts.


### Contact
In case of questions on how to run the ETAS models, please contact [Shyam Nandan](mailto:snandan@ethz.ch)
