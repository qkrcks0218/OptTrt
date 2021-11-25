# Replication Files for Park et al. (2021)

This Github repository contains replication files for [Park et al. (2021)](https://arxiv.org/abs/2111.09932 "OptTrt").

## Data 

Data folder contains replication files for the data analysis in Section 5 of the paper.

### Download the DHS data
* Access the data from https://dhsprogram.com/data/available-datasets.cfm.
* Download SNKR70FL.DTA , SNKR7HFL.DTA , SNKR7IFL.DTA , SNKR7ZFL.DTA , SNKR81FL.DTA and place the data sets in the following paths.
	* Data/DHS/2014/SNKR70FL.DTA
	* Data/DHS/2015/SNKR7HFL.DTA
	* Data/DHS/2016/SNKR7IFL.DTA
	* Data/DHS/2017/SNKR7ZFL.DTA
	* Data/DHS/2018/SNKR81FL.DTA
* Download "SNGE81FL.ZIP" from https://dhsprogram.com/data/dataset/Senegal_Continuous-DHS_2018.cfm?flag=1. 
* Unzip the zip file and place the 10  files including "SNGE81FL.shp" in "Data/GPS_2018" folder.

### Code 

* 1.Cleaning.R
	* The DHS data sets are cleaned. 
	* The training and test data sets are saved as "Sen_Training_1417.csv" and "Sen_Test_18.csv", respectively.
	* The test set including the identifier is saved as "TestSetID.csv".
* 2.Direct_NF.R
	* The outcome regressions and propensity scores are estimated. 
	* The results are saved as "NF\_Train\_Q\_B####.csv" files in "NF" folder.
	* The description is written in Section 3.3 of the main paper.
	* The code will take a lot of time, so it is recommended to use parallel computing.
	* To implement parallel computing, split BATCH=1,...,100 into multiple jobs.
* 3.Direct\_InitialPT.R.R
	* The initial values for the DC algorithm are evaluated, and are saved as "NF\_Train_Adj\_B####.csv" files in "NF" folder.
	* The description is written in Section A.3 of the supplementary material.
* 4.Direct\_CV.R
	* The empirical risks under each parameter are evaluated.
	* The results are saved as "Summary\_CV\_S###.csv" files in "CV" folder.
	* The description is written in Section A.4 of the supplementary material.
	* The code will take a lot of time, so it is recommended to use parallel computing.
	* To implement parallel computing, split ss.iter=1,...,10 and BATCH=1,...,3300 into multiple jobs.
* 5.Direct\_RULE.R
	* The direct OMARs for the training sets are obtained where parameters are chosen from CV.
	* The results are saved as "RULE\_S###\_Test.csv" files in "RULE" folder.
	* The code will take a lot of time, so it is recommended to use parallel computing.
	* To implement parallel computing, split ss.iter=1,...,10, nf.iter=1,...,100 into multiple jobs.
* 6.Direct\_Testset.R
	* The direct OMARs for the test sets are obtained where parameters are chosen from CV.
	* The results are saved as "RULE\_S###\_Test.csv" files in "RULE" folder.
* 7.Indirect.R
	* The indirect OMARs for the test sets using the Lasso and random forest are obtained.
	* The results are saved as "RULE\_LR\_Test.csv" files in "NF\_LR" folder.
* 8.GeoRegion.R
	* The simulation results are merged with GPS file "Data/GPS\_2018/SNGE81FL.shp"
	* The summary is saved as "RULE\_Mergged.csv".
* 9.Summary.R
	* The simulation results are summarized and the plots are drawn.
* 10.Check\_Assumptions.R
	* The assumptions are assessed as in Section 5.2.

### Folder

* CV: CV folder contains cross-validation results.
* DHS: DHS folder contains 2014-2018 Senegal DHS data sets.
* GPS\_2018: GPS_2018 folder contains the GPS information of 2018 Senegal DHS data set.
* NF: NF folder contains the estimates of the nuisance functions.
* NF\_LF: NF\_LR folder contains the indirect OMAR estimates for the test sets using the Lasso and random forest. 
* RULE: RULE folder contains the direct OMAR estimation results and the direct OMAR estiamtes for the test sets.
* Plot: Plot folder contains the graphical summaries of the data analysis and the assessment of the assumptions.

## Simulation

Simulation folder contains replication files for the simulation in Section 4 of the paper.

### Code

* 1.DataGenerate.R: The data sets used for the simulation are generated. 
* 2.Direct_NF.R: The outcome regressions and propensity scores are estimated. 
* 3.Direct_InitialPT.R: The initial values for the DC algorithm are evaluated.
* 4.Direct_CV.R: The empirical risks under each parameter are evaluated.
* 5.Direct_RULE.R: The direct OMAR estimates for the training sets are obtained where parameters are chosen from CV.
* 6.Direct_Testset.R: The direct OMAR estimates for the test sets are obtained where parameters are chosen from CV.
* 7.Indirect.R: The indirect OMAR estimates for the test sets using the Lasso and random forest are obtained.
* 8.Summary.R: The simulation results are summarized.

### Folder
	
* CV: CV folder contains cross-validation results.
* Data: Data folder contains generated data under the simulation models.
* NF: NF folder contains the estimates of the nuisance functions.
* NF_LR: NF_LR folder contains the indirect OMAR estimates for the test sets using the Lasso and random forest. 
* RULE: RULE folder contains the direct OMAR estimation reesults and the directc OMAR estiamtes for the test sets.
* Plot: Plot folder contains the graphical summaries of the simulation.

## Additional Files

* MySL.R contains functions used for implementing superlearner algorithm and estimating the nuisance functions.
* OptTrt_Sim_Function.R contains functions used for the OMAR using the direct and indirect methods.

## References
Chan Park, Guanhua Chen, Menggang Yu, and Hyunseung Kang (2021) **Optimal Allocation of Water and Sanitation Facilities To Prevent Communicable Diarrheal Diseases in Senegal Under Partial Interference**, _arXiv:2111.09932_ [[link](https://arxiv.org/abs/2111.09932 "OptTrt")]