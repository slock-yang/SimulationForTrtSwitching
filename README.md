# SimulationForTrtSwitching

Both settings' codes are written in the "RankPreserve_fun.R", where "sim_data" and "sim2_data" functions are introduced separately for scenario 1 and 2. Function "sim" is also introduced for estimating the treatment effects from the generated data applying recensoring, deleting, ITT analysis, aalen model and our proposed doubly robust estimator.

To conduct simulations, firstly, use "sim_data" and "sim2_data" to generate the data and store it in the assigned directory which is through an input parameter, "newdir". Then, use "sim" to extract the data (with the input parameter "newdir") from the directory and estimate.

 A sample is displayed in the "Simulation" section of file "RankPreserve.R".
 
 ## Algorithm
 
 The used computation codes are all stored in the directory "Function", where codes using "Rcpp" are stored in directory "integral". Though several parallel solutions are listed there, we typically use the "integral.R" which is a wrap of Rcpp codes in directory "integral".
