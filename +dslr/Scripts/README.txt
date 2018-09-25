Check Documentation.mlx for detailed explanations and examples

This folder contains the following scripts

	CameraController.m
	DarkAnalysisPlots.m
	DarkEvolution.m
	DarkEvolutionPlot.m
	DCNew.m
	designMatrixEstimation.m
	FlatCapsBlocks.m
	FlatRuns.m
	FlatsPlot.m
	meanCheck.m
	temperatureNoise.m

###############################################################

Some scripts are obsolete due to later iterations. Below are some descriptions for the scripts in this folder.


CameraController.m - downloaded from MathWorks. This script interfaces with the digiCam controller webserver and offers much controlover camera functionality.

DarkAnalysisPlots.m - accepts outputs of DCNew.m to plot median value and variance vs. exposure time, channel-by-channel. It also plots linear fits for both.

DarkEvolution.m - captures images at fixed exposure time for extended periods of time to analyze evolution of dark current as sensor heats up.

DarkEvolutionPlot.m - plots outputs of DarkEvolution.m, median exp. val and std vs. exposure time.

DCNew.m - the newer alternative to DarkGlobalCalc.m. Measures variance and mean pixel-by-pixel through stacks of dark images. 

designMatrixEstimation.m - sets up model for median and variance and estimates gain, read noise, and dark current using a design matrix.

FlatCapsBlocks.m - captures flats, takes median and variance through constant exp time stack, then divides those images into squares for easier handling.

FlatRuns.m - combines FlatCapsBlocks.m and designMatrixEstimation.m. Outputs gain, DC, RN for each iteration.

FlatsPlot.m - scatter plots output of FlatCapsBlocks. More useful as a diagnostic tool to check whether any exposure runs are abnormal in behavior.

meanCheck.m - checks mean for each channel. diagnostic to determine channel ordering

temperatureNoise.m - choose four images to compare low and high ISO frames at 2 different temperatures for visual inspection


############################################################