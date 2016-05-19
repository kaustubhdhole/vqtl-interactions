# vqtl-interactions
  QTL mapper || vQTL mapper || QTL-QTL interactions || vQTL-vQTL interactions
				
				Mapping vQTL(& QTL)
	To map vQTL(& QTL) you need to have the following files in the folder “Map_QTL” (or "vqtl-interactions"):
		a.	getphenogeno.R
		b.	mapping.R
		c.	getpeaks.R
		d.	BF_test.R
	System Requirements: 
		a.	For best results “R” version 3.1 or above
		b.	The following R packages must be installed	
			a)	rqtl
			b)	parallel
			c)	doparallel
			d)	foreach
			You may need to run the command 
			install.packages(package_name) 
			in R to install each one of them.

	Once the above requirements are met, follow the below steps to perform mapping vQTL(&QTL):
		Step 1: Place your desired phenogeno input file in the folder “Map_QTL” and rename the input file to “phenogeno.csv” (A sample file is supplemented for illustration purposes.)
		Step 2: Open “mapping. R” in a text editor and set the 2 variables no_of_env and nperm according to requirements.
		Step 3: Open “get_peaks.R” in a text editor and set the desired LODBF cutoff defined by the variable and enter   the chromosome starting and ending locations. 
		Step 4: Open R and change the working directory of R to Map_QTL.
		Step 5: Execute the following two commands in R:
			source(“mapping.R”)
			source(“get_peaks.R”)

			Check your output in the file “DataFile1.csv”.

	Note: After executing the command source(“mapping.R”), the folder Map_QTL/vQTL will contain a csv file named according to the environments. The BFLOD values for all the markers are present in them. The BFLOD plots will be generated in the Map_QTL/vQTL/Plots folder.
	The above codes also map QTL simultaneously. For QTL output you may want to check the folder “Map_QTL/QTL”. 
					***



				Deming Deviations
	To calculate the Deming deviations, you need to have the following files in the folder Map_QTL/bins:
		a.	DemingDeviations.R
		b.	Bins.csv
		c.	First_markers_of_bins.csv
		d.	QTL_and_vQTL.csv
	The R package “mcr” must be installed. 
	Once the requirements are met, perform the following steps:
		Step 1: Change the working directory of R to “Map_QTL\Bins”
		Step 2: Run the following script:
			source(“DemingDeviations.R”)
	Your output will be available in the file “DataFile2.csv” in the Maps_QTL\Bins folder.
					***


				Interactions
	To get the 2QTL and 2vQTL interactions, you need to have the following files in the folder Map_QTL/interactions:
		a.	interaction.R
		b.	vQTL_vQTL.py 
	Once the requirements are met, perform the following steps:
		Step 1: After performing peak extraction from getpeaks.R, place the 2 folders “QTL” and “vQTL” in the folder Map_QTL/interactions
		Step 2: Change the working directory of R to “Map_QTL\Bins”
		Step 3: Execute the Python script vQTL_vQTL.py
			vQTL_vQTL interactions will be generated in the folder 2vQTL.
			(Two folders 2QTL and 2vQTL are given for illustrative purposes.)
		Step 4: To collate the data, set the environments in interaction.R and run the R script
			source(“interaction.R”)
			Your output will be available in the file “DataFile3.csv” in the Maps_QTL\interactions folder.
					***

