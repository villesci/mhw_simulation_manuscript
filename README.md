# Marine Heatwave Simulation Analysis

Code and data to replicate results in Villeneuve and White 2024.

There are two major R markdown scripts that replicate analyses and figures from our paper:

* mhw_sim_analysis.Rmd - Bulk of analysis. Can knit entire document, or run select chunks. Some chunks are very time and computationally intensive, so we provide .RDs data files that can be read in at checkpoints to speed rendering of the document. These data files can be replicated by running chunks that are not evaluated in the markdown. A rendered .html file of the markdown is also available. 

* conceptual_figure.Rmd - Script to produce conceptual figure panels in the manuscript. A rendered .html file of the markdown is also available. 

Additonally, there are two R scripts containing functions from other works used in our analysis. 
* hour_fxns_heatwaveR - modified functions from Schlegel et al. 2018, the heatwaveR package. Functions altered to allow for inputs of data in hourly format. More information within mhw_sim_analysis.Rmd and manuscript text.

* Thermal_landscape_functions.R - modified functions from Rezende et al. 2020 to create thermal tolerance landscapes and produce estimates of mortality using dynamic survival model. More information within mhw_sim_analysis.Rmd and manuscript text.

Finally, there are two folders.

* data - .RDs files that are saved during computationally intensive sections of the analysis, so that they can be read in independent of running their origin code. These files are not necessary to replicate our results.

* result_figs - folder for storing figures from analysis. 
