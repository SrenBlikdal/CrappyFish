CrappyFish data analysis<br/>

List of scripts:\
1_Download_data.sh\ 
  - Download the data from Novogene's servers\

2_FastQC_Trimming_Alignment_MethylationCalling_pbs_batch_job_thinnode.sh\
  -QC with FastQC, trimming with adapterremoval, alingnment/methylation calling with bismark\

3.0_methylKit.R\
  -Importing the methylation call files into MethylKit\
  -Filtering for 10X and 5X coverage in all samples\
  -Finding significant deferentially methylated cytosines\
  -Finding significant deferentially methylated regions\
  
Output files generated are MethylKit_19samples_10X_5X_1X_210901.RData and MethylKit_19samples_10X_5X_210901.RData

3.1_PCA.R\
  -Make PCA of 5X and 10X coverage\
  -Plot using ggplot 

3.2_Manhattan.R
  -Make Sunset Manhattan plots and QQ-plots of 5X and 10X coverage\

