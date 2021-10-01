CrappyFish data analysis<br/>

List of scripts and output files:

1_Download_data.sh\ 
  - Download the data from Novogene's servers

2.0 FastQC AdapterRemoval and Bismark\
  - QC with FastQC, trimming with adapterremoval, alingnment/methylation calling with bismark\
Output files: "SAMPLE.pair1.truncated.gz_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz"

3.0 methylKit.R\
  - Importing the methylation call files into MethylKit\
  - Filtering for 10X and 5X coverage in all samples\
  - Finding significant deferentially methylated cytosines\
  - Finding significant deferentially methylated regions\
Output files: MethylKit_19samples_10X_5X_1X_210901.RData and MethylKit_19samples_10X_5X_210901.RData

3.1 PCA.R \
  - Make PCA of 5X and 10X coverage\
  - Plot using ggplot\  

3.2 Manhattan.R \
  - Make Sunset Manhattan plots and QQ-plots of 5X and 10X coverage\

3.3 Counting DMCs and DMRs.R\
  - Count the number of DMCs and DMRs\

4.0 Make BED for all covered CpGs, DMCs, and DMRs.R\
  - Convert the methylkit files to bedfiles for BEDtools\

4.1 Divide genome into features\
  - Base unix commands for downloading the reference genome and dividing into features\ 
  - Count the overlaps between the feature types and the CpGs, DMCs and DMRs\ 

4.2 Feature figure and G-test.r\
  - Make barchart with information from 4.1\
  - Use G-test for finding significant enriched or depleated feature types\

5.1 GO.R\
  - Download GO and KEGG terms for the genes found using 5.0\
  - Perform ORA\

6.0 Download nanopore data > call methylation\ 
  - Download the nanopore data from ERDA\
  - Index and map with minimap2\
  - Call methylation\
  - fix annoing bug with positions\

6.1 Count Nanopolish CpGs\
  - Count the number of CpG sites methylation called at different thresholds\

6.2 Exstract ROIs from methylKit\
  - Exstract methylation information specificly for ROIs\

6.3 Heatmap comparing ONT and WGBS\
  - Make heatmap for comparing methylation calls\ 

7.0 ORA overlap.R\
  - Find overlap of significant GO terms in different studies\

8.0 R Computerome\
How to launch r scripts at computerome2\

Im a molecular biologist and geneticist by training, and bioinformatiction by need.\
Good coding practice is not always applied, feel free to add suggestions and comments 
