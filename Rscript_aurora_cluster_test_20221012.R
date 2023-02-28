# Get the working directory to double check you are where you are supposed to be!
getwd()


# Install renv 
install.packages("renv")
# Initiate renv to create a project specific library of your packages. 
renv::init()
renv::snapshot()


# Install devtools and BiocManager and premessa
install.packages("devtools")
install.packages("BiocManager")

BiocManager::install("flowCore")


library(devtools)
install_github("ParkerICI/premessa")

# Check your antigen parameters and remove uneccessary channels here. Data will be exported in "renamed". Move to working directory and delete original fcs files.
# Copy "parameter" and "most common" from shinyapp for panel file "fsccolname" and "antigen"
library(premessa)
paneleditor_GUI()







# Load packages  - or check box in the packages list opposite.
library(readxl) 
library(CATALYST)
library(cowplot)
library(flowCore)
library(scater)
library(SingleCellExperiment)
library(openxlsx)
library(uwot)
library(rmarkdown)
library(knitr)

# From here I would switch to the markdown file. I keep these script notes for my reference and for annotations of the script.

# SCRIPTS FOR MARKDOWN FILE

# To load the metadata file from your working directory. If you ever need to re-load then unload first using remove(md) 
md <- "xxx.xlsx" 

#read excel reads the metadata file
md <- read_excel(md)                                           

# shows the dataframe of the metadata file
head(data.frame(md)) 


# this defines the flowset which contains the Cytof fcs files and also performs transformation on the data needed for analysis
fs <- read.flowSet(md$file_name, transformation = F, truncate_max_range = F)

# load the panel xls file from the working directory
panel <- "xxxx.xlsx" 
panel <- read_excel(panel)                                         
head(data.frame(panel))  

# extra checks to make sure column names match and panel names match 

all(panel$fcs_colname %in% colnames(fs))
# True

setdiff(fs@frames$filename.fcs@parameters@data$desc,panel$antigen)
# Character (0)

#define factors and levels in the metadata so here just "control" for condition as only looking at control samples
md$condition <- factor(md$condition, levels = c("transplant"))      
md$date <- factor(md$timepoint, levels = c("none"))
md$patient_id <- factor(md$patient_id, levels = c("ILC002", "ILC003", "ILC004", "ILC005", "ILC006", "ILC007",
                                                  "ILC008", "ILC009", "ILC011"))

md$sample_id <- factor(md$sample_id, levels = md$sample_id[order(md$condition)])   

# prepdata include factors of interest - need to specify if not defaults. sce is single cell experiment and defines what is used in later analysis ie the flowset, the metadata and the panel 
# ?prepdata - to get info on usage
# for Aurora work need to change the cofactor to 150
# for Aurora flow fcs files need to add FACS = TRUE
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 150, FACS = TRUE, features = panel$fcs_colname, md_cols = list(file = "file_name", id = "sample_id", 
                                                                            factors = c("condition", "timepoint", "batch", "patient_id")))

# Data checks prior to cluster analysis
# to count the numbers of cells in the sce

n_cells(sce)

# to plot the counts in the sce

plotCounts(sce, color_by = "patient_id")

# non-redundancy score (see paper for explanation)

plotNRS(sce, features = type_markers(sce), color_by = "condition")


sce <- cluster(sce, features = "type", 
               xdim = 10, ydim = 10, maxK = 20, 
               verbose = FALSE, seed = 1) 


# plots the MDS (see paper for explanation)

CATALYST::pbMDS(sce, color_by = "condition", size_by = TRUE) 

# plots expression heatmap  - Heatmap of aggregated marker expressions. Try changing scale to look at
# effect of scaling on expression values
# - When `scale = "first"`, the specified assay data will be scaled between 0 and 1 using lower (`q`) and upper (`1-q`) quantiles as boundaries, where `q = 0.01` by default. This way, while losing information on absolut expression values, marker expressions will stay comparable, and visualization is improved in cases where the expression range varies greatly between markers. 

# - When `scale = "last"`, assay data will be aggregated first and scaled subsequently. Thus, each marker's value range will be [0,1]. While all comparability between markers is lost, such scaling will improve seeing differences across, e.g., samples or clusters.

# - When `scale = "never"`, no scaling (and quantile trimming) is applied. The resulting heatmap will thus display *raw* pseudobulk data (e.g., median expression values).


plotExprHeatmap(sce, bin_anno = TRUE, row_anno = TRUE, scale="never")

# Clustering
# set seed for clustering - random seed selection

set.seed(1234)

# run clustering on the sce - type_marker means it will include any markers defined as "type" - all in this case (see panel) - maxK is the maximum number of clusters you want the algorithm to create from the data

sce <- cluster(sce, features = type_markers(sce),    
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

# plot of the cluster heat map to show expression level of all the markers in each cluster. 
# You can change k=metat20 to any number below 20 to show smaller clustering.
# Dont worry about the error that comes up here "deprecated" means that this command is going to
# be phased out soon; I'll need to re-write the code based on the new plot types....but not tonight :)

plotClusterHeatmap(sce, hm2 = NULL, k = "meta20", m = NULL,               
                   cluster_anno = TRUE, draw_freqs = TRUE)  

# plot of the cluster heat map to show expression of a specific marker across all samples and clusters. 
# You can change "CD294" to any of the antigens listed in the panel. 
# You will need to use this to identify the clusters later. 
# Have a go changing the antigen to see if you can identify the T-cells (CD3), the B-cells (CD19 and CD20), the Monocytes (CD14), gamma-delta T-cells (CD3 + TCRgd), Regulatory T-cells (CD3 + Foxp3 + CD4 + CD25). 
# Also you will need ot click zoom to open the plot - this is where jupyter lab was better :(

plotClusterHeatmap(sce, hm2 = "CD294", k = "meta20", draw_freqs = TRUE)

# TSNE and UMAP  are two different dimension reduction analysese so you can visualise the clustering. 
# Just to explain, the clustering is an algorithm that labels each cell with a cluster number.
# Whereas the dimension reduction is an algorithm to visualise that cluster definition. 
# Sometimes people get confused and think that a TSNE plot is a type of clustering but it isnt.

# I couldnt run the UMAP  as had package "uwot"  missing  - not sure how this happened. 

# use this to install if needed.
# install.packages("uwot")


set.seed(1234)                                                                      
sce <- runDR(sce, dr = "TSNE", cells = 500, features = "type") 
sce <- runDR(sce, dr = "UMAP", cells = 1e3, features = "type")


# packages required for further analysis and plots

library("ggplot2")
library("scatterplot3d")


# plots of the dimension reduction - you can change "TSNE" here for "UMAP" etc. 
# Also you can change what you color by and facet_wrap by eg. different antigen to color the plot by.
# Ive put some examples below.

plotDR(sce, "TSNE", color_by = "meta20") + facet_wrap("condition") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))     

plotDR(sce, "UMAP", color_by = "meta20") + facet_wrap("patient_id") +          
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))     

plotDR(sce, "TSNE", color_by = "CD127") + facet_wrap("condition") +          
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))  


plotDR(sce, "TSNE", color_by = "GranB")  + facet_wrap("condition") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))    

# This is a way to see the dimension reduction plots colored by expression of particular antigens rather than colored by cluster number (meta20). 
plotDR(sce, "UMAP", color_by = "GranB")  
+ guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))    


# This code allows you to see the TSNE and UMAP on the same graph to compare visualisation
p1 <- plotDR(sce, "TSNE", color_by = "meta20") +                              
  theme(legend.position = "none")                                           
p2 <- plotDR(sce, "UMAP", color_by = "meta20")                                
lgd <- get_legend(p2 +                                                        
                    guides(color =  guide_legend(ncol = 2, override.aes = list(size = 3))))   
p2 <- p2 + theme(legend.position = "none")                                    
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))   


# everything from here on is about labeling the clusters as cell types so dont worry about this for now.
# SKIP code from line 174 to 194.
merging_table1 <- "merging2.xlsx"   
merging_table1 <- read_excel(merging_table1)          
head(data.frame(merging_table1)) 

merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("NK_Cytotoxic", "Contaminant_cells", "NK_immature", "ILC", "ILC2_precursor")) 

sce <- mergeClusters(sce, k = "meta10",                                  
                     table = merging_table1, id = "merging2")     

plotDR(sce, "UMAP", color_by = "merging2")

plotClusterHeatmap(sce, k = "meta10", m = "merging2")

plotClusterHeatmap(sce, k = "merging5")

plotDR(sce, "TSNE", color_by = "merging5") + facet_wrap("condition") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))  

# Some code to show cluster population abundances
FDR_cutoff <- 0.05
plotAbundances(sce, k = "meta20", by = "cluster_id")
plotAbundances(sce, k = "meta20", by = "sample_id")
plotAbundances(sce, k = "meta20", by = "cluster_id", shape = "timepoint")

# Here is code to export cluster abundance 
# R-script to get cluster frequencies in dataframe; x = single cell experiment, 
 # k = cluster code (meta20, or if relabeled then merging id.)

ns <- table(
  cluster_id = cluster_ids(sce, k= "meta20"), 
  sample_id = sample_ids(sce))
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)


dff <- as.data.frame(cluster_ids(sce))



# export dataframe to csv or xlsx - make sure to define the directory you want it saved in.
# and then create a file name - Ive called it mydata.xlsx here - see below. Use getwd() to help.
getwd()

# write.csv (df, “path/mydata.csv”, row.names = TRUE
        
write.xlsx (df, "/Users/rebeccapayne/OneDrive - Newcastle University/Students/2021_Mres/Cytof_tutorial_Rstudio/ForRuthandMorven/mydata.xlsx", row.names = TRUE)

# check the xls file is exported and you are done!!!! Whoop!!           


#define new sce with subset of data
newSCEday30 <- filterSCE(sce, date == "day30")

sub <- filterSCE(sce, patient_id == "ILC009_BMT019")


## Create new sce (sub1) with monocyte clustes removed

sub1 <- filterSCE(sce, cluster_id %in% c(2, 1, 4, 5, 6, 7, 8, 9, 12, 17, 20, 11, 13, 14, 19, 18, 16), k = "meta20")


# FlowSOM cluster of sub1

sub1 <- cluster(sub1, features = "type", 
                xdim = 10, ydim = 10, maxK = 20, 
                verbose = FALSE, seed = 1) 

# TSNE Dimension reduction of sub1
set.seed(1234)                                                                      
sub1 <- runDR(sub1, dr = "TSNE", cells = 500, features = "type")


plotDR(sub1, "TSNE", color_by = "meta20") + facet_wrap("sample_id") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))   




# To get counts of cells in each cluster and for each sample_id in dataframe
n_cells <- table(sample = sce$sample_id, cluster = cluster_ids(sce, "merging3"))

# To export to directory (check directory first just in case)
getwd()
write.xlsx(n_cells, "n_cells.xlsx")




# plot counts
pC <- plotCounts(sce, group_by = "patient_id")
pC$data


# To extract clusters from a sce (called sub here) and export as fcs to directory. If want to extract one patient or timepoint then us filterSCE first

sub <- filterSCE(sce, patient_id == "ILC009_BMT019")

extractClusters(sub, k = "meta20" , clusters = NULL, as = c("fcs"), out_dir = "/Users/rebeccapayne/OneDrive - Newcastle University/ILC_study/Cytof_FCS_files/WBC/Matched_FCSfiles/Matched_fcs_LymphocyteGate/MATCHED_FCS",
                verbose = TRUE)




extractClusters(sub1, k = "meta20" , clusters = NULL, as = c("fcs"), out_dir = "/Users/rebeccapayne/OneDrive - Newcastle University/ILC_study/Cytof_FCS_files/WBC/Matched_FCSfiles/Matched_fcs_LymphocyteGate/MATCHED_FCS",
                verbose = TRUE)




# to get list of all packages used in session
sessionInfo()




