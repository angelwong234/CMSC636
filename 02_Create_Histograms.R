#=====================================================================================
#                                                                [CMSC 636] Drake Wong
#  Histogram Script
#
#=====================================================================================
library(tidyverse)

# Read data
type="BMDC"
file_location = "20_04_18 10 categories large/"
loc = paste0(file_location, type, ".csv")
data <-read_csv(loc, col_names=TRUE, col_types = cols(.default = "c"))

# Convert all columns to double, except for rownames
data <- data %>%
  select(-'X1') %>%
  mutate_if(is.character, as.double) %>%
  na_if(0)

# Extract gene expression values into data.frame
gene_expression<- as.data.frame(data) %>% 
  select(2:ncol(data))

# Set directory to store histograms
dir <- paste0(file_location, type)
new_path<-dir.create(file.path(dir)) 

# Loop through each gene. Create histogram of all samples for a particular gene.
for (i in 1:ncol(gene_expression)) {
  # Set filename in appropriate directory
  mypath <- file.path(dir,paste("myplot_",colnames(gene_expression)[i], ".jpg", sep = ""))
  jpeg(file=mypath, width=1000, height=1000)
  print(i) 
  # Creation of the plot without axis labels, tick marks, or axis numbers.
  trial<-ggplot(gene_expression, aes(x=gene_expression[,i])) +
    geom_histogram(binwidth=1)+
    theme_bw()+
    theme(axis.title = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  plot(trial)
  dev.off()
}

