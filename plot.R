
library(dplyr)
library(ggplot2)
library(readxl)
rm(list=ls())


## set the name of the input EXCEL file here
srcF = "../input/Input filr.xlsx"

## set this variable to be SVG or PDF depending upon what type of output you want
output_type = "PDF";



## read in the data from 'srcF'
rawd <- suppressMessages(read_excel(srcF, sheet="Table"))

## The first column is constant. It tells you the name of the Networks to display along the y-axis
## The data is in chunks of 4 columns: pvalue, gene ratio, module number, number of genes

## create a new data.frame that has the information cleaned up
## Each row of 'idx' represents one module as it is laid out in the 'rawd' data.frame
start_column <- grep("P-value", colnames(rawd))
end_column <- start_column + 3
idx <- data.frame(start=start_column, end=end_column)

d <- data.frame()
for(i in 1:nrow(idx)) {
	tmp <- rawd[, c(1,idx$start[i]:idx$end[i])]
	colnames(tmp) <- c("network", "pval", "geneRatio", "moduleNum", "numGenesInMod")
	tmp <- dplyr::select(tmp, network, moduleNum, numGenesInMod, geneRatio, pval)
	d <- rbind(d, tmp)
}
d <- d[complete.cases(d), ]
## We need the p-value to be represented as -log10(p-value)
d$negLog10P <- -log10(d$pval)


## Arrange the data by module number.
d <- arrange(d, moduleNum, network)

## make sure the networks are plotted in the same order as they appear in the input file
d$network <- factor(d$network, levels=rev(rawd$Network))


## Create a new text label for x-axis
d$xTxt <- paste0(d$moduleNum,"\n(",d$numGenesInMod,")")
tmp <- select(d, moduleNum, xTxt) %>% distinct() %>% arrange(moduleNum)
d$xTxt <- factor(d$xTxt, levels=tmp$xTxt)
rm(tmp)

## Convert module numbers to factors
d$moduleNum <- as.factor(d$moduleNum)


## Define how to color the dots
dot_color_scale <- scale_color_gradient(low="#4575b4", high="#d73027")


## Define overall theme here
overall_theme <- theme(
	axis.title=element_blank(),
	axis.text=element_text(color="black"),
	plot.margin = unit(c(0.5, 0.5, 2, 0.5), "lines")
)


## Define an annotation object to be the legend for the x-axis
## The 'x' and 'y' variables here define where the label is placed on the plot.
## The first "data point" of the plot starts at coordinates (1,1)
x_axis_label <- annotate(
	geom="text",
	x=0.2,y=-0.3,label="Module\n(No. of genes)",
	size=3,
	hjust=0
)

## You will need these for placing the x_axis_label object in the final plot
totNumModules <- length(unique(d$moduleNum))
totNumNetworks <- length(unique(d$network))



## Make the final plot
p1 <- ggplot(d) + 
	geom_point(aes(x=xTxt, y=network, color=negLog10P, size=geneRatio)) + 
	dot_color_scale + 
	x_axis_label + ## this is the x-axis label
	coord_cartesian(xlim=c(1,totNumModules), ylim=c(1,totNumNetworks), clip="off") + ## this is required for placing the x-axis table
	theme_bw() + 
	overall_theme



if(output_type == "SVG") {
	#-------------------------------------------------------------------------------
	## This is for SVG output
	svg(filename="../output/dfermin_plot.svg", width=11, height=8)
	plot(p1)
	dev.off()
	#-------------------------------------------------------------------------------
}


if(output_type == "PDF") {
	#-------------------------------------------------------------------------------
	## This is for PDF output
	pdf(file="../output/dfermin_plot.pdf", width=11, height=8)
	plot(p1)
	dev.off()
	#-------------------------------------------------------------------------------
}

