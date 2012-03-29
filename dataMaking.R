#Load in data
load("/Users/will/Documents/Stamatakis/data/butterflies/communities/yearlyUKBMS/workspace.RData")
pseudo.data$siteNo <- factor(pseudo.data$siteNo)

#Subset the data into phylogenetic partitions
pseudo.data$sci_name[pseudo.data$sci_name %in% tree$tip.label[1:11]] <- "A"
pseudo.data$sci_name[pseudo.data$sci_name %in% tree$tip.label[12:25]] <- "B"
pseudo.data$sci_name[pseudo.data$sci_name %in% tree$tip.label[26:32]] <- "C"
pseudo.data$sci_name[pseudo.data$sci_name %in% tree$tip.label[33:46]] <- "D"
pseudo.data$sci_name[pseudo.data$sci_name %in% tree$tip.label[46:55]] <- "E"
merged.data <- with(pseudo.data, data.frame(count = tapply(sindex, paste(sci_name, siteNo, year, sep="-"), sum)))
merged.data$sci_name <- rep("NA", nrow(merged.data))
merged.data$siteNo <- rep("NA", nrow(merged.data))
merged.data$year <- rep("NA", nrow(merged.data))
for(i in seq(nrow(merged.data))){
	if((i %% 94)==0) if((i %% 940)==0) cat("|") else cat(".")
	merged.data$sci_name[i] <- strsplit(rownames(merged.data)[i], "-", fixed=TRUE)[[1]][1]
	merged.data$siteNo[i] <- strsplit(rownames(merged.data)[i], "-", fixed=TRUE)[[1]][2]
	merged.data$year[i] <- strsplit(rownames(merged.data)[i], "-", fixed=TRUE)[[1]][3]
}
merged.data$sci_name <- factor(merged.data$sci_name)
merged.data$siteNo <- factor(merged.data$siteNo)
merged.data$year <- factor(merged.data$year)


#Function to make a usable dataset from a sample
make.output <- function(data, community=-1){
	
	#Replace all the names with something reasonable
	tree$tip.label[1:11]
	tree$tip.label[12:25]
	tree$tip.label[26:32]
	tree$tip.label[33:46]
	tree$tip.label[46:55]
	
	if(community > 0){
		community <- levels(data$siteNo)[community]
		data <- data[data$siteNo==community,]
	} else community <- "Everything"
	
	blanks <- names(which(with(data, tapply(count, sci_name, sum))==0)) #DODGY!!!
	data <- data[!data$sci_name %in% blanks,]
	output <- data.frame(abundance = data$count, species=data$sci_name, comm_name = data$siteNo, year = data$year)
	write.table(output, paste(community, ".csv", sep=""), sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#Write out the output
setwd("/Users/will/Documents/code/lottery/lottery/data")
for(comm in seq(levels(merged.data$siteNo))){
	make.output(merged.data, comm)
}
make.output(merged.data)
