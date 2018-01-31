

library(ggplot2)
library(phyloseq)
# source('http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')
library(ape)
# library(ICC)
library(reshape2)
library(matrixStats)
library(plyr)
library(gdata)
library(gtools)

# ******************************************************* #
#load microbiome data ####
# ******************************************************* #

if (dir.exists("/users/tg/jwillis/SLL")) { home_dir <- "/users/tg/jwillis/SLL"
} else if (dir.exists(sprintf("%s/SLL", getwd()))) { home_dir <- sprintf("%s/SLL", getwd())
} else if (dir.exists("~/Downloads/SLL")) { home_dir <- "~/Downloads/SLL"
}

# level <- 3
level <- 6

taxlevel <- sprintf('tax%s',level)

if (level==3) {
  
  taxfile <- sprintf("%s/Part_1/mothur/all_tax%d.csv", home_dir, level)
  data.tax <- read.delim(taxfile, row.names = 1, header=T)
  tax.table <- tax_table(do.call(rbind, lapply( rownames(data.tax), function(x) strsplit(as.character(x), ";")[[1]] )))
  phyla   <- as.character(tax.table[,'ta2'])
  classes <- as.character(tax.table[,'ta3'])
  for (i in 1:length(classes)) {
    if (classes[i] == 'unclassified') {
      rownames(data.tax)[i] <- paste('unclassified',substring(phyla[i],1,3),sep='.') #e.g. unclassified.Bac
    } else {
      rownames(data.tax)[i] <- classes[i]
    }
  }
  
  
} else if (level==6) {
  
  # taxfile <- sprintf("/Users/owner1/SLL/all_taxids%d.csv", level)
  taxfile <- sprintf("%s/Part_1/mothur/all_taxids%d.csv", home_dir, level)
  data.tax <- read.delim(taxfile, row.names = 1, header=T)
  tax.table <- tax_table(do.call(rbind, lapply( data.tax[,1], function(x) strsplit(as.character(x), ";")[[1]] )))
  data.tax <- data.tax[,2:length(data.tax)]
  
  genuses  <- as.character(tax.table[,'ta6'])
  for (i in 1:length(genuses)) {
    if (genuses[i] == 'unclassified') {
      rownames(data.tax)[i] <- paste('unclassified',substring(rownames(data.tax)[i],4),sep='') #e.g. unclassified.61432
    } else {
      rownames(data.tax)[i] <- genuses[i]
    }
  }

}

taxa_names(tax.table) <- rownames(data.tax)
data.tax[is.na(data.tax)] <- 0
# Filter out empty species
# data.tax <- data.tax[rowSums(data.tax) != 0,]

#ignore control samples
data.tax <- data.tax[,grepl("^Sdata\\d.*_SLL", colnames(data.tax))]
#change identifier to match those from the survey
colnames(data.tax) <- gsub("Sdata\\d.*_","",colnames(data.tax))


# ******************************************************* #
# Questionaire ####
# ******************************************************* #
# data.que <- read.delim("/Users/owner1/SLL/SLL_survey_data_repaired.csv", header=T)
data.que <- read.delim(sprintf("%s/Part_1/Surveys/SLL_survey_data_repaired.csv", home_dir), header=T)

# ignore samples of duplicated identifiers bc cannot determine which is correct or what the other should be changed to
dubs <- as.character(data.que[duplicated(data.que[,"Qsamples"]),"Qsamples"])
x <- 1
both_dubs <- c()
for (i in 1:length(rownames(data.que))) {
  if (data.que[i,"Qsamples"] %in% dubs) {
    both_dubs[x] <- as.numeric(rownames(data.que)[i])
    x <- x+1
  }
}
data.que <- data.que[-both_dubs,]
rownames(data.que) <- gsub('ll','LL', gsub("-","",as.vector(data.que[,"Qsamples"])))

#collect list of questions and remove question from top of data.que
full_questions <- data.que[1,]
data.que <- data.que[-1,]
data.que <- apply(data.que, 2, factor)
data.que <- as.data.frame(data.que) #because apply converts it to a matrix
colnames(data.que)[102] <- "Socioeconomic" #shortening super long question name

#need a better way to do this...make values null or something:
# data.que <- data.que[data.que[,'Q2'] != 'No Sabe/No Contesta',]
# data.que <- data.que[data.que[,'Q8'] != 'No Sabe/No Contesta',]
# data.que <- data.que[data.que[,'Q15'] != 'No Sabe/No Contesta',]
# data.que[data.que[,'Q2'] == 'No Sabe/No Contesta', 'Q2'] <- NA       #gender
# data.que[data.que[,'Q8'] == 'No Sabe/No Contesta', 'Q8'] <- NA       #ethnicity
# data.que[data.que[,'Q15'] == 'No Sabe/No Contesta', 'Q15'] <- NA     #zona municipal
# 
# data.que[data.que[,'Q29'] == 'No Sabe/No Contesta', 'Q29'] <- NA     #caries
# data.que[data.que[,'Q29.1'] == 'No Sabe/No Contesta', 'Q29.1'] <- NA #extr nervio
# data.que[data.que[,'Q29.2'] == 'No Sabe/No Contesta', 'Q29.2'] <- NA #dientes perdidos
# data.que[data.que[,'Q29.3'] == 'No Sabe/No Contesta', 'Q29.3'] <- NA #dientes reconstruidos
# 
# data.que[data.que[,'Q29']>10 ,'Q29'] <- 'No Sabe/No Contesta'
# data.que[data.que[,'Q29.1']>10 ,'Q29.1'] <- 'No Sabe/No Contesta'
# data.que[data.que[,'Q29.2']>10 ,'Q29.2'] <- 'No Sabe/No Contesta'
# data.que[data.que[,'Q29.3']>10 ,'Q29.3'] <- 'No Sabe/No Contesta'



# ******************************************************* #
# Construct phyloseq object ####
# ******************************************************* #
otu.table <- otu_table(data.tax, taxa_are_rows=T)
sample.data <- sample_data(data.que)

sll_tree <- read.tree(file = "/users/tg/jwillis/SLL/Part_1/SLL1.tree")
# st <- taxa_names(SLL)[ ! startsWith(taxa_names(SLL), "unclassified")]
# st[ ! st %in% sll_tree$tip.label]
sll_tree$tip.label[sll_tree$tip.label=="Rhodoferax"] <- "Albidiferax" # NCBI automatically changed this one


SLL <- phyloseq(otu.table, sample.data, tax.table)
SLL_with_tree <- phyloseq(otu.table, sample.data, tax.table, sll_tree)


# one sample had a pH of 4, which must be a typo since the measurable range was 5-10. So must remove it
sample_data(SLL)[ sample_data(SLL)[,"Q00.PH"]==4 ,"Q00.PH"] <- "No Sabe/No Contesta"

colnames(tax_table(SLL)) <- c("Kingdom","Phylum","Class","Order","Family","Genus")


# *********** #
# Plot boxes of diversity values in each school to show batch effect in first 5 schools ####
SLL_test <- SLL
rich_all <- estimate_richness(prune_taxa(taxa_sums(SLL_test)>0, SLL_test))

measure <- "Shannon"
# measure <- "Simpson"

div.by.school <- as.data.frame( cbind(rich_all[,measure], sample_data(SLL_test)[,"Qsamples.1"]) )
colnames(div.by.school) <- c("value","School_ID")

counts <- function(x) {
  ifelse(measure=="Shannon", y_add <<- 0.55, y_add <<- 0.12)
  return( data.frame(y=median(x) + y_add, label=length(x)) )
}

ggplot(div.by.school, aes(x=reorder(School_ID,-value,FUN=median), y=value, 
                          fill=reorder(School_ID,-value,FUN=median))) +
  geom_boxplot() + theme_minimal() +
  # theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold")) +
  # ggtitle(sprintf('%s per sample for %s', measure, "School_ID")) +
  xlab("School_ID") + ylab(measure) + guides(fill=FALSE) +#scale_fill_hue(name='School_ID') +#
  stat_summary(fun.data = counts, geom='text', size=6)
# *********** #



# remove samples from schools 1-5 because of batch effect during sequencing
for (i in 1:5) {
  good_school <- rownames(data.frame(sample_data(SLL)[ sample_data(SLL)[,'Qsamples.1'] != i]))
  SLL <- prune_samples(good_school, SLL)
  SLL_with_tree <- prune_samples(good_school, SLL_with_tree)
}

# remove taxa with 0 counts
SLL <- prune_taxa( names(rowSums(otu_table(SLL))[rowSums(otu_table(SLL)) != 0]), SLL )
SLL_with_tree <- prune_taxa( names(rowSums(otu_table(SLL_with_tree))[rowSums(otu_table(SLL_with_tree)) != 0]), SLL_with_tree )

# include columns for the various diversity measures ####
rich <- estimate_richness(prune_taxa(taxa_sums(SLL)>0, SLL))#, measures='Simpson')
sample_data(SLL)[,'Div.Observed'] <- rich$Observed
sample_data(SLL)[,'Div.Chao1'] <- rich$Chao1
sample_data(SLL)[,'Div.ACE'] <- rich$ACE
sample_data(SLL)[,'Div.Shannon'] <- rich$Shannon
sample_data(SLL)[,'Div.Simpson'] <- rich$Simpson
sample_data(SLL)[,'Div.InvSimpson'] <- rich$InvSimpson
sample_data(SLL)[,'Div.Fisher'] <- rich$Fisher

# include column for number of OTUs identified
check0 <- function(x) {
  return( sum( x != 0 ) )
}
num.otus <- apply(otu_table(SLL), 2, check0)
sample_data(SLL)[,'Num_OTUs'] <- num.otus




# include column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
library(picante)
faith <- pd(t(as.data.frame(otu_table(SLL_with_tree))), phy_tree(SLL_with_tree))#, include.root = F)
sample_data(SLL)[,'Faiths.PD'] <- faith$PD
sample_data(SLL)[,'Species_Richness'] <- faith$SR




# ******************************************************************** #
## Calculate Unifrac distances ####

# if necessary to calculate these values again, simply uncomment this section. 
# But takes a long time to calculate, better to read from files that were written from original calculation
# library(foreach)
# library(doParallel)
# 
# # This is required in order to register a parallel "backend" so the calculation can be run in parallel
# registerDoParallel(cores = 10)
# 
# weighted_Unifrac <- UniFrac(SLL_with_tree, weighted = T, parallel = T)
# unweighted_Unifrac <- UniFrac(SLL_with_tree, weighted = F, parallel = T)
# 
# write.csv(as.matrix(weighted_Unifrac), file = sprintf("%s/Part_1/SLL1_w_unifrac.csv", home_dir))
# write.csv(as.matrix(unweighted_Unifrac), file = sprintf("%s/Part_1/SLL1_uw_unifrac.csv", home_dir))

weighted_Unifrac <- read.csv(sprintf("%s/Part_1/SLL1_w_unifrac.csv", home_dir), row.names = 1)
diag(weighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
# weighted_Unifrac <- as.dist(weighted_Unifrac)
unweighted_Unifrac <- read.csv(sprintf("%s/Part_1/SLL1_uw_unifrac.csv", home_dir), row.names = 1)
diag(unweighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
# unweighted_Unifrac <- as.dist(unweighted_Unifrac)

sample_data(SLL)[ rownames(weighted_Unifrac), "Weighted_Unifrac"] <- rowMeans(weighted_Unifrac, na.rm = T)
sample_data(SLL)[ rownames(unweighted_Unifrac), "Unweighted_Unifrac"] <- rowMeans(unweighted_Unifrac, na.rm = T)
# ******************************************************************** #









# ********************************************************************************************************* #
###### Add column for age, based on date of birth ######
# ********************************************************************************************************* #
# using the date March 10th 2015 since samples were collected between February 3 and April 14 of 2015, so this is the middle date
# just a simpler method than using the exact dates of collections to get the ages in years. 
ages <- sapply(as.Date(as.matrix(sample_data(SLL)[, "Q1"])), function(x) as.numeric( difftime("2015-03-10", x, units = "days") / 365) )
# give those with ages less than 1 year as NA because these were simply typing errors by the participants that put 2014 or 2015 as their birthday
ages[ages < 1] <- NA

sample_data(SLL)[ , "Age"] <- ages




# ********************************************************************************************************* #
###### Add information about fungus in those samples for which we have data ######
# ********************************************************************************************************* #

fungus.data <- read.delim(sprintf("%s/Part_1/Fungus_data/MALDI-results_modified.csv", home_dir), header=T)
colnames(fungus.data)[2] <- "Plate.well"
colnames(fungus.data)[3] <- "Fungus"
fungus.data <- fungus.data[,2:ncol(fungus.data)]

well_IDs <- read.delim(sprintf("%s/all_sample_96well_locations.csv", home_dir), header=T)
well_IDs$Sample_ID <- gsub("-", "", well_IDs$Sample_ID)
#must flip well name to match fungus.data
well_IDs$Well <- paste0( substring(well_IDs$Well,nchar(as.character(well_IDs$Well))), substring(well_IDs$Well, 1, nchar(as.character(well_IDs$Well))-1) )
rownames(well_IDs) <- paste(well_IDs$Plate, well_IDs$Well, sep="/")

well_IDs[,"Fungi"] <- vector("numeric", nrow(well_IDs))

for (i in 1:nrow(fungus.data)){
  x <- as.character(fungus.data[i, "Plate.well"])
  if (grepl("[.]", x)) {
    # those samples that had colonies of multiple (potentially) of different fungal species
    if (well_IDs[ strsplit(x, "[.]")[[1]][1], "Fungi" ] == 0) {
      well_IDs[ strsplit(x, "[.]")[[1]][1], "Fungi" ] <- as.character(fungus.data[fungus.data[,1]==x, "Fungus"])
    } else {
      well_IDs[ strsplit(x, "[.]")[[1]][1], "Fungi" ] <- paste(well_IDs[ strsplit(x, "[.]")[[1]][1], "Fungi" ], as.character(fungus.data[fungus.data[,1]==x, "Fungus"]), sep=",")
    }

  } else {
    # those samples with colonies of only one supposed fungal species
    well_IDs[ strsplit(x, "[.]")[[1]][1], "Fungi" ] <- as.character(fungus.data[fungus.data[,1]==x, "Fungus"])
  }
}

# To exclude those 'test' wells that are all labeled "C+", and get just those that are kept in sample_data(SLL)
proper_well_IDs <- well_IDs[substr(well_IDs[,"Sample_ID"], 1,3)=="SLL","Sample_ID"]
proper_well_IDs <- proper_well_IDs[proper_well_IDs %in% rownames(sample_data(SLL))]

# Add column to sample_data of fungal species present in samples with data available
# ***Some have more than 1 species, separated by ","
#    If no data available for a sample, value is "0"
sample_data(SLL)[ proper_well_IDs, "Fungi" ] <- well_IDs[ well_IDs[,"Sample_ID"] %in% proper_well_IDs, "Fungi" ]

# ************* #
# Determine which samples have which fungi present
fungi.samples <- rownames(sample_data(SLL)[ sample_data(SLL)[, "Fungi"] != 0])

fungi_entries <- as.matrix(unique(sample_data(SLL)[fungi.samples, "Fungi"]))
diff_fungi <- vector()#mode = "character")

# make a vector of all the different fungi that appear in the sample_data
i = 1
for (f in fungi_entries) {
  if (grepl(',', f)) {
    f1 <- strsplit(f, ',')[[1]][1]
    f2 <- strsplit(f, ',')[[1]][2]
    if (! f1 %in% diff_fungi) {
      diff_fungi[i] <- f1
      i <- i + 1
    }
    if (! f2 %in% diff_fungi) {
      diff_fungi[i] <- f2
      i <- i + 1
    }
  } else {
    if (! f %in% diff_fungi) {
      diff_fungi[i] <- f
      i <- i + 1
    }
  }
}
fungi.list <- vector(mode="list", length=length(diff_fungi))
names(fungi.list) <- diff_fungi


for (s in fungi.samples) {
  if (grepl(',', sample_data(SLL)[s, "Fungi"] )) {
    f1 <- strsplit( as.character( sample_data(SLL)[s, "Fungi"] ), ',' )[[1]][1]
    f2 <- strsplit( as.character( sample_data(SLL)[s, "Fungi"] ), ',' )[[1]][2]
    fungi.list[[ f1 ]] <- c(fungi.list[[ f1 ]], s)
    fungi.list[[ f2 ]] <- c(fungi.list[[ f2 ]], s)
  } else {
    f <- as.character( sample_data(SLL)[s, "Fungi"] )
    fungi.list[[ f ]] <- c(fungi.list[[ f ]], s)
  }
}

# ************* #




# filter out OTUs with low counts
# wh0 = genefilter_sample(SLL, filterfun_sample(function(x) x > 0), A = 0.1 * nsamples(SLL))
# SLL1 = prune_taxa(wh0, SLL)

#normalize counts within each sample
# SLL1 = transform_sample_counts(SLL1, function(x) 100 * x/sum(x))
SLL1 = transform_sample_counts(SLL, function(x) 100 * x/sum(x))


# write tax_table to file
write.csv(tax_table(SLL1), sprintf("%s/Part_1/tax_table_SLL1.csv", home_dir))

# make file with list of sample names:
write(sort(sample_names(SLL1)), file = sprintf("%s/Part_1/sample_names.txt", home_dir))
s700 <- sort(sample_names(SLL1))[1:700]
s619 <- sort(sample_names(SLL1))[701:1319]
write(s700, file = sprintf("%s/Part_1/sample_names_s1-s700.txt", home_dir))
write(s619, file = sprintf("%s/Part_1/sample_names_s701-s1319.txt", home_dir))







# ******************************************************************** #
# Bray-Curtis and Canberra distances ####

library(vegan)

# Bray-Curtis distance combines absolute differences between features (more sensitive to a few large changes)
bray <- as.matrix( vegdist(t(otu_table(SLL1)), method = "bray") )
diag(bray) <- NA
# Canberra distance weights all differences equally (more sensitive to many small changes)
canberra <- as.matrix( vegdist(t(otu_table(SLL1)), method = "canberra") )
diag(canberra) <- NA

sample_data(SLL)[ rownames(bray), "Bray-Curtis"] <- rowMeans(bray, na.rm = T)
sample_data(SLL1)[ rownames(bray), "Bray-Curtis"] <- rowMeans(bray, na.rm = T)
sample_data(SLL)[ rownames(canberra), "Canberra"] <- rowMeans(canberra, na.rm = T)
sample_data(SLL1)[ rownames(canberra), "Canberra"] <- rowMeans(canberra, na.rm = T)


# ******************************************************************** #








# ## Add values for water hardness in Spain ###
# # values taken from here:
# #    http://www.sistemagua.com/informacion-sobre-el-agua/dureza-del-agua-en-espana/
# #    http://naturatips.com/agua/dureza-agua-provincias-espana/
# #    http://naturatips.com/agua/agua-en-espana/
# 
# dureza_by_prov <- list( "Almeria"=60,'Cádiz'=30,'Cordoba'=15,'Granada'=20,'Huelva'=35,'Jaen'=60,'Málaga'=55,'Sevilla'=25, # Andalucia
#                         'Huesca'=25,'Teruel'=30,'Zaragoza'=45, # Aragon
#                         'Asturias'=35, # Asturias (Oviedo)
#                         'Islas Baleares'=40,'Palma de Mallorca'=90, # Islas Baleares
#                         'Las Palmas de Gran Canaria'=25, 'Santa Cruz de Tenerife'=30, # Canarias
#                         'Cantabria'=10, # Cantabria (Santander)
#                         'Albacete'=48, 'Ciudad Real'=60, 'Cuenca'=40, 'Guadalajara'=35, 'Toledo'=55, # Castilla-La Mancha
#                         'Avila'=5,'Burgos'=10,'Leon'=10,'Palencia'=10,'Salamanca'=10,'Segovia'=5,'Soria'=25,'Valladolid'=20,'Zamora'=10, # Castilla y Leon
#                         'Barcelona'=20,'Barcelona_city'=60,'Girona'=20,'Lleida'=25,'Tarragona'=40, # Catalunya
#                         'Alicante'=35,'Alicante_city'=65,'Castellon'=55,'Valencia'=50,  # Comunidad Valenciana
#                         'Badajoz'=25, 'Cáceres'=10, # Extremadura
#                         'A Coruna'=10,'Lugo'=10,'Ourense'=10,'Pontevedra'=10, # Galicia ****No value for Pontevedra but every site says its the same****
#                         'La Rioja'=35, # La Rioja (Logroño)
#                         'Madrid'=5,
#                         'Murcia'=25,'Murcia_city'=60,
#                         'Navarra'=18, # Navarra (Pamplona)
#                         'Araba/Alava'=20,'Bizkaia'=20,'Gipuzkoa'=10, # Pais Vasco
#                         'Ceuta'=35, 'Melilla'=40 )
# 
# durezas <- list('soft'=0:10, 'intermediate'=11:25, 'hard'=26:40, 'very_hard'=41:100)
# 
# dureza_cities <- list( "Barcelona_city" = c(27,28, 41), "Barcelona" = c(26, 34), #Barcelona, Roda de Ter, Calella
#                        "Bizkaia" = c(11,17), #Santurtzi
#                        "Navarra" = c(32), #Altsasu
#                        "Madrid" = c(5,20,21, 33, 40), #Madrid, Cercedilla, Soto del Real
#                        "Málaga" = c(16,24, 31), "Sevilla" = c(6,18, 30), #Malaga, Torrox, Sevilla, El Saucejo
#                        "Palma de Mallorca" = c(29,39), "Islas Baleares" = c(19, 22), #Palma de Mallorca, Santanyí, Port d'Alcúdia
#                        "Murcia_city" = c(3,12), "Murcia" = c(25, 35), #Murcia_city, Bullas, La Paca
#                        "Valencia" = c(4,15, 10), "Castellon" = c(2,7, 8), #Valencia, Tavernes de Valldigna, L'Alcora, Nules
#                        "Pontevedra" = c(9,36), "Ourense" = c(37), #Vigo, Bande 
#                        "Burgos" = c(38), #Villarcayo
#                        "Zaragoza" = c(1,14, 13,23) ) #Zaragoza, Tauste
# 
# # add a column to sample_data for dureza, based on school ID. 
# # (27,28,41 will have the 'Barcelona_city' value while 34,26 will have the 'Barcelona' value)
# for (city in names(dureza_cities)) {
#   IDs <- dureza_cities[[city]]
#   samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% IDs, ])
#   
#   prov_dureza <- unlist( dureza_by_prov[city] )
#   for (d in names(durezas)) {
#     if (prov_dureza %in% durezas[d][[1]]) {
#       dureza_group <- d
#     }
#   }
#   
#   sample_data(SLL)[samps,"Water_hardness_group"] <- dureza_group
#   sample_data(SLL1)[samps,"Water_hardness_group"] <- dureza_group
#   
#   sample_data(SLL)[samps,"Water_hardness"] <- prov_dureza
#   sample_data(SLL1)[samps,"Water_hardness"] <- prov_dureza
# }













# Locations by school -----------------------------------------------------

# ******************* #
# list of cities with the school IDs from there
# cities <- list( Barcelona = c(26,27,28,34,41),
#                 Bilbao = c(11,17,32),
#                 Madrid = c(5,20,21,33,40),
#                 Malaga = c(16,24,31),
#                 Mallorca = c(19,22,29,39),
#                 Murcia = c(3,12,25,35),
#                 Sevilla = c(6,18,30),
#                 Valencia = c(2,4,7,8,10,15),
#                 Vigo = c(9,36,37,38),
#                 Zaragoza = c(1,13,14,23) )

comunidades <- list( "Andalucía" = c(16,24,31, 6,18,30), # Malaga, Torrox, El_Saucejo, Sevilla
                     "Aragón" = c(1,14, 13, 23), # Zaragoza, Tauste, Ejea de los Caballeros
                     "Castilla y León" = c(38), # Villarcayo
                     "Cataluña" = c(27,28, 41, 34, 26), # Barcelona, Hospitalet de Llorbregat, Calella, Roda_de_Ter
                     "Galicia" = c(9,36, 37), # Vigo, Bande
                     "Islas Baleares" = c(29,39, 19, 22), # Palma, Santanyi
                     "Comunidad de Madrid" = c(5,20, 21, 40, 33), # Madrid, Pozuelo de Alarcon, Soto_del_Real, Cercedilla
                     "Región de Murcia" = c(3,12, 25, 35), # Murcia, Bullas, La_Paca
                     "Comunidad Foral de Navarra" = c(32), # Altsasu
                     "País Vasco" = c(11,17), # Santurtzi
                     "Comunidad Valenciana" = c(4,15, 8, 2,7, 10)) # Valencia, Nules, L_Alcora, Tavernes_de_Valldigna

provinces <- list( "Barcelona" = c(27,28, 41, 34, 26),
                   "Baleares" = c(29,39, 19, 22),
                   "Castellón" = c(8, 2,7), "Valencia" = c(4,15, 10),
                   "Murcia" = c(3,12, 25, 35),
                   "Málaga" = c(16,24, 31), "Sevilla" = c(6,18, 30),
                   "Madrid" = c(5,20, 21, 40, 33),
                   "Pontevedra" = c(9,36), "Ourense" = c(37),
                   "Burgos" = c(38),
                   "Vizcaya" = c(11,17),
                   "Navarra" = c(32),
                   "Zaragoza" = c(1,14, 13, 23))

cities  <- list( "Barcelona" = c(27,28), "Hospitalet de Llobregat" = c(41), "Calella" = c(34), "Roda de Ter" = c(26),
                 "Santurtzi" = c(11,17),
                 "Altsasu" = c(32),
                 "Madrid" = c(5,20), "Pozuelo de Alarcon" = c(21), "Soto del Real" = c(40), "Cercedilla" = c(33),
                 "Málaga" = c(16,24), "Torrox" = c(31), "Sevilla" = c(6,18), "El Saucejo" = c(30),
                 "Palma de Mallorca" = c(29,39), "Santanyí" = c(19), "Port d'Alcúdia" = c(22),
                 "Murcia" = c(3,12), "Bullas" = c(25), "La Paca" = c(35),
                 "Valencia" = c(4,15), "Nules" = c(8), "L'Alcora" = c(2,7), "Tavernes de Valldigna" = c(10),
                 "Vigo" = c(9,36), "Bande" = c(37),
                 "Villarcayo" = c(38),
                 "Zaragoza" = c(1,14), "Tauste" = c(13), "Ejea de los Caballeros" = c(23) )

lat_lon <- list( "41.39 N, 2.17 E" = c(27,28), "41.37 N, 2.12 E" = c(41), "41.61 N, 2.65 E" = c(34), "41.98 N, 2.31 E" = c(26),
                 "43.33 N, 3.03 W" = c(11,17),
                 "42.90 N, 2.17 W" = c(32),
                 "40.42 N, 3.70 W" = c(5,20), "40.45 N, 3.81 W" = c(21), "40.75 N, 3.78 W" = c(40), "40.74 N, 4.06 W" = c(33),
                 "36.72 N, 4.42 W" = c(16,24), "36.75 N, 3.95 W" = c(31), "37.39 N, 5.98 W" = c(6,18), "37.07 N, 5.10 W" = c(30),
                 "39.57 N, 2.65 E" = c(29,39), "39.34 N, 3.12 E" = c(19), "39.84 N, 3.13 E" = c(22),
                 "37.99 N, 1.13 W" = c(3,12), "38.04 N, 1.67 W" = c(25), "37.63 N, 1.97 W" = c(35),
                 "39.47 N, 0.38 W" = c(4,15), "39.85 N, 0.16 W" = c(8), "40.07 N, 0.21 W" = c(2,7), "39.08 N, 0.25 W" = c(10),
                 "42.24 N, 8.72 W" = c(9,36), "42.03 N, 7.97 W" = c(37),
                 "42.94 N, 3.57 W" = c(38),
                 "41.65 N, 0.89 W" = c(1,14), "41.92 N, 1.26 W" = c(13), "42.13 N, 1.14 W" = c(23) )

# Add city, province and community to columns in sample_data
for (city in names(cities)) {
  school_IDs <- cities[[city]]
  samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% school_IDs, ])
  
  sample_data(SLL)[ samps, "City" ] <- city
  sample_data(SLL1)[ samps, "City" ] <- city
}

for (prov in names(provinces)) {
  school_IDs <- provinces[[prov]]
  samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% school_IDs, ])
  
  sample_data(SLL)[ samps, "Province" ] <- prov
  sample_data(SLL1)[ samps, "Province" ] <- prov
}

for (comm in names(comunidades)) {
  school_IDs <- comunidades[[comm]]
  samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% school_IDs, ])
  
  sample_data(SLL)[ samps, "Community" ] <- comm
  sample_data(SLL1)[ samps, "Community" ] <- comm
}



# ********************************************************************************************************* #
###### Add various water values for each city ######
# ********************************************************************************************************* #

water_data <- read.delim(sprintf("%s/Water_quality_data/water_qualities_by_city.csv", home_dir))
water_data <- water_data[1:30,]
rownames(water_data) <- water_data[,1]
rownames(water_data) <- gsub("Port d’Alcúdia", "Port d'Alcúdia", rownames(water_data))
rownames(water_data) <- gsub("L’Alcora", "L'Alcora", rownames(water_data))
water_data <- water_data[ , 2:22] # no need to include columns Taste, Smell, Color, CO2, SH2, Temperature

cont_water_data <- c("Conductivity","water_pH","Dry.residue.at.180.C","Dry.residue.at.110.C","Cl","F","HCO3","NO3",
                     "SO4","Na","K","Li","Ca","Mg","Sr","Hardness","Alcalinity")#,"CO3"
group_water_data <- c("Mineralization", "Composition", "Hardness_category")

for (city in rownames(water_data)) {
  print(city)
  school_IDs <- cities[[city]]
  samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% school_IDs, ])
  
  for (col in colnames(water_data)) {
    sample_data(SLL)[ samps, col ] <- water_data[ city, col ]
    sample_data(SLL1)[ samps, col ] <- water_data[ city, col ]
  }
}





# ********************************************************************************************************* #
###### Add population size for sample based on its city ######
# ********************************************************************************************************* #

city_pops <- list( "Barcelona" = 1605000, "Hospitalet de Llobregat" = 252171, "Calella" = 18226, "Roda de Ter" = 6122,
                   "Santurtzi" = 46284, 
                   "Altsasu" = 7490,
                   "Madrid" = 3142000, "Pozuelo de Alarcon" = 84558, "Soto del Real" = 8456, "Cercedilla" = 6781,
                   "Málaga" = 569130, "Torrox" = 15117, "Sevilla" = 693878, "El Saucejo" = 4399,
                   "Palma de Mallorca" = 400578, "Santanyí" = 11316, "Port d'Alcúdia" = 19763,
                   "Murcia" = 439889, "Bullas" = 11753, "La Paca" = 1318,
                   "Valencia" = 786189, "Nules" = 13442, "L'Alcora" = 10591, "Tavernes de Valldigna" = 17734,
                   "Vigo" = 294098, "Bande" = 1735, 
                   "Villarcayo" = 4826,
                   "Zaragoza" = 664953, "Tauste" = 6941, "Ejea de los Caballeros" = 16754 )

for (city in names(cities)) {
  IDs <- cities[[city]]
  samps <- rownames(sample_data(SLL1)[sample_data(SLL1)$Qsamples.1 %in% IDs, ])
  
  pop <- city_pops[[city]]
  
  sample_data(SLL)[samps,"Population"] <- pop
  sample_data(SLL1)[samps,"Population"] <- pop
}
# ********************************************************************************************************* #



# # ********************************************************************************************************* #
# #### Add columns for Ratios of common genera ####
# sample_data(SLL)[,"Neisseria_Prevotella"] <- as.numeric(otu_table(SLL1)["Neisseria",] / otu_table(SLL1)["Prevotella",])
# sample_data(SLL1)[,"Neisseria_Prevotella"] <- as.numeric(otu_table(SLL1)["Neisseria",] / otu_table(SLL1)["Prevotella",])
# 
# sample_data(SLL)[,"Neisseria_Veillonella"] <- as.numeric(otu_table(SLL1)["Neisseria",] / otu_table(SLL1)["Veillonella",])
# sample_data(SLL1)[,"Neisseria_Veillonella"] <- as.numeric(otu_table(SLL1)["Neisseria",] / otu_table(SLL1)["Veillonella",])
# 
# sample_data(SLL)[,"Haemophilus_Prevotella"] <- as.numeric(otu_table(SLL1)["Haemophilus",] / otu_table(SLL1)["Prevotella",])
# sample_data(SLL1)[,"Haemophilus_Prevotella"] <- as.numeric(otu_table(SLL1)["Haemophilus",] / otu_table(SLL1)["Prevotella",])
# 
# sample_data(SLL)[,"Haemophilus_Veillonella"] <- as.numeric(otu_table(SLL1)["Haemophilus",] / otu_table(SLL1)["Veillonella",])
# sample_data(SLL1)[,"Haemophilus_Veillonella"] <- as.numeric(otu_table(SLL1)["Haemophilus",] / otu_table(SLL1)["Veillonella",])
# 
# sample_data(SLL)[,"Streptococcus_Prevotella"] <- as.numeric(otu_table(SLL1)["Streptococcus",] / otu_table(SLL1)["Prevotella",])
# sample_data(SLL1)[,"Streptococcus_Prevotella"] <- as.numeric(otu_table(SLL1)["Streptococcus",] / otu_table(SLL1)["Prevotella",])
# 
# # this one they did in the Acharya paper
# sample_data(SLL)[,"Prevotella_Bacteroides"] <- as.numeric(otu_table(SLL1)["Prevotella",] / otu_table(SLL1)["Bacteroides",])
# sample_data(SLL1)[,"Prevotella_Bacteroides"] <- as.numeric(otu_table(SLL1)["Prevotella",] / otu_table(SLL1)["Bacteroides",])
# 
# # in case any genera were not present in some samples (as with Bacteroides), make value NA
# ratios <- c("Neisseria_Prevotella","Neisseria_Veillonella","Haemophilus_Prevotella","Haemophilus_Veillonella",
#             "Streptococcus_Prevotella","Prevotella_Bacteroides")
# for (r in ratios) {
#   sample_data(SLL)[ ! is.finite(as.numeric(as.matrix(sample_data(SLL))[,r])), r ] <- NA
#   sample_data(SLL1)[ ! is.finite(as.numeric(as.matrix(sample_data(SLL1))[,r])), r ] <- NA
# }
# # ********************************************************************************************************* #




# ********************************************************************************************************* #
###### Make biom-format object ######
# ********************************************************************************************************* #
# library(biomformat)
# 
# SLL.biom <- biomformat::make_biom(otu_table(SLL))
# biomformat::write_biom(SLL.biom, "/users/tg/jwillis/SLL/Part_1/PICRUSt/otu_table.biom")




# ********************************************************************************************************* #
region <- comunidades
# region <- provinces
# region <- cities

reg_name <- 'Comunidad Valenciana'
region_samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)$Qsamples.1 %in% region[reg_name][[1]], ])
SLL1.region <- prune_samples( region_samples, SLL1 )

region.top10 <- names(sort(taxa_sums(SLL1.region), decreasing = T))[1:10]
SLL1.top10 <- as.data.frame(otu_table(SLL1.region))[region.top10,]
# ******************* #

# ******************* #
# To look at the different types of water drunk at home
agua <- "Del Grifo (Filtrada)"
# agua <- "Del Grifo (No Filtrada)"
# agua <- "Embotellada"
# agua <- "No Tratada (Fuente, Pozo O Rio)"

SLL1.agua <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q27 == agua ]), SLL1)
agua.top10 <- names(sort(taxa_sums(SLL1.agua), decreasing = T))[1:10]
SLL1.top10 <- as.data.frame(otu_table(SLL1.agua))[agua.top10,]
# ******************* #


# top 10 most common OTUs
top10 <- names(sort(taxa_sums(SLL1), decreasing = T))[1:10]
SLL1.top10 <- as.data.frame(otu_table(SLL1))[top10,]


# ****************************************************************************************************************** #

# for whether or not to include a region name in plot titles:
# region_type <- "region"
region_type <- "all"

#for whether or not to include a water type in plot titles:
# agua_types <- "agua"
agua_types <- "all"








# ****************************************************************************************************************** #
# Get otu_tables at different taxa levels ####
# ****************************************************************************************************************** #
# tlev <- "Phylum"
# taxa <- unique(tax_table(SLL1)[,tlev])

t0 <- Sys.time()

tlev_otus <- list()
tlev_otus_rel <- list()

for (tlev in c("Family","Order","Class","Phylum")) {
  print(tlev)
  taxa <- unique(tax_table(SLL1)[,tlev])
  # get counts of all genera from otu_table(SLL) that are within each value of the given tlev
  tlev_otus[[tlev]] <- apply(otu_table(SLL), 2, function(x) sapply(taxa, function(y) sum(x[ rownames(tax_table(SLL)[ tax_table(SLL)[ , tlev ]==y, ]) ]) ) )
  rownames(tlev_otus[[tlev]]) <- taxa
  # get normalized values -- relative abundances
  tlev_otus_rel[[tlev]] <- apply(tlev_otus[[tlev]], 2, function(x) 100 * x/sum(x))
}
# # get counts of all genera from otu_table(SLL) that are within each value of the given tlev
# tlev_otus <- apply(otu_table(SLL), 2, function(x) sapply(taxa, function(y) sum(x[ rownames(tax_table(SLL)[ tax_table(SLL)[ , tlev ]==y, ]) ]) ) )
# rownames(tlev_otus) <- taxa
# 
# # get normalized values -- relative abundances
# tlev_otus_rel <- apply(tlev_otus, 2, function(x) 100 * x/sum(x))
print(Sys.time() - t0)










# ****************************************************************************************************************** #
# Get table of distribution of genera within each phylum for a given region (or overall) ####
# ****************************************************************************************************************** #
# tlev <- "Phylum"
# taxa <- unique(tax_table(SLL1)[,tlev])
uniq_gen <- unique(tax_table(SLL1)[,"Genus"])
uniq_phy <- unique(tax_table(SLL1)[,"Phylum"])

community <- "all"
# community <- "Andalucia"
# community <- as.character(as.matrix(unique(sample_data(SLL1)[,"Q14.2"])))[1]

if (community == "all") {
  which_samples <- rownames(sample_data(SLL1))
} else {
  which_samples <- rownames(sample_data(SLL1)[sample_data(SLL1)[, "Q14.2"] == community,])
}
community_otus <- otu_table(SLL1)[, which_samples]


gen_otu_table <- matrix(NA, nrow = length(uniq_gen), ncol = ncol(community_otus))
rownames(gen_otu_table) <- as.character(uniq_gen)
colnames(gen_otu_table) <- which_samples
for (ta in uniq_gen) {
  ta_genera <- rownames( tax_table(SLL1)[tax_table(SLL1)[ , "Genus"] == ta, ] )
  gen_otu_table[ta, ] <- colSums(community_otus[ta_genera, ])
}

phy_otu_table <- matrix(NA, nrow = length(uniq_phy), ncol = ncol(community_otus))
rownames(phy_otu_table) <- as.character(uniq_phy)
colnames(phy_otu_table) <- which_samples
for (ta in uniq_phy) {
  ta_genera <- rownames( tax_table(SLL1)[tax_table(SLL1)[ , "Phylum"] == ta, ] )
  phy_otu_table[ta, ] <- colSums(community_otus[ta_genera, ])
}



genus_abunds <- sort(rowMeans(gen_otu_table)[rowMeans(gen_otu_table)>0], decreasing = T)
phylum_abunds <- sort(rowMeans(phy_otu_table)[rowMeans(phy_otu_table)>0], decreasing = T)

# phylum_abunds <- unique(tax_table(SLL1)[tax_table(SLL1)[ , tlev] %in% names(genus_abunds), "Phylum"])

phy_list <- list()
for (p in names(phylum_abunds)) {
  # get all genera from this community with phylum "p"
  gen_in_p <- rownames(tax_table(SLL1)[ tax_table(SLL1)[,"Phylum"]==p & tax_table(SLL1)[,"Genus"] %in% names(genus_abunds), "Genus" ])
  # get only those that are not unclassified
  not_unclassified <- tryCatch( rownames(tax_table(SLL1)[tax_table(SLL1)[,"Genus"] %in% gen_in_p, "Phylum"]), error = function(x) NULL)
  # for those phyla that only had unclassified genera in this community
  if (is.null(not_unclassified)) {next}
  # then order the genera by abundance within this community
  not_unclassified <- sort(genus_abunds[not_unclassified], decreasing = T)
  # make elements for both the names and the abundances
  phy_list[[p]] <- names(not_unclassified)
  # round values for the sake of visualization
  phy_list[[ paste(p,'Abundance',sep = ' ')]] <- ifelse(not_unclassified<0.001, '< 0.001', round( not_unclassified, digits = 3))
  
}

phy_frame <- matrix('', nrow=max(sapply(phy_list, length)), ncol=length(phy_list))
colnames(phy_frame) <- names(phy_list)
for (p in names(phy_list)) {
  phy_frame[1:length(phy_list[[p]]), p] <- phy_list[[p]]
}

write.csv(phy_frame, sprintf("%s/Part_1/figures/Abundances/pie_charts/Abundance_tables/%s_abundances.csv",
                             home_dir,community), row.names = F)






# ******************************************************* #
# Pie charts / Donut charts of abundances ####
# ******************************************************* #
# # table_of_interest <- otu_table(SLL1)
# table_of_interest <- gen_otu_table
# 
# top10 <- names(sort(rowSums(table_of_interest), decreasing = T))[1:5]
# nottop10 <- names(sort(rowSums(table_of_interest), decreasing = T))[6:nrow(table_of_interest)]
# 
# SLL1.pie <- as.data.frame(table_of_interest)[top10,]
# SLL1.pie["Other",] <- colSums(table_of_interest[nottop10,])
# 
# 
# 
# #mean and sd of given taxa per sample
# tax.mps <- rowMeans(SLL1.pie)
# tax.sd <- rowSds(as.matrix(SLL1.pie))
# 
# #create frame that can be used for bar plots
# tax.frame <- as.data.frame(tax.mps)
# tax.frame[,2] <- tax.sd
# tax.frame[,3] <- rownames(SLL1.pie)
# colnames(tax.frame) <- c("mps","sd",tlev)
# 
# 
# #plot
# # ggplot(tax.frame, aes(x=factor(1), y=mps, fill=factor(Genus))) +
# ggplot(tax.frame, aes(x=factor(1), y=mps, fill=factor(reorder(tax.frame[,tlev], -mps)))) +
#   geom_bar(width=1, stat="identity") + coord_polar(theta = "y", start = 1.5708)+#, direction = -1) +
#   scale_fill_brewer(palette = "Spectral") + theme_void() +
#   labs(fill=tlev) + ggtitle(community) + theme(plot.title = element_text(hjust=0.5))

# **************************** # 
##### From the SO answer here: https://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot

#' x      numeric vector for each slice
#' group  vector identifying the group for each slice
#' labels vector of labels for individual slices
#' col    colors for each group
#' radius radius for inner and outer pie (usually in [0,1])

donuts <- function(x, group = 1, labels_out = NA, labels_in = NA, col = NULL, ti = NULL, radius = c(.7, 1)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]
  
  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    # rev(Vectorize(adjustcolor)(x, alpha.f = al))
    Vectorize(adjustcolor)(x, alpha.f = al)
  })
  
  plot.new()
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L], #init.angle = 30,
      col = unlist(col.sub), labels = labels_out, cex=1.75,
      main = ti, cex.main = 2)
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L], #init.angle = 30,
      col = unlist(col.main), labels = NA)#labels_in, cex=0.6)
  
  legend("topleft", legend = labels_in, col = names(col.main), pch = 15, bty = 'n', cex = 1.9)
  # text(x = c(0.1, -.25, 0, .35, .5, 0.55), y = c(0.15, 0, -.3, -.25, -.1, -0.025), 
  #      labels = unique(phyla_for_pie$phylum), col = 'white', cex = 1.2)
}
# **************************** # 
# **************************** # 

community <- "all"
# community <- "Islas Baleares"
# community <- as.character(as.matrix(unique(sample_data(SLL1)[,"Q14.2"])))[1]

tlev <- "Phylum"

if (community == "all") {
  which_samples <- rownames(sample_data(SLL1))
  pie_title <- sprintf("5 most common %s members and their most common Genera", tlev)
  # pie_title <- sprintf("5 Phyla más comunes y sus Géneros más comunes")
} else {
  which_samples <- rownames(sample_data(SLL1)[sample_data(SLL1)[, "Q14.2"] == community,])
  pie_title <- sprintf("5 most common %s members and their most common Genera\n%s", tlev, community)
  # pie_title <- sprintf("5 Phyla más comunes y sus Géneros más comunes\n%s", community)
}



# so need to make a table with 24 rows...
#   top 5 phyla + other...
#     top 3 genera in each phylum + other
topPhyla <- c(names(sort(rowSums(tlev_otus_rel[[tlev]][ , which_samples]), decreasing = T))[1:5], "Other")
topGenera <- list()
for (p in topPhyla) {
  # take top 3 genera per phylum for the top 4 phyla, then only the top 2 for the 5th and "Other" 
  # since they are small and get cluttered on the pie chart
  
  ifelse(p=="Other", {
    gen_in_phy <- rownames(tax_table(SLL1)[ ! tax_table(SLL1)[,tlev] %in% topPhyla[1:5],]);
    # topGenera[[p]] <- c(names(sort(rowSums(otu_table(SLL1)[gen_in_phy, which_samples]), decreasing = T))[1:2], "Other", NA)
    topGenera[[p]] <- c("Other", NA, NA, NA)
  }, {
    gen_in_phy <- rownames(tax_table(SLL1)[tax_table(SLL1)[,tlev]==p,]);
    ifelse(p %in% topPhyla[2:5],
           topGenera[[p]] <- c(names(sort(rowSums(otu_table(SLL1)[gen_in_phy, which_samples]), decreasing = T))[1:2], "Other", NA), 
           topGenera[[p]] <- c(names(sort(rowSums(otu_table(SLL1)[gen_in_phy, which_samples]), decreasing = T))[1:3], "Other"))
  })
}

# get list of abundances
topGen.frame <- as.data.frame(topGenera)
generaAbunds <- list()
for (p in topPhyla) {
  gens <- list()
  for (g in topGen.frame[,p]) {
    ifelse(g == "Other", {
      ifelse(p == "Other", {
        # if genus and phylum are both "Other"
        gen_in_phy <- rownames(tax_table(SLL1)[ ! tax_table(SLL1)[,tlev] %in% topPhyla[1:5],])
        g_in_p_other <- gen_in_phy[ ! gen_in_phy %in% topGenera[["Other"]] ]
        gens[[g]] <- sum(rowMeans(otu_table(SLL1)[g_in_p_other, which_samples]))
      }, {
        # if genus is "Other"
        gen_in_phy <- rownames(tax_table(SLL1)[ tax_table(SLL1)[,tlev] == p])
        g_in_p_other <- gen_in_phy[ ! gen_in_phy %in% topGenera[[p]] ]
        gens[[g]] <- sum(rowMeans(otu_table(SLL1)[g_in_p_other, which_samples]))
      })
    }, {
      # if neither genus nor phylum are "Other"
      gens[[g]] <- rowMeans(otu_table(SLL1)[g, which_samples])
    })
  }
  # then add the vector of genus abundances to each list element (a phylum name)
  generaAbunds[[p]] <- unname(unlist(gens))
}

# to determine structure for donut chart
if (tlev == "Phylum") {
  phylum_struc <- c(1L,1L,1L,1L,2L,2L,2L,3L,3L,3L,4L,4L,4L,5L,5L,5L,6L)
  genus_struc <- c(1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,11L,12L,13L,14L,15L,16L,17L)
  rownames_val <- -17L
}


# make data.frame for plot
phyla_for_pie <- structure(list(phylum = structure(phylum_struc,
                                                    .Label = topPhyla,#rep(topPhyla, each=4),
                                                    class = "factor"),
                                genus = structure(genus_struc,
                                                  .Label = unname(unlist(topGenera))[ ! is.na(unname(unlist(topGenera))) ],
                                                  class = "factor"),
                                abunds = unname(unlist(generaAbunds))),
                           .Names = c("phylum", "genus", "abunds"),
                           row.names = c(NA, rownames_val), 
                           class = "data.frame")
                                
phyla_for_pie$total <- with(phyla_for_pie, ave(abunds, phylum, FUN = sum))

# plot donut
with(phyla_for_pie,
     donuts(abunds, phylum, sprintf('%s: %s%%', genus, round(abunds,2)), unique(phylum),
            col = c('cyan2','red','orange','green','dodgerblue2','grey'),
            ti = pie_title)
)
# ****************************************************************************************************************** #







# top 10 most common OTUs
top10 <- names(sort(taxa_sums(SLL1), decreasing = T))[1:10]
SLL1.top10 <- as.data.frame(otu_table(SLL1))[top10,]


# ****************************************************************************************************************** #
# Bar charts of abundances ####
# ****************************************************************************************************************** #

#mean and sd of given taxa per sample
tax.mps <- rowMeans(SLL1.top10)
tax.sd <- rowSds(as.matrix(SLL1.top10))

#percentage of samples in which given taxa appears
tax.freqs <- round(rowSums(SLL1.top10 != 0) / length(colnames(SLL1.top10)) * 100, digits=2)
tax.freqs <- paste(as.character(tax.freqs),'%',sep='')


#create frame that can be used for bar plots
tax.frame <- as.data.frame(tax.mps)
tax.frame[,2] <- tax.sd
tax.frame[,3] <- tax.freqs
tax.frame[,4] <- rownames(SLL1.top10)
colnames(tax.frame) <- c("mps","sd","freqs","Genus")


if (region_type == "all") {
  if (agua_types == "all") {
    ti <- 'Relative abundance for most common OTU'
  } else if (agua_types == "agua") {
    ti <- sprintf('Relative abundance for most common OTU\nAgua %s', agua)
  }
} else if (region_type == "region") {
  ti <- sprintf('Relative abundance for most common OTU\n%s', reg_name)
}

#plot
ggplot(tax.frame, aes(x=reorder(tax.frame[,"Genus"], -mps), y=mps, 
                      fill=reorder(tax.frame[,"Genus"], -mps))) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mps-sd, ymax=mps+sd), width=0.2) +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
  ggtitle(ti) +
  xlab("Genus") + ylab('Mean normalized abundance per sample (as %)') + scale_fill_hue(name="Genus") +
  geom_text(aes(x=reorder(tax.frame[,"Genus"], -mps), label=freqs), vjust=-1)



### mean abundance for these top 10, and sd
mean(colSums(SLL1.top10))
sd(colSums(SLL1.top10))



# *********************************************************************** #
#frequency at which a given taxa is the max in a sample
maxes <- vector("integer", length(rownames(SLL1.top10)))
# names(maxes) <- rownames(SLL1.top10)
for (n in 1:length(rownames(SLL1.top10))) {
  for (s in colnames(SLL1.top10)) {
    if (SLL1.top10[n,s] == max(SLL1.top10[,s])) {
      maxes[[n]] <- maxes[[n]] + 1
    }
  }
}
maxes <- round(as.numeric(maxes) / length(colnames(SLL1.top10)) *100, digits=2)
# maxes.percent <- paste(as.character(maxes),'%','')

max.frame <- as.data.frame(maxes)
max.frame[,2] <- rownames(SLL1.top10)
max.frame[,3] <- tax.mps
colnames(max.frame) <- c("maxes", "Genus", "Means")


if (region_type == "all") {
  if (agua_types == "all") {
    ti <- 'Frequency as the most abundant OTU in a sample'
  } else if (agua_types == "agua") {
    ti <- sprintf('Frequency as the most abundant OTU in a sample\nAgua %s', agua)
  }
} else if (region_type == "region") {
  ti <- sprintf('Frequency as the most abundant OTU in a sample\n%s', reg_name)
}

# ggplot(max.frame, aes(x=Genus, y=maxes, fill=Genus)) +
ggplot(max.frame, aes(x=reorder(max.frame[,"Genus"], -Means), y=maxes, 
                      fill=reorder(max.frame[,"Genus"], -Means))) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
  ggtitle(ti) +
  xlab("Genus") + ylab('Frequency (%)') + scale_fill_hue(name="Genus")
# *********************************************************************** #





















# ****************************************************************************************************************** #
# Box plots of abundances within samples ####
# ****************************************************************************************************************** #

if (region_type == "all") {
  if (agua_types == "all") {
    ti <- 'Percent of given OTU per sample'
  } else if (agua_types == "agua") {
    ti <- sprintf('Percent of given OTU per sample\nAgua %s', agua)
  }
} else if (region_type == "region") {
  ti <- sprintf('Percent of given OTU per sample\n%s', reg_name)
}

ggplot(melt(t(SLL1.top10)), aes(x=factor(Var2, levels=rownames(SLL1.top10)), y=value, 
                                   fill=factor(Var2, levels=rownames(SLL1.top10)))) +
  geom_boxplot() + theme_minimal() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=14), plot.title = element_text(hjust=0.5)) + 
  # ggtitle(ti) + xlab("Genus") + ylab('% per sample') + scale_fill_hue(name="Genus")
  xlab("Genus") + ylab("") + scale_fill_hue(guide=F) # for the paper
# ****************************************************************************************************************** #





# ****************************************************************************************************************** #
#stacked bars of each of the top10 for each sample
# SLL.df <- as.data.frame(otu_table(SLL1))[rev(top10),]
SLL.df <- SLL1.top10
# rownames(SLL.df) <- rev(rownames(SLL1.top10))

tno <- matrix(0, nrow=length(rownames(SLL.df)), ncol=length(colnames(SLL.df)))
colnames(tno) <- colnames(SLL.df)
rownames(tno) <- rownames(SLL.df)

tax_ranks <- function(x) {
  rank( SLL.df[top10[10], ] ) #orders all bars by relative abundance of Propionibacterium for given sample
}
tno <- apply(SLL.df,1,tax_ranks)
# colnames(tno) <- rev(colnames(tno))
tno <- t(tno)

# tno.m <- cbind( melt(as.matrix(SLL.df[ rev(rownames(SLL.df)), ])), melt(tno[ rev(rownames(tno)), ]) )
tno.m <- cbind( melt(as.matrix(SLL.df)), melt(tno) )
tno.m <- tno.m[,c(1,2,3,6)]
colnames(tno.m)[4] <- 'ranks'


if (region_type == "all") {
  ti <- 'Total normalized abundances per sample for top 10 OTUs'
} else if (region_type == "region") {
  ti <- sprintf('Total normalized abundances per sample for top 10 OTUs\n%s', reg_name)
}

#all top 10 stacked, ordered by values for given taxa
ggplot(tno.m, aes(x=reorder(Var2,-value), y=value, fill=Var1)) +#factor(Var1, levels=top10))) +
  geom_bar(stat='identity') +
  # theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + ggtitle(ti) + 
  theme(axis.title=element_blank(), axis.text.y=element_text(size=14)) + # for paper
  guides(fill = guide_legend(title.theme=element_text(size=15, angle=0, face="bold"), 
                             label.theme=element_text(size=13, angle=0)) ) + # for paper
  xlab('Samples') + ylab('Abundances') + scale_fill_hue(name="Genus")


#all top 10 in separate plots, all ordered by values of most abundant taxa
ggplot(tno.m, aes(x=reorder(Var2,-ranks), y=value, fill=Var1)) +
  geom_bar(stat='identity') +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
  facet_wrap(~Var1) +
  ggtitle(ti) +
  xlab('Samples') + ylab('Abundances') + guides(fill=F) 
# ****************************************************************************************************************** #






# # ****************************************************************************************************************** #
# # Tables of abundances for a given taxonomic level ####
# # ****************************************************************************************************************** #
# 
# level_list <- c(2:6)
# # level_list <- 3
# 
# for (table_level in level_list) {
#   
#   print(sprintf("table_level = %s", table_level))
#   
#   unique_tax <- unique(tax_table(SLL1)[,c(1:table_level)])
#   full_taxa_names <- apply( unique_tax, 1, paste, collapse=';')
#   
#   otu_new_level <- matrix(0, nrow=length(full_taxa_names), ncol=ncol(otu_table(SLL1)))
#   rownames(otu_new_level) <- full_taxa_names
#   colnames(otu_new_level) <- colnames(otu_table(SLL1))
#   
#   
#   for (t_full in full_taxa_names) {
#     #break down full name separating at ';'
#     ta <- unlist(strsplit(t_full, ';'))
#     #makes matrix of character vectors, columns are genera (rownames of tax_table, which include unique identifiers for unclassified genera)
#     x <- apply(tax_table(SLL1)[, c(1:table_level)], 1, as.character)
#     #column will be all TRUEs at genera that match all columns of tax_table up to table_level
#     wheretrue <- (x == ta)
#     
#     #if only one match in wheretrue, wheretrue loses its "matrix" class and becomes "logical" class
#     # logical class then loses the colname (which is the genus name used as unique identifier, required for calculations in colSums() below)
#     # ***This occurs for all "unclassified" entries at table_level, as well as any with only one entry for the preceding level, so this must 
#     #    be done to determine the particular unclassified genus to avoid including all unclassified genera each time "unclassified" is encountered
#     if ( class(wheretrue[,colSums(wheretrue)==table_level]) == "logical" ) {
#       # in this case, must narrow down tax_table to those rows (should be only 1 row) matching ta, take rowname (unique genus)
#       new_table <- tax_table(SLL1)
#       for (i in c(1:table_level)) {
#         new_table <- new_table[new_table[,i] == ta[i], ]
#       }
#       #here seq_rows should be the genus name for the 1 row matching all columns of tax_table up to table_level
#       seq_rows <- rownames(new_table)
#     } else {
#       #for all other OTUs at table_level, will collect those genus names matching all columns of tax_table up to table_level
#       seq_rows <- colnames(wheretrue[,colSums(wheretrue)==table_level])
#     }
#     
#     #sum of abundances per sample for those genera matching all columns of tax_table up to table_level
#     seg <- colSums( otu_table(SLL1)[ seq_rows, ] )
#     otu_new_level[t_full, ] <- seg
#   }
#   
#   
#   taxa_stats <- matrix(0, nrow=length(full_taxa_names), ncol=(nclusters*2)+2)
#   rownames(taxa_stats) <- full_taxa_names
#   colnames(taxa_stats) <- c("mean", "sd", c("S1_mean", "S1_sd", "S2_mean", "S2_sd"))
#   
#   # to obtain stats for samples within each cluster
#   stom1 <- rownames(samples.by.cluster)[samples.by.cluster$otu.cluster=="1"]
#   stom2 <- rownames(samples.by.cluster)[samples.by.cluster$otu.cluster=="2"]
#   
#   taxa_stats[ rownames(otu_new_level), "mean"] <- rowMeans(otu_new_level)
#   taxa_stats[ rownames(otu_new_level), "sd"] <- rowSds(otu_new_level)
#   taxa_stats[ rownames(otu_new_level), "S1_mean"] <- rowMeans(otu_new_level[ , stom1])
#   taxa_stats[ rownames(otu_new_level), "S1_sd"] <- rowSds(otu_new_level[ , stom1])
#   taxa_stats[ rownames(otu_new_level), "S2_mean"] <- rowMeans(otu_new_level[ , stom2])
#   taxa_stats[ rownames(otu_new_level), "S2_sd"] <- rowSds(otu_new_level[ , stom2])
#   
#   sorted_names <- names(rev(sort(taxa_stats[,1])))
#   sorted_taxa_stats <- taxa_stats[sorted_names, ] #sorted by overall mean
#   
#   fname = sprintf("%s/Part_1/with_Pedro/taxa_stats_ta%s.csv", home_dir, table_level)
#   write.csv(sorted_taxa_stats, file=fname)
#   
#   # ********* #
#   # Then make stats about unclassified at each level, to show for how much of the data we can be informative
#   unclass_stats <- matrix(0, nrow=1, ncol=(nclusters*2)+2)
#   rownames(unclass_stats) <- sprintf("ta%s",table_level)
#   colnames(unclass_stats) <- c("mean", "sd", c("S1_mean", "S1_sd", "S2_mean", "S2_sd"))
#   
#   unclassifieds <- rownames(otu_new_level)[grep("unclassified",rownames(otu_new_level))]
#   
#   if (length(unclassifieds) == 1) {
#     #because otu_new_level will no longer be matrix, now is numeric, does not retain columns, cant compute colSums()
#     print(c(table_level, unclassifieds))
#     unclass_stats[1,"mean"] <- mean(otu_new_level[unclassifieds,])
#     unclass_stats[1,"sd"] <- sd(otu_new_level[unclassifieds,])
#     unclass_stats[1,"S1_mean"] <- mean(otu_new_level[unclassifieds, stom1])
#     unclass_stats[1,"S1_sd"] <- sd(otu_new_level[unclassifieds, stom1])
#     unclass_stats[1,"S2_mean"] <- mean(otu_new_level[unclassifieds, stom2])
#     unclass_stats[1,"S2_sd"] <- sd(otu_new_level[unclassifieds, stom2])
#   } else {
#     unclass_stats[1,"mean"] <- mean(colSums(otu_new_level[unclassifieds,]))
#     unclass_stats[1,"sd"] <- sd(colSums(otu_new_level[unclassifieds,]))
#     unclass_stats[1,"S1_mean"] <- mean(colSums(otu_new_level[unclassifieds, stom1]))
#     unclass_stats[1,"S1_sd"] <- sd(colSums(otu_new_level[unclassifieds, stom1]))
#     unclass_stats[1,"S2_mean"] <- mean(colSums(otu_new_level[unclassifieds, stom2]))
#     unclass_stats[1,"S2_sd"] <- sd(colSums(otu_new_level[unclassifieds, stom2]))
#   }
#   
#   uncl_csv <- read.delim(sprintf("%s/Part_1/with_Pedro/taxa_stats_unclassified.csv", home_dir), header=T, row.names=1, sep=',')
#   uncl_csv[rownames(unclass_stats),] <- unclass_stats
#   
#   fname = sprintf("%s/Part_1/with_Pedro/taxa_stats_unclassified.csv", home_dir)
#   write.csv(uncl_csv, file=fname)
# }
# 
# # ****************************************************************************************************************** #

## Average number of OTUs per found per sample
mean(colSums(otu_table(SLL1) != 0))








# ******************************************************* #
# Get phylo object of only those OTUs present in at least 75% of samples ####

tlev <- "Genus"; otutab <- otu_table(SLL)
# tlev <- "Phylum"; otutab <- tlev_otus[[tlev]]

P_A_table <- apply( otutab, 2, function(x) ifelse(x==0, 0, 1) )


min75 <- genefilter_sample(SLL1, filterfun_sample(function(x) x > 0), A = 0.75 * nsamples(SLL1))
SLL.min75 <- prune_taxa(min75, SLL1)


divs <- as.numeric(sample_data(SLL1)[,"Div.Shannon"][[1]])
quints <- quantile(divs, probs=seq(0,1,0.2))
# get samples by diversity quartile
lowdivs <- rownames(sample_data(SLL1)[ sample_data(SLL1)[,"Div.Shannon"] < quints['20%'], ] )

# slllow <- prune_samples(rownames(sample_data(SLL)[sample_data(SLL)[,"Diversity_group_Div.Shannon"]=="Low", ]), SLL)
slllow <- prune_samples(lowdivs, SLL1)
min75_low <- genefilter_sample(SLL1, filterfun_sample(function(x) x > 0), A = 0.75 * nsamples(slllow))
SLL.min75_low <- prune_taxa(min75_low, slllow)


means.min75 <- sort(rowMeans(otu_table(SLL1)[taxa_names(SLL.min75),]), decreasing = T)
pres.min75 <- sort(rowSums(P_A_table[taxa_names(SLL.min75),]) / 1319 * 100, decreasing = T)
tab.min75 <- as.data.frame(means.min75)
tab.min75[,2] <- pres.min75[names(means.min75)]
tab.min75[,c(3:6)] <- tax_table(SLL.min75)[names(means.min75), c(2:5)]
colnames(tab.min75) <- c("Mean relative abundance", "% of samples in which present", "Phylum","Class","Order","Family")
write.csv(tab.min75, sprintf("/users/tg/jwillis/SLL/Part_1/core_OTUs.csv"))

mean(colSums(otu_table(SLL.min75)))
sd(colSums(otu_table(SLL.min75)))
# ******************************************************* #






# ****************************************************************************************************************** #
# ******************************************************* #
# Clustering and network analyses ####
# ******************************************************* #

# Based on the steps in the Bork Enterotype paper

library(cluster)
library(clusterSim)

# Create distance matrix using the Jensen-Shannon distance
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j] = JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix) -> resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

#####  clutsering samples based on all OTUs that came through the original filter
# gender <- "male"
# gender <- "female"
gender <- "both"
# gender <- "CORE"

if (gender == "male") {
  males <- rownames(data.frame(sample_data(SLL)[ sample_data(SLL)[,'Q2'] == 1]))
  ttc <- prune_samples(males, SLL1)
} else if (gender == "female") {
  females <- rownames(data.frame(sample_data(SLL)[ sample_data(SLL)[,'Q2'] == 0]))
  ttc <- prune_samples(females, SLL1)
} else if (gender == "both") {
  # tax.to.clust <- otu_table(SLL1)[pred.core, ]
  ttc <- SLL1
} else if (gender == "CORE") {
  # In Hisayama, they clustered using only those "CORE" OTUs that had an average relative abundance of at least 1%
  core <- taxa_names(SLL.min75)[ rowMeans(otu_table(SLL.min75)) >= 1.0 ]
  ttc <- prune_taxa(core, SLL1)
}
tax.to.clust <- otu_table(ttc)



# Distance matrix object
dist_meas <- "JSD"
# dist_meas <- "Weighted_Unifrac"
# dist_meas <- "Unweighted_Unifrac"

if (dist_meas == "JSD") {
  otu.jsd = dist.JSD(tax.to.clust)
} else if (dist_meas == "Weighted_Unifrac") {
  otu.jsd = as.dist(weighted_Unifrac)
} else if (dist_meas == "Unweighted_Unifrac") {
  otu.jsd = as.dist(unweighted_Unifrac)
}



# Cluster using the Partitioning Around Medoids (PAM) algorithm
pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
  clust = as.vector(pam(as.dist(x), k, diss=TRUE))
  return(clust$clustering)
}

# Determine optimal number of clusters using the Calinski-Harabasz (CH) Index
CH.index = NULL

for (k in 1:10) { 
  if (k == 1) {
    CH.index[k] = NA 
  } else {
    otu.cluster_temp = pam.clustering(otu.jsd, k)
    CH.index[k] = index.G1(t(tax.to.clust), otu.cluster_temp,  d = otu.jsd, centrotypes = "medoids")
  }
}
plot(CH.index, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")


nclusters <- match(max(CH.index, na.rm=T), CH.index) #use # of clusters with max CH index
# nclusters <- 3


otu.cluster = pam.clustering( otu.jsd, k = nclusters )
CH.index.final = index.G1(t(tax.to.clust), otu.cluster, d = otu.jsd, centrotypes = "medoids")

# silhouette measures how similar one sample is to the others in the same cluster versus those in the neighbor cluster
obs.silhouette=mean(silhouette(otu.cluster, otu.jsd)[,3])
cat(obs.silhouette)
# So this says we should have 2 clusters and the silhouette coefficient (0.142112) is similar to
# the enterotype paper (0.1899451), though ideally it should be close to 1.0
#   3 clusters silhouette = 0.1089003
#   4 clusters silhouette = 0.09387917


# Remove those genera for which the average abundance across all samples is below 0.01%
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  # print(percent)
  return(Matrix_1)
}

tax.to.clust <- noise.removal(tax.to.clust, percent=0.01)


if (nclusters == 2) { 
  clust.col <- c("indianred1", "turquoise3")
} else if (nclusters == 3) { 
  clust.col <- c("indianred1", "green4", "dodgerblue")
} else if (nclusters == 4) { 
  clust.col <- c("indianred1", "yellowgreen", "turquoise3", "purple1")
} else if (nclusters == 5) {
  clust.col <- c("indianred1", "#829311", "seagreen3", "dodgerblue", "orchid")
} else if (nclusters == 6) {
  clust.col <- c("indianred1", "gold3", "green4", "turquoise3", "dodgerblue", "magenta")
} else if (nclusters == 10) {
  clust.col <- c("black", "red", "cyan", "darkgreen", "yellow", "purple", "blue", "green", "orange", "pink", "white")
}




# ************************* #
## plot 1: between-class analysis (BCA) ####
library(ade4)
library(Rmisc)

obs.pca=dudi.pca(data.frame(t(tax.to.clust)), scannf=F, nf=k)

# # ************************* #
# # check eigenvalues / variance explained per PC
# library(factoextra)
# eig.val <- get_eigenvalue(obs.pca)
# 
# fviz_eig(obs.pca, addlabels = TRUE, ncp = 15)
# 
# # results for active variables (genera)
# pca.var <- get_pca_var(obs.pca)
# 
# 
# 
# # to determine quality of representation of variables in each PC
# library(corrplot)
# corrplot(t(pca.var$cos2), is.corr = F)
# fviz_cos2(obs.pca, choice = "var", axes = 1:3)
# #this shows that, along the first 3 PCs, Prevotella, Ralstonia, Treponema have the highest quality of representation
# 
# fcon <- fviz_contrib(obs.pca, choice = "var", axes = 1:3)
# # same 3 genera are top contributors, red line is expected average contribution (1/62 => 1.6%)
# fcon.signif <- rownames(fcon$data[fcon$data$contrib >= 1.5*mean(fcon$data$contrib), ])
# # fcon.signif <- rownames(fcon$data)
# fcon$data[order(fcon$data$contrib, decreasing = T),][1:20,]
# fcon$data[order(fcon$data$contrib, decreasing = T),][20:40,]
# 
# 
# # plot variables in correlation circle
# fviz_pca_var(obs.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = T)
# 
# # can color variables by kmeans grouping
# res.km <- kmeans(pca.var$coord, centers = 3, nstart = 25)
# grp <- as.factor(res.km$cluster)
# fviz_pca_var(obs.pca, col.var = grp, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = T)
# 
# 
# 
# # results for individuals
# pca.ind <- get_pca_ind(obs.pca)
# fviz_pca_ind(obs.pca, col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              pointsize = "contrib", label="none")
# fviz_pca_ind(obs.pca, col.ind = as.factor(otu.cluster), palette = c("indianred1", "turquoise3"),
#              label="none", legend.title="Stomatotype")
# 
# 
# 
# # biplot
# clcol <- ifelse(sample_data(SLL1)[,"Stomatotype"]==1, "indianred1", "turquoise3")
# # phyla <- as.character(tax_table(SLL1)[ fcon.signif, "Phylum" ])
# phyla <- as.character(tax_table(SLL1)[ rownames(fcon$data), "Phylum" ])
# plist <- list(Actinobacteria="darkorange", Bacteroidetes="violet", Firmicutes="darkgreen", 
#               Fusobacteria="red", Proteobacteria="black", TM="purple", 
#               Spirochaetes="grey", SR="indianred", unclassified="brown", Tenericutes="darkcyan")
# phyla.col <- as.character(plist[phyla])#as.numeric(as.factor(phyla))
# 
# surv <- as.data.frame(as.matrix(sample_data(SLL1)))
# sigdifgen <- c("Haemophilus","Neisseria","Porphyromonas","Gemella","Granulicatella",
#                "unclassified.61432","unclassified.1524964","Aggregatibacter",
#                "Streptococcus","Prevotella","Veillonella","Actinomyces","Leptotrichia",
#                "TM_genus_incertae_sedis","Campylobacter","unclassified.186928","Atopobium",
#                "Megasphaera","Oribacterium","Selenomonas","Solobacterium","Eubacterium",
#                "Tannerella","unclassified.42214","Catonella","Olsenella","Butyrivibrio",
#                "unclassified.1394757","Mogibacterium","Slackia","Peptostreptococcus")
# ns1 <- nrow(sample_data(SLL)[ sample_data(SLL)[,"Stomatotype"]==1, ])
# ns2 <- nrow(sample_data(SLL)[ sample_data(SLL)[,"Stomatotype"]==2, ])
# Stom_w_labels <- as.factor(ifelse(surv$Stomatotype==1, 
#                                   sprintf("Stom 1 (n = %s)", ns1),
#                                   sprintf("Stom 2 (n = %s)", ns2)))
# 
# 
# fviz_pca_biplot(obs.pca, #axes = c(9,10),
#                 geom.ind = "point", #fill.ind = Stom_w_labels,#col.ind = surv$Stomatotype, #col.ind = "black",
#                 pointsize=2, habillage = Stom_w_labels,#surv$Stomatotype,# col.ind = "black",
#                 # palette = c("indianred1", "turquoise3", 1,2,3,4,5,6),
#                 label = "var", alpha.var=0.5, #arrowsize=factor(phyla),
#                 col.var = "black",#as.factor(phyla),
#                 # legend.title = list(fill="Stomatotype", color="Phylum"),
#                 legend.title = "Stomatotype",
#                 pointshape=16,# legend.labels = c("Stom 1 n=831", "Stom 2 n=488"),
#                 repel = T, labelsize=5,
#                 # select.var=list(name = fcon.signif),
#                 # select.var=list(name = c(fcon.signif, "Streptococcus","Haemophilus","Neisseria")),
#                 select.var=list(name = sigdifgen),
#                 ggtheme = theme_classic()) +
#   # ggpubr::fill_palette("ngp")+
#   # ggpubr::color_palette("ngp") +
#   # ggpubr::fill_palette(c("indianred1", "turquoise3")) +
#   # ggpubr::color_palette(phyla.col) +
#   xlim(-13,13) + ylim(-13,13)
# 
# set_palette(fp, "jco")
# 
# 
# #alt version of biplot:
# sdg_contribs <- sigdifgen[ sigdifgen %in% rownames(pca.var$contrib)]
# # sdg_contribs <- fcon.signif
# phyla <- as.character(tax_table(SLL1)[ sdg_contribs, "Phylum" ])
# phyla.col <- as.character(plist[phyla])
# 
# fpi <- fviz_pca_ind(obs.pca, col.ind = as.factor(otu.cluster), 
#                     palette = c("indianred1", "turquoise3"), pointsize=2,
#                     label="none", legend.title="Stomatotype") + 
#   xlim(-13,13) + ylim(-13,13)
# 
# fviz_add(fpi, pca.var$contrib[sdg_contribs, ], color=phyla.col, geom = "arrow",
#          repel = T, labelsize = 6, linetype = "solid") #+
# #  legend(12, 0, factor(phyla))
# 
#  # ************************* #


obs.bet=bca(obs.pca, fac=as.factor(otu.cluster), scannf=F, nf=k-1)
if (nclusters > 2) {
  #with only 2 clusters, apparently unable to produce figure bc does not create enough columns in obs.bet$ls
  s.class(obs.bet$ls, fac=as.factor(otu.cluster), grid=F,sub="Between-class analysis", col=clust.col)
}

# determine main contributor of a given cluster (heretofore referred to as ****stomatotypes****)
stomato.leaders <- NULL
for (i in 1:nclusters){
  stomato.leaders[i] <- colnames(obs.bet$tab)[obs.bet$tab[i,]==max(obs.bet$tab[i,])]
}

# check boxplots of rel abund of each stomatotype
samples.by.cluster <- as.data.frame(otu.cluster)
samples.by.cluster[,2] <- names(otu.cluster) #adding a column for corresponding stomatotype with a given sample

#add stomatotype to sample_data
if (gender == "both") {
  sample_data(SLL)[,'Stomatotype'] <- samples.by.cluster[,1]
  sample_data(SLL1)[,'Stomatotype'] <- samples.by.cluster[,1]
} else if (gender == "CORE") {
  sample_data(SLL)[,'Stomatotype_CORE'] <- samples.by.cluster[,1]
  sample_data(SLL1)[,'Stomatotype_CORE'] <- samples.by.cluster[,1]
}



stomato1.m <- cbind( melt(t(tax.to.clust[stomato.leaders[1], ])), melt(samples.by.cluster) )
stomato1.m <- stomato1.m[,c(1,2,3,6)]
colnames(stomato1.m)[4] <- 'cluster'
stomato2.m <- cbind( melt(t(tax.to.clust[stomato.leaders[2], ])), melt(samples.by.cluster) )
stomato2.m <- stomato2.m[,c(1,2,3,6)]
colnames(stomato2.m)[4] <- 'cluster'
if (nclusters > 2) {
  stomato3.m <- cbind( melt(t(tax.to.clust[stomato.leaders[3], ])), melt(samples.by.cluster) )
  stomato3.m <- stomato3.m[,c(1,2,3,6)]
  colnames(stomato3.m)[4] <- 'cluster'
} 
if (nclusters > 3) {
  stomato4.m <- cbind( melt(t(tax.to.clust[stomato.leaders[4], ])), melt(samples.by.cluster) )
  stomato4.m <- stomato4.m[,c(1,2,3,6)]
  colnames(stomato4.m)[4] <- 'cluster'
} 
if (nclusters > 4) {
  stomato5.m <- cbind( melt(t(tax.to.clust[stomato.leaders[5], ])), melt(samples.by.cluster) )
  stomato5.m <- stomato5.m[,c(1,2,3,6)]
  colnames(stomato5.m)[4] <- 'cluster'
}

# to give the same scale in each boxplot
ymax <- max(tax.to.clust[stomato.leaders,])

p1 <- ggplot(stomato1.m, aes(x=factor(cluster, levels=c(1:nclusters)), y=value, fill=factor(cluster, levels=c(1:nclusters)))) +
  geom_boxplot() +
  ggtitle(stomato.leaders[1]) + ylim(0, ymax) +
  xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
p2 <- ggplot(stomato2.m, aes(x=factor(cluster, levels=c(1:nclusters)), y=value, fill=factor(cluster, levels=c(1:nclusters)))) +
  geom_boxplot() +
  ggtitle(stomato.leaders[2]) + ylim(0, ymax) +
  xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
if (nclusters > 2) {
  p3 <- ggplot(stomato3.m, aes(x=factor(cluster, levels=c(1:nclusters)), y=value, fill=factor(cluster, levels=c(1:nclusters)))) +
    geom_boxplot() +
    ggtitle(stomato.leaders[3]) + ylim(0, ymax) +
    xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
} 
if (nclusters > 3) {
  p4 <- ggplot(stomato4.m, aes(x=factor(cluster, levels=c(1:nclusters)), y=value, fill=factor(cluster, levels=c(1:nclusters)))) +
    geom_boxplot() +
    ggtitle(stomato.leaders[4]) + ylim(0, ymax) +
    xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
} 
if (nclusters > 4) {
  p5 <- ggplot(stomato5.m, aes(x=factor(cluster, levels=c(1:nclusters)), y=value, fill=factor(cluster, levels=c(1:nclusters)))) +
    geom_boxplot() +
    ggtitle(stomato.leaders[5]) + ylim(0, ymax) +
    xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
}

if (nclusters == 2) {
  multiplot(p1, p2, cols=2)
} else if (nclusters == 3) {
  multiplot(p1, p2, p3, cols=3)
} else if (nclusters == 4) {
  multiplot(p1, p2, p3, p4, cols=4)
} else if (nclusters == 5) {
  multiplot(p1, p2, p3, p4, p5, cols=5)
}


# ************************* #
# plot 2: principal coordinates analysis (PCoA) ####
# *** In the tutorial by Bork created after their paper was published, they said that 
#     PCoA is actually better than BCA for these analyses....
#     ...But HOW to determine main contributors of a given cluster as in the BCA above ?????
obs.pcoa=dudi.pco(otu.jsd, scannf=F, nf=3)
if (length(table(sample_data(SLL)[,"Stomatotype"])) == 2) {
  ns1 <- nrow(sample_data(SLL)[ sample_data(SLL)[,"Stomatotype"]==1, ])
  ns2 <- nrow(sample_data(SLL)[ sample_data(SLL)[,"Stomatotype"]==2, ])
  labs <- c(sprintf("Stom %s \n n=%s", levels(as.factor(otu.cluster))[1],ns1), 
            sprintf("Stom %s \n n=%s", levels(as.factor(otu.cluster))[2],ns2))
} else {
  labs <- levels(as.factor(otu.cluster))
}

s.class(obs.pcoa$li, fac=as.factor(otu.cluster), grid=F, col=clust.col, clabel=1.5,
        label=labs)#, sub="Principal coordiante analysis")
# text(x=-0.5, y=0.5, labels=sprintf("n=%s", ns1), cex=1.7, font=2)
# text(x=1, y=0.5, labels=sprintf("n=%s", ns2), cex=1.7, font=2)
# obs.bet2=bca(obs.pcoa, fac=as.factor(otu.cluster), scannf=F, nf=k-1) 
# s.class(obs.bet2$ls, fac=as.factor(otu.cluster), grid=F,sub="Between-class analysis", col=c(3,2,4,6,5))



# ************************* #
# plot 3: 3d PCA plots ####
library(pca3d)
library(scales)
otu.prcomp <- prcomp(otu.jsd)
pca3d(otu.prcomp, group=as.factor(otu.cluster), show.centroids = T, show.group.labels = T, palette = hue_pal()(nclusters), 
      show.ellipses = T, ellipse.ci=0.75, show.plane = F)
# snapshotPCA3d(file=sprintf("%s/Part_1/figures/PCA/3D_PCA_%s_clust.png", home_dir, nclusters))

# ************************* #





# ****************************************************************************************************************** #
# ************************************************ #
# Boxplots of top15 within each stomatotype ####
# ************************************************ #

top15 <- names(sort(taxa_sums(ttc), decreasing = T))[1:15]
SLL.top <- cbind( melt(t(as.data.frame(tax.to.clust)[top15,])), melt(samples.by.cluster))


SLL.top <- SLL.top[, c(1,2,3,6)]
colnames(SLL.top) <- c("sample","OTU","value","cluster")
# to order the boxes by overall abundance
SLL.top$OTU <- as.character(SLL.top$OTU)
SLL.top$OTU <- factor(SLL.top$OTU, levels = top15)

cluster_labels <- c('1'=sprintf('Stomatotype 1 (n=%s)',sum(samples.by.cluster$otu.cluster==1)),
                    '2'=sprintf('Stomatotype 2 (n=%s)',sum(samples.by.cluster$otu.cluster==2)),
                    '3'=sprintf('Stomatotype 3 (n=%s)',sum(samples.by.cluster$otu.cluster==3)),
                    '4'=sprintf('Stomatotype 4 (n=%s)',sum(samples.by.cluster$otu.cluster==4)),
                    '5'=sprintf('Stomatotype 5 (n=%s)',sum(samples.by.cluster$otu.cluster==5)),
                    '6'=sprintf('Stomatotype 6 (n=%s)',sum(samples.by.cluster$otu.cluster==6)),
                    '7'=sprintf('Stomatotype 7 (n=%s)',sum(samples.by.cluster$otu.cluster==7)),
                    '8'=sprintf('Stomatotype 8 (n=%s)',sum(samples.by.cluster$otu.cluster==8)),
                    '9'=sprintf('Stomatotype 9 (n=%s)',sum(samples.by.cluster$otu.cluster==9)),
                    '10'=sprintf('Stomatotype 10 (n=%s)',sum(samples.by.cluster$otu.cluster==10)))

# ggplot(SLL.top, aes(x=reorder(factor(Var2, levels=top15), -value, FUN=median), y=value, 
# fill=reorder(factor(Var2, levels=top15), -value, FUN=median))) +
ggplot(SLL.top, aes(x=OTU, y=value, # This way keeps the order of top15 instead of ordering by median
                   fill=OTU)) +
  geom_boxplot() +
  facet_wrap(~cluster, ncol=1, labeller=as_labeller(cluster_labels)) +
  # theme(axis.text.x=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust=0.5)) +
  ggtitle(sprintf('Percent of given %s per sample by indicated stomatotype', "Genus")) +
  xlab("Genus") + ylab('% per sample') + scale_fill_hue(name="Genus") + guides(fill=F)




# ********** #
# to do paired boxes for top5
top5 <- names(sort(taxa_sums(ttc), decreasing = T))[1:5]
SLL.top5 <- cbind( melt(t(as.data.frame(tax.to.clust)[top5,])), melt(samples.by.cluster))


SLL.top5 <- SLL.top5[, c(1,2,3,6)]
colnames(SLL.top5) <- c("sample","Genus","value","Stomatotype")
# to order the boxes by overall abundance
SLL.top5$Genus <- as.character(SLL.top5$Genus)
SLL.top5$Genus <- factor(SLL.top5$Genus, levels = top5)
SLL.top5$Stomatotype <- factor(SLL.top5$Stomatotype)

ggplot(SLL.top5, aes(x=Genus, y=value, fill=Stomatotype)) +
  geom_boxplot() + theme_minimal() +
  theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=14),
        axis.title.y=element_blank(), axis.title.x=element_text(size=14, face="bold")) +
  guides(fill = guide_legend(title.theme=element_text(size=13, angle=0, face="bold"), 
                             label.theme=element_text(size=16, angle=0)) ) + 
  stat_compare_means(method = "t.test", label = "p.format", cex=5, color="red3", 
                     label.y=c(61, 36, 42, 39, 37))
  # theme(axis.text = element_text(size = 14), plot.title = element_text(hjust=0.5)) +
  # ggtitle(sprintf('Percent of given Genus per sample by indicated Stomatotype'))
# ********** #




# t test to see if there is a difference in abundances between stomatotypes 1 and 2
top.tax <- names(sort(taxa_sums(ttc), decreasing = T))
for (g in top.tax) {
  t.1 <- t.test(as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster==1,])])), 
                as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster %in% 2:10,])])),
                alternative = "greater")$p.value
  t.2 <- t.test(as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster %in% c(1,3:10),])])),
                as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster==2,])])),
                alternative = "less")$p.value
  if (nclusters >= 3) {
    t.3 <- t.test(as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster %in% c(1,2,4:10),])])), 
                  as.numeric(as.matrix(otu_table(ttc)[ g, rownames(samples.by.cluster[samples.by.cluster$otu.cluster==3,])])),
                  alternative = "less")$p.value
  } else{ t.3 <- 1 }
  
  if (is.nan(t.1)==F & t.1*length(top.tax) < 0.05) { # pval * length(top.tax)   as bonferroni correction
    print(c(g, t.1, "Stomatotype 1"))
  } else if (is.nan(t.2)==F & t.2*length(top.tax) < 0.05) {
    print(c(g, t.2, "Stomatotype 2"))
  } else if (nclusters>=3 & is.nan(t.3)==F & t.3*length(top.tax) < 0.05) {
    print(c(g, t.3, "Stomatotype 3"))
    # } else {
    #   print(g)
  }
}







# ************************************************ #
# Co-occurence network of OTUs ####
# tutorial here: https://statnet.org/trac/raw-attachment/wiki/Resources/introToSNAinR_sunbelt_2012_tutorial.pdf
# ************************************************ #

clus_n <- 0
if (clus_n==0) {
  clus = tax.to.clust
} else {
  clus <- tax.to.clust[ , rownames(sample_data(SLL1)[ sample_data(SLL1)[,"Stomatotype"] == clus_n, ]) ]
}
# must make the otu_table a as.numeric so it can be read by cor.test()
otu.subset <- apply(clus, 2, as.numeric) # try with c2 after, eventually make if elses for potential additional clusters, as above
rownames(otu.subset) <- rownames(clus)

res.matrix <- matrix(NA, nrow=length(rownames(clus)), ncol=length(rownames(clus)))
colnames(res.matrix) <- rownames(clus)
rownames(res.matrix) <- rownames(clus)
ps.matrix <- matrix(NA, nrow=length(rownames(clus)), ncol=length(rownames(clus)))
colnames(ps.matrix) <- rownames(clus)
rownames(ps.matrix) <- rownames(clus)

for (i in rownames(clus)) {
  for (j in rownames(clus)) {
    correl <- cor.test(otu.subset[i,], otu.subset[j,], na.rm=T)
    res.matrix[i,j] <- correl$estimate
    ps.matrix[i,j] <- correl$p.value
  }
}

ps.matrix.adj <- apply(ps.matrix, 2, p.adjust, method='bonferroni', n=length(rownames(clus)))


# p < 0.05 and cor_coeff > 0.6 for the graph in sna below****
# from soil bacteria paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4649028/

library(sna)
library(network)

# First must prepare a matrix that will serve as the coordinates for vertices in a graph of the network
# Square matrix where values are 1 if meets requirement of cor_coeff and p-val, 0 if not

top20 <- names(sort(taxa_sums(ttc), decreasing = T))[1:20] 
sna.gr <- matrix(0, nrow=length(top20), ncol=length(top20))
rownames(sna.gr) <- top20
colnames(sna.gr) <- top20
# Make another matrix 

min.cor <- 0.25
for (i in rownames(sna.gr)) {
  for (j in colnames(sna.gr)) {
    if (abs(res.matrix[i,j]) > min.cor & ps.matrix.adj[i,j] < 0.05) {
      if (i != j) {
        sna.gr[i,j] <- 1
      }
    }
  }
}

# Keep only those OTUs that have at least 1 edge (high enough cor coeff and signif)
sna.gr <- sna.gr[rowSums(sna.gr)>=1,colSums(sna.gr)>=1]
sna.net <- as.network(sna.gr, directed=FALSE)

deg <- degree(sna.net, gmode="graph") # Indegree for MIDs
phyla <- as.character(tax_table(SLL1)[ colnames(as.sociomatrix(sna.net)), "Phylum" ])
plist <- list(Actinobacteria="purple", Bacteroidetes="green", Firmicutes="yellow", Fusobacteria="cyan",
              Proteobacteria="darkorange", TM="black", Spirochaetes="grey", SR="indianred")
phyla.col <- as.character(plist[phyla])#as.numeric(as.factor(phyla))

# use as.edgelist(sna.net) to plot edge width and color based on correlation coefficient
edges <- as.edgelist(sna.net)
e.cor <- vector(length=nrow(edges))
for (i in 1:nrow(edges)) {
  snd <- attr(edges, "vnames")[ edges[i, ][1] ]
  rec <- attr(edges, "vnames")[ edges[i, ][2] ]
  e.cor[i] <- res.matrix[snd,rec]
}
e.colors <- ifelse(e.cor>0,"red","blue")
e.sizes <- (abs(e.cor*(1/min.cor)))^3


coords <- gplot(sna.net, gmode="graph", vertex.col=phyla.col, #vertex.cex=(deg)^1.5/5,
                displaylabels=TRUE, label.cex=1.1, label.pos=3, label.border=F, boxed.labels=T,
                edge.lwd=e.sizes, edge.col=e.colors) # mult by (1/min.cor) to ensure each is at least 1 before squaring 

if (clus_n==1 | clus_n==2) {
  # title(sprintf("Co-occurrence network of genera within stomatotype %s", clus_n))
  title(sprintf("Stomatotype %s", clus_n)) # for paper
  if (clus_n == 1) { # this is so as to combine the 2 stomatotypes for the figure in the paper
    legend(x=max(coords[,"x"])-2,y=max(coords[,"y"])+1.25, unique(phyla), 
           fill=unique(phyla.col), bty="n", title="Phylum", cex=1.5)
  } else if (clus_n == 2) {
    legend(x=min(coords[,"x"])-0.75,y=max(coords[,"y"])+1.5, c("Signif (+) cor", "Signif (-) cor"), 
           fill=c("red","blue"), bty="n", title="Edge color", cex=1.5)
  }
  
} else {
  # title(sprintf("Co-occurrence network of genera within all samples")) # dont use for paper figure
  # legend(x=min(coords[,"x"]),y=min(coords[,"y"])+2, unique(phyla), 
  #        fill=unique(phyla.col), bty="n", title="Phylum")
  legend(x=max(coords[,"x"])-2.5,y=max(coords[,"y"])+1, unique(phyla), 
         fill=unique(phyla.col), bty="n", title="Phylum", cex=1.5)
  legend(x=min(coords[,"x"]),y=max(coords[,"y"])+1, c("Signif (+) cor", "Signif (-) cor"), 
         fill=c("red","blue"), bty="n", title="Edge color", cex=1.5)
}








# ****************************************************************************************************************** #
# ********************************************************** #
# Heatmap of subjects, by area/day colored by stomatotype ####
# ********************************************************** #

SD <- sample_data(ttc)

school <- sort(as.numeric(as.matrix(unique(SD[,"Qsamples.1"]))))
l <- as.character(min(as.numeric(school)):max(as.numeric(school)))
indiv <- mapply(function(i) nrow(SD[SD[,"Qsamples.1"]==i]), l)

ph   <- sort(as.matrix(unique(SD[,"Q00.PH"])))

# subjs <- sort(as.matrix(unique(SD[,"subject"])))
# areas <- sort(as.matrix(unique(SD[,"area"]))) #must sort bc first appearance of area in data table is A4
# days  <- sort(as.matrix(unique(SD[,"day"])))
# samps <- sort(as.matrix(unique(SD[,"samplingMethod"])))
# paste_s_a <- function(x) { paste(x, areas, sep="_") }
# paste_d_s <- function(x) { paste(x, samps, sep="_") }
# 
# dermatomap <- matrix(NA, nrow = length(subjs)*length(areas), ncol = length(days)*length(samps))
# rownames(dermatomap) <- as.character(sapply(subjs, paste_s_a))
# colnames(dermatomap) <- as.character(sapply(days,  paste_d_s))

stomatomap <- matrix(NA, nrow=length(school), ncol=max(indiv))
rownames(stomatomap) <- as.character(school)
colnames(stomatomap) <- 1:max(indiv)

dwc <- cbind(SD, as.data.frame(otu.cluster))
# dwc <- transform(dwc, AgeCuts=cut(dwc$Age, 4))
# dwc <- transform(dwc, Dilution=gsub("/", ".", dwc$Dilution))
# dwc <- transform(dwc, Dilution=ifelse(dwc$Dilution=="0", "Zero", dwc$Dilution))


for (i in rownames(stomatomap)) {
  for (j in colnames(stomatomap)) {
    if (as.numeric(j) > indiv[i]) {
      stomatomap[i,j] <- 0
    } else {
      stomatomap[i,j] <- dwc[dwc$Qsamples.1==i, "otu.cluster"][as.numeric(j)]
    }
  }
}
rownames(stomatomap) <- paste('s', rownames(stomatomap), sep='_')
colnames(stomatomap) <- paste('n', colnames(stomatomap), sep='_')



row.ord <- hclust( dist( stomatomap, method="euclidean" ), method = "centroid" )$order
col.ord <- hclust( dist( t(stomatomap), method="euclidean" ), method = "centroid" )$order
# stom.m <- melt(stomatomap[ row.ord, rev(col.ord)])
stom.m <- melt(stomatomap)
stom.m[stom.m[,"value"] == 0, "value"] <- "No sample"



if (nclusters == 2) { 
  co <- colorRampPalette(c("indianred1", "turquoise3", "white"))(3)
  frameco <- "red"
} else if (nclusters == 3) { 
  co <- colorRampPalette(c("indianred1", "green4", "dodgerblue", "white"))(4)
  frameco <- "black"
} else if (nclusters == 4) { 
  # co <- colorRampPalette(c("lightblue", "red", "green", "yellow", "white"))(5)
  co <- colorRampPalette(c("indianred1", "yellowgreen", "turquoise3", "purple1", "white"))(5)
  frameco <- "black"
} else if (nclusters == 5) {
  co <- colorRampPalette(c("indianred1", "gold3", "green3", "dodgerblue", "orchid", "white"))(6)
  frameco <- "black"
} else if (nclusters == 6) {
  co <- colorRampPalette(c("indianred1", "gold3", "green4", "turquoise3", "dodgerblue", "magenta", "white"))(7)
  frameco <- "black"
} else if (nclusters == 10) {
  co <- colorRampPalette(c("black", "red", "cyan", "darkgreen", "yellow", "purple", "blue", "green", "orange", "pink", "white"))(11)
  frameco <- "black"
}

######## Horizontal heatmap without framing cells based on particular measures:
ggplot(data = stom.m, aes(x=Var1, y=Var2, fill=as.character(stom.m$value))) + 
  geom_tile(color="black") + 
  geom_tile(color="white", show.legend=F) + #maintains borders only in the legend
  scale_fill_manual("Stomatotype", values=setNames(co, levels(stom.m$value)), guide="legend") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust=0.5),
        axis.text.y = element_text(angle = 0, size = 10, hjust=1, vjust=0.5)) +
  guides(color = guide_legend()) +
  ggtitle(sprintf('Stomatotypes')) +
  xlab('School ID') + ylab('Sample # within given school')














# # ****************************************************************************************************************** #
# #Create heatmap of correlations with ICC
# # ****************************************************************************************************************** #
# 
# toremove <- c(1,2,13:34,41,43,45:50,99,101)
# #excel col:  (A,B, M-AH,AO,AQ,AS-AX,CU,CW )
# data.mix <- cbind(t(otu_table(SLL1)), sample_data(SLL1)[,-toremove])
# otus.mix <- rownames(otu_table(SLL1))
# questions.mix <- colnames(sample_data(SLL1))[-toremove]
# 
# cor.matrix <- matrix(NA, nrow=length(otus.mix), ncol=length(questions.mix))
# colnames(cor.matrix) <- questions.mix
# rownames(cor.matrix) <- otus.mix
# 
# for (i in otus.mix) {
#   for (j in questions.mix) {
#     cor.matrix[i,j] <- ICCbare(j, i, data.mix)
#   }
# }
# 
# 
# #good max/min cors within samples
# maxes <- apply(cor.matrix,2,max)
# goodmax <- maxes[maxes > 0.03]
# mins <- apply(cor.matrix,2,min)
# goodmin <- mins[mins < -0.03]
# goodcor <- unique(c(names(goodmax), names(goodmin)))
# 
# #good max/min cors within taxa
# tmaxes <- apply(cor.matrix,1,max)
# tgoodmax <- tmaxes[tmaxes > 0.03]
# tmins <- apply(cor.matrix,1,min)
# tgoodmin <- tmins[tmins < -0.03]
# tgoodcor <- unique(c(names(tgoodmax), names(tgoodmin)))
# good.cor.matrix <- cor.matrix[tgoodcor,goodcor]
# 
# 
# 
# 
# # cor.m <- melt(cor.matrix)
# # orderednames <- names(sort(rowSums(otu_table(SLL1))))
# # cor.m <- melt(good.cor.matrix[orderednames, ]) #to order names in figure by abundance
# 
# #cluster rows and cols for heatmap
# data <- scale(good.cor.matrix)
# ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
# data2 <- scale(t(good.cor.matrix))
# ord2 <- hclust( dist(data2, method = "euclidean"), method = "ward.D" )$order
# cor.m <- melt(good.cor.matrix[ord,ord2])
# 
# ## heatmap with ggplot
# ggplot(data = cor.m, aes(x=Var2, y=Var1, fill=value)) + 
#   geom_tile(color="white") + 
#   scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="ICC") +
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1)) +
#   ggtitle(sprintf('ICC correlations of survey with taxa at level %d', level)) +
#   xlab('') + ylab('')
# #   scale_x_discrete(limits=colnames(good.cor.matrix)) +
# #   scale_y_discrete(limits=rownames(good.cor.matrix))
# 
# #strongest correlation here...birth city of father (Q7) with Cardiobacterium (tax.2717)
# #also notable:
# #  Spiro...Treponema (tax.157) with leche (Q22.1), cafe (Q28.1), alcohol (Q28.9), Socioeconomic
# #  Fusobacterium (tax.848) with leche (Q22.1), cafe (Q28.1), alcohol (Q28.9), piercing (Q31), Socioeconomic
# #Socioeconomic seems to have the most high correlations...high code meant poorer school
# #  Actino...unclassified (tax.38987) with city/country questions
# #  As well as Firmic..Clostrid..Lachno..unclass (tax.186928), Fuso..Strep (tax.34104), many others
# 
# 
# #  cafe (28.1), alcohol (28.9 have many high correlations)
# 
# ## normal heatmap
# # heatmap(cor.matrix)
# # heatmap(good.cor.matrix)








# ****************************************************************************************************************** #
# Diversity within (alpha) and between (beta) samples ####
# ****************************************************************************************************************** #

tlev <- "Genus"; otutab <- otu_table(SLL); otutab_rel <- otu_table(SLL1)
# tlev <- "Phylum"; otutab <- tlev_otus[[tlev]]; otutab_rel <- tlev_otus_rel[[tlev]]

if (region_type == "all") {
  phy_obj <- SLL
  phy_obj_rel <- SLL1
} else if (region_type == "region") {
  phy_obj <- SLL.region
  phy_obj_rel <- SLL1.region
  otutab <- otutab[ , region_samples ]
  otutab_rel <- otutab_rel[ , region_samples ]
}


# ********************************************* #
# Will divide samples into quartiles based on their diversity values
# We have a number of different diversity measures, so will take a look at each here:

# div.estimate <- "Div.Observed"
# div.estimate <- "Div.Chao1"
# div.estimate <- "Div.ACE"
div.estimate <- "Div.Shannon"
# div.estimate <- "Div.Simpson"
# div.estimate <- "Div.InvSimpson"
# div.estimate <- "Div.Fisher"
# div.estimate <- "Weighted_Unifrac"
# div.estimate <- "Unweighted_Unifrac"
# div.estimate <- "Faiths.PD"
# div.estimate <- "Species_Richness"

divs <- as.numeric(sample_data(phy_obj_rel)[,div.estimate][[1]])
quarts <- quantile(divs, probs=seq(0,1,0.25))

# get samples by diversity quartile
q1.samples <- rownames(sample_data(phy_obj)[ sample_data(phy_obj)[,div.estimate] < quarts['25%'], ] )
q2.samples <- rownames(sample_data(phy_obj)[ (quarts['25%'] <= sample_data(phy_obj)[,div.estimate]) & (sample_data(phy_obj)[,div.estimate] < quarts['50%']), ] )
q3.samples <- rownames(sample_data(phy_obj)[ (quarts['50%'] <= sample_data(phy_obj)[,div.estimate]) & (sample_data(phy_obj)[,div.estimate] < quarts['75%']), ] )
q4.samples <- rownames(sample_data(phy_obj)[ sample_data(phy_obj)[,div.estimate] >= quarts['75%'], ] )

# add diversity group to sample_data 
sample_data(SLL)[q1.samples, sprintf("Diversity_group_%s", div.estimate)] <- "Low"
sample_data(SLL1)[q1.samples, sprintf("Diversity_group_%s", div.estimate)] <- "Low"
sample_data(SLL)[c(q2.samples,q3.samples), sprintf("Diversity_group_%s", div.estimate)] <- "Average"
sample_data(SLL1)[c(q2.samples,q3.samples), sprintf("Diversity_group_%s", div.estimate)] <- "Average"
sample_data(SLL)[q4.samples, sprintf("Diversity_group_%s", div.estimate)] <- "High"
sample_data(SLL1)[q4.samples, sprintf("Diversity_group_%s", div.estimate)] <- "High"

# update our objects after adding the Diversity_group column
if (region_type == "all") {
  phy_obj <- SLL
  phy_obj_rel <- SLL1
} else if (region_type == "region") {
  phy_obj <- SLL.region
  phy_obj_rel <- SLL1.region
}

# Get otu_tables for samples of each quantile
q1 <- otutab[ ,rownames(sample_data(phy_obj)[ sample_data(phy_obj)[,div.estimate] < quarts['25%'], ] ) ]
q2 <- otutab[ ,rownames(sample_data(phy_obj)[ (quarts['25%'] <= sample_data(phy_obj)[,div.estimate]) & (sample_data(phy_obj)[,div.estimate] < quarts['50%']), ] ) ]
q3 <- otutab[ ,rownames(sample_data(phy_obj)[ (quarts['50%'] <= sample_data(phy_obj)[,div.estimate]) & (sample_data(phy_obj)[,div.estimate] < quarts['75%']), ] ) ]
q4 <- otutab[ ,rownames(sample_data(phy_obj)[ sample_data(phy_obj)[,div.estimate] >= quarts['75%'], ] ) ]

q1.rel <- otutab_rel[ ,rownames(sample_data(phy_obj_rel)[ sample_data(phy_obj_rel)[,div.estimate] < quarts['25%'], ] ) ]
q2.rel <- otutab_rel[ ,rownames(sample_data(phy_obj_rel)[ (quarts['25%'] <= sample_data(phy_obj_rel)[,div.estimate]) & (sample_data(phy_obj_rel)[,div.estimate] < quarts['50%']), ] ) ]
q3.rel <- otutab_rel[ ,rownames(sample_data(phy_obj_rel)[ (quarts['50%'] <= sample_data(phy_obj_rel)[,div.estimate]) & (sample_data(phy_obj_rel)[,div.estimate] < quarts['75%']), ] ) ]
q4.rel <- otutab_rel[ ,rownames(sample_data(phy_obj_rel)[ sample_data(phy_obj_rel)[,div.estimate] >= quarts['75%'], ] ) ]

# low <- otutab_rel[ , sample_data(SLL1)$Diversity_group == "Low"]
# average <- otutab_rel[ , sample_data(SLL1)$Diversity_group == "Average"]
# high <- otutab_rel[ , sample_data(SLL1)$Diversity_group == "High"]

# 
# # boxplots of genera in the different quartiles
# # first check for genera in which there are significant differences between diversity groups
# signif_ttests <- function(x , n) {
#   tp1 <- t.test( as.numeric(low[x, ]), c(average[x,], high[x,]) )$p.value * n
#   tp2 <- t.test( as.numeric(average[x, ]), as.numeric(low[x,]) )$p.value * n
#   tp3 <- t.test( as.numeric(average[x, ]), as.numeric(high[x,]) )$p.value * n
#   tp4 <- t.test( as.numeric(high[x, ]), c(average[x,], low[x,]) )$p.value * n
#   
#   # return( tp1$p.value<0.05 | tp2$p.value<0.05 | tp3$p.value<0.05 | tp4$p.value<0.05 )
#   # if (min(tp1, tp2, tp3, tp4) < 0.005) print(c(x, min(tp1, tp2, tp3, tp4)))
#   return( min(tp1, tp2, tp3, tp4) < 0.05 )
# }
# 
# signif_genera <- sapply(rownames(otutab_rel), signif_ttests, length(rownames(otutab_rel)))
# signif_genera <- names(signif_genera[signif_genera])
# signif_genera <- signif_genera[ ! startsWith(signif_genera, "unclassified") ]
# signif_genera <- sapply( signif_genera, function(x) max(otutab_rel[x,]) > 5 )
# signif_genera <- names(signif_genera[signif_genera])

# check genera with signif differenecs amongst diversity groups with kruskal wallis test
signif_kruskal <- function(x, n) {
  dg <- sprintf("Diversity_group_%s", div.estimate)
  kp <- kruskal.test(as.numeric(otutab_rel[x,]), 
                     as.factor( as.character(as.matrix(sample_data(SLL1))[ , dg]) ))$p.value * n
  # if (kp < 0.05) print(c(x, kp, kp/n))
  return(kp)
}

kruskal_signif_genera <- sapply(rownames(otutab_rel), signif_kruskal, length(rownames(otutab_rel)))
# get only those significant genera
kruskal_signif_genera <- kruskal_signif_genera[kruskal_signif_genera<0.05]
# # remove the unclassified genera
# kruskal_signif_genera <- kruskal_signif_genera[ ! startsWith(names(kruskal_signif_genera), "unclassified") ]
# get the names of those genera that have an abundance of at least 5% in at least one sample
kruskal_signif_names <- sapply( names(kruskal_signif_genera), function(x) max(otutab_rel[x,]) > 5 )
kruskal_signif_genera <- kruskal_signif_genera[kruskal_signif_names]
# add stars for level of significance
sig_stars <- sapply(kruskal_signif_genera, function(x) ifelse(x<0.0005, '***',
                                                              ifelse(x<0.005, '**',
                                                                     ifelse(x<0.05, '*', NA))))


# object for ggplot
po_by_div <- cbind(sample_data(phy_obj_rel)[,sprintf("Diversity_group_%s", div.estimate)], sample_names(phy_obj_rel), t(otutab_rel)[,names(kruskal_signif_genera)])#[,1:5])
colnames(po_by_div)[2] <- "Sample_names"
# if only one significant OTU, name will not be kept as 3rd argument of cbind above will be numeric instead of matrix
if (length(kruskal_signif_genera) == 1) colnames(po_by_div)[3] <- names(kruskal_signif_genera)
po_by_div.m <- melt(po_by_div, id=1:2, measure=3:ncol(po_by_div))
colnames(po_by_div.m)[1] <- "Div_group" # to keep a general name that can be called by ggplot. Div.estimate will still be used in title



#***make a function for the x value for geom_jitter to say -0.1 for low, +0.1 for high*****


ymaxes <- sapply(names(kruskal_signif_genera), function(x) max(as.numeric(as.matrix(otutab_rel[x,])))+3)

ggplot(po_by_div.m, aes(x=variable, y=value, color=factor(Div_group, levels=c("Low", "Average", "High")))) +
  geom_boxplot(outlier.colour = NULL, outlier.alpha = 0.4) +
  ggtitle(sprintf("%s by %s", div.estimate, tlev)) +
  xlab(tlev) + ylab("Relative abundance") + scale_color_manual(name=div.estimate, values=c('blue','goldenrod','darkgreen')) +
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust=0.5)) +
  annotate("text", x=names(kruskal_signif_genera), y=ymaxes, label=sig_stars, color="red", size=6)




## histograms of different values ####

hv <- "Div.Shannon"
# hv <- "Weighted_Unifrac"
# hv <- "Q00.PH"

if (hv=="Div.Shannon") {
  gti <- "Shannon alpha diversity"; xl <- "Shannon"; nbin <- 20;
} else if (hv=="Weighted_Unifrac") {
  gti <- "Weighted UniFrac distance (beta diversity)"; xl <- "Weighted UniFrac"; nbin <- 20;
} else if (hv=="Q00.PH") {
  gti <- "pH of donor saliva"; xl <- "Saliva pH"; nbin <- 10;
}

hdf <- as.data.frame(as.matrix(sample_data(SLL))[, hv])
colnames(hdf) <- "meas"
hdf$meas <- as.numeric(as.character(hdf$meas))

ggplot(hdf, aes(meas)) +
  geom_histogram(bins = nbin, col="black", fill="goldenrod") + 
  ggtitle(gti) + theme_minimal() + xlab(xl) +
  theme(plot.title = element_text(hjust=0.5, size=20), axis.text = element_text(size = 15),
        axis.title = element_text(size=15))
  



# **************************************************************************************************** #
# Add values for BMI categories ####
# **************************************************************************************************** #

bmi.divs <- as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))
#ignore the very obvious outliers
bmi.divs <- bmi.divs[bmi.divs < 50]
bmi.quarts <- quantile(bmi.divs, probs=seq(0,1,0.25), na.rm = T)

# get samples by BMI quartile
q1.samples <- rownames(sample_data(SLL1)[ as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < bmi.quarts['25%'], ] )
q2.samples <- rownames(sample_data(SLL1)[ (bmi.quarts['25%'] <= as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) & (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < bmi.quarts['50%']), ] )
q3.samples <- rownames(sample_data(SLL1)[ (bmi.quarts['50%'] <= as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) & (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < bmi.quarts['75%']), ] )
q4.samples <- rownames(sample_data(SLL1)[ as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) >= bmi.quarts['75%'], ] )

# add BMI group to sample_data 
sample_data(SLL)[q1.samples, "BMI_group"] <- "Low"
sample_data(SLL1)[q1.samples, "BMI_group"] <- "Low"
sample_data(SLL)[c(q2.samples,q3.samples), "BMI_group"] <- "Average"
sample_data(SLL1)[c(q2.samples,q3.samples), "BMI_group"] <- "Average"
sample_data(SLL)[q4.samples, "BMI_group"] <- "High"
sample_data(SLL1)[q4.samples, "BMI_group"] <- "High"
sample_data(SLL)[ is.na(sample_data(SLL)[, "BMI_group"]), "BMI_group"] <- "Unknown"
sample_data(SLL1)[ is.na(sample_data(SLL1)[, "BMI_group"]), "BMI_group"] <- "Unknown"


# get samples by official BMI category
#   ***Obese is not including those with BMI >= 50 since those are outliers and almost surely due to input error (only 3 samples)
underweight <- rownames(sample_data(SLL1)[ (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < 18.5) & (! is.na(as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) ), ] )
normal      <- rownames(sample_data(SLL1)[ (18.5 <= as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) & (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < 25) & (! is.na(as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) ), ] )
overweight  <- rownames(sample_data(SLL1)[ (25.0 <= as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) & (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < 30) & (! is.na(as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) ), ] )
obese       <- rownames(sample_data(SLL1)[ (30.0 <= as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) & (as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])) < 50) & (! is.na(as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) ), ] )
unknown     <- rownames(sample_data(SLL1)[  is.na(as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"]))) , ] )

# add official BMI indicators to sample_data
sample_data(SLL)[underweight, "BMI_official"] <- "Underweight"
sample_data(SLL1)[underweight, "BMI_official"] <- "Underweight"
sample_data(SLL)[c(normal), "BMI_official"] <- "Normal"
sample_data(SLL1)[c(normal), "BMI_official"] <- "Normal"
sample_data(SLL)[overweight, "BMI_official"] <- "Overweight"
sample_data(SLL1)[overweight, "BMI_official"] <- "Overweight"
sample_data(SLL)[obese, "BMI_official"] <- "Obese"
sample_data(SLL1)[obese, "BMI_official"] <- "Obese"
sample_data(SLL)[ unknown, "BMI_official"] <- "Unknown"
sample_data(SLL1)[ unknown, "BMI_official"] <- "Unknown"

# histogram of BMIs
hist( as.numeric(as.matrix(sample_data(SLL1)[,"Q3.1"])),
      main = sprintf("Histogram of BMI"),
      xlab = "BMI", col = "darkgoldenrod",
      breaks = seq(0,200,5))
#
  
  
  

# ****************************************************************************************************************** #
# Heatmap of correlations for continuous/binary responses ####
# ****************************************************************************************************************** #
#numerical columns
only_cont <- c("Q3","Q16","Q17","Q19.2","Q19.4","Q19.7","Q19.10","Q28","Q28.1","Q28.2","Q28.3","Q28.4",
               "Q28.5","Q28.6","Q28.7","Q28.8","Q28.9","Q28.10","Q29","Q29.1","Q29.2","Q29.3","Q00.PH",
               "Socioeconomic","Q3.1","Q22.2","Q23.2","Q24.2","Q25.2","Q26.2","Q39.2","Q40.2","Q28.12",
               "Q28.14","Q28.16","Q28.18","Q28.20","Q28.22","Q28.24","Q28.26","Q28.28","Q28.30","Q28.32",
               "Q16.1","Div.Observed","Div.Chao1","Div.ACE","Div.Shannon","Div.Simpson","Div.InvSimpson",
               "Div.Fisher","Faiths.PD","Species_Richness","Num_OTUs",#"Stomatotype","Stomatotype_CORE",
               "Population","Age","Weighted_Unifrac","Unweighted_Unifrac",cont_water_data)#,"Water_hardness","Neisseria_Prevotella","Neisseria_Veillonella",
               #"Haemophilus_Prevotella","Haemophilus_Veillonella","Streptococcus_Prevotella","Prevotella_Bacteroides",

cont_bin <- c("Q2","Q3","Q4","Q12","Q12.1","Q12.2","Q12.3","Q12.4","Q12.5","Q12.8","Q13","Q13.1","Q13.2","Q13.3","Q13.4","Q13.5","Q13.6","Q13.8",
              "Q16","Q17","Q18","Q19","Q19.1","Q19.2","Q19.3","Q19.4","Q19.5","Q19.7","Q19.8","Q19.10","Q20.1","Q21","Q28","Q28.1","Q28.2","Q28.3",
              "Q28.4","Q28.5","Q28.6","Q28.7","Q28.8","Q28.9","Q28.10","Q29","Q29.1","Q29.2","Q29.3","Q30","Q31","Q32","Q33","Q35","Q36","Q37",
              "Q38","Q41","Q47","Q44","Q45","Q46","Q48","Q49","Q50","Q51","Q52","Q53","Q54","Q00.PH","Socioeconomic","Qsamples.1","Q3.1","Q22.1",
              "Q22.2","Q23.1","Q23.2","Q24.1","Q24.2","Q25.1","Q25.2","Q26.1","Q26.2","Q39.1","Q39.2","Q40.1","Q40.2","Q28.11","Q28.12","Q28.13",
              "Q28.14","Q28.15","Q28.16","Q28.17","Q28.18","Q28.19","Q28.20","Q28.21","Q28.22","Q28.23","Q28.24","Q28.25","Q28.26","Q28.27","Q28.28",
              "Q28.29","Q28.30","Q28.31","Q28.32","Q16.1","Div.Observed","Div.Chao1","Div.ACE","Div.Shannon","Div.Simpson","Div.InvSimpson",
              "Div.Fisher","Faiths.PD","Species_Richness","Num_OTUs","Stomatotype","Stomatotype_CORE","Population","Age","Weighted_Unifrac",
              "Unweighted_Unifrac",cont_water_data)#,"Water_hardness","Neisseria_Prevotella","Neisseria_Veillonella","Haemophilus_Prevotella","Haemophilus_Veillonella","Streptococcus_Prevotella","Prevotella_Bacteroides"

# cont_bin <- c(3,4,5,16:21,24:31,33,36:44,46,47,49,51,52,59:77,79:82,85,88:98,100,102,106:144,149:155,  156,  158)
# #excel col:  (C,D,E, P:V,  X:AE,AG,AJ:AR,AT,AU,AW,AY,AZ,BG:BY,CA:CD,CG,CJ:CT,CV,  CX, DB:EN ,div-stats,#otus,Stomatotype)


tlev <- "Genus"; otutab_rel <- otu_table(SLL1)
# tlev <- "Order"; otutab_rel <- tlev_otus_rel[[tlev]]
# tlev <- "Phylum"; otutab_rel <- tlev_otus_rel[[tlev]]

st <- "all"
# st <- "students"
# st <- "teachers"
# st <- "no-bottled"
# st <- "only-bottled"
# st <- "Stomatotype_1"
# st <- "Stomatotype_2"

# determine which region to look at
if (region_type == "all") {
  reg_s <- sample_names(SLL1)
} else if (region_type == "region") {
  reg_s <- region_samples
}

if (st == "students") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Q1.1"]=="Student" ]
} else if (st == "teachers") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Q1.1"]=="Teacher" ]
} else if (st == "no-bottled") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Q27"]!="Embotellada" ]
} else if (st == "only-bottled") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Q27"]=="Embotellada" ]
} else if (st == "Stomatotype_1") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Stomatotype"]==1 ]
} else if (st == "Stomatotype_2") {
  reg_s <- reg_s[ sample_data(SLL1)[,"Stomatotype"]==2 ]
}

### determine which genera to include (all, those present in at least 1/3 of samples, or those not present in at least 1/3 of samples)
# gs <- transform_sample_counts(gs.filt, function(x) 100 * x/sum(x))
# gs <- SLL1
# gs <- prune_taxa(!wh.cor, SLL1)



data.cont <- cbind(t(otutab_rel[ , reg_s]), sample_data(SLL1)[reg_s, only_cont])

data.cont[data.cont=='No Sabe/No Contesta'] <- NA
data.cont <- apply(data.cont, 2, as.numeric)

otus.cont <- rownames(otutab_rel)
# questions.cont <- colnames(sample_data(gs))[only_cont]
questions.cont <- only_cont

os_qs.cont <- c(otus.cont,questions.cont)

# result <- cor(data.cont, use = 'pairwise.complete.obs')

#instead try cor.test:
res.matrix <- matrix(NA, nrow=length(os_qs.cont), ncol=length(os_qs.cont))
colnames(res.matrix) <- os_qs.cont
rownames(res.matrix) <- os_qs.cont
ps.matrix <- matrix(NA, nrow=length(os_qs.cont), ncol=length(os_qs.cont))
colnames(ps.matrix) <- os_qs.cont
rownames(ps.matrix) <- os_qs.cont


for (i in os_qs.cont) {
  for (j in os_qs.cont) {
    
    # # do pearson if gen vs gen, spearman otherwise??
    # if (i %in% questions.cont | j %in% questions.cont) {
    #   correl <- cor.test(data.cont[,j], data.cont[,i], na.rm=T, method="spearman")
    # } else {
    #   correl <- cor.test(data.cont[,j], data.cont[,i], na.rm=T, method="pearson")
    # }
    
    correl <- cor.test(data.cont[,j], data.cont[,i], na.rm=T, method="pearson")
    # correl <- cor.test(data.cont[,j], data.cont[,i], na.rm=T, method="spearman")
    res.matrix[i,j] <- correl$estimate
    ps.matrix[i,j] <- correl$p.value
  }
}

# data.cnew <- data.cont[data.cont[,"Q29"] < 10, ]
# for (i in otus.cont) { 
#   correl <- cor.test(data.cnew[,"Q29"], data.cnew[,i], na.rm=T, method="pearson") 
#   if (correl$p.value * length(otus.cont) < 0.05) {
#     print(c(i, correl$p.value, correl$estimate))
#   }
# }

ps.matrix.adj <- apply(ps.matrix, 2, p.adjust, method='bonferroni', n=length(os_qs.cont))
#*****Is this too strict since there are more comparisons if we include otus vs otus and ques vs ques?***
diag(ps.matrix.adj) <- 1
#For some questions that have all 0s (occurs when doing cities alone)
res.matrix[is.na(res.matrix)] <- 0
ps.matrix.adj[is.na(ps.matrix.adj)] <- 1

pluses <- matrix('', nrow=length(rownames(res.matrix)), ncol=length(colnames(res.matrix)))
colnames(pluses) <- colnames(res.matrix)
rownames(pluses) <- rownames(res.matrix)
pluses[ps.matrix.adj < 0.05] <- '+'


# vs <- 'all'
# vs <- 'between'
vs <- 'otus'
# vs <- 'questions'
# vs <- "water_vals"
# vs <- "gen_v_water"


# ****************************** #
if (vs=='all') { #for all vs all
  #at least 1 good min p within samples
  mins <- apply(ps.matrix.adj,2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  tmins <- apply(ps.matrix.adj,1,min)
  tgoodmin <- tmins[tmins < 0.05]
  res.all <- res.matrix[names(tgoodmin),names(goodmin)]
  plus.all <- pluses[names(tgoodmin),names(goodmin)]
  
  cor.m.cont <- cbind( melt(res.all), melt(plus.all) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- subset( cor.m.cont[ upper.tri(res.all), ] )
  if (region_type == "all") {
    title <- sprintf('All vs All cor() %s', tlev)
  } else if (region_type == "region") {
    title <- sprintf('All vs All cor() %s\n%s', tlev, reg_name)
  }
  
  # ****************************** #  
  
} else if (vs=='between') { #for questions vs otus only
  res.Que_Otu <- res.matrix[questions.cont,otus.cont]
  plus.Que_Otu <- pluses[questions.cont,otus.cont]
  
  #at least 1 good min p within samples
  mins <- apply(ps.matrix.adj[questions.cont,otus.cont],2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  tmins <- apply(ps.matrix.adj[questions.cont,otus.cont],1,min)
  tgoodmin <- tmins[tmins < 0.05]
  res.Que_Otu <- res.Que_Otu[names(tgoodmin),names(goodmin)]
  plus.Que_Otu <- plus.Que_Otu[names(tgoodmin),names(goodmin)]
  
  ord <- hclust( dist( res.Que_Otu, method="euclidean" ), method = "ward.D" )$order
  ord2 <- hclust( dist( t(res.Que_Otu), method="euclidean" ), method = "ward.D" )$order
  cor.m.cont <- cbind( melt( res.Que_Otu[ord,ord2] ), melt( plus.Que_Otu[ord,ord2] ) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- cor.m.cont ## not using lower triangle here since cols dont match rows as in the other "vs" instances
  if (region_type == "all") {
    title <- sprintf('cor() for %s', tlev)
  } else if (region_type == "region") {
    title <- sprintf('cor() for %s\n%s', tlev, reg_name)
  }
  
  # **** #
  # filter those OTUs appearing in fewer than a given number of samples
  min_samples <- 20
  good_cols <- sapply(colnames(res.Que_Otu), function(x) sum(otutab_rel[x, reg_s]!=0) >= min_samples)
  good.otus <- colnames(res.Que_Otu)[ good_cols ]
  
  res.Que_Otu.signif <- matrix(NA, nrow=nrow(res.Que_Otu), ncol=length(good.otus))
  colnames(res.Que_Otu.signif) <- good.otus
  rownames(res.Que_Otu.signif) <- rownames(res.Que_Otu)
  
  for (i in rownames(res.Que_Otu)) {
    for (j in good.otus) {
      if (pluses[i,j] == '+') {
        res.Que_Otu.signif[i,j] <- res.matrix[i,j]
      } else {
        res.Que_Otu.signif[i,j] <- ''
      }
    }
  }
  
  # remove those variables that no longer have any significant correlations after removing the rare OTUs
  good_rows <- sapply(rownames(res.Que_Otu.signif), function(x) sum(res.Que_Otu.signif[ x, ] != '') > 0)
  good_rows <- good_rows[ good_rows ]
  
  res.Que_Otu.signif <- res.Que_Otu.signif[ names(good_rows), ]
  
  write.csv(res.Que_Otu.signif, sprintf('%s/Part_1/figures/Correlations/%s/%s/PAPER_%s_%s_Que-OTUs_%s-minSamps_signif_correlations.csv', 
                                        home_dir, tlev, st, tlev, st, min_samples))
  if (tlev=="Genus" & st=="all") {
    write.csv(res.Que_Otu.signif, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_5.csv', home_dir))
  }
  # **** #
  
  # ****************************** #
  
} else if (vs=='otus') {
  res.Otu_Otu <- res.matrix[otus.cont,otus.cont]
  plus.Otu_Otu <- pluses[otus.cont,otus.cont]
  
  # **** #
  # filter those OTUs appearing in fewer than a given number of samples
  min_samples <- 20
  good_rows <- sapply(otus.cont, function(x) sum(otutab_rel[x, reg_s]!=0) >= min_samples)
  good.otus <- otus.cont[ good_rows ]
  
  res.Otu_Otu.signif <- matrix(NA, nrow=length(good.otus), ncol=length(good.otus))
  colnames(res.Otu_Otu.signif) <- good.otus
  rownames(res.Otu_Otu.signif) <- good.otus
  
  for (i in good.otus) {
    for (j in good.otus) {
      if (pluses[i,j] == '+') {
        res.Otu_Otu.signif[i,j] <- res.matrix[i,j]
      } else {
        res.Otu_Otu.signif[i,j] <- ''
      }
    }
  }
  
  write.csv(res.Otu_Otu.signif, sprintf('%s/Part_1/figures/Correlations/%s/%s/PAPER_%s_%s_OTUs-only_%s-minSamps_signif_correlations.csv', 
                                home_dir, tlev, st, tlev, st, min_samples))
  if (tlev=="Genus" & st=="all") {
    write.csv(res.Otu_Otu.signif, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_2.csv', home_dir))
  }
  # **** #
  
  # **** #
  # for heatmap figure, filter out OTUs without at least 1 count in 33% of samples
  min33 <- genefilter_sample(SLL1, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(SLL1))
  min33 <- names(min33[min33])
  
  #at least 1 good min p within samples
  # mins <- apply(ps.matrix.adj[otus.cont,otus.cont],2,min)
  mins <- apply(ps.matrix.adj[min33,min33],2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  # tmins <- apply(ps.matrix.adj[otus.cont,otus.cont],1,min)
  tmins <- apply(ps.matrix.adj[min33,min33],1,min)
  tgoodmin <- tmins[tmins < 0.05]

  
  if (st %in% c("Stomatotype_1", "Stomatotype_2")) {
    ord <- otu.row.ord
    ord2 <- otu.col.ord
    res.Otu_Otu <- res.Otu_Otu[ord,ord2]
    plus.Otu_Otu <- plus.Otu_Otu[ord,ord2]
    
  } else {
    res.Otu_Otu <- res.Otu_Otu[names(tgoodmin),names(goodmin)]
    plus.Otu_Otu <- plus.Otu_Otu[names(tgoodmin),names(goodmin)]
    ord <- hclust( dist( res.Otu_Otu, method = "euclidean" ), method = "ward.D" )$order
    ord2 <- hclust( dist( t(res.Otu_Otu), method = "euclidean" ), method = "ward.D" )$order
  }
  
  cor.m.cont <- cbind( melt( res.Otu_Otu[ord,ord2] ), melt( plus.Otu_Otu[ord,ord2] ) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- subset( cor.m.cont[ upper.tri(res.Otu_Otu[ord,ord2]), ] )
  
  ifelse(st == "Stomatotype_1", title_extra <<- sprintf(" - %s", st), 
         ifelse(st == "Stomatotype_2", title_extra <<- sprintf(" - %s", st),
                title_extra <<- "") )
  if (region_type == "all") {
    title <- sprintf('OTUs cor() %s%s', tlev, title_extra)
  } else if (region_type == "region") {
    title <- sprintf('OTUs cor() %s%s\n%s', tlev, title_extra, reg_name)
  }
  
  res.OTUs.toprint <- res.Otu_Otu
  res.OTUs.toprint[ plus.Otu_Otu != '+' ] <- ''
  write.csv(res.OTUs.toprint[rev(ord),rev(ord2)], sprintf('%s/Part_1/figures/Correlations/%s/%s/%s_%s_OTUs_third-of-samples_signif_correlations.csv', 
                                       home_dir, tlev, st, tlev, st))
  write.csv(res.Otu_Otu[rev(ord),rev(ord2)], sprintf('%s/Part_1/figures/Correlations/%s/%s/%s_%s_OTUs_third-of-samples_all-values_correlations.csv', 
                                                          home_dir, tlev, st, tlev, st))
  
  if (tlev=="Genus" & st=="all") {
    write.csv(res.OTUs.toprint[min33,min33][rev(ord),ord2], sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_2_67.csv', home_dir))
  }
  
  # ****************************** #
  
} else if (vs=='questions') {
  res.Que_Que <- res.matrix[questions.cont,questions.cont]
  plus.Que_Que <- pluses[questions.cont,questions.cont]
  #at least 1 good min p within samples
  mins <- apply(ps.matrix.adj[questions.cont,questions.cont],2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  tmins <- apply(ps.matrix.adj[questions.cont,questions.cont],1,min)
  tgoodmin <- tmins[tmins < 0.05]
  res.Que_Que <- res.Que_Que[names(tgoodmin),names(goodmin)]
  plus.Que_Que <- plus.Que_Que[names(tgoodmin),names(goodmin)]
  
  ord <- hclust( dist( res.Que_Que, method = "euclidean" ), method = "ward.D" )$order
  ord2 <- hclust( dist( t(res.Que_Que), method = "euclidean" ), method = "ward.D" )$order
  cor.m.cont <- cbind( melt( res.Que_Que[ord,ord2] ), melt( plus.Que_Que[ord,ord2] ) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- subset( cor.m.cont[ upper.tri(res.Que_Que[ord,ord2]), ] )
  if (region_type == "all") {
    title <- sprintf('Questions cor() %s', tlev)
  } else if (region_type == "region") {
    title <- sprintf('Questions cor() %s\n%s', tlev, reg_name)
  }
  
  res.Que_Que.signif <- matrix(NA, nrow=nrow(res.Que_Que), ncol=ncol(res.Que_Que))
  colnames(res.Que_Que.signif) <- colnames(res.Que_Que)
  rownames(res.Que_Que.signif) <- rownames(res.Que_Que)
  
  for (i in rownames(res.Que_Que)) {
    for (j in colnames(res.Que_Que)) {
      if (pluses[i,j] == '+') {
        res.Que_Que.signif[i,j] <- res.matrix[i,j]
      } else {
        res.Que_Que.signif[i,j] <- ''
      }
    }
  }
  
  # good_rows <- sapply(rownames(res.Que_Otu.signif), function(x) sum(res.Que_Otu.signif[ x, ] != '') > 0)
  # good_rows <- good_rows[ good_rows ]
  # 
  # res.Que_Otu.signif <- res.Que_Otu.signif[ names(good_rows), ]
  
  write.csv(res.Que_Que.signif, sprintf('%s/Part_1/figures/Correlations/%s/%s/PAPER_%s_%s_Que-Que_%s-minSamps_signif_correlations.csv', 
                                 home_dir, tlev, st, tlev, st, min_samples))
  if (tlev=="Genus" & st=="all") {
    write.csv(res.Que_Que.signif, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_6.csv', home_dir))
  }
  
  # ****************************** #
  
} else if (vs=='water_vals') {
  res.Water_vals <- res.matrix[cont_water_data,cont_water_data]
  plus.Water_vals <- pluses[cont_water_data,cont_water_data]
  #at least 1 good min p within samples
  mins <- apply(ps.matrix.adj[cont_water_data,cont_water_data],2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  tmins <- apply(ps.matrix.adj[cont_water_data,cont_water_data],1,min)
  tgoodmin <- tmins[tmins < 0.05]
  res.Water_vals <- res.Water_vals[names(tgoodmin),names(goodmin)]
  plus.Water_vals <- plus.Water_vals[names(tgoodmin),names(goodmin)]
  
  ord <- hclust( dist( res.Water_vals, method = "euclidean" ), method = "ward.D" )$order
  ord2 <- hclust( dist( t(res.Water_vals), method = "euclidean" ), method = "ward.D" )$order
  cor.m.cont <- cbind( melt( res.Water_vals[ord,ord2] ), melt( plus.Water_vals[ord,ord2] ) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- subset( cor.m.cont[ upper.tri(res.Water_vals[ord,ord2]), ] )
  if (region_type == "all") {
    title <- sprintf('Water_vals cor() %s', tlev)
  } else if (region_type == "region") {
    title <- sprintf('Water_vals cor() %s\n%s', tlev, reg_name)
  }
  
  # ****************************** #
  
} else if (vs=='gen_v_water') {
  res.Water_vals <- res.matrix[otus.cont,cont_water_data]
  plus.Water_vals <- pluses[otus.cont,cont_water_data]
  # res.Water_vals <- res.matrix[os_qs.cont,cont_water_data]
  # plus.Water_vals <- pluses[os_qs.cont,cont_water_data]
  
  #at least 1 good min p within samples
  mins <- apply(ps.matrix.adj[otus.cont,cont_water_data],2,min)
  goodmin <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
  #at least 1 good min p within taxa
  tmins <- apply(ps.matrix.adj[otus.cont,cont_water_data],1,min)
  tgoodmin <- tmins[tmins < 0.05]
  res.Water_vals <- res.Water_vals[names(tgoodmin),names(goodmin)]
  plus.Water_vals <- plus.Water_vals[names(tgoodmin),names(goodmin)]
  
  # filter those OTUs appearing in fewer than a given number of samples
  min_samples <- 20
  keep_OTUs <- sapply(rownames(res.Water_vals), function(x) sum(otutab_rel[x, reg_s]!=0) >= min_samples)
  res.Water_vals <- res.Water_vals[keep_OTUs, ]
  plus.Water_vals <- plus.Water_vals[keep_OTUs, ]
  
  ord <- hclust( dist( res.Water_vals, method = "euclidean" ), method = "ward.D" )$order
  ord2 <- hclust( dist( t(res.Water_vals), method = "euclidean" ), method = "ward.D" )$order
  cor.m.cont <- cbind( melt( res.Water_vals[ord,ord2] ), melt( plus.Water_vals[ord,ord2] ) )
  cor.m.cont <- cor.m.cont[,c(1,2,3,6)]
  colnames(cor.m.cont)[4] <- 'signif'
  cor.m.lower <- cor.m.cont ## not using lower triangle here since cols dont match rows as in the other "vs" instances
  if (region_type == "all") {
    title <- sprintf('%s: gen_v_water cor() %s', st, tlev)
  } else if (region_type == "region") {
    title <- sprintf('%s: gen_v_water cor() %s\n%s', st, tlev, reg_name)
  }
  
  res.water.toprint <- res.Water_vals
  res.water.toprint[ plus.Water_vals != '+' ] <- ''
  write.csv(res.water.toprint, sprintf('%s/Part_1/figures/Correlations/%s/%s/PAPER_%s_%s_%s-minSamps_gen-v-water_signif_correlations.csv', 
                                       home_dir, tlev, st, tlev, st, min_samples))
  if (tlev=="Genus" & st=="no-bottled") {
    write.csv(res.water.toprint, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_7.csv', home_dir))
  }
}

# ************* #
# To get heatmaps of Stomatotype 1 and 2 samples for OTUs, 
# use the following to keep row and col orders from full dataset clustering
if (st == "all" & vs == "otus") {
  otu.row.ord <- ord
  otu.col.ord <- ord2
  otu.row.ord <- rownames(res.Otu_Otu)[ord]
  otu.col.ord <- colnames(res.Otu_Otu)[ord2]
}
# ************* #

ggplot(data = cor.m.cont, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color="white", data=cor.m.lower) + 
  geom_text(aes(label=signif), size = 4.5, data=cor.m.lower) + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Pearson") +
  theme_minimal() + 
  # theme(axis.text = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
        plot.title = element_text(hjust=0.5)) +
#   scale_x_discrete(limits=rownames(cor.cont), labels=c(otus.cont)) +
#   scale_y_discrete(limits=colnames(cor.cont), labels=c(questions.cont)) +
  xlab('') + ylab('') #+ ggtitle(title)



# ******** #
# A table with the correlation coefficient only for those values that are statistically significant
#    (have a '+' in the pluses matrix)
# ******** #

# filter those OTUs appearing in fewer than a given number of samples
min_samples <- 20
good_rows <- sapply(os_qs.cont, function(x) ifelse(x %in% rownames(otutab_rel),
                                                   sum(otutab_rel[x, reg_s]!=0) >= min_samples,
                                                   TRUE))
good.os_qs.cont <- os_qs.cont[ good_rows ]

res.signif <- matrix(NA, nrow=length(good.os_qs.cont), ncol=length(good.os_qs.cont))
colnames(res.signif) <- good.os_qs.cont
rownames(res.signif) <- good.os_qs.cont

for (i in good.os_qs.cont) {
  for (j in good.os_qs.cont) {
    if (pluses[i,j] == '+') {
      res.signif[i,j] <- res.matrix[i,j]
    } else {
      res.signif[i,j] <- ''
    }
  }
}

write.csv(res.signif, sprintf('%s/Part_1/figures/Correlations/%s/%s/%s_%s_%s-minSamps_signif_correlations.csv', home_dir, tlev, st, tlev, st, min_samples))
write.csv(pluses, sprintf('%s/Part_1/figures/Correlations/%s/%s/%s_%s_%s-minSamps_signif_correlations_pluses.csv', home_dir, tlev, st, tlev, st, min_samples))

# if (nrow(otu_table(gs))==67) {
#   write.csv(res.signif, sprintf('%s/Part_1/figures/Correlations/signif_correlations.csv', home_dir))
#   write.csv(pluses, sprintf('%s/Part_1/figures/Correlations/signif_correlations_pluses.csv', home_dir))
# } else if (nrow(otu_table(gs))==265) {
#   write.csv(res.signif, sprintf('%s/Part_1/figures/Correlations/signif_correlations_rare_genera.csv', home_dir))
#   write.csv(pluses, sprintf('%s/Part_1/figures/Correlations/signif_correlations_pluses_rare_genera.csv', home_dir))
# }




# ************************************************************************************************************** #
# **************** #
#### Make heatmap of OTUs for each stomatotype group and for all samples, color those squares that change correlation significance ####
cor.otus.both <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/all/%s_all_OTUs_third-of-samples_signif_correlations.csv', 
                                    home_dir, tlev, tlev), header=T, row.names=1, sep=',')
cor.otus.stom1 <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/Stomatotype_1/%s_Stomatotype_1_OTUs_third-of-samples_signif_correlations.csv', 
                                    home_dir, tlev, tlev), header=T, row.names=1, sep=',')
cor.otus.stom2 <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/Stomatotype_2/%s_Stomatotype_2_OTUs_third-of-samples_signif_correlations.csv', 
                                    home_dir, tlev, tlev), header=T, row.names=1, sep=',')

cor.all.both <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/all/%s_all_OTUs_third-of-samples_all-values_correlations.csv', 
                                    home_dir, tlev, tlev), header=T, row.names=1, sep=',')
cor.all.stom1 <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/Stomatotype_1/%s_Stomatotype_1_OTUs_third-of-samples_all-values_correlations.csv', 
                                     home_dir, tlev, tlev), header=T, row.names=1, sep=',')
cor.all.stom2 <- read.delim(sprintf('%s/Part_1/figures/Correlations/%s/Stomatotype_2/%s_Stomatotype_2_OTUs_third-of-samples_all-values_correlations.csv', 
                                     home_dir, tlev, tlev), header=T, row.names=1, sep=',')
# **************** #

# diff.both.stom1 <- matrix(NA, nrow=nrow(cor.otus.both), ncol=ncol(cor.otus.both))
# diff.both.stom2 <- matrix(NA, nrow=nrow(cor.otus.both), ncol=ncol(cor.otus.both))
# diff.stom1.stom2 <- matrix(NA, nrow=nrow(cor.otus.both), ncol=ncol(cor.otus.both))
# rownames(diff.both.stom1) <- rownames(diff.both.stom2) <- rownames(diff.stom1.stom2) <- rownames(cor.otus.both)
# colnames(diff.both.stom1) <- colnames(diff.both.stom2) <- colnames(diff.stom1.stom2) <- colnames(cor.otus.both)
# 
# for (i in rownames(cor.otus.both)) {
#   for (j in colnames(cor.otus.both)) {
#     
#     bo <- cor.otus.both[i,j]
#     s1 <- cor.otus.stom1[i,j]
#     s2 <- cor.otus.stom2[i,j]
#     
#     if ( !(is.na(bo) | is.na(s1)) & ((bo > 0 & s1 < 0) | (bo < 0 & s1 > 0)) ) {
#       diff.both.stom1[i,j] <- 1#"green" # neither na, opposite signs
#     } else if ( is.na(bo) & !is.na(s1) & s1 > 0 ) { diff.both.stom1[i,j] <- 2#"pink" # s1 becomes pos cor
#     } else if ( is.na(bo) & !is.na(s1) & s1 < 0 ) { diff.both.stom1[i,j] <- 3#"purple" # s1 becomes neg cor
#     } else if ( ! is.na(bo) & is.na(s1) ) { diff.both.stom1[i,j] <- 4#"gray" # s1 no longer signif
#     } else { diff.both.stom1[i,j] <- 0}#"white" }
#     
#     if ( !(is.na(bo) | is.na(s2)) & ((bo > 0 & s2 < 0) | (bo < 0 & s2 > 0)) ) { 
#       diff.both.stom2[i,j] <- 1#"green" # neither na, opposite signs
#     } else if ( is.na(bo) & !is.na(s2) & s2 > 0 ) { diff.both.stom2[i,j] <- 2#"pink" # s2 becomes pos cor
#     } else if ( is.na(bo) & !is.na(s2) & s2 < 0 ) { diff.both.stom2[i,j] <- 3#"purple" # s2 becomes neg cor
#     } else if ( ! is.na(bo) & is.na(s2) ) { diff.both.stom2[i,j] <- 4#"grey" # s2 no longer signif
#     } else { diff.both.stom2[i,j] <- 0}#"white" }
#     
#     if ( !(is.na(s1) | is.na(s2)) & ((s1 > 0 & s2 < 0) | (s1 < 0 & s2 > 0)) ) { 
#       diff.stom1.stom2[i,j] <- "green" # neither na, opposite signs
#     } else if ( is.na(s1) & !is.na(s2) & s2 > 0 ) { diff.stom1.stom2[i,j] <- "pink" # s1 insignif, s2 pos cor
#     } else if ( is.na(s1) & !is.na(s2) & s2 < 0 ) { diff.stom1.stom2[i,j] <- "purple" # s1 insignif, s2 neg cor
#     } else if ( (!is.na(s1) & is.na(s2)) | (!is.na(s2) & is.na(s1)) ) { diff.stom1.stom2[i,j] <- "grey" # signif in 1, insignif in other
#     } else { diff.stom1.stom2[i,j] <- "white" }
#   }
# }
# 
# cosig <- as.matrix(cor.all.stom2)
# cosig <- cosig[rev(rownames(cosig)), rev(colnames(cosig))]
# cor.plus <- as.matrix(cor.otus.stom2)
# cor.plus[!is.na(cor.plus)] <- '+'
# cor.plus[is.na(cor.plus)] <- ''
# cor.plus <- cor.plus[rev(rownames(cor.plus)), rev(colnames(cor.plus))]
# 
# diff.both.stom1 <- diff.both.stom1[rev(rownames(diff.both.stom1)), rev(colnames(diff.both.stom1))]
# diff.both.stom2 <- diff.both.stom2[rev(rownames(diff.both.stom2)), rev(colnames(diff.both.stom2))]
# diff.stom1.stom2 <- diff.stom1.stom2[rev(rownames(diff.stom1.stom2)), rev(colnames(diff.stom1.stom2))]
# 
# 
# cor.m.otus <- cbind( melt( cosig ), melt( cor.plus ), melt( diff.both.stom2 ) )
# cor.m.otus <- cor.m.otus[,c(1,2,3,6,9)]
# colnames(cor.m.otus)[4] <- 'signif'
# colnames(cor.m.otus)[5] <- 'borderCol'
# cor.m.otus.lower <- subset( cor.m.otus[ upper.tri(cosig), ] )
# 
# ggplot(data = cor.m.otus, aes(x=Var2, y=Var1, fill=value)) + 
#   # geom_tile(aes(color=as.factor(borderCol), width=0.7, height=0.7), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==3,], size=1, color="grey") + 
#   geom_tile(color="white", data=cor.m.otus.lower) + 
#   geom_tile(aes(width=0.9, height=0.9), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==0,], size=0.7, color="white") + 
#   geom_tile(aes(width=0.9, height=0.9), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==1,], size=0.7, color="green") + 
#   geom_tile(aes(width=0.9, height=0.9), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==2,], size=0.7, color="cyan") + 
#   geom_tile(aes(width=0.9, height=0.9), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==3,], size=0.7, color="purple") + 
#   geom_tile(aes(width=0.9, height=0.9), data=cor.m.otus.lower[cor.m.otus.lower$borderCol==4,], size=0.7, color="grey") + 
#   
#   # scale_fill_discrete(colours = c("white", "green", "pink","purple","gray"), values = c(0,1,2,3,4)) +
#   # scale_fill_manual(values=c("white", "green", "cyan", "purple", "gray")) +
#   geom_text(aes(label=signif), size = 4.5, data=cor.m.otus.lower) + 
#   scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Cor") +
#   # guides(fill=F)+
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
#         axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
#         plot.title = element_text(hjust=0.5)) +
#   ggtitle('title') + xlab('') + ylab('')
# 
# 
# 
# Or better to check for sites that changed same in both (green), only in s1 (pink), only in s2 (purple), etc
# 
# 
# # check specific correlations
# color <- "purple"
# # tableToCheck <- diff.both.stom1
# # tableToCheck <- diff.both.stom2
# tableToCheck <- diff.stom1.stom2
# 
# greens <- apply(tableToCheck, 2, function(x) color %in% x)
# greens <- greens[greens]
# green.list <- list()
# for (g in names(greens)) {
#   goods <- sapply(tableToCheck[g,], function(x) x==color)
#   green.list[[g]] <- goods[goods & !is.na(goods)]
# }

cor.test( as.numeric(otu_table(SLL1)["Neisseria",]), 
          as.numeric(otu_table(SLL1)["Haemophilus",]))

cor.test( as.numeric(otu_table(SLL1)["Neisseria", sample_names(SLL1)[ sample_data(SLL1)[,"Stomatotype"]==1 ]]), 
          as.numeric(otu_table(SLL1)["Haemophilus", sample_names(SLL1)[ sample_data(SLL1)[,"Stomatotype"]==1 ]]))

cor.test( as.numeric(otu_table(SLL1)["Neisseria", sample_names(SLL1)[ sample_data(SLL1)[,"Stomatotype"]==2 ]]), 
          as.numeric(otu_table(SLL1)["Haemophilus", sample_names(SLL1)[ sample_data(SLL1)[,"Stomatotype"]==2 ]]))
# **************** #



t.diff.both.stom1 <- matrix(NA, nrow=nrow(cor.all.both), ncol=ncol(cor.all.both))
t.diff.both.stom2 <- matrix(NA, nrow=nrow(cor.all.both), ncol=ncol(cor.all.both))
t.diff.stom1.stom2 <- matrix(NA, nrow=nrow(cor.all.both), ncol=ncol(cor.all.both))
rownames(t.diff.both.stom1) <- rownames(t.diff.both.stom2) <- rownames(t.diff.stom1.stom2) <- rownames(cor.all.both)
colnames(t.diff.both.stom1) <- colnames(t.diff.both.stom2) <- colnames(t.diff.stom1.stom2) <- colnames(cor.all.both)

min.diff <- 0.20

for (i in rownames(cor.all.both)) {
  for (j in colnames(cor.all.both)) {
    
    bo <- cor.all.both[i,j]
    s1 <- cor.all.stom1[i,j]
    s2 <- cor.all.stom2[i,j]
    
    bo.s1 <- s1 - bo
    bo.s2 <- s2 - bo
    s1.s2 <- s1 - s2
    
    if (bo.s1 >= min.diff) { t.diff.both.stom1[i,j] <- 1}#; print(c(i,j, bo.s1));}
    else if (bo.s1 <= -min.diff) {t.diff.both.stom1[i,j] <- 2}#; print(c(i,j, bo.s1));}
    else t.diff.both.stom1[i,j] <- 0
    if (bo.s2 >= min.diff) t.diff.both.stom2[i,j] <- 1
    else if (bo.s2 <= -min.diff) t.diff.both.stom2[i,j] <- 2
    else t.diff.both.stom2[i,j] <- 0
    if (s1.s2 >= min.diff) t.diff.stom1.stom2[i,j] <- 1
    else if (s1.s2 <= -min.diff) t.diff.stom1.stom2[i,j] <- 2
    else t.diff.stom1.stom2[i,j] <- 0

    
  }
}

# stom <- "Stomatotype_1"
# # stom <- "Stomatotype_2"
# # stom <- "both_Stomatotypes"
# 
# ds <- "bo.s1"
# # ds <- "bo.s2"
# # ds <- "s1.s2"
# 
# if (stom == "both_Stomatotypes") {
#   ca <- cor.all.both
#   cp <- cor.otus.both
# } else if (stom == "Stomatotype_1") {
#   ca <- cor.all.stom1
#   cp <- cor.otus.stom1
# } else if (stom == "Stomatotype_2") {
#   ca <- cor.all.stom2
#   cp <- cor.otus.stom2
# }
# 
# if (ds == "bo.s1") { td <- t.diff.both.stom1[rev(rownames(t.diff.both.stom1)), rev(colnames(t.diff.both.stom1))]
# } else if (ds == "bo.s2") { td <- t.diff.both.stom2[rev(rownames(t.diff.both.stom2)), rev(colnames(t.diff.both.stom2))]
# } else if (ds == "s1.s2") { td <- t.diff.stom1.stom2[rev(rownames(t.diff.stom1.stom2)), rev(colnames(t.diff.stom1.stom2))]
# }
# 
# cosig <- as.matrix(ca)
# cosig <- cosig[rev(rownames(cosig)), rev(colnames(cosig))]
# cor.plus <- as.matrix(cp)
# cor.plus[!is.na(cor.plus)] <- '+'
# cor.plus[is.na(cor.plus)] <- ''
# cor.plus <- cor.plus[rev(rownames(cor.plus)), rev(colnames(cor.plus))]
# 
# cor.m.otus <- cbind( melt( cosig ), melt( cor.plus ), melt( td ) )
# cor.m.otus <- cor.m.otus[,c(1,2,3,6,9)]
# colnames(cor.m.otus)[4] <- 'signif'
# colnames(cor.m.otus)[5] <- 'borderCol'
# cor.m.otus.lower <- subset( cor.m.otus[ upper.tri(cosig), ] )
# 
# ggplot(data = cor.m.otus, aes(x=Var2, y=Var1, fill=value)) + 
#   geom_tile(color="white", data=cor.m.otus.lower) + 
#   geom_tile( data=cor.m.otus.lower[cor.m.otus.lower$borderCol==0,], size=0.5, color="white") + #, aes(width=0.9, height=0.9)) 
#   geom_tile( data=cor.m.otus.lower[cor.m.otus.lower$borderCol==1,], size=0.5, color="red") + #, aes(width=0.9, height=0.9)) 
#   geom_tile( data=cor.m.otus.lower[cor.m.otus.lower$borderCol==2,], size=0.5, color="blue") + #, aes(width=0.9, height=0.9)) 
#   geom_text(aes(label=signif), size = 4.5, data=cor.m.otus.lower) + 
#   scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Cor") +
#   # guides(fill=F)+
#   theme_minimal() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
#         axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
#         plot.title = element_text(hjust=0.5)) +
#   ggtitle(sprintf('OTUs cor() Genus - %s\nFull dataset and %s differ by at least %s highlighted', stom,stom,min.diff)) +
#   xlab('') + ylab('')



library(gridExtra)

cosig <- as.matrix(cor.all.stom1)
cosig <- cosig[rev(rownames(cosig)), rev(colnames(cosig))]
cor.plus <- as.matrix(cor.otus.stom1)
cor.plus[!is.na(cor.plus)] <- '+'
cor.plus[is.na(cor.plus)] <- ''
cor.plus <- cor.plus[rev(rownames(cor.plus)), rev(colnames(cor.plus))]

cor.m.otus.1 <- cbind( melt( cosig ), melt( cor.plus ), melt( t.diff.stom1.stom2[rev(rownames(t.diff.stom1.stom2)), rev(colnames(t.diff.stom1.stom2))] ) )
cor.m.otus.1 <- cor.m.otus.1[,c(1,2,3,6,9)]
colnames(cor.m.otus.1)[4] <- 'signif'
colnames(cor.m.otus.1)[5] <- 'borderCol'
cor.m.otus.1.lower <- subset( cor.m.otus.1[ upper.tri(cosig), ] )


cosig <- as.matrix(cor.all.stom2)
cosig <- cosig[rev(rownames(cosig)), rev(colnames(cosig))]
cor.plus <- as.matrix(cor.otus.stom2)
cor.plus[!is.na(cor.plus)] <- '+'
cor.plus[is.na(cor.plus)] <- ''
cor.plus <- cor.plus[rev(rownames(cor.plus)), rev(colnames(cor.plus))]

cor.m.otus.2 <- cbind( melt( cosig ), melt( cor.plus ), melt( t.diff.stom1.stom2[rev(rownames(t.diff.stom1.stom2)), rev(colnames(t.diff.stom1.stom2))] ) )
cor.m.otus.2 <- cor.m.otus.2[,c(1,2,3,6,9)]
colnames(cor.m.otus.2)[4] <- 'signif'
colnames(cor.m.otus.2)[5] <- 'borderCol'
cor.m.otus.2.lower <- subset( cor.m.otus.2[ upper.tri(cosig), ] )


o1 <- ggplot(data = cor.m.otus.1, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color="white", data=cor.m.otus.1.lower) + 
  geom_tile( data=cor.m.otus.1.lower[cor.m.otus.1.lower$borderCol==0,], size=0.7, color="white") + #, aes(width=0.9, height=0.9)) 
  geom_tile( data=cor.m.otus.1.lower[cor.m.otus.1.lower$borderCol==1,], size=0.7, color="red") + #, aes(width=0.9, height=0.9)) 
  geom_tile( data=cor.m.otus.1.lower[cor.m.otus.1.lower$borderCol==2,], size=0.7, color="blue") + #, aes(width=0.9, height=0.9)) 
  geom_text(aes(label=signif), size = 4.5, data=cor.m.otus.1.lower) + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Pearson") +
  # guides(fill=F)+
  theme_minimal() + 
  # theme(axis.text.x = element_blank(),
  #       axis.text.y = element_blank(),
  #       plot.title = element_text(hjust=0.5, size=17)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 10, hjust = 1),
        plot.title = element_text(hjust=0.5)) +
  # guides(fill = guide_legend(title.theme=element_text(size=15, angle=0, face="bold"), 
  #                            label.theme=element_text(size=13, angle=0)) ) +
  # ggtitle(sprintf('OTUs cor() Genus - Stomatotype_1\nStomatotype_1 and Stomatotype_2 differ by at least %s highlighted', min.diff)) +
  ggtitle("Stomatotype 1") +
  xlab('') + ylab('')

o2 <- ggplot(data = cor.m.otus.2, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color="white", data=cor.m.otus.2.lower) + 
  geom_tile( data=cor.m.otus.2.lower[cor.m.otus.2.lower$borderCol==0,], size=0.7, color="white") + #, aes(width=0.9, height=0.9)) 
  geom_tile( data=cor.m.otus.2.lower[cor.m.otus.2.lower$borderCol==1,], size=0.7, color="blue") + #, aes(width=0.9, height=0.9)) 
  geom_tile( data=cor.m.otus.2.lower[cor.m.otus.2.lower$borderCol==2,], size=0.7, color="red") + #, aes(width=0.9, height=0.9)) 
  geom_text(aes(label=signif), size = 4.5, data=cor.m.otus.2.lower) + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, space="Lab",name="Pearson") +
  # guides(fill=F)+
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust=0.5, size=17)) +
  # ggtitle(sprintf('OTUs cor() Genus - Stomatotype_2\nStomatotype_1 and Stomatotype_2 differ by at least %s highlighted', min.diff)) +
  ggtitle("Stomatotype 2") +
  xlab('') + ylab('')

grid.arrange(o1, o2, nrow=1, ncol=2)

# ************************************************************************************************************** #








# ************* #
# to look at tax table entries for genera with signif pos or neg correlations to a certain question:
ques <- "Stomatotype"

signif_gen <- names( plus.Que_Otu[ ques, plus.Que_Otu[ques,]=='+' ] )
pos_signif_gen <- signif_gen[ res.Que_Otu[ques,signif_gen] > 0 ]
neg_signif_gen <- signif_gen[ res.Que_Otu[ques,signif_gen] < 0 ]

tax_table(SLL)[pos_signif_gen]
tax_table(SLL)[neg_signif_gen]



# ************* #

# to look at which questions/genera have few signif cors
for (i in names(goodmin)) {
  if (length(colnames(pluses)[pluses[i,]=='+']) < 2) {
    print(c(i, colnames(pluses)[pluses[i,]=='+']))
  }
  # print(c(i, length(res.Que_Que[i, names(plus.Que_Que[i,plus.Que_Que[i,]=='+'])])))
  # print(res.Que_Que[i, names(plus.Que_Que[i,plus.Que_Que[i,]=='+'])])
  # print("")
}

# to look at which questions/genera have ZERO signif cors
for (i in rownames(res.matrix)) {
  if (!i %in% names(goodmin)) {
    print(i)
  }
}








# Significant positive (or negative) correlations with Streptococcus (or whichever genus)
rownames( plus.Otu_Otu[ (plus.Otu_Otu["Streptococcus", ] == '+') & (res.Otu_Otu["Streptococcus",] > 0), ] )


## Which OTUs or questions have significant correlations with the indicated OTU/Question
oq <- "Streptococcus"
res.matrix[colnames(ps.matrix.adj)[pluses[oq,]=="+"], oq]
ps.matrix.adj[colnames(ps.matrix.adj)[pluses[oq,]=="+"], oq]
colnames(ps.matrix.adj)[pluses[oq,]=="+"]





# Explore those OTUs with significant correlations
tax.cor.signif <- tax.table[ colnames(res.matrix[otus.cont,otus.cont]), ] #res.matrix[otus.cont,otus.cont] === res.Otu_Otu
for (i in 1:length(rownames(tax.cor.signif))) {
  print( paste(tax.cor.signif[i,1],tax.cor.signif[i,2],tax.cor.signif[i,3],tax.cor.signif[i,4],tax.cor.signif[i,5],rownames(tax.cor.signif)[i], sep=";") )
}









# ****************************************************************************************************************** #
###### t-tests ######
# ****************************************************************************************************************** #

# questions with categories
group_qs <- c("Q2","Q8","Q15","Q15.1","Q18","Q19","Q19.1","Q19.3","Q19.5","Q19.6","Q19.8","Q19.9","Q20.1","Q21","Q27",
              "Q30","Q31","Q32","Q33","Q34","Q35","Q36","Q37","Q38","Q41","Q43","Q44","Q45","Q46","Q47","Q48",
              "Q49","Q50","Q51","Q52","Q53","Q54","Q1.1","Q14.1","Q14.2","Q22.1","Q23.1","Q24.1","Q25.1","Q26.1",
              "Q39.1","Q40.1","Q28.11","Q28.13","Q28.15","Q28.17","Q28.19","Q28.21","Q28.23","Q28.25","Q28.27","Q28.29",
              "Q28.31","Stomatotype","Stomatotype_CORE","Diversity_group_Div.Shannon","Diversity_group_Div.Simpson",
              "Diversity_group_Weighted_Unifrac","Diversity_group_Unweighted_Unifrac","Diversity_group_Faiths.PD",
              "Diversity_group_Species_Richness","BMI_group","BMI_official",group_water_data)#,"Water_hardness_group",
group_qs <- as.matrix(sample_data(SLL1)[, group_qs])

#numerical columns
only_cont <- c("Q3","Q16","Q17","Q19.2","Q19.4","Q19.7","Q19.10","Q28","Q28.1","Q28.2","Q28.3","Q28.4",
               "Q28.5","Q28.6","Q28.7","Q28.8","Q28.9","Q28.10","Q29","Q29.1","Q29.2","Q29.3","Q00.PH",
               "Socioeconomic","Q3.1","Q22.2","Q23.2","Q24.2","Q25.2","Q26.2","Q39.2","Q40.2","Q28.12",
               "Q28.14","Q28.16","Q28.18","Q28.20","Q28.22","Q28.24","Q28.26","Q28.28","Q28.30","Q28.32",
               "Q16.1","Div.Observed","Div.Chao1","Div.ACE","Div.Shannon","Div.Simpson","Div.InvSimpson",
               "Div.Fisher","Faiths.PD","Species_Richness","Num_OTUs","Stomatotype","Stomatotype_CORE",
               "Population","Age","Weighted_Unifrac","Unweighted_Unifrac",cont_water_data)#,"Water_hardness","Neisseria_Prevotella","Neisseria_Veillonella",
               #"Haemophilus_Prevotella","Haemophilus_Veillonella","Streptococcus_Prevotella","Prevotella_Bacteroides",


tlev <- "Genus"; otutab <- otu_table(SLL); otutab_rel <- otu_table(SLL1)
# tlev <- "Phylum"; otutab <- tlev_otus[[tlev]]; otutab_rel <- tlev_otus_rel[[tlev]]
  
  
# wh.ttest <- genefilter_sample(SLL, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(SLL))
# slltt.raw <- prune_taxa(wh.ttest, SLL)
# slltt.norm <- prune_taxa(wh.ttest, SLL1)
# 
# numerical_data.raw <- cbind(t(otu_table(slltt.raw)), sample_data(slltt.raw)[,only_cont])
# numerical_data.norm <- cbind(t(otu_table(slltt.norm)), sample_data(slltt.norm)[,only_cont])

numerical_data.raw <- cbind( t(otutab), sample_data(SLL)[,only_cont])
numerical_data.norm <- cbind( t(otutab_rel), sample_data(SLL1)[,only_cont])
num_rows <- rownames(numerical_data.raw)
numerical_data.raw[numerical_data.raw=='No Sabe/No Contesta'] <- NA
numerical_data.norm[numerical_data.norm=='No Sabe/No Contesta'] <- NA
numerical_data.raw <- apply(numerical_data.raw, 2, as.numeric)
numerical_data.norm <- apply(numerical_data.norm, 2, as.numeric)
rownames(numerical_data.raw) <- num_rows
rownames(numerical_data.norm) <- num_rows


# for column names, get Q--category
get_dim_names <- function(x) {
  categories <- unique(group_qs[,x])[ ! is.na(unique(group_qs[,x])) ]
  categories <- categories[categories != "No Sabe/No Contesta"]
  return( sapply(categories, function(y) paste(x, y, sep = '--' ) ) )
}

tt_cols <- unname(unlist( sapply(colnames(group_qs), get_dim_names) ))
tt_rows <- colnames(numerical_data.raw)

tt.p.greater <-  matrix(NA, nrow=length(tt_rows), ncol=length(tt_cols))
colnames(tt.p.greater) <- tt_cols
rownames(tt.p.greater) <- tt_rows
tt.p.less <-  matrix(NA, nrow=length(tt_rows), ncol=length(tt_cols))
colnames(tt.p.less) <- tt_cols
rownames(tt.p.less) <- tt_rows
fc.matrix <-  matrix(NA, nrow=length(tt_rows), ncol=length(tt_cols))
colnames(fc.matrix) <- tt_cols
rownames(fc.matrix) <- tt_rows

min.p.val <- 0.05
min.FC <- 1.5

# ************************** #
for (n in tt_rows) {
  for (q in colnames(group_qs)) {
    if (n == q) {
      next
    }
    
    # different categories for the given question
    q_cats <- unique(group_qs[ , q ])[ ! is.na(unique(group_qs[ , q ])) ]
    q_cats <- q_cats[q_cats!="No Sabe/No Contesta"] # dont care about these
    
    for (q_cat in q_cats) {
      if (class(group_qs[group_qs[,q] == q_cat, ]) == "character") {
        # for when there is only one sample with the given value, cannot get rowname, since it converts from matrix to character vector
        q_cat_rows <- rownames( group_qs[group_qs[,q] == q_cat, , drop = F] )
        next
        
        # } else if (q_cat == "No Sabe/No Contesta") {
        #   next # dont care about these
        
      } else {
        #get vectors of samples with given 
        q_cat_rows <- rownames( group_qs[group_qs[,q] == q_cat, ] )
        q_cat_rows <- q_cat_rows[! is.na(q_cat_rows)]
        q_cat_nots <- rownames( group_qs[group_qs[,q] %in% q_cats[q_cats != q_cat], ] )
        q_cat_nots <- q_cat_nots[! is.na(q_cat_nots)]
      }
      
      if (sum( ! is.na(numerical_data.norm[ q_cat_rows, n ])) <= 1 | sum( ! is.na(numerical_data.norm[ q_cat_nots, n ])) <= 1 ) {
        # if only 1 or none of samples has q_cat or non-q_cat, skip since t-test will say not enough observations
        next
      }

      # run t-test and store p-value
      tp.greater <- t.test(numerical_data.norm[ q_cat_rows, n ],
                           numerical_data.norm[ q_cat_nots, n ],
                           alternative = "greater")$p.value
      
      tp.less <- t.test(numerical_data.norm[ q_cat_rows, n ],
                        numerical_data.norm[ q_cat_nots, n ],
                        alternative = "less")$p.value
      
      # for fold change, will use non-normalized values because the normalization significantly 
      # changes absoluate values column-wise, affecting fold change values
      fc_vals <- foldchange(mean(numerical_data.raw[ q_cat_rows, n ], na.rm = T) ,
                            mean(numerical_data.raw[ q_cat_nots, n ], na.rm = T) )
      
      cname <- paste(q, q_cat, sep = '--')
      tt.p.greater[n, cname] <- tp.greater
      tt.p.less[n, cname] <- tp.less
      fc.matrix[n, cname] <- foldchange2logratio(fc_vals, base=2)
    }
  }
}
# ************************** #

# for those values that were skipped because of lack of data or whatever
tt.p.greater[ is.na(tt.p.greater) ] <- 1
tt.p.less[ is.na(tt.p.less) ] <- 1

# adjust for multiple testing
tt.p.greater.adj <- apply(tt.p.greater, 2, p.adjust, method='bonferroni', n=ncol(numerical_data.norm))
tt.p.less.adj <- apply(tt.p.less, 2, p.adjust, method='bonferroni', n=ncol(numerical_data.norm))

#at least 1 good p value within question response
mins.greater <- apply(tt.p.greater.adj,2,min)
goodPs.greater <- mins.greater[mins.greater < min.p.val] #-log10(0.05) == 1.3013
mins.less <- apply(tt.p.less.adj,2,min)
goodPs.less <- mins.less[mins.less < min.p.val] #-log10(0.05) == 1.3013 
#at least 1 good p within phenotype
tmins.greater <- apply(tt.p.greater.adj,1,min)
tgoodPs.greater <- tmins.greater[tmins.greater < min.p.val]
tmins.less <- apply(tt.p.less.adj,1,min)
tgoodPs.less <- tmins.less[tmins.less < min.p.val]

# a matrix with only those rows and columns that had a least 1 significant p-val
# use drop=F for those cases in which only 1 row or column had at least 1 significant p-value...maintains structure as a matrix
good.tt.p.greater.adj <- tt.p.greater.adj[names(tgoodPs.greater), names(goodPs.greater), drop = F]
good.tt.p.less.adj <- tt.p.less.adj[names(tgoodPs.less), names(goodPs.less), drop = F]

# any non-significant p-vals will be a blank entry
only.good.tt.p.greater.adj <- good.tt.p.greater.adj
only.good.tt.p.greater.adj[only.good.tt.p.greater.adj >= min.p.val] <- ''
only.good.tt.p.less.adj <- good.tt.p.less.adj
only.good.tt.p.less.adj[only.good.tt.p.less.adj >= min.p.val] <- ''

# dim(only.good.tt.p.less.adj)
# dim(only.good.tt.p.greater.adj)
print(c("less", dim(only.good.tt.p.less.adj)))
print(c("greater", dim(only.good.tt.p.greater.adj)))

write.csv(only.good.tt.p.greater.adj, sprintf("%s/Part_1/figures/ttests/%s/Abundance/%s_ttest_signif_greater.csv", home_dir, tlev, tlev) )
write.csv(only.good.tt.p.less.adj, sprintf("%s/Part_1/figures/ttests/%s/Abundance/%s_ttest_signif_less.csv", home_dir, tlev, tlev) )
write.csv(fc.matrix, sprintf("%s/Part_1/figures/ttests/%s/Abundance/%s_foldchange.csv", home_dir, tlev, tlev) )

# ************************** #



## get volcano plots
full_pvals_greater <- tt.p.greater.adj
full_pvals_less <- tt.p.less.adj

#make melted table with questions, phenotypes, pvals, FC
pv_fc.m <- cbind(melt(full_pvals_greater), melt(full_pvals_less), melt(fc.matrix))
pv_fc.m <- pv_fc.m[,c(1,2,3,6,9)]
colnames(pv_fc.m) <- c("Cont_data", "Categ_question", 'p_greater', 'p_less', 'FC')

pv_fc.m.final <- pv_fc.m[, c("Cont_data", "Categ_question", "FC")]
pv_fc.m.final["p_val"] <- ifelse( pv_fc.m.final$FC >= 0, pv_fc.m$p_greater, pv_fc.m$p_less )
# pv_fc.m.final["combo"] <- paste(pv_fc.m.final$Phenotype, pv_fc.m.final$Question, sep = "_")
pv_fc.m.final["group"] <- ifelse( pv_fc.m.final$p_val < min.p.val & abs(pv_fc.m.final$FC) > min.FC, "Signif&FC", 
                                  ifelse( abs(pv_fc.m.final$FC) > min.FC, "FC", 
                                          ifelse( pv_fc.m.final$p_val < min.p.val, "Signif", "NotSignif" )) )
# some FC values are Inf or -Inf if one side had a mean of 0
pv_fc.m.final <- pv_fc.m.final[ ! is.na(pv_fc.m.final$FC) & is.finite(pv_fc.m.final$FC), ]

# for x, y limits
max_FC <- max(pv_fc.m.final$FC)
min_FC <- min(pv_fc.m.final$FC)
max_log10_p <- max(-log10(pv_fc.m.final[ is.finite(-log10(pv_fc.m.final$p_val)), "p_val"])) #bc some p_val are 0, and thus Inf -log10(p)

vp <- ggplot(data = pv_fc.m.final, aes(x=FC, y=-log10(p_val), color=group)) + 
  geom_point() +
  scale_color_manual(name = "Significance",
                     values = c("FC" = "orange",
                                "NotSignif" = "black",
                                "Signif" = "red",
                                "Signif&FC" = "blue"),
                     labels = c("FC", "NotSignif", "Signif","Signif&FC")) +
  xlim(c( (-1)*max(abs(min_FC),max_FC), max(abs(min_FC),max_FC) )) +
  ylim(c(0, max_log10_p)) +
  # ylim(c(0, max(-log10( pv_fc.m.final[ pv_fc.m.final$group=="Signif&FC", "p_val"] )) + 1)) +
  geom_vline(xintercept = c(-1*min.FC, min.FC), linetype = "dotted") +
  geom_hline(yintercept = -log10(min.p.val), linetype = "dotted") +
  xlab("log2( fold change )") + ylab("-log10( p-value )") +
  # geom_text(aes(x=vp_labels_x, y=vp_labels_y, label=vp_labels), size = 3.5, color="black") +
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle( sprintf('Volcano plot for categorical questions vs continuous data') )
vp

# sort pv_fc.m.final by group labels, then by FC, to print them to a file
sorted_volcano <- pv_fc.m.final[ rev(order(pv_fc.m.final$group, -pv_fc.m.final$FC)), ]
sorted_volcano <- sorted_volcano[abs(sorted_volcano$FC) > min.FC | sorted_volcano$p_val < min.p.val, ]
write.csv(sorted_volcano, sprintf("%s/Part_1/figures/ttests/%s/Abundance/%s_sorted_volcano.csv", home_dir, tlev, tlev) )
# ****************************************************************************************************************** #





# ****************************************************************************************************************** #
# Make table of Presence/Absense for each genus in each sample ####
# ****************************************************************************************************************** #

tlev <- "Genus"; otutab <- otu_table(SLL)
# tlev <- "Phylum"; otutab <- tlev_otus[[tlev]]

P_A_table <- apply( otutab, 2, function(x) ifelse(x==0, 0, 1) )

# the vector cont_bin is found in the section "Heatmap of correlations for continuous/binary responses"


pa_g <- matrix(NA, nrow = nrow(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ]), ncol = length(cont_bin))
rownames(pa_g) <- rownames(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ])
colnames(pa_g) <- cont_bin
pa_l <- matrix(NA, nrow = nrow(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ]), ncol = length(cont_bin))
rownames(pa_l) <- rownames(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ])
colnames(pa_l) <- cont_bin
pa_fc <- matrix(NA, nrow = nrow(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ]), ncol = length(cont_bin))
rownames(pa_fc) <- rownames(P_A_table[ ! rowSums(P_A_table) %in% c(1,1319), ])
colnames(pa_fc) <- cont_bin

for (gen in rownames(pa_g)) {
  pres <- colnames( P_A_table[ , P_A_table[ gen, ]==1 ] )
  abse <- colnames( P_A_table[ , P_A_table[ gen, ]==0 ] )
  print(gen)
  
  for (que in colnames(pa_g)) {
    # print(c(pres, que))
    pq <- as.numeric( as.matrix(sample_data(SLL))[ pres, que ] )
    aq <- as.numeric( as.matrix(sample_data(SLL))[ abse, que ] )
    pq <- pq[ !is.na(pq) ]
    aq <- aq[ !is.na(aq) ]
    
    if (length(pq) < 2 | length(aq) < 2) {
      print(c(gen, que, length(pq), length(aq)))
      next
    }
    
    tt.g <- t.test( pq, aq, alternative = "greater" )
    tt.l <- t.test( pq, aq, alternative = "less" )
    fc_vals <- foldchange( mean(pq), mean(aq) )
    
    pa_g[ gen, que ] <- tt.g$p.value
    pa_l[ gen, que ] <- tt.l$p.value
    pa_fc[gen, que ] <- foldchange2logratio(fc_vals, base=2)
  }
}

# for those values that were skipped because too few samples for a given instance
pa_g[ is.na(pa_g) ] <- 1
pa_l[ is.na(pa_l) ] <- 1

# adjust for multiple testing
pa_g.adj <- apply(pa_g, 2, p.adjust, method='bonferroni', n=nrow(pa_g))
pa_l.adj <- apply(pa_l, 2, p.adjust, method='bonferroni', n=nrow(pa_l))

#at least 1 good p value within question response
mins.greater <- apply(pa_g.adj,2,min)
goodPs.greater <- mins.greater[mins.greater < min.p.val] #-log10(0.05) == 1.3013
mins.less <- apply(pa_l.adj,2,min)
goodPs.less <- mins.less[mins.less < min.p.val] #-log10(0.05) == 1.3013 
#at least 1 good p within phenotype
tmins.greater <- apply(pa_g.adj,1,min)
tgoodPs.greater <- tmins.greater[tmins.greater < min.p.val]
tmins.less <- apply(pa_l.adj,1,min)
tgoodPs.less <- tmins.less[tmins.less < min.p.val]

# a matrix with only those rows and columns that had a least 1 significant p-val
# use drop=F for those cases in which only 1 row or column had at least 1 significant p-value...maintains structure as a matrix
good.pa_g.adj <- pa_g.adj[names(tgoodPs.greater), names(goodPs.greater), drop = F]
good.pa_l.adj <- pa_l.adj[names(tgoodPs.less), names(goodPs.less), drop = F]

# any non-significant p-vals will be a blank entry
only.good.pa_g.adj <- good.pa_g.adj
only.good.pa_g.adj[only.good.pa_g.adj >= min.p.val] <- ''
only.good.pa_l.adj <- good.pa_l.adj
only.good.pa_l.adj[only.good.pa_l.adj >= min.p.val] <- ''

# dim(only.good.pa_l.adj)
# dim(only.good.pa_g.adj)
print(c("less", dim(only.good.pa_l.adj)))
print(c("greater", dim(only.good.pa_g.adj)))

write.csv(only.good.pa_g.adj, sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_signif_greater.csv", home_dir, tlev, tlev) )
write.csv(only.good.pa_l.adj, sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_signif_less.csv", home_dir, tlev, tlev) )
write.csv(pa_fc, sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_foldchange.csv", home_dir, tlev, tlev) )


####### Make volcano plot 
#make melted table with questions, genera, pvals, FC
pv_fc.m <- cbind(melt(pa_g.adj), melt(pa_l.adj), melt(pa_fc))
pv_fc.m <- pv_fc.m[,c(1,2,3,6,9)]
colnames(pv_fc.m) <- c(tlev, "Cont_data", 'p_greater', 'p_less', 'FC')

pv_fc.m.final <- pv_fc.m[, c(tlev, "Cont_data", "FC")]
pv_fc.m.final["p_val"] <- ifelse( pv_fc.m.final$FC >= 0, pv_fc.m$p_greater, pv_fc.m$p_less )
# pv_fc.m.final["combo"] <- paste(pv_fc.m.final$Phenotype, pv_fc.m.final$Question, sep = "_")
pv_fc.m.final["group"] <- ifelse( pv_fc.m.final$p_val < min.p.val & abs(pv_fc.m.final$FC) > min.FC, "Signif&FC", 
                                  ifelse( abs(pv_fc.m.final$FC) > min.FC, "FC", 
                                          ifelse( pv_fc.m.final$p_val < min.p.val, "Signif", "NotSignif" )) )
# some FC values are Inf or -Inf if one side had a mean of 0
pv_fc.m.final <- pv_fc.m.final[ ! is.na(pv_fc.m.final$FC) & is.finite(pv_fc.m.final$FC), ]

# for x, y limits
max_FC <- max(pv_fc.m.final$FC)
min_FC <- min(pv_fc.m.final$FC)
max_log10_p <- max(-log10(pv_fc.m.final[ is.finite(-log10(pv_fc.m.final$p_val)), "p_val"])) #bc some p_val are 0, and thus Inf -log10(p)

vp <- ggplot(data = pv_fc.m.final, aes(x=FC, y=-log10(p_val), color=group)) + 
  geom_point() +
  scale_color_manual(name = "Significance",
                     values = c("FC" = "orange",
                                "NotSignif" = "black",
                                "Signif" = "red",
                                "Signif&FC" = "blue"),
                     labels = c("FC", "NotSignif", "Signif","Signif&FC")) +
  xlim(c( (-1)*max(abs(min_FC),max_FC), max(abs(min_FC),max_FC) )) +
  ylim(c(0, max_log10_p)) +
  # ylim(c(0, max(-log10( pv_fc.m.final[ pv_fc.m.final$group=="Signif&FC", "p_val"] )) + 1)) +
  geom_vline(xintercept = c(-1*min.FC, min.FC), linetype = "dotted") +
  geom_hline(yintercept = -log10(min.p.val), linetype = "dotted") +
  xlab("log2( fold change )") + ylab("-log10( p-value )") +
  # geom_text(aes(x=vp_labels_x, y=vp_labels_y, label=vp_labels), size = 3.5, color="black") +
  theme(plot.title = element_text(hjust=0.5)) + 
  ggtitle( sprintf('Volcano plot for presence/absenece of %s vs continuous data', tlev) )
vp

# sort pv_fc.m.final by group labels, then by FC, to print them to a file
presence_volcano <- pv_fc.m.final[ rev(order(pv_fc.m.final$group, -pv_fc.m.final$FC)), ]
presence_volcano <- presence_volcano[abs(presence_volcano$FC) > min.FC | presence_volcano$p_val < min.p.val, ]
write.csv(presence_volcano, sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_volcano.csv", home_dir, tlev, tlev) )

# ****************************************************************************************************************** #










# *********************************************************** #
# Check # of schools in which each genus actually appears  ####
# *********************************************************** #
tlev <- "Genus"; otutab <- otu_table(SLL1)
# tlev <- "Phylum"; otutab <- tlev_otus_rel[[tlev]]

P_A_table <- apply( otutab, 2, function(x) ifelse(x==0, 0, 1) )


get_num_schools <- function(x) {
  samps <- colnames(P_A_table)[ P_A_table[ x, ]==1 ]
  as.character(as.matrix(unique(sample_data(SLL)[ samps, "Qsamples.1"])))
}

num_schools_tlev <- sapply( rownames(P_A_table), get_num_schools )
# order by length
num_schools_tlev <- num_schools_tlev[ order(sapply(num_schools_tlev,length)) ]#,decreasing=T) ]
dir.create(sprintf("%s/Part_1/figures/ttests/%s/Presence", home_dir, tlev), showWarnings = F)
capture.output(num_schools_tlev, file = sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_num_schools.csv", home_dir, tlev, tlev) )

# get the region names in which each genus appears
get_region_tlev <- function(region) {
  lapply(num_schools_tlev, function(x) {
    temp_reg <- lapply(region, function(y) as.numeric(x) %in% y)
    names(temp_reg[ lapply(temp_reg, sum) > 0])
  })
}
num_cities_tlev <- get_region_tlev(cities)
num_provinces_tlev <- get_region_tlev(provinces)
num_communities_tlev <- get_region_tlev(comunidades)

# Now will make a table of presence/absence of each genus in each region 
tlev_by_city <- matrix('', nrow = length(num_cities_tlev), ncol = length(cities))
rownames(tlev_by_city) <- names(num_cities_tlev)
colnames(tlev_by_city) <- names(cities)
tlev_by_province <- matrix('', nrow = length(num_provinces_tlev), ncol = length(provinces))
rownames(tlev_by_province) <- names(num_provinces_tlev)
colnames(tlev_by_province) <- names(provinces)
tlev_by_community <- matrix('', nrow = length(num_communities_tlev), ncol = length(comunidades))
rownames(tlev_by_community) <- names(num_communities_tlev)
colnames(tlev_by_community) <- names(comunidades)

for (g in names(num_schools_tlev)) {
  tlev_by_city[ g, num_cities_tlev[[g]] ] <- "+"
  tlev_by_province[ g, num_provinces_tlev[[g]] ] <- "+"
  tlev_by_community[ g, num_communities_tlev[[g]] ] <- "+"
}

write.csv(tlev_by_city, file = sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_by_city.csv", home_dir, tlev, tlev) )
write.csv(tlev_by_province, file = sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_by_province.csv", home_dir, tlev, tlev) )
write.csv(tlev_by_community, file = sprintf("%s/Part_1/figures/ttests/%s/Presence/%s_presence_by_community.csv", home_dir, tlev, tlev) )

## get the number of samples in which each OTU is present, and the mean abundance in those samples in which it is present
num_samples_present <- sapply(rownames(P_A_table), function(x) sum(P_A_table[x, ]))
percent_samples_present <- sapply(rownames(P_A_table), function(x) 100 * sum(P_A_table[x, ]) / ncol(P_A_table) )
mean_abund <- sapply(rownames(P_A_table), function(x) mean(otutab[x, otutab[x,]!=0]) )

# get sample names in which a given OTU is present
which_samples_present <- colnames(P_A_table[,P_A_table["Amaricoccus",]==1])









# ****************************************************************************************************************** #
# Test various differentiation statistics - Kruskal-Wallis ####
# ****************************************************************************************************************** #


#kruskal-wallis to test for any separations within a given survey response
#if one response has clear difference, can group responses binarily, then do wilcoxon
# toremove <- c(1,2,13:34,41,43,45:50,99,101)
#excel col:  (A,B, M-AH,AO,AQ,AS-AX,CU,CW )


tlev <- "Genus"; otutab_rel <- otu_table(SLL1)
# tlev <- "Class"; otutab_rel <- tlev_otus_rel[[tlev]]
# tlev <- "Phylum"; otutab_rel <- tlev_otus_rel[[tlev]]


st <- "all"
# st <- "students"
# st <- "teachers"
# st <- "no-bottled"
# st <- "only-bottled"

t0 <- Sys.time()

st_samps <- sample_names(SLL1)
if (st == "students") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q1.1"]=="Student" ]
} else if (st == "teachers") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q1.1"]=="Teacher" ]
} else if (st == "no-bottled") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q27"]!="Embotellada" ]
} else if (st == "only-bottled") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q27"]=="Embotellada" ]
}


# questions with categories
group_qs <- c("Q2","Q8","Q15","Q15.1","Q18","Q19","Q19.1","Q19.3","Q19.5","Q19.6","Q19.8","Q19.9","Q20.1","Q21","Q27",
              "Q30","Q31","Q32","Q33","Q34","Q35","Q36","Q37","Q38","Q41","Q43","Q44","Q45","Q46","Q47","Q48",
              "Q49","Q50","Q51","Q52","Q53","Q54","Q1.1","Q14.1","Q14.2","Q22.1","Q23.1","Q24.1","Q25.1","Q26.1",
              "Q39.1","Q40.1","Q28.11","Q28.13","Q28.15","Q28.17","Q28.19","Q28.21","Q28.23","Q28.25","Q28.27","Q28.29",
              "Q28.31","Stomatotype","Stomatotype_CORE","Diversity_group_Div.Shannon","Diversity_group_Div.Simpson",
              "Diversity_group_Weighted_Unifrac","Diversity_group_Unweighted_Unifrac","Diversity_group_Faiths.PD",
              "Diversity_group_Species_Richness","BMI_group","BMI_official",group_water_data)#,"Water_hardness_group"
if (st %in% c("students", "teachers")) group_qs <- group_qs[group_qs != "Q1.1"]
group_qs <- as.matrix(sample_data(SLL1)[ st_samps, group_qs])
# group_qs[group_qs=='No Sabe/No Contesta'] <- NA


data.mix <- cbind(t(otutab_rel[ , st_samps]), sample_data(SLL1)[ st_samps, only_cont])
data.mix[data.mix=='No Sabe/No Contesta'] <- NA
data.mix <- as.data.frame(apply(data.mix, 2, factor))


kw.stat <-  matrix(NA, nrow=ncol(data.mix), ncol=ncol(group_qs))
colnames(kw.stat) <- colnames(group_qs)
rownames(kw.stat) <- colnames(data.mix)
kw.p <-  matrix(NA, nrow=ncol(data.mix), ncol=ncol(group_qs))
colnames(kw.p) <- colnames(group_qs)
rownames(kw.p) <- colnames(data.mix)

#store kw$statistic, kw$p.value
for (i in colnames(data.mix)) {
  # print(i)
  for (j in colnames(group_qs)) {
    # print(c(i,j))
    if (length(table(sample_data(SLL1)[st_samps, j])) == 1) {
      kw.stat[i,j] <- 0
      kw.p[i,j] <- 1
    } else if (length(table(sample_data(SLL1)[st_samps, j])) == 2 & min(table(sample_data(SLL1)[st_samps, j])) == 1) {
      kw.stat[i,j] <- 0
      kw.p[i,j] <- 1
    } else {
      kw <- kruskal.test(as.numeric(data.mix[,i]), as.factor(group_qs[,j]), na.rm=T )#na.action='na.exclude'
      kw.stat[i,j] <- kw$statistic
      kw.p[i,j] <- kw$p.value
    }
  }
}



kw.p.adj <- apply(kw.p, 2, p.adjust, method='bonferroni', n=ncol(data.mix))
# for those values that are given as NaN because the given genus is not present in any samples with particular categories of the given question
kw.p.adj[ is.nan(kw.p.adj) ] <- 1
# kw.p.adj <- apply(kw.p.adj, 2, function(x) {-log10(x)}) #only necessary if plotting heatmap of Ps

#at least 1 good p value within samples
# maxes <- apply(kw.p.adj,2,max)
# goodPs <- maxes[maxes > 1.3013] #-log10(0.05) == 1.3013 
# #at least 1 good p within taxa
# tmaxes <- apply(kw.p.adj,1,max)
# tgoodPs <- tmaxes[tmaxes > 1.3013]

#at least 1 good p value within samples
mins <- apply(kw.p.adj,2,min)
goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
#at least 1 good p within taxa
tmins <- apply(kw.p.adj,1,min)
tgoodPs <- tmins[tmins < 0.05]

good.kw.p.adj <- kw.p.adj[names(tgoodPs),names(goodPs)]
good.kw.stat <- kw.stat[names(tgoodPs),names(goodPs)]

# filter those OTUs appearing in fewer than a given number of samples
min_samples <- 20
good_rows <- sapply(rownames(good.kw.p.adj), function(x) ifelse(x %in% rownames(otutab_rel), 
                                                                sum(otutab_rel[x, st_samps]!=0) >= min_samples,
                                                                TRUE))
good.kw.p.adj <- good.kw.p.adj[ good_rows, ]
good.kw.stat  <- good.kw.stat[ good_rows, ]


#clustering order
# test <- kw.p[!rownames(kw.p) %in% 'Streptococcus',]
# ord <- hclust( dist(good.kw.p.adj, method = "euclidean"), method = "ward.D" )$order
# ord2 <- hclust( dist(t(good.kw.p.adj), method = "euclidean"), method = "ward.D" )$order
# kw.p.m <- melt(good.kw.p.adj[ord,ord2])
ord <- hclust( dist(scale(good.kw.stat), method = "euclidean"), method = "ward.D" )$order
ord2 <- hclust( dist(scale(t(good.kw.stat)), method = "euclidean"), method = "ward.D" )$order
kw.p.m <- melt(good.kw.stat[ord,ord2])

pluses <- matrix('', nrow=length(rownames(good.kw.p.adj)), ncol=length(colnames(good.kw.p.adj)))
colnames(pluses) <- colnames(good.kw.p.adj)
rownames(pluses) <- rownames(good.kw.p.adj)
# pluses[good.kw.p.adj > 1.3013] <- '+'
pluses[good.kw.p.adj < 0.05] <- '+'
# kw.p.m <- cbind(melt(good.kw.p.adj), melt(pluses))
kw.p.m <- cbind(melt(good.kw.stat), melt(pluses))
kw.p.m <- kw.p.m[,c(1,2,3,6)]
colnames(kw.p.m)[4] <- 'signif'



ggplot(data = kw.p.m, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color="white") + 
  geom_text(aes(label=signif), size = 5) +
  scale_fill_gradient2(low="royalblue", high="red", mid="white", midpoint=max(good.kw.stat)/2, name="KW stat") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        plot.title = element_text(hjust=0.5)) +
  ggtitle(sprintf('KW statistic for %s', tlev)) +
  xlab('') + ylab('')
# ******************** #



# ******** #
# A table with the KW statistic only for those values that are statistically significant
#    (have a '+' in the pluses matrix)
# ******** #
kw.signif <- matrix(NA, nrow=nrow(good.kw.p.adj), ncol=ncol(good.kw.p.adj))
colnames(kw.signif) <- colnames(good.kw.p.adj)
rownames(kw.signif) <- rownames(good.kw.p.adj)

# verify significant differences with anova to filter out signif differences in "No sabe" response groups
kw_adf <- cbind(t(otutab_rel[ , st_samps]), sample_data(SLL1)[ st_samps, ])
kw_adf <- as.data.frame(as.matrix(kw_adf))
for (o in c(rownames(otutab_rel[ , st_samps]), only_cont)) {
  if (! o %in% c("Stomatotype","Stomatotype_CORE")) { # in this case treat Stomatotype as a categorical variable
    kw_adf[,o] <- as.numeric(as.matrix(kw_adf[,o]))
  }
}

for (i in rownames(kw.signif)) {
  for (j in colnames(kw.signif)) {
    # run pairwise anova tests to check where the significant differences lie
    res.aov <- aov(formula = as.matrix(kw_adf[, i]) ~ as.matrix(kw_adf[, j]), data = kw_adf)
    thsd <- TukeyHSD(res.aov)
    
    group_diffs <- rownames(thsd$`as.matrix(kw_adf[, j])`)
    sig_groups <- group_diffs[ thsd$`as.matrix(kw_adf[, j])`[group_diffs, "p adj"] < 0.05 ]
    
    # check if signif in Kruskal-Wallis test and in at least one anova test without a "No Sabe" response
    if (pluses[i,j] == '+' & sum(!grepl("No Sabe/No Contesta", sig_groups)) > 0) {
      kw.signif[i,j] <- good.kw.p.adj[i,j]
    } else {
      # check if anova wasnt signif because one group had only 0s as cont value
      j.1 <- names(table(kw_adf[, j])[1])
      j.2 <- names(table(kw_adf[, j])[2])
      j.1.s <- sum(kw_adf[ kw_adf[, j]==j.1, i ])
      j.2.s <- sum(kw_adf[ kw_adf[, j]==j.2, i ])
      
      if ( pluses[i,j] == '+' & length(table(kw_adf[, j]))==2 & 0 %in% c(j.1.s, j.2.s) ) {
        kw.signif[i,j] <- good.kw.p.adj[i,j]
      } else {
        kw.signif[i,j] <- ''
      }
    }
  }
}

keep_cols <- colnames(kw.signif)[ sapply(colnames(kw.signif), function(x) sum( unique(kw.signif[ , x]) != "") > 0) ]
keep_rows <- rownames(kw.signif)[ sapply(rownames(kw.signif), function(x) sum( unique(kw.signif[ x, ]) != "") > 0) ]
good.kw.signif <- kw.signif[ keep_rows, keep_cols]

print(Sys.time() - t0)

# write.csv(kw.signif, sprintf('%s/Part_1/figures/Kruskal-Wallis/%s/%s/%s_%s_%s-minSamps_signif_kw_stat.csv', home_dir, tlev, st, tlev, st, min_samples))
write.csv(good.kw.signif, sprintf('%s/Part_1/figures/Kruskal-Wallis/%s/%s/PAPER_%s_%s_anova_%s-minSamps_signif_kw_stat.csv', home_dir, tlev, st, tlev, st, min_samples))
if (tlev=="Genus") {
  write.csv(good.kw.signif, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_3.csv', home_dir))
}
# ***************************************************************** #


# some specific t-test checks ####
t.test( as.numeric(as.matrix(tlev_otus_rel[["Phylum"]]["Chloroflexi", rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Student",])])),
        as.numeric(as.matrix(tlev_otus_rel[["Phylum"]]["Chloroflexi", rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Teacher",])])))

t.test( as.numeric(as.matrix(otu_table(SLL1)["Burkholderia", rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Teacher",])])),
        as.numeric(as.matrix(otu_table(SLL1)["Burkholderia", rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Student",])])))

t.test( as.numeric(as.matrix(sample_data(SLL1)[rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q53"]==1,]), "Q20.1"])),
        as.numeric(as.matrix(sample_data(SLL1)[rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q53"]==0,]), "Q20.1"])))



# ******************** #
# boxplots of particular group_q vs cont data/otu abundance ####
group_col <- "Q1.1"
cont_col <- "Alloscardovia"
kw.box <- as.data.frame( cbind(as.matrix(data.mix[, cont_col]), group_qs[, group_col]) )
colnames(kw.box) <- c("cont","group")
# kw.box <- kw.box[ kw.box[ , "group"] != "No Sabe/No Contesta", ]

ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
                   fill=reorder(group,-as.numeric(as.character(cont)),median))) +
  geom_boxplot() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(group_col) + ylab(cont_col) + scale_fill_hue(name=group_col) #+ ylim(10,50)
# ******************** #



# ************************************* #
# To perform one-way anova ####

st <- "all"
# st <- "students"
# st <- "teachers"
# st <- "no-bottled"
# st <- "only-bottled"

st_samps <- sample_names(SLL1)
if (st == "students") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q1.1"]=="Student" ]
} else if (st == "teachers") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q1.1"]=="Teacher" ]
} else if (st == "no-bottled") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q27"]!="Embotellada" ]
} else if (st == "only-bottled") {
  st_samps <- st_samps[ sample_data(SLL1)[,"Q27"]=="Embotellada" ]
}

adf <- cbind(t(otutab_rel[ , st_samps]), sample_data(SLL1)[ st_samps, ])
adf <- as.data.frame(as.matrix(adf))
for (o in c(rownames(otutab_rel[ , st_samps]), only_cont)) {
  if (! o %in% c("Stomatotype","Stomatotype_CORE")) { # in this case treat Stomatotype as a categorical variable
    adf[,o] <- as.numeric(as.matrix(adf[,o]))
  }
}

res.aov <- aov(formula = as.matrix(adf[, cont_col]) ~ as.matrix(adf[, group_col]), data = adf)
# summary(res.aov)
TukeyHSD(res.aov)




library(dplyr)
library(rlang)
# look at means of cont_col for each group in group_col
group_by_(adf, group_col) %>%
  summarise(
    count = n(),
    mean = mean(UQ(sym(cont_col)), na.rm = TRUE),
    sd = sd(UQ(sym(cont_col)), na.rm = TRUE)
  )
# ************************************* #



table(kw.box[,2])
# run a t-test for those binary group questions (e.g. yes/no)
t.test( as.numeric(as.character(kw.box[kw.box[,2]=="1",1])), 
        as.numeric(as.character(kw.box[kw.box[,2]=="0",1])))
# run a t-test for those binary group questions (e.g. yes/no)
t.test( as.numeric(as.matrix(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Teacher", "Q20.1"])), 
        as.numeric(as.matrix(sample_data(SLL1)[sample_data(SLL1)[,"Q1.1"]=="Student", "Q20.1"])))
# how many samples have a given group for the given question
dim(kw.box[kw.box[,"group"]=="1",])
# get mean value of OTU for samples in indicated group
mean(otu_table(SLL1)["Burkholderia",sample_data(SLL1)[,"Q1.1"]=="Teacher"])
# get mean value of continuous question for samples in indicated group
mean(as.numeric(as.matrix(sample_data(SLL1))[sample_data(SLL1)[,"Q1.1"]=="Teacher", "Q28.1"]), na.rm = T)
mean(as.numeric(as.matrix(otu_table(SLL1))["Nevskia", sample_data(SLL1)[,"Q15.1"]=="Pueblo/Campo"]))
# get counts of responses for indicated group
table(as.numeric(as.matrix(sample_data(SLL1))[sample_data(SLL1)[,"Q1.1"]=="Student", "Q28.1"]))

# with Q19.9:
# Pyramidobacter -- pajaro, peces
# Rheinheimera -- gallina
# Erysipelotrichaceae_incertae_sedis -- Erizo
# Lactovum -- Agaporni, periquito
# Syntrophococcus -- pajaro, tortuga
# Negativicoccus -- pajaro, peces
# Bergeriella -- Agaporni, periquito
# Zoogloea -- Erizo
# Desulfovibrio -- Pajaro, peces
# Mycobacterium -- Erizo
# *** None of these is actually interesting because in all of these cases, there were only 2 with the indicated category of animal, 
#     and in each case, only 1 of the 2 samples had a non-0 abundance of the indicated genus

People with girl/boyfriend -- higher Num_OTUs
Cold alcoholic drinks -- higher Num_OTUs
No sweets -- more Anaerovorax
Stomatotype 1 -- higher Div.Simpson

Teachers -- higher unweighted unifrac
# ****************************************************************************************************************** #

# check how many samples contain genera of a given hgher level OTU: ####
sort(sapply(rownames(tax_table(SLL1)[tax_table(SLL1)[,"Family"]=="Rhodocyclaceae",]), function(x) num_samples_present[[x]]))
# check mean_abund of genera of a given hgher level OTU:
sort(sapply(rownames(tax_table(SLL1)[tax_table(SLL1)[,"Order"]=="Rhodospirillales",]), function(x) mean_abund[[x]]))







# ****************************************************************************************************************** #
# *********************************************************** #
# Chi-squared test ####
# *********************************************************** #

library(graphics)
library(vcd)
library(corrplot)

tlev <- "Genus"; otutab <- otu_table(SLL)
# tlev <- "Phylum"; otutab <- tlev_otus[[tlev]]

P_A_table <- apply( otutab, 2, function(x) ifelse(x==0, 0, 1) )
# take only those columns for OTUs which are not present in 100% of samples
diff_pres <- sapply( rownames(P_A_table), function(x) dim(table(P_A_table[ x, ])) != 1)
P_A_table <- P_A_table[ diff_pres, ]

# questions with categories
group_Qs <- c("Q2","Q8","Q15","Q15.1","Q18","Q19","Q19.1","Q19.3","Q19.5","Q19.6","Q19.8","Q19.9","Q20.1","Q21","Q27",
              "Q30","Q31","Q32","Q33","Q34","Q35","Q36","Q37","Q38","Q41","Q43","Q45","Q46","Q48",
              "Q49","Q50","Q51","Q52","Q53","Q54","Q1.1","Q14.1","Q14.2","Q22.1","Q23.1","Q24.1","Q25.1","Q26.1",
              "Q39.1","Q40.1","Q28.11","Q28.13","Q28.15","Q28.17","Q28.19","Q28.21","Q28.23","Q28.25","Q28.27","Q28.29",
              "Q28.31","Stomatotype","Stomatotype_CORE","Diversity_group_Div.Shannon","Diversity_group_Div.Simpson",
              "Diversity_group_Weighted_Unifrac","Diversity_group_Unweighted_Unifrac","Diversity_group_Faiths.PD",
              "Diversity_group_Species_Richness","BMI_group","BMI_official",group_water_data)#,"Water_hardness_group"
group_Qs <- as.matrix(sample_data(SLL1)[, group_Qs])
group_Qs <- cbind( group_Qs, t(P_A_table))

chi.matrix <- matrix(NA, nrow = ncol(group_Qs), ncol = ncol(group_Qs))
rownames(chi.matrix) <- colnames(chi.matrix) <- colnames(group_Qs)

for (i in colnames(group_Qs)) {
  for (j in colnames(group_Qs)) {
    # conting_tab <- table(as.character(as.matrix(sample_data(SLL)[, i])), 
    #                      as.character(as.matrix(sample_data(SLL)[, j])))
    conting_tab <- table(group_Qs[ , i], group_Qs[ , j])
    conting_tab <- conting_tab[rownames(conting_tab) != "No Sabe/No Contesta", 
                               colnames(conting_tab) != "No Sabe/No Contesta"]
    # conting_tab <- conting_tab[rownames(conting_tab)!="N.A.", colnames(conting_tab)!="N.A."]
    chi <- chisq.test(conting_tab)
    
    # since shouldnt use Chi-square when any expected value is below 5, use Fishers test here instead
    # if (length( chi$expected[chi$expected<5]) > 0) {
    #   print(c(i,j, length( chi$expected[chi$expected<5])))
    #   chi <- fisher.test(conting_tab, simulate.p.value = T)
    # }
    
    chi.matrix[i,j] <- chi$p.value
  }
}

# for those values that are given as NaN because at least one of the categories had all 0s 
chi.matrix[ is.nan(chi.matrix) ] <- 1
# adjust for multiple testing
chi.matrix.adj <- apply(chi.matrix, 2, p.adjust, method='bonferroni', n=ncol(group_Qs))
diag(chi.matrix.adj) <- 1 #because p for the diagonal will always be 0, and want to remove rows/cols without signif p values

#at least 1 good p value within samples
mins <- apply(chi.matrix.adj,2,min)
goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013 
#at least 1 good p within taxa
tmins <- apply(chi.matrix.adj,1,min)
tgoodPs <- tmins[tmins < 0.05]

good.chi.matrix.adj <- chi.matrix.adj[names(tgoodPs),names(goodPs)]
only.good.chi.matrix.adj <- good.chi.matrix.adj
only.good.chi.matrix.adj[only.good.chi.matrix.adj >= 0.05] <- ''

pluses <- matrix('', nrow=length(rownames(chi.matrix.adj)), ncol=length(colnames(chi.matrix.adj)))
colnames(pluses) <- colnames(chi.matrix.adj)
rownames(pluses) <- rownames(chi.matrix.adj)
pluses[chi.matrix.adj < 0.05] <- '+'

write.csv(only.good.chi.matrix.adj, sprintf("%s/Part_1/figures/Chi-squared/%s/signif_chi-squared.csv", home_dir, tlev))
if (tlev=="Genus") {
  write.csv(only.good.chi.matrix.adj, sprintf('%s/Part_1/SLL1_paper/supp_material/Supp_table_4.csv', home_dir))
}




rows <- "Q53"
cols <- "Q20.1"
# contingency <- table(as.character(as.matrix(sample_data(SLL)[ , rows ])), 
#                      as.character(as.matrix(sample_data(SLL)[ , cols ])),
#                      dnn = c(rows,cols))
contingency <- table(group_Qs[ , rows ], group_Qs[ , cols ], dnn = c(rows,cols) )
contingency <- contingency[rownames(contingency) != "No Sabe/No Contesta", 
                           colnames(contingency) != "No Sabe/No Contesta"]

# order diversity groups logically
if (startsWith(rows, "Diversity_group")) {contingency <- contingency[ c("Low","Average","High"), ]
} else if (startsWith(cols, "Diversity_group")) {contingency <- contingency[ , c("Low","Average","High") ]}
# order hardness groups logically
if (rows == "Hardness_category") {contingency <- contingency[ c("Muy blanda","Blanda","Dura","Muy dura","Extremadamente dura"), ]
} else if (cols == "Hardness_category") {contingency <- contingency[ , c("Muy blanda","Blanda","Dura","Muy dura","Extremadamente dura") ]}
# name location groups briefly
if (rows == "Q15") {rownames(contingency) <- c("Afueras/Campo","Urbano","Pueblo")
} else if (cols == "Q15") {colnames(contingency) <- c("Afueras/Campo","Urbano","Pueblo") }
# name ethnic groups briefly
if (rows == "Q8") {rownames(contingency) <- c("Africano","Arabe","Asiatico","Caucasico","Gitano","Nativo Americano")
} else if (cols == "Q8") {colnames(contingency) <- c("Africano","Arabe","Asiatico","Caucasico","Gitano","Nativo Americano") }
# order mineralization groups logically
if (rows == "Mineralization") {contingency <- contingency <- contingency[ c("Muy debil","Debil","Media","Fuerte","Oligometalica"), ]
} else if (cols == "Mineralization") {contingency <- contingency[ , c("Muy debil","Debil","Media","Fuerte","Oligometalica") ] }

if (nrow(contingency) > ncol(contingency)) contingency <- t(contingency)

chi <- chisq.test(contingency)

assoc(contingency, shade = T, main = sprintf("Association plot"))
# ****************************************************************************************************************** #








# # ****************************************************************************************************************** #
# Test various differentiation statistics - Wilcoxon ####
# # ****************************************************************************************************************** #
# 
# q.wil <- "Q5"
# o.wil <- "Streptococcus"
# cats <- levels(data.mix[ , q.wil])
# 
# wil.stat <-  matrix(NA, nrow=1, ncol=length(cats))
# colnames(wil.stat) <- cats
# rownames(wil.stat) <- "Streptococcus"
# wil.p <-  matrix(NA, nrow=1, ncol=length(cats))
# colnames(wil.p) <- cats
# rownames(wil.p) <- "Streptococcus"
# 
# for (i in cats) {
#   group <- data.mix[,q.wil]==i
#   wil.res <- wilcox.test(x = as.numeric(group), y = as.numeric(data.mix[,'Streptococcus']), na.rm =T)
#   wil.stat[1,i] <- wil.res$statistic
#   wil.p[1,i] <- wil.res$p.value
# }
# 
# #check size of groups with good scores
# new.wil <- wil.stat[, wil.stat > 1]
# # ***** FIX THIS *****






# ****************************************************************************************************************** #
# Boxplots for questions with significant kruskal against Diversity ####
# ****************************************************************************************************************** #
toremove <- c(1,2,13:34,41,43,45:50,99,101)
#excel col:  (A,B, M-AH,AO,AQ,AS-AX,CU,CW )
data.mix.box <- cbind(t(otu_table(gs)), sample_data(gs)[,-toremove])
otus.mix.box <- rownames(otu_table(gs))
questions.mix.box <- colnames(sample_data(gs))[-toremove]


qs <- vector('character')
j=1
measure <- 'Div.Shannon'
# measure <- 'Streptococcus'
for (i in questions.mix.box) {
  k <- kruskal.test(data.mix.box[,measure], data.mix.box[,i], na.rm=T )$p.value
  if (k<0.05) {
    qs[j] <- i
    print(c(k, i))
    j <- j+1
  }
}

# If measure is a diversity measure (not an OTU), must take only the indicated diversity measure
# and ignore the others because melt will group them all together under the "variable" column.
# If measure is an OTU, will simply ignore all the diversity measures.
divs_n_num <- c("Div.Observed", "Div.Chao1", "Div.ACE", "Div.Shannon", "Div.Simpson", "Div.InvSimpson", "Div.Fisher", "Num_OTUs")
other_divs <- divs_n_num[divs_n_num != measure]
qs <- qs[-which(qs %in% divs_n_num)]


#boxplots for questions that vary with Diversity significantly
q <- 'Qsamples.1'
data.m <- melt(data.mix.box[,c(measure,qs)])


# For boxes of the specific question
dat <- data.m[,c("value", q)]
dat <- dat[dat[,2]!="No Sabe/No Contesta",]
colnames(dat)[2] <- 'groups'
# counts <- count(dat[2])
counts <- function(x) {
  return( data.frame(y=median(x)+0.08, label=length(x)) )
}

ggplot(dat, aes(x=reorder(groups,-value,FUN=median), y=value, fill=reorder(groups,-value,FUN=median))) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
  ggtitle(sprintf('%s per sample for %s', measure, q)) +
  xlab(q) + ylab(measure) + scale_fill_hue(name='Levels') +#guides(fill=FALSE)+
  stat_summary(fun.data = counts, geom='text')


# ****************************************************************************************************************** #




# ******************* #
#wilcoxon to test for difference between binary respsonses

#list of questions for which kruskal test gave signif difference and can by separated to two categories
#Can create new cols in excel to separate (e.g. Q5,Q6,Q7 -> Barcelona or not, Q15 -> Urbano or not)

#Include (from colnames(good.kw.stat)): 2, 5,6,7,(14.1), 15, 19,19.1, 20.1,21,30,38,48,49,51,53
#Change these from ranges to binary (may have to play with the divider a bit): 22.1, 25.1, 39.1, 28.11, 28.12, 28.16,28.18,28.21, 24,29
#Think more about 27, 42,43, 4(greater/less than a given weight?), 8, 14.1

#Some of the signif from Kruskal will be ignored here because cannot separate into two cats (4,17,Socioeconomic,PH,16.1)...could do same as with weight?
#Some have been converted into binary responses already but were not signif with Kruskal (which?)

binary_qs <- c(3,)
#excel col:   (C,)
#...(* under those that are also ranges)


# ******************* #















# ***************** #
# Check if diet varies between ciudad and campo samples ******* ####
# ****Can use this for any question vs question comparison involving categorical values

# onlycont <- c(3,4,5,16:21,24:31,33,36:44,46,47,49,51,52,59:77,79:82,85,88:98,100,102,106:144,149:155,  156)
#  excel col:  (C,D,E, P:V,  X:AE,AG,AJ:AR,AT,AU,AW,AY,AZ,BG:BY,CA:CD,CG,CJ:CT,CV,  CX, DB:EN ,div-stats,#otus)
for.kw.qs <- c(59:77, 92:98, 107:143)
kw.qs <- colnames(sample_data(gs))[for.kw.qs]
kw.stat.questions <- matrix(NA, nrow=1, ncol=length(kw.qs))
colnames(kw.stat.questions) <- kw.qs
kw.p.questions <- matrix(NA, nrow=1, ncol=length(kw.qs))
colnames(kw.p.questions) <- kw.qs

for (i in kw.qs) {
  kw <- kruskal.test(data.mix[,i], data.mix[,"Q15"], na.rm=T)
  kw.stat.questions[1,i] <- kw$statistic
  kw.p.questions[1,i] <- kw$p.value
}

kw.p.questions.adj <- apply(kw.p.questions, 2, p.adjust, method='bonferroni', n=length(kw.qs))

# significant diffs amongst the categories:
kw.p.questions.adj[kw.p.questions.adj < 0.05]


ciudad_samples = rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q15"]=="Casco Urbano, Alto Nivel De Habitantes Por Metro Cuadrado"])
pueblo_samples = rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q15"]=="Pueblo, Menor Nivel Habitantes Por Metro Cuadrado"])
afueras_samples = rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q15"]=="Afueras/Campo"])
pue_afue_samples = c(rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q15"]=="Afueras/Campo"]),
                     rownames(sample_data(SLL1)[sample_data(SLL1)[,"Q15"]=="Pueblo, Menor Nivel Habitantes Por Metro Cuadrado"]))

mean(as.numeric(as.matrix(sample_data(SLL1)[ciudad_samples, "Q22.2"] )), na.rm=T)
mean(as.numeric(as.matrix(sample_data(SLL1)[pueblo_samples, "Q22.2"] )), na.rm=T)
mean(as.numeric(as.matrix(sample_data(SLL1)[afueras_samples, "Q22.2"] )), na.rm=T)
mean(as.numeric(as.matrix(sample_data(SLL1)[pue_afue_samples, "Q22.2"] )), na.rm=T)

# ***************** #









# ***************** #
# Check for differences between samples of the different cities or between a city and the rest of the samples: ####
# cities <- list( Barcelona = c(26,27,28,34,41),
#                 Bilbao = c(11,17,32),
#                 Madrid = c(5,20,21,33,40),
#                 Malaga = c(16,24,31),
#                 Mallorca = c(19,22,29,39),
#                 Murcia = c(3,12,25,35),
#                 Sevilla = c(6,18,30),
#                 Valencia = c(2,4,7,8,10,15),
#                 Vigo = c(9,36,37,38),
#                 Zaragoza = c(1,13,14,23) )

cities <- list( Barcelona = c(27,28,41), Calella = c(34), Roda_de_Ter = c(26),
                Santurtzi = c(11,17), Altsasu = c(32),
                Madrid = c(5,20,21), Soto_del_Real = c(40), Cercedilla = c(33),
                Malaga = c(16,24), Torrox = c(31),
                Palma = c(29,39), Santanyi = c(19), Port_d_Alcudia = c(22),
                Murcia = c(3,12), Bullas = c(25), La_Paca = c(35),
                Sevilla = c(6,18), El_Saucejo = c(30),
                Valencia = c(4,15), Nules = c(8), L_Alcora = c(2,7), Tavernes_de_Valldigna = c(10),
                Vigo = c(9,36), Bande = c(37), Villarcayo = c(38),
                Zaragoza = c(1,14), Tauste = c(13,23) )

city1 <- 'Tauste'
city2 <- 'Madrid'
SLL1.city1 <- prune_samples( rownames(sample_data(SLL1)[ sample_data(SLL1)$Qsamples.1 %in% cities[city1][[1]], ]), SLL1 )
SLL1.city2 <- prune_samples( rownames(sample_data(SLL1)[ sample_data(SLL1)$Qsamples.1 %in% cities[city2][[1]], ]), SLL1 )
SLL1.city.others <- prune_samples( rownames(sample_data(SLL1)[ !sample_data(SLL1)$Qsamples.1 %in% cities[city1][[1]], ]), SLL1 )


top.tax <- names(sort(taxa_sums(SLL1), decreasing = T))
for (g in top.tax) {
  t.1 <- t.test( as.numeric(as.matrix( otu_table(SLL1.city1)[ g, ]) ),
                 as.numeric(as.matrix( otu_table(SLL1.city.others)[ g, ]) ),
                 alternative = "greater")$p.value
  if (is.nan(t.1)==F & t.1*length(top.tax) < 0.05) {
    print(c(g,t.1))
  }
}
# ***************** #


# Check for differences between samples of the different water types consumed at home: ####
SLL1.filtrada <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q27 == "Del Grifo (Filtrada)" ]), SLL1)
SLL1.nofiltrada <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q27 == "Del Grifo (No Filtrada)" ]), SLL1)
SLL1.embotellada <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q27 == "Embotellada" ]), SLL1)
SLL1.notratada <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q27 == "No Tratada (Fuente, Pozo O Rio)" ]), SLL1)


top.tax <- names(sort(taxa_sums(SLL1), decreasing = T))
for (g in top.tax) {
  t.1 <- t.test( as.numeric(as.matrix( otu_table(SLL1.notratada)[ g, ]) ),
                 as.numeric(as.matrix( otu_table(SLL1.embotellada)[ g, ]) ),
                 alternative = "greater")$p.value
  if (is.nan(t.1)==F & t.1*length(top.tax) < 0.05) {
    print(c(g,t.1))
  }
}
# ***************** #


# Check for differences between samples of the different frequencies of dulces consumed: ####
SLL1.n <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "Nunca O Menos De 1 Vez Al Mes" ]), SLL1)
SLL1.1m <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "1-3 Veces Al Mes" ]), SLL1)
SLL1.1s <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "1-2 Veces A La Semana" ]), SLL1)
SLL1.3s <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "3-4 Veces A La Semana" ]), SLL1)
SLL1.5s <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "5-6 Veces A La Semana" ]), SLL1)
SLL1.1d <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "1 Vez Al Dia" ]), SLL1)
SLL1.2d <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "2 Veces Al Dia" ]), SLL1)
SLL1.m2d <- prune_samples(rownames(sample_data(SLL1)[ sample_data(SLL1)$Q24 == "Mas De 2 Veces Al Dia" ]), SLL1)


top.tax <- names(sort(taxa_sums(SLL1), decreasing = T))
for (g in top.tax) {
  t.1 <- t.test( as.numeric(as.matrix( otu_table(SLL1.1m)[ g, ]) ),
                 as.numeric(as.matrix( otu_table(SLL1.m2d)[ g, ]) ),
                 alternative = "greater")$p.value
  if (is.nan(t.1)==F & t.1*length(top.tax) < 0.05) {
    print(c(g,t.1))
  }
}
# ***************** #













# ****************************************************************************************************************** #
# Diversity within samples ####
# ****************************************************************************************************************** #
#Prune only those taxa with a count of 0, cannot normalize yet...plot_richness accepts only integers
if (region_type == "all") {
  SLL.diversity = prune_taxa(taxa_sums(SLL)>0, SLL)
} else if (region_type == "region") {
  SLL.region <- prune_samples( region_samples, SLL )
  SLL.diversity = prune_taxa(taxa_sums(SLL)>0, SLL.region)
}


alpha <- 'Simpson'
alpha <- 'Shannon'

plot_by <- 'sample'
# plot_by <- 'question'
# plot_by <- 'question_merge'


if (plot_by=='sample') {
  plot_richness(SLL.diversity, measures = alpha, 
                x=reorder(sample_names(SLL.diversity), estimate_richness(SLL.diversity)[,alpha])) +
    theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
    ggtitle(sprintf('Diversity index per sample tax%s', level)) +
    xlab('Samples') + ylab(sprintf('%s Alpha diversity', alpha))
  ## try plotting abundances ordered by richness to see trends here.....

} else if (plot_by=='question') {
  plot_richness(SLL.diversity, measures = alpha, x='Q29.4') +
    ggtitle(sprintf('%s Diversity index by tax%s', alpha,level))
  
} else if (plot_by=='question_merge') {
  #will need to make new column for whether sample is urbano or other
  #using function getVariable() %in% 'urbano...' ==> see 'plot_richness-examples from phyloseq
  
  SLLmerge <- merge_samples(SLL.diversity, 'Q8')
  sample_data(SLLmerge)[ ,'Q8'] <- factor(sample_names(SLLmerge))
}









# ****************************************************************************************************************** #
#here make plot similar to the abundances but for those samples with diversity below 0.7 or over 0.9
# ****************************************************************************************************************** #

# div.estimate <- "Div.Observed"
# div.estimate <- "Div.Chao1"
# div.estimate <- "Div.ACE"
# div.estimate <- "Div.Shannon"
div.estimate <- "Div.Simpson"
# div.estimate <- "Div.InvSimpson"
# div.estimate <- "Div.Fisher"

highlow <- 'low'
# highlow <- 'high'


if (highlow=='low') {
  
  if (level==3) {
    div.samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[,div.estimate] < 0.75 ])
  } else if (level==6) {
    div.samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[,div.estimate] < 0.85 ])
  }

} else if (highlow=='high') {
  
  if (level==3) {
    div.samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[,div.estimate] > 0.85 ])
  } else if (level==6) {
    div.samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[,div.estimate] > 0.9 ])
  }
}
SLL.div <- prune_samples(div.samples, SLL1)

top10 <- names( sort(rowSums(otu_table(SLL.div)), T)[1:10] )
# SLL.div <- prune_taxa(top10, SLL.div)
otu_table(SLL.div) <- otu_table(SLL.div)[top10,]


# ****************************************************************************************************************** #
#mean and sd of given taxa per sample
tax.mps.div <- rowMeans(otu_table(SLL.div))
tax.sd.div <- rowSds(otu_table(SLL.div))

#percentage of samples in which given taxa appears
tax.freqs.div <- round(rowSums(otu_table(SLL.div) != 0) / length(colnames(otu_table(SLL.div)))*100, digits=2)
tax.freqs.div <- paste(as.character(tax.freqs.div),'%','')

#create frame that can be used for bar plots
tax.frame.div <- as.data.frame(tax.mps.div)
tax.frame.div[,2] <- tax.sd.div
tax.frame.div[,3] <- tax.freqs.div
tax.frame.div[,4] <- rownames(otu_table(SLL.div))
colnames(tax.frame.div) <- c("mps","sd","freqs",taxlevel)


num <- length(div.samples)
tot <- length(rownames(sample_data(SLL1)))
ggplot(tax.frame.div, aes(x=reorder(tax.frame.div[,taxlevel], -mps), y=mps, 
                      fill=reorder(tax.frame.div[,taxlevel], -mps))) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mps-sd, ymax=mps+sd), width=0.2) +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
  ggtitle(sprintf('Relative abundance for common taxa in %s diversity samples (%d/%d)',highlow,num,tot)) +
  xlab(taxlevel) + ylab('Mean normalized abundance per sample (as %)') + scale_fill_hue(name=taxlevel) +
  geom_text(aes(x=reorder(tax.frame.div[,taxlevel], -mps), label=freqs), vjust=-1)


# ****************************************************************************************************************** #
#stacked bars of each of the top10 for each sample
divranks <- matrix(0, nrow=length(rownames(otu_table(SLL.div))), ncol=length(colnames(otu_table(SLL.div))))
colnames(divranks) <- colnames(otu_table(SLL.div))
rownames(divranks) <- rownames(otu_table(SLL.div))
tax_ranks <- function(x) {
  if (taxlevel=='tax3') {
    rank( as.data.frame(otu_table(SLL.div)['Bacilli', ]) )
  } else if (taxlevel=='tax6') {
    rank( as.data.frame(otu_table(SLL.div)['Streptococcus', ]) )
  }
}
divranks <-apply(otu_table(SLL.div),1,tax_ranks)
divranks <- t(divranks)

otu.m <- cbind( melt(otu_table(SLL.div)), melt(divranks) )
otu.m <- otu.m[,c(1,2,3,6)]
colnames(otu.m)[4] <- 'ranks'

# #all top 10 stacked, ordered by values for given taxa
# ggplot(otu.m, aes(x=reorder(Var2,-value), y=value, fill=Var1)) +
#   geom_bar(stat='identity') +
#   theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
#   ggtitle(sprintf('Total normalized abundances per sample for top 10 taxa (%s diversity)',highlow)) +
#   xlab('Samples') + ylab('Abundances') + scale_fill_hue(name=taxlevel)


#all top 10 in separate plots, all ordered by values of most abundant taxa
ggplot(otu.m, aes(x=reorder(Var2,-ranks), y=value, fill=Var1)) +
  geom_bar(stat='identity') +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
  facet_wrap(~Var1) +
  ggtitle(sprintf('Total normalized abundances per sample for top 10 taxa at level %s (%s diversity)',level,highlow)) +
  xlab('Samples') + ylab('Abundances') + guides(fill=F) 










# ****************************************************************************************************************** #
# Euclidean distances bt teachers and their students vs distance between all teachers ####
# ****************************************************************************************************************** #
teachers <- data.frame( sample_data(SLL1)[ sample_data(SLL1)[,'Q1.1']=='Teacher', c('Q1.1','Qsamples.1')] )
students <- data.frame( sample_data(SLL1)[ sample_data(SLL1)[,'Q1.1']=='Student', c('Q1.1','Qsamples.1')] )

# ***************** #
#teachers vs students at the same school
tvs.euc <- matrix(0, nrow=1, ncol=22) #for the 22 teacher samples that have otu data available
tvsind <- 1
for (i in 1:41) {
  if (as.character(i) %in% teachers[,'Qsamples.1']) {
    stud_samps <- rownames(students[students[,'Qsamples.1']==i, ])
    teach_samps <- rownames(teachers[teachers[,'Qsamples.1']==i, ])

    teach_otus <- data.frame(otu_table(SLL1)[,teach_samps])
    stud_otus  <- data.frame(otu_table(SLL1)[,stud_samps])
    
    
    for (ts in 1:ncol(teach_otus)) {
      d_total <- 0
      for (j in 1:ncol(stud_otus)) {
        temp <- rbind(t(teach_otus[,ts]), t(stud_otus[,j])) #distance caluclated between the rows of the object
        d <- dist(temp)
        # print(dist(temp))
        d_total <- d_total + d
      }
      d_avg <- d_total/ncol(stud_otus)
      print(c(d_avg, d_total, ncol(stud_otus), ncol(teach_otus)))
      tvs.euc[tvsind] <- d_avg
      tvsind <- tvsind + 1
    }
    # d_avg <- d_total/(ncol(stud_otus) * ncol(teach_otus))
    # print(c(d_avg, d_total, ncol(stud_otus), ncol(teach_otus)))
    
#     tvs.euc[tvsind] <- d_avg
#     tvsind <- tvsind + 1
  }
}
avg_tvs <- mean(tvs.euc)

# ***************** #
#teachers vs all other teachers
teach_only <- data.frame(otu_table(SLL1)[,rownames(teachers)])

tvt.euc <- matrix(0, nrow=1, ncol=22) # 22 teachers available when excluding schools 1-5
tvtind <- 1
for (i in 1:ncol(teach_only)) {
  d_total <- 0
  for (j in 1:ncol(teach_only)) {
    if (i!=j) {
      temp <- rbind(t(teach_only[,i]), t(teach_only[,j])) #distance caluclated between the rows of the object
      d <- dist(temp)
      # print(dist(temp))
      d_total <- d_total + d
    }
  }
  d_avg <- d_total/(ncol(teach_only)-1)
  print(c(d_avg, d_total, ncol(teach_only)))
  tvt.euc[tvtind] <- d_avg
  tvtind <- tvtind + 1
}

avg_tvt <- mean(tvt.euc)


# ********************** #
#students vs students of same school
stud_within <- data.frame(otu_table(SLL1)[,rownames(students)])

svws.euc <- matrix(0, nrow=1, ncol=36) #36 bc first 5 schools removed
svwsind <- 1
for (i in 1:41) {
  if (as.character(i) %in% students[,'Qsamples.1']) {
    # stud_samps_w <- rownames(stud_within[stud_within[,'Qsamples.1']==i, ])
    stud_samps_w <- rownames(students[students[,'Qsamples.1']==i, ])
    
    stud_otus_w  <- data.frame(otu_table(SLL1)[,stud_samps_w])
    
    d_total_school <- 0
    for (x in 1:ncol(stud_otus_w)) {
      d_total_student <- 0
      for (y in 1:ncol(stud_otus_w)) {
        if (x!=y) {
          temp <- rbind(t(stud_otus_w[,x]), t(stud_otus_w[,y]))
          d <- dist(temp)
          d_total_student <- d_total_student + d
        }
      }
      d_avg_stud <- d_total_student/(ncol(stud_otus_w)-1) # -1 bc not including student vs self
      print(c(d_avg_stud, d_total_student, ncol(stud_otus_w), i))
      d_total_school <- d_total_school + d_avg_stud
    }
    d_avg_school <- d_total_school / ncol(stud_otus_w) #num students in given school
    svws.euc[svwsind] <- d_avg_school
    svwsind <- svwsind + 1
  }
  
}
avg_svws <- mean(svws.euc)



# ********************** #
#students vs all other students
stud_all <- data.frame(otu_table(SLL1)[,rownames(students)])

svas.euc <- matrix(0, nrow=1, ncol=1473)
svasind <- 1
for (i in 1:ncol(stud_all)) {
  d_total <- 0
  for (j in 1:ncol(stud_all)) {
    temp <- rbind(t(stud_all[,i]), t(stud_all[,j])) #distance caluclated between the rows of the object
    d <- dist(temp)
    # print(dist(temp))
    d_total <- d_total + d
  }
  d_avg <- d_total/ncol(stud_all)
  print(c(d_avg, d_total, ncol(stud_all),i))
  svas.euc[svasind] <- d_avg
  svasind <- svasind + 1
}

avg_svas <- mean(svas.euc)


# ***************** #
euc.diff1 <- t.test(tvs.euc, tvt.euc)
euc.diff2 <- t.test(tvs.euc, svas.euc)
euc.diff3 <- t.test(tvt.euc, svas.euc)





















# ****************************************************************************************************************** #
# ****************************************************************************************************************** #
### Hisayama paper comparisons: ####

if (region_type == "all") {
  phy_obj <- SLL
  phy_obj_rel <- SLL1
} else if (region_type == "region") {
  phy_obj <- SLL.region
  phy_obj_rel <- SLL1.region
}

# they found 72 OTUs present in >= 75% of samples ("common" OTUs)
wh.in75 <- genefilter_sample(phy_obj, filterfun_sample(function(x) x > 0), A = 0.75 * nsamples(phy_obj))
summary(wh.in75)
# we had 44
phy_obj.75 <- prune_taxa(wh.in75, phy_obj)
phy_obj.75rel <- transform_sample_counts(phy_obj.75, function(x) 100 * x/sum(x))


# ********************************************* #
# what percentage of full microbiome is made up of our "common" OTUs
phy_obj.sums <- apply( otu_table(phy_obj), 2, sum )
phy_obj_rel.sums <- apply( otu_table(phy_obj)[ rownames(otu_table(phy_obj.75)), ], 2, sum )
full_percents <- phy_obj_rel.sums / phy_obj.sums
mean_per <- 100 * mean(full_percents)
sd_per <- 100 * sd(full_percents)
# so on average, these 44 genera that are present in at least 75% of samples make up 99.4% of the abundance


# ********************************************* #
# Their measure of diversity was phylogenetic diversity (PD), they divided into quintiles
# We have a number of different diversity measures, so will take a look at each here:

# div.estimate <- "Div.Observed"
# div.estimate <- "Div.Chao1"
# div.estimate <- "Div.ACE"
div.estimate <- "Div.Shannon"
# div.estimate <- "Div.Simpson"
# div.estimate <- "Div.InvSimpson"
# div.estimate <- "Div.Fisher"

# divs <- lapply( sample_data(phy_obj_rel)[,div.estimate], as.numeric )
divs <- as.numeric(sample_data(phy_obj_rel)[,div.estimate][[1]])
quints <- quantile(divs, probs=seq(0,1,0.20))

q1 <- otu_table(phy_obj)[ ,rownames(sample_data(phy_obj_rel)[ sample_data(phy_obj_rel)[,div.estimate] < quints['20%'], ] ) ]
q2 <- otu_table(phy_obj)[ ,rownames(sample_data(phy_obj_rel)[ (quints['20%'] <= sample_data(phy_obj_rel)[,div.estimate]) & (sample_data(phy_obj_rel)[,div.estimate] < quints['40%']), ] ) ]
q3 <- otu_table(phy_obj)[ ,rownames(sample_data(phy_obj_rel)[ (quints['40%'] <= sample_data(phy_obj_rel)[,div.estimate]) & (sample_data(phy_obj_rel)[,div.estimate] < quints['60%']), ] ) ]
q4 <- otu_table(phy_obj)[ ,rownames(sample_data(phy_obj_rel)[ (quints['60%'] <= sample_data(phy_obj_rel)[,div.estimate]) & (sample_data(phy_obj_rel)[,div.estimate] < quints['80%']), ] ) ]
q5 <- otu_table(phy_obj)[ ,rownames(sample_data(phy_obj_rel)[ sample_data(phy_obj_rel)[,div.estimate] >= quints['80%'], ] ) ]


# filter out "uncommon" OTUs from each quintile, change to a data.frame bc only way to add rows later
wh.q5.in75 <- genefilter_sample(q5, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(q5))
q5.1 <- prune_taxa(wh.q5.in75, q5)
q5.1 <- as.data.frame(transform_sample_counts(q5.1, function(x) 100 * x/sum(x)))
wh.q4.in75 <- genefilter_sample(q4, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(q4))
q4.1 <- prune_taxa(wh.q4.in75, q4)
q4.1 <- as.data.frame(transform_sample_counts(q4.1, function(x) 100 * x/sum(x)))
wh.q3.in75 <- genefilter_sample(q3, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(q3))
q3.1 <- prune_taxa(wh.q3.in75, q3)
q3.1 <- as.data.frame(transform_sample_counts(q3.1, function(x) 100 * x/sum(x)))
wh.q2.in75 <- genefilter_sample(q2, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(q2))
q2.1 <- prune_taxa(wh.q2.in75, q2)
q2.1 <- as.data.frame(transform_sample_counts(q2.1, function(x) 100 * x/sum(x)))
wh.q1.in75 <- genefilter_sample(q1, filterfun_sample(function(x) x > 0), A = 0.33 * nsamples(q1))
q1.1 <- prune_taxa(wh.q1.in75, q1)
q1.1 <- as.data.frame(transform_sample_counts(q1.1, function(x) 100 * x/sum(x)))




# names of OTUs only found in a given quintile or above
q1.orhigher.sort <- sort(rownames(q1.1))
q2.orhigher.sort <- sort(rownames(q2.1)[ !rownames(q2.1) %in% rownames(q1.1) ])
q3.orhigher.sort <- sort(rownames(q3.1)[ !rownames(q3.1) %in% c(rownames(q1.1), rownames(q2.1)) ])
q4.orhigher.sort <- sort(rownames(q4.1)[ !rownames(q4.1) %in% c(rownames(q1.1), rownames(q2.1), rownames(q3.1)) ])
q5.orhigher.sort <- sort(rownames(q5.1)[ !rownames(q5.1) %in% c(rownames(q1.1), rownames(q2.1), rownames(q3.1), rownames(q4.1)) ])

# keep names grouped by quintile, and order the names by indicated taxa level
tl <- "ta2"
q1.orhigher <- rownames(sort(tax_table(phy_obj_rel)[q1.orhigher.sort, tl]))
q2.orhigher <- rownames(sort(tax_table(phy_obj_rel)[q2.orhigher.sort, tl]))
q3.orhigher <- rownames(sort(tax_table(phy_obj_rel)[q3.orhigher.sort, tl]))
q4.orhigher <- rownames(sort(tax_table(phy_obj_rel)[q4.orhigher.sort, tl]))
q5.orhigher <- rownames(sort(tax_table(phy_obj_rel)[q5.orhigher.sort, tl]))

# add rows with 0s for those OTUs that don't appear in a given quintile
q1.1[ c(q2.orhigher,q3.orhigher,q4.orhigher,q5.orhigher), ] <- 0
q2.1[ c(q3.orhigher,q4.orhigher,q5.orhigher), ] <- 0
q3.1[ c(q4.orhigher,q5.orhigher), ] <- 0
q4.1[ q5.orhigher, ] <- 0

# matrix of detection rate within each quintile
detection <- matrix(NA, nrow=length(rownames(q1.1)), ncol=5)
colnames(detection) <- c(sprintf("Q1\n<%s\n(n=%s)", round(quints['20%'],3), dim(q1.1)[2]),
                         sprintf("Q2\n<%s\n(n=%s)", round(quints['40%'],3), dim(q2.1)[2]),
                         sprintf("Q3\n<%s\n(n=%s)", round(quints['60%'],3), dim(q3.1)[2]),
                         sprintf("Q4\n<%s\n(n=%s)", round(quints['80%'],3), dim(q4.1)[2]),
                         sprintf("Q5\n<%s\n(n=%s)", round(quints['100%'],3), dim(q5.1)[2]) )
rownames(detection) <- c(q1.orhigher, q2.orhigher, q3.orhigher, q4.orhigher, q5.orhigher)


quint.otus <- c(q1.1, q2.1, q3.1, q4.1, q5.1)
for (i in rownames(detection)) {
  detection[ i, 1 ] <- 100 * as.numeric(summary(as.numeric(q1.1[i,]) != 0)['TRUE']) / dim(q1.1)[2]
  detection[ i, 2 ] <- 100 * as.numeric(summary(as.numeric(q2.1[i,]) != 0)['TRUE']) / dim(q2.1)[2]
  detection[ i, 3 ] <- 100 * as.numeric(summary(as.numeric(q3.1[i,]) != 0)['TRUE']) / dim(q3.1)[2]
  detection[ i, 4 ] <- 100 * as.numeric(summary(as.numeric(q4.1[i,]) != 0)['TRUE']) / dim(q4.1)[2]
  detection[ i, 5 ] <- 100 * as.numeric(summary(as.numeric(q5.1[i,]) != 0)['TRUE']) / dim(q5.1)[2]
}
detection[is.na(detection)] <- 0 # bc any with a detection rate of 0 will not have a column for 'TRUE' in the summary

det.m <- melt(detection)

# To print phylum name alongside genus name
get_phylum <- function(x) {
  sprintf("%s (%s)", x, tax_table(phy_obj_rel)[x, tl])
}

if (region_type == "all") {
  ti <- sprintf("Detection rate of 'common' OTUs among each quintile\n(%s)", div.estimate)
} else if (region_type == "region") {
  ti <- sprintf("Detection rate of 'common' OTUs among each quintile\n%s (%s)", reg_name, div.estimate)
}

ggplot(data = det.m, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=50, name="Detection rate (%)") +
  scale_y_discrete(limits=rev(levels(det.m$Var1)), labels=get_phylum(rev(levels(det.m$Var1))))+#, position="right") +
  theme(plot.title = element_text(hjust=0.5)) +#, legend.position = "left") +
  ggtitle(ti) + xlab('') + ylab('')








# ********************************************* #
# check # OTUs with mean rel abund >= 1% in the lowest quintile
# In Hisayama they call these the "predominant core" OTUs
q1.mps <- rowSums(q1.1)/length(colnames(q1.1))
q1.sd <- rowSds(as.matrix(q1.1))
pred.core <- rownames(q1.1[q1.mps >= 0.5, ])

# What % does these OTUs make up from total
pred.core.sum <- apply(otu_table(phy_obj)[pred.core, ], 2, sum)
pred.core.per <- pred.core.sum / colSums(otu_table(phy_obj))
pred.core.mean <- 100 * mean(pred.core.per)
pred.core.sd <- 100 * sd(pred.core.per)


























# ****************************************************************************************************************** #
###### water hardness in Spain significance ######
# ****************************************************************************************************************** #

t.test(as.numeric(as.matrix(sample_data(SLL1)[sample_data(SLL1)[,"Stomatotype"]==1, "Water_hardness"])),
       as.numeric(as.matrix(sample_data(SLL1)[sample_data(SLL1)[,"Stomatotype"]==2, "Water_hardness"])) )

kruskal.test(as.numeric(as.matrix(sample_data(SLL1)[,"Stomatotype"])), 
             as.factor(as.matrix(sample_data(SLL1)[,"Water_hardness_group"])))

kruskal.test(as.numeric(as.matrix(sample_data(SLL1)[,"Water_hardness"])), 
             as.factor(as.matrix(sample_data(SLL1)[,"Stomatotype"])))


# sto_water <- sample_data(SLL1)[,c("Stomatotype","Water_hardness")]
# sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] != "Embotellada", c("Stomatotype","Water_hardness")]
# sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] != "Embotellada" & 
#                                  sample_data(SLL1)[ , "Q27"] != "No Sabe/No Contesta", c("Stomatotype","Water_hardness")]
# sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] == "No Tratada (Fuente, Pozo O Rio)", c("Stomatotype","Water_hardness")]
# sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] == "Del Grifo (No Filtrada)", c("Stomatotype","Water_hardness")]
# sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] == "Del Grifo (Filtrada)", c("Stomatotype","Water_hardness")]
sto_water <- sample_data(SLL1)[ sample_data(SLL1)[ , "Q27"] == "No Tratada (Fuente, Pozo O Rio)" | 
                                 sample_data(SLL1)[ , "Q27"] == "Del Grifo (No Filtrada)",
                               c("Stomatotype","Water_hardness")]

t.test(as.numeric(as.matrix(sto_water[sto_water[,"Stomatotype"]==1, "Water_hardness"])),
       as.numeric(as.matrix(sto_water[sto_water[,"Stomatotype"]==2, "Water_hardness"])) )

ggplot(sto_water, aes(x=as.character(Stomatotype), y=as.numeric(Water_hardness), fill=as.character(Stomatotype))) +
  geom_boxplot() +
  ggtitle('Stomatotype vs water hardness') + xlab("Stomatotype") + ylab('Water hardness') + scale_fill_hue(name="Stomatotype")










# ****************************************************************************************************************** #
###### Plotting on map of Spain ######
# ****************************************************************************************************************** #
library(maptools)
library(RColorBrewer)
library(plotrix)
# library(classInt)

# esp.shp <- readShapeSpatial("/users/tg/jwillis/SLL/ESP_adm_shp/ESP_adm0.shp")
esp.shp.com <- readShapeSpatial( sprintf("%s/ESP_adm_shp/ESP_adm1.shp", home_dir) )
esp.shp.prov <- readShapeSpatial( sprintf("%s/ESP_adm_shp/ESP_adm2.shp", home_dir) )
# esp.shp <- readShapeSpatial("/users/tg/jwillis/SLL/ESP_adm_shp/ESP_adm3.shp")
esp.shp.city <- readShapeSpatial( sprintf("%s/ESP_adm_shp/ESP_adm4.shp", home_dir) )
# plot(esp.shp)




plot_map <- function(g, fills, samp_data=sample_data(SLL1), otu_data=otu_table(SLL1)) {
  
  if (fills == "Provinces") { 
    shape <- esp.shp.prov
    shape <- shape[shape$NAME_2 != "Santa Cruz de Tenerife" & shape$NAME_2 != "Las Palmas", ]
    name <- shape$NAME_2
    loc_list <- provinces
  } else if (fills == "Autonomous Communities") { 
    shape <- esp.shp.com
    shape <- shape[shape$NAME_1 != "Islas Canarias", ]
    name <- shape$NAME_1
    loc_list <- comunidades
  }
  
  
  cont_vars <- c("Q00.PH","Div.Simpson","Div.Shannon","Div.Fisher","Faiths.PD","Species_Richness","Q22.2",
                 "Water_hardness","Weighted_Unifrac","Unweighted_Unifrac","Socioeconomic",cont_water_data)
  
  if (startsWith(g, "Stomatotype")) {
    # % of each autonomous community comprised of each stomatotype
    # By city according to school location (not by postal code of home address bc that was input by the students so may be unreliable):
    ## Apparently they are essentially the same anyway, with only 4 samples in differing locations
    ##    Will use this method for accuracy and because easier to match names to the shp file since I 
    ##    can define the names in this script, instead of relying on names imported from sample data
    
    shape$Stomatotype1 <- numeric(length=length(name))
    shape$Stomatotype2 <- numeric(length=length(name))
    
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        loc_count <- nrow( samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ] )
        s1_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
                                               rownames(samp_data)[samp_data[ , "Stomatotype"] == 1 ]) ] )
        s2_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
                                               rownames(samp_data)[samp_data[ , "Stomatotype"] == 2 ]) ] )
        s1_perc <- s1_count/loc_count
        s2_perc <- s2_count/loc_count
        print(c(loc, loc_count, s1_perc, s2_perc))
        
        attributes(shape)[["data"]][ which(name == loc), "Stomatotype1" ] <- s1_perc
        attributes(shape)[["data"]][ which(name == loc), "Stomatotype2" ] <- s2_perc
        
      } else {
        attributes(shape)[["data"]][ which(name == loc), "Stomatotype1" ] <- -1
        attributes(shape)[["data"]][ which(name == loc), "Stomatotype2" ] <- -1
      }
    }
    
    # colors <- c('black', rev(brewer.pal(11, "RdBu")))
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(0, 1, length.out=11))
    legend_labels = c("No data","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100") # leglabs(round(brks, digits=1))
    legend.title <- sprintf('%% of regions with\n%s', g)
    plot_title <- sprintf("Distribution of %s\namongst %s", g, fills)
    leg_x <- 1.75
  
  # } else if (g %in% c('soft','intermediate','hard','very_hard') ) {
  #   
  # 
  #   shape$hgroup <- numeric(length=length(name))
  #   shape$other_groups <- numeric(length=length(name))
  # 
  #   for (loc in name) {
  #     if (loc %in% names(loc_list)) {
  #       loc_count <- nrow( samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ] )
  #       hg_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
  #                                              rownames(samp_data)[samp_data[ , "Water_hardness_group"] == g ]) ] )
  #       og_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
  #                                              rownames(samp_data)[samp_data[ , "Water_hardness_group"] != g ]) ] )
  #       hg_perc <- hg_count/loc_count
  #       og_perc <- og_count/loc_count
  #       print(c(loc, loc_count, hg_perc, og_perc))
  # 
  #       attributes(shape)[["data"]][ which(name == loc), "hgroup" ] <- hg_perc
  #       attributes(shape)[["data"]][ which(name == loc), "other_groups" ] <- og_perc
  # 
  #     } else {
  #       attributes(shape)[["data"]][ which(name == loc), "hgroup" ] <- -1
  #       attributes(shape)[["data"]][ which(name == loc), "other_groups" ] <- -1
  #     }
  #   }
  # 
  #   # colors <- c('black', rev(brewer.pal(11, "RdBu")))
  #   colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
  #   brks <- c(-1, seq(0, 1, length.out=11))
  #   legend_labels = c("No data","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100") # leglabs(round(brks, digits=1))
  #   legend.title <- sprintf('%% of regions with\n%s water', g)
  #   plot_title <- sprintf("Distribution of %s water\namongst %s", g, fills)
  
  } else if (g %in% rownames(otu_data)) {
    
    ## For GENUS distributions
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        g_mean <- mean( otu_data[ g, samp_data$Qsamples.1 %in% loc_list[loc][[1]] ] )
        g_sd <- sd( otu_data[ g, samp_data$Qsamples.1 %in% loc_list[loc][[1]] ] )
        print(c(loc, g_mean, g_sd))
        attributes(shape)[["data"]][ which(name == loc), g ] <- g_mean
      } else {
        attributes(shape)[["data"]][ which(name == loc), g ] <- -1
      }
    }
    
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(sort(unique(attributes(shape)[["data"]][ , g]))[2], max(attributes(shape)[["data"]][ , g]), length.out=11))
    legend_labels = leglabs(round(brks, digits=2), over=">")
    legend_labels[1] <- "No data"
    legend_labels[2] <- sprintf("< %s", round(brks[3], digits=3))
    legend.title <- sprintf("Mean rel abund\nof %s", g)
    plot_title <- sprintf("Mean Relative Abundance of %s\namongst %s", g, fills)
    leg_x <- 1.65
    
  } else if (g %in% cont_vars) {
    
    ## For distributions of continuous variables
    print(c("Region", "# samples", "mean", "sd"))
    
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        cont_vals <- as.numeric(as.matrix(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], g ] ))
        cont_mean <- mean( cont_vals, na.rm=T )
        cont_sd <- sd( cont_vals, na.rm=T )
        print(c(loc, length(cont_vals), cont_mean, cont_sd))
        attributes(shape)[["data"]][ which(name == loc), g ] <- cont_mean
      } else {
        attributes(shape)[["data"]][ which(name == loc), g ] <- -1
      }
    }
    
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(sort(unique(attributes(shape)[["data"]][ , g]))[2], max(attributes(shape)[["data"]][ , g]), length.out=11))
    legend_labels = leglabs(round(brks, digits=3), over=">")
    legend_labels[1] <- "No data"
    legend_labels[2] <- sprintf("< %s", round(brks[3], digits=3))
    if (g == "Q22.2") {
      legend.title <- sprintf("Mean milk consumption\nwithin region")
      plot_title <- sprintf("Mean milk consumption amongst %s", fills)
    } else {
      legend.title <- sprintf("Mean %s\nwithin region", g)
      plot_title <- sprintf("Mean %s amongst %s", g, fills)
      leg_x <- 1.1
    }
    
  } else if (g %in% sort(as.character(as.matrix(unique(samp_data[,"Q27"]))))[-4]) {
    
    ## For different types of WATER consumed at home
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        loc_count <- nrow( samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ] )
        agua_count <- tryCatch(nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
                                                       rownames(samp_data[ samp_data[, "Q27"]==g ])) ] ),
                               error = function(x) 0) # In case there is an area/agua combo for which intersect returns no samples, agua_count becomes 0
        
        agua_perc <- (agua_count / loc_count)*100
        print(c(loc, loc_count, agua_perc))
        attributes(shape)[["data"]][ which(name == loc), g ] <- agua_perc
      } else {
        attributes(shape)[["data"]][ which(name == loc), g ] <- -1
      }
    }
    
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(0, 75, length.out=11))
    # brks <- c(-1, seq(sort(unique(attributes(shape)[["data"]][ , g]))[2], max(attributes(shape)[["data"]][ , g]), length.out=9))
    legend_labels = leglabs(round(brks, digits=2), over=">")
    legend_labels[1] <- "No data"
    legend_labels[2] <- sprintf("< %s", round(brks[3], digits=3))
    legend.title <- sprintf("%% of region that drinks\nwater %s", g)
    plot_title <- sprintf("Distribution of Agua %s\namongst %s", g, fills)
    
  } else {
    
    ## For distribution of "yes" responses in yes/no questions
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        ###  for g==Q19.2
        # g_mean <- mean( as.numeric(as.matrix(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], g ] )) )
        # g_sd <- sd( as.numeric(as.matrix(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], g ] )) )
        # print(c(loc, g_mean, g_sd))
        # attributes(shape)[["data"]][ which(name == loc), g ] <- g_mean
        
        ###  for g==Q19.1
        loc_count <- nrow( samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ] )
        g_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$Qsamples.1 %in% loc_list[loc][[1]], ]),
                                                       rownames(samp_data[ samp_data[, g]==1 ])) ] )
        g_perc <- (g_count / loc_count)*100
        print(c(loc, loc_count, g_perc))
        attributes(shape)[["data"]][ which(name == loc), g ] <- g_perc
      } else {
        attributes(shape)[["data"]][ which(name == loc), g ] <- -1
      }
    }
    
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(sort(unique(attributes(shape)[["data"]][ , g]))[2], max(attributes(shape)[["data"]][ , g]), length.out=11))
    legend_labels = leglabs(round(brks, digits=2), over=">")
    legend_labels[1] <- "No data"
    legend_labels[2] <- sprintf("< %s", round(brks[3], digits=3))
    
    if (g == "Q19.1") {
      legend.title <- sprintf("%% of region that owns a dog", g)
      plot_title <- sprintf("Distribution of dogs\namongst %s", fills)
    } else if (g == "Q21") {
      legend.title <- sprintf("%% of region with a\nsmoker in the house", g)
      plot_title <- sprintf("Distribution of homes with smokers\namongst %s", fills)
    }
    
  }
  
  # to shade with stripes those regions with no data
  angles <- c( 45, rep(0, 10) )
  densities <- c( 30, rep(250,10) )
  
  column <- attributes(shape)[["data"]][ , g]
  color_intervals <- colors[ findInterval(column, brks, rightmost.closed = T) ]
  angle_intervals <- angles[ findInterval(column, brks, rightmost.closed = T) ]
  density_intervals <- densities[ findInterval(column, brks, rightmost.closed = T)]
  coord <- coordinates(shape)[column != -1, ] #so we dont print the names of those areas for which there is no data
  printNames <- name[column != -1]            #  ^^^
  
  
  plot(shape, col=color_intervals, angle=angle_intervals, density=density_intervals)
  # boxed.labels(x=coord[,1], y=coord[,2], labels=printNames, bg="white", border=FALSE, cex=1)
  legend(x=leg_x, y=38.8, legend=legend_labels, fill=colors, bty="n", title=legend.title, 
         angle=angles, density=densities, cex = 1.45)
  # title(main=plot_title, adj = 0.5, line=-3)
  title(main=g, adj = 0.5, line=-2)
}


# fills <- "Provinces"
fills <- "Autonomous Communities"

# hardnesses <- c('soft','intermediate','hard','very_hard')
# g <- hardnesses[1] # this takes the categorical value
# g <- "Water_hardness" # this takes the numerical value
g <- "Stomatotype1"
# g <- "Veillonella"
# g <- "Q00.PH"
# g <- "Q22.2"
# g <- "Q19.1"
###[1] Del Grifo (Filtrada)  [2] Del Grifo (No Filtrada)  [3] Embotellada  [4] No Tratada (Fuente, Pozo O Rio):
# g <- sort(as.character(as.matrix(unique(sample_data(SLL1)[,"Q27"]))))[-4][1]
# g <- "Div.Shannon"
g <- "Weighted_Unifrac"
# g <- "Faiths.PD"
# g <- cont_water_data[18] # there are 18 values here

plot_map(g, fills)

# ****************************************************************************************************************** #








# ****************************************************************************************************************** #
###### Prepare table of sample attributes for Sequence Reads Archive ######
# ****************************************************************************************************************** #
attribs <- c("sample_name","bioproject_accession","organism","ecotype","collection_date","env_biome","env_feature",
             "env_material","geo_loc_name","host","isol_growth_condt","lat_lon","chem_administration",
             "host_age","host_body_mass_index","host_body_product","host_height","host_occupation","host_sex",
             "host_tissue_sampled","host_tot_mass","isolation_source","samp_collect_device","samp_store_temp")

# collection months
feb <- c("07","08","10","12","15","19","22","25","26","27","28","29","35","39")
mar <- c("06","09","11","16","17","18","20","21","24","30","31","32","33","36","37","38","40")
apr <- c("13","14","23","34","41")

# to get geo_loc_name
gln <- function(x) {
  sn <- as.numeric(as.matrix(sample_data(SLL1))[x,"Qsamples.1"])
  cit <- sapply(names(cities), function(y) sn %in% cities[[y]])
  c <- gsub("á","a", names(cit[cit]))
  c <- gsub("í","i", c)
  c <- gsub("ú","u", c)
  c <- gsub("ó","o", c)
  sprintf("Spain: %s", c)
}

samp_lat_lon <- function(x) {
  sn <- as.numeric(as.matrix(sample_data(SLL1))[x,"Qsamples.1"])
  ll <- sapply(names(lat_lon), function(y) sn %in% lat_lon[[y]])
  names(ll[ll])
}

attr.ages <- as.matrix(sample_data(SLL1))[,"Age"]
attr.bmi <- as.matrix(sample_data(SLL1))[,"Q3.1"]
attr.height <- as.matrix(sample_data(SLL1))[,"Q4"]
attr.occup <- as.matrix(sample_data(SLL1))[,"Q1.1"]
attr.sex <- as.matrix(sample_data(SLL1))[,"Q2"]
attr.mass <- as.matrix(sample_data(SLL1))[,"Q3"]

s700 <- sort(sample_names(SLL1))[1:700]
s619 <- sort(sample_names(SLL1))[701:1319]

for (lab in c("s1-s700","s701-s1319")) {
  if (lab=="s1-s700") {
    samp_vec <- s700
  } else {
    samp_vec <- s619
  }
  
  attr.table <- matrix(NA, nrow = length(samp_vec), ncol = length(attribs))
  rownames(attr.table) <- samp_vec
  colnames(attr.table) <- attribs
  
  for (samp in samp_vec) {
    
    attr.table[ samp, "sample_name" ] <- samp
    attr.table[ samp, "bioproject_accession" ] <- "PRJNA427101"
    attr.table[ samp, "organism" ] <- "Human oral metagenome"
    attr.table[ samp, "ecotype" ] <- "Spain"
    
    attr.table[ samp, "collection_date" ] <- ifelse( substr(samp,4,5) %in% feb, "Feb-2015", ifelse(substr(samp,4,5) %in% feb, "Mar-2015", "Apr-2015") )
    attr.table[ samp, "env_biome" ] <- "Human saliva"
    attr.table[ samp, "env_feature" ] <- "Human saliva"
    attr.table[ samp, "env_material" ] <- "Human saliva"
    attr.table[ samp, "geo_loc_name" ] <- gln(samp)
    attr.table[ samp, "host" ] <- "Homo sapiens"
    attr.table[ samp, "isol_growth_condt" ] <- "not applicable"
    attr.table[ samp, "lat_lon" ] <- samp_lat_lon(samp)
    
    # attr.table[ samp, "biotic_relationship" ] <- "from host"
    attr.table[ samp, "chem_administration" ] <- "Phosphate buffered saline (PBS)"
    attr.table[ samp, "host_age" ] <- ifelse( is.na(as.numeric(attr.ages[samp])), NA, sprintf("%s years", as.numeric(attr.ages[samp])) )
    attr.table[ samp, "host_body_mass_index" ] <- ifelse( is.na(as.numeric(attr.bmi[samp])), NA, as.numeric(attr.bmi[samp]) )
    attr.table[ samp, "host_body_product" ] <- "Saliva"
    attr.table[ samp, "host_height" ] <- ifelse( is.na(as.numeric(attr.height[samp])), NA, sprintf("%s cm", as.numeric(attr.height[samp])) )
    attr.table[ samp, "host_occupation" ] <- as.character(attr.occup[samp])
    attr.table[ samp, "host_sex" ] <- ifelse(as.character(attr.sex[samp])==0, "Female", ifelse(as.character(attr.sex[samp])==1, "Male", NA))
    attr.table[ samp, "host_tissue_sampled" ] <- "Saliva"
    attr.table[ samp, "host_tot_mass" ] <- ifelse( is.na(as.numeric(attr.mass[samp])), NA, sprintf("%s kg", as.numeric(attr.mass[samp])) )
    attr.table[ samp, "isolation_source" ] <- "Saliva"
    attr.table[ samp, "samp_collect_device" ] <- "Phosphate buffered saline (PBS)"
    attr.table[ samp, "samp_store_temp" ] <- "-20 C"
  }
  
  write.csv(attr.table, sprintf('%s/Part_1/SLL1_paper/SRA_submission/attr_tab_%s.csv', home_dir, lab), row.names = F)
}




# ****************************************************************************************************************** #
###### Prepare table of metadata for files for Sequence Reads Archive ######
# ****************************************************************************************************************** #
meta_labs <- c("bioproject_accession","biosample_accession","library_ID","title","library_strategy","library_source","library_selection",
             "library_layout","platform","instrument_model","design_description","filetype","filename","filename2")

# ****** #
# read sample accession numbers from the files provided by SRA after submitting the attributes info above
s700 <- sort(sample_names(SLL1))[1:700]
s619 <- sort(sample_names(SLL1))[701:1319]

samp_acc_s1_s700 <- read.delim(sprintf("%s/Part_1/SLL1_paper/SRA_submission/Sample_accessions_s1-s700.tsv", home_dir))
samp_acc_s701_s1319 <- read.delim(sprintf("%s/Part_1/SLL1_paper/SRA_submission/Sample_accessions_s701-s1319.tsv", home_dir))

get_samp_acc <- function(x) {
  if (x %in% s700) {
    samp_acc <- as.character( samp_acc_s1_s700[ samp_acc_s1_s700[,"sample_name"]==x, "accession" ] )
  } else {
    samp_acc <- as.character( samp_acc_s701_s1319[ samp_acc_s701_s1319[,"sample_name"]==x, "accession" ] )
  }
  return(samp_acc)
}

# ****** #
# read file containing the fastq filenames for each sample
filenames <- read.delim(sprintf("%s/Part_1/SLL1_paper/SRA_submission/filenames.csv", home_dir), row.names = 1)

# ****** #
for (lab in c("s1-s700","s701-s1319")) {
  if (lab=="s1-s700") {
    samp_vec <- s700
  } else {
    samp_vec <- s619
  }
  
  meta.table <- matrix(NA, nrow = length(samp_vec), ncol = length(meta_labs))
  rownames(meta.table) <- samp_vec
  colnames(meta.table) <- meta_labs
  
  for (samp in samp_vec) {
    
    meta.table[ samp, "bioproject_accession" ] <- "PRJNA427101"
    meta.table[ samp, "biosample_accession" ] <- get_samp_acc(samp)
    meta.table[ samp, "library_ID" ] <- samp
    meta.table[ samp, "title" ] <- sprintf("16S rRNA sequencing reads from sample %s", samp)
    
    meta.table[ samp, "library_strategy" ] <- "AMPLICON"
    meta.table[ samp, "library_source" ] <- "METAGENOMIC"
    meta.table[ samp, "library_selection" ] <- "PCR"
    meta.table[ samp, "library_layout" ] <- "paired"
    meta.table[ samp, "platform" ] <- "ILLUMINA"
    meta.table[ samp, "instrument_model" ] <- "Illumina MiSeq"
    
    meta.table[ samp, "design_description" ] <- "16S rRNA amplicon sequencing of the V3-V4 regions"
    meta.table[ samp, "filetype" ] <- "fastq"
    meta.table[ samp, "filename" ] <- as.character(filenames[samp, "File1"])
    meta.table[ samp, "filename2" ] <- as.character(filenames[samp, "File2"])
    
  }

  write.csv(meta.table, sprintf('%s/Part_1/SLL1_paper/SRA_submission/SRA_metadata_%s.csv', home_dir, lab), row.names = F)
}
# ****** #












# ****************************************************************************************************************** #
###### Fungal Analyses ######
# ****************************************************************************************************************** #
truecont <- c(4,5,36,37,41,43,46,49,52,59:73,100,102,106:144,149:155,  156)
#excel col:  (D,E,AJ,AK,AO,AQ,AT,AW,AZ,BG:BU,CV,  CX, DB:EN ,div-stats,#otus)

if (exists("samples.by.cluster")) {
  # meas.cont <- cbind(t(otu_table(SLL1)), sample_data(SLL1)[,truecont], samples.by.cluster$otu.cluster)
  meas.cont <- cbind(t(otu_table(SLL1)), sample_data(SLL1)[,only_cont], samples.by.cluster$otu.cluster)
  colnames(meas.cont)[dim(meas.cont)[2]] <- "Stomatotype"
} else {
  # meas.cont <- cbind(t(otu_table(SLL1)), sample_data(SLL1)[,truecont])
  meas.cont <- cbind(t(otu_table(SLL1)), sample_data(SLL1)[,only_cont])
}
meas.cont[meas.cont=='No Sabe/No Contesta'] <- NA
meas.cont <- apply(meas.cont, 2, as.numeric)
rownames(meas.cont) <- rownames(sample_data(SLL1))


meas <- "Q3.1"
meas_matrix <- matrix( nrow=length(fungi.list), ncol=max(sapply(fungi.list, length)) )
rownames(meas_matrix) <- names(fungi.list)
i <- 0

for (fung in fungi.list) {
  i <- i + 1
  fname <- names(fungi.list[i])
  
  # SLL1.fung <- as.data.frame(otu_table(SLL1)[meas,fung])
  # meas_matrix[fname,1:ncol(SLL1.fung)] <- as.numeric(SLL1.fung)
  
  SLL1.fung <- meas.cont[fung,meas]
  meas_matrix[fname,1:length(SLL1.fung)] <- SLL1.fung
}
mm.melt <- melt(t(meas_matrix))

# labels for the plot
ti <- sprintf('Relative abundance of %s in\nsamples for which given fungus is present', meas)
if (meas %in% taxa_names(SLL1)) {
  yl <- '% per sample'
} else {
  yl <- meas
}

counts <- function(x) {
  return( data.frame(y=max(x)+(0.05*max(x)), label=length(x)) )
}

# makes box plots of the values of meas (can be an OTU or a continuous survey variable)
# boxes are ordered by the median values
# numbers above boxes indicate the number of samples in which the given fungus are present
ggplot(mm.melt, aes(x=reorder(Var2, -value, FUN=median, na.rm=T), y=value, fill=reorder(Var2, -value, FUN=median, na.rm=T))) +
  geom_boxplot() +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) + 
  ggtitle(ti) +
  xlab("Fungus") + ylab(yl) + scale_fill_hue(name="Fungi") +
  stat_summary(fun.data = counts, geom='text')




## Run t-tests

# which alternative hypothesis for t-test
gl <- "greater"
# gl <- "less"

## Run correlations between the presence of particular fungi and all taxa
res.fungi <- matrix(NA, nrow=ncol(meas.cont), ncol=length(fungi.list))
colnames(res.fungi) <- names(fungi.list)
rownames(res.fungi) <- colnames(meas.cont)
ps.fungi <- matrix(NA, nrow=ncol(meas.cont), ncol=length(fungi.list))
colnames(ps.fungi) <- names(fungi.list)
rownames(ps.fungi) <- colnames(meas.cont)

for (i in colnames(meas.cont)) {
  for (j in names(fungi.list)) {
    
    if (length( fungi.list[[j]] ) > 1 ) {
      #check those samples containing a given fungus against the other samples
      # ... there may be some samples appearing in each group since some samples are listed under more than one fungus
      #....... should these be REMOVED???? ...
      #...try both ways
      this_fungi <- unname(unlist( fungi.list[j] ))
      other_fungi <- unname(unlist( fungi.list[names(fungi.list) != j] ))
      
      if (summary(! is.na(meas.cont[this_fungi, i]))["TRUE"] > 1 & summary(! is.na(meas.cont[other_fungi, i]))["TRUE"] > 1) {
        # because for some questions, there may be NAs and if only 1 or 0 samples have non-NA, cant run the t-test
        tt <- t.test(meas.cont[this_fungi, i], meas.cont[other_fungi, i], alternative = gl, na.rm=T)
        
        res.fungi[i,j] <- tt$statistic
        ps.fungi[i,j] <- tt$p.value
      }
    }
  }
}

# adjust for multiple testing
ps.fungi.adj <- apply(ps.fungi, 2, p.adjust, method='bonferroni', n=ncol(meas.cont))

#at least 1 good p value within question response
mins <- apply(ps.fungi.adj,2,min, na.rm=T)
goodPs <- mins[mins < min.p.val] #-log10(0.05) == 1.3013
#at least 1 good p within phenotype
tmins <- apply(ps.fungi.adj,1,min, na.rm=T)
tgoodPs <- tmins[tmins < min.p.val]

# a matrix with only those rows and columns that had a least 1 significant p-val
# use drop=F for those cases in which only 1 row or column had at least 1 significant p-value...maintains structure as a matrix
good.ps.fungi.adj <- ps.fungi.adj[names(tgoodPs), names(goodPs), drop = F]

# any non-significant p-vals will be a blank entry
only.good.ps.fungi.adj <- good.ps.fungi.adj
only.good.ps.fungi.adj[only.good.ps.fungi.adj >= min.p.val] <- ''

write.csv( only.good.ps.fungi.adj, sprintf('%s/Part_1/figures/Fungi/ttest_fungi.%s.signif.csv', home_dir, gl) )









# ****************************************************************************************************************** #
###### Prep for METAGENassist ######
# ****************************************************************************************************************** #

# get abundances/counts in a csv file
# must not quote values in csv because METAGENassist will not read file properly
write.csv(t(otu_table(SLL1)), sprintf('%s/Part_1/METAGENassist/SLL1_rel_abund.csv', home_dir), quote = F)
write.csv(t(otu_table(SLL)), sprintf('%s/Part_1/METAGENassist/SLL_total_counts.csv', home_dir), quote = F)

# ********* #
# metadata file must be qualitative variables only
qual_only <- c(3,6:33,35,38:40,42,44,45,47,48,51,53:58,74:98,100,102:106,108:148,           157:161)
#             (C,F:AG,AI,AL:AN,AP,AR,AS,AU,AV,AY,BA:BF,BV:CT,CV,  CX:DB,  DD:ER, Fungi,Stomatotype,Diversity_group,Water_hardness,Water_hardness_group)

# get columns from metadata_table to change 0/1 to yes/no
#  must shift down 1 from what appears in the csv file since rownames are included as column 1
yes_no_qs <- c(12:29,31:35,37,39,46:49,51:55,60:70,seq(77,111,2))
#             ( M:AD,AF:AJ,AL,AN,AU:AX,AZ:BC,BI:BS,BZ:DH(alternating))

# have to change commas to semicolons in the metadata file in order to save it with commas as separators and
#   not have to use quotes when printing to file to maintain values containing commas
#   METAGENassist will not read file properly if values are quoted
conv_com_to_semicol <- function(x) {
  gsub(',', ' ', x)
}

metadata_table <- as.matrix(sample_data(SLL1)[ , qual_only])
metadata_table[metadata_table == 'No Sabe/No Contesta'] <- NA
metadata_table[,yes_no_qs] <- sapply(metadata_table[,yes_no_qs], function(x) ifelse(x==0,"No",x))
metadata_table[,yes_no_qs] <- sapply(metadata_table[,yes_no_qs], function(x) ifelse(x==1,"Yes",x))
metadata_table[,"Q2"] <- sapply(metadata_table[,"Q2"], function(x) ifelse(x==0,"Female",x))
metadata_table[,"Q2"] <- sapply(metadata_table[,"Q2"], function(x) ifelse(x==1,"Male",x))
Fungi_present <- sapply(metadata_table[,"Fungi"], function(x) ifelse(x==0,"No","Yes"))
metadata_table <- cbind( metadata_table, Fungi_present)
# metadata_table[metadata_table == 0] <- "Zero"
# metadata_table[metadata_table == 1] <- "One"
# metadata_table[metadata_table == 2] <- "Two"
# metadata_table[metadata_table == 3] <- "Three"
# metadata_table[metadata_table == 4] <- "Four"
# metadata_table[metadata_table == 5] <- "Five"
# metadata_table[metadata_table == 6] <- "Six"
# metadata_table[metadata_table == 7] <- "Seven"
metadata_table <- apply(metadata_table, 2, conv_com_to_semicol)
write.csv(metadata_table, sprintf('%s/Part_1/METAGENassist/metadata.csv', home_dir), quote = F)
# ********* #

# ********* #
## Then do these by separating into each stomatotype
s1_samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[, "Stomatotype"] == 1, ])
s2_samples <- rownames(sample_data(SLL1)[ sample_data(SLL1)[, "Stomatotype"] == 2, ])

write.csv(t(otu_table(SLL)[,s1_samples]), sprintf('%s/Part_1/METAGENassist/SLL_total_counts_stomatotype1.csv', home_dir), quote = F)
write.csv(t(otu_table(SLL)[,s2_samples]), sprintf('%s/Part_1/METAGENassist/SLL_total_counts_stomatotype2.csv', home_dir), quote = F)


# ********* #
# the phenotypes database from METAGENassist is very large, with maaaany species, 
# so I will narrow it down and save a smaller version with just those genera that appear in our study
large_phenos <- read.delim(sprintf('%s/Part_1/METAGENassist/biography_cards.all.inference.tsv', home_dir), row.names = 1)
all_phenos_table <- large_phenos[ large_phenos$GENUS %in% taxa_names(SLL), ]

rm(large_phenos)

# all_phenos_table <- drop.levels(all_phenos_table)
all_phenos_table <- as.matrix(all_phenos_table)
all_phenos_table[is.na(all_phenos_table)] <- "NULL"
all_phenos_table <- as.data.frame(all_phenos_table)


### Keep only one row for each genus, selecting the most common response in each phenotype column amongst the rows for that genus
new_phenos_table <- matrix(NA, nrow=length(unique(all_phenos_table$GENUS)), ncol=ncol(all_phenos_table))
rownames(new_phenos_table) <- unique(all_phenos_table$GENUS)
colnames(new_phenos_table) <- colnames(all_phenos_table)


# *************************** #
get_phenotypes <- function(x, genus) {
  
  min.percent.species <- 0.75 #threshold percentage for phenotypes to appear in superior taxonomic level
  
  # check if at least 1 entry has multiple values separated by ';'
  # some fields, like "DISEASES" or "METABOLISM" may have multiple responses, will collect all, and take any of those appearing in 90% of species
  if (sum( grepl(';', x)) > 0) {
    
    # get vector of all responses
    all_Ps <- unique( unlist( strsplit(as.character(x), ';') ) )
    # check for those appearing in at least 90% of species
    good_Ps <- sapply(all_Ps, function(P, vec) ifelse(sum( ifelse(grepl(';',vec), 
                                                                  grepl(paste0(P,';'), vec) | grepl(paste0(';',P),vec), 
                                                                  grepl(P,vec))) / length(vec) > min.percent.species, 
                                                      P, "NULL"), x)
    good_Ps <- unique(unname(good_Ps))
    good_Ps <- good_Ps[good_Ps!="NULL"]
    good_Ps <- good_Ps[good_Ps!=""]
    
    # combine them again by ';' so can be handled as before
    final_Ps <- paste(good_Ps, collapse = ';')
    return(final_Ps)
    
  } else {
    # if not, check for most common phenotype for given category
    most_common <- names(sort(table(x), decreasing = T)[1])
    
    if ( sum( x == most_common) / length(x) > min.percent.species) {
      # most common phenotype must appear in at least 90% of species
      return(most_common)
    } else {
      return("UNKNOWN")
    }
  }
}
# *************************** #

for (genus in rownames(new_phenos_table)) {
  genus_table <- all_phenos_table[all_phenos_table$GENUS == genus, ]
  # # to remove those undetermined species (e.g. "Microbacterium sp. EB93)
  # genus_table <- genus_table[ ! startsWith(genus_table$SPECIES, 'sp.'), ]
  # # to remove those multiple strains for a given species (e.g. "Acetobacter pasteurianus IFO 3283-07")
  # genus_table <- genus_table[ sapply(sapply(rownames(genus_table), strsplit, ' '), length) == 2, ]

  # new_phenos_table[genus, ] <- unname(unlist(as.character( apply( genus_table, 2, function(x) names(sort(table(x), decreasing = T)[1])) )))
  new_phenos_table[genus, ] <- unname(unlist(as.character( apply( genus_table, 2, get_phenotypes, genus) )))
}

### First remove species column, since we only have genus data
new_phenos_table <- new_phenos_table[ , colnames(new_phenos_table) != "SPECIES"]
new_phenos_table[new_phenos_table == "NULL"] <- "UNKNOWN"

write.csv(new_phenos_table, sprintf("%s/Part_1/METAGENassist/phenotypes_all.csv", home_dir))





# Create large list of data tables
allTables <- list()
allTables$raw <- list()
allTables$norm <- list()


for (p in colnames(new_phenos_table)) {
  
  # get all the unique values for each phenotype (some are in a list separated by ';' or ',')
  p_cats <- unique(new_phenos_table[ , p])
  p_cats <- sort(unique(unlist( sapply(p_cats, function(x) strsplit(x, ';')) )))
  p_cats <- sort(unique(unlist( sapply(p_cats, function(x) strsplit(x, ', ')) )))
  p_cats <- sort(unique(unlist( sapply(p_cats, function(x) strsplit(x, ',')) )))
  p_cats[p_cats == "NULL"] <- "UNKNOWN"
  
  allTables$raw[[p]] <- matrix(0, nrow=ncol(otu_table(SLL)), ncol=length(p_cats))
  allTables$norm[[p]] <- matrix(0, nrow=ncol(otu_table(SLL)), ncol=length(p_cats))
  rownames(allTables$raw[[p]]) <- colnames(otu_table(SLL))
  rownames(allTables$norm[[p]]) <- colnames(otu_table(SLL))
  colnames(allTables$raw[[p]]) <- p_cats
  colnames(allTables$norm[[p]]) <- p_cats
}


acceptable_genera <- rownames(new_phenos_table)[rownames(new_phenos_table) %in% taxa_names(SLL)] #those both in our data and the METAGENassist database

for ( phenotype in colnames(new_phenos_table) ) {
# for ( phenotype in c("MOTILITY", "GRAM") ) {
  print(phenotype)
  for ( p_cat in colnames( allTables$raw[[phenotype]] ) ) {
    
    # first check if p_cat matches exactly, then check if it is part a list (some separated by ';' and some by ', ')
    # cannot do a simple grepl of p_cat since it may be part of other potential p_cats 
    # (e.g. "active" is one of the listed "DISEASES", but so is "active rheumatic fever", so must grep for "active;" or one of these below)
    genera <- acceptable_genera[ new_phenos_table[acceptable_genera, phenotype] == p_cat | 
                                   grepl( sprintf(';%s', p_cat), new_phenos_table[acceptable_genera, phenotype] ) |
                                   grepl( sprintf('%s;', p_cat), new_phenos_table[acceptable_genera, phenotype] ) |
                                   grepl( sprintf(', %s', p_cat), new_phenos_table[acceptable_genera, phenotype] ) |
                                   grepl( sprintf('%s,', p_cat), new_phenos_table[acceptable_genera, phenotype] )]
    
    # if no genera matching the given p_cat, give 0s for that p_cat for all samples, this column will then be removed outside of this for loop
    ifelse( length(genera)>0, {
      allTables$raw[[phenotype]][ , p_cat] <- apply( otu_table(SLL)[genera,], 2, sum )
      allTables$norm[[phenotype]][ , p_cat] <- apply( otu_table(SLL1)[genera,], 2, sum )
    }, {
      allTables$raw[[phenotype]][ , p_cat] <- numeric(length = length(sample_names(SLL)))
      allTables$norm[[phenotype]][ , p_cat] <- numeric(length = length(sample_names(SLL)))
    } )

  }
  # keep only those columns for phenotype categories with non-0 sums
  allTables$raw[[phenotype]] <- as.matrix( allTables$raw[[phenotype]][ , colSums(as.data.frame(allTables$raw[[phenotype]]))>0 ] )
  allTables$norm[[phenotype]] <- as.matrix( allTables$norm[[phenotype]][ , colSums(as.data.frame(allTables$norm[[phenotype]]))>0 ] )
}

# save allTables into a file
saveRDS( allTables, sprintf('%s/Part_1/METAGENassist/allTables.rds', home_dir) )



