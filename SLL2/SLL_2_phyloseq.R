
library(ggplot2)
library(phyloseq)
# source('http://bioconductor.org/biocLite.R')
# biocLite('phyloseq')
library(ape)
library(reshape2)
library(picante)
library(matrixStats)
library(cluster)
library(clusterSim)
library(ade4)
library(Rmisc)
library(scales)
library(factoextra)
library(microbiome)

library(knitr)
library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(vegan)

library(robCompositions)







# ******************************************************* #
#load microbiome data ####
# ******************************************************* #

# depending on where the script is located...
if (dir.exists("/users/tg/jwillis/SLL")) {
  home_dir <- "/users/tg/jwillis/SLL"
} else if (dir.exists("/home/jwillis/gpfs/projects/bsc40/current/jwillis/SLL")) {
  home_dir <- "/home/jwillis/gpfs/projects/bsc40/current/jwillis/SLL"
} else if (dir.exists("/gpfs/projects/bsc40/current/jwillis/SLL")) {
  home_dir <- "/gpfs/projects/bsc40/current/jwillis/SLL"
} else if (dir.exists("/home/jesse/gpfs/projects/bsc40/current/jwillis/SLL")) {
  home_dir <- "/home/jesse/gpfs/projects/bsc40/current/jwillis/SLL"
}else if (dir.exists(sprintf("%s/SLL", getwd()))) {
  home_dir <- sprintf("%s/SLL", getwd())
} else if (dir.exists("~/Downloads/SLL")) {
  home_dir <- "~/Downloads/SLL"
} else {
  home_dir <- getwd()
}

if (dir.exists("/users/tg/jwillis/SLL") || dir.exists("/home/jwillis/gpfs/projects/bsc40/current/jwillis/SLL") ||
    dir.exists("/home/jesse/gpfs/projects/bsc40/current/jwillis/SLL") || dir.exists("/gpfs/projects/bsc40/current/jwillis/SLL") ||
    dir.exists(sprintf("%s/SLL", getwd())) || dir.exists("~/Downloads/SLL")) {
  p2_dir    <- sprintf("%s/Part_2", home_dir)
  moth2_dir  <- sprintf("%s/Part_2/mothur", home_dir)
  dada2_dir <- sprintf("%s/Part_2/DADA2", home_dir)
  surv2_dir  <- sprintf("%s/Part_2/Surveys", home_dir)
  water_dir <- sprintf("%s/Water_quality_data", home_dir)
  esp_dir   <- sprintf("%s/ESP_adm_shp", home_dir)
} else {
  p2_dir <- moth2_dir <- dada2_dir <- surv2_dir <- water_dir <- esp_dir <- home_dir
}








# ***************************************************************************************** ####
# Objects to load ####




gloms     <- readRDS(sprintf("%s/R_objects/gloms.rds", p2_dir))
gloms_rel <- readRDS(sprintf("%s/R_objects/gloms_rel.rds", p2_dir))
gloms_clr <- readRDS(sprintf("%s/R_objects/gloms_clr.rds", p2_dir))
taxTables.both <- readRDS(sprintf("%s/R_objects/taxTables.both.rds", p2_dir))




# **************************************************************** #
# function to print all species of a taxa within a given level
taxa_print <- function(tl, tax, phy=SLL2, tl.limit=NULL) {
  ttp <- phy@tax_table[ phy@tax_table[ , tl ] %in% tax, ]
  
  if ( ! is.null(tl.limit)) {
    ttp <- unique(ttp[ , 1:match(tl.limit, colnames(ttp))])
  }
  
  return(ttp)
}

# function to print mean abundances of all species of a taxa within a given level
taxa_abunds <- function(tl, tax, phy=SLL2_rel) {
  species <- rownames(phy@tax_table[ phy@tax_table[ , tl ] == tax, ])
  rowMeans(phy@otu_table[ species, ])
}
# **************************************************************** #





SLL2      <- readRDS(sprintf("%s/R_objects/SLL2.rds", p2_dir))
SLL2_rel  <- readRDS(sprintf("%s/R_objects/SLL2_rel.rds", p2_dir))
SLL2.meta <- readRDS(sprintf("%s/R_objects/SLL2.meta.rds", p2_dir))
meta.healthy <- readRDS(sprintf("%s/R_objects/meta.healthy.rds", p2_dir))
phy.healthy <- phyloseq(otu_table(prune_samples(rownames(meta.healthy), SLL2)),
                        sample_data(meta.healthy),
                        phy_tree(prune_samples(rownames(meta.healthy), SLL2)),
                        tax_table(prune_samples(rownames(meta.healthy), SLL2)))

# list of the subSample groups for each of the 100 iterations for each variable:
all_subSamps <- readRDS(file = sprintf("%s/R_objects/all_subSamps.rds", p2_dir))


cont_water_data <- c("Conductivity","water_pH","Dry.residue.at.180.C","Dry.residue.at.110.C","Cl","F","HCO3","NO3",
                     "SO4","Na","K","Li","Ca","Mg","Sr","Hardness","Alcalinity")#,"CO3"
group_water_data <- c("Mineralization", "Composition", "Hardness_category")


additional_diseases <- c("Diabetes","Hypertension","Cholesterol","Depression","Anxiety","Headaches",
                         "Lactose_intolerant","Gastritis","Intestinal_issues","Anemia","Sinusitis",
                         "Fibrosis_carrier","Thyroid_issue","Hypothyroidism","Cancer","Transplant",
                         "Immune_issues","Skin_issues","Lung_issues","Circulatory_issues","Kidney_issues",
                         "Tonsil_issues","Samp_collect_issues","Birth_control","Antihistamines")


# # **************** #
groupQs <- c( "Gender","Age_groups",
              # "Family_participants","Grandparent","Parent","Sibling","Grandchild","Child","Partner","Gender",
              # "Country_of_birth","Province_of_birth","City_of_birth","Country_of_birth.mother","Province_of_birt.mother",
              # "City_of_birth.mother","Country_of_birth.father","Province_of_birth.father","City_of_birth.father",
              "Family_unit","Sibling_unit","Twin_unit","Partner_unit","Parent_Child_unit",
              # "Mother_Child_unit","Father_Child_unit",
              "Grandparent_Grandchild_unit",
              "Ethnicity.Caucasian","Ethnicity.Asian","Ethnicity.African","Ethnicity.Arab","Ethnicity.Gypsy",
              "Ethnicity.Native_American","Ethnicity.No_response","Education.mother","Education.father","Education",
              "Occupation.mother","Occupation.father","Municipal_zone",
              "Moisture_in_home","Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds",
              "Pets.Reptiles_amphibians","Pets.Fish","Pets.Type.Small_furry_animals","Pets.Rabbits","Pets.Rodents",
              "Pets.Mammals","Pets.Horses",
              "Smoker","Water_type_home","Braces","Braces.binary","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement",
              "Mouth_wounds","Reason_dental_visit",
              "Chronic_disorder","Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome","Eating_disorder",
              "Other_disorder_binary","Medications","Antibiotics","Analgesics","Vitamin_supplements",
              "Other_medications_binary","Asthma","Wheezing","How_do_you_feel","Do_you_feel_well",additional_diseases,
              "DS_family","CF_family","Celiac_family",
              "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen",
              "Allergy.Animals","Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
              "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary",
              "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
              "Kissing_partner","BMI_group","BMI_official",group_water_data,
              "MALDI.Yeast_detected","MALDI.Mold_detected","MALDI.Bacteria_detected",#fung.specs.Qs, fung.bact.Qs,
              "Full_MALDI.Candida","Full_MALDI.Candida_albicans",
              "Full_MALDI.Candida_dubliniensis","Full_MALDI.Candida_glabrata",
              "Full_MALDI.Candida_parapsilosis","Full_MALDI.Candida_guillermondii",
              "Full_MALDI.Candida_intermedia","Full_MALDI.Candida_krusei",
              "Full_MALDI.Candida_lusitaniae","Full_MALDI.Cryptococcus_spp",
              "Full_MALDI.Debaryomyces_hansenii","Full_MALDI.Rhodotorula_mucilaginosa",
              "Consumption.Milk.Binary","Consumption.Yogurt.Binary","Consumption.Sweets.Binary",
              "Consumption.Chewing_gum.Binary","Consumption.Nuts.Binary","Drinks.Decaf_coffee.Binary",
              "Drinks.Coffee.Binary","Drinks.Tea.Binary","Drinks.Infusion.Binary","Drinks.Soda.Binary",
              "Drinks.Soda_sugarless.Binary","Drinks.Soda_decaf.Binary","Drinks.RedBull.Binary",
              "Drinks.Other_sugary_drinks.Binary","Drinks.Alcohol_cold.Binary","Drinks.Alcohol_hot.Binary",
              "Brushing.Binary","Floss.Binary","Last_dental_visit.Binary","City","Province","Community","seqGroup",
              # "Diversity_group_Div.Shannon","Diversity_group_Div.Simpson","Diversity_group_Weighted_Unifrac",
              # "Diversity_group_Unweighted_Unifrac","Diversity_group_Faiths.PD","Diversity_group_Species_Richness",
              # "Diversity_group_Bray.Curtis","Diversity_group_Canberra",
              "Stomatotype_JSD","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
              "Stomatotype_VAW_GUnifrac","Stomatotype_a0_GUnifrac","Stomatotype_a05_GUnifrac",
              "Stomatotype_Bray","Stomatotype_Jaccard","Stomatotype_Canberra","Stomatotype_Aitchison",
              "Stomatotype_JSD.3",
              "Stomatotype_Weighted_Unifrac.3","Stomatotype_Weighted_Unifrac.4","Stomatotype_VAW_GUnifrac.3")
# "Stomatotype_Weighted_Unifrac.2","Stomatotype_Weighted_Unifrac.4" )


#numerical columns
only_cont <- c( "Weight","Height","BMI","Number_inhabitants","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats",
                "Pets.Number.Small_furry_animals","Pets.Number.Birds","Pets.Number.Reptiles_amphibians",
                "Pets.Number.Fish","Pets.Number.Rabbits","Pets.Number.Rodents","Pets.Number.Mammals","Pets.Number.Horses",
                "Number_smokers_home","Consumption.Milk","Consumption.Yogurt","Consumption.Sweets",
                "Consumption.Chewing_gum","Consumption.Nuts",
                "Drinks.Decaf_coffee","Drinks.Coffee","Drinks.Tea","Drinks.Infusion","Drinks.Soda",
                "Drinks.Soda_sugarless","Drinks.Soda_decaf","Drinks.RedBull","Drinks.Other_sugary_drinks",
                "Drinks.Alcohol_cold","Drinks.Alcohol_hot",
                "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth",
                "Brushing","Floss","Last_dental_visit","pH",cont_water_data,#"FQ.Diet.Isotonic_drinks_frequency",
                "MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies",
                "Div.Observed","Div.Chao1","Div.ACE","Div.Shannon","Div.Simpson","Div.InvSimpson",#"Div.Fisher",
                "Faiths.PD","Species_Richness","Num_OTUs","Gene_counts","JSD","Weighted_Unifrac","Unweighted_Unifrac",
                "Bray.Curtis","Canberra","Age","Latitude","Longitude","Population" )



# ***************************************************************************************** ####










# data.tax <- read.delim(sprintf("%s/mothur_otu_table.csv", moth2_dir, level), row.names = 2, header=T)

# ******************************************************************** #
# Get tax tax and otu tables
taxTables.orig <- readRDS(sprintf("%s/taxTables.rds", dada2_dir))
otus.orig <- readRDS(sprintf("%s/otus.rds", dada2_dir))
otus_rel.orig <- readRDS(sprintf("%s/otus_rel.rds", dada2_dir))

taxTables.added <- readRDS(sprintf("%s/added_samples/taxTables.rds", dada2_dir))
otus.added <- readRDS(sprintf("%s/added_samples/otus.rds", dada2_dir))
otus_rel.added <- readRDS(sprintf("%s/added_samples/otus_rel.rds", dada2_dir))



# ignore the mock community samples which are not useful here
mocks <- c("HM-782D-analysis1","HM-782D-analysis2","HM-782D-analysis3","HM-782D-analysis4","HM-782D-analysis5",
           "HM-782D-analysis6","HM-783D-analysis1","HM-783D-analysis2","HM-783D-analysis3","HM-783D-analysis4",
           "HM-783D-analysis5","HM-783D-analysis6","MZ","NTC1-analysis1","NTC1-analysis2","NTC1-analysis3",
           "NTC1-analysis4","NTC1-analysis5","NTC1-analysis6","NTC2-analysis1","NTC2-analysis2",
           "NTC2-analysis3","NTC2-analysis4","NTC2-analysis5","NTC2-analysis6")

mocks.added <- c("HM-782D.1","HM-783D.1","NTC1.1","HM-782D.2","HM-783D.2","NTC1.2","HM-782D.3","HM-783D.3","NTC1.3")

for (tl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  otus.orig[[ tl ]] <- otus.orig[[ tl ]][ , ! colnames(otus.orig[[ tl ]]) %in% mocks ]
  otus_rel.orig[[ tl ]] <- otus_rel.orig[[ tl ]][ , ! colnames(otus_rel.orig[[ tl ]]) %in% mocks ]
  
  otus.added[[ tl ]] <- otus.added[[ tl ]][ , ! colnames(otus.added[[ tl ]]) %in% mocks.added ]
  otus_rel.added[[ tl ]] <- otus_rel.added[[ tl ]][ , ! colnames(otus_rel.added[[ tl ]]) %in% mocks.added ]
}





# ******************************************************************** #
# get questionnaire table ####
# ************************** #

que <- read.delim(sprintf("%s/SLL_survey_part2_total_repaired.csv", surv2_dir), row.names = 1, stringsAsFactors = F)
# keep only those rows for which we have 16S data
# que <- que[ rownames(que) %in% colnames(otus.orig[["Species"]]), ]
que <- que[ rownames(que) %in% c( colnames(otus.orig[["Species"]]), colnames(otus.added[["Species"]])), ]


# Update the column names so they are understandable
qIDs <- as.matrix( read.delim(sprintf("%s/question_ID_translations.txt", surv2_dir), sep = "=", header = F) )
colnames(que)[ colnames(que) == qIDs[,1] ] <- qIDs[,2]


# ************************** #
# have to fix a few clear typos with family IDs that I discovered
# or some for which the relative gave samples later, so the first person did not put any ID
que["10-05","Parent"] <- 1
que["10-05","Parent_ID"] <- "15-13"

que["10-13","Parent"] <- 1
que["10-13","Parent_ID"] <- "15-12"

que["10-16","Parent"] <- 1
que["10-16","Parent_ID"] <- "15-11"

que["10-18","Parent"] <- 1
que["10-18","Parent_ID"] <- "15-09"

que["10-19","Parent"] <- 1
que["10-19","Parent_ID"] <- "15-10"

que["10-22","Parent"] <- 1
que["10-22","Parent_ID"] <- "15-05"

que["10-52","Parent"] <- 1
que["10-52","Parent_ID"] <- "15-14"

que["10-38","Sibling_ID"] <- "15-01"

que["12-01","Partner"] <- 1
que["12-01","Partner_ID"] <- "12-09"

que["16-10","Parent"] <- 1
que["16-10","Parent_ID"] <- "16-44"

que["27-15","Parent"] <- 1
que["27-15","Parent_ID"] <- "27-83,27-84"

que["30-34","Sibling"] <- 0
que["30-34","Sibling_ID"] <- 0

que["31-15","Parent"] <- 1
que["31-15","Parent_ID"] <- "31-98"

que["31-23","Parent"] <- 1
que["31-23","Parent_ID"] <- "31-72"

que["31-24","Parent"] <- 1
que["31-24","Parent_ID"] <- "31-86,31-87"
que["31-86","Partner"] <- 1
que["31-86","Partner_ID"] <- "31-87"
que["31-87","Partner"] <- 1
que["31-87","Partner_ID"] <- "31-86"

que["31-36","Parent"] <- 1
que["31-36","Parent_ID"] <- "31-81"

que["31-44","Parent"] <- 1
que["31-44","Parent_ID"] <- "31-66"

que["31-45","Parent"] <- 1
que["31-45","Parent_ID"] <- "31-97"
que["31-97","Parent"] <- 0
que["31-97","Parent_ID"] <- 0
que["31-97","Child"] <- 1
que["31-97","Child_ID"] <- "31-45"

que["31-56","Parent"] <- 1
que["31-56","Parent_ID"] <- "31-0A"

que["31-80","Partner_ID"] <- "31-78"

que["35-14","Parent"] <- 1
que["35-14","Parent_ID"] <- "36-05"

que["35-15","Parent"] <- 1
que["35-15","Parent_ID"] <- "36-06"

que["46-12","Child"] <- 1
que["46-12","Child_ID"] <- "46-24"

que["48-06","Sibling"] <- 1
que["48-06","Sibling_ID"] <- "48-01"

que["49-25","Sibling"] <- 1
que["49-25","Sibling_ID"] <- "49-11"

que["52-17","Parent"] <- 0
que["52-17","Parent_ID"] <- 0
que["52-17","Child"] <- 1
que["52-17","Child_ID"] <- "52-52"

# 49-11, 59-17 and 60-28 were apparently not sequenced yet, but their sibling, child and parent (respectively) were

# ************************** #

# add columns for Mother and Mother_ID and for Father and Father_ID

que[ , "Mother" ] <- sapply(rownames(que), function(x) {
  if (que[x,"Parent_ID"]==0) {
    return(0)
  } else {
    pid <- que[x,"Parent_ID"]
    if (pid=="No Sabe/No Contesta")
      return(0)
    
    if (grepl(",",pid)) {
      pid <- strsplit(pid, ",")[[1]]
    }
    
    if (sum(que[pid, "Gender"]==0) > 0)
      return(1)
    else
      return(0)
  }
})

que[ , "Mother_ID" ] <- sapply(rownames(que), function(x) {
  if (que[x,"Parent_ID"]==0) {
    return(0)
  } else {
    pid <- que[x,"Parent_ID"]
    if (pid=="No Sabe/No Contesta")
      return(0)
    
    if (grepl(",",pid)) {
      pid <- strsplit(pid, ",")[[1]]
    }
    
    if (sum(que[pid, "Gender"]==0) > 0)
      return(pid[que[pid, "Gender"]==0])
    else
      return(0)
  }
})

que[ , "Father" ] <- sapply(rownames(que), function(x) {
  if (que[x,"Parent_ID"]==0) {
    return(0)
  } else {
    pid <- que[x,"Parent_ID"]
    if (pid=="No Sabe/No Contesta")
      return(0)
    
    if (grepl(",",pid)) {
      pid <- strsplit(pid, ",")[[1]]
    }
    
    if (sum(que[pid, "Gender"]==1) > 0)
      return(1)
    else
      return(0)
  }
})

que[ , "Father_ID" ] <- sapply(rownames(que), function(x) {
  if (que[x,"Parent_ID"]==0) {
    return(0)
  } else {
    pid <- que[x,"Parent_ID"]
    if (pid=="No Sabe/No Contesta")
      return(0)
    
    if (grepl(",",pid)) {
      pid <- strsplit(pid, ",")[[1]]
    }
    
    if (sum(que[pid, "Gender"]==1) > 0)
      return(pid[que[pid, "Gender"]==1])
    else
      return(0)
  }
})


# ************************** #

# Update format of some columns that should be numerical, and make columns for binary (yes/no) ####
consumption_freqs <- list( "Nunca O Menos De 1 Vez Al Mes" = 0, "1-3 Veces Al Mes" = 1,
                           "1-2 Veces A La Semana" = 2, "3-4 Veces A La Semana" = 3,
                           "5-6 Veces A La Semana" = 4, "1 Vez Al Dia" = 5,
                           "2 Veces Al Dia" = 6, "Mas De 2 Veces Al Dia" = 7,
                           "No Sabe/No Contesta" = NA )
for (col in c("Consumption.Milk","Consumption.Yogurt","Consumption.Sweets","Consumption.Chewing_gum","Consumption.Nuts")) {
  que[ , col ] <- sapply(as.character(que[, col]), function(x) consumption_freqs[[x]])
  que[ , sprintf("%s.Binary", col) ] <- sapply( que[, col], function(x) ifelse(x==0, "No", "Yes"))
}

# For drinks questions
for (col in c("Drinks.Decaf_coffee","Drinks.Coffee","Drinks.Tea","Drinks.Infusion","Drinks.Soda",
              "Drinks.Soda_sugarless","Drinks.Soda_decaf","Drinks.RedBull","Drinks.Other_sugary_drinks",
              "Drinks.Alcohol_cold","Drinks.Alcohol_hot")) {
  que[ , col ] <- sapply(as.character(que[, col]), function(x) ifelse(x=="No Sabe/No Contesta", NA, as.numeric(x)))
  que[ , sprintf("%s.Binary", col) ] <- sapply( que[, col], function(x) ifelse(x==0, "No", "Yes"))
}

# For brushing
brush_freqs <- list( "Nunca" = 0, "1 Vez Al Dia" = 1, "2 Veces Al Dia" = 2, 
                     "Mas De 2 Veces Al Dia" = 3, "No Sabe/No Contesta" = NA)
que[ , "Brushing" ] <- sapply(as.character(que[ , "Brushing" ]), function(x) brush_freqs[[x]])
que[ , "Brushing.Binary" ] <- sapply( que[, "Brushing"], function(x) ifelse(x==0, "No", "Yes"))

# For flossing
floss_freqs <- list( "Nunca O Menos De 1 Vez Al Mes" = 0, "1-3 Veces Al Mes" = 1, "1-3 Veces A La Semana" = 2, 
                     "4-6  Veces A La Semana" = 3, "4-6 Veces A La Semana" = 3, # extra space in some of these
                     "Cada Dia" = 4, "No Sabe/No Contesta" = NA)
que[ , "Floss" ] <- sapply(as.character(que[ , "Floss" ]), function(x) floss_freqs[[x]])
que[ , "Floss.Binary" ] <- sapply( que[, "Floss"], function(x) ifelse(x==0, "No", "Yes"))

# For last_dental_visit
dentist_freqs <- list( "Menos De 6 Meses" = 1, "Entres 6 Meses Y 1 Ano" = 2,
                     "Entre 1 Y 2 Anos" = 3, "Entre 2 Y 3 Anos" = 4,
                     "Entre 3 Y 5 Anos" = 5, "Entres 3 Y 5 Anos" = 5, # extra "s" in some of these
                     "Mas De 5 Anos" = 6, "Nunca" = 7)
que[ , "Last_dental_visit" ] <- sapply(as.character(que[ , "Last_dental_visit" ]), function(x) dentist_freqs[[x]])
que[ , "Last_dental_visit.Binary" ] <- sapply( que[, "Last_dental_visit"], function(x) ifelse(x==7, "No", "Yes"))


# Fix spaces in Reason_dental_visit so they match answers in SLL1
que[ , "Reason_dental_visit"] <- gsub('/ ', '/', as.character(que[ , "Reason_dental_visit"]))



# ************************** #
# change "10 O Mas" to 10 in Number_inhabitants
que$Number_inhabitants <- as.numeric(sapply(que$Number_inhabitants, function(x) 
  ifelse(x=="10 O Mas", 10, ifelse(x=="No Sabe/No Contesta", NA, x))))

# change "5 O Mas" to 5 in Number_smokers_home
que$Number_smokers_home <- as.numeric(sapply(que$Number_smokers_home, function(x) 
  ifelse(x=="5 O Mas", 5, ifelse(x=="No Sabe/No Contesta", NA, x))))

# change pH values which were likely typos to the values they should be, make numeric
que$pH[ que$pH %in% c(7.55, 7.6) ] <- "7.5"
que$pH[ que$pH %in% c(6.6) ] <- "6.5"
que$pH <- as.numeric(que$pH)

# add variables for pH labels
que$pH_label <- sapply(que$pH, function(x) ifelse(is.na(x), NA, 
                                                  ifelse(x < 7, "Acidic", 
                                                         ifelse(x < 7.5, "Neutral", "Alkaline"))))
# and wider ranging labels
que$pH_label.wide <- sapply(que$pH, function(x) ifelse(is.na(x), NA, 
                                                       ifelse(x < 6.5, "Highly_Acidic", 
                                                              ifelse(x < 7, "Acidic", 
                                                                     ifelse(x < 7.5, "Neutral", 
                                                                            ifelse(x < 8, "Highly_Alkaline", "Alkaline"))))))


# remove "No Sabe/No Contesta" from these to make them numeric
for (num in c("Weight","Height","Years_in_home",
              "Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
              "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish",
              "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth")) {
  que[ , num ] <- as.numeric( que[ , num ] )
}

c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
  "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish","Number_smokers_home",
  "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")
# ************************** #


# there are 2 samples for which we have 16S but cannot use:
# 42-14 we cannot use because did not have the consent form I think
# 55-55 seems to be mislabelled since we only had 1-13 for site 55. 
#   - Might be sample 52-55 since we have a questionnaire for it, but apparently no 16S, but cant be sure so wont use it
for (tl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  otus.orig[[ tl ]] <- otus.orig[[ tl ]][ , colnames(otus.orig[[ tl ]]) %in% rownames(que) ]
  otus_rel.orig[[ tl ]] <- otus_rel.orig[[ tl ]][ , colnames(otus_rel.orig[[ tl ]]) %in% rownames(que) ]
  
  otus.added[[ tl ]] <- otus.added[[ tl ]][ , colnames(otus.added[[ tl ]]) %in% rownames(que) ]
  otus_rel.added[[ tl ]] <- otus_rel.added[[ tl ]][ , colnames(otus_rel.added[[ tl ]]) %in% rownames(que) ]
}

# ******************************************************************** #









# # # first get list of genus names in order to prepare the tax tree
# # write( sort(unique( rownames(otus[["Genus"]]) )), file = sprintf("%s/tree/taxa_names_genus.txt", p2_dir) )
# # # also get a list of species names. For those unclassified at species level, get genus name only
# # # sp.ge <- unique(unname(sapply(rownames(otus[["Species"]]), function(x) ifelse(startsWith(strsplit(x,' ')[[1]][2], 'unclassified'), strsplit(x,' ')[[1]][1], paste(strsplit(x,' ')[[1]],collapse = '_')))))
# # sp.ge <- unique(unname(sapply(rownames(otus[["Species"]]), function(x) ifelse(startsWith(strsplit(x,' ')[[1]][2], 'unclassified'), strsplit(x,' ')[[1]][1], x))))
# # write( sort(sp.ge), file = sprintf("%s/tree/taxa_names_species.txt", p2_dir) )
# # # then use the common tree tool from NCBI to get a tree: https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
# # # upload the genus names file, then save tree as a phylip tree
# 
# 
# # get tree object
# # # wont use this tree, since it is based only on genus level
# # ocd.tree <- read.tree(file = "/users/tg/jwillis/SLL/OCD/OCD.tree")
# 
# # read.tree() function has a problem with spaces in the names (as in the <Genus species> values)
# # so first must read the tree file in as text, and put it all into one character object
# sll2.TreeNoSpaces <- paste(as.character(as.matrix(read.delim(sprintf("%s/tree/SLL2.tree_species",p2_dir),header = F))), 
#                            collapse = '')
# # then replace the spaces with "_" (checked to make sure there were no genus names that had any "_" already, 
# # since this would cause a problem when changing back to spaces)
# sll2.TreeNoSpaces <- gsub(' ', '_', sll2.TreeNoSpaces)
# # get tree
# sll2.tree <- read.tree(text = sll2.TreeNoSpaces)
# # now change the "_" back to spaces so that they match the taxa names in all the other objects
# sll2.tree$tip.label <- gsub('_', ' ', sll2.tree$tip.label)
# 
# 
# # # unfortunately any genera with also had unclassified species will be ignored here, so will artificially add tips for these
# # genus.nodes <- sll2.tree$node.label[sll2.tree$node.label %in% rownames(otus[["Genus"]])]
# # # get genera that have an unclassified species present but not in tree
# # all.species <- unique( rownames(otus[["Species"]]) )
# # genus.nodes <- genus.nodes[ ! sapply(genus.nodes, 
# #                                      function(x) all.species[ sapply(all.species, function(y) strsplit(y, ' ')[[1]][1]) == strsplit(x, ' ')[[1]][1] & 
# #                                                                 startsWith(sapply(all.species, function(z) strsplit(z, ' ')[[1]][2]), "unclassified")]) %in% sll2.tree$tip.label ]
# # # add these to the sll2.tree$tip.label
# # sll2.tree$tip.label <- c(sll2.tree$tip.label, genus.nodes)
# 
# 
# # Finally, must change those tip.labels with only genus names, to include the proper unclassified value for species,
# # so that these names will also match the taxa names present in all other objects
# # Because any tip.labels that do not match taxa names in otu_table, sample_data, and tax_table will be 
# # ignored when creating phyloseq object
# all.species <- unique( rownames(otus[["Species"]]) )
# 
# sll2.tree$tip.label <- unname(sapply(sll2.tree$tip.label, function(x) {
#   ifelse(length(strsplit(x, ' ')[[1]]) == 1,
#          all.species[ sapply(all.species, function(y) strsplit(y, ' ')[[1]][1]) == strsplit(x, ' ')[[1]][1] & 
#                         startsWith(sapply(all.species, function(z) strsplit(z, ' ')[[1]][2]), "unclassified") ],
#          x)
# }))
# 
# # ************************ #
# sll2_genus.tree <- read.tree(file = sprintf("%s/tree/SLL2.tree_genus", p2_dir))
# 
# # NCBI automatically changed these, so will change them back to the names that appear in our dataset
# sll2_genus.tree$tip.label[sll2_genus.tree$tip.label=="Floricoccus"] <- "Anthococcus"
# 
# # ************************ #


# ******************************************************************** #






# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####



# ******************************************************* #
# Phylogenetic tree ####
# ******************************************************* #

# Generated during the dada2 pipeline
# fitGTR <- readRDS(sprintf("%s/R_objects/SLL2.fitGTR.rds", p2_dir))
# fitGTR.added <- readRDS(sprintf("%s/R_objects/SLL2_added_samples.fitGTR.rds", p2_dir))

fitGTR <- readRDS(sprintf("%s/R_objects/SLL2_all.fitGTR.rds", p2_dir))

# ******************************************************************** #


# Create phyloseq objects ####
SLL2 <- phyloseq(otu_table(otus.orig[["Species"]], taxa_are_rows = T), 
                 sample_data(que), 
                 tax_table(taxTables.orig[["Species"]]))#,
# phy_tree(fitGTR$tree))


# remove samples with under 1000 total counts
SLL2 <- prune_samples( sample_sums(SLL2) >= 1000, SLL2)




# Create phyloseq object for added samples
SLL2.added <- phyloseq(otu_table(otus.added[["Species"]], taxa_are_rows = T), 
                       sample_data(que), 
                       tax_table(taxTables.added[["Species"]]))#,
# phy_tree(fitGTR.added$tree))


# remove samples with under 1000 total counts
SLL2.added <- prune_samples( sample_sums(SLL2.added) >= 1000, SLL2.added)



# merge these objects to make one phyloseq object
SLL2 <- merge_phyloseq(SLL2, SLL2.added)
# then add the tree used for both datasets
SLL2@phy_tree <- phy_tree(fitGTR$tree)



# change taxa and sample names to replace '-' with '.' because many functions seem to prefer it this way
taxa_names(SLL2) <- gsub('-', '\\.', taxa_names(SLL2))

sample_names(SLL2) <- gsub('-', '\\.', sample_names(SLL2))
sample_names(SLL2) <- paste("SLL", sample_names(SLL2), sep = '.')






# ******************************************************************** #


# #### check number of species in each data set originally found by DADA2, and after filtering
# os <- length(rownames(otus[["Species"]])) # 425 - originally
# ds <- length(unique(rownames(otus[["Species"]])[ ! grepl("unclassified", rownames(otus[["Species"]])) ])) # 156 - known at Species level from DADA2
# fs <- length(unique(SLL2@tax_table[,"Species"][ ! grepl("unclassified", SLL2@tax_table[,"Species"]) ])) # 66 - after filtering
# 
# og <- length(rownames(gloms[["Genus"]])) # 172 - originally
# dg <- length(unique(rownames(gloms[["Genus"]])[ ! grepl("unclassified", rownames(gloms[["Genus"]])) ])) # 127 - known at Genus level from DADA2
# fg <- length(unique(SLL2@tax_table[,"Genus"][ ! grepl("unclassified", SLL2@tax_table[,"Genus"]) ])) # 34 - after filtering
# 
# check_nums <- as.data.frame(matrix(c(os,ds,fs, og,dg,fg), nrow=3))
# rownames(check_nums) <- c("Unique", "DADA2", "Filtered")
# colnames(check_nums) <- c("Species", "Genus")
# 
# check_nums.m <- melt(check_nums)
# categ <- factor(rep(c("Unique", "DADA2", "Filtered"), 2))
# check_nums.m[,"category"] <- factor(categ, levels = c("Unique", "DADA2", "Filtered"))
# 
# ggplot(check_nums.m, aes(x = variable, y = value)) +
#   geom_bar(aes(fill = category), position = "dodge", stat="identity") +
#   geom_vline(xintercept = 1.5) +
#   ggtitle("Remove taxa not seen in at least 25 samples") +
#   theme(axis.title = element_text(size=14), axis.text = element_text(size = 14),
#         legend.text = element_text(size=14), legend.title = element_text(size=14)) +
#   geom_text(aes(x=c(0.66,1,1.33, 1.66,2,2.33),
#                 c(471,231,145, 250,186,106)),
#             label=c(os,ds,fs, og,dg,fg))



# ******************************************************************** #

# # now do same filtering for the otu tables at all levels, since these will also be used throughout the analyses
# for (tl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
#   
#   keep.tl <- sapply(rownames(otus[[ tl ]]), function(x) sum(otus[[ tl ]][x,] > 0) > 25)
#   keep.tl <- names(keep.tl[ keep.tl ])
#   otus[[ tl ]] <- otus[[ tl ]][ keep.tl, ]
#   otus_rel[[ tl ]] <- otus_rel[[ tl ]][ keep.tl, ]
#   taxTables[[ tl ]] <- taxTables[[ tl ]][ keep.tl, ]
# }



# ******************************************************************** #

# Add school ID to sample_data
SLL2@sam_data[, "School_ID"] <- sapply( sample_names(SLL2), function(x) as.integer( strsplit(x,'\\.')[[1]][2] ))
SLL2@sam_data[, "School_ID.letter"] <- sapply( sample_names(SLL2), function(x) 
  sprintf("S.%s", as.integer( strsplit(x,'\\.')[[1]][2] )) )

# *********** #
# Plot boxes of diversity values in each school to show batch effect in first 5 schools ####
SLL_test <- SLL2
rich_all <- estimate_richness(prune_taxa(taxa_sums(SLL_test)>0, SLL_test),
                              measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))

measure <- "Shannon"
# measure <- "Simpson"
# measure <- "Veillonella"

# div.by.school <- as.data.frame( cbind(t(tagl@otu_table)[ , measure ], tagl@sam_data[,"School_ID"]) )
div.by.school <- as.data.frame( cbind(rich_all[,measure], sample_data(SLL_test)[,"School_ID"]) )
# div.by.school <- as.data.frame( cbind(SLL_test@sam_data[,"Unweighted_Unifrac"], sample_data(SLL_test)[,"School_ID"]) )
colnames(div.by.school) <- c("value","School_ID")

counts <- function(x) {
  ifelse(measure=="Shannon", y_add <<- 0.55, y_add <<- 0.12)
  return( data.frame(y=median(x) + y_add, label=length(x)) )
}

ggplot(div.by.school, aes(x=reorder(School_ID,-value,FUN=median), y=value, 
                          fill=reorder(School_ID,-value,FUN=median))) +
  geom_boxplot() + theme_minimal() +
  # theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90), axis.title=element_text(size=16, face="bold")) +
  # ggtitle(sprintf('%s per sample for %s', measure, "School_ID")) +
  xlab("School_ID") + ylab(measure) + guides(fill=FALSE) +#scale_fill_hue(name='School_ID') +#
  stat_summary(fun.data = counts, geom='text', size=4, angle=90)
# *********** #






# ******************************************************************** #

# include columns for the various diversity measures ####
rich <- estimate_richness(prune_taxa(taxa_sums(SLL2)>0, SLL2),
                          measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))#, measures='Simpson')
rich[ is.nan(as.matrix(rich)) ] <- NA

SLL2@sam_data[,'Div.Observed'] <- rich$Observed
SLL2@sam_data[,'Div.Chao1'] <- rich$Chao1
SLL2@sam_data[,'Div.ACE'] <- rich$ACE
SLL2@sam_data[,'Div.Shannon'] <- rich$Shannon
SLL2@sam_data[,'Div.Simpson'] <- rich$Simpson
SLL2@sam_data[,'Div.InvSimpson'] <- rich$InvSimpson
# SLL2@sam_data[,'Div.Fisher'] <- rich$Fisher


# include column for number of OTUs identified
SLL2@sam_data[,'Num_OTUs'] <- apply(SLL2@otu_table, 2, function(x) sum( x != 0))


# include column for number of 16S gene counts per sample
SLL2@sam_data[,'Gene_counts'] <- colSums(SLL2@otu_table)




# include column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
faith <- pd(t(as.data.frame(SLL2@otu_table)), SLL2@phy_tree, include.root = F)
faith$PD[ is.na(faith$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
SLL2@sam_data[,'Faiths.PD'] <- faith$PD
SLL2@sam_data[,'Species_Richness'] <- faith$SR



# ********************************************** #
# Check distributions of a div measures ####
#   get info from this tutorial: https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html
par(mfrow = c(3, 2))

#Then plot each metric.
hist(SLL2.meta$Div.Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(SLL2.meta$Faiths.PD, main="Faiths PD", xlab="", breaks=15)
hist(SLL2.meta$Div.Simpson, main="Simpson diversity", xlab="", breaks=10)
hist(SLL2.meta$Div.InvSimpson, main="Inverse Simpson", xlab="", breaks=15)
hist(SLL2.meta$Div.Chao1, main="Chao richness", xlab="", breaks=15)
hist(SLL2.meta$Species_Richness, main="Species Richness", xlab="", breaks=15)

par(mfrow = c(1, 1))

# test for normalcy of div measures with Shapiro-Wilk test of normality
shapiro.test(SLL2.meta$Div.Shannon)
shapiro.test(SLL2.meta$Faiths.PD)
shapiro.test(SLL2.meta$Div.Simpson)
shapiro.test(SLL2.meta$Div.InvSimpson)
shapiro.test(SLL2.meta$Div.Chao1)
shapiro.test(SLL2.meta$Species_Richness)

# this suggests that none of these measures is normally distributed, should use non-parametric tests later.
# although the Shapiro-Wilk test is quite sensitive, and some of these measures are at least roughly normal

# ********************************************** #





# ******************************************************************** #
# remove taxa not seen at least once in at least 25 samples (2.3% of the 1085 samples).
# This helps protect against an OTU with small mean & trivially large Coefficient of Variation.
# from /users/tg/jwillis/SLL/Papers/Tools/phyloseq-article-source-files-figs-code-03/phyloseq_plos1_2012-source-doc.html
SLL2 <- filter_taxa(SLL2, function(x) sum(x > 0) > 25, prune = TRUE)
# SLL2 <- filter_taxa(SLL2, function(x) sum(x > 0) > (0.05 * length(x)), prune = TRUE)


# ******************************************************************** #










# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####







# ******************************************************************** #
# give roots to trees that are unrooted

library(QsRutils)

SLL2 <- root_phyloseq_tree(SLL2)

# ******************************************************************** #




# ******************************************************************** #
## Calculate Unifrac distances ####

# # if necessary to calculate these values again, simply uncomment this section.
# # But takes a long time to calculate, better to read from files that were written from original calculation
library(foreach)
library(doParallel)
library(GUniFrac)
# 
# # This is required in order to register a parallel "backend" so the calculation can be run in parallel
# registerDoParallel(cores = 10)
# 
# weighted_Unifrac <- UniFrac(SLL2, weighted = T, parallel = T)
# unweighted_Unifrac <- UniFrac(SLL2, weighted = F, parallel = T)
# saveRDS(as.data.frame(as.matrix(weighted_Unifrac)), file = sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))
# saveRDS(as.data.frame(as.matrix(unweighted_Unifrac)), file = sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))
# 
# guni <- GUniFrac(t(SLL2@otu_table), SLL2@phy_tree)
# guni.VAW <- guni$unifracs[ , , "d_VAW"]
# guni.a0  <- guni$unifracs[ , , "d_0"]
# guni.a05 <- guni$unifracs[ , , "d_0.5"]
# saveRDS(as.data.frame(guni.VAW), file = sprintf("%s/R_objects/beta_diversities/SLL2_VAW_Gunifrac.rds", p2_dir))
# saveRDS(as.data.frame(guni.a0), file = sprintf("%s/R_objects/beta_diversities/SLL2_a0_Gunifrac.rds", p2_dir))
# saveRDS(as.data.frame(guni.a05), file = sprintf("%s/R_objects/beta_diversities/SLL2_a05_Gunifrac.rds", p2_dir))


weighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))
diag(weighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
unweighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))
diag(unweighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
guni.VAW <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_VAW_Gunifrac.rds", p2_dir))
diag(guni.VAW) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
guni.a0 <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_a0_Gunifrac.rds", p2_dir))
diag(guni.a0) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
guni.a05 <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_a05_Gunifrac.rds", p2_dir))
diag(guni.a05) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored



SLL2@sam_data[ rownames(weighted_Unifrac), "Weighted_Unifrac"] <- rowMeans(weighted_Unifrac, na.rm = T)
SLL2@sam_data[ rownames(unweighted_Unifrac), "Unweighted_Unifrac"] <- rowMeans(unweighted_Unifrac, na.rm = T)
SLL2@sam_data[ rownames(guni.VAW), "VAW_GUnifrac"] <- rowMeans(guni.VAW, na.rm = T)
SLL2@sam_data[ rownames(guni.a0), "a0_GUnifrac"] <- rowMeans(guni.a0, na.rm = T)
SLL2@sam_data[ rownames(guni.a05), "a05_GUnifrac"] <- rowMeans(guni.a05, na.rm = T)
# ******************************************************************** #




# ******************************************************************** #
## Calculate JSD distances ####

# # if necessary to calculate these values again, simply uncomment this section.
# # But takes a long time to calculate, better to read from files that were written from original calculation
# library(foreach)
# library(doParallel)
# 
# # This is required in order to register a parallel "backend" so the calculation can be run in parallel
# registerDoParallel(cores = 10)
# 
# jsd <- distance(SLL2, method = "jsd", parallel = T)
# saveRDS(as.data.frame(as.matrix(jsd)), file = sprintf("%s/R_objects/beta_diversities/SLL2_jsd.rds", p2_dir))

# *********************** #
jsd <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_jsd.rds", p2_dir))
diag(jsd) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
# jsd <- as.dist(jsd)

SLL2@sam_data[ rownames(jsd), "JSD"] <- rowMeans(jsd, na.rm = T)
# ******************************************************************** #






# ******************************************************************** #
# Bray-Curtis and Canberra distances ####

# Bray-Curtis distance combines absolute differences between features (more sensitive to a few large changes)
bray <- as.matrix( vegdist(t(SLL2@otu_table), method = "bray") )
diag(bray) <- NA

# Jaccard is similar to Bray-curtis, but is metric, where BC is semimetric
jaccard <- as.matrix( vegdist(decostand(as.data.frame(t(SLL2@otu_table)), method="pa"), method = "jaccard") )
diag(jaccard) <- NA

# Canberra distance weights all differences equally (more sensitive to many small changes)
canberra <- as.matrix( vegdist(t(SLL2@otu_table), method = "canberra") )
diag(canberra) <- NA

SLL2@sam_data[ rownames(bray), "Bray.Curtis"] <- rowMeans(bray, na.rm = T)
SLL2@sam_data[ rownames(jaccard), "Jaccard"] <- rowMeans(jaccard, na.rm = T)
SLL2@sam_data[ rownames(canberra), "Canberra"] <- rowMeans(canberra, na.rm = T)


# saveRDS(bray, file = sprintf("%s/R_objects/beta_diversities/SLL2_bray.rds", p2_dir))
# saveRDS(jaccard, file = sprintf("%s/R_objects/beta_diversities/SLL2_jaccard.rds", p2_dir))
# saveRDS(canberra, file = sprintf("%s/R_objects/beta_diversities/SLL2_canberra.rds", p2_dir))


# ******************************************************************** #








# ******************************************************************** #
# Compositionally appropriate calculations: ####
# R_Block_1
# filter the dataset
f <- codaSeq.filter(SLL2@otu_table, 
                    min.reads=1000, # filter out samples with fewer reads
                    min.prop=0.001, # filter out taxa with lower minimum abundance
                    min.occurrence=0.05, # filter out taxa not appearing in at least this prop of samps
                    samples.by.row=FALSE)
# replace 0 values with an estimate
f.n0 <- cmultRepl(t(f), method="CZM", label=0)




# generate the CLR values for plotting later
f.clr <- codaSeq.clr(f.n0)



# ******************************************************************** #
# Aitchison distance ####
aitch <- as.matrix( aDist(f.n0) )
diag(aitch) <- NA

# saveRDS(aitch, file = sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))


SLL2@sam_data[ rownames(aitch), "Aitchison"] <- rowMeans(aitch, na.rm = T)
# ******************************************************************** #






# then get phyloseq object with centered log ratio values
SLL2_clr <- phyloseq(otu_table(t(f.clr), taxa_are_rows = T),
                       sample_data(SLL2),
                       phy_tree(SLL2),
                       tax_table(SLL2))
# ******************************************************************** #









# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####








# ********************************************************************************************************* #
###### Adjust values of some columns ######
# ********************************************************************************************************* #
SLL2@sam_data$Gender <- unname(sapply(as.character(SLL2@sam_data$Gender), function(x) ifelse(x==1, "M", ifelse(x==0, "F", x))))


yes.nos <- c("Family_participants","Grandparent","Parent","Sibling","Grandchild","Child","Partner",
             "Ethnicity.Caucasian","Ethnicity.Asian","Ethnicity.African","Ethnicity.Arab","Ethnicity.Gypsy",
             "Ethnicity.Native_American","Ethnicity.No_response",
             "Moisture_in_home","Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds",
             "Pets.Reptiles_amphibians","Pets.Fish",
             "Smoker","Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement",
             "Chronic_disorder","Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome",
             "Eating_disorder","Other_disorder_binary","Other_disorder",
             "Medications","Antibiotics","Analgesics","Vitamin_supplements","Other_medications_binary",
             "Asthma","Wheezing",
             "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen","Allergy.Animals",
             "Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
             "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary","Allergy.other",
             "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
             "Kissing_partner")

for (col in yes.nos) {
  SLL2@sam_data[,col] <- unname(sapply(as.character(as.matrix(SLL2@sam_data)[,col]), 
                                        function(x) ifelse(x==1, "Yes", ifelse(x==0, "No", x))))
}

# add column for simple yes/no for Braces
SLL2@sam_data[,"Braces.binary"] <- sapply( as.character(as.matrix(SLL2@sam_data)[, "Braces"]),
                                           function(x) ifelse(x=="No", "No", 
                                                              ifelse(x=="No Sabe/No Contesta", NA,
                                                                     "Yes")) )
# and for Mouth_wounds
SLL2@sam_data[,"Mouth_wounds.binary"] <- sapply( as.character(as.matrix(SLL2@sam_data)[, "Mouth_wounds"]),
                                                 function(x) ifelse(x=="0", "No", 
                                                                    ifelse(x=="No Sabe/No Contesta", NA,
                                                                           "Yes")) )

# Add a column also for How_do_you_feel in a binary format: Do_you_feel_well
SLL2@sam_data[,"Do_you_feel_well"] <- sapply( as.character(as.matrix(SLL2@sam_data)[, "How_do_you_feel"]),
                                              function(x) ifelse(x=="Bien/ Normal", "Yes", 
                                                                 ifelse(x=="No Sabe/No Contesta", NA,
                                                                        "No")) )



# ****************************************************************************************************************** #
# add column for ages ####

ages <- as.numeric(sapply(as.matrix(SLL2@sam_data)[, "Birthdate"], function(x) 
  ifelse(nchar(x)==10,
         ifelse( is.na(as.Date(x, "%d-%m-%Y")),
                 ifelse( is.na(as.Date(strsplit(x,'-')[[1]][3], "%Y")),
                         NA,
                         as.numeric( difftime("2017-03-07", as.Date(strsplit(x,'-')[[1]][3], "%Y"), units = "days") / 365)),
                 as.numeric( difftime("2017-03-07", as.Date(x, "%d-%m-%Y"), units = "days") / 365)),
         NA)))
# give those with ages less than 1 year as NA because these were simply typing errors by the participants that put 2014 or 2015 as their birthday
ages[ is.na(ages) ] <- "No Sabe/No Contesta"

SLL2@sam_data[ , "Age"] <- as.numeric(ages)


# column for ages groups
SLL2@sam_data[ , "Age_groups"] <- sapply(SLL2@sam_data$Age, function(x) ifelse(x < 13, "Child",
                                                                               ifelse(x < 20, "Teen", 
                                                                                      ifelse(x < 60, "Adult", "Senior"))))
# make a factor to order the groups logically
SLL2@sam_data$Age_groups <- factor(SLL2@sam_data$Age_groups, levels = c("Child","Teen","Adult","Senior"))



# ****************************************************************************************************************** #
# add columns for various family units ####
# ************************************ #
get_famIDs <- function(samps, sam_table, relations, gender=NULL, twins=FALSE) {
  
  # *************** #
  get_fi <- function(x, sam_table, relations, gender, twins) {
    fi <- as.matrix(sam_table[x, relations])
    
    # # get specified gender if indicated
    # if (! is.null(gender)) {
    #   if (gender=="Mother") {
    #     
    #   }
    # }
    
    # get only twins if indicated
    if (twins) {
      fi <- ifelse(as.matrix(sam_table[x, "Sibling"]) %in% c("Gemelo","Melliza"), fi, 0)
    }
    
    # ignore those with no response or "No Sabe/No Contesta"
    fi <- as.character(fi[ ! fi %in% c(0, "No Sabe/No Contesta") ])
    
    # separate those with multiple IDs in a given category to treat them as individuals
    # print(fi)
    fi <- as.character(unlist(sapply(fi, function(x) unlist(strsplit(x,",")))))
    # change format to match sample names
    fi <- unname(unlist(sapply(fi, function(x) paste("SLL", gsub("-","\\.",x), sep="."))))
    # include self in vector so that itll be included later in calculations
    fi <- sort(c(fi, x))
    # remove any samples that are not in rownames (these were samples that have not been sequenced)
    fi <- fi[ fi %in% rownames(sam_table) ]
  }
  # *************** #
  famIDs <- sapply(samps, get_fi, sam_table, relations, gender, twins)
  # *************** #
  
  # keep only those with values
  famIDs <- famIDs[ unlist(lapply(famIDs, function(x) length(x) > 1)) ]
  
  
  # character vector of all those samples involved in a family unit
  all.fam_members <- sort(unique(unname(unlist(famIDs))))
  
  # make list of unique family units, including all children, parents, siblings, partners, etc that are connected
  famUnits <- unique(sapply(all.fam_members, function(af) 
    list(sort(unique(unname(unlist(famIDs[sapply(names(famIDs), function(x) 
      af %in% famIDs[[x]])])))))))
  # give arbitrary names to each unit
  names(famUnits) <- paste0("Family_", 1:length(famUnits))
  
  # character vector of family unit numbers
  famUnits_vector <- sapply(samps, function(samp) 
    ifelse(samp %in% all.fam_members,
           names(famUnits)[ unlist(lapply(famUnits, function(un) samp %in% un)) ],
           "None"))
  
  return(famUnits_vector)
}
# ************************************ #

# add column to phyloseq object
SLL2@sam_data[ , "Family_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2),
                                               c("Grandparent_ID","Parent_ID","Sibling_ID","Grandchild_ID",
                                                 "Child_ID","Partner_ID"))
SLL2@sam_data[ , "Sibling_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), "Sibling_ID")
SLL2@sam_data[ , "Twin_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), "Sibling_ID", twins = TRUE)
SLL2@sam_data[ , "Partner_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), "Partner_ID")
SLL2@sam_data[ , "Parent_Child_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), c("Parent_ID","Child_ID"))
# SLL2@sam_data[ , "Mother_Child_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), c("Mother_ID","Child_ID"))
# SLL2@sam_data[ , "Father_Child_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), c("Father_ID","Child_ID"))
SLL2@sam_data[ , "Grandparent_Grandchild_unit" ] <- get_famIDs(sample_names(SLL2), meta(SLL2), 
                                                               c("Grandparent_ID","Grandchild_ID"))

**** Mother_Child_unit and Father_Child_unit are messed up because a person can be both a parent and a child 
in those cases where there is a grandparent unit ****




# ****************************************************************************************************************** #

# first have to make data.frame version of sam_data
SLL2.meta <- meta(SLL2)

for (disorder in c("Downs_Syndrome","Cystic_fibrosis","Celiac")) {
  
  # get disorder samples
  disSamps <- rownames( SLL2.meta[ SLL2.meta[ , disorder ] == "Yes", ] )
  
  # get family member samples
  dis.famUnits <- unique( SLL2.meta[disSamps, "Family_unit"] )
  dis.famUnits <- dis.famUnits[ dis.famUnits != "None" ]
  disFams <- rownames( SLL2.meta[ SLL2.meta[, "Family_unit"] %in% dis.famUnits, ])
  
  disFamVar <- sprintf("%s_family", ifelse(disorder=="Downs_Syndrome", "DS", 
                                           ifelse(disorder=="Cystic_fibrosis","CF",
                                                  "Celiac")))
  # have to make a variable for being in a DS/CF family
  SLL2@sam_data[ , disFamVar ] <- "No"
  SLL2@sam_data[ disFams, disFamVar ] <- "Yes"
  SLL2.meta <- meta(SLL2)
  
}
# ****************************************************************************************************************** #




  




# ****************************************************************************************************************** #
# add columns for BMI, BMI classifications ####
# SLL2@sam_data[ SLL2@sam_data[,"Weight"]=="No Sabe/No Contesta" ,"Weight"] <- NA
# SLL2@sam_data[ SLL2@sam_data[,"Height"]=="No Sabe/No Contesta" ,"Height"] <- NA

SLL2@sam_data[ , "BMI" ] <- sapply(sample_names(SLL2), function(x)
  as.numeric(as.matrix(SLL2@sam_data[x,"Weight"])) / ( as.numeric(as.matrix(SLL2@sam_data[x,"Height"])) / 100)**2 )


# get samples by BMI quartile
bmi.quarts <- quantile(SLL2@sam_data$BMI, probs=seq(0,1,0.25), na.rm = T)

q1.samples <- rownames(SLL2@sam_data[ (SLL2@sam_data$BMI < bmi.quarts['25%']) & (! is.na(SLL2@sam_data$BMI)), ] )
q2.samples <- rownames(SLL2@sam_data[ (bmi.quarts['25%'] <= SLL2@sam_data$BMI) & (SLL2@sam_data$BMI < bmi.quarts['50%']) & (! is.na(SLL2@sam_data$BMI)), ] )
q3.samples <- rownames(SLL2@sam_data[ (bmi.quarts['50%'] <= SLL2@sam_data$BMI) & (SLL2@sam_data$BMI < bmi.quarts['75%']) & (! is.na(SLL2@sam_data$BMI)), ] )
q4.samples <- rownames(SLL2@sam_data[ (SLL2@sam_data$BMI >= bmi.quarts['75%']) & (! is.na(SLL2@sam_data$BMI)), ] )

# add BMI group to sample_data
SLL2@sam_data[q1.samples, "BMI_group"]                           <- "Low"
SLL2@sam_data[c(q2.samples,q3.samples), "BMI_group"]             <- "Average"
SLL2@sam_data[q4.samples, "BMI_group"]                           <- "High"
SLL2@sam_data[ is.na(SLL2@sam_data[, "BMI_group"]), "BMI_group"] <- NA

# get samples by official BMI category
underweight <- rownames(SLL2@sam_data[ (SLL2@sam_data$BMI < 18.5) & (! is.na(SLL2@sam_data$BMI) ), ] )
normal      <- rownames(SLL2@sam_data[ (18.5 <= SLL2@sam_data$BMI) & (SLL2@sam_data$BMI < 25) & (! is.na(SLL2@sam_data$BMI) ), ] )
overweight  <- rownames(SLL2@sam_data[ (25.0 <= SLL2@sam_data$BMI) & (SLL2@sam_data$BMI < 30) & (! is.na(SLL2@sam_data$BMI) ), ] )
obese       <- rownames(SLL2@sam_data[ (30.0 <= SLL2@sam_data$BMI) & (SLL2@sam_data$BMI < 50) & (! is.na(SLL2@sam_data$BMI) ), ] )
unknown     <- rownames(SLL2@sam_data[  is.na(SLL2@sam_data$BMI) , ] )

# add official BMI indicators to sample_data
SLL2@sam_data[underweight, "BMI_official"] <- "Underweight"
SLL2@sam_data[normal, "BMI_official"]      <- "Normal"
SLL2@sam_data[overweight, "BMI_official"]  <- "Overweight"
SLL2@sam_data[obese, "BMI_official"]       <- "Obese"
SLL2@sam_data[unknown, "BMI_official"]     <- NA


# histogram of BMIs
hist( as.numeric(as.matrix(SLL2@sam_data[,"BMI"])),
      main = sprintf("Histogram of BMI"),
      xlab = "BMI", col = "darkgoldenrod")
#







# ****************************************************************************************************************** #

# Locations by school -----------------------------------------------------

# ******************* #
# list of cities with the school IDs from there

comunidades.sll2 <- list( "Andalucía" = c(22, 23,24,25,26,29, 27,28),
                          "Aragón" = c(47, 48,50),
                          "Cantabria" = c(46),
                          "Cataluña" = c(1,2,54,55,58,59, 3, 51, 53, 52, 49, 60),
                          "Galicia" = c(35, 36, 37, 38,39, 40,41),
                          "Islas Baleares" = c(4,5,6,7,8,9),
                          "Comunidad de Madrid" = c(30,33,57, 31, 32, 34),
                          "Murcia" = c(17, 18,19,21, 20),
                          "País Vasco" = c(42, 43,44, 45),
                          "Comunidad Valenciana" = c(11,12, 16, 10,14,15, 13))

provinces.sll2 <- list ( "Barcelona" = c(1,2,54,55,58,59, 3, 51, 53), "Tarragona" = c(52), "Girona" = c(49), "Lleida" = c(60),
                         "Baleares" = c(4,5,6,7,8,9),
                         "Valencia" = c(11,12, 16), "Castellón" = c(10,14,15, 13),
                         "Murcia" = c(17, 18,19,21, 20),
                         "Málaga" = c(22, 23,24,25,26,29), "Sevilla" = c(27,28),
                         "Madrid" = c(30,33,57, 31, 32), "Toledo" = c(34),
                         "A Coruña" = c(35, 36, 37), "Pontevedra" = c(38,39, 40,41),
                         "Vizcaya" = c(42, 43,44), "Guipúzcoa" = c(45),
                         "Cantabria" = c(46),
                         "Zaragoza" = c(47, 48,50))

cities.sll2 <- list( "Barcelona" = c(1,2,54,55,58,59), "Roda de Ter" = c(3), "Sant Feliu de Llobregat" = c(51), "San Cugat del Vallés" = c(53), "Amposta" = c(52), "Figueres" = c(49), "Lleida" = c(60),
                     "Palma de Mallorca" = c(4,5,6,7,8,9),
                     "L'Alcora" = c(10,14,15), "Nules" = c(13), "Valencia" = c(11,12), "Benetússer" = c(16),
                     "La Paca" = c(17), "Murcia" = c(18,19,21), "Bullas" = c(20),
                     "Alhaurín de la Torre" = c(22), "Málaga" = c(23,24,25,26,29), "Sevilla" = c(27,28),
                     "Madrid" = c(30,33,57), "Colmenar Viejo" = c(31), "Moralzarzal" = c(32), "Yepes" = c(34),
                     "Os Campons" = c(35), "Cambre" = c(36), "Ponteceso" = c(37), "Ponteareas" = c(38,39), "Vigo" = c(40,41),
                     "Santurtzi" = c(42), "Bilbao" = c(43,44), "Villarreal de Urrechua" = c(45), 
                     "Santander" = c(46),
                     "Tauste" = c(47), "Zaragoza" = c(48,50))

lat_lon.sll2 <- list( "41.39 N, 2.17 E" = c(1,2,54,55,58,59), "41.98 N, 2.31 E" = c(3), "41.40 N, 2.05 E" = c(51), "41.45 N, 2.08 E" = c(53), "40.71 N, 0.58 E" = c(52), "42.27 N, 2.96 E" = c(49), "41.62 N, 0.62 E" = c(60),
                      "39.57 N, 2.65 E" = c(4,5,6,7,8,9),
                      "40.07 N, 0.21 W" = c(10,14,15), "39.85 N, 0.16 W" = c(13), "39.47 N, 0.38 W" = c(11,12), "39.42 N, 0.40 W" = c(16),
                      "37.63 N, 1.97 W" = c(17), "37.99 N, 1.13 W" = c(18,19,21), "38.04 N, 1.67 W" = c(20),
                      "36.66 N, 4.56 W" = c(22), "36.72 N, 4.42 W" = c(23,24,25,26,29), "37.39 N, 5.98 W" = c(27,28),
                      "40.42 N, 3.70 W" = c(30,33,57), "40.66 N, 3.77 W" = c(31), "40.68 N, 3.97 W" = c(32), "39.90 N, 3.63 W" = c(34),
                      "43.28 N, 8.36 W" = c(35), "43.29 N, 8.32 W" = c(36), "43.24 N, 8.90 W" = c(37), "42.18 N, 8.51 W" = c(38,39), "42.24 N, 8.72 W" = c(40,41),
                      "43.33 N, 3.03 W" = c(42), "43.26 N, 2.94 W" = c(43,44), "43.09 N, 2.32 W" = c(45),
                      "43.46 N, 3.81 W" = c(46),
                      "41.92 N, 1.26 W" = c(47), "41.65 N, 0.89 W" = c(48,50))


# Add city, province and community to columns in sample_data
for (city in names(cities.sll2)) {
  IDs <- cities.sll2[[city]]
  sams <- sample_names( subset_samples(SLL2, School_ID %in% IDs) )
  SLL2@sam_data[ sams, "City" ] <- city
}

for (prov in names(provinces.sll2)) {
  IDs <- provinces.sll2[[prov]]
  sams <- sample_names( subset_samples(SLL2, School_ID %in% IDs) )
  SLL2@sam_data[ sams, "Province" ] <- prov
}

for (comm in names(comunidades.sll2)) {
  IDs <- comunidades.sll2[[comm]]
  sams <- sample_names( subset_samples(SLL2, School_ID %in% IDs) )
  SLL2@sam_data[ sams, "Community" ] <- comm
}

# Add columns for latitude and longitude
for (lat_lon in names(lat_lon.sll2)) {
  IDs <- lat_lon.sll2[[lat_lon]]
  sams <- sample_names( subset_samples(SLL2, School_ID %in% IDs) )
  
  lat <- strsplit(lat_lon, ', ')[[1]][1]
  lat <- ifelse(endsWith(lat, 'N'),
                as.numeric(strsplit(lat, ' ')[[1]][1]),
                as.numeric(strsplit(lat, ' ')[[1]][1]) * -1)
  
  lon <- strsplit(lat_lon, ', ')[[1]][2]
  lon <- ifelse(endsWith(lon, 'E'),
                as.numeric(strsplit(lon, ' ')[[1]][1]),
                as.numeric(strsplit(lon, ' ')[[1]][1]) * -1)
  
  SLL2@sam_data[ sams, "Latitude" ] <- lat
  SLL2@sam_data[ sams, "Longitude" ] <- lon
}






# ********************************************************************************************************* #
# /gpfs/projects/bsc40/current/jwillis/SLL/Part_2/annick/Rutas_v1.3.xls
schools.sll2 <- list("Escola Virolai"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.01") ],
                     "Fabrica de Sol"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.02") ],
                     "IES Miquel Marti i Pol"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.03") ],
                     "Institut Frances"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.04") ],
                     "Cafe a tres bandas"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.05") ],
                     "Collegi Sant Vicenc de Paul"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.06") ],
                     "Caixaforum"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.07") ],
                     "Collegi Sant Pere"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.08") ],
                     "IES Josep Maria Llompart"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.09") ],
                     "IES Ximen dUrrea"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.10") ],
                     "IES Misericordia"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.11") ],
                     "Octubre centro de cultura contemporanea"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.12") ],
                     "IES Gilabert de Centelles"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.13") ],
                     "Colegio La Salle"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.14") ],
                     "Padres de Alcora y Ximen"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.15") ],
                     "Nuestra Senora del Socorro"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.16") ],
                     "IES Pedanias Altas"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.17") ],
                     "Colegio La Merced Fuensanta"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.18") ],
                     "IES Florida Blanca"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.19") ],
                     "IES Los Cantos"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.20") ],
                     "Museo de la Ciencia y el Agua"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.21") ],
                     "Colegio El Pinar"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.22") ],
                     "Bar La Pinta Craft Beer"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.23") ],
                     "Down Malaga"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.24") ],
                     "Aula de Mayores de la Universidad de Malaga"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.25") ],
                     "Celiacos Malaga"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.26") ],
                     "IES Torreblanca"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.27") ],
                     "Centro Civico Las Sirenas"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.28") ],
                     "Meeting Illumina - Hotel Barcelo"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.29") ],
                     "IES El Espinillo"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.30") ],
                     "IES Rosa Chacel"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.31") ],
                     "IES Carmen Martin Gaite"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.32") ],
                     "Medialab Prado"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.33") ],
                     "IES Carpetania"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.34") ],
                     "IES David Bujan"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.35") ],
                     "Ayuntamiento de Cambre"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.36") ],
                     "IES de Ponteceso"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.37") ],
                     "CPR Plurilingue Santiago"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.38") |
                                                                      startsWith(sample_names(SLL2), "SLL.39")],
                     "Colegio Nino Jesus de Praga"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.40") ],
                     "Instalaciones de Down Vigo"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.41") ],
                     "Colegio San Jose Carmelitas"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.42") ],
                     "Down Bilbao"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.43") ],
                     "Asociacion de Celiacos de Euskadi"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.44") ],
                     "Urretxu Zumarraga Ikastola"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.45") ],
                     "COCEMFE-Cantabria"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.46") ],
                     "IES Rio Arba de Tauste"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.47") ],
                     "IES Goya"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.48") ],
                     "INS Deulofeu de Figueres"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.49") ],
                     "Wetlab Laboratorios Cesar en Etopia"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.50") ],
                     "IES Olorda"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.51") ],
                     "Institut Ramon Berenguer IV"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.52") ],
                     "NS Leonardo da Vinci Sant Cugat del Valles"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.53") ],
                     "CRG"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.54") ], 
                     "DIY bio BCN"=sample_names(SLL2)[ startsWith(sample_names(SLL2), "SLL.55") ])

# first add column of school name to metadata
for (school in names(schools.sll2)) {
  SLL2@sam_data[ schools.sll2[[ school ]], "School_name" ] <- school
}
# ********************************************************************************************************* #








# ********************************************************************************************************* #
###### Add various water values for each city ######
# ********************************************************************************************************* #

water_data <- read.delim(sprintf("%s/water_qualities_by_city.csv", water_dir), row.names = 1)
rownames(water_data) <- gsub("Port d’Alcúdia", "Port d'Alcúdia", rownames(water_data))
rownames(water_data) <- gsub("L’Alcora", "L'Alcora", rownames(water_data))
water_data <- water_data[ , 1:21] # no need to include columns Taste, Smell, Color, CO2, SH2, Temperature

cont_water_data <- c("Conductivity","water_pH","Dry.residue.at.180.C","Dry.residue.at.110.C","Cl","F","HCO3","NO3",
                     "SO4","Na","K","Li","Ca","Mg","Sr","Hardness","Alcalinity")#,"CO3"
group_water_data <- c("Mineralization", "Composition", "Hardness_category")

for (city in names(cities.sll2)) {
  print(city)
  school_IDs <- cities.sll2[[city]]
  samps <- rownames(SLL2@sam_data[SLL2@sam_data$School_ID %in% school_IDs, ])
  
  for (col in colnames(water_data)) {
    SLL2@sam_data[ samps, col ] <- water_data[ city, col ]
  }
}


# ********************************************************************************************************* #
###### Add population size for sample based on its city ######
# ********************************************************************************************************* #

city_pops <- list( "Barcelona" = 1605000, "Roda de Ter" = 6122, "Sant Feliu de Llobregat" = 44086, "San Cugat del Vallés" = 88921, "Amposta" = 20654, "Figueres" = 45726, "Lleida" = 138144,
                   "Palma de Mallorca" = 400578, 
                   "L'Alcora" = 10591, "Nules" = 13442, "Valencia" = 786189, "Benetússer" = 14505,
                   "La Paca" = 1318, "Murcia" = 439889, "Bullas" = 11753,
                   "Alhaurín de la Torre" = 38794, "Málaga" = 569130, "Sevilla" = 693878, 
                   "Madrid" = 3142000, "Colmenar Viejo" = 48020, "Moralzarzal" = 12372, "Yepes" = 5066,
                   "Os Campons" = 258, "Cambre" = 24141, "Ponteceso" = 5703, "Ponteareas" = 22963, "Vigo" = 294098,
                   "Santurtzi" = 46284, "Bilbao" = 345122, "Villarreal de Urrechua" = 6805,
                   "Santander" = 172656,
                   "Tauste" = 6941, "Zaragoza" = 664953 )

for (city in names(cities.sll2)) {
  IDs <- cities.sll2[[city]]
  samps <- rownames(SLL2@sam_data[SLL2@sam_data$School_ID %in% IDs, ])
  
  pop <- city_pops[[city]]
  
  SLL2@sam_data[samps,"Population"] <- pop
}
# ********************************************************************************************************* #













# ****************************************************************************************************************** #
# Add columns for other diseases ####
med_cols <- c("Other_disorder","Other_medical_matters", "Other_medications")

# DIABETES
samps.diabetes <- sample_names(subset_samples(SLL2, grepl('diabet',Other_disorder,ignore.case=T) | grepl('diabet',Other_medical_matters,ignore.case=T) | 
                                                grepl('insulin',Other_medications,ignore.case=T) ))
samps.diabetes <- samps.diabetes[ ! samps.diabetes %in% c("13.43") ] # "Abuela paterna diabetica"...not interesting here
SLL2@sam_data[ , "Diabetes"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.diabetes, 'Yes', 'No'))

# HYPERTENSION
samps.hiperten <- sample_names(subset_samples(SLL2, grepl('tension',Other_disorder,ignore.case=T) | grepl('hipertenso',Other_disorder,ignore.case=T) | 
                                                grepl('tens',Other_medical_matters,ignore.case=T) | grepl('pressure',Other_medical_matters,ignore.case=T) | 
                                                grepl('tensio',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Hypertension"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.hiperten, 'Yes', 'No'))

# CHOLESTEROL
samps.colest <- sample_names(subset_samples(SLL2, grepl('colester',Other_disorder,ignore.case=T) | grepl('colester',Other_medical_matters,ignore.case=T) | 
                                              grepl('colester',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Cholesterol"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.colest, 'Yes', 'No'))

# DEPRESSION
samps.depres <- sample_names(subset_samples(SLL2, grepl('depres',Other_disorder,ignore.case=T) | grepl('depres',Other_medical_matters,ignore.case=T) | 
                                              grepl('depres',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Depression"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.depres, 'Yes', 'No'))

# ANXIETY
samps.anxi <- sample_names(subset_samples(SLL2, grepl('ansiedad',Other_disorder,ignore.case=T) | grepl('ansiedad',Other_medical_matters,ignore.case=T) | 
                                            grepl('ansiol',Other_medications,ignore.case=T) | grepl('loraz',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Anxiety"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.anxi, 'Yes', 'No'))

# HEADACHES
samps.migr <- sample_names(subset_samples(SLL2, grepl('jaqueca',Other_disorder,ignore.case=T) | grepl('migran',Other_disorder,ignore.case=T) | 
                                            grepl('jaqueca',Other_medical_matters,ignore.case=T) | grepl('migran',Other_medical_matters,ignore.case=T) | 
                                            grepl('migran',Other_medications,ignore.case=T) | grepl('dolor de cabeza',Other_medications,ignore.case=T) ))#| 
# grepl('ibuprofeno',Other_medications,ignore.case=T) | grepl('iboprufeno',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Headaches"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.migr, 'Yes', 'No'))

# LACTOSE INTOLERANT
samps.lact <- sample_names(subset_samples(SLL2, grepl('lactos',Other_disorder,ignore.case=T) | 
                                            grepl('lactos',Other_medical_matters,ignore.case=T) | grepl('lacteo',Other_medical_matters,ignore.case=T) | 
                                            grepl('lactos',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Lactose_intolerant"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.lact, 'Yes', 'No'))

# GASTRITIS
samps.gast <- sample_names(subset_samples(SLL2, grepl('gastritis',Other_disorder,ignore.case=T) | grepl('gastritis',Other_medical_matters,ignore.case=T) | 
                                            grepl('abdomin',Other_medical_matters,ignore.case=T) | grepl('gastritis',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Gastritis"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.gast, 'Yes', 'No'))

# INTESTINAL ISSUES
samps.intes <- sample_names(subset_samples(SLL2, grepl('colitis',Other_disorder,ignore.case=T) | grepl('colon',Other_disorder,ignore.case=T) |
                                             grepl('colon',Other_medical_matters,ignore.case=T) | grepl('crohn',Other_disorder,ignore.case=T) |
                                             grepl('chron',Other_disorder,ignore.case=T) | grepl('crohn',Other_medical_matters,ignore.case=T) ))
SLL2@sam_data[ , "Intestinal_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.intes, 'Yes', 'No'))

# ANEMIA
samps.anem <- sample_names(subset_samples(SLL2, grepl('anemia',Other_disorder,ignore.case=T) | grepl('anemia',Other_medical_matters,ignore.case=T) | 
                                            grepl('bajo en hierro',Other_medical_matters,ignore.case=T) ))
SLL2@sam_data[ , "Anemia"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.anem, 'Yes', 'No'))

# SINUSITIS
samps.sinu <- sample_names(subset_samples(SLL2, grepl('sinusitis',Other_disorder,ignore.case=T) | grepl('sinusitis',Other_medical_matters,ignore.case=T) | 
                                            grepl('sinusitis',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Sinusitis"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.sinu, 'Yes', 'No'))

# FIBROSIS CARRIER
samps.fqcar <- sample_names(subset_samples(SLL2, grepl('portador',Other_disorder,ignore.case=T) | grepl('portador',Other_medical_matters,ignore.case=T) | 
                                            grepl('portador',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Fibrosis_carrier"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.fqcar, 'Yes', 'No'))

# THYROID ISSUES
samps.tyro <- sample_names(subset_samples(SLL2, grepl('tiroid',Other_disorder,ignore.case=T) | grepl('tiroid',Other_medical_matters,ignore.case=T) | 
                                            grepl('tiroid',Other_medications,ignore.case=T) ))
samps.hypotyro <- sample_names(subset_samples(SLL2, grepl('hipotiroid',Other_disorder,ignore.case=T) | grepl('tirox',Other_medications,ignore.case=T) ))
samps.hypertyro <- sample_names(subset_samples(SLL2, grepl('hipertiroid',Other_disorder,ignore.case=T) | grepl('tirodril',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Thyroid_issue"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.tyro, 'Yes', 'No'))
SLL2@sam_data[ , "Hypothyroidism"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.hypotyro, 'Yes', 'No'))

# CANCER
samps.canc <- sample_names(subset_samples(SLL2, grepl('cancer',Other_disorder,ignore.case=T) | grepl('cancer',Other_medical_matters,ignore.case=T) | 
                                            grepl('tumor',Other_medical_matters,ignore.case=T) | grepl('carcino',Other_disorder,ignore.case=T) ))
SLL2@sam_data[ , "Cancer"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.canc, 'Yes', 'No'))

# TRANSPLANT
samps.trans <- sample_names(subset_samples(SLL2, grepl('plant',Other_disorder,ignore.case=T) | grepl('plant',Other_medical_matters,ignore.case=T) | 
                                            grepl('plant',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Transplant"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.trans, 'Yes', 'No'))

# IMMUNE ISSUES
samps.immu <- sample_names(subset_samples(SLL2, grepl('inmun',Other_disorder,ignore.case=T) | grepl('inmun',Other_medical_matters,ignore.case=T) | 
                                             grepl('inmun',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Immune_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.immu, 'Yes', 'No'))

# SKIN ISSUES
samps.skin <- sample_names(subset_samples(SLL2, grepl('dermatitis',Other_disorder,ignore.case=T) | grepl('psorias',Other_disorder,ignore.case=T) |
                                            grepl('dermatitis',Other_medical_matters,ignore.case=T) | grepl('eccema',Other_medical_matters,ignore.case=T) | 
                                            grepl('exema',Other_medical_matters,ignore.case=T) | grepl('atopica',Other_medical_matters,ignore.case=T) | 
                                            grepl('psorias',Other_medications,ignore.case=T) | grepl('dercutane',Other_medical_matters,ignore.case=T) |
                                            grepl('cutan',Other_medications,ignore.case=T) | grepl('granos',Other_medications,ignore.case=T) |
                                            grepl('acne',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Skin_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.skin, 'Yes', 'No'))

# LUNG ISSUES
samps.lung <- sample_names(subset_samples(SLL2, grepl('pulm',Other_disorder,ignore.case=T) | grepl('pulm',Other_medical_matters,ignore.case=T) | 
                                            grepl('pulm',Other_medications,ignore.case=T) | grepl('asma',Other_disorder,ignore.case=T) |
                                            grepl('asma',Other_medications,ignore.case=T) | grepl('bronquitis',Other_disorder,ignore.case=T) |
                                            grepl('bronquitis',Other_medical_matters,ignore.case=T) | grepl('neumonia',Other_medical_matters,ignore.case=T) ))
SLL2@sam_data[ , "Lung_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.lung, 'Yes', 'No'))

# CIRCULATORY ISSUES
samps.circ <- sample_names(subset_samples(SLL2, grepl('artrosis',Other_disorder,ignore.case=T) | grepl('artrosis',Other_medical_matters,ignore.case=T) |
                                            grepl('bloqueo',Other_disorder,ignore.case=T) | grepl('bloqueo',Other_medical_matters,ignore.case=T) | 
                                            grepl('cardi',Other_disorder,ignore.case=T) | grepl('cardi',Other_medical_matters,ignore.case=T) | 
                                            grepl('coraz',Other_medical_matters,ignore.case=T) | grepl('arritmia',Other_medical_matters,ignore.case=T) |
                                            grepl('aortica',Other_medical_matters,ignore.case=T)))
SLL2@sam_data[ , "Circulatory_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.circ, 'Yes', 'No'))

# KDINEY ISSUES
samps.kidn <- sample_names(subset_samples(SLL2, grepl('renal',Other_disorder,ignore.case=T) | grepl('renal',Other_medical_matters,ignore.case=T) |
                                            grepl('rinon',Other_disorder,ignore.case=T) | grepl('rinon',Other_medical_matters,ignore.case=T) | 
                                            grepl('nefr',Other_disorder,ignore.case=T) | grepl('orina',Other_medical_matters,ignore.case=T) |
                                            grepl('renitis',Other_disorder,ignore.case=T) ))
SLL2@sam_data[ , "Kidney_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.kidn, 'Yes', 'No'))

# TONSILS
samps.tons <- sample_names(subset_samples(SLL2, grepl('amigdal',Other_medical_matters,ignore.case=T) | grepl('vegetaciones',Other_medical_matters,ignore.case=T) | 
                                            grepl('angina',Other_medical_matters,ignore.case=T) ))
SLL2@sam_data[ , "Tonsil_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.tons, 'Yes', 'No'))


# PANCREAS
samps.panc <- sample_names(subset_samples(SLL2, grepl('pancrea',Other_medications,ignore.case=T) | 
                                            grepl('kreon',Other_medications,ignore.case=T) | 
                                            grepl('pancrea',Other_disorder,ignore.case=T) ))
SLL2@sam_data[ , "Pancreas_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.panc, 'Yes', 'No'))


# SAMPLE COLLECTION ISSUES
samps.coll <- sample_names(subset_samples(SLL2, grepl('Justo Despues Del Recreo',Other_medical_matters,ignore.case=T) | 
                                            grepl('Poca Cantidad',Other_medical_matters,ignore.case=T) | 
                                            grepl('Se Trago El Pbs',Other_medical_matters,ignore.case=T) ))
SLL2@sam_data[ , "Samp_collect_issues"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.coll, 'Yes', 'No'))

# BIRTH CONTROL
samps.bc <- sample_names(subset_samples(SLL2, grepl('concep',Other_medical_matters,ignore.case=T) | grepl('concep',Other_medications,ignore.case=T) | 
                                          grepl('pildora',Other_medications,ignore.case=T) | grepl('Anticnceptvo',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Birth_control"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.bc, 'Yes', 'No'))

# ANTIHISTAMINES
samps.hist <- sample_names(subset_samples(SLL2, grepl('stamin',Other_medications,ignore.case=T) ))
SLL2@sam_data[ , "Antihistamines"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.hist, 'Yes', 'No'))


additional_diseases <- c("Diabetes","Hypertension","Cholesterol","Depression","Anxiety","Headaches",
                         "Lactose_intolerant","Gastritis","Intestinal_issues","Anemia","Sinusitis",
                         "Fibrosis_carrier","Thyroid_issue","Hypothyroidism","Cancer","Transplant",
                         "Immune_issues","Skin_issues","Lung_issues","Circulatory_issues","Kidney_issues",
                         "Tonsil_issues","Samp_collect_issues","Birth_control","Antihistamines")

# # check to see if any missing diabetes samples
# samps.to.check <- samps.diabetes
# samps.to.check <- samps.hiperten
# samps.to.check <- samps.migr
# samps.to.check <- samps.gast
# samps.to.check <- samps.canc
# samps.to.check <- samps.lung
# 
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.to.check ], "Other_disorder"])
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.to.check ], "Other_medical_matters"])
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.to.check ], "Other_medications"])
# 
# samps.not <- sort(unique(c(samps.diabetes,samps.hiperten,samps.colest,samps.depres,samps.anxi,samps.migr,samps.lact,
#                            samps.gast,samps.intes,samps.anem,samps.sinu,samps.fqcar,samps.tyro,samps.canc,samps.trans,
#                            samps.immu,samps.skin,samps.lung,samps.circ,samps.kidn,samps.tons,samps.coll,samps.bc,samps.hist)))
# 
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.not ], "Other_disorder"])
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.not ], "Other_medical_matters"])
# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.not ], "Other_medications"])
# 
# Hierro
# Epilepsia
# Artritis
# Herpes  Liquen Plano
# Sobrepeso Obesidad







# ****************************************************************************************************************** #

# Add columns for types of other small furry animals ####

# RABBITS
samps.rabb <- sample_names(subset_samples(SLL2, grepl('conejo',Pets.Type.Small_furry_animals,ignore.case=T) ))
SLL2@sam_data[ , "Pets.Rabbits"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.rabb, 'Yes', 'No'))
SLL2@sam_data[ , "Pets.Number.Rabbits"] <- sapply(sample_names(SLL2), function(x) 
  ifelse(SLL2@sam_data[x, "Pets.Rabbits"]=='Yes', 
         as.numeric(as.matrix(SLL2@sam_data[x, "Pets.Number.Small_furry_animals"])), 
         0))

# RODENTS
samps.rodent <- sample_names(subset_samples(SLL2, grepl('ardilla',Pets.Type.Small_furry_animals,ignore.case=T) | 
                                              grepl('cobayo',Pets.Type.Small_furry_animals,ignore.case=T) |
                                              grepl('chinchilla',Pets.Type.Small_furry_animals,ignore.case=T) |
                                              grepl('hamster',Pets.Type.Small_furry_animals,ignore.case=T) |
                                              grepl('jerbo',Pets.Type.Small_furry_animals,ignore.case=T) | 
                                              grepl('rata',Pets.Type.Small_furry_animals,ignore.case=T) ))
SLL2@sam_data[ , "Pets.Rodents"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.rodent, 'Yes', 'No'))
SLL2@sam_data[ , "Pets.Number.Rodents"] <- sapply(sample_names(SLL2), function(x) 
  ifelse(SLL2@sam_data[x, "Pets.Rodents"]=='Yes', 
         as.numeric(as.matrix(SLL2@sam_data[x, "Pets.Number.Small_furry_animals"])), 
         0))

# MAMMALS
samps.mamm <- sample_names(subset_samples(SLL2, grepl('erizo',Pets.Type.Small_furry_animals,ignore.case=T) | 
                                            grepl('huron',Pets.Type.Small_furry_animals,ignore.case=T) |
                                            grepl('mamifero',Pets.Type.Small_furry_animals,ignore.case=T) ))
samps.mamm <- sort(unique(c(samps.mamm, samps.rabb, samps.rodent)))
SLL2@sam_data[ , "Pets.Mammals"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.mamm, 'Yes', 'No'))
SLL2@sam_data[ , "Pets.Number.Mammals"] <- sapply(sample_names(SLL2), function(x) 
  ifelse(SLL2@sam_data[x, "Pets.Mammals"]=='Yes', 
         as.numeric(as.matrix(SLL2@sam_data[x, "Pets.Number.Small_furry_animals"])), 
         0))

# HORSES
samps.horse <- sample_names(subset_samples(SLL2, grepl('caballo',Pets.Type.Small_furry_animals,ignore.case=T) ))
SLL2@sam_data[ , "Pets.Horses"] <- sapply(sample_names(SLL2), function(x) ifelse(x %in% samps.horse, 'Yes', 'No'))
SLL2@sam_data[ , "Pets.Number.Horses"] <- sapply(sample_names(SLL2), function(x)
  ifelse(SLL2@sam_data[x, "Pets.Horses"]=='Yes',
         as.numeric(as.matrix(SLL2@sam_data[x, "Pets.Number.Small_furry_animals"])),
         0))


# table(SLL2@sam_data[sample_names(SLL2)[ ! sample_names(SLL2) %in% samps.rabb ], "Pets.Type.Small_furry_animals"])

# ****************************************************************************************************************** #









# ********************************************************************************************************* #
###### Add data for Fungal identification ######
# ********************************************************************************************************* #
fung <- read.csv(sprintf("%s/2017_SLL2_DATABASE_FINAL_SENTon20190205.csv", p2_dir), row.names = 2, stringsAsFactors = F)
rownames(fung) <- sapply(rownames(fung), function(x) paste0("SLL.", gsub("-","\\.",x)))
fung <- fung[ sample_names(SLL2)[ sample_names(SLL2) %in% rownames(fung) ], ]
fung[ is.na(fung) ] <- "None"

# fix strange numbers
fung$X..Yeast.colonies[fung$X..Yeast.colonies==">100"] <- 100
fung$X..Yeast.colonies[fung$X..Yeast.colonies==">1000"] <- 1000

# change " + " to " and " so it can work with phyloseq column names
fung$Yeast.species <- gsub(" \\+ ", " AND ", fung$Yeast.species)
fung$Bacteria.detected.by.MALDI <- gsub(" \\+ ", " AND ", fung$Bacteria.detected.by.MALDI)

# remove strange parts of some yeast names
fung$Yeast.species <- gsub("Candida_guilliermondii\\[ana\\] \\(Pichia_guilliermondii\\[teleo\\]#\\)", "Candida guillermondii", fung$Yeast.species)
fung$Yeast.species <- gsub("Candida_krusei\\[ana\\] \\(Issatchenkia_orientalis\\[teleo\\]#\\)", "Candida krusei", fung$Yeast.species)
fung$Yeast.species <- gsub("Candida_kefyr\\[ana\\] \\(Kluyveromyces_marxianus_ssp_marxianus\\[teleo\\]#\\)", "Candida kefyr", fung$Yeast.species)
fung$Yeast.species <- gsub("Candida_lusitaniae\\[ana\\] \\(Clavispora_lusitaniae\\[teleo\\]#\\)", "Candida lusitaniae", fung$Yeast.species)
fung$Yeast.species <- gsub("Kloeckera_apiculata\\[ana\\] \\(Hanseniaspora_uvarum\\[teleo\\]#\\)", "Kloeckera apiculata", fung$Yeast.species)

# fix names with typos
fung$Yeast.species[ fung$Yeast.species=="Sarocladiun strictum" ] <- "Sarocladium strictum"
fung$Yeast.species <- gsub("Candidad", "Candida", fung$Yeast.species)
fung$Yeast.species[ fung$Yeast.species=="Debaryomyces hansenii " ] <- "Debaryomyces hansenii"

fung$Bacteria.detected.by.MALDI <- gsub("fluorescens ", "fluorescens", fung$Bacteria.detected.by.MALDI)
fung$Bacteria.detected.by.MALDI[ fung$Bacteria.detected.by.MALDI=="Ralstonia insidiosa " ] <- "Ralstonia insidiosa"

# vector of all yeast species detected (and combos of yeast species)
fung.specs <- sort(unique(c(fung$Yeast.species,
                            unname(unlist(sapply(names(table(fung$Yeast.species)), function(x) strsplit(x," AND ")))))))
fung.specs <- fung.specs[ ! fung.specs %in% c("","None") ]

# vector of all bacteria species detected (and combos of bacteria species)
fung.bact <- sort(unique(c(fung$Bacteria.detected.by.MALDI,
                           unname(unlist(sapply(names(table(fung$Bacteria.detected.by.MALDI)), function(x) strsplit(x," AND ")))))))
fung.bact <- fung.bact[ ! fung.bact %in% c("","None") ]

# Add columns to SLL2@sam_data
SLL2@sam_data[ rownames(fung), "MALDI.Num_Yeast_Colonies" ] <- as.numeric(fung$X..Yeast.colonies)
SLL2@sam_data[ rownames(fung), "MALDI.Num_Mold_Colonies" ]  <- as.numeric(fung$X..Mold.colonies)
SLL2@sam_data[ rownames(fung), "MALDI.Yeast_detected" ]     <- ifelse(fung$X..Yeast.colonies %in% c("None","0"), "No","Yes")
SLL2@sam_data[ rownames(fung), "MALDI.Mold_detected" ]      <- ifelse(fung$X..Mold.colonies %in% c("None","0"), "No","Yes")
SLL2@sam_data[ rownames(fung), "MALDI.Bacteria_detected" ]  <- chartr("yn", "YN", gsub("None","no",fung$Bacteria))

for (fs in fung.specs) {
  SLL2@sam_data[ rownames(fung), sprintf("MALDI.Yeast.%s", gsub(" ","_",fs)) ] <- ifelse(fung$Yeast.species == fs, "Yes", "No")
}

for (fb in fung.bact) {
  SLL2@sam_data[ rownames(fung), sprintf("MALDI.Bacteria.%s", gsub(" ","_",fb)) ] <- ifelse(fung$Bacteria.detected.by.MALDI == fb, "Yes", "No")
}

fung.specs.Qs <- sprintf("MALDI.Yeast.%s", gsub(" ","_",fung.specs))
fung.bact.Qs <- sprintf("MALDI.Bacteria.%s", gsub(" ","_",fung.bact))


# ********************************************************************************************************* #




# ******************************************************************** #
# More appropriate versions of some MALDI variables ####
maldi_vars <- c("MALDI.Yeast.Candida_albicans","MALDI.Yeast.Candida_albicans_AND_Candida_dubliniensis",
                "MALDI.Yeast.Candida_albicans_AND_Candida_glabrata","MALDI.Yeast.Candida_albicans_AND_Candida_guillermondii",
                "MALDI.Yeast.Candida_albicans_AND_Candida_intermedia","MALDI.Yeast.Candida_albicans_AND_Candida_parapsilosis",
                "MALDI.Yeast.Candida_albicans_AND_Candida_parapsilosis_AND_Trichosporon_spp.",
                "MALDI.Yeast.Candida_albicans_AND_Candida_tropicalis","MALDI.Yeast.Candida_albicans_AND_Debaryomyces_hansenii",
                "MALDI.Yeast.Candida_boidinii","MALDI.Yeast.Candida_dubliniensis",
                "MALDI.Yeast.Candida_glabrata","MALDI.Yeast.Candida_glabrata_AND_Candida_albicans",
                "MALDI.Yeast.Candida_glabrata_AND_Candida_krusei","MALDI.Yeast.Candida_guillermondii",
                "MALDI.Yeast.Candida_guillermondii_AND_Candida_albicans","MALDI.Yeast.Candida_intermedia",
                "MALDI.Yeast.Candida_intermedia_AND_Candida_dubliniensis_AND_Candida_albicans","MALDI.Yeast.Candida_kefyr",
                "MALDI.Yeast.Candida_krusei","MALDI.Yeast.Candida_lusitaniae",
                "MALDI.Yeast.Candida_lusitaniae_AND_Candida_parapsilosis","MALDI.Yeast.Candida_parapsilosis",
                "MALDI.Yeast.Candida_parapsilosis_AND_Candida_albicans","MALDI.Yeast.Candida_parapsilosis_AND_Candida_tropicalis",
                "MALDI.Yeast.Candida_spp.","MALDI.Yeast.Candida_tropicalis","MALDI.Yeast.Candida_zeylanoides",
                "MALDI.Yeast.Cryptococcus_spp.","MALDI.Yeast.Debaryomyces_hansenii",
                "MALDI.Yeast.Debaryomyces_hansenii_AND_Candida_zeylanoides","MALDI.Yeast.Kloeckera_apiculata",
                "MALDI.Yeast.Kloeckera_apiculata_AND_Candida_intermedia_AND_Candida_boidinii",
                "MALDI.Yeast.Rhodotorula_mucilaginosa","MALDI.Yeast.Rhodotorula_mucilaginosa_AND_Candida_parapsilosis",
                "MALDI.Yeast.Rhodotorula_mucilaginosa_AND_Cryptococcus_spp.","MALDI.Yeast.Sarocladium_strictum",
                "MALDI.Yeast.Trichosporon_spp.")

# first must update the MALDI values to reflect presence of particular yeast properly:
yeast_of_interest <- c("Candida","Candida_albicans","Candida_dubliniensis","Candida_glabrata","Candida_parapsilosis",
                       #"Candida_boidinii", # too few - no need to run any tests
                       "Candida_guillermondii",
                       #"Candida_tropicalis",
                       "Candida_intermedia","Candida_krusei","Candida_lusitaniae",
                       #"Candida_kefyr","Candida_zeylanoides","Trichosporon_spp",
                       "Cryptococcus_spp","Debaryomyces_hansenii",
                       #"Kloeckera_apiculata",,"Sarocladium_strictum"
                       "Rhodotorula_mucilaginosa")

for (yeast in yeast_of_interest) {
  
  yVars <- maldi_vars[ grepl(yeast, maldi_vars) ]
  
  SLL2@sam_data[ , sprintf("Full_MALDI.%s", yeast) ] <- sapply(sample_names(SLL2), function(x) 
    ifelse( sum(SLL2@sam_data[x, yVars]=="Yes") > 0, "Yes", "No" ))
  
  # print(yeast)
  # print(table(SLL2.meta[ , sprintf("Full_MALDI.%s", yeast) ]))
  # print(table(mTab.ds.h[ , sprintf("Full_MALDI.%s", yeast) ]))
  # print("")
  
}
# ********************************************************************************************************* #









# ********************************************************************************************************* #

# add column for which sequencing batch samples are from
sample_names(SLL2.added) <- gsub('-', '\\.', sample_names(SLL2.added))
sample_names(SLL2.added) <- paste("SLL", sample_names(SLL2.added), sep = '.')

SLL2@sam_data[ , "seqGroup"] <- sapply(sample_names(SLL2), function(x) 
  ifelse(x %in% sample_names(SLL2.added), "Two", "One"))


# ********************************************************************************************************* #

# Add a column that has all 1 value to force all to be same shape by default
SLL2@sam_data[ , "Project"] <- rep("SLL2", nrow(SLL2@sam_data))


# get rid of "No Sabe/No Contesta"
SLL2@sam_data[ SLL2@sam_data == "No Sabe/No Contesta"] <- NA

# ********************************************************************************************************* #


# save phyloseq object
saveRDS(SLL2, file = sprintf("%s/R_objects/SLL2.phyloseq.rds", p2_dir))




# Normalize counts within each sample ####
SLL2_rel = transform_sample_counts(SLL2, function(x) 100 * x/sum(x))













# ****************************************************************************************************************** #
# **************************************************************** #
# ||| Get values of taxa at each level after filtering as above ####
# **************************************************************** #

# ********************** #
get_gloms <- function(phy, tl, tTabs, clr=F) {
  
  # agglomerate taxa at given level
  tg <- tax_glom(phy, taxrank = tl)@otu_table
  # fix strange name if necessary
  if (sum(c("Family_XI","Family_XII") %in% phy@tax_table[ , tl ]) > 0) {
    # ***************************** #
    rownames(tg) <- unname(sapply( rownames(tg), function(x) {
      if (phy@tax_table[x, tl] %in% c("Family_XI","Family_XI")) {
        # since this taxon at the family level appears for both orders Clostridiales and Bacillales
        sprintf("%s.%s", phy@tax_table[x, tl], phy@tax_table[x, "Order"])
      } else {
        phy@tax_table[x, tl]
      }
    } ))
    # get appropriate name for row of agglomerated table
    rownames(tg) <- unname(sapply( rownames(tg), function(x) {
      if (x %in% c("Family_XI.Clostridiales", "Family_XI.Bacillales")) {
        names(tTabs[[ tl ]][ , tl][tTabs[[ tl ]][ , tl] == "Family_XI" & 
                                     tTabs[[ tl ]][ , "Order"] == strsplit(x,'\\.')[[1]][2]])
      } else if (x %in% c("Family_XII.Clostridiales", "Family_XII.Bacillales")) {
        names(tTabs[[ tl ]][ , tl][tTabs[[ tl ]][ , tl] == "Family_XII" & 
                                     tTabs[[ tl ]][ , "Order"] == strsplit(x,'\\.')[[1]][2]])
      } else {
        names(tTabs[[ tl ]][ , tl][tTabs[[ tl ]][ , tl] == x])
      }
    } ))
    # ***************************** #
  } else {
    rownames(tg) <- unname(sapply( rownames(tg), function(x) phy@tax_table[x , tl] ))
  }
  
  # ***************************** #
  if (clr == T) {
    tg.f <- codaSeq.filter(tg, 
                           min.reads=1000, # filter out samples with fewer reads
                           min.prop=0.001, # filter out taxa with lower minimum abundance
                           min.occurrence=0.05, # filter out taxa not appearing in at least this prop of samps
                           samples.by.row=FALSE)
    # replace 0 values with an estimate
    tg.f.n0 <- cmultRepl(t(tg.f), method="CZM", label=0)
    # generate the CLR values for plotting later
    tg.f.clr <- t(codaSeq.clr(tg.f.n0))
    
  } else {
    return(tg)
  }
  
}
# ********************** #


# first make a combined version of the taxTables tables
taxTables.both <- list()
for (tl in rev(c("Phylum", "Class", "Order", "Family", "Genus", "Species"))) {
  taxTables.both[[ tl ]] <- unique(rbind(taxTables.orig[[ tl ]], taxTables.added[[ tl ]], stringsAsFactors = F))
  taxTables.both[[ tl ]] <- taxTables.both[[ tl ]][ sort(rownames(taxTables.both[[ tl ]])), ]
}


gloms <- list()
gloms_rel <- list()
gloms_clr <- list()

for (tl in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  print(tl)
  
  
  
  gloms[[ tl ]] <- get_gloms(SLL2, tl, taxTables.both)
  gloms_rel[[ tl ]] <- get_gloms(SLL2_rel, tl, taxTables.both)
  gloms_clr[[ tl ]] <- get_gloms(SLL2, tl, taxTables.both, clr = T)
  
}

saveRDS(gloms,     sprintf("%s/R_objects/gloms.rds", p2_dir))
saveRDS(gloms_rel, sprintf("%s/R_objects/gloms_rel.rds", p2_dir))
saveRDS(gloms_clr, sprintf("%s/R_objects/gloms_clr.rds", p2_dir))
saveRDS(taxTables.both, sprintf("%s/R_objects/taxTables.both.rds", p2_dir))


# ************************************************************ #
# get matrix version of phyloseq objects (since those can be tricky to work with sometimes)
SLL2.counts <- abundances(SLL2)
SLL2.abunds <- abundances(SLL2, transform = "compositional")*100
SLL2.abun_clr <- abundances(SLL2, transform = "clr")
SLL2.meta   <- meta(SLL2)


# ****************************************************************************************************************** #








# ******************************************************************** #
# alternative method of Aitchison ####
#   from here: https://www.ucl.ac.uk/~ucfbpve/papers/VermeeschSedGeol2016/
library(provenance)

# package requires this object to be of class "compositional", cant use matrix or dataframe
f.n0.prov <- as.compositional(f.n0)
# this produces exact same distance matrix as aDist above, but is now of classes "diss" and "dist",
#    allows use of other functions in provenance package
aitch.prov <- diss(f.n0.prov, method="aitchison")


# ********** #
MDS.aitch.prov <- MDS(aitch.prov, classical=TRUE)
plot(MDS.aitch.prov, xaxt='s', yaxt='s')

PCA.aitch.prov <- PCA(f.n0.prov)
plot(PCA.aitch.prov)


# ********** #
# comp10 <- names(sort(colSums(f.n0), decreasing = T))[1:10]
# fit.aitch.prov <- envfit(MDS.aitch.prov, f.n0[ , comp10])
comp10 <- names(sort(rowSums(gloms_clr$Genus), decreasing = T))[1:10]
fit.aitch.prov <- envfit(MDS.aitch.prov, t(gloms_clr$Genus)[ , comp10])

plot(MDS.aitch.prov$points, xaxt='s', yaxt='s', type="n")
points(MDS.aitch.prov$points, pch = 16, col = c("cyan2", "purple")[as.numeric(as.factor(SLL2.meta$seqGroup))])
# text(MDS.aitch.prov$points, col = as.numeric(as.factor(SLL2.meta$Chronic_disorder))+2)
plot(fit.aitch.prov, col="black", cex=1.75)


# ********** #
fitArrows <- as.data.frame(fit.aitch.prov$vectors$arrows)

# arrowLengthMult <- sqrt(fit.aitch.prov$vectors$r) *0.25
arrowLengthMult <- fit.aitch.prov$vectors$r * 50

plot_ordination(SLL2, MDS.aitch.prov, type="sites", color="seqGroup") + 
  #scale_colour_manual(values=c("2w"="green", "8w"="red", "1yr"="blue")) + 
  theme_bw() + 
  stat_ellipse() +
  ggtitle("Aitchison") +
  geom_segment(data = fitArrows, arrow=arrow(length=unit(1/2, "picas")), color = "black",
               # x*(-1) so that goes in same direction as this plot
               aes(x=0, xend=(fitArrows$Dim1*arrowLengthMult*(-1)),
                   y=0, yend=(fitArrows$Dim2*arrowLengthMult)) ) +
  annotate("text", x = (fitArrows$Dim1*arrowLengthMult*(-1))-0.02, y = (fitArrows$Dim2*arrowLengthMult),
           label = rownames(fitArrows), fontface=2, color="black", size=6)



# ********** #
fit.aitch.prov.all <- envfit(MDS.aitch.prov, t(gloms_clr$Genus))

# ******************************************************************** #













# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####










# ****************************************************************************************************************** #

# ******************************************************* #
# Pie charts / Donut charts of abundances ####
# ******************************************************* #

# ******************************************************** #

min_abund <- 0.8
i.tax <- "Phylum"
o.tax <- "Genus"
# o.tax <- "Species"

samps <- sample_names(SLL2)
# first open the png device to prepare it to save an image at the given location
png(filename = sprintf("%s/figures/Abundances/donut_charts/All.donut_chart.png", p2_dir), 
    width = 1505, height = 951, pointsize = 16)
# then plot the image, which will get saved at the above location
plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel)
# then turn off png device so it can be repeated for next image
dev.off()


# INDIVIDUALS
for (s in sample_names(SLL2)) {
  plot_donut(i.tax, o.tax, s, min_abund, gloms_rel, to.save = T, group = "individuals", fname = s)
}

# DISORDERS
# for (disorder in c("Celiaca","Fibrosis quistica","Gingivitis-Periodontitis","Sindrome de Down","Trastorno alimenticio"))
samps <- sample_names(subset_samples(SLL2, Q39.1==1)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Celiaca")
samps <- sample_names(subset_samples(SLL2, Q39.2==1)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Fibrosis_quistica")
samps <- sample_names(subset_samples(SLL2, Q39.3==1)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Gingivitis-Periodontitis")
samps <- sample_names(subset_samples(SLL2, Q39.4==1)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Sindrome_de_Down")
samps <- sample_names(subset_samples(SLL2, Q39.5==1)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Trastorno_alimenticio")
samps <- sample_names(subset_samples(SLL2, Q39==0)); plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, 
                                                                  group = "disorders", fname = "Sano")

# REGIONS
for (comm in names(comunidades.sll2)) {
  samps <- sample_names( subset_samples(SLL2, Community == comm) )
  plot_donut(i.tax, o.tax, samps, min_abund, gloms_rel, to.save = T, group = "regions", fname = comm)
}

# ******************************************************************** #


# ****************************************************************************************************************** #












# ****************************************************************************************************************** #
# STACKED BAR charts for comparison between groups ####
# tlev <- "Species"; otutab_rel <- SLL2_rel@otu_table
tlev <- "Genus"; otutab_rel <- gloms_rel[[ tlev ]]
# tlev <- "Phylum"; otutab_rel <- gloms_rel[[ tlev ]]


# top 10 most common OTUs
top15 <- rev(names(sort(rowSums(otutab_rel), decreasing = T))[1:15])
SLL2.top15 <- as.data.frame(otutab_rel)[top15,]


teens  <- sample_names(subset_samples(SLL2, Age<19))
adult  <- sample_names(subset_samples(SLL2, Age>=19 & Age<56))
senior <- sample_names(subset_samples(SLL2, Age>=56))


# DISORDERS
# for (disorder in c("Celiaca","Fibrosis quistica","Gingivitis-Periodontitis","Sindrome de Down","Trastorno alimenticio"))
samps.disorders <- list()
samps.disorders[[ "Celíaco" ]] <- sample_names(subset_samples(SLL2, Celiac==1))
samps.disorders[[ "Fibrosis Quística" ]] <- sample_names(subset_samples(SLL2, Cystic_fibrosis==1))
samps.disorders[[ "Gingivitis-Periodontitis" ]] <- sample_names(subset_samples(SLL2, Gingivitis_periodontitis==1))
samps.disorders[[ "Síndrome de Down" ]] <- sample_names(subset_samples(SLL2, Downs_Syndrome==1))
samps.disorders[[ "Celíaco" ]] <- samps.diabetes
samps.disorders[[ "Hipertensión" ]] <- samps.hiperten
samps.disorders[[ "Migrañas" ]] <- samps.migr
samps.disorders[[ "Trastorno alimenticio" ]] <- sample_names(subset_samples(SLL2, Eating_disorder==1))
samps.disorders[[ "Sin trastorno crónico" ]] <- sample_names(subset_samples(SLL2, Chronic_disorder==0))

samps.celiac <- sample_names(subset_samples(SLL2, Celiac=="Yes"))
samps.cf <- sample_names(subset_samples(SLL2, Cystic_fibrosis==1))
samps.gp <- sample_names(subset_samples(SLL2, Gingivitis_periodontitis=="Yes"))
samps.down <- sample_names(subset_samples(SLL2, Downs_Syndrome==1))
samps.eating <- sample_names(subset_samples(SLL2, Eating_disorder==1))
samps.noDisorder <- sample_names(subset_samples(SLL2, Chronic_disorder==0))


# REGION
samps.regions <- sapply(unique(SLL2.meta$Community), function(com) rownames(SLL2.meta[SLL2.meta$Community == com,] ))
samps.regions.healthy <- lapply(samps.regions, function(x) x[x %in% samps.noDisorder])



# ******************************************************************************************************************* #




orient <- "vertical"
# orient <- "horizontal"

plot_stacked_bars(samps.to.plot = teens, lab = "Adolescentes", orientation = orient)
plot_stacked_bars(samps.to.plot = adult, lab = "Adultos", orientation = orient)
plot_stacked_bars(samps.to.plot = senior, lab = "Mayores", orientation = orient)

plot_stacked_bars(samps.to.plot = sample_names(SLL2), lab = "Todas las Muestras", orientation = orient)

plot_stacked_bars(plotType = "Ages", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.celiac, lab = "Celíaco", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.cf, lab = "Fibrosis Quística", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.down, lab = "Síndrome de Down", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.gp, lab = "Gingivitis-Periodontitis", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.diabetes, lab = "Diabetes", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.hiperten, lab = "Hipertensión", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.migr, lab = "Migrañas", orientation = orient)
plot_stacked_bars(plotType = "Ages", samps.to.plot = samps.noDisorder, lab = "El resto de población analizada", orientation = orient)

plot_stacked_bars(samps.to.plot = samps.celiac, lab = "Celíaco", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.cf, lab = "Fibrosis Quística", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.down, lab = "Síndrome de Down", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.gp, lab = "Gingivitis-Periodontitis", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.diabetes, lab = "Diabetes", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.hiperten, lab = "Hipertensión", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.migr, lab = "Migrañas", orientation = orient)
plot_stacked_bars(samps.to.plot = samps.noDisorder, lab = "El resto de población analizada", orientation = orient)

plot_stacked_bars(plotType = "All_disorders", orientation = orient)
plot_stacked_bars(plotType = "All_regions", orientation = orient)
plot_stacked_bars(plotType = "All_regions_healthy", orientation = orient)

# for all individuals
for (s in sample_names(SLL2)) {
  plot_stacked_bars(samps.to.plot = s, lab = s, orientation = orient)
}

# for all regions
for (reg in names(samps.regions)) {
  plot_stacked_bars(samps.to.plot = samps.regions[[ reg ]], lab = reg, orientation = orient)
}

# for healthy samples in all regions
for (reg in names(samps.regions.healthy)) {
  plot_stacked_bars(samps.to.plot = samps.regions.healthy[[ reg ]], lab = sprintf("%s-healthy",reg), orientation = orient)
}



# ********** #
SLL2.top15.indiv <- SLL2.top15[ , s]
SLL2.top15.disorder <- SLL2.top15[ , samps.disorder]
SLL2.top15.region <- SLL2.top15[ , samps.region]

tax.mps <- cbind(rowMeans(SLL2.top15.region), rowMeans(SLL2.top15.disorder), SLL2.top15.indiv)
colnames(tax.mps) <- c(comm,"Celiaco","Tu microbioma oral")
tax.sd <- cbind(rowSds(as.matrix(SLL2.top15.region)), rowSds(as.matrix(SLL2.top15.disorder)), rep(0,15))

tax.mps.m <- reshape2::melt(tax.mps)
colnames(tax.mps.m) <- c("Genus","Samples","value")

# # ********** #
# samps.to.plot <- teens; lab <- "Adolescentes"
# samps.to.plot <- adult; lab <- "Adultos"
# samps.to.plot <- senior; lab <- "Mayores"
# 
# SLL2.top15.toplot <- rowMeans(SLL2.top15[ , samps.to.plot])
# SLL2.top15.toplot <- c(100-sum(SLL2.top15.toplot), SLL2.top15.toplot)
# names(SLL2.top15.toplot)[1] <- "Otros"
# tax.mps.m <- melt(as.matrix(SLL2.top15.toplot))
# colnames(tax.mps.m) <- c("Genus","Samples","value")
# tax.mps.m$Samples <- lab
# 
# 
# 
# 
# ggplot(tax.mps.m, aes(x=Samples, y=value, fill=Genus)) +
#   geom_bar(aes(x=Samples, y=value, fill=Genus), stat="identity", color="black") +
#   coord_flip() + theme_minimal() + #ylim(0,100) + #scale_y_log10() +
#   scale_fill_manual(values=colPal, name="Género", guide=F) +
#   # guides(fill = guide_legend(reverse=T)) +
#   theme(plot.title = element_text(hjust=0.5, size=20), axis.title.x = element_text(size=17),
#         axis.text.y = element_text(size=17), axis.text.x = element_text(size=15),
#         legend.title = element_text(size=17), legend.text = element_text(size=15)) +
#   ggtitle(sprintf('Abundancia relativa de los géneros más comunes')) +
#   xlab('') + ylab('Abundancia media normalizada por muestra (%)') + #scale_fill_hue(name=tlev)
#   ggsave(sprintf("%s/figures/Abundances/stacked_bars/%s.png", p2_dir, gsub(' ','_',iconv(lab,to='ASCII//TRANSLIT'))),
#          width = 16.82, height = 8.41, device = "png")


# ******************************************************************************************* #

DS.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Downs_Syndrome=="Yes"])
mC.DS <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
mC.DS <- mC.DS[ SLL2.meta[mC.DS, "Downs_Syndrome"]=="No" ]

plot_stacked_bars.disorders("Downs_Syndrome", DS.samps, mC.DS)



CF.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Cystic_fibrosis=="Yes"])
mC.CF <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
mC.CF <- mC.CF[ SLL2.meta[mC.CF, "Cystic_fibrosis"]=="No" ]

plot_stacked_bars.disorders("Cystic_fibrosis", CF.samps, mC.CF)
# ******************************************************************************************* #










# ******************************************************************************************* #
# get tables of significant differences of particular taxa between particular groups
library(gtools)
print_diffs <- function(int, other, tl="Genus", reg=NA) {
  
  top15 <- rev(names(sort(rowSums(gloms_rel[[ tl ]]), decreasing = T))[1:15])
  
  prints <- list()
  for (gen in top15) {
    vals <- as.numeric(as.matrix(gloms_rel[[ tl ]])[gen, int ])
    other_vals <- as.numeric(as.matrix(gloms_rel[[ tl ]])[gen, other ])
    
    # tt <- t.test(vals, other_vals)
    tt <- wilcox.test(vals, other_vals)
    fcl <- foldchange2logratio( foldchange = foldchange(mean(vals), mean(other_vals)) )
    fc <- logratio2foldchange(fcl)
    
    if (tt$p.value*15 < 0.05 & abs(fc) > 1.5) {
      if ( is.na(reg)) {
        # prints[[gen]] <- c(p_val=tt$p.value*15, interest=unname(tt$estimate[1]), other=unname(tt$estimate[2]), fc=fc, fcLog=fcl)
        prints[[gen]] <- c(p_val=tt$p.value*15, interest=mean(vals), other=mean(other_vals), fc=fc, fcLog=fcl)
      } else {
        # prints[[gen]] <- c(p_val=tt$p.value*15, interest=unname(tt$estimate[1]), other=unname(tt$estimate[2]), fc=fc, fcLog=fcl, Region=reg)
        prints[[gen]] <- c(p_val=tt$p.value*15, interest=mean(vals), other=mean(other_vals), fc=fc, fcLog=fcl, Region=reg)
      }
    }
  }
  return(t(as.data.frame(prints)))
}
# ******************************************************************************************* #



# Celiac
print_diffs(samps.celiac, sample_names(SLL2)[!sample_names(SLL2) %in% samps.celiac])

# Celiac - ages
dt <- samps.celiac[samps.celiac %in% teens]; dnt <- samps.celiac[!samps.celiac %in% teens]
da <- samps.celiac[samps.celiac %in% adult]; dna <- samps.celiac[!samps.celiac %in% adult]
ds <- samps.celiac[samps.celiac %in% senior]; dns <- samps.celiac[!samps.celiac %in% senior]

print_diffs(dt, dnt)
print_diffs(da, dna)
print_diffs(ds, dns)



# FQ
print_diffs(samps.cf, sample_names(SLL2)[!sample_names(SLL2) %in% samps.cf])

# FQ - ages
dt <- samps.cf[samps.cf %in% teens]; dnt <- samps.cf[!samps.cf %in% teens]
da <- samps.cf[samps.cf %in% adult]; dna <- samps.cf[!samps.cf %in% adult]
ds <- samps.cf[samps.cf %in% senior]; dns <- samps.cf[!samps.cf %in% senior]

print_diffs(dt, dnt)
print_diffs(da, dna)
# print_diffs(ds, dns)



# Down
print_diffs(samps.down, sample_names(SLL2)[!sample_names(SLL2) %in% samps.down])

# Down - ages
dt <- samps.down[samps.down %in% teens]; dnt <- samps.down[!samps.down %in% teens]
da <- samps.down[samps.down %in% adult]; dna <- samps.down[!samps.down %in% adult]

print_diffs(dt, dnt)
print_diffs(da, dna)




# healthy vs disorders
print_diffs(samps.noDisorder, sample_names(SLL2)[!sample_names(SLL2) %in% samps.noDisorder])


# healthy - ages
dt <- samps.noDisorder[samps.noDisorder %in% teens]; dnt <- samps.noDisorder[!samps.noDisorder %in% teens]
da <- samps.noDisorder[samps.noDisorder %in% adult]; dna <- samps.noDisorder[!samps.noDisorder %in% adult]
ds <- samps.noDisorder[samps.noDisorder %in% senior]; dns <- samps.noDisorder[!samps.noDisorder %in% senior]

print_diffs(dt, dnt)
print_diffs(da, dna)
print_diffs(ds, dns)

print_diffs(dt, da)
print_diffs(dt, ds)
print_diffs(da, ds)



# comparing regions
reg_diffs <- list()
for (reg in sort(names(samps.regions.healthy))) {
  mat <- print_diffs(samps.regions.healthy[[ reg ]], unlist(samps.regions.healthy[names(samps.regions.healthy) != reg]), reg = reg)
  if (nrow(mat)>0) reg_diffs[[reg]] <- mat
}
reg_diffs <- as.data.frame(do.call(rbind, reg_diffs))
reg_diffs[sort(rownames(reg_diffs)),]


# plot_bar(phyloseq(otu_table(tax.mps, taxa_are_rows = T),tax_table(taxTables[["Genus"]])), fill = "Genus") +
#   coord_flip() +
#   xlab('') + ylab('Mean normalized abundance per sample (as %)') 


# ****************************************************************************************************************** #



























tlev <- "Species"; otutab_rel <- SLL2_rel@otu_table
# tlev <- "Genus"; otutab_rel <- gloms_rel[[ tlev ]]
# tlev <- "Phylum"; otutab_rel <- gloms_rel[[ tlev ]]


# top 10 most common OTUs
top10 <- names(sort(rowSums(otutab_rel), decreasing = T))[1:10]
SLL2.top10 <- as.data.frame(otutab_rel[top10,])


# ****************************************************************************************************************** #
# Bar charts of abundances ####
# ****************************************************************************************************************** #

#mean and sd of given taxa per sample
tax.mps <- rowMeans(SLL2.top10)
tax.sd <- rowSds(as.matrix(SLL2.top10))

#percentage of samples in which given taxa appears
tax.freqs <- round(rowSums(SLL2.top10 != 0) / length(colnames(SLL2.top10)) * 100, digits=2)
tax.freqs <- paste(as.character(tax.freqs),'%',sep='')


#create frame that can be used for bar plots
tax.frame <- as.data.frame(tax.mps)
tax.frame[,2] <- tax.sd
tax.frame[,3] <- tax.freqs
tax.frame[,4] <- rownames(SLL2.top10)
colnames(tax.frame) <- c("mps","sd","freqs",tlev)

#plot
ggplot(tax.frame, aes(x=reorder(tax.frame[,tlev], -mps), y=mps, 
                      fill=reorder(tax.frame[,tlev], -mps))) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mps-sd, ymax=mps+sd), width=0.2) +
  theme(axis.text.x=element_blank(), plot.title = element_text(hjust=0.5)) +
  ggtitle(sprintf('Relative abundance for most common %s', tlev)) +
  xlab(tlev) + ylab('Mean normalized abundance per sample (as %)') + scale_fill_hue(name=tlev) +
  geom_text(aes(x=reorder(tax.frame[,tlev], -mps), label=freqs), vjust=-1)

### mean abundance for these top 10, and sd
mean(colSums(SLL2.top10)) # 58.78991 (76.70266 for top 10 at Genus level)
sd(colSums(SLL2.top10)) # 6.918236 (10.1354 for top 10 at Genus level)




# ****************************************************************************************************************** #
# Box plots of abundances within samples ####
# ****************************************************************************************************************** #

ggplot(reshape2::melt(t(SLL2.top10)), aes(x=factor(Var2, levels=rownames(SLL2.top10)), y=value, 
                                fill=factor(Var2, levels=rownames(SLL2.top10)))) +
  geom_boxplot(notch = T) + theme_minimal() +
  theme(axis.text.x=element_text(angle = 45, vjust = 0.75, size=14), axis.text.y=element_text(size=14), 
        axis.title.x=element_text(size=14), plot.title = element_text(hjust=0.5)) + 
  ggtitle(sprintf('Percent of given %s per sample', tlev)) +
  xlab(tlev) + ylab("") + scale_fill_hue(guide=F) # for the paper

# ****************************************************************************************************************** #




## Average number of OTUs per found per sample
mean(colSums(SLL2_rel@otu_table != 0)) # 87.08981
mean(colSums(tax_glom(SLL2_rel, "Genus")@otu_table != 0)) # (54.76396 at genus level)




















# ****************************************************************************************************************** #
# Diversity within (alpha) and between (beta) samples ####
# ****************************************************************************************************************** #

tlev <- "Species"; otutab <- SLL2@otu_table; otutab_rel <- SLL2_rel@otu_table
# tlev <- "Genus"; otutab <- gloms[[ tlev ]]; otutab_rel <- gloms_rel[[ tlev ]]

# ********************************************* #
# Will divide samples into quartiles based on their diversity values
# We have a number of different diversity measures, so will take a look at each here:

# first create columns in sample_data for each diversity value
for (div.estimate in c("Div.Shannon","Div.Simpson","Weighted_Unifrac","Unweighted_Unifrac","Faiths.PD",
                       "Species_Richness","Bray.Curtis","Canberra")) {
  divs <- SLL2.meta[ , div.estimate]
  quarts <- quantile(divs, probs=seq(0,1,0.25))
  
  # get samples by diversity quartile
  q1.samples <- rownames(SLL2.meta[ SLL2.meta[,div.estimate] < quarts['25%'], ] )
  q2.samples <- rownames(SLL2.meta[ (quarts['25%'] <= SLL2.meta[,div.estimate]) & (SLL2.meta[,div.estimate] < quarts['50%']), ] )
  q3.samples <- rownames(SLL2.meta[ (quarts['50%'] <= SLL2.meta[,div.estimate]) & (SLL2.meta[,div.estimate] < quarts['75%']), ] )
  q4.samples <- rownames(SLL2.meta[ SLL2.meta[,div.estimate] >= quarts['75%'], ] )
  
  # add diversity group to sample_data 
  SLL2@sam_data[q1.samples, sprintf("Diversity_group_%s", div.estimate)] <- "Low"
  SLL2_rel@sam_data[q1.samples, sprintf("Diversity_group_%s", div.estimate)] <- "Low"
  SLL2@sam_data[c(q2.samples,q3.samples), sprintf("Diversity_group_%s", div.estimate)] <- "Average"
  SLL2_rel@sam_data[c(q2.samples,q3.samples), sprintf("Diversity_group_%s", div.estimate)] <- "Average"
  SLL2@sam_data[q4.samples, sprintf("Diversity_group_%s", div.estimate)] <- "High"
  SLL2_rel@sam_data[q4.samples, sprintf("Diversity_group_%s", div.estimate)] <- "High"
}

SLL2.meta <- meta(SLL2)













# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####












# ****************************************************************************************************************** #
# ******************************************************* #
# Clustering and network analyses ####
# ******************************************************* #

# Based on the steps in the Bork Enterotype paper

# # Create distance matrix using the Jensen-Shannon distance
# dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
#   KLD <- function(x,y) sum(x *log(x/y))
#   JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
#   matrixColSize <- length(colnames(inMatrix))
#   matrixRowSize <- length(rownames(inMatrix))
#   colnames <- colnames(inMatrix)
#   resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
#   
#   inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
#   
#   for(i in 1:matrixColSize) {
#     for(j in 1:matrixColSize) { 
#       resultsMatrix[i,j] = JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))
#     }
#   }
#   colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
#   as.dist(resultsMatrix) -> resultsMatrix
#   attr(resultsMatrix, "method") <- "dist"
#   return(resultsMatrix) 
# }






# # get prediction strength objects for each distance measure ####
# library(fpc)
# ps.dists.temp <- list()
# ps.dists.temp[[ "JSD" ]] <- prediction.strength(as.dist(jsd), cutoff = 0.75)
# ps.dists.temp[[ "Weighted_Unifrac" ]] <- prediction.strength(as.dist(weighted_Unifrac), cutoff = 0.75)
# ps.dists.temp[[ "Unweighted_Unifrac" ]] <- prediction.strength(as.dist(unweighted_Unifrac), cutoff = 0.75)
# ps.dists.temp[[ "VAW_GUnifrac" ]] <- prediction.strength(as.dist(guni.VAW), cutoff = 0.75)
# ps.dists.temp[[ "a0_GUnifrac" ]] <- prediction.strength(as.dist(guni.a0), cutoff = 0.75)
# ps.dists.temp[[ "a05_GUnifrac" ]] <- prediction.strength(as.dist(guni.a05), cutoff = 0.75)
# ps.dists.temp[[ "Bray" ]] <- prediction.strength(as.dist(bray), cutoff = 0.75)
# ps.dists.temp[[ "Jaccard" ]] <- prediction.strength(as.dist(jaccard), cutoff = 0.75)
# ps.dists.temp[[ "Canberra" ]] <- prediction.strength(as.dist(canberra), cutoff = 0.75)
# ps.dists.temp[[ "Aitchison" ]] <- prediction.strength(as.dist(aitch), cutoff = 0.75)
# 
# saveRDS(ps.dists.temp, sprintf("%s/R_objects/ps.SLL2.dists.temp.rds", p2_dir))


ps.dists <- readRDS(sprintf("%s/R_objects/ps.SLL2.dists.temp.rds", p2_dir))


# ************************************************************************************************ #




tl <- "Species"; tagl <- SLL2_rel
# tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2_rel))

# dist_meas <- "JSD"
# dist_meas <- "JSD_alt"
# dist_meas <- "Weighted_Unifrac"
# dist_meas <- "Unweighted_Unifrac"
# dist_meas <- "Bray"

dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
           "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
           "Bray","Jaccard","Canberra","Aitchison")

for (dist_meas in dists) {
  
  glomTab <- gloms_rel
  
  print(dist_meas)
  clu <- get_clusters( tagl, tl, dist_meas, ps.dists, glomTab, saveClusPlots = T,
                       clusPlotsDir = sprintf("%s/figures/Stomatotypes-cluster_measures", p2_dir, dist_meas))
  clus_var_name <- clu[[ "clus_var_name" ]]
  cluster_full <- clu[[ "cluster_full" ]]
  nclusters <- clu[[ "nclusters" ]]
  tax.to.clust <- clu[[ "tax.to.clust" ]]
  dist.matrix <- clu[[ "dist.matrix" ]]
  clus.meas.norm <- clu[[ "clus.meas.norm" ]]
  
  #add stomatotype to sample_data
  #first must remove the current values of the column if they are there and being rewritten, 
  # because can cause problems when adding factors with different NA values
  SLL2@sam_data <- SLL2@sam_data[ , colnames(SLL2@sam_data) != sprintf("%s_%s",clus_var_name, dist_meas)]
  SLL2_rel@sam_data <- SLL2_rel@sam_data[ , colnames(SLL2_rel@sam_data) != sprintf("%s_%s",clus_var_name, dist_meas)]
  
  # must specify row order in the case that not all samples were included for the clustering, 
  # since the order in cluster_full.stool will be different
  SLL2@sam_data[ names(cluster_full), sprintf("%s_%s",clus_var_name, dist_meas) ] <- as.factor(cluster_full)
  SLL2_rel@sam_data[ names(cluster_full), sprintf("%s_%s",clus_var_name, dist_meas) ] <- as.factor(cluster_full)
}


# a few extras just to check them out (since the PS values suggest maybe other nclusters are better)
for (dis in c("JSD.3", "Weighted_Unifrac.3", "Weighted_Unifrac.4", "VAW_GUnifrac.3")) {
# for (dis in c("JSD.3", "Weighted_Unifrac.2", "Weighted_Unifrac.4")) {
  
  glomTab <- gloms_rel
  
  print(dis)
  dist_meas <- strsplit(dis, '\\.')[[1]][1]
  nc <- as.integer(strsplit(dis, '\\.')[[1]][2])
  
  clu <- get_clusters( tagl, tl, dist_meas, ps.dists, glomTab, nClustForce = nc )
  
  # remove current values if necessary
  SLL2@sam_data <- SLL2@sam_data[ , colnames(SLL2@sam_data) != sprintf("%s_%s", clu[["clus_var_name"]], dis)]
  SLL2_rel@sam_data <- SLL2_rel@sam_data[ , colnames(SLL2_rel@sam_data) != sprintf("%s_%s", clu[["clus_var_name"]], dis)]
  
  # add stomatotype values
  SLL2@sam_data[ names(clu[["cluster_full"]]), sprintf("%s_%s",clu[["clus_var_name"]], dis) ] <- as.factor(clu[["cluster_full"]])
  SLL2_rel@sam_data[ names(clu[["cluster_full"]]), sprintf("%s_%s",clu[["clus_var_name"]], dis) ] <- as.factor(clu[["cluster_full"]])
}

SLL2.meta <- meta(SLL2)


saveRDS(SLL2,      sprintf("%s/R_objects/SLL2.rds", p2_dir))
saveRDS(SLL2_rel,  sprintf("%s/R_objects/SLL2_rel.rds", p2_dir))
saveRDS(SLL2.meta, sprintf("%s/R_objects/SLL2.meta.rds", p2_dir))









# ******************************************************* #
# Get lists of taxa that are significantly greater in each Stomatotype of each type ####
# ******************************************************* #


tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2))

sig_tax.dists <- list()
for(dist_meas in c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
                   "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
                   "Bray","Jaccard","Canberra")) {
  print(dist_meas)
  sig_tax.dists[[ dist_meas ]] <- get_sig_tax(tagl, dist_meas, ps.dists, gloms_rel)
}

# add a list for 3 clusters with JSD since that seems potentially more suitable at the moment
# as well as 3 and 4 for Weighted_Unifrac
sig_tax.dists[[ "JSD.3" ]] <- get_sig_tax(tagl, "JSD", ps.dists, gloms_rel, ncf = 3)
# sig_tax.dists[[ "Weighted_Unifrac.2" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ncf = 2)
sig_tax.dists[[ "Weighted_Unifrac.3" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ps.dists, gloms_rel, ncf = 3)
sig_tax.dists[[ "Weighted_Unifrac.4" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ps.dists, gloms_rel, ncf = 4)
sig_tax.dists[[ "VAW_GUnifrac.3" ]] <- get_sig_tax(tagl, "VAW_GUnifrac", ps.dists, gloms_rel, ncf = 3)


sig_tax.dists[[ "Aitchison" ]] <- get_sig_tax(phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)), 
                                              "Aitchison", ps.dists, gloms_clr)
# ******************************************************* #






# Do the same, but then filter the taxa in each Stomatotype in order to determine which
#   is the most significant Stomatotype for a given taxon in the case that a taxon
#   appears as significant in more than 1 Stomatotype. 
#   (Still possible to have multiple Stomatotypes if not a significant difference between them)
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2))

sig_tax.anova.dists <- list()
for(dist_meas in c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
                   "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
                   "Bray","Jaccard","Canberra")) {
  print(dist_meas)
  sig_tax.anova.dists[[ dist_meas ]] <- get_sig_tax(tagl, dist_meas, ps.dists, gloms_rel, filt_anova = TRUE)
}

# add a list for 3 clusters with JSD since that seems potentially more suitable at the moment
# as well as 3 and 4 for Weighted_Unifrac
sig_tax.anova.dists[[ "JSD.3" ]] <- get_sig_tax(tagl, "JSD", ps.dists, gloms_rel, ncf = 3, filt_anova = TRUE)
# sig_tax.anova.dists[[ "Weighted_Unifrac.2" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ncf = 2, filt_anova = TRUE)
sig_tax.anova.dists[[ "Weighted_Unifrac.3" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ps.dists, gloms_rel, ncf = 3, filt_anova = TRUE)
sig_tax.anova.dists[[ "Weighted_Unifrac.4" ]] <- get_sig_tax(tagl, "Weighted_Unifrac", ps.dists, gloms_rel, ncf = 4, filt_anova = TRUE)
sig_tax.anova.dists[[ "VAW_GUnifrac.3" ]] <- get_sig_tax(tagl, "VAW_GUnifrac", ps.dists, gloms_rel, ncf = 3, filt_anova = TRUE)


sig_tax.anova.dists[[ "Aitchison" ]] <- get_sig_tax(phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)), 
                                                    "Aitchison", ps.dists, gloms_clr, filt_anova = TRUE)
# ******************************************************* #

for (na in names(sig_tax.anova.dists)) {
  cat(c(na, "\t", length(unlist(sig_tax.anova.dists[[na]])), length(unique(unlist(sig_tax.anova.dists[[na]]))),"\n"))
}
# according to this only JSD.3 has genera that still appear in more than 1 Stomatotype
#   JSD.3 => Tannerella in Stomatotypes 1 and 2 (3 total)
#   the p-value in the anova between Stomatotypes 1 and 2 for the given genus was ~0.98

# The comments below were true before filtering out samples with very low counts...only affected Unweighted Unifrac (reasonable)
# # # according to this, only Unweighted_Unifrac and JSD.3 have genera that still appear in more than 1 Stomatotype
# # #   Unweighted_Unifrac => Solobacterium in Stomatotypes 1 and 2 (3 total)
# # #   JSD.3 => Tannerella in Stomatotypes 1 and 2 (3 total)
# # # In both cases, the p-value in the anova between Stomatotypes 1 and 2 for the given genus was ~0.99











# ********************************************************************************************** #
# get list of consensus clusters based on different distance measures ####

# ******************************************************* #
get_stom_dists <- function(sig_tax, dis_list, minPer = 0.33) {
  sto.dist <- list()
  
  for (dis in dis_list) {
    # select first dis for comparison to all others
    sig.tax <- sig_tax[[ dis ]]
    sto.dist[[ dis ]] <- sapply(dists[ dis != dists ], function(dis2) {
      # for each of the remaining dis, compare to the given dis
      sig.tax2 <- sig_tax[[ dis2 ]]
      
      # make table in which rows are the Stomatotypes for a given dis, columns are the other dis
      #   values are the stomatotype for a given dis that is most similar to given row (similarity of at least minPer)
      st.mat <- sapply(names(sig.tax), function(st) {
        sapply(names(sig.tax2), function(st2) {
          # get percentage of stomatotype from one measure that matches stomatotype of another measure
          length(sig.tax[[st]][ sig.tax[[st]] %in% sig.tax2[[st2]] ]) / length(sig.tax[[st]])
        })
      })
      apply(st.mat, 2, function(x) ifelse(max(x) > minPer,
                                          names(which(x == max(x))),
                                          NA))
      # {NA;print(c(max(x),dis,dis2))}))
    })
  }
  
  return(sto.dist)
}
# ******************************************************* #

dists <- c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")
# dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
#            "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
#            "Bray","Jaccard","Canberra",
#            "JSD.3","Weighted_Unifrac.3","Weighted_Unifrac.4","VAW_GUnifrac.3")
#            # "JSD.3","Weighted_Unifrac.2","Weighted_Unifrac.4")

stoms.dists <- get_stom_dists(sig_tax.dists, dists)
stoms.anova.dists <- get_stom_dists(sig_tax.anova.dists, dists)











# ******************************************************* #

# # Since all 5 measures suggest 2 Stomatotypes, will now create 2 consensus Stomatotypes,
# #   and any samples that fall within non-corresponding stomatotypes (as determined in above section)
# #   will be labeled NA
# 
# s1.samps <- sample_names(subset_samples(SLL2, Stomatotype_JSD=="1" 
#                                         # & Stomatotype_JSD_alt=="1"
#                                         # & Stomatotype_Weighted_Unifrac=="2"
#                                         & (Stomatotype_Weighted_Unifrac=="1" | Stomatotype_Weighted_Unifrac=="2")
#                                         # & Stomatotype_Weighted_Unifrac.2=="2"
#                                         # & Stomatotype_Bray=="1"
#                                         # & Stomatotype_Canberra=="1" 
#                                         # & Stomatotype_Unweighted_Unifrac=="1"
#                                         ))
# 
# s2.samps <- sample_names(subset_samples(SLL2, Stomatotype_JSD=="2" 
#                                         # & Stomatotype_JSD_alt=="2" 
#                                         & Stomatotype_Weighted_Unifrac=="3"
#                                         # & Stomatotype_Weighted_Unifrac.2=="1"
#                                         # & Stomatotype_Bray=="2"
#                                         # & Stomatotype_Canberra=="2"
#                                         # & Stomatotype_Unweighted_Unifrac=="3"
#                                         ))
# 
# noStom.samps <- sample_names(SLL2)[ ! sample_names(SLL2) %in% c(s1.samps,s2.samps)]
# 
# 
# stoms <- sapply(sample_names(SLL2), function(x) {
#   if (x %in% s1.samps) 1
#   else if (x %in% s2.samps) 2
#   else "None"
# })
# table(stoms)





# ******************************************************* #

# Will also now create 3 consensus Stomatotypes, since JSD and Weighted_Unifrac may actually suggest 3

# s1.samps <- sample_names(subset_samples(SLL2, Stomatotype_JSD.3=="1"
#                                         # & Stomatotype_Weighted_Unifrac=="1"
#                                         & Stomatotype_Weighted_Unifrac.3=="1"
#                                         # & Stomatotype_Weighted_Unifrac.4=="1"
# ))
# 
# s2.samps <- sample_names(subset_samples(SLL2, Stomatotype_JSD.3=="2"
#                                         # & Stomatotype_Weighted_Unifrac=="2"
#                                         & Stomatotype_Weighted_Unifrac.3=="2"
#                                         # & Stomatotype_Weighted_Unifrac.4=="2"
# ))
# 
# s3.samps <- sample_names(subset_samples(SLL2, Stomatotype_JSD.3=="3"
#                                         # & Stomatotype_Weighted_Unifrac=="3"
#                                         & Stomatotype_Weighted_Unifrac.3=="3"
#                                         # & Stomatotype_Weighted_Unifrac.4=="4"
# ))
# 
# noStom.samps <- sample_names(SLL2)[ ! sample_names(SLL2) %in% c(s1.samps, s2.samps, s3.samps)]



s1.samps <- sample_names(subset_samples(SLL2, Stomatotype_Aitchison=="1"
                                        # & Stomatotype_Weighted_Unifrac=="1"
                                        & Stomatotype_Bray=="1"
                                        # & Stomatotype_Weighted_Unifrac.4=="1"
))

s2.samps <- sample_names(subset_samples(SLL2, Stomatotype_Aitchison=="2"
                                        # & Stomatotype_Weighted_Unifrac=="2"
                                        & Stomatotype_Weighted_Unifrac=="2"
                                        & Stomatotype_Bray=="2"
                                        # & Stomatotype_Weighted_Unifrac.4=="2"
))

noStom.samps <- sample_names(SLL2)[ ! sample_names(SLL2) %in% c(s1.samps, s2.samps)]


stoms <- sapply(sample_names(SLL2), function(x) {
  if (x %in% s1.samps) 1
  else if (x %in% s2.samps) 2
  # else if (x %in% s3.samps) 3
  else "None"
})
table(stoms)














# ************************************************* #
# Add to list of significantly different genera for consensus Stomatotype samples ####
# just using a modified version of the get_sig_tax() function above

# include "None" samples - do not follow consensus of clustering by indicated methods
stoms.consensus <- sapply(unique(stoms), function(x) names(stoms)[ stoms == x ])
# tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2))
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2))
top.tax <- names(sort(taxa_sums(tagl), decreasing = T))

sig_tax.dists[[ "consensus" ]] <- sapply(names(stoms.consensus), function(st) {
  other_st <- names(stoms.consensus)[ st != names(stoms.consensus) ]
  
  sig <- sapply(top.tax, function(tax) {
    wilcox.test(as.numeric(as.matrix(otu_table(tagl)[ tax, stoms.consensus[[ st ]] ])), 
                as.numeric(as.matrix(otu_table(tagl)[ tax, unlist(stoms.consensus[ other_st ]) ])),
                alternative = "greater")$p.value
  })
  sig.adj <- p.adjust(sig, method = "fdr")
  names(sig.adj[sig.adj<0.05])
})

# also filtered version with anova
sig_tax.anova.dists[[ "consensus" ]] <- filter_sig_tax_w_anova(sig_tax.dists[[ "consensus" ]], 
                                                               stoms, 
                                                               tagl)




# include only samples that follow consensus of clustering by indicated methods
stoms.consensus_only <- sapply(unique(stoms[stoms != "None"]), function(x) names(stoms)[ stoms == x ])
# tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2))
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2))
top.tax <- names(sort(taxa_sums(tagl), decreasing = T))

sig_tax.dists[[ "consensus_only" ]] <- sapply(names(stoms.consensus_only), function(st) {
  other_st <- names(stoms.consensus_only)[ st != names(stoms.consensus_only) ]
  
  sig <- sapply(top.tax, function(tax) {
    wilcox.test(as.numeric(as.matrix(otu_table(tagl)[ tax, stoms.consensus_only[[ st ]] ])), 
                as.numeric(as.matrix(otu_table(tagl)[ tax, unlist(stoms.consensus_only[ other_st ]) ])),
                alternative = "greater")$p.value
  })
  sig.adj <- p.adjust(sig, method = "fdr")
  names(sig.adj[sig.adj<0.05])
})
# also filtered version with anova
sig_tax.anova.dists[[ "consensus_only" ]] <- filter_sig_tax_w_anova(sig_tax.dists[[ "consensus_only" ]], 
                                                                    stoms[stoms != "None"], 
                                                                    prune_samples(names(stoms[stoms != "None"]), tagl))
# ************************************************* #

# just checking again
for (na in names(sig_tax.anova.dists)) {
  print(c(na, length(unlist(sig_tax.anova.dists[[na]])), length(unique(unlist(sig_tax.anova.dists[[na]])))))
}
# neither "consensus" nor "consensus_only" have genera that still appear in more than 1 Stomatotype












# ***************************************************************************************** #
# Plot 1: Drivers ####
library(gridExtra)
library(grid)

plot_drivers <- function(tagl, tl, dist_meas, ncf=NULL) {
  
  clu <- get_clusters( tagl, tl, dist_meas, ps.dists, gloms_clr, nClustForce = ncf, CH.by.tagl = FALSE )
  
  if (clu[[ "nclusters" ]] == 2) { 
    clust.col <- c("indianred1", "turquoise3")
  } else if (clu[[ "nclusters" ]] == 3) { 
    clust.col <- c("indianred1", "green4", "dodgerblue")
  } else if (clu[[ "nclusters" ]] == 4) { 
    clust.col <- c("indianred1", "yellowgreen", "turquoise3", "purple1")
  } else if (clu[[ "nclusters" ]] == 5) {
    clust.col <- c("indianred1", "#829311", "seagreen3", "dodgerblue", "orchid")
  }
  
  # ************************* #
  ## between-class analysis (BCA)
  k <- 10
  pca = dudi.pca(t(clu[[ "tax.to.clust" ]]), scannf=F, nf=k)
  bet = bca(pca, fac=as.factor(clu[[ "cluster_full" ]]), scannf=F, nf=k-1)
  
  # ************************* #
  stomato.leaders <- NULL
  for (i in 1:clu[[ "nclusters" ]]){
    stomato.leaders[i] <- colnames(bet$tab)[bet$tab[i,]==max(bet$tab[i,])]
  }
  
  # ************************* #
  # stomato.leaders <- c("Porphyromonas","Prevotella","Pseudomonas")
  stomato.melts <- list()
  for (s in stomato.leaders) {
    stomato.melts[[ s ]] <- cbind( reshape2::melt(t(clu[[ "tax.to.clust" ]][ s , ])), 
                                   reshape2::melt(as.matrix(clu[[ "cluster_full" ]])) )
    stomato.melts[[ s ]] <- stomato.melts[[ s ]][,c(1,2,3,6)]
    colnames(stomato.melts[[ s ]])[4] <- 'cluster_full'
  }
  
  # ************************* #
  
  # to give the same scale in each boxplot
  ymax <- max(clu[[ "tax.to.clust" ]][stomato.leaders,])
  
  stomato.box.leaders <- list()
  for (s in stomato.leaders) {
    stomato.box.leaders[[ s ]] <- ggplot(stomato.melts[[ s ]], 
                                         aes(x=factor(cluster_full, levels=c(1:clu[[ "nclusters" ]])), y=value, 
                                             fill=factor(cluster_full, levels=c(1:clu[[ "nclusters" ]])))) +
      geom_boxplot(notch = T) +
      # geom_violin(scale = "count", trim = F, draw_quantiles = c(0.25, 0.5, 0.75)) +
      # geom_point() +
      ggtitle( s ) + ylim(0, ymax) +
      xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
  }
  
  grid.arrange(grobs = stomato.box.leaders,
               top = textGrob(sprintf("Drivers at %s level for %s distances", tl, dist_meas),
                              gp = gpar(fontsize=15)),
               ncol = length( stomato.box.leaders ))
}
# ******************************* #

# tl <- "Species"; tagl <- SLL2_rel
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2_rel))

plot_drivers(tagl, tl, "JSD")
plot_drivers(tagl, tl, "JSD", ncf = 3)
plot_drivers(tagl, tl, "Weighted_Unifrac")
# plot_drivers(tagl, tl, "Weighted_Unifrac", ncf = 2)
plot_drivers(tagl, tl, "Weighted_Unifrac", ncf = 3)
plot_drivers(tagl, tl, "Weighted_Unifrac", ncf = 4)
plot_drivers(tagl, tl, "VAW_GUnifrac")
plot_drivers(tagl, tl, "VAW_GUnifrac", ncf = 3)

# ***************************************************************************************** #












# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####







# ***************************************************************************************** #
# Plot 2: principal coordinates analysis (PCoA) ####
library(scales)
library(viridis)

# tl <- "Species"; tagl <- SLL2_rel
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2_rel))

pcoas <- list()
for (dis in c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
              "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
              "Bray","Jaccard","Canberra","Aitchison")) {
  print(dis)
  clu <- get_clusters( tagl, tl, dist_meas = dis, ps.dists )
  # pcoas[[ dis ]] <- dudi.pco(clu[[ "dist.matrix" ]], scannf=F, nf=3)
  pcoas[[ dis ]] <- ape::pcoa(clu[[ "dist.matrix" ]])
}

pcoas[[ "JSD.3" ]] <- pcoas[[ "JSD" ]]
# pcoas[[ "Weighted_Unifrac.2" ]] <- pcoas[[ "Weighted_Unifrac" ]]
pcoas[[ "Weighted_Unifrac.3" ]] <- pcoas[[ "Weighted_Unifrac" ]]
pcoas[[ "Weighted_Unifrac.4" ]] <- pcoas[[ "Weighted_Unifrac" ]]
pcoas[[ "VAW_GUnifrac.3" ]] <- pcoas[[ "VAW_GUnifrac" ]]


col.to.plot <- sprintf("Stomatotype_%s", dist_meas)
ade4::s.class(pcoas[[ dist_meas ]]$vectors, 
              fac=SLL2.meta[ rownames(pcoas[[ dist_meas ]]$vectors), col.to.plot ],
              grid=F, clabel=1.5, 
              col=hue_pal()(length(unique(SLL2.meta[ rownames(pcoas[[ dist_meas ]]$vectors), col.to.plot ]))), 
              sub = sprintf("PCoA separated by %s", col.to.plot ))



# col.to.plot <- "Sibling_unit"
# fa <- as.factor(SLL2.meta[ rownames(pcoas[[ dist_meas ]]$li), col.to.plot ])
# ade4::s.class(pcoas[[ dist_meas ]]$li[ rownames(pcoas[[ dist_meas ]]$li)[fa != "None"],  ], 
#               fac=fa[fa != "None"],
#               grid=F, clabel=1.5, 
#               col=hue_pal()(length(unique(fa[ fa != "None"]))), 
#               sub = sprintf("PCoA separated by %s", col.to.plot ))

# ***************************************************************************************** #




# make PCoAs for stoms.consensus_only samples
consOnly.samps <- unname( unlist(stoms.consensus_only) )



# # *********************************************************** #
# dists.consensus_only <- list()
# 
# phy <- prune_samples(consOnly.samps, SLL2)
# 
# # Unifracs
# dists.consensus_only[[ "Weighted_Unifrac" ]] <- UniFrac(phy, weighted = T, parallel = T)
# dists.consensus_only[[ "Unweighted_Unifrac" ]] <- UniFrac(phy, weighted = F, parallel = T)
# # bray and jaccard
# dists.consensus_only[[ "Bray" ]] <- vegdist(t(phy@otu_table), method = "bray")
# dists.consensus_only[[ "Jaccard" ]] <- vegdist(decostand(as.data.frame(t(phy@otu_table)), method="pa"), method = "jaccard")
# # Aitchison
# f.CO <- codaSeq.filter(phy@otu_table[ , consOnly.samps ],
#                        min.reads=950, # filter out samples with fewer reads
#                        min.prop=0.001, # filter out taxa with lower minimum abundance
#                        min.occurrence=0.05, # filter out taxa not appearing in at least this prop of samps
#                        samples.by.row=FALSE)
# # replace 0 values with an estimate
# f.n0.CO <- cmultRepl(t(f.CO), method="CZM", label=0)
# dists.consensus_only[[ "Aitchison" ]] <- aDist(f.n0.CO)
# 
# saveRDS(dists.consensus_only, sprintf("%s/R_objects/dists.consensus_only.rds", p2_dir))

dists.consensus_only <- readRDS(sprintf("%s/R_objects/dists.consensus_only.rds", p2_dir))
# *********************************************************** #



# *********************************************************** #
pcoas.consensus_only <- list()

for (dist_meas in names(dists.consensus_only)) {
  pcoas.consensus_only[[ dist_meas ]] <- ape::pcoa(dists.consensus_only[[ dist_meas]])
}

# *********************************************************** #


dist_meas <- "Aitchison"
dist_meas <- "Weighted_Unifrac"
dist_meas <- "Unweighted_Unifrac"
dist_meas <- "Bray"
dist_meas <- "Jaccard"

col.to.plot <- sprintf("Stomatotype_%s", dist_meas)

ade4::s.class(pcoas.consensus_only[[ dist_meas ]]$vectors, 
              fac=as.factor(unname(unlist(sapply(names(stoms.consensus_only), function(x) rep(x, length(stoms.consensus[[ x]])))))),
              grid=F, clabel=1.5, 
              col=hue_pal()(length(stoms.consensus)), 
              sub = sprintf("PCoA separated by %s", col.to.plot ))
# ***************************************************************************************** #



















# ***************************************************************************************** #
# Plot GRADIENTS of abundances of particular taxa or groups of taxa ####
library(RColorBrewer)

# tl <- "Species"; tagl <- SLL2_rel
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2_rel))

# ************************************* #

plot_gradients(tl, tagl, "JSD", "Stom 1 sig")
plot_gradients(tl, tagl, "JSD.3", "Stom 1 sig")
plot_gradients(tl, tagl, "Weighted_Unifrac", "Stom 1 sig")
# plot_gradients(tl, tagl, "Weighted_Unifrac.2", "Stom 1 sig")
plot_gradients(tl, tagl, "Weighted_Unifrac.3", "Stom 1 sig")
plot_gradients(tl, tagl, "Weighted_Unifrac.4", "Stom 1 sig")
plot_gradients(tl, tagl, "VAW_GUnifrac.3", "Stom 1 sig")




# ************************************* #
tl <- "Genus"
nTop <- 5

# ************************************* #
# Aitchison
plot_gradients(tl, phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)), 
               "Aitchison", pcoas, gloms_clr, SLL2.meta,
               "Stom_1_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)), 
               "Aitchison", pcoas, gloms_clr, SLL2.meta, 
               "Stom_2_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop=nTop))

plot_gradients(tl, phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)),
               "Aitchison", pcoas, gloms_clr, SLL2.meta, "Stom 1 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_clr[[ tl ]], taxa_are_rows = T), sample_data(SLL2)), 
               "Aitchison", pcoas, gloms_clr, SLL2.meta, "Stom 2 sig")



# ************************************* #
# Weighted_Unifrac
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Weighted_Unifrac", pcoas, gloms_rel, SLL2.meta, 
               "Stom_1_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Weighted_Unifrac", gloms_rel, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Weighted_Unifrac", pcoas, gloms_rel, SLL2.meta, 
               "Stom_2_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Weighted_Unifrac", gloms_rel, SLL2.meta, nTop=nTop))

plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)),
               "Weighted_Unifrac", pcoas, gloms_rel, SLL2.meta, "Stom 1 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Weighted_Unifrac", pcoas, gloms_rel, SLL2.meta, "Stom 2 sig")



# ************************************* #
# Unweighted_Unifrac
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Unweighted_Unifrac", pcoas, gloms_rel, SLL2.meta, 
               "Stom_1_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Unweighted_Unifrac", gloms_rel, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Unweighted_Unifrac", pcoas, gloms_rel, SLL2.meta, 
               "Stom_2_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Unweighted_Unifrac", gloms_rel, SLL2.meta, nTop=nTop))

plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)),
               "Unweighted_Unifrac", pcoas, gloms_rel, SLL2.meta, "Stom 1 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Unweighted_Unifrac", pcoas, gloms_rel, SLL2.meta, "Stom 2 sig")



# ************************************* #
# Bray
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Bray", pcoas, gloms_rel, SLL2.meta, 
               "Stom_1_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Bray", pcoas, gloms_rel, SLL2.meta, 
               "Stom_2_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Bray", pcoas, gloms_rel, SLL2.meta, 
               "Stom_3_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop=nTop))

plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)),
               "Bray", pcoas, gloms_rel, SLL2.meta, "Stom 1 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Bray", pcoas, gloms_rel, SLL2.meta, "Stom 2 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Bray", pcoas, gloms_rel, SLL2.meta, "Stom 3 sig")



# ************************************* #
# Jaccard
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Jaccard", pcoas, gloms_rel, SLL2.meta, 
               "Stom_1_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Jaccard", gloms_rel, SLL2.meta, nTop=nTop))
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Jaccard", pcoas, gloms_rel, SLL2.meta, 
               "Stom_2_drivers", save.grads = T, plot_dirs = "Driver_gradients",
               drivers = get_drivers("Genus", "Jaccard", gloms_rel, SLL2.meta, nTop=nTop))

plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)),
               "Jaccard", pcoas, gloms_rel, SLL2.meta, "Stom 1 sig")
plot_gradients(tl, phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2)), 
               "Jaccard", pcoas, gloms_rel, SLL2.meta, "Stom 2 sig")
# ***************************************************************************************** #














# ***************************************************************************************** #
# Check drivers and plot GRADIENTS ####

# ************************************* #


get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop = 5)

get_drivers("Genus", "Weighted_Unifrac", gloms_rel, SLL2.meta, nTop = 5)
get_drivers("Genus", "Unweighted_Unifrac", gloms_rel, SLL2.meta, nTop = 5)

get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop = 5)
get_drivers("Genus", "Jaccard", gloms_rel, SLL2.meta, nTop = 5)


library(venn)

nTop <- 10
venn(list("Aitchison_1"=names(get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop = nTop)[[1]]),
          "Weighted_Unifrac_2"=names(get_drivers("Genus", "Weighted_Unifrac", gloms_rel, SLL2.meta, nTop = nTop)[[2]]),
          "Jaccard_1"=names(get_drivers("Genus", "Jaccard", gloms_rel, SLL2.meta, nTop = nTop)[[1]]),
          "Unweighted_Unifrac_1"=names(get_drivers("Genus", "Unweighted_Unifrac", gloms_rel, SLL2.meta, nTop = nTop)[[1]]),
          "Bray_1"=names(get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop = nTop)[[1]])),
     ilabels = T, zcolor = "style", counts = T, sncs = 1.5, ilcs = 1.5)

venn(list("Aitchison_2"=names(get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop = nTop)[[2]]),
          "Weighted_Unifrac_1"=names(get_drivers("Genus", "Weighted_Unifrac", gloms_rel, SLL2.meta, nTop = nTop)[[1]]),
          "Jaccard_2"=names(get_drivers("Genus", "Jaccard", gloms_rel, SLL2.meta, nTop = nTop)[[2]]),
          "Unweighted_Unifrac_1"=names(get_drivers("Genus", "Unweighted_Unifrac", gloms_rel, SLL2.meta, nTop = nTop)[[1]]),
          "Bray_3"=names(get_drivers("Genus", "Bray", gloms_rel, SLL2.meta, nTop = nTop)[[3]])),
     ilabels = T, zcolor = "style", counts = T, sncs = 1.5, ilcs = 1.5)
# title("All samples - nonsynonymous SNPs", line = -1)
# ***************************************************************************************** #





























# ****************************************************************************************************************** #
# *********************************************************** #
# PERMANOVA based on various categorical variables using different distance matrices ####
# *********************************************************** #

# *********************************************************** #
get_permanova_ps <- function(groups, phy, useStrata=FALSE) {
  
  b_dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac",
               "VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac",
               "Bray","Jaccard","Canberra")
  bdist.codes <- structure(c("jsd","weighted_Unifrac","unweighted_Unifrac",
                             "guni.VAW","guni.a0","guni.a05",
                             "bray","jaccard","canberra"), .Names = b_dists)
  
  pnov.p <- pnov.R <- matrix(NA, nrow=length(b_dists), ncol = length(groups))
  rownames(pnov.p) <- rownames(pnov.R) <- b_dists
  colnames(pnov.p) <- colnames(pnov.R) <- groups
  
  for (b in rownames(pnov.p)) {
    print(b)
    
    for (g in colnames(pnov.p)) {
      print(sprintf("    %s",g))
      remove_NAs <- ! is.na(get_variable(phy, g))
      if (sum(remove_NAs) > 0) {
        # only run if not all NAs for given g
        phy.ad <- prune_samples(remove_NAs, phy)
        # in case there were NAs, must use only the rows and columns of the distance matrix that are not NAs here
        sn <- sample_names(phy.ad)
        # strataNames <- unname(sapply(sample_names(phy.ad), function(x) 
        #   ifelse(grepl('CTRL',x), x, strsplit(x,'\\.')[[1]][1])))
        
        if (length(table(phy.ad@sam_data[,g])) > 1) {
          # only run when there is more than 1 group present for the given categorical variable
          # if (sum(duplicated(strataNames)) > 0 & useStrata==TRUE) {
          #   # in sample groups that include the T0 and T3 pairs, use strata to constrain permutations
          #   ad <- adonis(as.formula(sprintf('distance(phy.ad, method = "%s") ~ %s', bdist.codes[b], g)),
          #                data=as(sample_data(phy.ad), "data.frame"), 
          #                strata = strataNames)
          # } else {
          #   ad <- adonis(as.formula(sprintf('distance(phy.ad, method = "%s") ~ %s', bdist.codes[b], g)),
          #                data=as(sample_data(phy.ad), "data.frame"))
          # }
          ad <- adonis(as.formula(sprintf('as.dist(%s[sn, sn]) ~ %s', bdist.codes[b], g)),
                       data=as(sample_data(phy.ad), "data.frame"))
          
          pnov.p[b, g] <- ad$aov.tab$`Pr(>F)`[1]
          pnov.R[b, g] <- ad$aov.tab$R2[1]
          
          # ad <- adonis2(as.formula(sprintf('distance(phy.ad, method = "%s") ~ %s', bdist.codes[b], g)),
          #               data=as(sample_data(phy.ad), "data.frame"))
          # pnov.p[b, g] <- ad$`Pr(>F)`[1]
          # pnov.R[b, g] <- ad$R2[1]
        }
      }
    }
  }
  
  pnov.p.adj <- apply(pnov.p, 2, p.adjust, method='bonferroni')
  
  pnov.p.adj[ is.nan(pnov.p.adj) ] <- 1
  pnov.p.adj[ is.na(pnov.p.adj) ] <- 1
  
  #at least 1 good p value within samples
  mins <- apply(pnov.p.adj, 2, min)
  goodPs <- mins[mins < 0.05] #-log10(0.05) == 1.3013
  #at least 1 good p within taxa
  tmins <- apply(pnov.p.adj, 1, min)
  tgoodPs <- tmins[tmins < 0.05]
  
  # get tables of only significant values
  pnov.p.adj.toWrite <- pnov.p.adj[ names(tgoodPs), names(goodPs) ]
  pnov.p.adj.toWrite[ pnov.p.adj.toWrite >= 0.05 ] <- ""
  
  pnov.R.toWrite <- pnov.R[ names(tgoodPs), names(goodPs) ]
  pnov.R.toWrite[ pnov.p.adj.toWrite == "" ] <- ""
  
  # write.csv(pnov.p.adj.toWrite,
  #           sprintf("%s/permanova/signif_permanova.csv", p1_dir))
  
  return(list(R2 = t(pnov.R), R2.toWrite = t(pnov.R.toWrite), p = t(pnov.p.adj), p.toWrite = t(pnov.p.adj.toWrite)))
}
# *********************************************************** #



groupQs.permanova <- c("Gender","Age_groups","Municipal_zone","Moisture_in_home",
                       "Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds",
                       "Pets.Reptiles_amphibians","Pets.Fish","Pets.Rabbits","Pets.Rodents","Pets.Mammals",
                       "Smoker","Water_type_home","Braces","Mouth_piercing","Fluoride_toothpaste",
                       "Fluoride_supplement","Mouth_wounds","Reason_dental_visit",
                       "Chronic_disorder","Celiac","Cystic_fibrosis","Gingivitis_periodontitis",
                       "Downs_Syndrome","Eating_disorder","Other_disorder_binary",
                       "Medications","Antibiotics","Analgesics","Vitamin_supplements",
                       "Other_medications_binary","Asthma","Wheezing","How_do_you_feel","Do_you_feel_well",additional_diseases,
                       "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen",
                       "Allergy.Animals","Allergy.Sun","Allergy.Medications","Allergy.Nickel",
                       "Allergy.Stings","Allergy.Latex","Allergy.Anisakis","Allergy.Seasonal",
                       "Allergy.other_binary","Bite_nails","Hair_in_mouth","Chew_pens",
                       "Wash_hands_before_eat","Wash_hands_after_bathroom","Kissing_partner",
                       "BMI_group","BMI_official",group_water_data,
                       "MALDI.Yeast_detected","MALDI.Mold_detected",
                       "MALDI.Bacteria_detected","MALDI.Yeast.Candida_albicans",
                       "MALDI.Yeast.Candida_guillermondii","MALDI.Yeast.Candida_parapsilosis",
                       "MALDI.Bacteria.Pseudomonas_putida","MALDI.Bacteria.Ralstonia_insidiosa",
                       "Consumption.Milk.Binary","Consumption.Yogurt.Binary","Consumption.Sweets.Binary",
                       "Consumption.Chewing_gum.Binary","Consumption.Nuts.Binary",
                       "Drinks.Decaf_coffee.Binary","Drinks.Coffee.Binary","Drinks.Tea.Binary",
                       "Drinks.Infusion.Binary","Drinks.Soda.Binary","Drinks.Soda_sugarless.Binary",
                       "Drinks.Soda_decaf.Binary","Drinks.RedBull.Binary","Drinks.Other_sugary_drinks.Binary",
                       "Drinks.Alcohol_cold.Binary","Drinks.Alcohol_hot.Binary",
                       "Brushing.Binary","Floss.Binary","Last_dental_visit.Binary",
                       "City","Province","Community",
                       "Diversity_group_Div.Shannon","Diversity_group_Div.Simpson",
                       "Diversity_group_Weighted_Unifrac","Diversity_group_Unweighted_Unifrac",
                       "Diversity_group_Faiths.PD","Diversity_group_Species_Richness",
                       "Diversity_group_Bray.Curtis","Diversity_group_Canberra",
                       "Stomatotype_JSD","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
                       "Stomatotype_VAW_GUnifrac","Stomatotype_a0_GUnifrac","Stomatotype_a05_GUnifrac",
                       "Stomatotype_Bray","Stomatotype_Jaccard","Stomatotype_Canberra",
                       "Stomatotype_Aitchison",
                       "Stomatotype_JSD.3","Stomatotype_Weighted_Unifrac.3",
                       "Stomatotype_Weighted_Unifrac.4","Stomatotype_VAW_GUnifrac.3")

# permanovas.2 <- get_permanova_ps(groupQs.permanova, SLL2)
# 
# # write.csv(permanovas.2[["p.toWrite"]], sprintf("%s/permanova/signif_permanova.pvals.csv", p2_dir))
# # write.csv(permanovas.2[["R2.toWrite"]], sprintf("%s/permanova/signif_permanova.R2.csv", p2_dir))
# 
# saveRDS(permanovas.2, file = sprintf("%s/R_objects/permanovas.2.rds", p2_dir))
# 
# # **************** #

permanovas.2 <- readRDS(sprintf("%s/R_objects/permanovas.2.rds", p2_dir))

# *********************************************************** #







# *********************************************************** #

# ... For families ####


weighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))
diag(weighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
unweighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))
diag(unweighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored

bray <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_bray.rds", p2_dir))
diag(bray) <- NA
jaccard <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_jaccard.rds", p2_dir))
diag(jaccard) <- NA

aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

full_ords <- list("Aitchison"=aitch, "Weighted_Unifrac"=as.matrix(weighted_Unifrac),
                  "Unweighted_Unifrac"=as.matrix(unweighted_Unifrac), "Bray"=bray,
                  "Jaccard"=jaccard)


# PERMANOVA
# fam.adonis_list <- list()
# for (dist_meas in c("weighted_Unifrac","unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "bray","jaccard","aitch")) {
#   
#   fam.adonis_list[[ dist_meas ]] <- list()
#   
#   for (unit_lab in c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
#                      "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
#                      "Grandparent_Grandchild_unit")) {
#     print(c(dist_meas, unit_lab))
#     fam.adonis_list[[ dist_meas ]][[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(%s) ~ %s', 
#                                                                               dist_meas, unit_lab)), 
#                                                            strata = SLL2.meta[,"seqGroup"],
#                                                            data=SLL2.meta)
#   }
# }
# 
# saveRDS(fam.adonis_list, sprintf("%s/R_objects/fam.adonis_list.rds", p2_dir))


fam.adonis_list <- readRDS(sprintf("%s/R_objects/fam.adonis_list.rds", p2_dir))

# for (unit_lab in c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
#                    "Parent_Child_unit","Mother_Child_unit","Father_Child_unit","Grandparent_Grandchild_unit")) {
#   adonis_list$Weighted_Unifrac[[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(weighted_Unifrac) ~ %s', unit_lab)), 
#                                                        data=as(sample_data(SLL2), "data.frame"))
#   adonis_list$Unweighted_Unifrac[[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(unweighted_Unifrac) ~ %s', unit_lab)), 
#                                                          data=as(sample_data(SLL2), "data.frame"))
#   adonis_list$JSD[[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(jsd) ~ %s', unit_lab)), 
#                                           data=as(sample_data(SLL2), "data.frame"))
#   adonis_list$Bray.Curtis[[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(bray) ~ %s', unit_lab)), 
#                                                   data=as(sample_data(SLL2), "data.frame"))
#   adonis_list$Canberra[[ unit_lab ]] <- adonis(as.formula(sprintf('as.dist(canberra) ~ %s', unit_lab)), 
#                                                data=as(sample_data(SLL2), "data.frame"))
# }
# adonis(as.dist(weighted_Unifrac) ~ Family_unit, data=as(sample_data(SLL2), "data.frame"))
# adonis(as.dist(unweighted_Unifrac) ~ Family_unit, data=as(sample_data(SLL2), "data.frame"))
# 
# adonis(as.dist(weighted_Unifrac) ~ Sibling_unit, data=as(sample_data(SLL2), "data.frame"))
# adonis(as.dist(unweighted_Unifrac) ~ Sibling_unit, data=as(sample_data(SLL2), "data.frame"))

# ANOSIM
# fam.anosim_list <- list()
# for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "Bray","Jaccard")) {
#   
#   fam.anosim_list[[ dist_meas ]] <- list()
#   
#   for (unit_lab in c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
#                      "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
#                      "Grandparent_Grandchild_unit")) {
#     print(c(dist_meas, unit_lab))
#     
#     
#     # fam.anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas)) ),
#     #                                                        as.character(as.matrix(SLL2@sam_data[, unit_lab])),
#     #                                                        strata = SLL2.meta[ , "seqGroup"])
#     
#     mTab <- SLL2.meta[ SLL2.meta[,unit_lab] != "None", ]# **** must remove the None ****
#     fam.anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( full_ords[[ dist_meas ]][rownames(mTab),rownames(mTab)] ),
#                                                        as.character(as.matrix(mTab[, unit_lab])))
#   }
# }
# 
# saveRDS(fam.anosim_list, sprintf("%s/R_objects/fam.anosim_list.rds", p2_dir))


# # do the same for classmates 
# #   for those, use only teens, non-Chronic disorder, and those that came from high schools
# for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "Bray","Jaccard")) {
#   print(dist_meas)
#   
#   mTabClass <- meta.healthy[ ! is.na(meta.healthy$School_name) & 
#                                ! meta.healthy$School_name %in% c("Fabrica de Sol",
#                                                                  "Octubre centro de cultura contemporanea",
#                                                                  "Institut Frances", # only 1 sample would be left here, just remove
#                                                                  "Wetlab Laboratorios Cesar en Etopia") &
#                                ! is.na(meta.healthy$Age_groups) &
#                                meta.healthy$Age_groups == "Teen", ]
#   
#   fam.anosim_list[[dist_meas]][[ "School_name" ]] <- anosim(as.dist(full_ords[[ dist_meas ]][rownames(mTabClass),rownames(mTabClass)] ),
#                                                             as.character(as.matrix(mTabClass[, "School_name"])))
# }



# anosim(as.dist(weighted_Unifrac), SLL2@sam_data$Family_unit)
# anosim(as.dist(unweighted_Unifrac), SLL2@sam_data$Family_unit)
# 
# anosim(as.dist(weighted_Unifrac), SLL2@sam_data$Sibling_unit)
# anosim(as.dist(unweighted_Unifrac), SLL2@sam_data$Sibling_unit)

fam.anosim_list <- readRDS(sprintf("%s/R_objects/fam.anosim_list.rds", p2_dir))




# plot differences in average distances for group members vs non group members
# ****************************** #
plot_famUnit_diffs <- function(mTab, dist_obj, dist_name, comp.method, plotType="box", fam.anosim=NULL, fam.adonis=NULL) {
  
  fam_units <- c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
                 "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
                 "Grandparent_Grandchild_unit",
                 "School_name")
  
  all_unit_labs <- list()
  unit_labs_pvals <- list()
  unit_counts <- list()
  num_units <- list()
  
  for (unit_lab in fam_units) {
    
    # ***************** #
    if (unit_lab == "School_name") {
      sam_tab <- mTab[ ! is.na(mTab$School_name) & 
                         ! mTab$School_name %in% c("Fabrica de Sol",
                                                   "Octubre centro de cultura contemporanea",
                                                   "Institut Frances", # only 1 sample would be left here, just remove
                                                   "Wetlab Laboratorios Cesar en Etopia") &
                         ! is.na(mTab$Age_groups) &
                         mTab$Age_groups == "Teen" &
                         mTab$Chronic_disorder == "No", ]
    } else {
      sam_tab <- mTab
    }
    # ***************** #
    
    
    # ***************** #
    # get vector of units of particular type
    units <- unique( sam_tab[ , unit_lab] )
    units <- units[ units != "None" ]
    
    # get counts of individuals in each category
    unit_counts[[ gsub("_unit","",unit_lab) ]] <- nrow(sam_tab[ sam_tab[ , unit_lab ] != "None", ])
    num_units[[ gsub("_unit","",unit_lab) ]] <- length(unique(sam_tab[ sam_tab[ , unit_lab ] != "None", unit_lab]))
    
    # get vector of dists for each unit - members of unit against same unit
    all_unit_labs[[ sprintf("%s.within", unit_lab) ]] <- unlist(sapply(units, function(x) {
      unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
      # make dist_obj as.dist() first, 
      #   in this case it removes NAs and duplicated values
      dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
    }))
    
    # get vector of dists for each unit - members of unit against those not in that unit 
    #  (or against only other units??)
    all_unit_labs[[ sprintf("%s.between", unit_lab) ]] <- unlist(sapply(units, function(x) {
      unit_samps     <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
      non.unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
      
      # make dist_obj as.matrix() first,
      #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
      # dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
      dists <- as.numeric( as.matrix(dist_obj[ unit_samps, non.unit_samps ] ))
    }))
    # ***************** #
    
    
    # with.bet.diff <- t.test(within_dist, between_dist)
    # print(with.bet.diff)
    
    if (comp.method == "ANOSIM") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.anosim[[ dist_name ]][[ unit_lab ]]$signif
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- fam.anosim[[ dist_name ]][[ unit_lab ]]$statistic
      
    } else if (comp.method == "PERMANOVA") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$`Pr(>F)`[ 1 ]
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- sqrt(fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$R2[ 1 ])
      unit_labs_pvals[[ sprintf("%s.F", unit_lab) ]]   <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$F.Model[ 1 ]
    }
    
  }
  # print(str(all_unit_labs))
  # print(str(unit_labs_pvals))
  
  names(all_unit_labs) <- gsub("School_name", "Classmates", names(all_unit_labs))
  names(unit_labs_pvals) <- gsub("School_name", "Classmates", names(unit_labs_pvals))
  names(unit_counts) <- gsub("School_name", "Classmates", names(unit_counts))
  names(num_units) <- gsub("School_name", "Classmates", names(num_units))
  fam_units <- gsub("School_name", "Classmates", fam_units)
  
  # ************************ #
  if (plotType == "bar") {
    with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
                           "labs" = rep(c("Same Family","Different Families"), length(fam_units)),
                           "sds"  = unlist(lapply(all_unit_labs, sd)),
                           "unit" = rep(gsub("_unit","",fam_units), each = 2))
    
  } else if (plotType == "box") {
    wb.labs <- sapply(names(all_unit_labs), function(x) {
      if (endsWith(x, "within")) rep("Same Family", length(all_unit_labs[[ x ]]))
      else if (endsWith(x, "between")) rep("Different Families", length(all_unit_labs[[ x ]]))
    })
    
    wb.units <- sapply(gsub("_unit","",fam_units), function(x) {
      aul.names <- names(all_unit_labs)[ startsWith(names(all_unit_labs), x) ]
      rep(x, sum(sapply(aul.names, function(y) length(all_unit_labs[[ y ]]))))
    })
    
    with.bet <- data.frame("vals" = unlist(all_unit_labs),
                           "labs" = unlist(wb.labs),
                           "unit" = unlist(wb.units))
  }
  
  with.bet$labs <- factor(with.bet$labs, levels = rev(c("Same Family","Different Families")))
  with.bet$unit <- sapply(as.character(with.bet$unit), function(x) 
    sprintf("**%s**<br># of samples = %s<br># of units = %s", x, unit_counts[[ x ]], num_units[[ x ]]))
  # with.bet$unit <- factor(with.bet$unit, levels = rev(gsub("_unit","",fam_units)))
  with.bet$unit <- factor(with.bet$unit, levels = 
                            sapply(rev(gsub("_unit","",fam_units)), function(x) 
                              sprintf("**%s**<br># of samples = %s<br># of units = %s", x, unit_counts[[ x ]], num_units[[ x ]])))
  
  # ************************ #
  
  # print(unlist(lapply(all_unit_labs, mean)))
  # with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
  #                        "labs" = rep(c("Same Family","Different Families"), length(fam_units)),
  #                        "sds"  = unlist(lapply(all_unit_labs, sd)),
  #                        "unit" = rep(gsub("_unit","",fam_units), each = 2))
  # 
  # with.bet$unit <- factor(with.bet$unit, levels = rev(gsub("_unit","",fam_units)))
  
  # with.bet <- data.frame("vals" = c(mean(within_dist), mean(between_dist)),
  #                        "labs" = c("Within","Between"),
  #                        "sds"  = c(sd(within_dist), sd(between_dist)),
  #                        "unit" = c(unit_lab, unit_lab))
                         # "labs"=c( rep("Within", length(within_dist)), rep("Between",length(between_dist)) ))
  # print(with.bet)
  # print(unit_labs_pvals)
  
  star_labels <- sapply(fam_units, function(x) {
    if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
      sprintf("R=%s **", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
      sprintf("R=%s *", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else {
      # sprintf("R=%s", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      ""
    }
  })
  # print(star_labels)
  
  if (plotType == "bar") {
    
    ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_bar(stat="identity", position = "dodge2") +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
      coord_flip() +
      # geom_text(aes(x=unit, y=vals+0.05))
      annotate("text", x=1:length(fam_units), y=with.bet$vals[with.bet$labs=="Different Families"]+0.05, 
               label = rev(star_labels), size=5, color="red") +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15), legend.position = "bottom") +
      xlab(sprintf("%s in family units", comp.method)) + ylab(dist_name)
    
  } else if (plotType == "box") {
    
    ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_boxplot(notch = T) +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      # geom_text(aes(x=unit, y=vals+0.05))
      annotate("text", x=1:length(fam_units), 
               y=max(with.bet$vals[ startsWith(as.character(with.bet$labs), "Same") ]),
               label = rev(star_labels), size=5, color="red") +
      theme_minimal() +
      theme(axis.text.y = ggtext::element_markdown(size=15), 
            axis.text.x = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15), legend.position = "bottom") +
      xlab(sprintf("%s in family units", comp.method)) + ylab(dist_name)
    
  }
  
}
# ****************************** #

plot_famUnit_diffs(SLL2.meta, aitch, "Aitchison", "ANOSIM", fam.anosim = fam.anosim_list, plotType = "box")
plot_famUnit_diffs(SLL2.meta, weighted_Unifrac, "Weighted_Unifrac", "ANOSIM", fam.anosim = fam.anosim_list, plotType = "box")


# plot_famUnit_diffs(SLL2.meta, weighted_Unifrac, "weighted_Unifrac", "PERMANOVA")
# plot_famUnit_diffs(SLL2.meta, weighted_Unifrac, "weighted_Unifrac", "ANOSIM")
# 
# plot_famUnit_diffs(SLL2.meta, unweighted_Unifrac, "unweighted_Unifrac", "PERMANOVA")
# plot_famUnit_diffs(SLL2.meta, unweighted_Unifrac, "unweighted_Unifrac", "ANOSIM")
# 
# plot_famUnit_diffs(SLL2.meta, jsd, "jsd", "PERMANOVA")
# plot_famUnit_diffs(SLL2.meta, jsd, "jsd", "ANOSIM")
# 
# plot_famUnit_diffs(SLL2.meta, aitch, "aitch", "PERMANOVA")
# plot_famUnit_diffs(SLL2.meta, aitch, "Aitchison", "ANOSIM", fam.anosim_list)
# 
# plot_famUnit_diffs(SLL2.meta, weighted_Unifrac, "Weighted_Unifrac", "ANOSIM", fam.anosim_list)

# *********************************************************** #

# ******************************************* #
library(ggpubr)
plot_within_unit_dists <- function(mTab, dist_obj, dist_name, fam.anosim, comp.method="ANOSIM") {
  
  fam_units <- c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
                 "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
                 "Grandparent_Grandchild_unit",
                 "School_name")
  
  all_unit_labs <- list()
  unit_labs_pvals <- list()
  unit_counts <- list()
  num_units <- list()
  
  for (unit_lab in fam_units) {
    
    # ***************** #
    if (unit_lab == "School_name") {
      sam_tab <- mTab[ ! is.na(mTab$School_name) & 
                         ! mTab$School_name %in% c("Fabrica de Sol",
                                                   "Octubre centro de cultura contemporanea",
                                                   "Institut Frances", # only 1 sample would be left here, just remove
                                                   "Wetlab Laboratorios Cesar en Etopia") &
                         ! is.na(mTab$Age_groups) &
                         mTab$Age_groups == "Teen" &
                         mTab$Chronic_disorder == "No", ]
    } else {
      sam_tab <- mTab
    }
    # ***************** #
    
    
    # ***************** #
    # get vector of units of particular type
    units <- unique( sam_tab[ , unit_lab] )
    units <- units[ units != "None" ]
    
    # get counts of individuals in each category
    unit_counts[[ gsub("_unit","",unit_lab) ]] <- nrow(sam_tab[ sam_tab[ , unit_lab ] != "None", ])
    num_units[[ gsub("_unit","",unit_lab) ]] <- length(unique(sam_tab[ sam_tab[ , unit_lab ] != "None", unit_lab]))
    
    # get vector of dists for each unit - members of unit against same unit
    all_unit_labs[[ unit_lab ]] <- unlist(sapply(units, function(x) {
      unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
      # make dist_obj as.dist() first, 
      #   in this case it removes NAs and duplicated values
      dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
    }))
    # ***************** #
    
    
    # with.bet.diff <- t.test(within_dist, between_dist)
    # print(with.bet.diff)
    
    if (comp.method == "ANOSIM") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.anosim[[ dist_name ]][[ unit_lab ]]$signif
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- fam.anosim[[ dist_name ]][[ unit_lab ]]$statistic
      
    } else if (comp.method == "PERMANOVA") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$`Pr(>F)`[ 1 ]
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- sqrt(fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$R2[ 1 ])
      unit_labs_pvals[[ sprintf("%s.F", unit_lab) ]]   <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$F.Model[ 1 ]
    }
    
  }
  
  
  names(all_unit_labs) <- gsub("School_name", "Classmates", names(all_unit_labs))
  names(unit_labs_pvals) <- gsub("School_name", "Classmates", names(unit_labs_pvals))
  names(unit_counts) <- gsub("School_name", "Classmates", names(unit_counts))
  names(num_units) <- gsub("School_name", "Classmates", names(num_units))
  fam_units <- gsub("School_name", "Classmates", fam_units)
  
  aul.melt <- reshape::melt.list(all_unit_labs)
  colnames(aul.melt) <- c("dist","units")
  
  res.aov <- aov(formula = as.numeric(as.matrix(aul.melt$dist)) ~ as.factor(as.matrix(aul.melt$units)), data = aul.melt)
  TukeyTab <- TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]
  
  group1s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][1])
  group2s <- sapply(rownames(TukeyTab), function(x) strsplit(x, "-")[[1]][2])
  cont.max <- max(aul.melt$dist, na.rm = T)
  # arbitrarily make the space between lines based on the range of values
  rango <- (max(aul.melt$dist, na.rm = T)-min(aul.melt$dist, na.rm = T)) / 10
  # cont.max will determine the y.position for signif lines
  cont.max <- cont.max + rango
  cont.max <- seq( cont.max, cont.max+((nrow(TukeyTab)-1)*rango), by=rango)
  stat.test <- data.frame(".y."=rep("cont", nrow(TukeyTab)), 
                          "group1"=group1s, "group2"=group2s,
                          "p"=rep(NA, nrow(TukeyTab)), 
                          "p.adj"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                          "p.format"=formatC(TukeyTab[,"p adj"], format="e", digits=3),
                          "p.signif"=symnum(TukeyTab[,"p adj"], 
                                            cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                            symbols = c("****", "***", "**", "*", "ns")),
                          "method"=rep("Tukey", nrow(TukeyTab)),
                          "y.position"=cont.max)
  
  # print(stat.test)
  
  # ******************* #
  aul.melt.sibs <- aul.melt[ aul.melt$units %in% c("Twin_unit","Sibling_unit"), ]
  
  # print(kruskal.test(as.numeric(aul.melt.sibs$dist), as.factor(aul.melt.sibs$units)))
  # print(wilcox.test(aul.melt.sibs$dist[ aul.melt.sibs$units == "Twin_unit"], aul.melt.sibs$dist[ aul.melt.sibs$units == "Sibling_unit"]))
  print(wilcox.test(dist ~ units, data=aul.melt.sibs))
  # print(t.test(aul.melt.sibs$dist[ aul.melt.sibs$units == "Twin_unit"], aul.melt.sibs$dist[ aul.melt.sibs$units == "Sibling_unit"]))
  # hist(aul.melt.sibs$dist[ aul.melt.sibs$units == "Sibling_unit"], breaks=6)
  
  # ******************* #
  ggplot(aul.melt, aes(x=units, y=dist)) +
    geom_boxplot(notch = T) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    stat_pvalue_manual(stat.test, label = "p = {p.adj}", inherit.aes=F)
}
# ******************************************* #

# prints Mann-Whitney test result for Twin vs Sibling units
plot_within_unit_dists(SLL2.meta, aitch, "Aitchison", fam.anosim_list)



hist(SLL2.meta[SLL2.meta$Twin_unit!="None", "Age"])
sort(SLL2.meta[SLL2.meta$Twin_unit!="None", "Age"])

hist(SLL2.meta[SLL2.meta$Sibling_unit!="None", "Age"])
sort(SLL2.meta[SLL2.meta$Sibling_unit!="None", "Age"])

mTabClass <- meta.healthy[ ! is.na(meta.healthy$School_name) &
                             ! meta.healthy$School_name %in% c("Fabrica de Sol",
                                                               "Octubre centro de cultura contemporanea",
                                                               "Institut Frances", # only 1 sample would be left here, just remove
                                                               "Wetlab Laboratorios Cesar en Etopia") &
                             ! is.na(meta.healthy$Age_groups) &
                             meta.healthy$Age_groups == "Teen", ]

# *********************************************************** #

get_drivers("Genus", "Twin_unit", gloms_clr, SLL2.meta[SLL2.meta$Twin_unit != "None",], nTop = 10)
topgen.Twin <- unname(unlist(lapply(get_drivers("Genus", "Twin_unit", gloms_clr, SLL2.meta[SLL2.meta$Twin_unit != "None",], nTop = 10), names)))
sort(topgen.Twin[ duplicated(topgen.Twin) ])

get_drivers("Genus", "Twin_unit", gloms_clr, SLL2.meta[SLL2.meta$Twin_unit != "None",], nTop = 10, strongest = F)
bottomgen.Twin <- unname(unlist(lapply(get_drivers("Genus", "Twin_unit", gloms_clr, SLL2.meta[SLL2.meta$Twin_unit != "None",], nTop = 10, strongest = F), names)))
sort(bottomgen.Twin[ duplicated(bottomgen.Twin) ])
# *********************************************************** #



# plot Age ranges for particular family unit types
hist(SLL2.meta$Age, main="All samples", xlab = "Age")

hist(SLL2.meta[ ! SLL2.meta$Family_participants %in% c("No", "No Sabe/No Contest", "None", NA), "Age"], 
     main="Family", xlab = "Age")

hist(SLL2.meta[ ! SLL2.meta$Sibling %in% c("No", "No Sabe/No Contest", "None", NA), "Age"], 
     main="Sibling", xlab = "Age")

hist(SLL2.meta[ SLL2.meta$Sibling %in% c("Gemelo","Melliza"), "Age"], 
     main="Twin", xlab = "Age")

hist(SLL2.meta[ ! SLL2.meta$Partner %in% c("No", "No Sabe/No Contest", "None", NA), "Age"], 
     main="Partner", xlab = "Age")

hist(SLL2.meta[ ! SLL2.meta$Parent_Child_unit %in% c("No", "No Sabe/No Contest", "None", NA), "Age"], 
     main="Parent_Child_unit", xlab = "Age")

hist(SLL2.meta[ ! SLL2.meta$Grandparent_Grandchild_unit %in% c("No", "No Sabe/No Contest", "None", NA), "Age"], 
     main="Grandparent_Grandchild_unit", xlab = "Age")

# *********************************************************** #
























































# samps <- "All"
# 
# if (samps == "All") {
#   samps <- sample_names(SLL2_rel)
#   ttc <- SLL2_rel
#   clus_var_name <- "Stomatotype"
# }
# 
# tax.to.clust <- ttc@otu_table
# 
# 
# 
# noise.removal <- function(dataframe, percent=0.01, top=NULL){
#   dataframe->Matrix
#   bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
#   Matrix_1 <- Matrix[bigones,]
#   # print(percent)
#   return(Matrix_1)
# }
# 
# tax.to.clust <- noise.removal(tax.to.clust, percent=0.01)
# 
# 
# 
# # Distance matrix object
# dist_meas <- "JSD"
# # dist_meas <- "Weighted_Unifrac"
# # dist_meas <- "Unweighted_Unifrac"
# 
# if (dist_meas == "JSD") {
#   jsd = dist.JSD(tax.to.clust)
# } else if (dist_meas == "Weighted_Unifrac") {
#   jsd = as.dist(weighted_Unifrac)
# } else if (dist_meas == "Unweighted_Unifrac") {
#   jsd = as.dist(unweighted_Unifrac)
# }
# 
# # Cluster using the Partitioning Around Medoids (PAM) algorithm
# pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
#   clust = as.vector(pam(as.dist(x), k, diss=TRUE))
#   return(clust$clustering)
# }
# 
# 
# 
# CH.index <- obs.silhouette <- NULL
# 
# for (k in 1:10) { 
#   if (k == 1) {
#     CH.index[k] = NA
#     obs.silhouette[k] <- NA
#   } else {
#     cluster_temp = pam.clustering(jsd, k)
#     CH.index[k] = index.G1(t(tax.to.clust), cluster_temp,  d = jsd, centrotypes = "medoids")
#     obs.silhouette[k] <- mean(silhouette(cluster_temp, jsd)[,3])
#   }
# }
# 
# 
# clus.meas <- cbind(CH.index, obs.silhouette)
# clus.meas <- apply(clus.meas, 2, function(x) x/max(x, na.rm = T))
# 
# ggplot(melt(clus.meas), aes(x=Var1, y=value, fill=Var2)) +
#   geom_bar(stat="identity", position = position_dodge(), width = 0.5) +
#   scale_x_continuous(breaks = seq(0,10,1)) + ggtitle('Optimal number of clusters by CH index and obs silhouette') +
#   xlab("k clusters") + ylab('CH index / obs silhouette') + scale_fill_hue(name="measure")
# 
# nc <- sapply(1:nrow(clus.meas), function(x) sum(clus.meas[x,"CH.index"], clus.meas[x,"obs.silhouette"]))
# nclusters <- match(max(nc, na.rm = T), nc)
# 
# 
# 
# cluster_vals = pam.clustering( jsd, k = nclusters )
# 
# # must add values for those samples not included when appropriate
# cluster_vals_full <- cluster_vals
# cluster_vals_full[ sample_names(SLL2_rel)[ ! sample_names(SLL2_rel) %in% samps ] ] <- NA
# 
# 
# if (nclusters == 2) { 
#   clust.col <- c("indianred1", "turquoise3")
# } else if (nclusters == 3) { 
#   clust.col <- c("indianred1", "green4", "dodgerblue")
# } else if (nclusters == 4) { 
#   clust.col <- c("indianred1", "yellowgreen", "turquoise3", "purple1")
# } else if (nclusters == 5) {
#   clust.col <- c("indianred1", "#829311", "seagreen3", "dodgerblue", "orchid")
# }
# 
# 
# 
# 
# # ************************* #
# ## plot 1: between-class analysis (BCA) ####
# 
# pca = dudi.pca(t(tax.to.clust), scannf=F, nf=k)
# bet = bca(pca, fac=as.factor(cluster_vals), scannf=F, nf=k-1)
# 
# if (nclusters > 2) {
#   #with only 2 clusters, apparently unable to produce figure bc does not create enough columns in obs.bet$ls
#   s.class(bet$ls, fac=as.factor(cluster_vals), grid=F,sub="Between-class analysis", col=clust.col)
#   # color.by <- "Enterotype"
#   # s.class(bet.stool$ls, fac=as.factor(as.matrix(OCD.stool@sam_data)[ , color.by]),
#   #         grid=F,sub="Between-class analysis", col=hue_pal()(nrow(unique(OCD.stool@sam_data[ , color.by]))))
# }
# 
# 
# 
# # ************************* #
# stomato.leaders <- NULL
# for (i in 1:nclusters){
#   stomato.leaders[i] <- colnames(bet$tab)[bet$tab[i,]==max(bet$tab[i,])]
# }
# 
# # ************************* #
# # check boxplots of rel abund of each stomatotype
# samples.by.stomato <- as.data.frame(cluster_vals)
# samples.by.stomato[,2] <- names(cluster_vals) #adding a column for corresponding enterotype with a given sample
# 
# #add stomatotype to sample_data
# # must specify row order in the case that not all samples were included for the clustering, since the order in cluster_vals_full will be different
# SLL2@sam_data[ names(cluster_vals_full), clus_var_name ] <- as.factor(cluster_vals_full)
# SLL2_rel@sam_data[ names(cluster_vals_full), clus_var_name ] <- as.factor(cluster_vals_full)
# 
# 
# # ************************* #
# 
# stomato.melts <- list()
# for (s in stomato.leaders) {
#   stomato.melts[[ s ]] <- cbind( melt(t(tax.to.clust[ s , ])), melt(samples.by.stomato) )
#   stomato.melts[[ s ]] <- stomato.melts[[ s ]][,c(1,2,3,6)]
#   colnames(stomato.melts[[ s ]])[4] <- 'cluster_vals'
# }
# 
# 
# # to give the same scale in each boxplot
# ymax <- max(tax.to.clust[stomato.leaders,])
# 
# stomato.box.leaders <- list()
# for (s in stomato.leaders) {
#   stomato.box.leaders[[ s ]] <- ggplot(stomato.melts[[ s ]], 
#                                       aes(x=factor(cluster_vals, levels=c(1:nclusters)), y=value, 
#                                           fill=factor(cluster_vals, levels=c(1:nclusters)))) +
#     geom_boxplot(notch = T) +
#     ggtitle( s ) + ylim(0, ymax) +
#     xlab('Stomatotype') + ylab('% per sample') + theme(legend.position="none")
# }
# 
# library(Rmisc)
# multiplot(plotlist = stomato.box.leaders, cols = length( stomato.box.leaders ))
# 
# 
# 
# # ************************* #
# # plot 2: principal coordinates analysis (PCoA) ####
# library(scales)
# pcoa = dudi.pco(jsd, scannf=F, nf=3)
# ade4::s.class(pcoa$li, fac=as.factor(as.matrix(SLL2@sam_data)[ samps, clus_var_name]), grid=F, 
#               col=hue_pal()(nrow(unique(SLL2@sam_data[ samps, clus_var_name]))), clabel=1.5,
#               sub = sprintf("PCoA separated by %s", clus_var_name))



# ************************* #
# plot 2.2: Correspondence Analysis and checking contributions of both samples and taxa ####

dist_meas <- "Weighted_Unifrac"

color.by <- sprintf("Stomatotype_%s", dist_meas)
coa.noOutliers <- dudi.coa(t(tax.to.clust[,!colnames(tax.to.clust) %in% c("SLL.59.05")]), scannf = F, nf = 10)
# coa <- dudi.coa(t(tax.to.clust), scannf = F, nf = 10)
facs <- as.factor(as.matrix(SLL2@sam_data)[ !sample_names(SLL2) %in% c("SLL.59.05"), color.by])
ade4::s.class(coa.noOutliers$li, fac = facs, grid=F, 
              col=hue_pal()(length(levels(facs))), clabel=1.5,
              sub = sprintf("PCoA separated by %s", color.by))
library(factoextra)
fviz_contrib(coa.noOutliers, top = 15, choice = "col")
fviz_contrib(coa.noOutliers, top = 15, choice = "row")
fviz_ca(coa.noOutliers)

# ************************* #





# ************************* #
# plot 3: 3d PCA plots ####
library(pca3d)

prcomp <- prcomp(jsd)
# pca3d(prcomp.stool, group=as.factor(cluster.stool), show.centroids = T, show.group.labels = T, palette = hue_pal()(nclusters.stool), 
#       show.ellipses = T, ellipse.ci=0.75, show.plane = F)
pca3d(prcomp, group=as.factor(as.matrix(SLL2@sam_data)[ , color.by]), show.centroids = T, 
      show.group.labels = T, palette = hue_pal()(nrow(unique(SLL2@sam_data[ , color.by]))), 
      show.ellipses = T, ellipse.ci=0.75, show.plane = F,
      title = sprintf("3D PCA separated by %s", clus_var_name))

# ************************* #






# ****************************************************************************************************************** #
# ************************************************ #
# Boxplots of top15 within each cluster ####
# ************************************************ #

samps <- sample_names(SLL2)

dis <- "Weighted_Unifrac"

tlev <- "Species"; otutab_rel <- SLL2_rel@otu_table[,samps]
# tlev <- "Genus"; otutab_rel <- gloms_rel[[ tlev ]][,samps]
# tlev <- "Phylum"; otutab_rel <- gloms_rel[[ tlev ]][,samps]

top15 <- names(sort(rowSums(otutab_rel[,samps]), decreasing = T))[1:15]
SLL2.top <- cbind( reshape2::melt(t(otutab_rel[top15[!is.na(top15)], samps])), 
                   reshape2::melt(as.matrix(tagl@sam_data[,sprintf("Stomatotype_%s",dis)])))


SLL2.top <- SLL2.top[, c(1,2,3,6)]
colnames(SLL2.top) <- c("sample","OTU","value","cluster_vals")
# to order the boxes by overall abundance
SLL2.top$OTU <- as.character(SLL2.top$OTU)
SLL2.top$OTU <- factor(SLL2.top$OTU, levels = top15)

samples.by.stomato <- as.character(as.matrix(tagl@sam_data[,sprintf("Stomatotype_%s",dis)]))

cluster_labels <- c('1'=sprintf('Stomatotype 1 (n=%s)',sum(samples.by.stomato==1)),
                    '2'=sprintf('Stomatotype 2 (n=%s)',sum(samples.by.stomato==2)),
                    '3'=sprintf('Stomatotype 3 (n=%s)',sum(samples.by.stomato==3)),
                    '4'=sprintf('Stomatotype 4 (n=%s)',sum(samples.by.stomato==4)),
                    '5'=sprintf('Stomatotype 5 (n=%s)',sum(samples.by.stomato==5)),
                    '6'=sprintf('Stomatotype 6 (n=%s)',sum(samples.by.stomato==6)),
                    '7'=sprintf('Stomatotype 7 (n=%s)',sum(samples.by.stomato==7)),
                    '8'=sprintf('Stomatotype 8 (n=%s)',sum(samples.by.stomato==8)),
                    '9'=sprintf('Stomatotype 9 (n=%s)',sum(samples.by.stomato==9)),
                    '10'=sprintf('Stomatotype 10 (n=%s)',sum(samples.by.stomato==10)))

# ggplot(OCD.top.stool, aes(x=reorder(factor(Var2, levels=top15.stool), -value, FUN=median), y=value, 
# fill=reorder(factor(Var2, levels=top15.stool), -value, FUN=median))) +
ggplot(SLL2.top, aes(x=OTU, y=value, # This way keeps the order of top15 instead of ordering by median
                          fill=OTU)) +
  geom_boxplot(notch = T) +
  facet_wrap(~cluster_vals, ncol=1, labeller=as_labeller(cluster_labels)) +
  # theme(axis.text.x=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust=0.5)) +
  ggtitle(sprintf('Percent of given %s per sample by indicated %s', tlev, sprintf("Stomatotype_%s",dis))) +
  xlab(tlev) + ylab('% per sample') + scale_fill_hue(name=tlev) + guides(fill=F)





# ********** #
# to do paired boxes for top5

tl <- "Species"; tagl <- SLL2_rel
tl <- "Genus"; tagl <- phyloseq(otu_table(gloms_rel[[ tl ]]), sample_data(SLL2_rel))

samps <- sample_names(tagl)

dis <- "Weighted_Unifrac"

library(ggpubr)
top5 <- names(sort(rowSums(tagl@otu_table[,samps]), decreasing = T))[1:5]
SLL2.top5 <- cbind( reshape2::melt(t(as.data.frame(tagl@otu_table)[top5,])), 
                    reshape2::melt(as.matrix(tagl@sam_data[,sprintf("Stomatotype_%s",dis)])))


SLL2.top5 <- SLL2.top5[, c(1,2,3,6)]
colnames(SLL2.top5) <- c("sample","Species","value","Stomatotype")
# to order the boxes by overall abundance
SLL2.top5$Species <- as.character(SLL2.top5$Species)
SLL2.top5$Species <- factor(SLL2.top5$Species, levels = top5)
SLL2.top5$Stomatotype <- factor(SLL2.top5$Stomatotype)

ggplot(SLL2.top5, aes(x=Species, y=value, fill=Stomatotype)) +
  geom_boxplot() + theme_minimal() +
  theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle = 30),
        axis.title.y=element_blank(), axis.title.x=element_text(size=14, face="bold")) +
  guides(fill = guide_legend(title.theme=element_text(size=13, angle=0, face="bold"), 
                             label.theme=element_text(size=16, angle=0)) ) + 
  stat_compare_means(method = "wilcox.test", label = "p.format", cex=5, color="red3", 
                     label.y=c(61, 37, 30, 25, 20))
# theme(axis.text = element_text(size = 14), plot.title = element_text(hjust=0.5)) +
# ggtitle(sprintf('Percent of given Species per sample by indicated Stomatotype'))
# ********** #










# ************************************************ #
# get list of taxa for each stomatotype that are significantly more abundant than in other types
# ************************************************ #

stomato.signif.greater <- list()
for (s in 1:nclusters) {
  stomato.signif.greater[[ sprintf("Stomatotype_%s",s) ]] <- c()
  
  for (tax in names(sort(taxa_sums(SLL2_rel), decreasing = T))) {
    
    tg <- t.test(as.numeric(as.matrix(SLL2_rel@otu_table[ tax, rownames(samples.by.stomato[samples.by.stomato$cluster_vals == s,])])), 
                 as.numeric(as.matrix(SLL2_rel@otu_table[ tax, rownames(samples.by.stomato[samples.by.stomato$cluster_vals %in% 
                                                                                                 (1:nclusters)[1:nclusters != s],])])),
                 alternative = "greater")$p.value
    
    if (tg < 0.05) {
      stomato.signif.greater[[ sprintf("Stomatotype_%s",s) ]] <- c( stomato.signif.greater[[ sprintf("Stomatotype_%s",s) ]], tax )
    }
  }
}


# ****************************************************************************************************************** #















# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####



# ****************************************************************************************************************** #
# PERMANOVA with adonis ####
# ****************************************************************************************************************** #


# ***************************************************************************************** #
# ***************************************************************************************** #

adonis_vars <- c("Gender","Age","Age_groups","Weight","Height","BMI","Municipal_zone","Number_inhabitants",
                 "Years_in_home","Moisture_in_home",
                 "Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds","Pets.Reptiles_amphibians",
                 "Pets.Fish","Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
                 # "Pets.Type.Small_furry_animals",
                 "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish","Pets.Rabbits",
                 "Pets.Number.Rabbits","Pets.Rodents","Pets.Number.Rodents","Pets.Mammals","Pets.Number.Mammals",
                 "Pets.Horses","Pets.Number.Horses",
                 "Smoker","Number_smokers_home",
                 "Consumption.Milk","Consumption.Yogurt","Consumption.Sweets","Consumption.Chewing_gum",
                 "Consumption.Nuts","Water_type_home",
                 "Drinks.Decaf_coffee","Drinks.Coffee","Drinks.Tea","Drinks.Infusion","Drinks.Soda",
                 "Drinks.Soda_sugarless","Drinks.Soda_decaf","Drinks.RedBull","Drinks.Other_sugary_drinks",
                 "Drinks.Alcohol_cold","Drinks.Alcohol_hot",
                 "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth",
                 "Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement","Brushing","Floss",
                 "Mouth_wounds","Last_dental_visit",
                 "Chronic_disorder","Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome",
                 "Eating_disorder","Medications","Antibiotics","Analgesics","Vitamin_supplements",
                 "Asthma","Wheezing","How_do_you_feel","Do_you_feel_well",
                 "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen","Allergy.Animals",
                 "Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
                 "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary",
                 "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
                 "Kissing_partner","pH",
                 "MALDI.Yeast_detected","MALDI.Mold_detected",
                 cont_water_data, group_water_data, "Population","Community","Province",additional_diseases,"seqGroup")

# adonis_vars.short <- c("Gender","Age","Age_groups","BMI","Municipal_zone","Moisture_in_home",
#                        "Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds","Pets.Reptiles_amphibians",
#                        "Smoker",
#                        "Consumption.Milk","Consumption.Yogurt","Consumption.Sweets","Consumption.Chewing_gum",
#                        "Consumption.Nuts","Water_type_home",
#                        "Drinks.Decaf_coffee","Drinks.Coffee","Drinks.Tea","Drinks.Infusion","Drinks.Soda",
#                        "Drinks.Soda_sugarless","Drinks.Soda_decaf","Drinks.RedBull","Drinks.Other_sugary_drinks",
#                        "Drinks.Alcohol_cold","Drinks.Alcohol_hot",
#                        "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth",
#                        "Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement","Brushing","Floss",
#                        "Mouth_wounds","Last_dental_visit",
#                        "Medications","Antibiotics","Analgesics","Vitamin_supplements",
#                        "Asthma","Wheezing","How_do_you_feel","Do_you_feel_well",
#                        "Allergy",
#                        "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
#                        "Kissing_partner","pH","Population")

covs <- c("Chronic_disorder","Gender","Age","Population")#



# # make a table without "No Sabe/No Contesta"
# metaTab <- SLL2.meta
# metaTab[ metaTab=="No Sabe/No Contesta"] <- NA
# # Add a column that has all 1 value to force all to be same shape by default
# metaTab$Project <- rep("SLL2", nrow(metaTab))
# metaTab$Weight                     <- as.numeric(metaTab$Weight)
# metaTab$Height                     <- as.numeric(metaTab$Height)
# metaTab$Dental.Fillings            <- as.numeric(metaTab$Dental.Fillings)
# metaTab$Dental.Nerve_extractions   <- as.numeric(metaTab$Dental.Nerve_extractions)
# metaTab$Dental.Lost_teeth          <- as.numeric(metaTab$Dental.Lost_teeth)
# metaTab$Dental.Reconstructed_teeth <- as.numeric(metaTab$Dental.Reconstructed_teeth)
# metaTab$pH[ metaTab$pH %in% c(7.55,7.6) ] <- "7.5"
# metaTab$pH <- as.numeric(metaTab$pH)  




# # *********************************************************** #
# adonis_list <- list()
# for (trait in adonis_vars) {
# # for (trait in adonis_vars.short) {
# # for (trait in c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats",
# #                 "Pets.Number.Small_furry_animals","Pets.Number.Birds","Pets.Number.Reptiles_amphibians",
# #                 "Pets.Number.Fish","Number_smokers_home","Dental.Fillings","Dental.Nerve_extractions",
# #                 "Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")) {
# 
#   print(trait)
#   adonis_list[[ trait ]] <- list()
# 
#   if (trait %in% c("Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome",
#                    "Eating_disorder",additional_diseases)) {
#     # ************* #
#     if (trait %in% c("Cystic_fibrosis","Transplant","Immune_issues")) {
#       cov_vars <- c("Antibiotics","Gender","Age","Population")#
#     } else {
#       cov_vars <- c("Gender","Age","Population")#
#     }
#     # ************* #
#   } else if (trait %in% c(cont_water_data, group_water_data, "Community","Province")) {
#     cov_vars <- c("Chronic_disorder","Gender","Age")#
#   } else {
#     cov_vars <- covs
#   }
# 
#   for (dist_meas in c("Weighted_Unifrac","Unweighted_Unifrac",
#                       # "JSD","VAW_GUnifrac","a0_GUnifrac","Canberra","a05_GUnifrac",
#                       "Bray","Jaccard", "Aitchison"
#   )) {
# 
#     # print(dist_meas)
# 
#     adonis_list[[ trait ]][[ dist_meas ]] <- get_adonis(trait, cov_vars, dist_meas, SLL2.meta)
#   }
# }
# # *********************************************************** #
# saveRDS(adonis_list, sprintf("%s/figures/Adonis_signifs/adonis_list.rds", p2_dir))
# # ******************************************** #

adonis_list <- readRDS(sprintf("%s/figures/Adonis_signifs/adonis_list.rds", p2_dir))
# *********************************************************** #

c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
  "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish","Number_smokers_home",
  "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")








# *********************************************************** #

# get the pvals for each trait with each dist_meas
adonis.pvals <- sapply(names(adonis_list), function(x) 
  lapply(adonis_list[[ x ]], function(y)
    y[x, "Pr(>F)"]))
# write only those pvals < 0.1
adonis.signifs <- adonis.pvals
# adonis.signifs[ adonis.signifs > 0.1 ] <- ""
adonis.signifs[ adonis.signifs > 0.05 ] <- ""

t(adonis.signifs[ , colSums(adonis.signifs!="") > 0 ])




# *********************************************************** #

# for each dist_meas, get pvals of covariates, to see what covaries with each trait
adonis.covs.pvals <- lapply(adonis_list, function(y) {
  cov.ps <- sapply(names(y), function(x) {
    # include any interaction terms, if they appear
    trait.interactions <- rownames(y[[x]])[ grepl(":", rownames(y[[x]])) ]
    cps <- y[[ x ]][ c(covs, trait.interactions), "Pr(>F)"]
    names(cps) <- c(covs, trait.interactions)
    return(cps)
  })
  
  if (class(cov.ps) == "matrix" ) {
    # in this case, all interaction terms were included for each dist_meas, 
    #   so no need to deal with added values
    return(cov.ps)
    
  } else {
    # give value of 1 for interaction terms that do not appear for a given dist_meas, 
    #   will simply be filtered out in next step, but this allows to make a data frame
    allcovs <- unique(unlist(lapply(cov.ps, names)))
    cov.ps.all <- lapply(cov.ps, function(x) {
      newX <- x
      # add 1 for missing interaction term, along with the name
      newX[ allcovs[ ! allcovs %in% names(x) ] ] <- 1
      # then return this vector with all terms in the same order each time, to make sure rows are correct in next step
      return( newX[ allcovs ] )
    })
    
    return( as.matrix( as.data.frame(cov.ps.all) ) )
  }
  
})
# *********************************** #

# then write just those which are significant, ignoring the values for same trait (if one of the covs)
adonis.covs.signifs <- lapply(names(adonis.covs.pvals), function(x) {
  covsTab <- adonis.covs.pvals[[ x ]]
  # covsTab[ covsTab > 0.1 ] <- ""
  covsTab[ covsTab > 0.05 ] <- ""
  
  if (x %in% covs) {
    covsTab[x, covsTab[x, ] != "" ] <- ""
  }
  return(covsTab)
}); names(adonis.covs.signifs) <- names(adonis.covs.pvals) # lapply does not give the names to the list items


# ***************************************************************************************** #















# ****************************************************************************************************************** #
# plot ordination with significant adonis variables ####
# *********************************************************** #






# # make a table without "No Sabe/No Contesta"
# metaTab <- SLL2.meta
# metaTab[ metaTab=="No Sabe/No Contesta"] <- NA
# # Add a column that has all 1 value to force all to be same shape by default
# metaTab$Project <- rep("SLL2", nrow(metaTab))
# metaTab$Weight                     <- as.numeric(metaTab$Weight)
# metaTab$Height                     <- as.numeric(metaTab$Height)
# metaTab$Dental.Fillings            <- as.numeric(metaTab$Dental.Fillings)
# metaTab$Dental.Nerve_extractions   <- as.numeric(metaTab$Dental.Nerve_extractions)
# metaTab$Dental.Lost_teeth          <- as.numeric(metaTab$Dental.Lost_teeth)
# metaTab$Dental.Reconstructed_teeth <- as.numeric(metaTab$Dental.Reconstructed_teeth)
# metaTab$pH[ metaTab$pH %in% c(7.55,7.6) ] <- "7.5"
# metaTab$pH <- as.numeric(metaTab$pH)  

# *********************************************************** #

for (trait in colnames(adonis.signifs)) {
  
  print(trait)
  
  for (dist_meas in rownames(adonis.signifs)) {
    
    # plot only those which have signif pvals
    if ( adonis.signifs[ dist_meas, trait ] != "") {
      
      # Case.control is default shape
      # if Case.control is not signif as covariate, check the others
      # if none other is signif, Case.control will remain the shape variable
      shapeBy <- "Project"
      # if (trait != "Chronic_disorder" & # dont want to change shape for Case.control
      #     #     adonis.covs.signifs[[ trait ]]["Case.control",dist_meas] %in% c("",NA) & # if Case.control not a covar
      #     #     ! adonis.covs.signifs[[ trait ]]["Outcome",dist_meas] %in% c("",NA)) # Outcome is a covar
      #     #   shapeBy <- "Outcome"
      #     # else if (trait != "Case.control" & # dont want to change shape for Case.control
      #     adonis.covs.signifs[[ trait ]]["Chronic_disorder",dist_meas] %in% c("",NA) & # if Case.control not a covar
      #     ! adonis.covs.signifs[[ trait ]]["Gender",dist_meas] %in% c("",NA)) # Gender is a covar
      #   shapeBy <- "Gender"
      
      # *********************************************************** #
      # make plot
      # if (trait == "Month.BAL")
      #   plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
      #                        adonis.pvals, adonis.covs.signifs, save.adoPlot = T, pcoaList = pcoas, ellipses = F)
      # else if (trait %in% c("Number_inhabitants","Dental.Fillings","Dental.Nerve_extractions",
      #                       "Dental.Lost_teeth","Dental.Reconstructed_teeth"))
      #   plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
      #                        adonis.pvals, adonis.covs.signifs, save.adoPlot = T, pcoaList = pcoas, log_of_color = T)
      # else
      #   plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
      #                        adonis.pvals, adonis.covs.signifs, save.adoPlot = T, pcoaList = pcoas)
      
      # *********************************************************** #
      # then also make plots using the ncMCE coordinates
      if (trait == "Month.BAL") {
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T, ellipses = F,
                             mce_tab = ncMCEs[[ dist_meas ]], mce_centered = "No")
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T, ellipses = F,
                             mce_tab = MCEs[[ dist_meas ]], mce_centered = "Yes")
        
      } else if (trait %in% c("Number_inhabitants","Dental.Fillings","Dental.Nerve_extractions",
                            "Dental.Lost_teeth","Dental.Reconstructed_teeth")) {
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T, log_of_color = T,
                             mce_tab = ncMCEs[[ dist_meas ]], mce_centered = "No")
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T, log_of_color = T,
                             mce_tab = MCEs[[ dist_meas ]], mce_centered = "Yes")
      } else {
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T,
                             mce_tab = ncMCEs[[ dist_meas ]], mce_centered = "No")
        plot_adonis_w_covars(dist_meas, SLL2, SLL2.meta, gloms_rel, shapeBy, trait, NA, 
                             adonis.pvals, adonis.covs.signifs, save.adoPlot = T,
                             mce_tab = MCEs[[ dist_meas ]], mce_centered = "Yes")
      }
      
      # *********************************************************** #
    }
  }
}

c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
  "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish","Number_smokers_home",
  "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")




# *********************************************************** #

# plot_adonis_w_covars("Jaccard", SLL2, metaTab, gloms_rel, 
#                      shapeBy = "Project", colorBy = "Age", 
#                      NA, adonis.pvals, adonis.covs.signifs,
#                      pcoaList = NULL)

plot_adonis_w_covars("Jaccard", SLL2, SLL2.meta, gloms_rel, 
                     shapeBy = "Project", colorBy = "Age", 
                     NA, adonis.pvals, adonis.covs.signifs,
                     pcoaList = pcoas)


# ****************************************************************************************************************** #































# ****************************************************************************************************************** #
# Run linear model ####



# ************************** #
# check for collinearity between continuous variables
fixed.continuous <- c("Age","BMI",
                      "Consumption.Milk","Consumption.Yogurt","Consumption.Sweets","Consumption.Chewing_gum",
                      "Consumption.Nuts",
                      "Drinks.Decaf_coffee","Drinks.Coffee","Drinks.Tea","Drinks.Infusion","Drinks.Soda",
                      "Drinks.Soda_sugarless","Drinks.Soda_decaf","Drinks.RedBull","Drinks.Other_sugary_drinks",
                      "Drinks.Alcohol_cold","Drinks.Alcohol_hot",
                      "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth",
                      "Brushing","Floss","Last_dental_visit",
                      "pH","Population")

collin.mat <- matrix("", nrow=length(fixed.continuous), ncol=length(fixed.continuous))
rownames(collin.mat) <- fixed.continuous
colnames(collin.mat) <- fixed.continuous


mta <- SLL2.meta[,fixed.continuous]
mta[ mta=="No Sabe/No Contesta" ] <- NA

for (f1 in fixed.continuous) {
  for (f2 in fixed.continuous) {
    ct <- cor.test(as.numeric(mta[,f1]), as.numeric(mta[,f2]))
    
    # if (ct$p.value*length(fixed.continuous) < 0.05 & f1 != f2) {
    #   collin.mat[ f1, f2 ] <- round(ct$p.value*length(fixed.continuous), 5)
    # }
    if (ct$p.value < 0.05 & f1 != f2) {
      collin.mat[ f1, f2 ] <- round(ct$p.value, 5)
    }
  }
}
# ************************** #

# ************************** #
# then check between categorical and continuous variables
fixed.categorical <- c("Gender","Municipal_zone","Moisture_in_home",
                       "Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds","Pets.Reptiles_amphibians",
                       "Smoker","Water_type_home",
                       "Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement",
                       "Mouth_wounds","Medications","Antibiotics","Analgesics","Vitamin_supplements",
                       "Asthma","Wheezing","How_do_you_feel","Do_you_feel_well",
                       "Allergy",
                       "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
                       "Kissing_partner")

collin.cat <- matrix("", nrow=length(fixed.continuous), ncol=length(fixed.categorical))
rownames(collin.cat) <- fixed.continuous
colnames(collin.cat) <- fixed.categorical

mta <- SLL2.meta[,c(fixed.continuous, fixed.categorical)]
mta[ mta=="No Sabe/No Contesta" ] <- NA

for (f1 in fixed.continuous) {
  for (f2 in fixed.categorical) {
    kt <- kruskal.test(as.numeric(mta[,f1]), as.factor(mta[,f2]))
    
    # if (ct$p.value*length(fixed.continuous) < 0.05 & f1 != f2) {
    #   collin.cat[ f1, f2 ] <- round(ct$p.value*length(fixed.continuous), 5)
    # }
    if (kt$p.value < 0.05 & f1 != f2) {
      collin.cat[ f1, f2 ] <- round(kt$p.value, 5)
    }
  }
}

# ************************** #

# ************************** #
# finally check between categorical variables

collin.chi <- matrix("", nrow=length(fixed.categorical), ncol=length(fixed.categorical))
rownames(collin.chi) <- fixed.categorical
colnames(collin.chi) <- fixed.categorical

mta <- SLL2.meta[,fixed.categorical]
mta[ mta=="No Sabe/No Contesta" ] <- NA

for (f1 in fixed.categorical) {
  for (f2 in fixed.categorical) {
    cht <- chisq.test(as.factor(mta[,f1]), as.factor(mta[,f2]))
    
    if (cht$p.value < 0.05 & f1 != f2) {
      collin.chi[ f1, f2 ] <- round(cht$p.value, 5)
    }
  }
}

# ************************** #








library(lme4)
library(car)
library(nnet)

div_vars <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
              "Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
              "Stomatotype_a05_GUnifrac","Stomatotype_Bray","Stomatotype_Jaccard",
              "Stomatotype_JSD","Stomatotype_VAW_GUnifrac","Stomatotype_a0_GUnifrac",
              "Stomatotype_Canberra","Stomatotype_Aitchison")

# ************************************************************************ #
# ************************************************************************ #




# # make a table without "No Sabe/No Contesta"
# metaTab <- SLL2.meta
# metaTab[ metaTab=="No Sabe/No Contesta"] <- NA
# # Add a column that has all 1 value to force all to be same shape by default
# metaTab$Project <- rep("SLL2", nrow(metaTab))
# metaTab$Weight                     <- as.numeric(metaTab$Weight)
# metaTab$Height                     <- as.numeric(metaTab$Height)
# metaTab$Dental.Fillings            <- as.numeric(metaTab$Dental.Fillings)
# metaTab$Dental.Nerve_extractions   <- as.numeric(metaTab$Dental.Nerve_extractions)
# metaTab$Dental.Lost_teeth          <- as.numeric(metaTab$Dental.Lost_teeth)
# metaTab$Dental.Reconstructed_teeth <- as.numeric(metaTab$Dental.Reconstructed_teeth)
# metaTab$pH[ metaTab$pH %in% c(7.55,7.6) ] <- "7.5"
# metaTab$pH <- as.numeric(metaTab$pH)  




# Anova.pvals <- list()
# 
# # these are the fixed effects that will be included each time
# fixeds <- c("Chronic_disorder","Gender","Age","Population","seqGroup")#,"School_ID.letter")#
# 
# # ******************************************** #
# for (nor in c("clr","rel")) {
# 
#   Anova.pvals[[ nor ]] <- list()
# 
#   for (tl in c("contVar","Phylum","Class","Order","Family","Genus")) {
# 
#     # ******************** #
#     if (tl == "contVar") {
#       # no need to run it twice, will just place "contVar" in the "rel" list
#       if ( nor=="rel") {
#         next
#       }
#       glomTab <- NULL
#       dv <- div_vars
#       # ******************** #
#     } else {
#       if (nor=="clr")
#         glomTab <- gloms_clr[[ tl ]]
#       else if (nor=="rel")
#         glomTab <- gloms_rel[[ tl ]]
#       dv <- NULL
#     }
#     # ******************** #
# 
#     print(sprintf("%s - %s", nor, tl))
#     Anova.pvals[[ nor ]][[ tl ]] <- list()
# 
# 
#     # then will go through each of these variables and add them to the list of fixed effects
#     # for (f in c(adonis_vars.short,cont_water_data,
#     for (f in c(adonis_vars,
#                 "Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
#                 "Stomatotype_a05_GUnifrac","Stomatotype_Bray","Stomatotype_Jaccard",
#                 "Stomatotype_JSD","Stomatotype_VAW_GUnifrac","Stomatotype_a0_GUnifrac",
#                 "Stomatotype_Canberra","Stomatotype_Aitchison")) {
#     # for (f in c("Cystic_fibrosis","Transplant","Immune_issues")) {
#       # for (f in c("Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome")) {#,"Stomatotype_Aitchison")) {
#     # for (f in c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats",
#     #             "Pets.Number.Small_furry_animals","Pets.Number.Birds","Pets.Number.Reptiles_amphibians",
#     #             "Pets.Number.Fish","Number_smokers_home","Dental.Fillings","Dental.Nerve_extractions",
#     #             "Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")) {
# 
#       if (tl == "contVar" & startsWith(f, "Stomatotype_"))
#         # only include Stomatotype groups as a fixed effect against taxa, not diversity variables
#         next
# 
#       print(sprintf("  %s", f))
# 
#       if (f %in% c("Municipal_zone",cont_water_data, group_water_data, "Community","Province")) {
#         # for these, dont include Population bc ought to be strongly correlated
#         Anova.pvals[[ nor ]][[ tl ]][[ f ]] <- get_lm( c(f, fixeds[ fixeds != "Population" ]), tl, SLL2.meta,
#                                                        glomTab, dv )
#       } else if (f %in% c("Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome",
#                         "Gingivitis_periodontitis","Eating_disorder",
#                         additional_diseases[ ! additional_diseases %in% c("Samp_collect_issues","Birth_control",
#                                                                           "Antihistamines")])) {
#         if (f %in% c("Cystic_fibrosis","Transplant","Immune_issues")) {
#           # Dont include Chronic_disorder here because these values are all "Yes" in these cases
#           # also include Antibiotics in this case
#           Anova.pvals[[ nor ]][[ tl ]][[ f ]] <- get_lm( c(f, "Antibiotics", fixeds[ fixeds != "Chronic_disorder" ]), tl, SLL2.meta,
#                                                          glomTab, dv )
#         } else {
#           # Dont include Chronic_disorder here because these values are all "Yes" in these cases
#           Anova.pvals[[ nor ]][[ tl ]][[ f ]] <- get_lm( c(f, fixeds[ fixeds != "Chronic_disorder" ]), tl, SLL2.meta,
#                                                          glomTab, dv )
#         }
#         
#       } else {
#         # otherwise just use the normal vector, when f is one of the fixeds, put that first
#         Anova.pvals[[ nor ]][[ tl ]][[ f ]] <- get_lm( c(f, fixeds[ fixeds != f ]), tl, SLL2.meta, glomTab, dv )
#       }
#     }
# 
#   }
# 
# }
# 
# 
# saveRDS(Anova.pvals, sprintf("%s/figures/LMM_signifs/Anova.pvals.rds", p2_dir))
# # ******************************************** #

# Anova.pvals <- readRDS(sprintf("%s/figures/LMM_signifs/Anova.pvals.rds", p2_dir))
Anova.pvals.unadjusted <- readRDS(sprintf("%s/figures/LMM_signifs/Anova.pvals_unadjusted.rds", p2_dir))

# # ******************************************** #
# Anovas for pH ####
# to get supplementary figure S7
# anova.pH <- get_lm( c("pH","Chronic_disorder","Gender","Age","Population","seqGroup"), 
#                     "contVar", SLL2.meta, NULL, 
#                     c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
#                     noRemove = T, rerun.nonSig = F)

anova.healthy_pH <- get_lm( c("pH","Gender","Age","Population","seqGroup"), 
                            "contVar", meta.healthy, NULL, 
                            c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
                            noRemove = T, rerun.nonSig = F)

# anova.pH_label.wide <- get_lm( c("pH_label.wide","Chronic_disorder","Gender","Age","Population","seqGroup"), 
#                                "contVar", SLL2.meta, NULL, 
#                                c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
#                                noRemove = T, rerun.nonSig = F)

anova.healthy_pH_label.wide <- get_lm( c("pH_label.wide","Gender","Age","Population","seqGroup"), 
                                       "contVar", meta.healthy, NULL, 
                                       c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
                                       noRemove = T, rerun.nonSig = F)
# ******* #
meta.pH_no_5 <- meta.healthy[ ! is.na(meta.healthy$pH) & meta.healthy$pH > 5.5, ]

anova.pH_no_5 <- get_lm( c("pH","Gender","Age","Population","seqGroup"), 
                         "contVar", meta.pH_no_5, NULL, 
                         c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
                         noRemove = T, rerun.nonSig = F)

anova.pH_no_5.label.wide <- get_lm( c("pH_label.wide","Gender","Age","Population","seqGroup"), 
                                    "contVar", meta.pH_no_5, NULL, 
                                    c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts"),
                                    noRemove = T, rerun.nonSig = F)
# ******* #


# subsampling.plots.box("pH","contVar","pH", 
#                       c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
#                       c(groupQs,"pH"), only_cont, 
#                       gloms_clr, SLL2.meta, SLL2, NULL, 
#                       ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(SLL2.meta),
#                       plotType = "box", plot_tukey = F)

subsampling.plots.box("pH","contVar","pH", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      c(groupQs,"pH"), only_cont, 
                      gloms_clr, meta.healthy, SLL2, NULL, 
                      ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(meta.healthy),
                      plotType = "box", plot_tukey = F)


# subsampling.plots.box("pH_label.wide","contVar","pH_label.wide", 
#                       c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
#                       c(groupQs,"pH_label.wide"), only_cont, 
#                       gloms_clr, SLL2.meta, SLL2, NULL, 
#                       ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(SLL2.meta),
#                       plotType = "box", plot_tukey = F, xAngle = 30 )

subsampling.plots.box("pH_label.wide","contVar","pH_label.wide", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      c(groupQs,"pH_label.wide"), only_cont, 
                      gloms_clr, meta.healthy, SLL2, NULL, 
                      ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(meta.healthy),
                      plotType = "box", plot_tukey = F, xAngle = 30 )



subsampling.plots.box("pH","contVar","pH", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      c(groupQs,"pH"), only_cont, 
                      gloms_clr, meta.pH_no_5, SLL2, NULL, 
                      ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(meta.pH_no_5),
                      plotType = "box", plot_tukey = F)


subsampling.plots.box("pH_label.wide","contVar","pH_label.wide", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      c(groupQs,"pH_label.wide"), only_cont, 
                      gloms_clr, meta.pH_no_5, SLL2, NULL, 
                      ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(meta.pH_no_5),
                      plotType = "box", plot_tukey = F, xAngle = 30 )



meta.pH_6_8 <- meta.healthy[ ! is.na(meta.healthy$pH) & meta.healthy$pH %in% c(6, 8), ]

subsampling.plots.box("pH","contVar","pH", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      c(groupQs,"pH"), only_cont, 
                      gloms_clr, meta.pH_6_8, SLL2, NULL, 
                      ignore_pvals = T, ignore_dstStruc = T, chosenSamps = rownames(meta.pH_6_8),
                      plotType = "box", plot_tukey = F)
# # ******************************************** #

# check taxa that are significant by pH:
Anova.pvals.unadjusted$clr$Genus$pH_no_5 <- get_lm( c("pH","Chronic_disorder","Gender","Age","Population","seqGroup"), 
                                                    "Genus", meta.pH_no_5, 
                                                    gloms_clr$Genus[ , rownames(meta.pH_no_5) ], NULL)

ph.taxa <- Anova.pvals$clr$Genus$pH_no_5
ph.taxa <- apply(ph.taxa, 2, as.numeric)
rownames(ph.taxa) <- rownames(Anova.pvals$clr$Genus$pH_no_5)

sort(ph.taxa[,"pH"])

group_vs_cont_box(meta.pH_no_5, rownames(meta.pH_no_5), "Genus", gloms_clr$Genus[ , rownames(meta.pH_no_5)], "pH", 
                  "Porphyromonas", c(groupQs,"pH"), only_cont, print_tukey = F, plot_tukey = F )
# # ******************************************** #




adjust_anova_pvals <- function(anovaPTab) {
  anovaPTab[ anovaPTab== "" ] <- 1
  # adjust pvalues
  pval.mat <- apply(anovaPTab, 2, p.adjust, method="fdr")
  
  # print only significant values
  pval.mat[ pval.mat > 0.05 ] <- ""
  pval.mat[ is.na(pval.mat) ] <- ""

  # remove rows with dependVar values that had no significant effects
  pval.mat <- as.data.frame(pval.mat)[ pval.mat[ , 1] != "", ]
  
  return( as.matrix(pval.mat) )
}



Anova.pvals <- list()

for (nor in names(Anova.pvals.unadjusted)) {
  Anova.pvals[[ nor ]] <- list()
  for (tl in names(Anova.pvals.unadjusted[[ nor ]])) {
    Anova.pvals[[ nor ]][[ tl ]] <- list()
    for (variable in names(Anova.pvals.unadjusted[[ nor ]][[ tl ]])) {
      if (nrow(Anova.pvals.unadjusted[[ nor ]][[ tl ]][[ variable ]]) %in% c(0,1)) 
        Anova.pvals[[ nor ]][[ tl ]][[ variable ]] <- Anova.pvals.unadjusted[[ nor ]][[ tl ]][[ variable ]]
      else
        Anova.pvals[[ nor ]][[ tl ]][[ variable ]] <- adjust_anova_pvals(Anova.pvals.unadjusted[[ nor ]][[ tl ]][[ variable ]])
    }
  }
}

# ****************************************************************************************************************** #





# # ******************************************** #
# Fungi vs Genera ####

# cand.anova <- list()
# for (i in names(all_subSamps$Full_MALDI.Candida)) {
#   print(i)
#   meta.cand <- SLL2.meta[ all_subSamps$Full_MALDI.Candida[[ i ]], ]
#   glom.cand <- gloms_clr$Genus[ , all_subSamps$Full_MALDI.Candida[[ i ]] ]
#   
#   cand.anova[[ i ]] <- get_lm( c("Full_MALDI.Candida","Gender","Age","Population"), "Genus", meta.cand, 
#                                glom.cand, NULL,
#                                noRemove = T, rerun.nonSig = F)
# }
# saveRDS(cand.anova, sprintf("%s/R_objects/cand.anova.rds", p2_dir))

cand.anova <- readRDS(sprintf("%s/R_objects/cand.anova.rds", p2_dir))


gen.pMeans.cand   <- sapply(rownames(cand.anova$`1`), function(x)
  mean(p.adjust(unlist(lapply(cand.anova, function(y)
    y[ gsub("-","\\.",x), "Full_MALDI.Candida" ]
  )), method = "fdr"))
)

gen.num_sig.cand   <- sapply(rownames(cand.anova$`1`), function(x)
  sum(p.adjust(unlist(lapply(cand.anova, function(y)
    y[ gsub("-","\\.",x), "Full_MALDI.Candida" ]
  )), method = "fdr") < 0.05 )
)


sort(gen.num_sig.cand)
sort(gen.pMeans.cand, decreasing = T)



# ******************** #
# cand.kw.p <- list()
# cand.kw.s <- list()
# for (i in names(all_subSamps$Full_MALDI.Candida)) {
#   print(i)
#   meta.cand <- SLL2.meta[ all_subSamps$Full_MALDI.Candida[[ i ]], ]
#   glom.cand <- gloms_clr$Genus[ , all_subSamps$Full_MALDI.Candida[[ i ]] ]
# 
#   cand.kw.p[[ i ]] <- sapply(rownames(gloms_clr$Genus), function(gen)
#     kruskal.test( as.numeric(glom.cand[ gen, ]),
#                   as.factor(meta.cand$Full_MALDI.Candida))$p.value
#   )
#   
#   cand.kw.s[[ i ]] <- sapply(rownames(gloms_clr$Genus), function(gen)
#     kruskal.test( as.numeric(glom.cand[ gen, ]),
#                   as.factor(meta.cand$Full_MALDI.Candida))$statistic
#   )
# }
# saveRDS(cand.kw.p, sprintf("%s/R_objects/cand.kw.p.rds", p2_dir))
# saveRDS(cand.kw.s, sprintf("%s/R_objects/cand.kw.s.rds", p2_dir))

cand.kw.p <- readRDS(sprintf("%s/R_objects/cand.kw.p.rds", p2_dir))
cand.kw.s <- readRDS(sprintf("%s/R_objects/cand.kw.s.rds", p2_dir))

gen.kw.pMeans.cand <- sapply(names(cand.kw.p$`1`), function(x)
  mean(p.adjust(unlist(sapply(cand.kw.p, function(y)
    y[ x ]
  )), method = "fdr"))
)

gen.kw.num_sig.cand <- sapply(names(cand.kw.p$`1`), function(x)
  sum(p.adjust(unlist(sapply(cand.kw.p, function(y)
    y[ x ]
  )), method = "fdr") < 0.05 )
)


sort(gen.kw.num_sig.cand)
sort(gen.kw.pMeans.cand, decreasing = T)

# ****************************************************************************************************************** #








# ****************************************************************************************************************** #
# ****************************************************************************************************************** #

c("Weight","Height","Years_in_home","Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
  "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish","Number_smokers_home",
  "Dental.Fillings","Dental.Nerve_extractions","Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")













c("Acinetobacter","Actinobacillus","Aggregatibacter","Alloprevotella",
  "Bacteroides","Bradyrhizobium","Brevundimonas","Capnocytophaga","Chryseobacterium",
  "Comamonas","Corynebacterium","Curvibacter","Desulfobulbus","Dialister","Hydrogenophaga",
  "Hyphomicrobium","Mesorhizobium","Mobilunculus","Moraxella","Oceanivirga","Oribacterium",
  "Peptococcus","Peptostreptococcus","Phocaeicola","Prevotella","Pseudomonas","Ralstonia",
  "Shuttleworthia","Simonsiella","Solobacterium","Sphingobacterium","Staphylococcus",
  "Stenotrophomonas","Streptococcus","Treponema","Variovorax")

c("Abiotrophia","Anaetoglobulus","Atopobium","Bulleidia","Cardiobacterium","Eikenella",
  "Fretibacterium","Fusobacterium","Gemella","Granulicatella","Haemophilus","Johnsonella",
  "Kingella","Leptotrichia","Megasphaera","Microbacterium","Mogibacterium","Mycoplasma","Neisseria",
  "Parvimonas","Porphyromonas","Pseudopropionibacterium","Rothia","Streptobacillus","Veillonella")



c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Spirochaetes","unclassified.P1")

c("Epsilonbacteraeota","Fusobacteria")


# metaTab <- SLL2.meta
# metaTab[ metaTab=="No Sabe/No Contesta"] <- NA
# # Add a column that has all 1 value to force all to be same shape by default
# metaTab$Project <- rep("SLL2", nrow(metaTab))
# metaTab$Weight                     <- as.numeric(metaTab$Weight)
# metaTab$Height                     <- as.numeric(metaTab$Height)
# metaTab$Dental.Fillings            <- as.numeric(metaTab$Dental.Fillings)
# metaTab$Dental.Nerve_extractions   <- as.numeric(metaTab$Dental.Nerve_extractions)
# metaTab$Dental.Lost_teeth          <- as.numeric(metaTab$Dental.Lost_teeth)
# metaTab$Dental.Reconstructed_teeth <- as.numeric(metaTab$Dental.Reconstructed_teeth)
# metaTab$pH[ metaTab$pH %in% c(7.55,7.6) ] <- "7.5"
# metaTab$pH <- as.numeric(metaTab$pH)  


# ************************************************************************ #
# Make plots for significant LMM results ####

lmCounter <- list("clr"=0, "rel"=0)

for (nor in names(Anova.pvals)) {
  
  # ******************** #
  if (nor=="clr")
    glomTab <- gloms_clr
  else if (nor=="rel")
    glomTab <- gloms_rel
  # ******************** #
  
  for (tl in names(Anova.pvals[[ nor ]])) {
    print(sprintf("%s - %s", nor, tl))
    
    # ******************** #
    if (tl == "contVar") {
      # arbitrarily make tlev Genus, this value wont be used in this case
      tlev <- "Genus"
      comparison <- "questions"
    } else {
      tlev <- tl
      comparison <- "TvsQ"
    }
    # ******************** #
    
    for (f in names(Anova.pvals[[ nor ]][[ tl ]])) {
    # for (f in c("Weight","Height","Years_in_home",
    #             "Pets.Number.Dogs","Pets.Number.Cats","Pets.Number.Small_furry_animals",
    #             "Pets.Number.Birds","Pets.Number.Reptiles_amphibians","Pets.Number.Fish",
    #             "Number_smokers_home","Dental.Fillings","Dental.Nerve_extractions",
    #             "Dental.Lost_teeth","Dental.Reconstructed_teeth","pH")) {
    # for (f in cont_water_data) {
    # for (f in c("Age","Brushing","Floss")) {
    # for (f in c("Consumption.Milk","Consumption.Yogurt","Consumption.Sweets","Consumption.Chewing_gum",
    #             "Consumption.Nuts","Last_dental_visit","pH")) {
    # for (f in c("Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome")) {
    # for (f in c("Braces","Water_type_home","Municipal_zone","How_do_you_feel")) {
      
      
      print(sprintf("  %s", f))
      
      for (div_or_tax in rownames(Anova.pvals[[ nor ]][[ tl ]][[ f ]])) {
        # print(sprintf("    %s", div_or_tax))
        
        # only plot for those div_or_tax that f was significant
        if (Anova.pvals[[ nor ]][[ tl ]][[ f ]][ div_or_tax, f ] != "") {
          
          lmCounter[[ nor ]] <- lmCounter[[ nor ]] + 1
          
          # ****************************************************************** #
          # have to change the characters that needed to be changed during lm tests above back to original format
          if (div_or_tax == "Absconditabacteriales_.SR1.")
            div_or_tax <- "Absconditabacteriales_(SR1)"
          else if ( tl != "contVar")
            div_or_tax <- gsub("unclassified-", "unclassified.", 
                               gsub("Family_XI-", "Family_XI.", 
                                    gsub("\\.", "-", div_or_tax) ))
          
          # ****************************************************************** #
          if (startsWith(f, "Stomatotype") ) {
            # do both gradient on PCoA plot...
            plot_gradients(tl, SLL2, gsub("Stomatotype_","",f), pcoas, glomTab, SLL2.meta, div_or_tax,
                           save.grads = T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f), 
                           anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
            # as well as boxplot separating groups
            # in this case f is the group_col and div_or_tax is the cont_col
            group_vs_cont_box(SLL2.meta, sample_names(SLL2_rel), tlev, glomTab[[ tlev ]], 
                              f, div_or_tax, groupQs, only_cont, print_tukey = F, plot_tukey = T,
                              save.boxes = T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                              anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
            # ********************** #
          } else if (f %in% c(groupQs,"Brushing","Floss","Consumption.Milk","Consumption.Yogurt",
                              "Consumption.Sweets","Consumption.Chewing_gum","Consumption.Nuts",
                              "Last_dental_visit","pH")) {
            if ( ! startsWith(div_or_tax, "Stomatotype")) {
              # for now, not plotting non taxa with Stomatotypes ... as will need to do an assoc plot
              #    or represent it in a PCoA plot, which already basically did with adonis tests
              group_vs_cont_box(SLL2.meta, sample_names(SLL2_rel), tlev, glomTab[[ tlev ]], 
                                f, div_or_tax, 
                                c(groupQs,"Brushing","Floss","Consumption.Milk","Consumption.Yogurt",
                                  "Consumption.Sweets","Consumption.Chewing_gum","Consumption.Nuts",
                                  "Last_dental_visit","pH"), 
                                only_cont, print_tukey = F, plot_tukey = T,
                                xAngle = ifelse(f %in% c("Braces","Water_type_home","Municipal_zone","How_do_you_feel",
                                                         "Community","Province","Composition","Hardness_category"), 
                                                15, 0), # angle x labels so they are readable
                                save.boxes = T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                                anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
            }
            # ********************** #
          } else if (startsWith(div_or_tax, "Stomatotype")) {
            # in this case, represent the continuous f variables with the gradient plot
            plot_gradients(tl, SLL2, gsub("Stomatotype_","",div_or_tax), pcoas, glomTab, SLL2.meta, f,
                           save.grads = T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                           anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
            # as well as boxplot separating groups
            # in this case div_or_tax is the group_col and f is the cont_col
            group_vs_cont_box(SLL2.meta, sample_names(SLL2_rel), tlev, glomTab[[ tlev ]], 
                              div_or_tax, f, groupQs, only_cont, print_tukey = F, plot_tukey = T,
                              save.boxes = T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                              anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
            # also use the adonis plot function bc better at displaying the range of continuous values
            colorLogs <- c("Neutrophils","Lymphocytes","IL.1β","IL.6","IL.8","IL.10","IL.17A","IL.18","IL.23")
            
            plot_adonis_w_covars(gsub("Stomatotype_","",div_or_tax), SLL2, SLL2.meta, glomTab, 
                                 shapeBy = div_or_tax, colorBy = f, tl = ifelse(tl=="contVar",NA,tlev), 
                                 anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]], pcoaList = pcoas, 
                                 log_of_color = ifelse(f %in% colorLogs, T, F), 
                                 save.adoPlot = T)
            # ********************** #
          } else {
            plot_data.cont(f, div_or_tax, comparison, tlev, sample_names(SLL2), only_cont, glomTab, SLL2.meta, 
                           save.scats=T, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
                           anova.pTab = Anova.pvals[[ nor ]][[ tl ]][[ f ]] )
          }
          # ****************************************************************** #
          
          
        }
      }
    }
  }
  
}


# ************************************************************************ #

f <- "Community"
tl <- tlev <- "Genus"

for (div_or_tax in cont_water_data) {
  
  print(div_or_tax)
  
  group_vs_cont_box(SLL2.meta, sample_names(SLL2_rel), tlev, glomTab[[ tlev ]], 
                    f, div_or_tax, 
                    c(groupQs,"Brushing","Floss","Consumption.Milk","Consumption.Yogurt",
                      "Consumption.Sweets","Consumption.Chewing_gum","Consumption.Nuts",
                      "Last_dental_visit","pH"), 
                    only_cont, print_tukey = F, plot_tukey = F,
                    xAngle = ifelse(f %in% c("Braces","Water_type_home","Municipal_zone","How_do_you_feel",
                                             "Community","Province","Composition","Hardness_category"), 
                                    15, 0), # angle x labels so they are readable
                    save.boxes = T, plot_dirs = sprintf("LMM_signifs/%s/%s/Water_data", nor, f))
  
}


# ************************************************************************ #
















# ****************************************************************************************************************** #
# ANOSIM for some categorical variables ####

# If you have very different group sizes, you may consider analysis of similarities (ANOSIM) instead of PERMANOVA. 
# This test does not assume equal group variances. However, it only allows simple 1 variable models with no 
# interactions and can only be used for categorical (AgeGroup), not continuous (ADG) variables. So, ANOSIM has a 
# lot of limitations and should only be used if you group sizes are very, very different, like 10 vs 100.

anosim_list <- list()
for (dist_meas in c("jsd","weighted_Unifrac","unweighted_Unifrac",
                    "guni.VAW","guni.a0","guni.a05",
                    "bray","jaccard","canberra","aitch")) {
  
  anosim_list[[ dist_meas ]] <- list()
  
  for (unit_lab in c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
                     "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
                     "Grandparent_Grandchild_unit")) {
    print(c(dist_meas, unit_lab))
    
    
    anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas)) ),
                                                       as.character(as.matrix(SLL2.meta[, unit_lab])),
                                                       strata = SLL2.meta[ , "seqGroup"])
    
    # mTab <- SLL2.meta[ SLL2.meta[,unit_lab] != "None", ]
    # anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas))[rownames(mTab),rownames(mTab)] ), 
    #                                                    as.character(as.matrix(mTab[, unit_lab])))
  }
}

saveRDS(anosim_list, sprintf("%s/R_objects/anosim_list.rds", p2_dir))

anosim_list <- readRDS(sprintf("%s/R_objects/anosim_list.rds", p2_dir))
# ****************************************************************************************************************** #

















# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# test plotting with non-centered minimum curvilinear embedding ####

# ********************************************************************** #
# function taken from here: 
#   https://sites.google.com/site/carlovittoriocannistraci/5-datasets-and-matlab-code/minimum-curvilinearity-ii-april-2012
library(igraph)
mce <- function(x, n, centring){
  
  # This code is performs just MCE and ncMCE (both SVD-based)
  
  #Given a distance or correlation matrix x, it performs Minimum Curvilinear 
  #Embedding (MCE) or non-centred MCE (ncMCE) (coded 24-March-2013 by 
  #Gregorio Alanis-Lobato and checked by Carlo Vittorio Cannistraci)
  
  #INPUT
  #   x => Distance (example: Euclidean) or distance-adjusted correlation matrix (example: x = 1 - Pearson_correlation)
  #   n => Dimension into which the data is to be embedded
  #   centring => 'yes' if x should be centred or 'no' if not
  #OUTPUT
  #   s => Coordinates of the samples being embedded into the reduced n-dimesional space
  
  # #Make sure the required library 'igraph' is installed and load it
  # if(require("igraph")){
  #   print("igraph has been loaded...");
  # } else{
  #   print("Trying to install igraph...");
  #   install.packages("igraph");
  #   if(require("igraph")){
  #     print("igraph has been installed and loaded...");
  #   } else{
  #     stop("Could not install igraph");
  #   }
  # }
  
  #Make sure the matrix is symmetric
  x <- pmax(x, t(x));
  
  #Create a graph object out of the adjacency matrix x
  g <- graph.adjacency(x, mode = "undirected", weighted = TRUE);
  
  #MC-kernel computation
  mst <- minimum.spanning.tree(g);
  kernel <- shortest.paths(mst);
  
  #Kernel centring
  if(centring == "yes"){
    N <- nrow(kernel);
    J <- diag(N) - (1/N)*matrix(1, N, N); #Form the centring matrix J
    kernel <- (-0.5)*(J %*% kernel^2 %*% J);
  }
  
  #SVD-based Embedding
  res <- svd(kernel);
  L <- diag(res$d);
  V <- res$v;
  
  sqrtL <- sqrt(L[1:n, 1:n]);
  V <- V[, 1:n];
  
  s <- t(sqrtL %*% t(V));
  
  return(s);
}
# ********************************************************************** #

ncMCEs <- list()

ncMCEs[[ "Aitchison" ]]          <- mce(aitch, 3, "no")
ncMCEs[[ "Weighted_Unifrac" ]]   <- mce(as.matrix(weighted_Unifrac), 3, "no")
ncMCEs[[ "Unweighted_Unifrac" ]] <- mce(as.matrix(unweighted_Unifrac), 3, "no")
ncMCEs[[ "Bray" ]]               <- mce(bray, 3, "no")
ncMCEs[[ "Jaccard" ]]            <- mce(jaccard, 3, "no")


MCEs <- list()

MCEs[[ "Aitchison" ]]          <- mce(aitch, 3, "yes")
MCEs[[ "Weighted_Unifrac" ]]   <- mce(as.matrix(weighted_Unifrac), 3, "yes")
MCEs[[ "Unweighted_Unifrac" ]] <- mce(as.matrix(unweighted_Unifrac), 3, "yes")
MCEs[[ "Bray" ]]               <- mce(bray, 3, "yes")
MCEs[[ "Jaccard" ]]            <- mce(jaccard, 3, "yes")




























# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####
# Adonis and ANOVA for disorders with matched control subsamplings ####


# for (i in names(disSubsTests$default$Cystic_fibrosis)) {
#   print(i)
#   soi <- disSubsTests$default$Cystic_fibrosis[[ i ]]$samples
#   mta <- disSubsTests$default$Cystic_fibrosis[[ i ]]$mTab
#   
#   disSubsTests$default$Cystic_fibrosis[[ i ]]$Anova[[ "Species" ]] <- get_lm( c("Cystic_fibrosis", "Antibiotics","Gender","Age","Population"), 
#                                                                               "Species", mta, gloms_clr$Species[ , soi ], 
#                                                                               NULL, noRemove = T, rerun.nonSig = F)
# }

# # ************************************************************************ #
# # ************************************************************************ #
# 
# disSubsTests <- list()
# 
# # for (sampMode in c("default","sameComs","equalSize")) {
# for (sampMode in c("default","sameComs")) {
#   
#   disSubsTests[[ sampMode ]] <- list()
#   
#   for (disorder in c("Cystic_fibrosis","Downs_Syndrome")) {
#     
#     disSubsTests[[ sampMode ]][[ disorder ]] <- run_full_subsampling_calcs.disorder(disorder, SLL2, SLL2.meta, 100, sampMode)
#   }
# }
# 
# saveRDS(disSubsTests, file = sprintf("%s/R_objects/disSubsTests.rds", p2_dir))
# # ************************************************************************ #

disSubsTests <- readRDS(sprintf("%s/R_objects/disSubsTests.rds", p2_dir))
# ************************************************************************ #


# disSubsTests.2 <- list()
# disSubsTests.2[[ "default" ]] <- list()
# 
# for (disorder in c("Cystic_fibrosis","Downs_Syndrome")) {
# 
#   disSubsTests.2[[ "default" ]][[ disorder ]] <- list()
# 
#   # ******************************** #
#   if (disorder == "Cystic_fibrosis") {
#     cov_vars <- c("Antibiotics","Gender","Age","Population")#
#   } else {
#     cov_vars <- c("Gender","Age","Population")#
#   }
#   # ******************************** #
# 
#   for (i in names(disSubsTests$default[[ disorder ]])) {
# 
#     disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]] <- list()
# 
#     disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "samples" ]] <- disSubsTests[[ "default" ]][[ disorder ]][[ i ]][[ "samples" ]]
#     disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "mTab" ]] <- disSubsTests[[ "default" ]][[ disorder ]][[ i ]][[ "mTab" ]]
#     disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "Adonis" ]] <- disSubsTests[[ "default" ]][[ disorder ]][[ i ]][[ "Adonis" ]]
# 
#     samples.to_use <- disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "samples" ]]
#     mTab.to_use <- disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "mTab" ]][ samples.to_use, ]
#     # ************************************************************ #
#     # run linear model for subsamples
#     disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "Anova" ]] <- list()
# 
#     for (tl in c("contVar","Phylum","Class","Order","Family","Genus")) {
# 
#       # ******************** #
#       if (tl == "contVar") {
#         glomTab <- NULL
#         dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
#                 "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
#                 "Stomatotype_Bray","Stomatotype_Jaccard")
#         # ******************** #
#       } else {
#         glomTab <- gloms_clr
#         dv <- NULL
#       }
#       # ******************** #
# 
#       disSubsTests.2[[ "default" ]][[ disorder ]][[ i ]][[ "Anova" ]][[ tl ]] <- get_lm( c(disorder, cov_vars),
#                                                  tl, mTab.to_use, glomTab[[ tl ]][ , samples.to_use ],
#                                                  dv, noRemove = T, rerun.nonSig = F)
#     }
#     # ************************************************************ #
#   }
# 
# }
# saveRDS(disSubsTests.2, file = sprintf("%s/R_objects/disSubsTests.2.rds", p2_dir))


# ************************************************************************ #

gen.pMeans.ds.default <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Downs_Syndrome"]
    }
    )))
  )
phy.pMeans.ds.default <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )))
)
conts.pMeans.ds.default <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )))
)


gen.num_sig.ds.default <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )) < 0.05)
)
phy.num_sig.ds.default <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )) < 0.05)
)
conts.num_sig.ds.default <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )) < 0.05)
)

# gen.pMeans.ds.sameComs <- sapply(rownames(gloms_clr$Genus), function(x)
#   mean(unlist(lapply(disSubsTests$sameComs$Downs_Syndrome, function(y) {
#     y$Anova$Genus[gsub("-","\\.",x),"Downs_Syndrome"]
#   }
#   )))
# )
# gen.pMeans.ds.equalSize <- sapply(rownames(gloms_clr$Genus), function(x)
#   mean(unlist(lapply(disSubsTests$equalSize$Downs_Syndrome, function(y) {
#     y$Anova$Genus[gsub("-","\\.",x),"Downs_Syndrome"]
#   }
#   )))
# )




gen.pMeans.cf.default <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)
phy.pMeans.cf.default <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)
conts.pMeans.cf.default <- sapply(rownames(disSubsTests$default$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)


gen.num_sig.cf.default <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)
phy.num_sig.cf.default <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)
conts.num_sig.cf.default <- sapply(rownames(disSubsTests$default$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)





# ************************************************************************ #
# mean and num sig for each adonis test ####
ado.Pmeans.ds <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                        function(x) mean(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) 
                          y$Adonis[[ x ]]["Downs_Syndrome","Pr(>F)"]))))

ado.Psds.ds <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                      function(x) sd(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) 
                        y$Adonis[[ x ]]["Downs_Syndrome","Pr(>F)"]))))

ado.num_sig.ds <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                         function(x) sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) 
                           y$Adonis[[ x ]]["Downs_Syndrome","Pr(>F)"])) < 0.05))




ado.Pmeans.cf <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                        function(x) mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) 
                          y$Adonis[[ x ]]["Cystic_fibrosis","Pr(>F)"]))))

ado.Psds.cf <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                      function(x) sd(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) 
                        y$Adonis[[ x ]]["Cystic_fibrosis","Pr(>F)"]))))

ado.num_sig.cf <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                         function(x) sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) 
                           y$Adonis[[ x ]]["Cystic_fibrosis","Pr(>F)"])) < 0.05))


CF.mC <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
CF.mC <- CF.mC[ SLL2.meta[CF.mC, "Cystic_fibrosis"]=="No" ]
CF.samps <- rownames(SLL2.meta[ SLL2.meta$Cystic_fibrosis=="Yes", ])

phy.cf.h <- prune_samples(c(CF.samps, CF.mC), SLL2)
mTab.cf.h <- SLL2.meta[ c(CF.samps, CF.mC), ]

mTab.cf.h$CF_ab <- sapply(rownames(mTab.cf.h), function(x)
  ifelse(mTab.cf.h[x, "Cystic_fibrosis"]=="Yes" & mTab.cf.h[x, "Antibiotics"]=="Yes", "CF.Y_Ab.Y",
         ifelse(mTab.cf.h[x, "Cystic_fibrosis"]=="Yes" & mTab.cf.h[x, "Antibiotics"]=="No", "CF.Y_Ab.N",
                ifelse(mTab.cf.h[x, "Cystic_fibrosis"]=="No" & mTab.cf.h[x, "Antibiotics"]=="Yes", "CF.N_Ab.Y", "CF.N_Ab.N"))))

# ordObj.cf.h <- subsampling_ordination_objects(c(CF.samps, CF.mC), phy.cf.h, distsOnly=T)
# saveRDS(ordObj.cf.h, file = sprintf("%s/R_objects/ordObj.cf.h.rds", p2_dir))
ordObj.cf.h <- readRDS(sprintf("%s/R_objects/ordObj.cf.h.rds", p2_dir))

ordObj.cf.h.pcoa <- list()
ordObj.cf.h.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.cf.h$Aitchison)
ordObj.cf.h.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.cf.h$Weighted_Unifrac)
ordObj.cf.h.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.cf.h$Unweighted_Unifrac)
ordObj.cf.h.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.cf.h$Bray)
ordObj.cf.h.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.cf.h$Jaccard)

plot_adonis_w_covars("Jaccard", phy.cf.h, mTab.cf.h, gloms_clr, 
                     shapeBy = "CF_ab", colorBy = "CF_ab", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.cf.h.pcoa)

plot_adonis_w_covars("Aitchison", phy.cf.h, mTab.cf.h, gloms_clr, 
                     shapeBy = "CF_ab", colorBy = "CF_ab", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.cf.h.pcoa)



plot_adonis_w_covars("Aitchison", phy.cf.h, mTab.cf.h, gloms_clr, 
                     shapeBy = "Cystic_fibrosis", colorBy = "Community", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.cf.h.pcoa)

plot_adonis_w_covars("Jaccard", phy.cf.h, mTab.cf.h, gloms_clr, 
                     shapeBy = "Cystic_fibrosis", colorBy = "Community", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.cf.h.pcoa)












# ************************************************************************ #
# Check how often Antibiotics was significant for CF tests ####


gen.pMeans.ab_in_cf <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Antibiotics"]
  }
  )))
)
phy.pMeans.ab_in_cf <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Antibiotics"]
  }
  )))
)
conts.pMeans.ab_in_cf <- sapply(rownames(disSubsTests$default$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Antibiotics"]
  }
  )))
)


gen.num_sig.ab_in_cf <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Antibiotics"]
  }
  )) < 0.05)
)
phy.num_sig.ab_in_cf <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Antibiotics"]
  }
  )) < 0.05)
)
conts.num_sig.ab_in_cf <- sapply(rownames(disSubsTests$default$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Antibiotics"]
  }
  )) < 0.05)
)


# ************************************************************************ #













# ************************************************************************ #
# those genera that had a mean pval < 0.05 and were < 0.05 in at least 70/100 subsamplings
siggen <- sort(gen.pMeans.ds.default[ gen.pMeans.ds.default < 0.05 & gen.num_sig.ds.default > 70 ])
# table with mean pval next to frequency of being significant among 100 subsamplings
cbind( gen.pMeans.ds.default[ names(siggen) ], gen.num_sig.ds.default[ names(siggen) ] )
# ************************************************************************ #


ds.89 <- disSubsTests$sameComs$Downs_Syndrome$`89`$samples
ds.7  <- disSubsTests$sameComs$Downs_Syndrome$`7`$samples

ds.35 <- disSubsTests$sameComs$Downs_Syndrome$`35`$samples
ds.63 <- disSubsTests$sameComs$Downs_Syndrome$`63`$samples

group_vs_cont_box(SLL2.meta, ds.89, "Genus", gloms_clr$Genus, 
                  "Downs_Syndrome", "Kingella", groupQs, only_cont, print_tukey = F, plot_tukey = T,
                  save.boxes = F)

group_vs_cont_box(SLL2.meta, ds.7, "Genus", gloms_clr$Genus, 
                  "Downs_Syndrome", "Kingella", groupQs, only_cont, print_tukey = F, plot_tukey = T,
                  save.boxes = F)

group_vs_cont_box(SLL2.meta, ds.35, "Genus", gloms_clr$Genus, 
                  "Downs_Syndrome", "Kingella", groupQs, only_cont, print_tukey = F, plot_tukey = T,
                  save.boxes = F)

group_vs_cont_box(SLL2.meta, ds.63, "Genus", gloms_clr$Genus, 
                  "Downs_Syndrome", "Kingella", groupQs, only_cont, print_tukey = F, plot_tukey = T,
                  save.boxes = F)
# ************************************************************************ #































# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# Adonis and ANOVA for DS FAMILY MEMBERS with matched control subsamplings ####
# famSubsTests <- list()
# 
# for (disorder in c("Downs_Syndrome","Cystic_fibrosis")) {
#   
#   # get disorder samples
#   disSamps <- rownames( SLL2.meta[ SLL2.meta[ , disorder ] == "Yes", ] )
#   
#   # get family member samples 
#   dis.famUnits <- unique( SLL2.meta[disSamps, "Family_unit"] )
#   dis.famUnits <- dis.famUnits[ dis.famUnits != "None" ]
#   disFams <- rownames( SLL2.meta[ SLL2.meta[, "Family_unit"] %in% dis.famUnits, ])
#   
#   
#   
#   # take matched controls, fam members, remove the DS samples
#   matchCont <- lapply(disSubsTests$default[[ disorder ]], function(x) {
#     mC <- x$samples
#     mC <- unique(c( mC, disFams ))
#     mC <- mC[ SLL2.meta[ mC, disorder ] == "No"]
#     })
#   
#   # run tests using these same matched controls with the fam members now instead of disorder
#   famSubsTests[[ disorder ]] <- run_full_subsampling_calcs.disorder(disorder, SLL2, SLL2.meta, 100,
#                                                                     sampMode = "family", chosenControls = matchCont)
#   
# }
#
# saveRDS(famSubsTests, file = sprintf("%s/R_objects/famSubsTests.rds", p2_dir))
# # ************************************************************************ #

famSubsTests <- readRDS(sprintf("%s/R_objects/famSubsTests.rds", p2_dir))
# ************************************************************************ #


# write questionnaires to tables with only fam samples to see who the family members are, etc
dsf <- SLL2.meta[ SLL2.meta$DS_family=="Yes",]
DS.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Downs_Syndrome=="Yes"])
# order the samples to have DS first, then their families after
dsf <- dsf[c(rownames(dsf)[ rownames(dsf) %in% DS.samps], rownames(dsf)[ ! rownames(dsf) %in% DS.samps]), ]
write.csv(dsf, file = sprintf("%s/Paper_DS/SLL_survey_part2.DS_family_only.csv", p2_dir))



cff <- SLL2.meta[ SLL2.meta$CF_family=="Yes",]
CF.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Cystic_fibrosis=="Yes"])
# order the samples to have DS first, then their families after
cff <- cff[c(rownames(cff)[ rownames(cff) %in% CF.samps], rownames(cff)[ ! rownames(cff) %in% CF.samps]), ]
write.csv(cff, file = sprintf("%s/Paper_CF/SLL_survey_part2.CF_family_only.csv", p2_dir))






gen.pMeans.ds.fam <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"DS_family"]
  }
  )))
)
phy.pMeans.ds.fam <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"DS_family"]
  }
  )))
)
conts.pMeans.ds.fam <- sapply(rownames(famSubsTests$Downs_Syndrome$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"DS_family"]
  }
  )))
)


gen.num_sig.ds.fam <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"DS_family"]
  }
  )) < 0.05)
)
phy.num_sig.ds.fam <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"DS_family"]
  }
  )) < 0.05)
)
conts.num_sig.ds.fam <- sapply(rownames(famSubsTests$Downs_Syndrome$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"DS_family"]
  }
  )) < 0.05)
)




adonis.pMeans.ds.fam <- sapply(names(famSubsTests$Downs_Syndrome$`1`$Adonis), function(x)
  mean(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Adonis[[ x ]][ "DS_family", "Pr(>F)"]
  }
  )))
)

adonis.num_sig.ds.fam <- sapply(names(famSubsTests$Downs_Syndrome$`1`$Adonis), function(x)
  sum(unlist(lapply(famSubsTests$Downs_Syndrome, function(y) {
    y$Adonis[[ x ]][ "DS_family", "Pr(>F)"]
  }
  )) < 0.05)
)



# ****************************************************************************************************************** #



gen.pMeans.cf.fam <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"CF_family"]
  }
  )))
)
phy.pMeans.cf.fam <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"CF_family"]
  }
  )))
)
conts.pMeans.cf.fam <- sapply(rownames(famSubsTests$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"CF_family"]
  }
  )))
)


gen.num_sig.cf.fam <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"CF_family"]
  }
  )) < 0.05)
)
phy.num_sig.cf.fam <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"CF_family"]
  }
  )) < 0.05)
)
conts.num_sig.cf.fam <- sapply(rownames(famSubsTests$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"CF_family"]
  }
  )) < 0.05)
)




adonis.pMeans.cf.fam <- sapply(names(famSubsTests$Cystic_fibrosis$`1`$Adonis), function(x)
  mean(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Adonis[[ x ]][ "CF_family", "Pr(>F)"]
  }
  )))
)

adonis.num_sig.cf.fam <- sapply(names(famSubsTests$Cystic_fibrosis$`1`$Adonis), function(x)
  sum(unlist(lapply(famSubsTests$Cystic_fibrosis, function(y) {
    y$Adonis[[ x ]][ "CF_family", "Pr(>F)"]
  }
  )) < 0.05)
)



# ****************************************************************************************************************** #


group_vs_cont_box(SLL2.meta, unique(unlist(matchCont)), tlev, gloms_clr[[ tlev ]], 
                  "DS_family", "pH", c(groupQs,"DS_family"), only_cont, print_tukey = F, plot_tukey = T)

group_vs_cont_box(SLL2.meta, unique(unlist(matchCont)), tlev, gloms_clr[[ tlev ]], 
                  "DS_family", "Div.Shannon", c(groupQs,"DS_family"), only_cont, print_tukey = F, plot_tukey = T)



group_vs_cont_box(SLL2.meta, unique(unlist(matchCont)), tlev, gloms_clr[[ tlev ]], 
                  "DS_family", "Kingella", c(groupQs,"DS_family"), only_cont, print_tukey = F, plot_tukey = T)



# ****************************************************************************************************************** #






new.dv <- c("MALDI.Yeast_detected","MALDI.Mold_detected","MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies",
            sprintf("Full_MALDI.%s", yeast_of_interest),
            "Allergy","pH","BMI")


# more_anovas.fam <- list()
# 
# for ( disorder in names(famSubsTests)) {
#   
#   # have to make a variable for being in a DS/CF family
#   SLL2@sam_data[ , sprintf("%s_family", ifelse(disorder=="Downs_Syndrome", "DS", "CF")) ] <- "No"
#   SLL2@sam_data[ disFams, sprintf("%s_family", ifelse(disorder=="Downs_Syndrome", "DS", "CF")) ] <- "Yes"
#   SLL2.meta <- meta(SLL2)
#   
#   # ******************************** #
#   if (disorder == "Cystic_fibrosis") {
#     cov_vars <- c("Antibiotics","Gender","Age","Population")#
#     dis.famName <- "CF_family"
#   } else {
#     cov_vars <- c("Gender","Age","Population")#
#     dis.famName <- "DS_family"
#   }
#   # ******************************** #
# 
#   more_anovas.fam[[ disorder ]] <- list()
# 
#   for (i in names(famSubsTests[[ disorder ]])) {
#     print(sprintf("%s  %s", disorder, i))
# 
#     more_anovas.fam[[ disorder ]][[ i ]] <- list()
# 
#     samps.sub <- famSubsTests[[ disorder ]][[ i ]]$samples
#     mTab.sub <- SLL2.meta[ samps.sub, ]
# 
#     for (tl in c("contVar")) {
# 
#       # ******************** #
#       if (tl == "contVar") {
#         glomTab <- NULL
#         dv <- new.dv
#       } else {
#         glomTab <- gloms_clr[[ tl ]][ , samps.sub]
#         dv <- NULL
#       }
#       # ******************** #
# 
#       more_anovas.fam[[ disorder ]][[ i ]][[ tl ]] <- get_lm( c(dis.famName, cov_vars), tl, mTab.sub, glomTab, dv,
#                                                           noRemove = T, rerun.nonSig = F )
# 
#       # ******************** #
#     }
#   }
# }
# saveRDS(more_anovas.fam, sprintf("%s/R_objects/more_anovas.fam.rds", p2_dir))
# # ******************************************************************** #
more_anovas.fam <- readRDS(sprintf("%s/R_objects/more_anovas.fam.rds", p2_dir))
# ******************************************************************** #



more_anovas.pMeans.cf.fam <- sapply(new.dv, function(x)
  mean(unlist(lapply(more_anovas.fam$Cystic_fibrosis, function(y) {
    y$contVar[x,"CF_family"]
  }
  )))
)

more_anovas.pMeans.ds.fam <- sapply(new.dv, function(x)
  mean(unlist(lapply(more_anovas.fam$Downs_Syndrome, function(y) {
    y$contVar[x,"DS_family"]
  }
  )))
)



more_anovas.num_sig.cf.fam <- sapply(new.dv, function(x)
  sum(unlist(lapply(more_anovas.fam$Cystic_fibrosis, function(y) {
    y$contVar[x,"CF_family"]
  }
  )) < 0.05)
)

more_anovas.num_sig.ds.fam <- sapply(new.dv, function(x)
  sum(unlist(lapply(more_anovas.fam$Downs_Syndrome, function(y) {
    y$contVar[x,"DS_family"]
  }
  )) < 0.05)
)
# ******************************************************************** #


















# ****************************************************************************************************************** ####
# Plots for subsamplings ####
# 
# 
# ds.conts <- disSubsTests$default$Downs_Syndrome$`1`$samples
# mTab.sub <- disSubsTests$default$Downs_Syndrome$`1`$mTab
# DS.samps <- ds.conts[ mTab.sub[ds.conts, "Downs_Syndrome"]=="Yes" ]
# cont.samps <- ds.conts[ mTab.sub[ds.conts, "Downs_Syndrome"]=="No" ]
# 
# 
# group_vs_cont_box(mTab.sub, c(DS.samps, cont.samps), tlev, glomTab[[ tlev ]], 
#                   "Downs_Syndrome", "Div.Shannon", groupQs, only_cont, print_tukey = F, plot_tukey = T,
#                   save.boxes = F, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
#                   anova.pTab = disSubsTests$default$Downs_Syndrome$`1`$Anova$contVar )


# ************************************************************************************** #

subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","Div.Shannon", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box")

subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","pH", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "box")
subsampling.plots.box("Downs_Syndrome","contVar","Age","pH", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "scatter", facetVar = "Downs_Syndrome")

subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","BMI", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "box")
subsampling.plots.box("Downs_Syndrome","contVar","Age","BMI", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "scatter", facetVar = "Downs_Syndrome")



# subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","Full_MALDI.Candida_parapsilosis", groupQs, only_cont,
#                       gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
#                       plotType = "gradient", dist_meas = "Aitchison", swapShapeColor=T)
subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","Full_MALDI.Candida_parapsilosis", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "assoc")

dsc <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
dsc.lm <- get_lm( c("Downs_Syndrome", "Gender","Age","Population"), "contVar", SLL2.meta[ dsc, ], 
                  NULL, new.dv, noRemove = T, rerun.nonSig = F )
dsc.lm[ ! is.na(dsc.lm) & dsc.lm > 0.05 ] <- ""
dsc.lm


subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","Full_MALDI.Candida_dubliniensis", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "assoc")


# subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","MALDI.Mold_detected", groupQs, only_cont,
#                       gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
#                       plotType = "assoc")


subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome","MALDI.Yeast_detected", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "assoc")



subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome",
                      c("MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies"), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "box")




# ************************************************************************************** #
# Downs_Syndrome ####
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ds <- gen.num_sig.ds.default[ gen.num_sig.ds.default >= 70]
sort(gen.pMeans.ds.default[ names(freqSig.ds) ])

gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
                 "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
                 "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ds)))
# as.data.frame(cbind(gen.num_sig.ds.default[ gensToCheck ], formatC(gen.pMeans.ds.default[ gensToCheck ],  format="f", digits=5)))

freq.mean.ps <- as.data.frame(cbind(gen.num_sig.ds.default[ gensToCheck ], gen.pMeans.ds.default[ gensToCheck ]))
colnames(freq.mean.ps) <- c("num_sig","meanP")
freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig)), ]
freq.mean.ps


subsampling.plots.box("Downs_Syndrome","Genus","Downs_Syndrome", names(rev(sort(freqSig.ds))), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)

subsampling.plots.box("Downs_Syndrome","Genus","Downs_Syndrome", 
                      c(names(rev(sort(freqSig.ds))),"Actinobacillus"),#"Porphyromonas","Streptococcus"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)


subsampling.plots.box("Downs_Syndrome","Phylum","Downs_Syndrome", 
                      c("Patescibacteria","Bacteroidetes","unclassified.P1","Proteobacteria"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)
subsampling.plots.box("Downs_Syndrome","Phylum","Age", 
                      c("Patescibacteria","Bacteroidetes","unclassified.P1","Proteobacteria"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "scatter", plot_tukey = F)


subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome", 
                      c("Div.Shannon","Species_Richness","Div.Simpson","Faiths.PD"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)


subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome", 
                      c("Stomatotype_Aitchison"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "assoc", plot_tukey = F)



subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome",c("pH","BMI"), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "box")

subsampling.plots.box("Downs_Syndrome","contVar","Downs_Syndrome",c("pH","Div.Shannon"), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="combo", additional.dst=more_anovas,
                      plotType = "box")






subsampling.plots.box("Downs_Syndrome","Genus","Downs_Syndrome",
                      c("Kingella","Alloprevotella","pH","Div.Shannon"), 
                      groupQs, only_cont, gloms_clr, SLL2.meta, SLL2, disSubsTests, 
                      dstStruc="combo", additional.dst=more_anovas, plotType = "box")




# combine plots for paper Figure 1:
library(ggpubr)

DS.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Downs_Syndrome=="Yes"])
mC.DS <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
mC.DS <- mC.DS[ SLL2.meta[mC.DS, "Downs_Syndrome"]=="No" ]

psb.DS <- plot_stacked_bars.disorders("Downs_Syndrome", DS.samps, mC.DS)
spb.DS <- subsampling.plots.box("Downs_Syndrome","Genus","Downs_Syndrome",
                                c("Kingella","Alloprevotella","pH","Div.Shannon"), 
                                groupQs, only_cont, gloms_clr, SLL2.meta, SLL2, disSubsTests, 
                                dstStruc="combo", additional.dst=more_anovas, plotType = "box")

ggarrange(psb.DS, spb.DS, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))








subsampling.plots.box("Downs_Syndrome","contVar","Age","pH", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "scatter", facetVar = "Downs_Syndrome")
subsampling.plots.box("Downs_Syndrome","contVar","Age","BMI", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "scatter")



table(sort( (unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) 
  x$Anova$Genus["Actinobacillus","Downs_Syndrome"])) )) < 0.05)
table(sort(p.adjust(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) 
  x$Anova$Genus["Actinobacillus","Downs_Syndrome"])), method="BH")) < 0.05)

mean(sort( (unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) 
  x$Anova$Genus["Actinobacillus","Downs_Syndrome"])) )) )
mean(sort(p.adjust(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) 
  x$Anova$Genus["Actinobacillus","Downs_Syndrome"])), method="BH")) )









# ******************************************************* #
gen.pMeans.water.ds <- sapply(cont_water_data, function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )))
  ))
gen.pMeans.water.ds.cov <- sapply(cont_water_data, function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )))
  ))


gen.num_sig.water.ds <- sapply(cont_water_data, function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )) < 0.05)
  ))
gen.num_sig.water.ds.cov <- sapply(cont_water_data, function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )) < 0.05)
  ))

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.water.ds <- gen.num_sig.water.ds[ gen.num_sig.water.ds > 90]
freqSig.water.ds <- freqSig.water.ds[ ! is.na(freqSig.water.ds) ]


# ******************************************************* #
gen.pMeans.more.ds <- sapply(rownames(more_anovas$Downs_Syndrome$`1`$contVar), function(x)
  mean(unlist(lapply(more_anovas$Downs_Syndrome, function(y) {
    y$contVar[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )))
)

gen.num_sig.more.ds <- sapply(rownames(more_anovas$Downs_Syndrome$`1`$contVar), function(x)
  sum(unlist(lapply(more_anovas$Downs_Syndrome, function(y) {
    y$contVar[gsub("-","\\.",x),"Downs_Syndrome"]
  }
  )) < 0.05, na.rm = T)
)

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.more.ds <- gen.num_sig.more.ds[ gen.num_sig.more.ds > 90]
freqSig.more.ds <- freqSig.more.ds[ ! is.na(freqSig.more.ds) ]

# ******************************************************* #






# ******************************************************* #
# water vals
conts.pMeans.water.ds <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )))
  ))
conts.pMeans.water.ds.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )))
  ))


conts.num_sig.water.ds <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )) < 0.05)
  ))
conts.num_sig.water.ds.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )) < 0.05)
  ))

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.water.ds <- conts.num_sig.water.ds[ conts.num_sig.water.ds > 90]
freqSig.water.ds <- freqSig.water.ds[ ! is.na(freqSig.water.ds) ]



# **************************** #

gen.num_sig.water.ds <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$Genus$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )) < 0.05)
  ))
gen.num_sig.water.ds.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$Genus$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )) < 0.05)
  ))

gen.pMeans.water.ds <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$Genus$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )))
  ))
gen.pMeans.water.ds.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Downs_Syndrome$`1`$Genus$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )))
  ))


siggen.num_sig.water.ds <- rownames(gen.num_sig.water.ds > 90)[rowSums(gen.num_sig.water.ds > 90) > 0 ]
gen.num_sig.water.ds[ siggen.num_sig.water.ds, ]
gen.pw.ds.sigs <- gen.pMeans.water.ds[ siggen.num_sig.water.ds, ]
gen.pw.ds.sigs <- signif(gen.pw.ds.sigs, digits=5)
gen.pw.ds.sigs[ gen.pw.ds.sigs > 0.05 ] <- ""
gen.pw.ds.sigs

# siggen.num_sig.water.ds.cov <- rownames(gen.num_sig.water.ds.cov > 90)[rowSums(gen.num_sig.water.ds.cov > 90) > 0 ]
gen.num_sig.water.ds.cov[ siggen.num_sig.water.ds, ]
gen.pw.ds.cov.sigs <- gen.pMeans.water.ds.cov[ siggen.num_sig.water.ds, ]
gen.pw.ds.cov.sigs <- signif(gen.pw.ds.cov.sigs, digits=5)
gen.pw.ds.cov.sigs[ gen.pw.ds.cov.sigs > 0.05 ] <- ""
gen.pw.ds.cov.sigs

# those water vals that were sig dif, but also sig by Downs_Syndrome
gpwdss.toCheck <- gen.pw.ds.sigs
gpwdss.toCheck[ ! (gen.pw.ds.sigs != "" & gen.pw.ds.cov.sigs != "") ] <- ""
gpwdss.toCheck

gpwdss.cov.toCheck <- gen.pw.ds.cov.sigs
gpwdss.cov.toCheck[ ! (gen.pw.ds.sigs != "" & gen.pw.ds.cov.sigs != "") ] <- ""
gpwdss.cov.toCheck


# subsampling.plots.box("Downs_Syndrome","Genus","Conductivity","Curvibacter", groupQs, only_cont,
#                       gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="water", additional.dst=water_anovas,
#                       plotType = "scatter", facetVar = "Downs_Syndrome")

for (gen in rownames(gpwdss.toCheck)) {
  for (wat in colnames(gpwdss.toCheck)) {
    
    if (gpwdss.toCheck[ gen, wat ] != "") {
      
      subsampling.plots.box("Downs_Syndrome","Genus", wat, gen, groupQs, only_cont,
                            gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="water", additional.dst=water_anovas,
                            plotType = "scatter", facetVar = "Downs_Syndrome",
                            save_plot = T, plot_dir = "Paper.Downs_Syndrome_plots/water_vals")
      
    }
  }
}





# ************************************************************************************** #



# # trying to get average values of all pvalues in anova tests, then plot this with all the controls used
# ds.conts.all <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
# mTab.sub.all <- SLL2.meta[ds.conts.all, ]
# DS.samps <- ds.conts.all[ mTab.sub.all[ds.conts.all, "Downs_Syndrome"]=="Yes" ]
# cont.samps.all <- ds.conts.all[ mTab.sub.all[ds.conts.all, "Downs_Syndrome"]=="No" ]
# 
# # here are the mean pvalues for Downs_Syndrome and the other 3 covars 
# #  problem is that when using na.rm=T maybe only like 3 instances in which it worked
# mps <- as.data.frame(sapply(colnames(disSubsTests$default$Downs_Syndrome$`1`$Anova$contVar), function(y) 
#   mean(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) {
#     val <- x$Anova$contVar["Div.Shannon", y]
#     ifelse(is.na(val), 1, val)
#     })) )))#, 
#     # na.rm=F)))
# colnames(mps) <- "Div.Shannon"
# mps <- as.matrix(t(mps))
# 
# 
# # only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
# gen.num_sig.ds.default[ gen.num_sig.ds.default > 95]
# 
# group_vs_cont_box(mTab.sub.all, c(DS.samps, cont.samps.all), tlev, glomTab[[ tlev ]], 
#                   "Downs_Syndrome", "Div.Shannon", groupQs, only_cont, print_tukey = F, plot_tukey = T)
#                   # save.boxes = F, plot_dirs = sprintf("LMM_signifs/%s/%s", nor, f),
#                   # anova.pTab = disSubsTests$default$Downs_Syndrome$`1`$Anova$contVar )

# ****************************************************************************************************************** #



# ************************************************************************************** #
# Cystic_fibrosis ####
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.cf <- gen.num_sig.cf.default[ gen.num_sig.cf.default >= 80]
sort(gen.pMeans.cf.default[ names(freqSig.cf) ])

gensToCheck <- unique(c("Porphyromonas","Treponema","Tannerella",
                        "Fusobacterium","Prevotella","Parvimonas","Campylobacter",
                        "Eikenella","Aggregatibacter",
                        "Haemophilus","Staphylococcus",
                        "Pseudomonas","Rothia","Brevundimonas","Streptococcus",
                        # "Burkholderia","Achromobacter",
                        names(freqSig.cf)))
# as.data.frame(cbind(gen.num_sig.cf.default[ gensToCheck ], formatC(gen.pMeans.cf.default[ gensToCheck ],  format="f", digits=5)))

freq.mean.ps.cf <- as.data.frame(cbind(gen.num_sig.cf.default[ gensToCheck ], gen.pMeans.cf.default[ gensToCheck ]))
colnames(freq.mean.ps.cf) <- c("num_sig","meanP")
# freq.mean.ps.cf <- freq.mean.ps.cf[ rev(order(freq.mean.ps.cf$num_sig, freq.mean.ps.cf$meanP, decreasing = c(F, T))), ]
freq.mean.ps.cf <- freq.mean.ps.cf[ order(freq.mean.ps.cf$meanP), ]
freq.mean.ps.cf





# combine plots for paper Figure 1:
library(ggpubr)

CF.samps <- sort(rownames(SLL2.meta)[SLL2.meta$Cystic_fibrosis=="Yes"])
mC.CF <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
mC.CF <- mC.CF[ SLL2.meta[mC.CF, "Cystic_fibrosis"]=="No" ]

psb.CF <- plot_stacked_bars.disorders("Cystic_fibrosis", CF.samps, mC.CF, ang = 15)
spb.CF <- subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis",
                                c("Rothia","Treponema","Div.Shannon","pH"),
                                # c("Rothia","Aggregatibacter","pH","Div.Shannon"), 
                                groupQs, only_cont, gloms_clr, SLL2.meta, SLL2, disSubsTests, 
                                dstStruc="combo", additional.dst=more_anovas, plotType = "box",
                                atex.s=25, atit.s=27, pts=23, sts=22)

ggarrange(psb.CF, spb.CF, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),
          font.label = list(size=19))








subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis",
                      c("Rothia","Treponema","Div.Shannon","pH"), 
                      groupQs, only_cont, gloms_clr, SLL2.meta, SLL2, disSubsTests, 
                      dstStruc="combo", additional.dst=more_anovas, plotType = "box")






# genera associated with lung infections
subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis", 
                      c("Pseudomonas","Rothia","Brevundimonas","Haemophilus","Staphylococcus","Streptococcus"),
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)


# genera associated with periodontitis
subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis", 
                      c("Porphyromonas","Treponema","Tannerella",
                        "Fusobacterium","Prevotella","Parvimonas","Campylobacter",
                        "Eikenella","Aggregatibacter",
                        "Peptostreptococcus","Bergeyella"),
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)



# check particular species within some of these to see if they correspond to correct ones
subsampling.plots.box("Cystic_fibrosis","Species","Cystic_fibrosis", 
                      c("Eikenella corrodens","Porphyromonas gingivalis",
                        "Fusobacterium nucleatum",
                        "Aggregatibacter aphrophilus",#"Aggregatibacter actinomycetemcomitans",
                        "Tannerella forsythia"),
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)

# check all Prevotella species
prev.specs <- rownames(taxa_print("Genus","Prevotella"))[ rownames(taxa_print("Genus","Prevotella")) %in% rownames(gloms_clr$Species) ]
subsampling.plots.box("Cystic_fibrosis","Species","Cystic_fibrosis", 
                      prev.specs,
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)

# check all Treponema species
trep.specs <- rownames(taxa_print("Genus","Treponema"))[ rownames(taxa_print("Genus","Treponema")) %in% rownames(gloms_clr$Species) ]
subsampling.plots.box("Cystic_fibrosis","Species","Cystic_fibrosis", 
                      trep.specs,
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)




# just those that were signif in at least 90 subsamplings
subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis", 
                      c(names(rev(sort(freqSig.cf))), "Fusobacterium", "Haemophilus"),
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F,
                      atex.s=20, atit.s=27, pts=25, sts=16)



# phyla that were signif
sort(phy.num_sig.cf.default)
subsampling.plots.box("Cystic_fibrosis","Phylum","Cystic_fibrosis", 
                      c("Firmicutes","unclassified.P1","Actinobacteria","Spirochaetes","Fusobacteria","Patescibacteria"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F,
                      atex.s=25, atit.s=27, pts=23, sts=20)

# a div values
sort(conts.num_sig.cf.default)
subsampling.plots.box("Cystic_fibrosis","contVar","Cystic_fibrosis", 
                      c("Div.Shannon","Faiths.PD","Species_Richness","Div.Simpson"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)







subsampling.plots.box("Cystic_fibrosis","Genus","Cystic_fibrosis", 
                      c("Pseudomonas","Burkholderia","Achromobacter"),
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, plotType = "box", plot_tukey = F)



# ******************************************************* #
# beta div values
vcd::assoc(table(disSubsTests$default$Cystic_fibrosis$`1`$mTab$Cystic_fibrosis, 
                 disSubsTests$default$Cystic_fibrosis$`1`$mTab$Stomatotype_Unweighted_Unifrac), 
           shade=T, main="Stomatotype_Unweighted_Unifrac (i=1)")
table(disSubsTests$default$Cystic_fibrosis$`1`$mTab$Cystic_fibrosis, 
      disSubsTests$default$Cystic_fibrosis$`1`$mTab$Stomatotype_Unweighted_Unifrac,
      dnn = c("Cystic_Fibrosis","Stomatotype_Unweighted_Unifrac"))




vcd::assoc(table(disSubsTests$default$Cystic_fibrosis$`2`$mTab$Cystic_fibrosis, 
                 disSubsTests$default$Cystic_fibrosis$`2`$mTab$Stomatotype_Jaccard), 
           shade=T, main="Stomatotype_Jaccard (i=2)")
table(disSubsTests$default$Cystic_fibrosis$`2`$mTab$Cystic_fibrosis, 
      disSubsTests$default$Cystic_fibrosis$`2`$mTab$Stomatotype_Jaccard,
      dnn = c("Cystic_Fibrosis","Stomatotype_Jaccard"))





vcd::assoc(table(disSubsTests$default$Cystic_fibrosis$`1`$mTab$Cystic_fibrosis, 
                 disSubsTests$default$Cystic_fibrosis$`1`$mTab$Stomatotype_Aitchison), 
           shade=T, main="Stomatotype_Aitchison (i=1)")
table(disSubsTests$default$Cystic_fibrosis$`1`$mTab$Cystic_fibrosis, 
      disSubsTests$default$Cystic_fibrosis$`1`$mTab$Stomatotype_Aitchison,
      dnn = c("Cystic_Fibrosis","Stomatotype_Aitchison"))

vcd::assoc(table(disSubsTests$default$Cystic_fibrosis$`2`$mTab$Cystic_fibrosis, 
                 disSubsTests$default$Cystic_fibrosis$`2`$mTab$Stomatotype_Aitchison), 
           shade=T, main="Stomatotype_Aitchison (i=2)")
table(disSubsTests$default$Cystic_fibrosis$`2`$mTab$Cystic_fibrosis, 
      disSubsTests$default$Cystic_fibrosis$`2`$mTab$Stomatotype_Aitchison,
      dnn = c("Cystic_Fibrosis","Stomatotype_Aitchison"))



# ******************************************************* #
gen.pMeans.more.cf <- sapply(rownames(more_anovas$Cystic_fibrosis$`1`$contVar), function(x)
  mean(unlist(lapply(more_anovas$Cystic_fibrosis, function(y) {
    y$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)

gen.num_sig.more.cf <- sapply(rownames(more_anovas$Cystic_fibrosis$`1`$contVar), function(x)
  sum(unlist(lapply(more_anovas$Cystic_fibrosis, function(y) {
    y$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05, na.rm = T)
)

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.more.cf <- gen.num_sig.more.cf[ gen.num_sig.more.cf > 90]
freqSig.more.cf <- freqSig.more.cf[ ! is.na(freqSig.more.cf) ]

subsampling.plots.box("Cystic_fibrosis","contVar","Cystic_fibrosis",c("pH","BMI"), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "box")

# subsampling.plots.box("Cystic_fibrosis","contVar","Age","pH", groupQs, only_cont,
#                       gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
#                       plotType = "scatter", facetVar = "Cystic_fibrosis")
subsampling.plots.box("Cystic_fibrosis","contVar","Age","BMI", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "scatter")


subsampling.plots.box("Cystic_fibrosis","contVar","Cystic_fibrosis","Full_MALDI.Candida_albicans", groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="additional", additional.dst=more_anovas,
                      plotType = "assoc")





# ******************************************************* #
# water vals
conts.pMeans.water.cf <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )))
  ))
conts.pMeans.water.cf.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$contVar$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Cystic_fibrosis"]
    }
    )))
  ))


conts.num_sig.water.cf <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )) < 0.05)
  ))
conts.num_sig.water.cf.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$contVar$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$contVar[[ cwd ]][gsub("-","\\.",x), "Cystic_fibrosis"]
    }
    )) < 0.05)
  ))

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.water.cf <- conts.num_sig.water.cf[ conts.num_sig.water.cf > 90]
freqSig.water.cf <- freqSig.water.cf[ ! is.na(freqSig.water.cf) ]



# **************************** #

gen.num_sig.water.cf <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$Genus$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )) < 0.05)
  ))
gen.num_sig.water.cf.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$Genus$Conductivity), function(x)
    sum(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), "Cystic_fibrosis"]
    }
    )) < 0.05)
  ))

gen.pMeans.water.cf <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$Genus$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), cwd]
    }
    )))
  ))
gen.pMeans.water.cf.cov <- sapply(c(cont_water_data, "Age"), function(cwd)
  sapply(rownames(water_anovas$Cystic_fibrosis$`1`$Genus$Conductivity), function(x)
    mean(unlist(lapply(water_anovas$Cystic_fibrosis, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Genus[[ cwd ]][gsub("-","\\.",x), "Cystic_fibrosis"]
    }
    )))
  ))


siggen.num_sig.water.cf <- rownames(gen.num_sig.water.cf > 90)[rowSums(gen.num_sig.water.cf > 90) > 0 ]
gen.num_sig.water.cf[ siggen.num_sig.water.cf, ]
gen.pw.cf.sigs <- gen.pMeans.water.cf[ siggen.num_sig.water.cf, ]
gen.pw.cf.sigs <- signif(gen.pw.cf.sigs, digits=5)
gen.pw.cf.sigs[ gen.pw.cf.sigs > 0.05 ] <- ""
gen.pw.cf.sigs

# siggen.num_sig.water.cf.cov <- rownames(gen.num_sig.water.cf.cov > 90)[rowSums(gen.num_sig.water.cf.cov > 90) > 0 ]
gen.num_sig.water.cf.cov[ siggen.num_sig.water.cf, ]
gen.pw.cf.cov.sigs <- gen.pMeans.water.cf.cov[ siggen.num_sig.water.cf, ]
gen.pw.cf.cov.sigs <- signif(gen.pw.cf.cov.sigs, digits=5)
gen.pw.cf.cov.sigs[ gen.pw.cf.cov.sigs > 0.05 ] <- ""
gen.pw.cf.cov.sigs

# those water vals that were sig dif, but also sig by Cystic_fibrosis
gpwcfs.toCheck <- gen.pw.cf.sigs
gpwcfs.toCheck[ ! (gen.pw.cf.sigs != "" & gen.pw.cf.cov.sigs != "") ] <- ""
gpwcfs.toCheck

gpwcfs.cov.toCheck <- gen.pw.cf.cov.sigs
gpwcfs.cov.toCheck[ ! (gen.pw.cf.sigs != "" & gen.pw.cf.cov.sigs != "") ] <- ""
gpwcfs.cov.toCheck


# subsampling.plots.box("Cystic_fibrosis","Genus","Conductivity","Curvibacter", groupQs, only_cont,
#                       gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="water", additional.dst=water_anovas,
#                       plotType = "scatter", facetVar = "Cystic_fibrosis")

for (gen in rownames(gpwcfs.toCheck)) {
  for (wat in colnames(gpwcfs.toCheck)) {
    
    if (gpwcfs.toCheck[ gen, wat ] != "") {
      
      subsampling.plots.box("Cystic_fibrosis","Genus", wat, gen, groupQs, only_cont,
                            gloms_clr, SLL2.meta, SLL2, disSubsTests, dstStruc="water", additional.dst=water_anovas,
                            plotType = "scatter", facetVar = "Cystic_fibrosis",
                            save_plot = T, plot_dir = "Paper.Cystic_fibrosis_plots/water_vals")
      
    }
  }
}

# ****************************************************************************************************************** #









# ****************************************************************************************************************** ####
# Test family similarities in disorder subsamplings ####


# *********************************************************************************** #
get_anosim.disFam_only <- function(disorder, mTab, phy) {
  
  # disorder samples
  disSamps <- rownames( mTab[ mTab[ , disorder ] == "Yes", ] )
  
  # get family member samples 
  dis.famUnits <- unique( mTab[disSamps, "Family_unit"] )
  dis.famUnits <- dis.famUnits[ dis.famUnits != "None" ]
  disFams <- rownames( mTab[ mTab[, "Family_unit"] %in% dis.famUnits, ])
  
  # objects with disorders and their family members only
  dis_with_Fam <- unique(c(disSamps, disFams))
  phy.Fams <- prune_samples( dis_with_Fam, phy)
  mTab.Fams <- mTab[ dis_with_Fam, ]
  
  ordObj <- subsampling_ordination_objects(dis_with_Fam, phy.Fams, distsOnly = T)
  
  
  anosim_list <- list()
  for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                      "Bray","Jaccard")) {
    
    anosim_list[[ dist_meas ]] <- list()
    
    for (unit_lab in c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
                       "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
                       "Grandparent_Grandchild_unit", disorder)) {
      if (disorder == "Downs_Syndrome" & unit_lab %in% c("Twin_unit","Partner_unit","Grandparent_Grandchild_unit"))
        next
      
      print(c(dist_meas, unit_lab))
      
      # # make variable to show unit+disorder
      # mTab.Fams[ , sprintf("%s_%s", unit_lab, disorder) ] <- sapply(rownames(mTab.Fams), function(x) 
      #   sprintf("%s-%s", mTab.Fams[x, unit_lab], mTab.Fams[x, disorder]))
      
      
      # anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( ordObj[[ dist_meas ]] ),
      #                                                    as.character(as.matrix(mTab.Fams[, sprintf("%s_%s", unit_lab, disorder)])))#,
      # strata = SLL2.meta[ , "seqGroup"])
      
      anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( ordObj[[ dist_meas ]] ),
                                                         as.character(as.matrix(mTab.Fams[, unit_lab])))#,
      
      # mTab <- SLL2.meta[ SLL2.meta[,unit_lab] != "None", ]
      # anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas))[rownames(mTab),rownames(mTab)] ), 
      #                                                    as.character(as.matrix(mTab[, unit_lab])))
    }
  }
  return(anosim_list)
}
# *********************************************************************************** #


# anosim_list <- get_anosim.disorders("Cystic_fibrosis", 1, disSubsTests, SLL2.meta, SLL2)
anosim.disFams <- list()
dis_with_Fam <- list()
phy.Fams <- list()
mTab.Fams <- list()
ordObj <- list()

for (disorder in c("Cystic_fibrosis","Downs_Syndrome")) {
  
  anosim.disFams[[ disorder ]] <- get_anosim.disFam_only(disorder, SLL2.meta, SLL2)
  
  # get disorder samples
  disSamps <- rownames( SLL2.meta[ SLL2.meta[ , disorder ] == "Yes", ] )
  
  # get family member samples 
  dis.famUnits <- unique( SLL2.meta[disSamps, "Family_unit"] )
  dis.famUnits <- dis.famUnits[ dis.famUnits != "None" ]
  disFams <- rownames( SLL2.meta[ SLL2.meta[, "Family_unit"] %in% dis.famUnits, ])
  
  # objects with disorders and their family members only
  dis_with_Fam[[ disorder ]] <- unique(c(disSamps, disFams))
  phy.Fams[[ disorder ]] <- prune_samples( dis_with_Fam[[ disorder ]], SLL2)
  mTab.Fams[[ disorder ]] <- SLL2.meta[ dis_with_Fam[[ disorder ]], ]
  
  ordObj[[ disorder ]] <- subsampling_ordination_objects(dis_with_Fam[[ disorder ]], phy.Fams[[ disorder ]], distsOnly = T)
}


# plot differences in average distances for group members vs non group members


# ****************************** #
plot_famUnit_diffs.disorders <- function(sam_tab, dist_obj, dist_name, comp.method, disorder, 
                                         fam.anosim=NULL, fam.adonis=NULL) {
  
  fam_units <- c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
                 "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
                 "Grandparent_Grandchild_unit")
  if (disorder == "Downs_Syndrome")
    fam_units <- c("Family_unit","Sibling_unit","Parent_Child_unit")
  
  all_unit_labs <- list()
  unit_labs_pvals <- list()
  
  # ****************************** #
  for (unit_lab in fam_units) {
    
    # get vector of units of particular type
    units <- unique( sam_tab[ , unit_lab] )
    units <- units[ units != "None" ]
    
    # get vector of dists for each unit - members of unit against same unit
    all_unit_labs[[ sprintf("%s.within", unit_lab) ]] <- unlist(sapply(units, function(x) {
      unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
      # make dist_obj as.dist() first, 
      #   in this case it removes NAs and duplicated values
      dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
    }))
    
    # get vector of dists for each unit - members of unit against those not in that unit 
    #  (or against only other units??)
    all_unit_labs[[ sprintf("%s.between", unit_lab) ]] <- unlist(sapply(units, function(x) {
      unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] != x, ] )
      # unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
      
      # make dist_obj as.matrix() first,
      #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
      dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
    }))
    
    
    # with.bet.diff <- t.test(within_dist, between_dist)
    # print(with.bet.diff)
    
    if (comp.method == "ANOSIM") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.anosim[[ dist_name ]][[ unit_lab ]]$signif
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- fam.anosim[[ dist_name ]][[ unit_lab ]]$statistic
      
    } else if (comp.method == "PERMANOVA") {
      unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$`Pr(>F)`[ 1 ]
      unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- sqrt(fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$R2[ 1 ])
      unit_labs_pvals[[ sprintf("%s.F", unit_lab) ]]   <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$F.Model[ 1 ]
    }
    
  }
  # ****************************** #
  # print(str(all_unit_labs))
  # print(str(unit_labs_pvals))
  
  
  # ****************************** #
  disSamps <- rownames(sam_tab[ sam_tab[ , disorder ] == "Yes", ] )
  non.disSamps <- rownames(sam_tab[ sam_tab[ , disorder ] == "No", ] )
  
  all_unit_labs[[ sprintf("disorder") ]] <- as.numeric( as.dist( dist_obj[ disSamps, disSamps ] ))
  all_unit_labs[[ sprintf("controls") ]] <- as.numeric( as.dist( dist_obj[ non.disSamps, non.disSamps ] ))
  all_unit_labs[[ sprintf("dis_vs_con") ]] <- as.numeric( as.matrix( dist_obj[ disSamps, non.disSamps ] ))
  
  if (comp.method == "ANOSIM") {
    unit_labs_pvals[[ sprintf("%s.pval", "disorder") ]] <- fam.anosim[[ dist_name ]][[ disorder ]]$signif
    unit_labs_pvals[[ sprintf("%s.R2", "disorder") ]]   <- fam.anosim[[ dist_name ]][[ disorder ]]$statistic
    
  } else if (comp.method == "PERMANOVA") {
    unit_labs_pvals[[ sprintf("%s.pval", "disorder") ]] <- fam.adonis[[ dist_name ]][[ disorder ]]$aov.tab$`Pr(>F)`[ 1 ]
    unit_labs_pvals[[ sprintf("%s.R2", "disorder") ]]   <- sqrt(fam.adonis[[ dist_name ]][[ disorder ]]$aov.tab$R2[ 1 ])
    unit_labs_pvals[[ sprintf("%s.F", "disorder") ]]   <- fam.adonis[[ dist_name ]][[ disorder ]]$aov.tab$F.Model[ 1 ]
  }
  # ****************************** #
  
  
  
  # print(unlist(lapply(all_unit_labs, mean)))
  with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
                         "labs" = c(rep(c("Same Family","Different Families"), length(fam_units)), 
                                    "disorder", "controls", "dis_vs_con"),
                         "sds"  = unlist(lapply(all_unit_labs, sd)),
                         "unit" = c(rep(gsub("_unit","",fam_units), each = 2), rep("disorder", 3)))
  
  with.bet$labs <- factor(with.bet$labs, levels = rev(c("Same Family","Different Families","disorder","dis_vs_con","controls")))
  with.bet$unit <- factor(with.bet$unit, levels = c("disorder",rev(gsub("_unit","",fam_units))))
  
  # with.bet <- data.frame("vals" = c(mean(within_dist), mean(between_dist)),
  #                        "labs" = c("Within","Between"),
  #                        "sds"  = c(sd(within_dist), sd(between_dist)),
  #                        "unit" = c(unit_lab, unit_lab))
  # "labs"=c( rep("Within", length(within_dist)), rep("Between",length(between_dist)) ))
  # print(with.bet)
  
  star_labels <- sapply(c(fam_units,"disorder"), function(x) {
    if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
      sprintf("R=%s **", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
      sprintf("R=%s *", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
    } else {
      # sprintf("R=%s", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
      ""
    }
  })
  # print(star_labels)
  
  ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
    geom_bar(stat="identity", position = "dodge2") +
    guides(fill=guide_legend(title=NULL, reverse = T)) +
    geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
    coord_flip() +
    # geom_text(aes(x=unit, y=vals+0.05))
    annotate("text", 
             x=1:length(c(fam_units,"disorder")), 
             y=with.bet$vals[with.bet$labs %in% c("Different Families","controls")]+0.05, 
             label = rev(star_labels), size=5) +
    theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
          legend.text = element_text(size=15)) +#, legend.position = "bottom") +
    xlab(sprintf("%s in family units", comp.method)) + ylab(dist_name)
}
# ****************************** #


plot_famUnit_diffs.disorders(mTab.Fams[[ disorder ]], as.matrix(ordObj[[ disorder ]]$Aitchison), 
                             "Aitchison", "ANOSIM", disorder, fam.anosim = anosim.disFams[[ disorder ]])

plot_famUnit_diffs.disorders(mTab.Fams[[ disorder ]], as.matrix(ordObj[[ disorder ]]$Weighted_Unifrac), 
                             "Weighted_Unifrac", "ANOSIM", disorder, fam.anosim = anosim.disFams[[ disorder ]])

plot_famUnit_diffs.disorders(mTab.Fams[[ disorder ]], as.matrix(ordObj[[ disorder ]]$Unweighted_Unifrac), 
                             "Unweighted_Unifrac", "ANOSIM", disorder, fam.anosim = anosim.disFams[[ disorder ]])

plot_famUnit_diffs.disorders(mTab.Fams[[ disorder ]], as.matrix(ordObj[[ disorder ]]$Bray), 
                             "Bray", "ANOSIM", disorder, fam.anosim = anosim.disFams[[ disorder ]])

plot_famUnit_diffs.disorders(mTab.Fams[[ disorder ]], as.matrix(ordObj[[ disorder ]]$Jaccard), 
                             "Jaccard", "ANOSIM", disorder, fam.anosim = anosim.disFams[[ disorder ]])
# *********************************************************************************** #












# *********************************************************************************** #
# Anosim of disorders vs controls ####

# to do so, must get the list of disorder samples and controls
# *********************************************************************************** #
get_anosim.disorders <- function(disorder, i, dst, mTab, phy, sampMode="default") {
  
  disSamps <- dst[[ sampMode ]][[ disorder ]][[ i ]]$samples
  
  phy.subs <- prune_samples( disSamps, phy)
  mTab.subs <- mTab[ disSamps, ]
  
  ordObj <- subsampling_ordination_objects(disSamps, phy.subs, distsOnly = T)
  
  
  anosim_list <- list()
  for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                      "Bray","Jaccard")) {
    
    anosim_list[[ dist_meas ]] <- anosim(as.dist( ordObj[[ dist_meas ]] ),
                                         as.character(as.matrix(mTab.subs[, disorder ])))#,
      # strata = SLL2.meta[ , "seqGroup"])
      
      # mTab <- SLL2.meta[ SLL2.meta[,unit_lab] != "None", ]
      # anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas))[rownames(mTab),rownames(mTab)] ), 
      #                                                    as.character(as.matrix(mTab[, unit_lab])))
  }
  return(anosim_list)
}

# *********************************************************************************** #

anosim.disorders <- list()

for (disorder in c("Cystic_fibrosis","Downs_Syndrome")) {
  
  anosim.disorders[[ disorder ]] <- list()
  
  for (i in names(disSubsTests$default[[ disorder ]])) {
    
    anosim.disorders[[ disorder ]][[ i ]] <- get_anosim.disorders(disorder, i, disSubsTests, SLL2.meta, SLL2)
  }
}




# *********************************************************************************** #



# # ****************************** #
# plot_famUnit_diffs.disorders3 <- function(sam_tab, dist_obj, dist_name, comp.method, disorder, fam.anosim=NULL, fam.adonis=NULL) {
#   
#   fam_units <- c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
#                  "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
#                  "Grandparent_Grandchild_unit")
#   all_unit_labs <- list()
#   unit_labs_pvals <- list()
#   
#   for (unit_lab in fam_units) {
#     
#     sam_tab[ , sprintf("%s_%s", unit_lab, disorder) ] <- sapply(rownames(sam_tab), function(x) 
#       sprintf("%s-%s", sam_tab[x, unit_lab], sam_tab[x, disorder]))
#     
#     
#     # get vector of units of particular type
#     units <- unique( sam_tab[ , unit_lab] )
#     units <- units[ units != "None" ]
#     
#     # ******* #
#     # get vector of dists for each unit - members of unit against same unit - with disorder
#     all_unit_labs[[ sprintf("%s.within-dis", unit_lab) ]] <- unlist(sapply(units, function(x) {
#       unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x &
#                                         sam_tab[ , disorder ] == "Yes", ] )
#       # make dist_obj as.dist() first, 
#       #   in this case it removes NAs and duplicated values
#       dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
#     }))
#     # ******* #
#     # get vector of dists for each unit - members of unit against same unit - without disorder
#     all_unit_labs[[ sprintf("%s.within-con", unit_lab) ]] <- unlist(sapply(units, function(x) {
#       unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x &
#                                         sam_tab[ , disorder ] == "No", ] )
#       # make dist_obj as.dist() first, 
#       #   in this case it removes NAs and duplicated values
#       dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
#     }))
#     
#     # ******* #
#     # get vector of dists for each unit - members of unit against those not in that unit - with disorder
#     #  (or against only other units??)
#     all_unit_labs[[ sprintf("%s.between-dis", unit_lab) ]] <- unlist(sapply(units, function(x) {
#       unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] != x &
#                                         sam_tab[ , disorder ] == "Yes", ] )
#       # unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
#       
#       # make dist_obj as.matrix() first,
#       #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
#       dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
#     }))
#     # ******* #
#     # get vector of dists for each unit - members of unit against those not in that unit - without disorder
#     #  (or against only other units??)
#     all_unit_labs[[ sprintf("%s.between-con", unit_lab) ]] <- unlist(sapply(units, function(x) {
#       unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] != x &
#                                         sam_tab[ , disorder ] == "No", ] )
#       # unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
#       
#       # make dist_obj as.matrix() first,
#       #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
#       dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
#     }))
#     # ******* #
#     
#     
#     # with.bet.diff <- t.test(within_dist, between_dist)
#     # print(with.bet.diff)
#     
#     if (comp.method == "ANOSIM") {
#       unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.anosim[[ dist_name ]][[ unit_lab ]]$signif
#       unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- fam.anosim[[ dist_name ]][[ unit_lab ]]$statistic
#       
#     } else if (comp.method == "PERMANOVA") {
#       unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$`Pr(>F)`[ 1 ]
#       unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- sqrt(fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$R2[ 1 ])
#       unit_labs_pvals[[ sprintf("%s.F", unit_lab) ]]   <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$F.Model[ 1 ]
#     }
#     
#   }
#   # print(str(all_unit_labs))
#   # print(str(unit_labs_pvals))
#   
#   # # get vector of units of particular type
#   # units <- unique( sam_tab[ , unit_lab] )
#   # units <- units[ units != "None" ]
#   # 
#   # # get vector of dists for each unit - members of unit against same unit
#   # within_dist <- unlist(sapply(units, function(x) {
#   #   unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
#   #   # make dist_obj as.dist() first, 
#   #   #   in this case it removes NAs and duplicated values
#   #   dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
#   # }))
#   # 
#   # # get vector of dists for each unit - members of unit against those not in that unit 
#   # #  (or against only other units??)
#   # between_dist <- unlist(sapply(units, function(x) {
#   #   unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] != x, ] )
#   #   # make dist_obj as.matrix() first, 
#   #   #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
#   #   dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
#   # }))
#   # # between_dist <- sapply(units, function(x) {
#   # #   unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
#   # #   dists <- as.numeric( dist_obj[ unit_samps, unit_samps ])
#   # # }, sam_tab, unit_lab, dist_obj)
#   # 
#   # with.bet.diff <- t.test(within_dist, between_dist)
#   # print(with.bet.diff)
#   
#   # print(unlist(lapply(all_unit_labs, mean)))
#   with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
#                          # "labs" = rep(c("Same Family","Different Families"), length(fam_units)),
#                          "labs" = rep(c("Same Family - dis","Same Family - con",
#                                         "Different Families - dis","Different Families - con"), length(fam_units)),
#                          "sds"  = unlist(lapply(all_unit_labs, sd)),
#                          "unit" = rep(gsub("_unit","",fam_units), each = 4))
#   
#   with.bet$unit <- factor(with.bet$unit, levels = rev(gsub("_unit","",fam_units)))
#   
#   # with.bet <- data.frame("vals" = c(mean(within_dist), mean(between_dist)),
#   #                        "labs" = c("Within","Between"),
#   #                        "sds"  = c(sd(within_dist), sd(between_dist)),
#   #                        "unit" = c(unit_lab, unit_lab))
#   # "labs"=c( rep("Within", length(within_dist)), rep("Between",length(between_dist)) ))
#   # print(with.bet)
#   
#   star_labels <- sapply(fam_units, function(x) {
#     if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
#       sprintf("R=%s **", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
#     } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
#       sprintf("R=%s *", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
#     } else {
#       # sprintf("R=%s", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       ""
#     }
#   })
#   # print(star_labels)
#   print(with.bet)
#   print(star_labels)
#   ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
#     geom_bar(stat="identity", position = "dodge2") +
#     guides(fill=guide_legend(title=NULL, reverse = T)) +
#     geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
#     coord_flip() +
#     # geom_text(aes(x=unit, y=vals+0.05))
#     annotate("text", x=1:length(fam_units), y=with.bet$vals[with.bet$labs=="Different Families - dis"]+0.05, 
#              label = rev(star_labels), size=5) +
#     theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
#           legend.text = element_text(size=15), legend.position = "bottom") +
#     xlab(sprintf("%s in family units", comp.method)) + ylab(dist_name)
# }
# # ****************************** #
# 
# # ****************************** #
# plot_famUnit_diffs.disorders2 <- function(sam_tab, dist_obj, dist_name, comp.method, disorder, fam.anosim=NULL, fam.adonis=NULL) {
#   
#   fam_units <- c("Family_unit","Sibling_unit","Twin_unit","Partner_unit",
#                  "Parent_Child_unit",#"Mother_Child_unit","Father_Child_unit",
#                  "Grandparent_Grandchild_unit")
#   all_unit_labs <- list()
#   unit_labs_pvals <- list()
#   
#   for (unit_lab in fam_units) {
#     
#     sam_tab[ , sprintf("%s_%s", unit_lab, disorder) ] <- sapply(rownames(sam_tab), function(x) 
#       sprintf("%s-%s", sam_tab[x, unit_lab], sam_tab[x, disorder]))
#     
#     
#     # get vector of units of particular type
#     units <- unique( sam_tab[ , unit_lab] )
#     units <- units[ units != "None" ]
#     
#     # ******* #
#     # get vector of dists for each unit - members of unit against same unit - with disorder
#     all_unit_labs[[ sprintf("%s.within-disOnly", unit_lab) ]] <- unlist(sapply(units, function(x) {
#       unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ])# &
#       # sam_tab[ , disorder ] == "Yes", ] )
#       # print(unit_samps)
#       tempDisObj <- as.matrix(dist_obj)[ unit_samps, unit_samps]
#       print(tempDisObj)
#       # make those values that are between 2 nonDisorder samples NA
#       #   also between same sample
#       print(unit_lab)
#       for (s1 in unit_samps) {
#         for (s2 in unit_samps) {
#           if (s1 == s2 | 
#               sum(sam_tab[ c(s1, s2), disorder] == "No") == 2) {
#             tempDisObj[s1, s2] <- NA
#           }
#         }
#       }
#       print(tempDisObj)
#       return( as.numeric(tempDisObj)[ ! is.na(as.numeric(tempDisObj))] )
#       # make dist_obj as.dist() first, 
#       #   in this case it removes NAs and duplicated values
#       # dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
#     }))
#     # print(all_unit_labs[[ sprintf("%s.within-disOnly", unit_lab) ]])
#     
#     
#     
#     # with.bet.diff <- t.test(within_dist, between_dist)
#     # print(with.bet.diff)
#     
#     if (comp.method == "ANOSIM") {
#       unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.anosim[[ dist_name ]][[ unit_lab ]]$signif
#       unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- fam.anosim[[ dist_name ]][[ unit_lab ]]$statistic
#       
#     } else if (comp.method == "PERMANOVA") {
#       unit_labs_pvals[[ sprintf("%s.pval", unit_lab) ]] <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$`Pr(>F)`[ 1 ]
#       unit_labs_pvals[[ sprintf("%s.R2", unit_lab) ]]   <- sqrt(fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$R2[ 1 ])
#       unit_labs_pvals[[ sprintf("%s.F", unit_lab) ]]   <- fam.adonis[[ dist_name ]][[ unit_lab ]]$aov.tab$F.Model[ 1 ]
#     }
#     
#   }
#   
#   # print(all_unit_labs[[ sprintf("%s.within-disOnly", "Family_unit") ]])
#   
#   # ******* #
#   disSamps <- rownames(sam_tab[ sam_tab[ , disorder ] == "Yes", ] )
#   nonDisSamps <- rownames(sam_tab[ sam_tab[ , disorder ] == "No", ] )
#   
#   # get vector of dists for given disorder - members of disorder against same disorder
#   all_unit_labs[[ sprintf("dis_vs_dis") ]] <- as.numeric( as.dist( dist_obj[ disSamps, disSamps ] ))
#   # ******* #
#   
#   # ******* #
#   # get vector of dists for given disorder - members of disorder against same disorder
#   all_unit_labs[[ sprintf("dis_vs_cont") ]]  <- as.numeric( as.matrix( dist_obj)[ disSamps, nonDisSamps ] )#c()
#   # print(disSamps)
#   # print(nonDisSamps)
#   # for (s1 in disSamps) {
#   #   for (s2 in nonDisSamps) {
#   #     all_unit_labs[[ sprintf("dis_vs_cont") ]] <- c(all_unit_labs[[ sprintf("dis_vs_cont") ]], 
#   #                                                    as.numeric( as.matrix( dist_obj)[ s1, s2 ] ))
#   #   }
#   # }
#   # print(all_unit_labs[[ sprintf("dis_vs_cont") ]])
#   # ******* #
#   # print(str(all_unit_labs))
#   # print(str(unit_labs_pvals))
#   
#   # # get vector of units of particular type
#   # units <- unique( sam_tab[ , unit_lab] )
#   # units <- units[ units != "None" ]
#   # 
#   # # get vector of dists for each unit - members of unit against same unit
#   # within_dist <- unlist(sapply(units, function(x) {
#   #   unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] == x, ] )
#   #   # make dist_obj as.dist() first, 
#   #   #   in this case it removes NAs and duplicated values
#   #   dists <- as.numeric( as.dist( dist_obj[ unit_samps, unit_samps ] ))
#   # }))
#   # 
#   # # get vector of dists for each unit - members of unit against those not in that unit 
#   # #  (or against only other units??)
#   # between_dist <- unlist(sapply(units, function(x) {
#   #   unit_samps <- rownames(sam_tab[ sam_tab[ , unit_lab ] != x, ] )
#   #   # make dist_obj as.matrix() first, 
#   #   #   in this case, unlike for within_dist, it keeps all values in this subtable, since rows and columns are different
#   #   dists <- as.numeric( as.matrix(dist_obj[ unit_samps, colnames(dist_obj)[ ! colnames(dist_obj) %in% unit_samps ] ] ))
#   # }))
#   # # between_dist <- sapply(units, function(x) {
#   # #   unit_samps <- rownames(sam_tab[ ! sam_tab[ , unit_lab ] %in% c(x,"None"), ] )
#   # #   dists <- as.numeric( dist_obj[ unit_samps, unit_samps ])
#   # # }, sam_tab, unit_lab, dist_obj)
#   # 
#   # with.bet.diff <- t.test(within_dist, between_dist)
#   # print(with.bet.diff)
#   
#   print(unlist(lapply(all_unit_labs, mean)))
#   with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
#                          "labs" = names(all_unit_labs),
#                          # "labs" = rep(c("Same Family - dis","Same Family - con",
#                          #                "Different Families - dis","Different Families - con"), length(fam_units)),
#                          "sds"  = unlist(lapply(all_unit_labs, sd)))
#   # with.bet <- data.frame("vals" = unlist(lapply(all_unit_labs, mean)),
#   #                        "labs" = rep(c("Same Family","Same Disorder"), length(fam_units)),
#   #                        # "labs" = rep(c("Same Family - dis","Same Family - con",
#   #                        #                "Different Families - dis","Different Families - con"), length(fam_units)),
#   #                        "sds"  = unlist(lapply(all_unit_labs, sd)),
#   #                        "unit" = rep(gsub("_unit","",fam_units), each = 2))
#   
#   # with.bet$unit <- factor(with.bet$unit, levels = rev(gsub("_unit","",fam_units)))
#   
#   # with.bet <- data.frame("vals" = c(mean(within_dist), mean(between_dist)),
#   #                        "labs" = c("Within","Between"),
#   #                        "sds"  = c(sd(within_dist), sd(between_dist)),
#   #                        "unit" = c(unit_lab, unit_lab))
#   # "labs"=c( rep("Within", length(within_dist)), rep("Between",length(between_dist)) ))
#   # print(with.bet)
#   
#   star_labels <- sapply(fam_units, function(x) {
#     if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.001) {
#       sprintf("R=%s **", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       # sprintf("F=%s**", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
#     } else if (unit_labs_pvals[[ sprintf("%s.pval", x) ]] < 0.05) {
#       sprintf("R=%s *", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       # sprintf("F=%s*", round(unit_labs_pvals[[ sprintf("%s.F", x) ]], digits = 3))
#     } else {
#       # sprintf("R=%s", round(unit_labs_pvals[[ sprintf("%s.R2", x) ]], digits = 3))
#       ""
#     }
#   })
#   # print(star_labels)
#   print(with.bet)
#   print(star_labels)
#   ggplot(with.bet, aes(x=labs, y=vals, fill=labs)) +
#     geom_bar(stat="identity", position = "dodge2") +
#     guides(fill=guide_legend(title=NULL, reverse = T)) +
#     geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
#     coord_flip() +
#     # geom_text(aes(x=unit, y=vals+0.05))
#     # annotate("text", x=1:length(fam_units), y=with.bet$vals[with.bet$labs=="Different Families - dis"]+0.05, 
#     #          label = rev(star_labels), size=5) +
#     theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
#           legend.text = element_text(size=15), legend.position = "bottom") +
#     xlab(sprintf("%s in family units", comp.method)) + ylab(dist_name)
# }
# # ****************************** #
# 
# plot_famUnit_diffs.disorders3(mTab.Fams, as.matrix(ordObj$Aitchison), "Aitchison", "ANOSIM", disorder, fam.anosim = anosim_list)
# plot_famUnit_diffs(mTab.Fams, as.matrix(ordObj$Bray), "Bray", "ANOSIM", fam.anosim = anosim_list)

# ****************************************************************************************************************** #











# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# ANOSIM for a list of disorders against each other ####
anosim_disorder <- c("Cystic_fibrosis","Downs_Syndrome","Diabetes","Lactose_intolerant",
                     "Circulatory_issues","Lung_issues","Kidney_issues","Thyroid_issue",
                     "Intestinal_issues","Sinusitis","Gastritis","Immune_issues")
# ****************************************************************************************************************** #



# ANOSIM for disorders

weighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))
diag(weighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
unweighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))
diag(unweighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored

bray <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_bray.rds", p2_dir))
diag(bray) <- NA
jaccard <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_jaccard.rds", p2_dir))
diag(jaccard) <- NA

aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))


# ************************************************************************* #
# dis.anosim_list <- list()
# for (dist_meas in c("weighted_Unifrac","unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "bray","jaccard","aitch")) {
# 
#   dis.anosim_list[[ dist_meas ]] <- list()
# 
#   for (disorder in anosim_disorder) {
#     print(c(dist_meas, disorder))
#     
#     # get vector of dists for each unit - members of unit against same unit
#     dS     <- rownames(SLL2.meta[ SLL2.meta[ , disorder ] == "Yes", ] )
#     odS    <- rownames(SLL2.meta[ SLL2.meta[ , "Chronic_disorder" ] == "Yes" & 
#                                     SLL2.meta[ , disorder ] == "No", ] )
#     non.dS <- rownames(SLL2.meta[ SLL2.meta[ , disorder ] == "No", ] )
#     hS     <- rownames(SLL2.meta[ SLL2.meta[ , "Chronic_disorder" ] == "No" &
#                                     SLL2.meta[ , disorder ] == "No", ] )
#     
#     dis.anosim_list[[ dist_meas ]][[ disorder ]] <- list()
#     
#     dis.anosim_list[[ dist_meas ]][[ disorder ]][[ "other.disSamps" ]] <- 
#       anosim(as.dist( eval(parse(text = dist_meas))[ c(dS, odS), c(dS, odS) ] ),
#              as.character(SLL2.meta[ c(dS, odS), disorder]),
#              strata = SLL2.meta[ c(dS, odS), "seqGroup"])
#     
#     dis.anosim_list[[ dist_meas ]][[ disorder ]][[ "all_non.disSamps" ]] <- 
#       anosim(as.dist( eval(parse(text = dist_meas))[ c(dS, non.dS), c(dS, non.dS) ] ),
#              as.character(SLL2.meta[ c(dS, non.dS), disorder]),
#              strata = SLL2.meta[ c(dS, non.dS), "seqGroup"])
#     
#     dis.anosim_list[[ dist_meas ]][[ disorder ]][[ "healthySamps" ]] <- 
#       anosim(as.dist( eval(parse(text = dist_meas))[ c(dS, hS), c(dS, hS) ] ),
#              as.character(SLL2.meta[ c(dS, hS), disorder]),
#              strata = SLL2.meta[ c(dS, hS), "seqGroup"])
# 
#   }
# }
# 
# saveRDS(dis.anosim_list, sprintf("%s/R_objects/dis.anosim_list.rds", p2_dir))

dis.anosim_list <- readRDS(sprintf("%s/R_objects/dis.anosim_list.rds", p2_dir))

# ************************************************************************* #


age.child  <- ifelse(SLL2.meta$Age_groups == "Child", "Yes", "No")
age.teen   <- ifelse(SLL2.meta$Age_groups == "Teen", "Yes", "No")
age.adult  <- ifelse(SLL2.meta$Age_groups == "Adult", "Yes", "No")
age.senior <- ifelse(SLL2.meta$Age_groups == "Senior", "Yes", "No")

no.na.samps <- rownames(SLL2.meta[ ! is.na(SLL2.meta$Age_groups), ])

# ************************************************************************* #
# age.anosim_list <- list()
# for (dist_meas in c("weighted_Unifrac","unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "bray","jaccard","aitch")) {
#   
#   print(dist_meas)
#   
#   age.anosim_list[[ dist_meas ]] <- list()
#   
#   age.anosim_list[[ dist_meas ]][[ "Children" ]]  <- anosim(as.dist( eval(parse(text = dist_meas))[no.na.samps, no.na.samps] ),
#                                                             age.child[ ! is.na(age.child) ],
#                                                             strata = SLL2.meta[ no.na.samps, "seqGroup"])
#   age.anosim_list[[ dist_meas ]][[ "Teens" ]]   <- anosim(as.dist( eval(parse(text = dist_meas))[no.na.samps, no.na.samps] ),
#                                                           age.teen[ ! is.na(age.teen) ],
#                                                           strata = SLL2.meta[ no.na.samps, "seqGroup"])
#   age.anosim_list[[ dist_meas ]][[ "Adults" ]]  <- anosim(as.dist( eval(parse(text = dist_meas))[no.na.samps, no.na.samps] ),
#                                                           age.adult[ ! is.na(age.adult) ],
#                                                           strata = SLL2.meta[ no.na.samps, "seqGroup"])
#   age.anosim_list[[ dist_meas ]][[ "Seniors" ]] <- anosim(as.dist( eval(parse(text = dist_meas))[no.na.samps, no.na.samps] ),
#                                                           age.senior[ ! is.na(age.senior) ],
#                                                           strata = SLL2.meta[ no.na.samps, "seqGroup"])
#   
# }
# 
# saveRDS(age.anosim_list, sprintf("%s/R_objects/age.anosim_list.rds", p2_dir))

age.anosim_list <- readRDS(sprintf("%s/R_objects/age.anosim_list.rds", p2_dir))


# ************************************************************************* #







# For disorders
plot_anosim_groups(SLL2.meta, aitch, "aitch", "ANOSIM", rev(anosim_disorder), "Disorders",
                   anoList = dis.anosim_list)
plot_anosim_groups(SLL2.meta, aitch, "aitch", "ANOSIM", rev(anosim_disorder), "Disorders",
                   anoList = dis.anosim_list, plotType = "box")

plot_anosim_groups(SLL2.meta, weighted_Unifrac, "weighted_Unifrac", "ANOSIM", rev(anosim_disorder),  "Disorders",
                   anoList = dis.anosim_list)

plot_anosim_groups(SLL2.meta, jaccard, "jaccard", "ANOSIM", rev(anosim_disorder),  "Disorders",
                   anoList = dis.anosim_list)

# ****************************** #



# For age_groups
mTab.ageGroups <- SLL2.meta[ no.na.samps, ]
mTab.ageGroups$Children <- ifelse(mTab.ageGroups$Age_groups == "Child", "Yes", "No")
mTab.ageGroups$Teens    <- ifelse(mTab.ageGroups$Age_groups == "Teen", "Yes", "No")
mTab.ageGroups$Adults   <- ifelse(mTab.ageGroups$Age_groups == "Adult", "Yes", "No")
mTab.ageGroups$Seniors  <- ifelse(mTab.ageGroups$Age_groups == "Senior", "Yes", "No")

plot_anosim_groups(mTab.ageGroups, aitch[ no.na.samps, no.na.samps ], "aitch", "ANOSIM", 
                   rev(c("Children","Teens","Adults","Seniors")), "Age_groups", 
                   anoList = age.anosim_list)
plot_anosim_groups(mTab.ageGroups, aitch[ no.na.samps, no.na.samps ], "aitch", "ANOSIM", 
                   rev(c("Children","Teens","Adults","Seniors")), "Age_groups", 
                   anoList = age.anosim_list, plotType = "box")

plot_anosim_groups(mTab.ageGroups, weighted_Unifrac[ no.na.samps, no.na.samps ], "weighted_Unifrac", "ANOSIM", 
                   rev(c("Children","Teens","Adults","Seniors")), "Age_groups", 
                   anoList = age.anosim_list)

plot_anosim_groups(mTab.ageGroups, jaccard[ no.na.samps, no.na.samps ], "jaccard", "ANOSIM", 
                   rev(c("Children","Teens","Adults","Seniors")), "Age_groups", 
                   anoList = age.anosim_list)

# ************************************************************************* #




# ANOSIM for community
healthySamps <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No", ])
meta.healthy <- SLL2.meta[ healthySamps, ]

# anosim.community <- list()
# anosim.city <- list()
# for (dist_meas in c("weighted_Unifrac","unweighted_Unifrac",
#                     # "guni.VAW","guni.a0","guni.a05","jsd","canberra",
#                     "bray","jaccard","aitch")) {
#   
#   anosim.community[[ dist_meas ]] <- anosim(as.dist( eval(parse(text = dist_meas))[healthySamps, healthySamps] ),
#                                             as.character(SLL2.meta[ healthySamps, "Community"]),
#                                             strata = SLL2.meta[ healthySamps, "seqGroup"])
#   
#   anosim.city[[ dist_meas ]] <- anosim(as.dist( eval(parse(text = dist_meas))[healthySamps, healthySamps] ),
#                                        as.character(SLL2.meta[ healthySamps, "City"]),
#                                        strata = SLL2.meta[ healthySamps, "seqGroup"])
# }
# saveRDS(anosim.community, file = sprintf("%s/R_objects/anosim.community.rds", p2_dir))
# saveRDS(anosim.city, file = sprintf("%s/R_objects/anosim.city.rds", p2_dir))

anosim.community <- readRDS(sprintf("%s/R_objects/anosim.community.rds", p2_dir))
anosim.city <- readRDS(sprintf("%s/R_objects/anosim.city.rds", p2_dir))


# ************************************************************************* #
library(ggpubr)
plot_location_anosim <- function(loc_type, dist_obj, dist_meas, sam_tab, plotType) {
  
  locs <- sort(unique(sam_tab[ , loc_type ]))
  
  loc.dists <- list()
  
  for (loc in locs) {
    
    loc.Samps <- rownames(sam_tab[ sam_tab[ , loc_type ] == loc, ] )
    other_loc.Samps <- rownames(sam_tab[ sam_tab[ , loc_type ] != loc, ] )
    
    loc.dists[[ sprintf("%s.within",loc) ]]    <- as.numeric( as.dist(dist_obj[ loc.Samps, loc.Samps ] ))
    loc.dists[[ sprintf("%s.other_loc",loc) ]] <- as.numeric( as.matrix(dist_obj[ loc.Samps, other_loc.Samps ] ))
    
  }
  
  # ****************************** #
  
  plot_labs <- c(sprintf("Same %s", loc_type), sprintf("Other %s", loc_type))
  unit_reps <- 2
  plot_xLab <- sprintf("%s amongst %s", "ANOSIM", ifelse(loc_type=="Community","communities",
                                                         ifelse(loc_type=="Province","provinces","cities")))
  
  if (plotType == "bar") {
    with.bet <- data.frame("vals" = unlist(lapply(loc.dists, mean)),
                           "labs" = rep(plot_labs, length(locs)),
                           "sds"  = unlist(lapply(loc.dists, sd)),
                           "unit" = rep(locs, each = unit_reps))
    
  } else if (plotType == "box") {
    wb.labs <- sapply(names(loc.dists), function(x) {
      if (endsWith(x, "within")) rep(plot_labs[1], length(loc.dists[[ x ]]))
      else if (endsWith(x, "other_loc") | endsWith(x, "other_age")) rep(plot_labs[2], length(loc.dists[[ x ]]))
      # else if (endsWith(x, "non_dis")) rep(plot_labs[3], length(loc.dists[[ x ]]))
      # else if (endsWith(x, "healthy")) rep(plot_labs[4], length(loc.dists[[ x ]]))
    })
    
    wb.units <- sapply(locs, function(x) {
      aul.names <- names(loc.dists)[ startsWith(names(loc.dists), x) ]
      rep(x, sum(sapply(aul.names, function(y) length(loc.dists[[ y ]]))))
    })
    
    with.bet <- data.frame("vals" = unlist(loc.dists),
                           "labs" = unlist(wb.labs),
                           "unit" = unlist(wb.units))
  }
  
  with.bet$labs <- factor(with.bet$labs, levels = rev(plot_labs))
  with.bet$unit <- factor(with.bet$unit, levels = rev(locs))
  # ****************************** #
  
  
  if (plotType == "bar") {
    ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_bar(stat="identity", position = "dodge2") +
      geom_errorbar(aes(ymin=vals-sds, ymax=vals+sds), position = position_dodge2(width = 0.2, padding = 0.8)) +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      # annotate("text",
      #          x=1:length(c(cols_to_check)),
      #          y=with.bet$vals[with.bet$labs %in% c("Healthy","Other Age Groups")]+0.05,
      #          label = rev(star_labels), size=5) +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      xlab(plot_xLab) + ylab(dist_meas)
    
  } else if (plotType == "box") {
    
    p.x <- ifelse(dist_meas=="Aitchison", 60, 
                  ifelse(dist_meas=="Weighted_Unifrac", 0.4,
                         ifelse(dist_meas=="Unweighted_Unifrac", 0.6,
                                ifelse(dist_meas %in% c("Bray","Jaccard"), 0.8, 0.5))))
    
    boxes <- ggplot(with.bet, aes(x=unit, y=vals, fill=labs)) +
      geom_boxplot(notch = T) +
      guides(fill=guide_legend(title=NULL, reverse = T)) +
      coord_flip() +
      # annotate("text",
      #          x=1:length(c(cols_to_check)),
      #          # y=max(with.bet$vals[with.bet$labs %in% c("Healthy","Other Age Groups")]),
      #          y=max(with.bet$vals[ startsWith(as.character(with.bet$labs), "Same") ]),
      #          label = rev(star_labels), size=6,
      #          col="red") +
      theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
            legend.text = element_text(size=15)) +#, legend.position = "bottom") +
      stat_compare_means(aes(group=labs, label=paste0("p = ",..p.format..)), color="red", hide.ns=T, method="kruskal",
                         # label.x = c(5, 4, 3, 2, 1), 
                         label.y =rep(p.x, length(locs))) +
      xlab(plot_xLab) + ylab(dist_meas)
    
    # return(stat_compare_means(aes(group=labs, label=paste0("p = ",..p.format..)), color="red", hide.ns=T, method="kruskal",
    #                           # label.x = c(5, 4, 3, 2, 1), 
    #                           label.y =rep(p.x, length(locs))))
    boxes
    return(list("loc.dists"=loc.dists, "boxes"=boxes))
  }
  # ****************************** #
}
# ************************************************************************* #

plot_location_anosim("Community", aitch, "Aitchison", SLL2.meta[healthySamps, ], "box")
plot_location_anosim("Community", weighted_Unifrac, "Weighted_Unifrac", SLL2.meta[healthySamps, ], "box")
plot_location_anosim("Community", unweighted_Unifrac, "Unweighted_Unifrac", SLL2.meta[healthySamps, ], "box")
plot_location_anosim("Community", bray, "Bray", SLL2.meta[healthySamps, ], "box")
plot_location_anosim("Community", jaccard, "Jaccard", SLL2.meta[healthySamps, ], "box")


# Cantabria is most consistently similar community, while Pais Vasco is the only that is not significantly more similar to itself than to others
#   may suggest less diversity overall in Cantabria and much more in Pais Vasco - lets check
group_vs_cont_box(SLL2.meta[healthySamps,], healthySamps, "Genus", gloms_clr$Genus, "Community", "Div.Shannon", 
                  groupQs, only_cont, xAngle=45)
# Indeed Pais Vasco is a bit higher than other communities, and Cantabria is much lower (though only 9 healthySamps here)
table(SLL2.meta[healthySamps, "Community"])


# ************************************************************************* #
city.aitch <- plot_location_anosim("City", aitch, "Aitchison", SLL2.meta[healthySamps, ], "box")

# check for connections between similarity within a community and the water values
# Strange cities include Alhaurí de la Torre (not sig), Figueres, Valencia, Villarreal de Urrechua (not sig)
# Cities with very consistent similarity include Colmenar Viejo, L'Alcora, La Paca, Moralzarzal, Palma de Mallorca, Roda de Ter, Santander

group_vs_cont_box(SLL2.meta[healthySamps,], healthySamps, "Genus", gloms_clr$Genus, "City", "Div.Shannon", 
                  groupQs, only_cont, xAngle=45, plot_tukey = F)
# Villarreal de Urrechua has high Shannon and Alhaurin de la Torre somewhat, but these 2 were not signif in the anosim plot
#   Figueres and Valencia were signif in the anosim plot and have somewhat high Shannon,  but close to median
# those cities that remain consistent do not however have particularly low Shannon, they are generally around the average

group_vs_cont_box(SLL2.meta[healthySamps,], healthySamps, "Genus", gloms_clr$Genus, "City", "Div.Simpson", 
                  groupQs, only_cont, xAngle=45, plot_tukey = F)

# So this trend of high diversity == less similarty within a location does not hold as well at the City level as it does at the Community level

table(SLL2.meta[healthySamps, "City"])

# ************************************************************************************************ #














# ************************************************************************************************ #
# Check correlations between distances and differences in waterVals ####


# *************************************************************** #
dists_and_water_difs <- function(waterVal, sam_tab, dist_obj) {
  
  # *************** #
  # first have to make a table of sample vs sample to get the difference in the water value between those
  water_difs <- matrix(NA, nrow = nrow(sam_tab), ncol = nrow(sam_tab))
  rownames(water_difs) <- colnames(water_difs) <- rownames(sam_tab)
  
  for (samp1 in rownames(water_difs)) {
    for (samp2 in rownames(water_difs)) {
      # must be absolute value of difference
      water_difs[ samp1, samp2 ] <- abs(sam_tab[samp1, waterVal] - sam_tab[samp2, waterVal])
    }
  }
  # *************** #
  
  # *************** #
  loc.dists <- list()
  loc.wds   <- list()
  # cities <- sort(unique(sam_tab[, "City"]))
  # cities <- cities[ nrow(sam_tab[ sam_tab[,"City"]==cities, ]) > 5 ]
  for (loc in sort(unique(sam_tab[, "City"]))) {
    
    loc.Samps <- rownames(sam_tab[ sam_tab[ , "City" ] == loc, ] )
    other_loc.Samps <- rownames(sam_tab[ sam_tab[ , "City" ] != loc, ] )
    
    # loc.dists[[ sprintf("%s.within",loc) ]]    <- as.numeric( as.dist(dist_obj[ loc.Samps, loc.Samps ] ))
    loc.dists[[ loc ]] <- as.numeric( as.matrix(dist_obj[ loc.Samps, other_loc.Samps ] ))
    loc.wds[[ loc ]]   <- as.numeric( water_difs[ loc.Samps, other_loc.Samps ] )
    
  }
  # *************** #
  return(list("loc.dists"=loc.dists, "loc.wds"=loc.wds))
}
# *************************************************************** #


# *************************************************************** #

healthy.nonBottle <- healthySamps[ ! is.na(SLL2.meta[healthySamps, "Water_type_home"]) &
                                     SLL2.meta[healthySamps, "Water_type_home"] != "Embotellada" ]
meta.healthy_nonBottle <- SLL2.meta[ healthy.nonBottle, ]


# dwd_cors <- list()
# for (dist_meas in c("aitch","weighted_Unifrac","unweighted_Unifrac","bray","jaccard")) {
#   print(dist_meas)
#   dwd_cors[[ dist_meas ]] <- list()
#   
#   for (waterVal in cont_water_data) {
#     print(sprintf("   %s", waterVal))
#     dwd_cors[[ dist_meas ]][[ waterVal ]] <- list()
#     
#     # first get distance values between samples from each location vs other locations, 
#     #    and differences in given water values between those samples
#     dwd <- dists_and_water_difs(waterVal, SLL2.meta[healthy.nonBottle, ], eval(parse(text = dist_meas)))
#     
#     for (corType in c("pearson","spearman","kendall")) {
#       # then do indicated correlations between distance values and water differences
#       corTab <- t(sapply(sort(unique(SLL2.meta[healthy.nonBottle, "City"])), function(loc) {
#         ct <- cor.test(dwd$loc.dists[[ loc ]], dwd$loc.wds[[ loc ]], method = corType)
#         # c(round(ct$p.value,5), round(ct$estimate,5))
#         c(ct$p.value, ct$estimate)
#       }))
#       colnames(corTab) <- c("pval","cor")
#       dwd_cors[[ dist_meas ]][[ waterVal ]][[ corType ]] <- corTab
#     }
#     
#   }
#   # saveRDS(dwd_cors, file = sprintf("%s/R_objects/dwd_cors.rds", p2_dir))
# }
# saveRDS(dwd_cors, file = sprintf("%s/R_objects/dwd_cors.rds", p2_dir))
# # *************************************************************** #

dwd_cors <- readRDS(file = sprintf("%s/R_objects/dwd_cors.rds", p2_dir))




interesting.cities <- c("Alhaurín de la Torre","Figueres","Valencia","Villarreal de Urrechua",
                        "Colmenar Viejo","L'Alcora","La Paca","Moralzarzal","Palma de Mallorca","Roda de Ter","Santander")
lapply(dwd_cors$aitch$Conductivity, function(x) x[interesting.cities,])


# *************************************************************** #
# those cities that had a correlation value of at least 0.2 for indicated waterVal
sig.dwd_cors <- list()
for (dist_meas in names(dwd_cors)) {
  sig.dwd_cors[[ dist_meas ]] <- list()
  
  for (waterVal in cont_water_data) {
    sig.dwd_cors[[ dist_meas ]][[ waterVal ]] <- lapply(dwd_cors[[ dist_meas ]][[ waterVal ]], 
                                                        function(x) as.matrix(as.data.frame(x)[abs(x[,"cor"])>0.175,]))
  }
}

# *************************************************************** #
# and those that are not statistically significant
non_sig.dwd_cors <- list()
for (dist_meas in names(dwd_cors)) {
  non_sig.dwd_cors[[ dist_meas ]] <- list()
  
  for (waterVal in cont_water_data) {
    non_sig.dwd_cors[[ dist_meas ]][[ waterVal ]] <- lapply(dwd_cors[[ dist_meas ]][[ waterVal ]], 
                                                            function(x) as.matrix(as.data.frame(x)[abs(x[,"pval"])>0.05,]))
  }
}

# ****************************************************************************************************************** #
























# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####
# Testing effects with age in DS samples ####


# ****************************************************************************************************************** #
# 1) Are DS more similar to older bins of healthy population than younger bins? ####
#      Since DS is supposed to be similar to a premature aging disease in some aspects, would suggest more similar to older


DS.samps <- rownames(SLL2.meta[ SLL2.meta$Downs_Syndrome=="Yes", ])
healthyPop <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder=="No", ])
# ignore the 2 samples without Age values
healthyPop <- healthyPop[ ! is.na(SLL2.meta[ healthyPop, "Age_groups"]) ]

phy.ds.h <- prune_samples(c(DS.samps, healthyPop), SLL2)
mTab.ds.h <- SLL2.meta[ c(DS.samps, healthyPop), ]

# ordObj.ds.h <- subsampling_ordination_objects(c(DS.samps, healthyPop), phy.ds.h, distsOnly=T)
# saveRDS(ordObj.ds.h, file = sprintf("%s/R_objects/ordObj.ds.h.rds", p2_dir))
ordObj.ds.h <- readRDS(sprintf("%s/R_objects/ordObj.ds.h.rds", p2_dir))

# ******************************** #


# ******************************** #
# plot mean dists between DS and each age group in healthyPop
# ******************************** #
get_dists_to_ageGroups <- function(ordObj, dist_meas, mTab, dss, hp) {
  
  # Aitchison distances on scale of 0-100, others 0-1, so will have to adjust Aitchison to be comparable
  mult <- 1
  if (dist_meas=="Aitchison")
    mult <- 100
  
  Child  <- hp[ mTab[ hp, "Age_groups"] == "Child" ]
  Teen   <- hp[ mTab[ hp, "Age_groups"] == "Teen" ]
  Adult  <- hp[ mTab[ hp, "Age_groups"] == "Adult" ]
  Senior <- hp[ mTab[ hp, "Age_groups"] == "Senior" ]
  
  return( list( "Child"  = as.numeric(ordObj[dss, Child]) / mult, 
                "Teen"   = as.numeric(ordObj[dss, Teen]) / mult, 
                "Adult"  = as.numeric(ordObj[dss, Adult]) / mult, 
                "Senior" = as.numeric(ordObj[dss, Senior]) / mult ))
}
# ******************************** #

ag.dists <- list()
ag.dists[[ "Aitchison" ]]          <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Aitchison), "Aitchison", mTab.ds.h, DS.samps, healthyPop)
ag.dists[[ "Weighted_Unifrac" ]]   <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Weighted_Unifrac), "Weighted_Unifrac", mTab.ds.h, DS.samps, healthyPop)
ag.dists[[ "Unweighted_Unifrac" ]] <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Unweighted_Unifrac), "Unweighted_Unifrac", mTab.ds.h, DS.samps, healthyPop)
ag.dists[[ "Bray" ]]               <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Bray), "Bray", mTab.ds.h, DS.samps, healthyPop)
ag.dists[[ "Jaccard" ]]            <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Jaccard), "Jaccard", mTab.ds.h, DS.samps, healthyPop)

ag.dists.m <- reshape2::melt(ag.dists)
colnames(ag.dists.m) <- c("value","Age_group","Distance")
ag.dists.m$Age_group <- factor(ag.dists.m$Age_group, levels=rev(c("Child","Teen","Adult","Senior")))
ag.dists.m$Distance <- factor(ag.dists.m$Distance, levels=rev(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                                                                "Bray","Jaccard")))

# for adult DS only:
DS.adult <- DS.samps[ mTab.ds.h[ DS.samps, "Age_groups"]=="Adult" ]
ag.dists.adult <- list()
ag.dists.adult[[ "Aitchison" ]]          <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Aitchison), "Aitchison", mTab.ds.h, DS.adult, healthyPop)
ag.dists.adult[[ "Weighted_Unifrac" ]]   <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Weighted_Unifrac), "Weighted_Unifrac", mTab.ds.h, DS.adult, healthyPop)
ag.dists.adult[[ "Unweighted_Unifrac" ]] <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Unweighted_Unifrac), "Unweighted_Unifrac", mTab.ds.h, DS.adult, healthyPop)
ag.dists.adult[[ "Bray" ]]               <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Bray), "Bray", mTab.ds.h, DS.adult, healthyPop)
ag.dists.adult[[ "Jaccard" ]]            <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Jaccard), "Jaccard", mTab.ds.h, DS.adult, healthyPop)

ag.dists.adult.m <- reshape2::melt(ag.dists.adult)
colnames(ag.dists.adult.m) <- c("value","Age_group","Distance")
ag.dists.adult.m$Age_group <- factor(ag.dists.adult.m$Age_group, levels=rev(c("Child","Teen","Adult","Senior")))
ag.dists.adult.m$Distance <- factor(ag.dists.adult.m$Distance, levels=rev(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                                                                            "Bray","Jaccard")))

# for DS youths only
DS.youth <- DS.samps[ mTab.ds.h[ DS.samps, "Age_groups"]!="Adult" ]
ag.dists.youth <- list()
ag.dists.youth[[ "Aitchison" ]]          <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Aitchison), "Aitchison", mTab.ds.h, DS.youth, healthyPop)
ag.dists.youth[[ "Weighted_Unifrac" ]]   <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Weighted_Unifrac), "Weighted_Unifrac", mTab.ds.h, DS.youth, healthyPop)
ag.dists.youth[[ "Unweighted_Unifrac" ]] <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Unweighted_Unifrac), "Unweighted_Unifrac", mTab.ds.h, DS.youth, healthyPop)
ag.dists.youth[[ "Bray" ]]               <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Bray), "Bray", mTab.ds.h, DS.youth, healthyPop)
ag.dists.youth[[ "Jaccard" ]]            <- get_dists_to_ageGroups(as.matrix(ordObj.ds.h$Jaccard), "Jaccard", mTab.ds.h, DS.youth, healthyPop)

ag.dists.youth.m <- reshape2::melt(ag.dists.youth)
colnames(ag.dists.youth.m) <- c("value","Age_group","Distance")
ag.dists.youth.m$Age_group <- factor(ag.dists.youth.m$Age_group, levels=rev(c("Child","Teen","Adult","Senior")))
ag.dists.youth.m$Distance <- factor(ag.dists.youth.m$Distance, levels=rev(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                                                                            "Bray","Jaccard")))




library(ggpubr)
ggplot(ag.dists.m, aes(x=Distance, y=value, fill=Age_group)) +
  geom_boxplot(notch = T) + coord_flip() + ylim(0, 1.1) +
  guides(fill=guide_legend(title=NULL, reverse = T)) +
  stat_compare_means(aes(group=Age_group, label=..p.signif..), color="red", hide.ns=T, method="kruskal",
                     label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
  ggtitle("DS vs Age Groups") +
  theme_classic() + ylab(NULL) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size=15), legend.position = "bottom")

ggplot(ag.dists.adult.m, aes(x=Distance, y=value, fill=Age_group)) +
  geom_boxplot(notch = T) + coord_flip() + ylim(0, 1.1) +
  guides(fill=guide_legend(title=NULL, reverse = T)) +
  stat_compare_means(aes(group=Age_group, label=..p.signif..), color="red", hide.ns=T, method="kruskal",
                     label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
  ggtitle("(a) DS vs Age Groups - DS adults only") +
  theme_classic() + ylab(NULL) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size=15), legend.position = "bottom")

ggplot(ag.dists.youth.m, aes(x=Distance, y=value, fill=Age_group)) +
  geom_boxplot(notch = T) + coord_flip() + ylim(0, 1.1) +
  guides(fill=guide_legend(title=NULL, reverse = T)) +
  stat_compare_means(aes(group=Age_group, label=..p.signif..), color="red", hide.ns=T, method="kruskal",
                     label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
  ggtitle("(b) DS vs Age Groups - DS youths only") +
  theme_classic() + ylab(NULL) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size=15), legend.position = "bottom")

# ******************************** #


tukey_kruskal <- function(dist_meas, agList, sigOnly=F) {
  mag <- reshape2::melt(agList[[ dist_meas ]])
  mag$L1 <- factor(mag$L1)
  ra <- aov(value~L1, data = mag)
  if (sigOnly==T)
    print(TukeyHSD(ra)$L1[ TukeyHSD(ra)$L1[,"p adj"]<0.05, ])
  else
    print(TukeyHSD(ra))
  print(kruskal.test(value~L1, data = mag))
}

tukey_kruskal("Aitchison", ag.dists)
tukey_kruskal("Weighted_Unifrac", ag.dists)
tukey_kruskal("Unweighted_Unifrac", ag.dists)
tukey_kruskal("Bray", ag.dists)
tukey_kruskal("Jaccard", ag.dists)


# ******************************** #
# cov_vars <- c("Age_groups","Gender","Population","seqGroup")
# anova.ds.age <- list()
# for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
#                     "Bray","Jaccard")) {
#   
#   anova.ds.age[[ dist_meas ]] <- get_lm( c(f, fixeds[ fixeds != f ]), tl, SLL2.meta, glomTab, dv )# get_adonis("", cov_vars, dist_meas, mTab.ds.h)
# }

# ****************************************************************************************************************** #
# Anosim for age groups among healthy samples ####

anosim.ageGroups.H <- list()
for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                    "Bray","Jaccard")) {
  print(dist_meas)
  
  anosim.ageGroups.H[[ dist_meas ]] <- anosim(as.dist( as.matrix(ordObj.ds.h[[ dist_meas ]])[healthyPop, healthyPop] ),
                                              as.character(as.matrix(mTab.ds.h[healthyPop, "Age_groups" ])))#,
  # strata = SLL2.meta[ , "seqGroup"])
  
  # mTab <- SLL2.meta[ SLL2.meta[,unit_lab] != "None", ]
  # anosim_list[[ dist_meas ]][[ unit_lab ]] <- anosim(as.dist( eval(parse(text = dist_meas))[rownames(mTab),rownames(mTab)] ), 
  #                                                    as.character(as.matrix(mTab[, unit_lab])))
}

ag.c <- healthyPop[ mTab.ds.h[ healthyPop, "Age_groups"] == "Child" ]
ag.t <- healthyPop[ mTab.ds.h[ healthyPop, "Age_groups"] == "Teen" ]
ag.a <- healthyPop[ mTab.ds.h[ healthyPop, "Age_groups"] == "Adult" ]
ag.s <- healthyPop[ mTab.ds.h[ healthyPop, "Age_groups"] == "Senior" ]

diffs <- list()
for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
  diffs[[ dist_meas ]] <- list()
  
  diffs[[ dist_meas ]][[ "Child" ]] <- list("vs_Child"  = as.numeric(as.dist( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.c, ag.c ])),
                                            "vs_Teen"   = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.c, ag.t ]),
                                            "vs_Adult"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.c, ag.a ]),
                                            "vs_Senior" = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.c, ag.s ]))
  
  diffs[[ dist_meas ]][[ "Teen" ]]  <- list("vs_Child"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.t, ag.c ]),
                                            "vs_Teen"   = as.numeric(as.dist( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.t, ag.t ])),
                                            "vs_Adult"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.t, ag.a ]),
                                            "vs_Senior" = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.t, ag.s ]))
  
  diffs[[ dist_meas ]][[ "Adult" ]] <- list("vs_Child"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.a, ag.c ]),
                                            "vs_Teen"   = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.a, ag.t ]),
                                            "vs_Adult"  = as.numeric(as.dist( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.a, ag.a ])),
                                            "vs_Senior" = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.a, ag.s ]))
  
  diffs[[ dist_meas ]][[ "Senior" ]] <- list("vs_Child"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.s, ag.c ]),
                                             "vs_Teen"   = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.s, ag.t ]),
                                             "vs_Adult"  = as.numeric( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.s, ag.a ]),
                                             "vs_Senior" = as.numeric(as.dist( as.matrix(ordObj.ds.h[[ dist_meas ]])[ ag.s, ag.s ])))
  
}

dist_meas <- "Weighted_Unifrac"

dm <- melt(diffs[[ dist_meas ]])
dm$L2 <- factor(dm$L2, levels = c("vs_Child","vs_Teen","vs_Adult","vs_Senior"))
dm$L1 <- factor(dm$L1, levels = c("Child","Teen","Adult","Senior"))

ggplot( dm, aes(x=L1, y=value, fill=L2)) +
  geom_boxplot(notch = T) +
  stat_compare_means(aes(group=L2, label=..p.signif..), hide.ns=T, method="anova") +
  ggtitle(dist_meas)




# ****************************************************************************************************************** #






# ****************************************************************************************************************** #
# 2) Are differences across age are more extreme in DS than healthy? ####
#      i.e. can we see more rapid aging effects present in the DS microbiome?

# list of ageDiff vs dist correlations for all DS samples:
dif.dis.DS <- list()
for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
  
  # for each sample, column for diff bt that samps age and another samp, other column for distance bt those 2 samps
  age_diffs <- as.numeric(sapply(DS.samps, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[DS.samps[DS.samps!=x], "Age"] )))
  dist_vals <- as.numeric(sapply(DS.samps, function(x) as.numeric(as.matrix(ordObj.ds.h[[ dist_meas ]])[x, DS.samps[DS.samps!=x] ])) )
  
  dif.dis.DS[[ dist_meas ]] <- list("vals"=cbind(age_diffs, dist_vals),
                                    "corTest"=cor.test(age_diffs, dist_vals))
  
}



# list of ageDiff vs dist correlations for all healthyPop samples:
dif.dis.healthyPop <- list()
for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
  
  age_diffs <- as.numeric(sapply(healthyPop, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[healthyPop[healthyPop!=x], "Age"] )))
  dist_vals <- as.numeric(sapply(healthyPop, function(x) as.numeric(as.matrix(ordObj.ds.h[[ dist_meas ]])[x, healthyPop[healthyPop!=x] ])) )
  
  dif.dis.healthyPop[[ dist_meas ]] <- list("vals"=cbind(age_diffs, dist_vals),
                                            "corTest"=cor.test(age_diffs, dist_vals))
}


# running correlation for all controls used in subsampling, only those up to age 33, to have same trajectory as DS
dsc.33 <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Downs_Syndrome"]=="No" ]
dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Age"] < 34 ]

dif.dis.matched33 <- list()
for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
  
  age_diffs <- as.numeric(sapply(dsc.33, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc.33[dsc.33!=x], "Age"] )))
  dist_vals <- as.numeric(sapply(dsc.33, function(x) as.numeric(as.matrix(ordObj.ds.h[[ dist_meas ]])[x, dsc.33[dsc.33!=x] ])) )
  
  dif.dis.matched33[[ dist_meas ]] <- list("vals"=cbind(age_diffs, dist_vals),
                                           "corTest"=cor.test(age_diffs, dist_vals))
}


# ************************************************************************************** #



# for each sample, column for diff bt that samps age and another samp, other column for distance bt those 2 samps
dif.dis.DS <- cbind(as.numeric(sapply(DS.samps, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[DS.samps[DS.samps!=x], "Age"] ))),
                    as.numeric(sapply(DS.samps, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, DS.samps[DS.samps!=x] ])) ))
colnames(dif.dis.DS) <- c("age_diffs","dists")
saveRDS(dif.dis.DS, file = sprintf("%s/R_objects/dif.dis.DS.rds", p2_dir))

plot(dif.dis.DS)
cor.test(dif.dis.DS[,"age_diffs"], dif.dis.DS[,"dists"])


# running correlation for all healthyPop
dif.dis.healthyPop <- cbind(as.numeric(sapply(healthyPop, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[healthyPop[healthyPop!=x], "Age"] ))),
                            as.numeric(sapply(healthyPop, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, healthyPop[healthyPop!=x] ])) ))
colnames(dif.dis.healthyPop) <- c("age_diffs","dists")
saveRDS(dif.dis.healthyPop, file = sprintf("%s/R_objects/dif.dis.healthyPop.rds", p2_dir))

plot(dif.dis.healthyPop)
cor.test(dif.dis.healthyPop[,"age_diffs"], dif.dis.healthyPop[,"dists"])


# correlations for all healthPop under 34
hp33 <- healthyPop[ mTab.ds.h[healthyPop, "Age"] < 34 ]
dif.dis.hp33 <- cbind(as.numeric(sapply(hp33, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[hp33[hp33!=x], "Age"] ))),
                      as.numeric(sapply(hp33, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, hp33[hp33!=x] ])) ))
colnames(dif.dis.hp33) <- c("age_diffs","dists")
saveRDS(dif.dis.hp33, file = sprintf("%s/R_objects/dif.dis.hp33.rds", p2_dir))

plot(dif.dis.hp33)
cor.test(dif.dis.hp33[,"age_diffs"], dif.dis.hp33[,"dists"])





dsc <- disSubsTests$default$Downs_Syndrome$`1`$samples
dsc <- dsc[ mTab.ds.h[dsc, "Downs_Syndrome"]=="No" ]
dif.dis <- cbind(as.numeric(sapply(dsc, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc[dsc!=x], "Age"] ))),
                 as.numeric(sapply(dsc, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, dsc[dsc!=x] ])) ))
colnames(dif.dis) <- c("age_diffs","dists")

plot(dif.dis)
cor.test(dif.dis[,"age_diffs"], dif.dis[,"dists"])





# running correlation for controls in each subsampling separately:
dif.dis.conts <- list()
for (i in names(disSubsTests$default$Downs_Syndrome)) {
  dsc <- disSubsTests$default$Downs_Syndrome[[ i ]]$samples
  dsc <- dsc[ mTab.ds.h[dsc, "Downs_Syndrome"]=="No" ]
  
  dif.dis <- cbind(as.numeric(sapply(dsc, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc[dsc!=x], "Age"] ))),
                   as.numeric(sapply(dsc, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, dsc[dsc!=x] ])) ))
  colnames(dif.dis) <- c("age_diffs","dists")
  
  dif.dis.conts[[ i ]] <- cor.test(dif.dis[,"age_diffs"], dif.dis[,"dists"])
  
}
saveRDS(dif.dis.conts, file = sprintf("%s/R_objects/dif.dis.conts.rds", p2_dir))

mean(unlist(lapply(dif.dis.conts, function(x) x$p.value)))
mean(unlist(lapply(dif.dis.conts, function(x) x$estimate)))


hist(-log(unlist(lapply(dif.dis.conts, function(x) x$p.value))))
hist(unlist(lapply(dif.dis.conts, function(x) x$estimate)))

plot(-log(unlist(lapply(dif.dis.conts, function(x) x$p.value))),
     unlist(lapply(dif.dis.conts, function(x) x$estimate)))



# running correlation for controls in each subsampling separately, only those up to age 33, to have same trajectory as DS:
dif.dis.conts.33 <- list()
for (i in names(disSubsTests$default$Downs_Syndrome)) {
  dsc.33 <- disSubsTests$default$Downs_Syndrome[[ i ]]$samples
  dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Downs_Syndrome"]=="No" ]
  dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Age"] < 34 ]
  
  dif.dis.33 <- cbind(as.numeric(sapply(dsc.33, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc.33[dsc.33!=x], "Age"] ))),
                   as.numeric(sapply(dsc.33, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, dsc.33[dsc.33!=x] ])) ))
  colnames(dif.dis.33) <- c("age_diffs","dists")
  
  dif.dis.conts.33[[ i ]] <- cor.test(dif.dis.33[,"age_diffs"], dif.dis.33[,"dists"])
  
}
saveRDS(dif.dis.conts.33, file = sprintf("%s/R_objects/dif.dis.conts.33.rds", p2_dir))

mean(unlist(lapply(dif.dis.conts.33, function(x) x$p.value)))
mean(unlist(lapply(dif.dis.conts.33, function(x) x$estimate)))







# running correlation for all controls used in subsampling
dsc <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
dsc <- dsc[ mTab.ds.h[dsc, "Downs_Syndrome"]=="No" ]
# dsc <- dsc[ mTab.ds.h[dsc, "Age"] < 34 ]
dif.dis <- cbind(as.numeric(sapply(dsc, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc[dsc!=x], "Age"] ))),
                 as.numeric(sapply(dsc, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, dsc[dsc!=x] ])) ))
colnames(dif.dis) <- c("age_diffs","dists")

plot(dif.dis)
cor.test(dif.dis[,"age_diffs"], dif.dis[,"dists"])






# running correlation for all controls used in subsampling, only those up to age 33, to have same trajectory as DS
dsc.33 <- unique(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(x) x$samples)))
dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Downs_Syndrome"]=="No" ]
dsc.33 <- dsc.33[ mTab.ds.h[dsc.33, "Age"] < 34 ]
dif.dis.33 <- cbind(as.numeric(sapply(dsc.33, function(x) abs(mTab.ds.h[x, "Age"]-mTab.ds.h[dsc.33[dsc.33!=x], "Age"] ))),
                    as.numeric(sapply(dsc.33, function(x) as.numeric(as.matrix(ordObj.ds.h$Aitchison)[x, dsc.33[dsc.33!=x] ])) ))
colnames(dif.dis.33) <- c("age_diffs","dists")
saveRDS(dif.dis.33, file = sprintf("%s/R_objects/dif.dis.33.rds", p2_dir))

plot(dif.dis.33)
cor.test(dif.dis.33[,"age_diffs"], dif.dis.33[,"dists"])





# ****************************************************************************************************************** #
# Scatterplot of healthyPop and DS age_diffs vs distances, combined ####

library(ggnewscale)

dif.dis.DS <- readRDS(file = sprintf("%s/R_objects/dif.dis.DS.rds", p2_dir))
dif.dis.healthyPop <- readRDS(file = sprintf("%s/R_objects/dif.dis.healthyPop.rds", p2_dir))

# combine into a single table
dda <- as.data.frame( rbind(dif.dis.healthyPop, dif.dis.DS) )
dda$samps <- c( rep("non-DS", nrow(dif.dis.healthyPop)), rep("DS", nrow(dif.dis.DS)) )

ggplot(dda, aes(x=age_diffs, y=dists, color=samps)) +
  geom_point(aes(shape = samps, color=samps), size=2.5) +
  scale_shape_manual(values=c(16,1), guide=F) +
  scale_color_manual(values=c("#CB6767","#73BDD3")) +
  labs(x="Age Difference", y="Aitchison Distance", color="") +
  new_scale("color") +
  geom_smooth(method = lm, aes(color=samps)) +
  scale_color_manual(values=c("#bf4342","#2f6690"), guide=F) +
  theme_classic() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=15),
        legend.title = element_text(size=17), legend.text = element_text(size=15))




# try also with only matched controls under 34
dif.dis.33 <- readRDS(file = sprintf("%s/R_objects/dif.dis.33.rds", p2_dir))

# combine into a single table
dd33 <- as.data.frame( rbind(dif.dis.33, dif.dis.DS) )
dd33$samps <- c( rep("non-DS", nrow(dif.dis.33)), rep("DS", nrow(dif.dis.DS)) )

ggplot(dd33, aes(x=age_diffs, y=dists)) +
  geom_point(aes(shape = samps, color=samps), size=2.5) +
  scale_shape_manual(values=c(16,1), guide=F) +
  scale_color_manual(values=c("#CB6767","#73BDD3")) +
  labs(x="Age Difference", y="Aitchison Distance", color="") +
  new_scale("color") +
  geom_smooth(method = lm, aes(color=samps)) +
  scale_color_manual(values=c("#bf4342","#2f6690"), guide=F) +
  theme_classic() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=15),
        legend.title = element_text(size=17), legend.text = element_text(size=15))





# combine plots for paper Figure 2:
library(ggpubr)
vsAgeGrp <- ggplot(ag.dists.m, aes(x=Distance, y=value, fill=Age_group)) +
  geom_boxplot(notch = T) + coord_flip() + ylim(0, 1.1) +
  guides(fill=guide_legend(title=NULL, reverse = T)) +
  stat_compare_means(aes(group=Age_group, label=..p.signif..), color="red", hide.ns=T, method="kruskal",
                     label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
  ggtitle("DS vs Age Groups") +
  theme_classic() + ylab(NULL) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size=15), legend.position = "bottom")

ageDif.DS <- ggplot(dd33, aes(x=age_diffs, y=dists)) +
  geom_point(aes(shape = samps, color=samps), size=2.5) +
  scale_shape_manual(values=c(16,1), guide=F) +
  scale_color_manual(values=c("#CB6767","#73BDD3")) +
  labs(x="Age Difference", y="Aitchison Distance", color="") +
  new_scale("color") +
  geom_smooth(method = lm, aes(color=samps)) +
  scale_color_manual(values=c("#bf4342","#2f6690"), guide=F) +
  theme_classic() +
  theme(axis.title = element_text(size=17), axis.text = element_text(size=15),
        legend.title = element_text(size=17), legend.text = element_text(size=15))

ggarrange(vsAgeGrp, ageDif.DS, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))



# or only healthyPop under 34
dif.dis.hp33 <- readRDS(file = sprintf("%s/R_objects/dif.dis.hp33.rds", p2_dir))





# ****************************************************************************************************************** #








# ****************************************************************************************************************** #

# Make bins of ages by either 5 or 10 years ####
mTab.ds.h$Age_bins <- sapply(mTab.ds.h$Age, function(x)
  ifelse(x < 10, "0_10",
         ifelse(x < 15, "10_15",
                ifelse(x < 20, "15_20",
                       ifelse(x < 30, "20_30",
                              ifelse(x < 40, "30_40",
                                     ifelse(x < 50, "40_50",
                                            ifelse(x < 60, "50_60",
                                                   ifelse(x < 70, "60_70", "70_+")))))))) )

# ******************************** #
get_dists_to_ageBins <- function(ordObj, dist_meas, mTab, dss, hp) {
  
  # Aitchison distances on scale of 0-100, others 0-1, so will have to adjust Aitchison to be comparable
  mult <- 1
  if (dist_meas=="Aitchison")
    mult <- 100
  
  b.0_10  <- hp[ mTab[ hp, "Age_bins"] == "0_10" ]
  b.10_15 <- hp[ mTab[ hp, "Age_bins"] == "10_15" ]
  b.15_20 <- hp[ mTab[ hp, "Age_bins"] == "15_20" ]
  b.20_30 <- hp[ mTab[ hp, "Age_bins"] == "20_30" ]
  b.30_40 <- hp[ mTab[ hp, "Age_bins"] == "30_40" ]
  b.40_50 <- hp[ mTab[ hp, "Age_bins"] == "40_50" ]
  b.50_60 <- hp[ mTab[ hp, "Age_bins"] == "50_60" ]
  b.60_70 <- hp[ mTab[ hp, "Age_bins"] == "60_70" ]
  b.70    <- hp[ mTab[ hp, "Age_bins"] == "70_+" ]
  
  return( list( "0_10"  = as.numeric(ordObj[dss, b.0_10]) / mult, 
                "10_15" = as.numeric(ordObj[dss, b.10_15]) / mult, 
                "15_20" = as.numeric(ordObj[dss, b.15_20]) / mult, 
                "20_30" = as.numeric(ordObj[dss, b.20_30]) / mult,
                "30_40" = as.numeric(ordObj[dss, b.30_40]) / mult,
                "40_50" = as.numeric(ordObj[dss, b.40_50]) / mult,
                "50_60" = as.numeric(ordObj[dss, b.50_60]) / mult,
                "60_70" = as.numeric(ordObj[dss, b.60_70]) / mult,
                "70_+"   = as.numeric(ordObj[dss, b.70]) / mult ))
}
# ******************************** #


ab.dists <- list()
ab.dists[[ "Aitchison" ]]          <- get_dists_to_ageBins(as.matrix(ordObj.ds.h$Aitchison), "Aitchison", mTab.ds.h, DS.samps, healthyPop)
ab.dists[[ "Weighted_Unifrac" ]]   <- get_dists_to_ageBins(as.matrix(ordObj.ds.h$Weighted_Unifrac), "Weighted_Unifrac", mTab.ds.h, DS.samps, healthyPop)
ab.dists[[ "Unweighted_Unifrac" ]] <- get_dists_to_ageBins(as.matrix(ordObj.ds.h$Unweighted_Unifrac), "Unweighted_Unifrac", mTab.ds.h, DS.samps, healthyPop)
ab.dists[[ "Bray" ]]               <- get_dists_to_ageBins(as.matrix(ordObj.ds.h$Bray), "Bray", mTab.ds.h, DS.samps, healthyPop)
ab.dists[[ "Jaccard" ]]            <- get_dists_to_ageBins(as.matrix(ordObj.ds.h$Jaccard), "Jaccard", mTab.ds.h, DS.samps, healthyPop)

ab.dists.m <- reshape2::melt(ab.dists)
colnames(ab.dists.m) <- c("value","Age_bin","Distance")
ab.dists.m$Age_bin <- factor(ab.dists.m$Age_bin, levels=rev(c("0_10","10_15","15_20","20_30","30_40",
                                                              "40_50","50_60","60_70","70_+")))
ab.dists.m$Distance <- factor(ab.dists.m$Distance, levels=rev(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac",
                                                                "Bray","Jaccard")))

library(ggpubr)
ggplot(ab.dists.m, aes(x=Distance, y=value, fill=Age_bin)) +
  geom_boxplot(notch = T) + coord_flip() + ylim(0, 1.1) +
  guides(fill=guide_legend(title=NULL, reverse = T)) +
  stat_compare_means(aes(group=Age_bin, label=..p.signif..), color="red", hide.ns=T, method="kruskal",
                     label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
  ggtitle("(c) DS vs Age Bins") +
  theme_classic() + ylab(NULL) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
        plot.title = element_text(size = 18),
        legend.text = element_text(size=15), legend.position = "bottom")



tukey_kruskal("Aitchison", ab.dists, sigOnly = T)
tukey_kruskal("Weighted_Unifrac", ab.dists, sigOnly = T)
tukey_kruskal("Unweighted_Unifrac", ab.dists, sigOnly = T)
tukey_kruskal("Bray", ab.dists, sigOnly = T)
tukey_kruskal("Jaccard", ab.dists, sigOnly = T)

# ******************************** #






ab.anosims <- list()

ab.anosims[[ "DS" ]] <- anosim(as.dist( as.matrix(ordObj.ds.h$Aitchison)[DS.samps, DS.samps] ),
       as.character(as.matrix(mTab.ds.h[DS.samps, "Age_bins"])) )

ab.anosims[[ "healthPop" ]] <- anosim(as.dist( as.matrix(ordObj.ds.h$Aitchison)[healthyPop, healthyPop] ),
       as.character(as.matrix(mTab.ds.h[healthyPop, "Age_bins"])) )


# ****************************** #
plot_age_anosims.DS <- function(sam_tab, dist_obj, dist_meas, comp.method, age_var, fam.anosim=NULL, fam.adonis=NULL) {
  
  bins <- sort(unique(sam_tab[ , age_var]))
  # Aitchison distances on scale of 0-100, others 0-1, so will have to adjust Aitchison to be comparable
  mult <- 1
  if (dist_meas=="Aitchison")
    mult <- 100
  
  all_bin_dists <- list()
  
  for (bin in bins) {
    bin.samps <- rownames(sam_tab[ sam_tab[ , age_var]==bin, ])
    
    all_bin_dists[[ bin ]] <- sapply(bins, function(x) {
      x.samps <- rownames(sam_tab[ sam_tab[ , age_var]==x, ])
      as.numeric(dist_obj[bin.samps, x.samps]) / mult
    })
    names(all_bin_dists[[ bin ]]) <- sprintf("vs_%s", bins)
  }
  
  abd.m <- melt(all_bin_dists)
  colnames(abd.m) <- c("value","vs","bin")
  # remove rows with value==0 since those are just sample vs itself
  abd.m <- abd.m[ abd.m$value > 0, ]
  abd.m$vs <- factor(abd.m$vs, levels=rev(sprintf("vs_%s", bins)))
  abd.m$bin <- factor(abd.m$bin, levels=rev(bins))
  
  ggplot(abd.m, aes(x=bin, y=value, fill=vs)) +
    geom_boxplot(notch = T) + coord_flip() +# ylim(0, 1.1) +
    guides(fill=guide_legend(title=NULL, reverse = T)) +
    stat_compare_means(aes(group=vs, label=..p.signif..), hide.ns=T, method="anova") +#,
                       # label.x = c(5, 4, 3, 2, 1), label.y = c(0.95, 1.05, 0.8, 0.58, 0.7)) +
    theme(axis.text = element_text(size=15), axis.title = element_text(size=17), 
          legend.text = element_text(size=15), legend.position = "bottom")

}
# ****************************** #

plot_age_anosims.DS(mTab.ds.h[ DS.samps, ], as.matrix(ordObj.ds.h$Aitchison), 
                    "Aitchison", "ANOSIM", "Age_groups", fam.anosim = ab.anosims$DS)

plot_age_anosims.DS(mTab.ds.h[ DS.samps, ], as.matrix(ordObj.ds.h$Aitchison), 
                    "Aitchison", "ANOSIM", "Age_bins", fam.anosim = ab.anosims$DS)

# ******************************** #




gen.num_sig.age.ds <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$Genus), function(x)
    sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Anova$Genus[gsub("-","\\.",x), "Age"]
    }
    )) < 0.05)
  )
gen.num_sig.age.ds.cov <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$Genus), function(x)
    sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
      # *************** fix this here after using noRemove=T in water_anovas ************* 
      y$Anova$Genus[gsub("-","\\.",x), "Downs_Syndrome"]
    }
    )) < 0.05)
  )

# 10 most frequently significant by age, how often were those significant by DS:
gen.num_sig.age.ds.cov[ tail(names(sort(gen.num_sig.age.ds)), n=10) ]


plot_data.cont("Age", "Alloprevotella", "TvsQ", "Genus", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")
plot_data.cont("Age", "Peptostreptococcus", "TvsQ", "Genus", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")
plot_data.cont("Age", "Haemophilus", "TvsQ", "Genus", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")





phy.num_sig.age.ds <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$Phylum), function(x)
  sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    # *************** fix this here after using noRemove=T in water_anovas ************* 
    y$Anova$Phylum[gsub("-","\\.",x), "Age"]
  }
  )) < 0.05)
)
phy.num_sig.age.ds.cov <- sapply(rownames(disSubsTests$default$Downs_Syndrome$`1`$Anova$Phylum), function(x)
  sum(unlist(lapply(disSubsTests$default$Downs_Syndrome, function(y) {
    # *************** fix this here after using noRemove=T in water_anovas ************* 
    y$Anova$Phylum[gsub("-","\\.",x), "Downs_Syndrome"]
  }
  )) < 0.05)
)

# 10 most frequently significant by age, how often were those significant by DS:
sort(phy.num_sig.age.ds)
phy.num_sig.age.ds.cov[ names(sort(phy.num_sig.age.ds)) ]



plot_data.cont("Age", "Proteobacteria", "TvsQ", "Phylum", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")
plot_data.cont("Age", "Bacteroidetes", "TvsQ", "Phylum", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")
plot_data.cont("Age", "unclassified.P1", "TvsQ", "Phylum", c(dsc.33,DS.samps), only_cont, gloms_clr, SLL2.meta, facetVar = "Downs_Syndrome")

# ****************************************************************************************************************** #






# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####
# Test effects of stomatotypes within disorder+control groups ####
stom_chis <- list()
for ( disorder in names(disSubsTests$default)) {
  stom_chis[[ disorder ]] <- list()
  
  for (i in names(disSubsTests$default[[ disorder ]])) {
    
    stom_chis[[ disorder ]][[ i ]] <- list()
    mTab.sub <- disSubsTests$default[[ disorder ]][[ i ]]$mTab
    
    for (dist_meas in c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard")) {
      stom_chis[[ disorder ]][[ i ]][[ dist_meas ]] <- chisq.test(mTab.sub[, disorder], mTab.sub[ , sprintf("Stomatotype_%s", dist_meas)])
    }
    
  }
}

stom_chis.pMeans.cf <- sapply(names(stom_chis$Cystic_fibrosis$`1`), function(x)
  mean(p.adjust(unlist(lapply(stom_chis$Cystic_fibrosis, function(y) y[[ x ]]$p.value)), method = "fdr") )
)

stom_chis.num_sig.cf <- sapply(names(stom_chis$Cystic_fibrosis$`1`), function(x)
  sum(p.adjust(unlist(lapply(stom_chis$Cystic_fibrosis, function(y) y[[ x ]]$p.value)), method = "fdr") < 0.05)
)



stom_chis.pMeans.ds <- sapply(names(stom_chis$Downs_Syndrome$`1`), function(x)
  mean(p.adjust(unlist(lapply(stom_chis$Downs_Syndrome, function(y) y[[ x ]]$p.value)), method = "fdr") )
)

stom_chis.num_sig.ds <- sapply(names(stom_chis$Downs_Syndrome$`1`), function(x)
  sum(p.adjust(unlist(lapply(stom_chis$Downs_Syndrome, function(y) y[[ x ]]$p.value)), method = "fdr") < 0.05)
)


# ****************************************************************************************************************** #
# ****************************************************************************************************************** #

# Test effects of water vals within disorder+control groups ####

dv.subs <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
             "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac",
             "Stomatotype_Unweighted_Unifrac","Stomatotype_Bray","Stomatotype_Jaccard")

# # water_anovas <- list()
# 
# for ( disorder in names(disSubsTests$default)) {
# 
#   # water_anovas[[ disorder ]] <- list()
#   fix.subs <- c(disorder, "Age", "Gender")
# 
#   for (i in names(disSubsTests$default[[ disorder ]])) {
#     print(sprintf("%s  %s", disorder, i))
# 
#     # water_anovas[[ disorder ]][[ i ]] <- list()
#     mTab.sub <- disSubsTests$default[[ disorder ]][[ i ]]$mTab
#     samps.sub <- disSubsTests$default[[ disorder ]][[ i ]]$samples
# 
#     for (tl in c("contVar","Phylum","Class","Order","Family","Genus")) {
#       # water_anovas[[ disorder ]][[ i ]][[ tl ]] <- list()
#       # ******************** #
#       if (tl == "contVar") {
#         glomTab <- NULL
#         dv <- dv.subs
#       } else {
#         glomTab <- gloms_clr[[ tl ]][ , samps.sub]
#         dv <- NULL
#       }
#       # ******************** #
# 
#       for (water_val in c(cont_water_data, "Age")) {
#         water_anovas[[ disorder ]][[ i ]][[ tl ]][[ water_val ]] <- get_lm( c(water_val, fix.subs[fix.subs!=water_val]), tl, mTab.sub,
#                                                                             glomTab, dv, noRemove = T, rerun.nonSig = F )
#       }
#       # ******************** #
#     }
#   }
# }
# saveRDS(water_anovas, sprintf("%s/R_objects/water_anovas.rds", p2_dir))
# # ******************************************************************** #
water_anovas <- readRDS(sprintf("%s/R_objects/water_anovas.rds", p2_dir))
# ******************************************************************** #


# ******************************************************************** #
# vars to add to DV ####
more_vars <- c("Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
               "Stomatotype_Bray","Stomatotype_Jaccard",
               "MALDI.Yeast_detected","MALDI.Mold_detected","MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies",
               "Allergy","pH","BMI",
               "MALDI.Yeast.Candida_albicans","MALDI.Yeast.Candida_albicans_AND_Candida_dubliniensis",
               "MALDI.Yeast.Candida_albicans_AND_Candida_glabrata","MALDI.Yeast.Candida_albicans_AND_Candida_guillermondii",
               "MALDI.Yeast.Candida_albicans_AND_Candida_intermedia","MALDI.Yeast.Candida_albicans_AND_Candida_parapsilosis",
               "MALDI.Yeast.Candida_albicans_AND_Candida_parapsilosis_AND_Trichosporon_spp.",
               "MALDI.Yeast.Candida_albicans_AND_Candida_tropicalis","MALDI.Yeast.Candida_albicans_AND_Debaryomyces_hansenii",
               "MALDI.Yeast.Candida_boidinii","MALDI.Yeast.Candida_dubliniensis",
               "MALDI.Yeast.Candida_glabrata","MALDI.Yeast.Candida_glabrata_AND_Candida_albicans",
               "MALDI.Yeast.Candida_glabrata_AND_Candida_krusei","MALDI.Yeast.Candida_guillermondii",
               "MALDI.Yeast.Candida_guillermondii_AND_Candida_albicans","MALDI.Yeast.Candida_intermedia",
               "MALDI.Yeast.Candida_intermedia_AND_Candida_dubliniensis_AND_Candida_albicans","MALDI.Yeast.Candida_kefyr",
               "MALDI.Yeast.Candida_krusei","MALDI.Yeast.Candida_lusitaniae",
               "MALDI.Yeast.Candida_lusitaniae_AND_Candida_parapsilosis","MALDI.Yeast.Candida_parapsilosis",
               "MALDI.Yeast.Candida_parapsilosis_AND_Candida_albicans","MALDI.Yeast.Candida_parapsilosis_AND_Candida_tropicalis",
               "MALDI.Yeast.Candida_spp.","MALDI.Yeast.Candida_tropicalis","MALDI.Yeast.Candida_zeylanoides",
               "MALDI.Yeast.Cryptococcus_spp.","MALDI.Yeast.Debaryomyces_hansenii",
               "MALDI.Yeast.Debaryomyces_hansenii_AND_Candida_zeylanoides","MALDI.Yeast.Kloeckera_apiculata",
               "MALDI.Yeast.Kloeckera_apiculata_AND_Candida_intermedia_AND_Candida_boidinii",
               "MALDI.Yeast.Rhodotorula_mucilaginosa","MALDI.Yeast.Rhodotorula_mucilaginosa_AND_Candida_parapsilosis",
               "MALDI.Yeast.Rhodotorula_mucilaginosa_AND_Cryptococcus_spp.","MALDI.Yeast.Sarocladium_strictum",
               "MALDI.Yeast.Trichosporon_spp.",
               cont_water_data)

# first must update the MALDI values to reflect presence of particular yeast properly:
yeast_of_interest <- c("Candida_albicans","Candida_dubliniensis","Candida_glabrata","Candida_parapsilosis",
                       #"Candida_boidinii", # too few - no need to run any tests
                       "Candida_guillermondii",
                       #"Candida_tropicalis",
                       "Candida_intermedia","Candida_krusei","Candida_lusitaniae",
                       #"Candida_kefyr","Candida_zeylanoides","Trichosporon_spp",
                       "Cryptococcus_spp","Debaryomyces_hansenii",
                       #"Kloeckera_apiculata",,"Sarocladium_strictum"
                       "Rhodotorula_mucilaginosa")

for (yeast in yeast_of_interest) {
  
  yVars <- more_vars[ grepl(yeast, more_vars) ]
  
  SLL2.meta[ , sprintf("Full_MALDI.%s", yeast) ] <- sapply(rownames(SLL2.meta), function(x) 
    ifelse( sum(SLL2.meta[x, yVars]=="Yes") > 0, "Yes", "No" ))
  
  mTab.ds.h[ , sprintf("Full_MALDI.%s", yeast) ] <- sapply(rownames(mTab.ds.h), function(x) 
    ifelse( sum(mTab.ds.h[x, yVars]=="Yes") > 0, "Yes", "No" ))
  
  # print(yeast)
  # print(table(SLL2.meta[ , sprintf("Full_MALDI.%s", yeast) ]))
  # print(table(mTab.ds.h[ , sprintf("Full_MALDI.%s", yeast) ]))
  # print("")
  
}


new.dv <- c("MALDI.Yeast_detected","MALDI.Mold_detected","MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies",
            sprintf("Full_MALDI.%s", yeast_of_interest),
            "Allergy","pH","BMI")


# more_anovas <- list()
# 
# for ( disorder in names(disSubsTests$default)) {
# 
#   # ******************************** #
#   if (disorder == "Cystic_fibrosis") {
#     cov_vars <- c("Antibiotics","Gender","Age","Population")#
#   } else {
#     cov_vars <- c("Gender","Age","Population")#
#   }
#   # ******************************** #
# 
#   more_anovas[[ disorder ]] <- list()
# 
#   for (i in names(disSubsTests$default[[ disorder ]])) {
#     print(sprintf("%s  %s", disorder, i))
# 
#     more_anovas[[ disorder ]][[ i ]] <- list()
# 
#     samps.sub <- disSubsTests$default[[ disorder ]][[ i ]]$samples
#     mTab.sub <- SLL2.meta[ samps.sub, ]
# 
#     for (tl in c("contVar")) {
# 
#       # ******************** #
#       if (tl == "contVar") {
#         glomTab <- NULL
#         dv <- new.dv
#       } else {
#         glomTab <- gloms_clr[[ tl ]][ , samps.sub]
#         dv <- NULL
#       }
#       # ******************** #
# 
#       more_anovas[[ disorder ]][[ i ]][[ tl ]] <- get_lm( c(disorder, cov_vars), tl, mTab.sub, glomTab, dv, 
#                                                           noRemove = T, rerun.nonSig = F )
# 
#       # ******************** #
#     }
#   }
# }
# saveRDS(more_anovas, sprintf("%s/R_objects/more_anovas.rds", p2_dir))
# # ******************************************************************** #
more_anovas <- readRDS(sprintf("%s/R_objects/more_anovas.rds", p2_dir))
# ******************************************************************** #



more_anovas.pMeans.cf <- sapply(new.dv, function(x)
  mean(unlist(lapply(more_anovas$Cystic_fibrosis, function(y) {
    y$contVar[x,"Cystic_fibrosis"]
  }
  )))
)

more_anovas.pMeans.ds <- sapply(new.dv, function(x)
  mean(unlist(lapply(more_anovas$Downs_Syndrome, function(y) {
    y$contVar[x,"Downs_Syndrome"]
  }
  )))
)



more_anovas.num_sig.cf <- sapply(new.dv, function(x)
  sum(unlist(lapply(more_anovas$Cystic_fibrosis, function(y) {
    y$contVar[x,"Cystic_fibrosis"]
  }
  )) < 0.05)
)

more_anovas.num_sig.ds <- sapply(new.dv, function(x)
  sum(unlist(lapply(more_anovas$Downs_Syndrome, function(y) {
    y$contVar[x,"Downs_Syndrome"]
  }
  )) < 0.05)
)
# ******************************************************************** #






# vars to add to fixed effects against tl and contVar
more_fixeds <- c("pH", # if signif for something interesting by DS vs H
                 "Consumption.Sweets","Consumption.Nuts")
# ****************************************************************************************************************** #
# ****************************************************************************************************************** #






































# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# Rerun tests for CF excluding Antibiotics as a covar, but use same matched control subsamplings ####
no_ab_subsTests <- list()

for (disorder in c("Cystic_fibrosis")) {
  
  # take disorder samples and matched controls
  matchCont <- lapply(disSubsTests$default[[ disorder ]], function(x) x$samples)
  
  # run tests using these same matched controls with the fam members now instead of disorder
  no_ab_subsTests[[ disorder ]] <- run_full_subsampling_calcs.disorder(disorder, SLL2, SLL2.meta, 100, sampMode = "useAntibiotics", 
                                                                       chosenControls = matchCont, useAntibiotics = F)
  
}

saveRDS(no_ab_subsTests, file = sprintf("%s/R_objects/no_ab_subsTests.rds", p2_dir))
# ************************************************************************ #

no_ab_subsTests <- readRDS(sprintf("%s/R_objects/no_ab_subsTests.rds", p2_dir))
# ************************************************************************ #



gen.pMeans.cf.no_ab <- sapply(rownames(gloms_clr$Genus), function(x)
  mean(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)
phy.pMeans.cf.no_ab <- sapply(rownames(gloms_clr$Phylum), function(x)
  mean(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)
conts.pMeans.cf.no_ab <- sapply(rownames(no_ab_subsTests$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  mean(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )))
)


gen.num_sig.cf.no_ab <- sapply(rownames(gloms_clr$Genus), function(x)
  sum(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$Genus[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)
phy.num_sig.cf.no_ab <- sapply(rownames(gloms_clr$Phylum), function(x)
  sum(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$Phylum[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)
conts.num_sig.cf.no_ab <- sapply(rownames(no_ab_subsTests$Cystic_fibrosis$`1`$Anova$contVar), function(x)
  sum(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Anova$contVar[gsub("-","\\.",x),"Cystic_fibrosis"]
  }
  )) < 0.05)
)




adonis.pMeans.cf.no_ab <- sapply(names(no_ab_subsTests$Cystic_fibrosis$`1`$Adonis), function(x)
  mean(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Adonis[[ x ]][ "Cystic_fibrosis", "Pr(>F)"]
  }
  )))
)

adonis.num_sig.cf.no_ab <- sapply(names(no_ab_subsTests$Cystic_fibrosis$`1`$Adonis), function(x)
  sum(unlist(lapply(no_ab_subsTests$Cystic_fibrosis, function(y) {
    y$Adonis[[ x ]][ "Cystic_fibrosis", "Pr(>F)"]
  }
  )) < 0.05)
)






# ****************************************************************************************************************** #










# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# Effects of antibiotic use within CF samples only ####

CF.samps <- rownames(SLL2.meta[ SLL2.meta$Cystic_fibrosis=="Yes", ])

# first must update the MALDI values to reflect presence of particular yeast properly:
yeast_of_interest <- c("Candida_albicans","Candida_dubliniensis","Candida_glabrata","Candida_parapsilosis",
                       #"Candida_boidinii", # too few - no need to run any tests
                       "Candida_guillermondii",
                       #"Candida_tropicalis",
                       "Candida_intermedia","Candida_krusei","Candida_lusitaniae",
                       #"Candida_kefyr","Candida_zeylanoides","Trichosporon_spp",
                       "Cryptococcus_spp","Debaryomyces_hansenii",
                       #"Kloeckera_apiculata",,"Sarocladium_strictum"
                       "Rhodotorula_mucilaginosa")

new.dv <- c("MALDI.Yeast_detected","MALDI.Mold_detected","MALDI.Num_Yeast_Colonies","MALDI.Num_Mold_Colonies",
            sprintf("Full_MALDI.%s", yeast_of_interest),
            "Allergy","pH","BMI","Height","Weight")

# make a new variable to test - pancreatic problems
mtab.cf <- SLL2.meta[CF.samps, ]
panc.liver <- sample_names(subset_samples(SLL2, grepl('pancrea',Other_medications,ignore.case=T) | 
                                            grepl('kreon',Other_medications,ignore.case=T) |
                                            grepl('urso',Other_medications,ignore.case=T) | 
                                            grepl('pancrea',Other_disorder,ignore.case=T) ))
panc.liver <- panc.liver[ panc.liver %in% CF.samps]
mtab.cf$Panc.Liver <- sapply(CF.samps, function(x) ifelse(x %in% panc.liver, 'Yes', 'No'))



ab_in_cf <- list()
diab_in_cf <- list()
pl_in_cf <- list()
transplant_in_cf <- list()
ca_in_cf <- list()

for (tl in c("contVar","moreVars","Phylum","Class","Order","Family","Genus")) {
  
  print(tl)
  
  # ******************** #
  if (tl == "contVar") {
    glomTab <- NULL
    tl.general <- tl
    dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
            "Stomatotype_Aitchison","Stomatotype_Weighted_Unifrac","Stomatotype_Unweighted_Unifrac",
            "Stomatotype_Bray","Stomatotype_Jaccard","pH","BMI")
    # ******************** #
  } else if (tl == "moreVars") {
    glomTab <- NULL
    tl.general <- "contVar"
    dv <- new.dv
    # ******************** #
  } else {
    glomTab <- gloms_clr
    tl.general <- tl
    dv <- NULL
  }
  # ******************** #
  
  ab_in_cf[[ tl ]] <- get_lm( c("Antibiotics","Gender","Age","Population"), 
                              tl.general, mtab.cf, glomTab[[ tl.general ]][ , CF.samps ], 
                              dv, noRemove = T, rerun.nonSig = F)
  
  # if (tl == "moreVars")
  #   next
  diab_in_cf[[ tl ]] <- get_lm( c("Diabetes","Gender","Age","Population"), 
                              tl, mtab.cf, glomTab[[ tl.general ]][ , CF.samps ], 
                              dv, noRemove = T, rerun.nonSig = F)
  
  pl_in_cf[[ tl ]] <- get_lm( c("Panc.Liver","Gender","Age","Population"), 
                              tl.general, mtab.cf, glomTab[[ tl.general ]][ , CF.samps ], 
                              dv, noRemove = T, rerun.nonSig = F)
  
  # if (tl == "moreVars")
  #   next
  transplant_in_cf[[ tl ]] <- get_lm( c("Transplant","Gender","Age","Population"), 
                                      tl, mtab.cf, glomTab[[ tl.general ]][ , CF.samps ], 
                                      dv, noRemove = T, rerun.nonSig = F)
  
  ca_in_cf[[ tl ]] <- get_lm( c("Full_MALDI.Candida_albicans","Gender","Age","Population"), 
                              tl.general, mtab.cf, glomTab[[ tl.general ]][ , CF.samps ], 
                              dv, noRemove = T, rerun.nonSig = F)
  
}



# ****************************************************************************************************************** #






# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####




# Table for github ####
s.anon <- SLL2.meta
rownames(s.anon) <- paste0("S_",1:nrow(s.anon))

toRemove_s.anon <- c("Birthdate", 
                     "Country_of_birth","Province_of_birth","City_of_birth","Country_of_birth.mother","Province_of_birth.mother",
                     "City_of_birth.mother","Country_of_birth.father","Province_of_birth.father","City_of_birth.father", 
                     "Education.mother","Education.father","Education",
                     "Occupation.mother","Occupation.mother_other","Occupation.father","Occupation.father_other",
                     "Postal_code", "Want_info","Email",
                     "FQ.Mutation","FQ.Lungs.Pseudomonas_aeruginosa","FQ.Lungs.Staphylococcus_aureus","FQ.Lungs.Haemophilus_sp",
                     "FQ.Lungs.Others","FQ.Lungs.Fungi_binary","FQ.Lungs.Fungi","FQ.Infection.Bronchopneumonia","FQ.Infection.Pneumonia",
                     "FQ.Infection.Bronchiolitis","FQ.Transplant.Lung_transplant","FQ.Transplant.Immunosuppressant",
                     "FQ.Transplant.Other_binary","FQ.Transplant.Other","FQ.Diet.Specialized","FQ.Diet.Reduced_carbs",
                     "FQ.Diet.More_fats","FQ.Diet.More_salt","FQ.Diet.Isotonic_drinks","FQ.Diet.Isotonic_drinks_frequency",
                     "FQ.Diet.Other_binary","FQ.Diet.Other","FQ.Medications.Corticoids","FQ.Medications.Pancreatic_enzymes",
                     "FQ.Medications.Insulin","FQ.Medications.Bronchodilators","FQ.Medications.Immunosuppressants",
                     "FQ.Medications.Others_binary","FQ.Medications.Others","FQ.Med_format.Oral","FQ.Med_format.Intravenous",
                     "FQ.Med_format.Inhaled","FQ.Hygiene.Wash_hands_often","FQ.Hygiene.Dont_share_bottles",
                     "FQ.Hygiene.Avoid_contact_other_FQ","FQ.Hygiene.Other_binary","FQ.Hygiene.Other",
                     "Weighted_Unifrac","Unweighted_Unifrac","VAW_GUnifrac","a0_GUnifrac","a05_GUnifrac","JSD","Bray.Curtis",
                     "Jaccard","Canberra","Aitchison")
s.anon <- s.anon[ , ! colnames(s.anon) %in% toRemove_s.anon ]

write.csv(s.anon, sprintf("%s/SLL2.metadata.csv", p2_dir))





# ****************************************************************************************************************** #



# ****************************************************************************************************************** #
###### Prepare table of sample attributes for Sequence Reads Archive ######
# ****************************************************************************************************************** #


attribs <- c("sample_name","bioproject_accession","organism","ecotype","collection_date","env_biome","env_feature",
             "env_material","geo_loc_name","host","isol_growth_condt","lat_lon","chem_administration",
             "host_age","host_body_mass_index","host_body_product","host_height","host_sex",
             "host_tissue_sampled","host_tot_mass","isolation_source","samp_collect_device","samp_store_temp")

# collection months
jan <- c("01","02","03","04")
feb <- c("06","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26")
mar <- c("27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50")
apr <- c("51","52","53","54","55","57","58","59","60")


# update rownames to anonymize 
SLL2.SRA_meta <- SLL2.meta
rownames(SLL2.SRA_meta) <- paste0("S_",1:nrow(SLL2.SRA_meta))

# for particular columns
attrs.date <- sapply(rownames(SLL2.meta), function(x) {
  sID <- strsplit(x,"\\.")[[1]][2]
  ifelse(sID %in% jan, "Jan-2017",
         ifelse(sID %in% feb, "Feb-2017",
                ifelse(sID %in% mar, "Mar-2017",
                       ifelse(sID %in% apr, "Apr-2017", "No_date"))))
})
names(attrs.date) <- paste0("S_", 1:length(attrs.date))

# names for the 2 separate tables
s.1.824    <- rownames(SLL2.SRA_meta)[1:824]
s.825.1648 <- rownames(SLL2.SRA_meta)[825:1648]

samp_lat_lon <- function(x) {
  lat <- sprintf("%s N", SLL2.SRA_meta[x, "Latitude"])
  lon <- sprintf("%s %s", abs(SLL2.SRA_meta[x, "Longitude"]),
                 ifelse(SLL2.SRA_meta[x, "Longitude"] < 0, "W", "E"))
  return(sprintf("%s, %s", lat, lon))
}

# since the age was somewhat estimated for some samples, and S_517 and S_518 now have all exactly the same values in this table
#  must adjust the age slightly for one so that the table can be accepted by SRA
SLL2.SRA_meta["S_517", "Age"]
SLL2.SRA_meta["S_962", "Age"]  <- 15.45
SLL2.SRA_meta["S_1222", "Age"] <- 14.45
SLL2.SRA_meta["S_1502", "Age"] <- 15.45
SLL2.SRA_meta["S_1580", "Age"] <- 14.45


# make tables
for (lab in c("s1-s824","s825-s1648")) {
  if (lab=="s1-s824") {
    samp_vec <- s.1.824
  } else {
    samp_vec <- s.825.1648
  }
  
  attr.table <- matrix(NA, nrow = length(samp_vec), ncol = length(attribs))
  rownames(attr.table) <- samp_vec
  colnames(attr.table) <- attribs
  
  for (samp in samp_vec) {
    
    attr.table[ samp, "sample_name" ] <- samp
    attr.table[ samp, "bioproject_accession" ] <- "PRJNA427101"
    attr.table[ samp, "organism" ] <- "Human oral microbiome"
    attr.table[ samp, "ecotype" ] <- "Spain"
    
    attr.table[ samp, "collection_date" ] <- attrs.date[ samp ]
    attr.table[ samp, "env_biome" ] <- "Human saliva"
    attr.table[ samp, "env_feature" ] <- "Human saliva"
    attr.table[ samp, "env_material" ] <- "Human saliva"
    attr.table[ samp, "geo_loc_name" ] <- sprintf("Spain: %s", 
                                                  gsub("á","a", 
                                                       gsub("é","e", 
                                                            gsub("í","i", 
                                                                 gsub("ú","u", SLL2.SRA_meta[samp, "City"])))) )
    attr.table[ samp, "host" ] <- "Homo sapiens"
    attr.table[ samp, "isol_growth_condt" ] <- "not applicable"
    attr.table[ samp, "lat_lon" ] <- samp_lat_lon(samp)
    
    # attr.table[ samp, "biotic_relationship" ] <- "from host"
    attr.table[ samp, "chem_administration" ] <- "Phosphate buffered saline (PBS)"
    attr.table[ samp, "host_age" ] <- ifelse( is.na(SLL2.SRA_meta[samp, "Age"]), 
                                              NA, 
                                              sprintf("%s years", round(SLL2.SRA_meta[samp, "Age"], 3)) )
    attr.table[ samp, "host_body_mass_index" ] <- ifelse( is.na(SLL2.SRA_meta[samp, "BMI"]), 
                                                          NA, 
                                                          round(SLL2.SRA_meta[samp, "BMI"], 2) )
    attr.table[ samp, "host_body_product" ] <- "Saliva"
    attr.table[ samp, "host_height" ] <- ifelse( is.na(SLL2.SRA_meta[samp, "Height"]), 
                                                 NA, 
                                                 sprintf("%s cm", SLL2.SRA_meta[samp, "Height"]) )
    attr.table[ samp, "host_sex" ] <- ifelse( is.na(SLL2.SRA_meta[samp, "Gender"]), 
                                              NA, 
                                              ifelse(SLL2.SRA_meta[samp, "Gender"] == "F", "Female","Male"))
    attr.table[ samp, "host_tissue_sampled" ] <- "Saliva"
    attr.table[ samp, "host_tot_mass" ] <- ifelse( is.na(SLL2.SRA_meta[samp, "Weight"]), 
                                                   NA, 
                                                   sprintf("%s kg", SLL2.SRA_meta[samp, "Weight"]) )
    attr.table[ samp, "isolation_source" ] <- "Saliva"
    attr.table[ samp, "samp_collect_device" ] <- "Phosphate buffered saline (PBS)"
    attr.table[ samp, "samp_store_temp" ] <- "-20 C"
  }
  
  write.csv(attr.table, sprintf('%s/SRA_submission/Biosamp.attributes.%s.csv', p2_dir, lab), row.names = F)
}











# ****************************************************************************************************************** #
###### Prepare table of metadata for files for Sequence Reads Archive ######
# ****************************************************************************************************************** #
# meta_labs <- c("bioproject_accession","biosample_accession","library_ID","title",
#                "library_strategy","library_source","library_selection",
#                "library_layout","platform","instrument_model","design_description",
#                "filetype","filename","filename2")

meta_labs <- c("sample_name","library_ID","title",
               "library_strategy","library_source","library_selection",
               "library_layout","platform","instrument_model","design_description",
               "filetype","filename","filename2")

# ****** #
# read sample accession numbers from the files provided by SRA after submitting the attributes info above
s.1.824    <- rownames(SLL2.SRA_meta)[1:824]
s.825.1648 <- rownames(SLL2.SRA_meta)[825:1648]

# get file with original sample names and anonymized names
s824_names <- rownames(SLL2.meta)[1:824]
names(s824_names) <- s.1.824
write.table(cbind(s824_names, names(s824_names) ), 
            sprintf("%s/SRA_submission/sample_names.s.1.824.csv", p2_dir), 
            row.names = F, col.names = F, sep = ",", quote = F)

s1648_names <- rownames(SLL2.meta)[825:1648]
names(s1648_names) <- s.825.1648
write.table(cbind(s1648_names, names(s1648_names) ), 
            sprintf("%s/SRA_submission/sample_names.s.825.1648.csv", p2_dir), 
            row.names = F, col.names = F, sep = ",", quote = F)


# samp_acc_s1_s700 <- read.delim(sprintf("%s/SRA_submission/Sample_accessions_s1-s700.tsv", p2_dir))
# samp_acc_s701_s1319 <- read.delim(sprintf("%s/SRA_submission/Sample_accessions_s701-s1319.tsv", p2_dir))
# 
# get_samp_acc <- function(x) {
#   if (x %in% s700) {
#     samp_acc <- as.character( samp_acc_s1_s700[ samp_acc_s1_s700[,"sample_name"]==x, "accession" ] )
#   } else {
#     samp_acc <- as.character( samp_acc_s701_s1319[ samp_acc_s701_s1319[,"sample_name"]==x, "accession" ] )
#   }
#   return(samp_acc)
# }

# ****** #
# read file containing the fastq filenames for each sample
filenames <- read.delim(sprintf("%s/SRA_submission/filenames.csv", p2_dir), row.names = 1)

# ****** #
for (lab in c("s1-s824","s825-s1648")) {
  if (lab=="s1-s824") {
    samp_vec <- s.1.824
  } else {
    samp_vec <- s.825.1648
  }
  
  meta.table <- matrix(NA, nrow = length(samp_vec), ncol = length(meta_labs))
  rownames(meta.table) <- samp_vec
  colnames(meta.table) <- meta_labs
  
  for (samp in samp_vec) {
    
    meta.table[ samp, "sample_name" ] <- samp
    # meta.table[ samp, "bioproject_accession" ] <- "PRJNA667146"
    # meta.table[ samp, "biosample_accession" ] <- get_samp_acc(samp)
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
  
  write.csv(meta.table, sprintf('%s/SRA_submission/SRA_metadata_%s.csv', p2_dir, lab), row.names = F)
}
# ****** #







# ****************************************************************************************************************** #
###### Table of MALDI-TOF results with anonymyzed sample names matching github table ######
# ****************************************************************************************************************** #

samps.anon <- paste0("S_",1:nrow(SLL2.meta))
names(samps.anon) <- rownames(SLL2.meta)

fung.anon <- fung
colnames(fung.anon) <- c("Center","Plate","Num_Yeast_colonies","Num_Mold_colonies","Bacteria",
                         "Bacteria_detected_by_MALDI","Yeast_species")
# first remove the 4 samples that are not used in the 1648 SLL2 samples
fung.anon <- fung.anon[ rownames(fung.anon) %in% names(samps.anon), ]
rownames(fung.anon) <- samps.anon[ rownames(fung.anon) ]

write.csv(fung.anon, sprintf("%s/SLL2.MALDI_results.csv", p2_dir))

# ****************************************************************************************************************** #
























# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# compare heatmap of correlations with water values to SLL1 ####


gen.cors <- gen.ps <- gen.cor.ps <- matrix(NA, nrow=nrow(gloms_clr$Genus), ncol=length(cont_water_data))
rownames(gen.cors) <- rownames(gen.ps) <- rownames(gen.cor.ps) <- rownames(gloms_clr$Genus)
colnames(gen.cors) <- colnames(gen.ps) <- colnames(gen.cor.ps) <- cont_water_data

for (gen in rownames(gloms_clr$Genus)) {
  for (wat in cont_water_data) {
    if (gen %in% rownames(Anova.pvals$clr$Genus[[ wat ]]) ) {
      cor.ps <- cor.test(gloms_clr$Genus[ gen, rownames(SLL2.meta) ], SLL2.meta[ , wat], method = "pearson")
      
      gen.cors[ gen, wat ] <- cor.ps$estimate
      gen.ps[ gen, wat ] <- Anova.pvals$clr$Genus[[ wat ]][ gen, wat ]
      gen.cor.ps[ gen, wat ] <- cor.ps$p.value
      
    } else {
      
      gen.cors[ gen, wat ] <- 0
      gen.ps[ gen, wat ] <- 1
      gen.cor.ps[ gen, wat ] <- 1
      
    }
  }
}









# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# test MDiNE package for comparing co-occurence networks between groups ####

library(mdine)


# they use table of counts, with a column at the end for the sum of counts of all OTUs not included
#  so here will use 20 most common genera as a start, so 21st column has to be sum of genera not in top20
top20 <- names(sort(rowMeans(gloms_rel$Genus), decreasing = T)[1:20])
Y <- as.data.frame( t(gloms$Genus[top20,]) )
Y$ref <- colSums(gloms$Genus[ rownames(gloms$Genus)[ ! rownames(gloms$Genus) %in% top20 ], ])

# then make covariate table with the variables of choice
network.covs <- c("Chronic_disorder","Gender","Age","Population","seqGroup")
network.covs <- c("Stomatotype_Aitchison","Gender","Age","Population","seqGroup")
network.covs <- c("Smoker","Gender","Age","Population","seqGroup")
netCovTab <- SLL2.meta[ , network.covs ]

# get a subsample to test first
netSamps <- sort( sample(rownames(SLL2.meta), 50) )

Y <- Y[ netSamps, ]
netCovTab <- netCovTab[ netSamps, ]


# make model matrix
X <- model.matrix(~Chronic_disorder, data = netCovTab)
X <- model.matrix(~Stomatotype_Aitchison, data = netCovTab)
X.seqGroup <- model.matrix(~seqGroup, data = netCovTab)
X.Smoker <- model.matrix(~Smoker, data = netCovTab)


# run mdine
md.fit <- mdine(Y=as.matrix(Y), X=X, Z=X[,2], mc.cores = 4, iter = 500)
md.fit.seqGroup <- mdine(Y=as.matrix(Y), X=X.seqGroup, Z=X.seqGroup[,2], mc.cores = 4, iter = 500)
md.fit.Smoker <- mdine(Y=as.matrix(Y), X=X.Smoker, Z=X.Smoker[,2], mc.cores = 4, iter = 500)

mdf0 <- md.fit$post_mean$invsigma0
mdf1 <- md.fit$post_mean$invsigma1
colnames(mdf0) <- colnames(mdf1) <- rownames(mdf0) <- rownames(mdf1) <- top20

plot_networks(md.fit)
plot_networks(md.fit.seqGroup)
plot_networks(md.fit.Smoker)



# Weighted adjacency matrices based on each precision matrix
adj <- ci2adj(md.fit, weighted = T)
adj


# Weighted adjacency matrices based on each precision matrix
ig0 <- adj2ig(adj$adj0)
igraph::plot.igraph(ig0)


ig1 <- adj2ig(adj$adj1)
igraph::plot.igraph(ig1)

# ****************************************************************************************************************** #


bin_vars <- c("Gender","Smoker","Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement","Mouth_wounds",
              "Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome","Gingivitis_periodontitis",
              "Eating_disorder","Medications","Antibiotics","Analgesics","Vitamin_supplements",
              "Asthma","Wheezing",
              "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen","Allergy.Animals",
              "Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
              "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary",
              "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
              "Kissing_partner",
              "MALDI.Yeast_detected","MALDI.Mold_detected",
              additional_diseases,"seqGroup")

bin_vars.short <- c("Gender","Smoker",#"Braces",
                    # "Mouth_piercing",#"Fluoride_toothpaste","Fluoride_supplement",#"Mouth_wounds",
                    "Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome",#"Gingivitis_periodontitis",
                    #"Eating_disorder","Medications",
                    "Antibiotics","Analgesics","Vitamin_supplements",
                    # "Asthma","Wheezing",
                    "Allergy",#"Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen","Allergy.Animals",
                    # "Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
                    # "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary",
                    # "Bite_nails","Hair_in_mouth","Chew_pens",
                    "Wash_hands_before_eat","Wash_hands_after_bathroom",
                    "Kissing_partner",
                    "MALDI.Yeast_detected","MALDI.Mold_detected",
                    # additional_diseases,
                    "seqGroup")

# they use table of counts, with a column at the end for the sum of counts of all OTUs not included
#  so here will use 20 most common genera as a start, so 21st column has to be sum of genera not in top20
top20 <- names(sort(rowMeans(gloms_rel$Genus), decreasing = T)[1:20])
Y <- as.data.frame( t(gloms$Genus[top20,]) )
Y$ref <- colSums(gloms$Genus[ rownames(gloms$Genus)[ ! rownames(gloms$Genus) %in% top20 ], ])

# then make covariate table with the variables of choice
network.covs <- c("Chronic_disorder","Gender","Age","Population","seqGroup")

# get a subsample to test first
netSamps <- sort( sample(rownames(SLL2.meta), 100) )

# ********************************************************************** #
# md.fits <- list()
# md.fits[[ "netSamps" ]] <- netSamps
# 
# for (binTrait in bin_vars.short) {
#   
#   print(binTrait)
#   
#   if (binTrait %in% c("Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome","Gingivitis_periodontitis",
#                              "Eating_disorder",additional_diseases)) {
#     nCovs <- network.covs[ network.covs != "Chronic_disorder" ]
#   } else {
#     nCovs <- network.covs[ network.covs != binTrait ]
#   }
#   
#   # make model matrix
#   X <- model.matrix(as.formula( sprintf("~%s", binTrait) ), data = SLL2.meta[ netSamps, c(binTrait, nCovs) ] )
#   
#   # run mdine
#   md.fits[[ binTrait ]] <- mdine(Y=as.matrix(Y[ netSamps, ]), X=X, Z=X[,2], mc.cores = 4, iter = 500)
#   
# }
# saveRDS(md.fits, sprintf("%s/R_objects/md.fits.rds", p2_dir))
# ********************************************************************** #
md.fits <- readRDS(sprintf("%s/R_objects/md.fits.rds", p2_dir))


binTrait <- "Gender"

plot_networks(md.fits[[ binTrait ]])

# Weighted adjacency matrices based on each precision matrix
adj <- ci2adj(md.fits[[ binTrait ]], weighted = T)
adj


# Weighted adjacency matrices based on each precision matrix
ig0 <- adj2ig(adj$adj0)
igraph::plot.igraph(ig0)


ig1 <- adj2ig(adj$adj1)
igraph::plot.igraph(ig1)




taxTables.both$Genus[top20, "Phylum"]

# ****************************************************************************************************************** #




# ****************************************************************************************************************** #
# instead of 20 most common, try top20 drivers of Aitchison clusters

aitch_drivers <- get_drivers("Genus", "Aitchison", gloms_clr, SLL2.meta, nTop = 10)
Y <- as.data.frame( t(gloms$Genus[ names(unlist(unname(aitch_drivers))), ]) )
Y$ref <- colSums(gloms$Genus[ rownames(gloms$Genus)[ ! rownames(gloms$Genus) %in% names(unlist(unname(aitch_drivers))) ], ])

# then make covariate table with the variables of choice
network.covs <- c("Chronic_disorder","Gender","Age","Population","seqGroup")
network.covs <- c("Stomatotype_Aitchison","Gender","Age","Population","seqGroup")
netCovTab <- SLL2.meta[ , network.covs ]

# get a subsample to test first
netSamps <- sort( sample(rownames(SLL2.meta), 50) )

Y <- Y[ netSamps, ]
netCovTab <- netCovTab[ netSamps, ]


# make model matrix
X <- model.matrix(~Chronic_disorder, data = netCovTab)
X <- model.matrix(~Stomatotype_Aitchison, data = netCovTab)


# run mdine
md.fit.aitch_drivers <- mdine(Y=as.matrix(Y), X=X, Z=X[,2], mc.cores = 4, iter = 500)

mdf0 <- md.fit.aitch_drivers$post_mean$invsigma0
mdf1 <- md.fit.aitch_drivers$post_mean$invsigma1
colnames(mdf0) <- colnames(mdf1) <- rownames(mdf0) <- rownames(mdf1) <- names(unlist(unname(aitch_drivers)))

plot_networks(md.fit.aitch_drivers)



# Weighted adjacency matrices based on each precision matrix
adj <- ci2adj(md.fit.aitch_drivers, weighted = T)
adj


# Weighted adjacency matrices based on each precision matrix
ig0 <- adj2ig(adj$adj0)
igraph::plot.igraph(ig0)

ig1 <- adj2ig(adj$adj1)
igraph::plot.igraph(ig1)


# ****************************************************************************************************************** #
# ****************************************************************************************************************** #



# Run MDiNE with disorder subsamplings ####

disorder <- trait <- "Cystic_fibrosis"
i <- 1

# Must get disorder samples and matched controls from the 'i'th subsampling that was already performed
disSamps <- disSubsTests$default[[ disorder ]][[ i ]]$samples

# ******************************************************* #
# they use table of counts, with a column at the end for the sum of counts of all OTUs not included
#  so here will use 50 most common genera as a start, so 51st column has to be sum of genera not in top50
top50 <- names(sort(rowMeans(gloms_rel$Genus), decreasing = T)[1:40])
# this list should be the same for all subsamplings to maintain consistency

# but the table of counts will of course only include the matched controls with the disorder samples
# disorder.samps <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
# controls <- get_samps_for_test("disorder", SLL2.meta, i, phy=SLL2, defaultDisorder = T)

# get tables with only those samples of interest
glomTab <- gloms$Genus[ , disSamps ] #[ , c(disorder.samps, controls) ]
meta.disorder <- SLL2.meta[ disSamps, ] #[ c(disorder.samps, controls), ]

Y <- as.data.frame( t(glomTab[top50,]) )
Y$ref <- colSums(glomTab[ rownames(glomTab)[ ! rownames(glomTab) %in% top50 ], ])
# ******************************************************* #




# ******************************************************* #
# then make covariate table with the variables of choice
# the first variable will be the trait of interest, and will be passed from the command line
# trait <- argsFromCommandLine[1]
# print(trait)

# will always include these 5
network.covs <- c("Chronic_disorder","Gender","Age","Population")
# remove correlated covariate when applicable
if (trait %in% c("Municipal_zone",cont_water_data, group_water_data, "Community","Province")) {
  network.covs <- c("Chronic_disorder","Gender","Age")
  
} else if (trait %in% c("Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome",
                        "Gingivitis_periodontitis","Eating_disorder",
                        additional_diseases[ ! additional_diseases %in% c("Samp_collect_issues","Birth_control",
                                                                          "Antihistamines")])) {
  network.covs <- c("Gender","Age","Population")
}
# then make vector with trait and covariates, removing duplicates if trait of interest is one of the standard covariates
network.covs <- c(trait, network.covs[ network.covs != trait ])
# ******************************************************* #

# ******************************************************* #
# get data.frame of just those variables of interest
netCovTab <- meta.disorder[ , network.covs ]

# cant use samples that were NA for given trait, must remove from both count and var tables
non.na.samps <- rownames(meta.disorder)[ ! is.na(meta.disorder[ , trait] )]
Y <- Y[ non.na.samps, ]
netCovTab <- netCovTab[ non.na.samps, ]
# ******************************************************* #

# ******************************************************* #
# make model matrix
X <- model.matrix(as.formula(sprintf("~%s", trait)), data = netCovTab)

print(dim(Y))
print(dim(X))
print(dim(netCovTab))
# ******************************************************* #


# ******************************************************* #
# run mdine
md.fit <- mdine(Y=as.matrix(Y), X=X, Z=X[,2], mc.cores = 4, chains = 4, iter = 500)
# ******************************************************* #
# save object containing result, as well as the samples that were selected for the test
fit_and_samps <- list("md.fit"=md.fit, "samples"=disSamps) #c(disorder.samps, controls))

# make directory if necessary
dir.create(sprintf("%s/R_objects/MDiNE/%s/%s", p2_dir, disorder, trait), showWarnings = F)
saveRDS(fit_and_samps, sprintf("%s/R_objects/MDiNE/%s/%s/md.fit.%s.%s.rds", p2_dir, disorder, trait, trait, i))
# ******************************************************* #

md.fit <- readRDS(sprintf("%s/R_objects/MDiNE/%s/%s/md.fit.%s.%s.rds", p2_dir, disorder, trait, trait, i))






# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####



















# ****************************************************************************************************************** #

# test SpiecEasi package for comparing co-occurence networks between groups ####

library(SpiecEasi)

# can run with a phyloseq object.
# will do so with one at the genus level
tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

se.mb.sll2 <- spiec.easi(tagl, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))

ig2.mb <- adj2igraph(getRefit(se.mb.sll2),  vertex.attr=list(name=taxa_names(tagl)))
plot_network(ig2.mb, tagl, type='taxa', color="Phylum")


# ************************************************** #
# can run with a phyloseq object.
# will do so with one at the genus level
tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

tagl.Yes <- subset_samples(tagl, Chronic_disorder=="Yes")
tagl.No  <- subset_samples(tagl, Chronic_disorder=="No")

se.mb.sll2.Yes <- spiec.easi(tagl.Yes, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
se.mb.sll2.No  <- spiec.easi(tagl.No, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))

ig2.mb.Yes <- adj2igraph(getRefit(se.mb.sll2.Yes),  vertex.attr=list(name=taxa_names(tagl.Yes)))
plot_network(ig2.mb.Yes, tagl.Yes, type='taxa', color="Phylum")

ig2.mb.No <- adj2igraph(getRefit(se.mb.sll2.No),  vertex.attr=list(name=taxa_names(tagl.No)))
plot_network(ig2.mb.No, tagl.No, type='taxa', color="Phylum")

# ************************************************** #

# can run with a phyloseq object.
# will do so with one at the genus level
tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

tagl.1 <- subset_samples(tagl, Stomatotype_Aitchison=="1")
tagl.2  <- subset_samples(tagl, Stomatotype_Aitchison=="2")

se.mb.sll2.1 <- spiec.easi(tagl.1, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
se.mb.sll2.2  <- spiec.easi(tagl.2, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))

ig2.mb.1 <- adj2igraph(getRefit(se.mb.sll2.1),  vertex.attr=list(name=taxa_names(tagl.1)))
plot_network(ig2.mb.1, tagl.1, type='taxa', color="Phylum")

ig2.mb.2 <- adj2igraph(getRefit(se.mb.sll2.2),  vertex.attr=list(name=taxa_names(tagl.2)))
plot_network(ig2.mb.2, tagl.2, type='taxa', color="Phylum")


# ************************************************** #
# SpiecEasi for Cystic Fibrosis ####

# can run with a phyloseq object.
# will do so with one at the genus level
tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

tagl.CF.Yes <- subset_samples(tagl, Cystic_fibrosis=="Yes")
CF.mC <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
# CF.mC <- disSubsTests$default$Cystic_fibrosis$`1`$samples
CF.mC <- CF.mC[ SLL2.meta[CF.mC, "Cystic_fibrosis"]=="No" ]
tagl.CF.all_No  <- subset_samples(tagl, sample_names(tagl) %in% CF.mC)

# # remove those taxa that have 0 counts in either Yes or No
# CF.non0s <- taxa_names(tagl.CF.all_No)[ taxa_sums(tagl.CF.all_No)>0 & taxa_sums(tagl.CF.Yes)>0 ]
# tagl.CF.all_No <- prune_taxa(CF.non0s, tagl.CF.all_No)
# tagl.CF.Yes <- prune_taxa(CF.non0s, tagl.CF.Yes)


se.mb.sll2.CF.Yes <- spiec.easi(tagl.CF.Yes, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
se.mb.sll2.CF.all_No  <- spiec.easi(tagl.CF.all_No, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))

ig2.mb.CF.Yes <- adj2igraph(getRefit(se.mb.sll2.CF.Yes),  vertex.attr=list(name=taxa_names(tagl.CF.Yes)), rmEmptyNodes = T)
plot_network(ig2.mb.CF.Yes, tagl.CF.Yes, type='taxa', color="Phylum", title = "Co-occurrence network: CF samples")

ig2.mb.CF.all_No <- adj2igraph(getRefit(se.mb.sll2.CF.all_No),  vertex.attr=list(name=taxa_names(tagl.CF.all_No)), rmEmptyNodes = T)
plot_network(ig2.mb.CF.all_No, tagl.CF.all_No, type='taxa', color="Phylum", title = "Co-occurrence network: all mC samples")


ig2.mb.CF.Yes.allNodes <- adj2igraph(getRefit(se.mb.sll2.CF.Yes),  vertex.attr=list(name=taxa_names(tagl.CF.Yes)), rmEmptyNodes = F)
ig2.mb.CF.all_No.allNodes <- adj2igraph(getRefit(se.mb.sll2.CF.all_No),  vertex.attr=list(name=taxa_names(tagl.CF.all_No)), rmEmptyNodes = F)


se.slr.sll2.CF.Yes <- spiec.easi(tagl.CF.Yes, method="slr", r=10, lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=30, ncores=2))
se.slr.sll2.CF.all_No  <- spiec.easi(tagl.CF.all_No, method="slr", r=10, lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=30, ncores=2))



# ************************************************** #
huge::huge.roc(se.mb.sll2.Yes$est$path, graph, verbose=FALSE)


# ************************************************** #









# ****************************************************************************************************************** ####

# ************************************************** #

# Calculate Hamming distances between networks ####
# get Hamming distances between CF and all matched control groups, as well as between all matched control groups 
library(nettools)
library(SpiecEasi)

tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

# first get objects for CF samples
samps.CF <- rownames(SLL2.meta[ SLL2.meta[, "Cystic_fibrosis"]=="Yes", ])
tagl.CF.Yes  <- subset_samples(tagl, sample_names(tagl) %in% samps.CF)
# then for all of the matched controls together
samps.mC <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
samps.mC <- samps.mC[ SLL2.meta[samps.mC, "Cystic_fibrosis"]=="No" ]
tagl.CF.all_No  <- subset_samples(tagl, sample_names(tagl) %in% samps.mC)

# # remove those taxa that dont have at least 100 counts across all CF + mC samples,
# #   then will make tagl.CF.Yes and talg.CF.all_No again below
# tagl <- prune_taxa(taxa_sums(tagl.CF.all_No) + taxa_sums(tagl.CF.Yes) > 100, tagl)

# ************************** #
# from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
# filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
#    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
library(seqtime)
minCount <- 15
minOcc   <- 20
g10 <- as.data.frame(gloms$Genus[ , c(samps.CF, samps.mC) ])
g10[ g10 <= minCount ] <- 0
# filterobj <- filterTaxonMatrix(gloms$Genus[ , c(samps.CF, samps.mC) ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
otus.f <- filterobj$mat

# replace the values that were made to be 0s for the filtering
anti.g10 <- as.data.frame(gloms$Genus[ , c(samps.CF, samps.mC) ])
anti.g10[ anti.g10 > minCount ] <- 0
# replace them into g10
g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
# first update the values that were made 0
upd.otus_f <- otus.f
upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(gloms$Genus[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], c(samps.CF, samps.mC) ])
# then update the summed row from the removed values
upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))

# finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
taxa.f <- tax_table(tagl)[setdiff(1:nrow(tax_table(tagl)),filterobj$filtered.indices),]
dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
taxa.f <- rbind(taxa.f, dummyTaxonomy)
rownames(taxa.f)[nrow(taxa.f)] <- "0"
rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"

tagl <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                 sample_data(SLL2)[ c(samps.CF, samps.mC), ],
                 tax_table(taxa.f))
# ************************** #


# first get objects for CF samples
samps.CF <- rownames(SLL2.meta[ SLL2.meta[, "Cystic_fibrosis"]=="Yes", ])
tagl.CF.Yes  <- subset_samples(tagl, sample_names(tagl) %in% samps.CF)
# se.mb.sll2.CF.Yes <- spiec.easi(tagl.CF.Yes, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.CF.Yes, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.CF.Yes.rds", p2_dir))
se.mb.sll2.CF.Yes <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.CF.Yes.rds", p2_dir))
igraphs.CF.Yes.allNodes <- adj2igraph(getRefit(se.mb.sll2.CF.Yes),  vertex.attr=list(name=taxa_names(tagl.CF.Yes)), rmEmptyNodes = F)
igraphs.CF.Yes.noEmpty <- adj2igraph(getRefit(se.mb.sll2.CF.Yes),  vertex.attr=list(name=taxa_names(tagl.CF.Yes)), rmEmptyNodes = T)
# plot_network(igraphs.CF.Yes.noEmpty, tagl.CF.Yes, type='taxa', color="Phylum", title = "Co-occurrence network: CF samples")



# then for all of the matched controls together
samps.mC <- unique(unlist(lapply(disSubsTests$default$Cystic_fibrosis, function(x) x$samples)))
samps.mC <- samps.mC[ SLL2.meta[samps.mC, "Cystic_fibrosis"]=="No" ]
tagl.CF.all_No  <- subset_samples(tagl, sample_names(tagl) %in% samps.mC)
# se.mb.sll2.CF.all_No <- spiec.easi(tagl.CF.all_No, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.CF.all_No, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.CF.all_No.rds", p2_dir))
se.mb.sll2.CF.all_No <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.CF.all_No.rds", p2_dir))
igraphs.CF.all_No.allNodes <- adj2igraph(getRefit(se.mb.sll2.CF.all_No),vertex.attr=list(name=taxa_names(tagl.CF.all_No)), rmEmptyNodes = F)
igraphs.CF.all_No.noEmpty <- adj2igraph(getRefit(se.mb.sll2.CF.all_No),vertex.attr=list(name=taxa_names(tagl.CF.all_No)), rmEmptyNodes = T)



# CF_mC.hammings <- list()
# SE_objects <- list()
# igraphs.allNodes <- list()
# igraphs.noEmpty <- list()
# 
# for (i in names(disSubsTests$default$Cystic_fibrosis)) {
#   
#   print(i)
#   
#   # get samples
#   samps <- disSubsTests$default$Cystic_fibrosis[[ i ]]$samples
#   samps.mC <- samps[ SLL2.meta[samps, "Cystic_fibrosis"]=="No" ]
#   
#   # make phyloseq objects
#   tagl.CF.No   <- subset_samples(tagl, sample_names(tagl) %in% samps.mC)
#   
#   # calculate spiec.easi object
#   se.mb.sll2.CF.No  <- spiec.easi(tagl.CF.No, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
#   SE_objects[[ i ]] <- se.mb.sll2.CF.No
#   
#   # get igraphs with all nodes for hamming calculations
#   ig2.mb.CF.No.allNodes <- adj2igraph(getRefit(se.mb.sll2.CF.No),  vertex.attr=list(name=taxa_names(tagl.CF.No)), rmEmptyNodes = F)
#   igraphs.allNodes[[ i ]] <- ig2.mb.CF.No.allNodes
#   
#   # get igraphs without empty nodes for plotting
#   ig2.mb.CF.No <- adj2igraph(getRefit(se.mb.sll2.CF.No),  vertex.attr=list(name=taxa_names(tagl.CF.No)), rmEmptyNodes = T)
#   igraphs.noEmpty[[ i ]] <- ig2.mb.CF.No
#   
#   # calculating Hamming distances between CF and mC networks
#   CF_mC.hammings[[ i ]] <- netdist(ig2.mb.CF.No.allNodes, igraphs.CF.Yes.allNodes, d = "Hamming")
#   
# }
# 
# saveRDS(CF_mC.hammings, file = sprintf("%s/R_objects/SpiecEasi/CF_mC.hammings.rds", p2_dir))
# saveRDS(SE_objects, file = sprintf("%s/R_objects/SpiecEasi/SE_objects.rds", p2_dir))
# saveRDS(igraphs.allNodes, file = sprintf("%s/R_objects/SpiecEasi/igraphs.allNodes.rds", p2_dir))
# saveRDS(igraphs.noEmpty, file = sprintf("%s/R_objects/SpiecEasi/igraphs.noEmpty.rds", p2_dir))

# ************************************************** #

CF_mC.hammings   <- readRDS(sprintf("%s/R_objects/SpiecEasi/CF_mC.hammings.rds", p2_dir))
SE_objects       <- readRDS(sprintf("%s/R_objects/SpiecEasi/SE_objects.rds", p2_dir))
igraphs.allNodes <- readRDS(sprintf("%s/R_objects/SpiecEasi/igraphs.allNodes.rds", p2_dir))
igraphs.noEmpty  <- readRDS(sprintf("%s/R_objects/SpiecEasi/igraphs.noEmpty.rds", p2_dir))

# ************************************************** #


# # get pairwise distances between networks of all the matched control subsets
# mC_mC.hammings <- list()
# 
# for (i in 1:99) {
#   for (j in 2:100) {
#     mC_mC.hammings[[ sprintf("%s_%s", i, j) ]] <- netdist(igraphs.allNodes[[ as.character(i) ]], 
#                                                           igraphs.allNodes[[ as.character(j) ]], 
#                                                           d = "Hamming")
#   }
# }
# 
# saveRDS(mC_mC.hammings, file = sprintf("%s/R_objects/SpiecEasi/mC_mC.hammings.rds", p2_dir))

mC_mC.hammings   <- readRDS(sprintf("%s/R_objects/SpiecEasi/mC_mC.hammings.rds", p2_dir))
# ************************************************** #

mean(unlist(CF_mC.hammings))
mean(unlist(mC_mC.hammings))

sd(unlist(CF_mC.hammings))
sd(unlist(mC_mC.hammings))

# dist between CF and all mC
netdist(ig2.mb.CF.all_No.allNodes, igraphs.CF.Yes.allNodes, d = "Hamming")


# plot a single mC network
i <- "1"

samps <- disSubsTests$default$Cystic_fibrosis[[ i ]]$samples
tagl.CF.No   <- subset_samples(tagl, sample_names(tagl) %in% samps[ SLL2.meta[samps, "Cystic_fibrosis"]=="No" ])

plot_network(igraphs.noEmpty[[ i ]], 
             tagl.CF.No, 
             type='taxa', color="Phylum", title = sprintf("Co-occurrence network: mC[ %s ] samples", i))


# boxplot of distances
hammings <- list("CF vs mC" = unlist(CF_mC.hammings), "mC vs mC" = unlist(mC_mC.hammings))
boxplot(hammings)
kruskal.test(hammings)

h.df <- data.frame(unlist(hammings), c(rep("CF vs mC", length(hammings$`CF vs mC`)), rep("mC vs mC", length(hammings$`mC vs mC`))))
colnames(h.df) <- c("Hamming", "comp")

library(ggpubr)
ggplot(h.df, aes(x=comp, y=Hamming, fill=comp)) +
  geom_boxplot() +
  ggtitle("(c) Hamming distances between networks") +
  scale_fill_manual(values=c("#E4BA4E","#73BDD3")) +
  theme_minimal() +
  theme(axis.title = element_text(size=30), axis.text = element_text(size=27, face = "bold"),
        plot.title = element_text(size=25, face = "bold")) +
  guides(fill=F) + xlab(NULL) +
  stat_compare_means(label.x = 1.2, cex=7.75)


# ************************************************** #
# ************************************************** #

# use package intergraph to convert igraphs to networks, from which I can pull out all the edge connections ####
library(intergraph)

edge.mat.CF <- as.matrix(asNetwork(igraphs.CF.Yes.allNodes))

# for each genus, get vector of those other genera with which it has an edge
gen_edges.CF <- list()
for (gen in rownames(edge.mat.CF)) {
  
  gen_edges.CF[[ gen ]] <- colnames(edge.mat.CF)[ edge.mat.CF[ gen, ] != 0 ]
  
}



# ******************** #

edge.mats.CF_no <- list()
gen_edges.CF_no <- list()
for (i in names(igraphs.allNodes)) {
  # get matrix of node connections
  edge.mats.CF_no[[ i ]] <- as.matrix(asNetwork(igraphs.allNodes[[ i ]]))
  
  # get genus connection vectors
  gen_edges.CF_no[[ i ]] <- list()
  for (gen in rownames(edge.mats.CF_no[[ i ]])) {
    
    gen_edges.CF_no[[ i ]][[ gen ]] <- colnames(edge.mats.CF_no[[ i ]])[ edge.mats.CF_no[[ i ]][ gen, ] != 0 ]
    
  }
}




# ******************** #
# first test with 1st subsamp 
#   - for each genus, get vector of connected genera that are same in both CF and mC
#   - vector of genera only connected in CF
#   - vector of genera only connected in mC

gen_connect <- list("both"=list(), "CF"=list(), "mC"=list())

for (i in names(gen_edges.CF_no)) {
  gen_connect$both[[ i ]] <- sapply(names(gen_edges.CF), function(gen) 
    unique(c(gen_edges.CF[[ gen ]][ gen_edges.CF[[ gen ]] %in% gen_edges.CF_no[[ i ]][[ gen ]] ],
      gen_edges.CF_no[[ i ]][[ gen ]][ gen_edges.CF_no[[ i ]][[ gen ]] %in% gen_edges.CF[[ gen ]] ])) )
  
  gen_connect$CF[[ i ]] <- sapply(names(gen_edges.CF), function(gen) 
    gen_edges.CF[[ gen ]][ ! gen_edges.CF[[ gen ]] %in% gen_edges.CF_no[[ i ]][[ gen ]] ])
  
  gen_connect$mC[[ i ]] <- sapply(names(gen_edges.CF), function(gen) 
    gen_edges.CF_no[[ i ]][[ gen ]][ ! gen_edges.CF_no[[ i ]][[ gen ]] %in% gen_edges.CF[[ gen ]] ])
}

# *************************************** #
gen_con_freqs <- function(gc_list, gen) {
  
  for (con in rev(names(gc_list))) {
    cat(sprintf("**** %s connections in %s %s ****\n", gen, con, ifelse(con=="both","","only")))
    print(sort(table(unlist(lapply(gc_list[[ con ]], function(x) x[[ gen ]])))))
  }
  
  
  # cat(sprintf("**** %s connections in mC only ****", gen))
  # print(sort(table(unlist(lapply(gc_list$mC, function(x) x[[ gen ]])))))
  # 
  # cat(sprintf("\n**** %s connections in CF only ****", gen))
  # print(sort(table(unlist(lapply(gc_list$CF, function(x) x[[ gen ]])))))
  # 
  # cat(sprintf("\n**** %s connections in both ****", gen))
  # print(sort(table(unlist(lapply(gc_list$both, function(x) x[[ gen ]])))))
  
}
# *************************************** #

gen_con_freqs(gen_connect, "Porphyromonas")

gen_con_freqs(gen_connect, "Pseudomonas")


CF-related: c("Pseudomonas","Rothia","Staphylococcus","Chryseobacterium","Microbacterium","Brevundimonas",
              "Stenotrophomonas","Streptococcus","Delftia","Comamonas","Scardovia","Mobiluncus","Sphingobacterium",
              "Mogibacterium", "Desulfobulbus")

- expect Desulfobulbus-Pseudomonas connection in CF?
- also Rothia-Pseudomonas
- maybe anti-correlation between Pseudomonas-Streptococcus if S. salivarius inhibiting P. aeruginosa
   - or positive correlation between Pseudomonas-Streptococcus if other Streptococcus producing 2,3-butanediol
   - in that case also Rothia-Streptococcus, since Rothia also uses it.
- Also Rothia-Stenotrophomonas from reduction of Fe3+ by S. maltophilia so R. mucilaginosa can use Fe2+



Periodontitis-related: c("Treponema","Tannerella","Porphyromonas","Fusobacterium","Parvimonas","Aggregatibacter",
                         "Peptostreptococcus","Clostrdiiales Family_XIII","Clostridiales_vadinBB60_group",
                         "Mogibacterium","Desulfobulbus","Stenotrophomonas")
- Staphylococcus and Rothia perhaps involved in early development of periodontitis, check for connection with all these



caries-related: c("Chryseobacterium","Comamonas","Scardovia","Mogibacterium","Aggregatibacter","Alloprevotella")

caries-free: c("Treponema","Fusobacterium","Bergeyella","Patescibacteria","Rothia","Stenotrophomonas","Delftia")




# ************************************************** #
# ************************************************** #

# use package intergraph to convert igraphs to networks, from which I can pull out all the edge connections
library(intergraph)

edge.mat.CF.all_No <- as.matrix(asNetwork(igraphs.CF.all_No.allNodes))

# for each genus, get vector of those other genera with which it has an edge
gen_edges.CF.all_No <- list()
for (gen in rownames(edge.mat.CF.all_No)) {
  
  gen_edges.CF.all_No[[ gen ]] <- colnames(edge.mat.CF.all_No)[ edge.mat.CF.all_No[ gen, ] != 0 ]
  
}



# ******************** #

gen_connect.to_all_No <- list()

gen_connect.to_all_No[[ "both" ]] <- sapply(names(gen_edges.CF), function(gen) 
    gen_edges.CF[[ gen ]][ gen_edges.CF[[ gen ]] %in% gen_edges.CF.all_No[[ gen ]] ])
  
gen_connect.to_all_No[[ "CF" ]] <- sapply(names(gen_edges.CF), function(gen) 
    gen_edges.CF[[ gen ]][ ! gen_edges.CF[[ gen ]] %in% gen_edges.CF.all_No[[ gen ]] ])
  
gen_connect.to_all_No[[ "mC" ]] <- sapply(names(gen_edges.CF), function(gen) 
    gen_edges.CF.all_No[[ gen ]][ ! gen_edges.CF.all_No[[ gen ]] %in% gen_edges.CF[[ gen ]] ])


# *************************************** #
gen_con_freqs.all_No <- function(gc_list, gen) {
  
  for (con in rev(names(gc_list))) {
    cat(sprintf("**** %s connections in %s %s ****\n", gen, con, ifelse(con=="both","","only")))
    print(sort(gc_list[[ con ]][[ gen ]]))
  }
  
  # cat(sprintf("**** %s connections in mC only ****\n", gen))
  # print(sort(gc_list$mC[[ gen ]]))
  # 
  # cat(sprintf("\n**** %s connections in CF only ****\n", gen))
  # print(sort(gc_list$CF[[ gen ]]))
  # 
  # cat(sprintf("\n**** %s connections in both ****\n", gen))
  # print(sort(gc_list$both[[ gen ]]))
  
}
# *************************************** #

gen_con_freqs.all_No(gen_connect.to_all_No, "Porphyromonas")







# ****************************************************************************************************************** #
# Make network plots ####
# Following this tutorial here: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html
# also check this one: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php


# **************************************************************************** #


CF.net   <- get_network_objects("CF", se.mb.sll2.CF.Yes, tagl.CF.Yes, taxTables.both$Genus)
mC.net   <- get_network_objects("mC", se.mb.sll2.CF.all_No, tagl.CF.all_No, taxTables.both$Genus)

i <- "1"
samps.mC_i <- disSubsTests$default$Cystic_fibrosis[[ i ]]$samples
samps.mC_i <- samps.mC_i[ SLL2.meta[samps.mC_i, "Cystic_fibrosis"]=="No" ]
mC_i.net <- get_network_objects("mC_i", SE_objects, subset_samples(tagl, sample_names(tagl) %in% samps.mC_i), 
                                taxTables.both$Genus, i=i)

# **************************************************************************** #


CF.net <- get_network_objects("CF", se.mb.sll2.CF.Yes, tagl.CF.Yes, taxTables.both$Genus, 
                              addLabels = T, Vprop = T, Eprop = T)



CF.net <- get_network_objects("CF", se.mb.sll2.CF.Yes, tagl.CF.Yes, taxTables.both$Genus, 
                              addLabels = T, Vprop = T, Eprop = T,
                              taxa_to_label = c("Alloprevotella","Brevundimonas","unclassified.G20","Solobacterium",
                                                "Aggregatibacter","Mesorhizobium","Actinobacillus",
                                                "unclassified.G33","Fretibacterium","Dialister","Treponema",
                                                "Gemella","Lactobacillus","Fusobacterium","Campylobacter","unclassified.G31",
                                                "Leptotrichia","Pseudomonas"))

c("Chryseobacterium","Microbacterium","Brevundimonas","Peptostreptococccus","Stenotrophomonas","unclassified.G20",
  "Alloprevotella","Streptococcus","Rothia","Staphylococcus","Delftia","Comamonas","unclassified.G43","Scardovia",
  "Treponema","Aggregatibacter","Desulfobulbus","Parvimonas","Bergeyella","unclassified.G19","Mobiluncus","Sphingobacterium",
  "unclassified.G33","Mogibacterium","Fusobacterium","Haemophilus")



# Final paper version *****
CF.net <- get_network_objects("CF", se.mb.sll2.CF.Yes, tagl.CF.Yes, taxTables.both$Genus, 
                              addLabels = T, Vprop = T, Eprop = T,
                              addLegend = F,
                              taxa_to_label = c("Alloprevotella","Brevundimonas", # 100 - mention periodontitis
                                                "Fusobacterium","Lactobacillus", # 100 - negative - lacto inhibits, is assoc with caries, can thrive in CF
                                                "Gemella","Lactobacillus", # 100 - negative gemella is inital cariogen, may pave way for lacto, then inhibited later as gem can also be a periopathogen
                                                "Aggregatibacter","Mesorhizobium","Actinobacillus", # ??? 100 - negative, and 99 - look into mesorhiz more
                                                "Campylobacter","Leptotrichia", # 100 - use of lactate by campylobacter in caries
                                                "Megasphaera","Kingella", # 100 - caries
                                                "Corynebacterium","Hyphomicrobium","Selenomonas", # 100, 100 - biofilms, caries, some CF
                                                "Veillonella","Granulicatella", # 97 - caries, gran is acidogenic, veil consumes lactate
                                                "Atopobium","Stomatobaculum", # 86 - caries, some CF
                                                "Streptococcus","Rothia", # 74 caries and particularly lactate
                                                "Chryseobacterium","Brevundimonas", # 70 CF + caries
                                                "Variovorax","Acinetobacter","Pseudomonas" # 79, 53 - connections to CF, virulence factors
                              ))





# **************************************** #



mC.net <- get_network_objects("mC", se.mb.sll2.CF.all_No, tagl.CF.all_No, taxTables.both$Genus,
                              addLabels = T, Vprop = T, Eprop = T)


mC.net <- get_network_objects("mC", se.mb.sll2.CF.all_No, tagl.CF.all_No, taxTables.both$Genus,
                              addLabels = T, Vprop = T, Eprop = T,
                              taxa_to_label = c("Chryseobacterium","Bradyrhizobium",
                                                "Alloprevotella","Peptostreptococcus","Peptococcus","unclassified.G20","Solobacterium",
                                                "Parvimonas","Catonella","Variovorax","Curvibacter","Porphyromonas","Prevotella","Ralstonia",
                                                "Streptococcus","Haemophilus","Granulicatella","Veillonella","Actinomyces",
                                                "Treponema","Filifactor","Dialister","F0058","unclassified.G10",
                                                "Bergeyella","Neisseria","Acinetobacter",
                                                "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Capnocytophaga",
                                                "Selenomonas","Bifidobacterium","Hyphomicrobium","Johnsonella","Lachnoanaerobaculum",
                                                "Ruminococcaceae_UCG-014","unclassified.G32","Atopobium",
                                                "Mesorhizobium",
                                                "Gemella","Lactobacillus","Fusobacterium","Campylobacter",
                                                "Leptotrichia","Pseudomonas"))

mC.net <- get_network_objects("mC", se.mb.sll2.CF.all_No, tagl.CF.all_No, taxTables.both$Genus,
                              addLabels = T, Vprop = T, Eprop = T,
                              taxa_to_label = c(
                                #"Chryseobacterium","Bradyrhizobium", # caries - against narrative
                                "unclassified.G20","Parvimonas", # 50 - periodontitis
                                "Streptococcus","Gemella", # 92 - caries - against narrative, periodontitis - with narrative
                                "Rothia", "Actinomyces","Granulicatella", # 58, 55 - caries - against narrative
                                "Gemella","Granulicatella", # 72 - caries - against narrative
                                "Treponema","Filifactor","Dialister", # 47, 33 - periodontitis
                                "Prevotella","Veillonella","Solobacterium","Atopobium", # 92, 61, 50 - periodontitis - look into more for the other 3
                                # "Bergeyella","Neisseria", # 39 - low Berg in CF, and its assoc with low caries, Neisseria helps cariogens, so their assoc in mC makes sense
                                "Actinomyces","Stomatobaculum","Rothia","Atopobium", # 62, 58, 46 - caries
                                "Ralstonia","Variovorax", # 96 - AHL-acylase, periodontitis
                                "Bradyrhizobium","Hyphomicrobium" # 86 - periodontitis?
                              ))


# **************************************************************************** #

# For figure 2 of the CF paper


par(mfrow = c(1,2))

CF.net <- get_network_objects("CF", se.mb.sll2.CF.Yes, tagl.CF.Yes, taxTables.both$Genus, 
                              addLabels = T, Vprop = T, Eprop = T,
                              addLegend = F,
                              taxa_to_label = c("Alloprevotella","Brevundimonas", # 100 - mention periodontitis
                                                "Fusobacterium","Lactobacillus", # 100 - negative - lacto inhibits, is assoc with caries, can thrive in CF
                                                "Gemella","Lactobacillus", # 100 - negative gemella is inital cariogen, may pave way for lacto, then inhibited later as gem can also be a periopathogen
                                                "Aggregatibacter","Mesorhizobium",#"Actinobacillus", # ??? 100 - negative, and 99 - look into mesorhiz more
                                                "Campylobacter","Leptotrichia", # 100 - use of lactate by campylobacter in caries
                                                "Megasphaera","Kingella", # 100 - caries
                                                "Corynebacterium","Selenomonas"#,"Hyphomicrobium", # 100, 100 - biofilms, caries, some CF
                                                # "Veillonella","Granulicatella", # 97 - caries, gran is acidogenic, veil consumes lactate
                                                # "Atopobium","Stomatobaculum", # 86 - caries, some CF
                                                # "Streptococcus","Rothia", # 74 caries and particularly lactate
                                                # "Chryseobacterium","Brevundimonas", # 70 CF + caries
                                                # "Variovorax","Acinetobacter","Pseudomonas" # 79, 53 - connections to CF, virulence factors
                              ))

mC.net <- get_network_objects("mC", se.mb.sll2.CF.all_No, tagl.CF.all_No, taxTables.both$Genus,
                              addLabels = T, Vprop = T, Eprop = T,
                              addLegend = F,
                              taxa_to_label = c(
                                #"Chryseobacterium","Bradyrhizobium", # caries - against narrative
                                "unclassified.G20","Parvimonas", # 50 - periodontitis
                                "Streptococcus","Gemella", # 92 - caries - against narrative, periodontitis - with narrative
                                "Rothia", "Actinomyces","Granulicatella", # 58, 55 - caries - against narrative
                                "Gemella","Granulicatella", # 72 - caries - against narrative
                                "Treponema","Filifactor","Dialister", # 47, 33 - periodontitis
                                "Prevotella","Veillonella","Solobacterium","Atopobium", # 92, 61, 50 - periodontitis - look into more for the other 3
                                # "Bergeyella","Neisseria", # 39 - low Berg in CF, and its assoc with low caries, Neisseria helps cariogens, so their assoc in mC makes sense
                                "Actinomyces","Stomatobaculum","Rothia","Atopobium", # 62, 58, 46 - caries
                                "Ralstonia","Variovorax", # 96 - AHL-acylase, periodontitis
                                "Bradyrhizobium","Hyphomicrobium" # 86 - periodontitis?
                              ))


# **************************************************************************** #





ggarrange(plot(CF.net$ig.no_weight, layout = CF.net$coords.fr,
               vertex.size = CF.net$vsize[ V(CF.net$ig.no_weight)$name ],
               edge.width = abs(E(CF.net$ig)$weight)*30,
               vertex.label.cex = 0.85,
               main = CF.net$net.title), 
          plot(mC.net$ig.no_weight, layout = mC.net$coords.fr,
               vertex.size = mC.net$vsize[ V(mC.net$ig.no_weight)$name ],
               edge.width = abs(E(mC.net$ig)$weight)*30,
               vertex.label.cex = 0.85,
               main = mC.net$net.title), 
          hamming_boxes, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))


# **************************************************************************** #
if (addLabels == FALSE) 
  netPlot <- plot(ig, layout=coords.fr, 
                  vertex.size = vs.plot, edge.width = ew.plot, 
                  vertex.label.cex = 0.001, 
                  main=net.title)

else if (addLabels == TRUE) 
  netPlot <- plot(ig.no_weight.taxLabs, layout=coords.fr, 
                  vertex.size = vs.plot, edge.width = ew.plot, 
                  vertex.label.cex = vertLabCex, 
                  main=net.title)

if (Vprop == TRUE) vs.plot <- vsize[ V(ig.no_weight)$name ]

if (Eprop == TRUE) ew.plot <- (abs(E(ig)$weight))*30














# **************************************************************************** #


n.c.CF_Yes <- symBeta(getOptBeta( se.mb.sll2.CF.Yes ))
betaMat.CF_Yes <- as.matrix( n.c.CF_Yes )

# to add abundance values to vertices
colnames(n.c.CF_Yes) <- rownames(n.c.CF_Yes) <- taxa_names(tagl.CF.Yes)#rownames(gloms$Genus)
# add log abundance as properties of vertex/nodes.
# vsize.CF_Yes <- log2(apply(gloms$Genus[ , samps.CF ], 1, mean)+1)
vsize.CF_Yes <- log2(apply(tagl.CF.Yes@otu_table, 1, mean)+1)

# get number of positive and negative edges
# We divide by two since an edge is represented by two entries in the matrix. 
pos_edge.CF_Yes <- length(betaMat.CF_Yes[ betaMat.CF_Yes > 0 ])/2
neg_edge.CF_Yes <- length(betaMat.CF_Yes[ betaMat.CF_Yes < 0 ])/2
tot_edge.CF_Yes <- length(betaMat.CF_Yes[ betaMat.CF_Yes != 0 ])/2



# prepare data for plotting
ig.CF_Yes <- graph.adjacency(n.c.CF_Yes, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# remove those nodes with no edges
ig.CF_Yes <- delete.vertices(ig.CF_Yes, V(ig.CF_Yes)[degree(ig.CF_Yes) == 0])
ig.CF_Yes <- delete.vertices(ig.CF_Yes, V(ig.CF_Yes)[name == "0"])
# # we can see all the attributes and weights
ig.CF_Yes

# one more object, to be able to check weights and genus names together
ig.CF_Yes.unchanged <- ig.CF_Yes
E(ig.CF_Yes.unchanged)[weight < 0]


# now color the edges based on their values positive is steelblue
# E(igraphs.CF.Yes.noEmpty)[weight > 0]$color<-"steelblue"
# # now color the edges based on their values
# E(igraphs.CF.Yes.noEmpty)[weight < 0]$color<-"orange"
E(ig.CF_Yes)[weight > 0]$color <- "#73BDD3" #"steelblue"
# now color the edges based on their values
E(ig.CF_Yes)[weight < 0]$color <- "#CB6767" # "orange"


# before making coordinates, must remove weights
ig.CF_Yes.no_weight <- ig.CF_Yes
E(ig.CF_Yes.no_weight)$weight <- 1

# make coordinates
set.seed(100)
coords.fr.CF_Yes <- layout_with_fr(ig.CF_Yes.no_weight)#igraphs.CF.Yes.noEmpty)

# color dots by phylum
# prepare colors for plot
colPal.cols <- c("#e38f87","#9abddf","#90C98A",
                 "#FFBC67","#a854a8","#ab8b67",
                 "darkblue","#bbbbbb","pink","#56423c",
                 "brown","black","gray","darkorange3")
names(colPal.cols) <- c("Firmicutes","Bacteroidetes","Proteobacteria",
                        "Fusobacteria","Actinobacteria","Epsilonbacteraeota",
                        "Patescibacteria","Spirochaetes","Synergistetes","unclassified.P1",
                        "Verrucomicrobia","Tenericutes","Cyanobacteria","Chloroflexi")

V(ig.CF_Yes)$name <- unname(taxTables.both$Genus[V(ig.CF_Yes)$name,"Phylum"])
# V(ig.CF_Yes)$color <- match(V(ig.CF_Yes)$name, unique(V(ig.CF_Yes)$name))
V(ig.CF_Yes)$color <- colPal.cols[ V(ig.CF_Yes)$name ]

# give same colors to no_weight graph, just in case
E(ig.CF_Yes.no_weight)$color <- E(ig.CF_Yes)$color
V(ig.CF_Yes.no_weight)$color <- V(ig.CF_Yes)$color



# (1) plot with no labels or proportional anything
plot(ig.CF_Yes, layout=coords.fr.CF_Yes, vertex.size = 6, 
     vertex.label.cex = 0.001, #edge.width = exp(abs(E(ig.CF_Yes)$weight))*3, 
     main="Co-occurrence network: CF samples")

# (2) plot with proportional vertex sizes
plot(ig.CF_Yes, layout=coords.fr.CF_Yes, vertex.size = vsize.CF_Yes[ V(ig.CF_Yes.no_weight)$name ], 
     vertex.label.cex = 0.001, #edge.width = exp(E(ig.CF_Yes)$weight),
     main="Co-occurrence network: CF samples")

# (3) plot with proportional edge widths
plot(ig.CF_Yes, layout=coords.fr.CF_Yes, vertex.size = 6, 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.CF_Yes)$weight))*30, 
     main="Co-occurrence network: CF samples")

# (4) plot with proportional edge widths and vertex sizes
plot(ig.CF_Yes, layout=coords.fr.CF_Yes, vertex.size = vsize.CF_Yes[ V(ig.CF_Yes.no_weight)$name ], 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.CF_Yes)$weight))*30, 
     main="Co-occurrence network: CF samples")

# (5) plot with genus labels, proportional edge widths and vertex sizes
plot(ig.CF_Yes.no_weight, layout=coords.fr.CF_Yes, vertex.size = vsize.CF_Yes[ V(ig.CF_Yes.no_weight)$name ], 
     vertex.label.cex = 0.75, edge.width = (abs(E(ig.CF_Yes)$weight))*30, 
     main="Co-occurrence network: CF samples")

# (6) plot with genus labels
plot(ig.CF_Yes.no_weight, layout=coords.fr.CF_Yes, vertex.size = 6, 
     vertex.label.cex = 0.75, #edge.width = exp(abs(E(ig.CF_Yes)$weight))*3, 
     main="Co-occurrence network: CF samples")

# plot_network(igraphs.CF.Yes.noEmpty, tagl.CF.Yes, type='taxa', color="Phylum", title = "Co-occurrence network: CF samples")

legend(x=-1.83, y=1, inset=0.2, title="Phylum", 
       gsub("unclassified.P1", "unclassified", sort(unique(V(ig.CF_Yes)$name))),
       fill = colPal.cols[ sort(unique(V(ig.CF_Yes)$name)) ], 
       cex=1.2, bty="n")
# add legend
********


# *************************************** #
# c("#CB6767","#73BDD3")
# c("#bf4342","#2f6690")


n.c.CF_all_No <- symBeta(getOptBeta( se.mb.sll2.CF.all_No ))
betaMat.CF_all_No <- as.matrix( n.c.CF_all_No )

# to add abundance values to vertices
colnames(n.c.CF_all_No) <- rownames(n.c.CF_all_No) <- taxa_names(tagl.CF.all_No)#rownames(gloms$Genus)
# add log abundance as properties of vertex/nodes.
# vsize.CF_all_No <- log2(apply(gloms$Genus[ , samps.mC ], 1, mean)+1)
vsize.CF_all_No <- log2(apply(tagl.CF.all_No@otu_table, 1, mean)+1)


# get number of positive and negative edges
# We divide by two since an edge is represented by two entries in the matrix. 
pos_edge.CF_all_No <- length(betaMat.CF_all_No[ betaMat.CF_all_No > 0 ])/2
neg_edge.CF_all_No <- length(betaMat.CF_all_No[ betaMat.CF_all_No < 0 ])/2
tot_edge.CF_all_No <- length(betaMat.CF_all_No[ betaMat.CF_all_No != 0 ])/2




# prepare data for plotting
ig.CF_all_No <- graph.adjacency(n.c.CF_all_No, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# remove those nodes with no edges
ig.CF_all_No <- delete.vertices(ig.CF_all_No, V(ig.CF_all_No)[degree(ig.CF_all_No) == 0])
ig.CF_all_No <- delete.vertices(ig.CF_all_No, V(ig.CF_all_No)[name == "0"])
# we can see all the attributes and weights
ig.CF_all_No

# one more object, to be able to check weights and genus names together
ig.CF_all_No.unchanged <- ig.CF_all_No
E(ig.CF_all_No.unchanged)[weight < 0]

# now color the edges based on their values positive is steelblue
# E(ig.CF_all_No)[weight > 0]$color <- "#CB6767" #"steelblue"
# # now color the edges based on their values
# E(ig.CF_all_No)[weight < 0]$color <- "#73BDD3" #"orange"
E(ig.CF_all_No)[weight > 0]$color <- "#73BDD3" #"steelblue"
# now color the edges based on their values
E(ig.CF_all_No)[weight < 0]$color <- "#CB6767" # "orange"



# before making coordinates, must remove weights
ig.CF_all_No.no_weight <- ig.CF_all_No
E(ig.CF_all_No.no_weight)$weight <- 1

# make coordinates
set.seed(104)
coords.fr.CF_all_No <- layout_with_fr(ig.CF_all_No.no_weight)



# color dots by phylum
# prepare colors for plot
colPal.cols <- c("#e38f87","#9abddf","#90C98A",
                 "#FFBC67","#a854a8","#ab8b67",
                 "darkblue","#bbbbbb","pink","#56423c",
                 "brown","black","gray","darkorange3")
names(colPal.cols) <- c("Firmicutes","Bacteroidetes","Proteobacteria",
                        "Fusobacteria","Actinobacteria","Epsilonbacteraeota",
                        "Patescibacteria","Spirochaetes","Synergistetes","unclassified.P1",
                        "Verrucomicrobia","Tenericutes","Cyanobacteria","Chloroflexi")

V(ig.CF_all_No)$name <- unname(taxTables.both$Genus[V(ig.CF_all_No)$name,"Phylum"])
# V(ig.CF_all_No)$color <- match(V(ig.CF_all_No)$name, unique(V(ig.CF_all_No)$name))
V(ig.CF_all_No)$color <- colPal.cols[ V(ig.CF_all_No)$name ]

# give same colors to no_weight graph, just in case
E(ig.CF_all_No.no_weight)$color <- E(ig.CF_all_No)$color
V(ig.CF_all_No.no_weight)$color <- V(ig.CF_all_No)$color



# (1) plot with no labels or proportional anything
plot(ig.CF_all_No, layout=coords.fr.CF_all_No, vertex.size = 6, 
     vertex.label.cex = 0.001, main="Co-occurrence network: matched control samples")

# (2) plot with proportional vertex sizes
plot(ig.CF_all_No, layout=coords.fr.CF_all_No, vertex.size = vsize.CF_all_No[ V(ig.CF_all_No.no_weight)$name ], 
     vertex.label.cex = 0.001, main="Co-occurrence network: matched control samples")

# (3) plot with proportional edge widths
plot(ig.CF_all_No, layout=coords.fr.CF_all_No, vertex.size = 6, 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.CF_all_No)$weight))*30, 
     main="Co-occurrence network: matched control samples")

# (4) plot with proportional edge widths and vertex sizes
plot(ig.CF_all_No, layout=coords.fr.CF_all_No, vertex.size = vsize.CF_all_No[ V(ig.CF_all_No.no_weight)$name ], 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.CF_all_No)$weight))*30, 
     main="Co-occurrence network: matched control samples")


# (5) plot with genus labels, proportional edge widths and vertex sizes
plot(ig.CF_all_No.no_weight, layout=coords.fr.CF_all_No, vertex.size = vsize.CF_all_No[ V(ig.CF_all_No.no_weight)$name ], 
     vertex.label.cex = 0.75, edge.width = (abs(E(ig.CF_all_No)$weight))*30, 
     main="Co-occurrence network: matched control samples")

# (6) plot with genus labels
plot(ig.CF_all_No.no_weight, layout=coords.fr.CF_all_No, vertex.size = 6, 
     vertex.label.cex = 0.75, #edge.width = exp(abs(E(ig.CF_all_No)$weight))*3, 
     main="Co-occurrence network: matched control samples")


legend(x=-1.83, y=1, inset=0.2, title="Phylum", 
       gsub("unclassified.P1", "unclassified", sort(unique(V(ig.CF_all_No)$name))),
       fill = colPal.cols[ sort(unique(V(ig.CF_all_No)$name)) ], 
       cex=1.2, bty="n")









# *************************************** #
# c("#CB6767","#73BDD3")
# c("#bf4342","#2f6690")

i <- "1"
n.c.mC_1 <- symBeta(getOptBeta( SE_objects[[ i ]] ))
betaMat.mC_1 <- as.matrix( n.c.mC_1 )

samps.mC_1 <- disSubsTests$default$Cystic_fibrosis[[ i ]]$samples
samps.mC_1 <- samps.mC_1[ SLL2.meta[samps.mC_1, "Cystic_fibrosis"]=="No" ]
tagl.mC_1  <- subset_samples(tagl, sample_names(tagl) %in% samps.mC_1)

# to add abundance values to vertices
colnames(n.c.mC_1) <- rownames(n.c.mC_1) <- taxa_names(tagl.mC_1)#rownames(gloms$Genus)
# add log abundance as properties of vertex/nodes.
# vsize.mC_1 <- log2(apply(gloms$Genus[ , samps.mC_1 ], 1, mean)+1)
vsize.mC_1 <- log2(apply(tagl.mC_1@otu_table, 1, mean)+1)

# get number of positive and negative edges
# We divide by two since an edge is represented by two entries in the matrix. 
pos_edge.mC_1 <- length(betaMat.mC_1[ betaMat.mC_1 > 0 ])/2
neg_edge.mC_1 <- length(betaMat.mC_1[ betaMat.mC_1 < 0 ])/2
tot_edge.mC_1 <- length(betaMat.mC_1[ betaMat.mC_1 != 0 ])/2




# prepare data for plotting
ig.mC_1 <- graph.adjacency(n.c.mC_1, mode='undirected', add.rownames = TRUE, weighted = TRUE)
# remove those nodes with no edges
ig.mC_1 <- delete.vertices(ig.mC_1, V(ig.mC_1)[degree(ig.mC_1) == 0])
ig.mC_1 <- delete.vertices(ig.mC_1, V(ig.mC_1)[name == "0"])
# we can see all the attributes and weights
ig.mC_1

# one more object, to be able to check weights and genus names together
ig.mC_1.unchanged <- ig.mC_1
E(ig.mC_1.unchanged)[weight < 0]

# now color the edges based on their values positive is steelblue
# E(ig.mC_1)[weight > 0]$color <- "#CB6767" #"steelblue"
# # now color the edges based on their values
# E(ig.mC_1)[weight < 0]$color <- "#73BDD3" #"orange"
E(ig.mC_1)[weight > 0]$color <- "#73BDD3" #"steelblue"
# now color the edges based on their values
E(ig.mC_1)[weight < 0]$color <- "#CB6767" # "orange"



# before making coordinates, must remove weights
ig.mC_1.no_weight <- ig.mC_1
E(ig.mC_1.no_weight)$weight <- 1
coords.fr.mC_1 <- layout_with_fr(ig.mC_1.no_weight)



# color dots by phylum
# prepare colors for plot
colPal.cols <- c("#e38f87","#9abddf","#90C98A",
                 "#FFBC67","#a854a8","#ab8b67",
                 "darkblue","#bbbbbb","pink","#56423c",
                 "brown","black","gray","darkorange3")
names(colPal.cols) <- c("Firmicutes","Bacteroidetes","Proteobacteria",
                        "Fusobacteria","Actinobacteria","Epsilonbacteraeota",
                        "Patescibacteria","Spirochaetes","Synergistetes","unclassified.P1",
                        "Verrucomicrobia","Tenericutes","Cyanobacteria","Chloroflexi")

V(ig.mC_1)$name <- unname(taxTables.both$Genus[V(ig.mC_1)$name,"Phylum"])
# V(ig.mC_1)$color <- match(V(ig.mC_1)$name, unique(V(ig.mC_1)$name))
V(ig.mC_1)$color <- colPal.cols[ V(ig.mC_1)$name ]

# give same colors to no_weight graph, just in case
E(ig.mC_1.no_weight)$color <- E(ig.mC_1)$color
V(ig.mC_1.no_weight)$color <- V(ig.mC_1)$color


# (1) plot with no labels or proportional anything
plot(ig.mC_1, layout=coords.fr.mC_1, vertex.size = 6, 
     vertex.label.cex = 0.001, main="Co-occurrence network: matched controls - subsampling #1")

# (2) plot with proportional vertex sizes
plot(ig.mC_1, layout=coords.fr.mC_1, vertex.size = vsize.mC_1[ V(ig.mC_1.no_weight)$name ], 
     vertex.label.cex = 0.001, main="Co-occurrence network: matched controls - subsampling #1")

# (3) plot with proportional edge widths
plot(ig.mC_1, layout=coords.fr.mC_1, vertex.size = 6, 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.mC_1)$weight))*30, 
     main="Co-occurrence network: matched controls - subsampling #1")

# (4) plot with proportional edge widths and vertex sizes
plot(ig.mC_1, layout=coords.fr.mC_1, vertex.size = vsize.mC_1[ V(ig.mC_1.no_weight)$name ], 
     vertex.label.cex = 0.001, edge.width = (abs(E(ig.mC_1)$weight))*30, 
     main="Co-occurrence network: matched controls - subsampling #1")

# (5) plot with genus labels, proportional edge widths and vertex sizes
plot(ig.mC_1.no_weight, layout=coords.fr.mC_1, vertex.size = vsize.mC_1[ V(ig.mC_1.no_weight)$name ], 
     vertex.label.cex = 0.75, edge.width = (abs(E(ig.mC_1)$weight))*30, 
     main="Co-occurrence network: matched controls - subsampling #1")


legend(x=-1.83, y=1, inset=0.2, title="Phylum", 
       gsub("unclassified.P1", "unclassified", sort(unique(V(ig.mC_1)$name))),
       fill = colPal.cols[ sort(unique(V(ig.mC_1)$name)) ], 
       cex=1.2, bty="n")







# # for tkplot
# spiec.graph.b <- ig.CF_all_No
# nodenames <- V(spiec.graph.b)$name
# # V(spiec.graph.b)$name <- getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
# E(spiec.graph.b)$arrow.size=5
# V(spiec.graph.b)$color="white"
# V(spiec.graph.b)$frame.color="black"
# tkplot(spiec.graph.b)

# ****************************************************************************************************************** ####















# ****************************************************************************************************************** ####

# ************************************************** #

# Networks for Age_groups ####

library(nettools)
library(SpiecEasi)

tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))
# keep only those samples that do not have any chronic disorder
tagl_healthy <- subset_samples(tagl, Chronic_disorder=="No")

# first get objects for each age group

# Child
samps.Child <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"]=="Child" &
                                     SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Child  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Child)

# Teen
samps.Teen <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"]=="Teen" &
                                    SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Teen  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Teen)

# Adult
samps.Adult <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"]=="Adult" &
                                     SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Adult  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Adult)

# Senior
samps.Senior <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                    SLL2.meta[, "Age_groups"]=="Senior" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Senior  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Senior)

# # remove those taxa that dont have at least 100 counts across all CF + mC samples,
# #   then will make tagl.CF.Yes and talg.CF.all_No again below
# tagl <- prune_taxa(taxa_sums(tagl.CF.all_No) + taxa_sums(tagl.CF.Yes) > 100, tagl)

# ************************** #
# from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
# filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
#    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
library(seqtime)
minCount <- 15
minOcc   <- 20
g10 <- as.data.frame(gloms$Genus[ , c(samps.Child, samps.Teen, samps.Adult, samps.Senior) ])
g10[ g10 <= minCount ] <- 0
# filterobj <- filterTaxonMatrix(gloms$Genus[ , c(samps.Child, samps.Teen, samps.Adult, samps.Senior) ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
otus.f <- filterobj$mat

# replace the values that were made to be 0s for the filtering
anti.g10 <- as.data.frame(gloms$Genus[ , c(samps.Child, samps.Teen, samps.Adult, samps.Senior) ])
anti.g10[ anti.g10 > minCount ] <- 0
# replace them into g10
g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
# first update the values that were made 0
upd.otus_f <- otus.f
upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(gloms$Genus[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], c(samps.Child, samps.Teen, samps.Adult, samps.Senior) ])
# then update the summed row from the removed values
upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))

# finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
taxa.f <- tax_table(tagl_healthy)[setdiff(1:nrow(tax_table(tagl_healthy)),filterobj$filtered.indices),]
dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
taxa.f <- rbind(taxa.f, dummyTaxonomy)
rownames(taxa.f)[nrow(taxa.f)] <- "0"
rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"

tagl_healthy <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                 sample_data(SLL2)[ c(samps.Child, samps.Teen, samps.Adult, samps.Senior), ],
                 tax_table(taxa.f))
# ************************** #


# first get objects for each age group

# Child
samps.Child <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"]=="Child" &
                                     SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Child  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Child)
# se.mb.sll2.Child <- spiec.easi(tagl.Child, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Child, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Child.rds", p2_dir))
se.mb.sll2.Child <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Child.rds", p2_dir))
igraphs.Child.allNodes <- adj2igraph(getRefit(se.mb.sll2.Child),  vertex.attr=list(name=taxa_names(tagl.Child)), rmEmptyNodes = F)
igraphs.Child.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Child),  vertex.attr=list(name=taxa_names(tagl.Child)), rmEmptyNodes = T)



# Teen
samps.Teen <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                    SLL2.meta[, "Age_groups"]=="Teen" &
                                    SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Teen  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Teen)
# se.mb.sll2.Teen <- spiec.easi(tagl.Teen, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Teen, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Teen.rds", p2_dir))
se.mb.sll2.Teen <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Teen.rds", p2_dir))
igraphs.Teen.allNodes <- adj2igraph(getRefit(se.mb.sll2.Teen),  vertex.attr=list(name=taxa_names(tagl.Teen)), rmEmptyNodes = F)
igraphs.Teen.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Teen),  vertex.attr=list(name=taxa_names(tagl.Teen)), rmEmptyNodes = T)



# Adult
samps.Adult <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"]=="Adult" &
                                     SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Adult  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Adult)
# se.mb.sll2.Adult <- spiec.easi(tagl.Adult, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Adult, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Adult.rds", p2_dir))
se.mb.sll2.Adult <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Adult.rds", p2_dir))
igraphs.Adult.allNodes <- adj2igraph(getRefit(se.mb.sll2.Adult),  vertex.attr=list(name=taxa_names(tagl.Adult)), rmEmptyNodes = F)
igraphs.Adult.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Adult),  vertex.attr=list(name=taxa_names(tagl.Adult)), rmEmptyNodes = T)


# Senior
samps.Senior <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                      SLL2.meta[, "Age_groups"]=="Senior" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Senior  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Senior)
# se.mb.sll2.Senior <- spiec.easi(tagl.Senior, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Senior, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Senior.rds", p2_dir))
se.mb.sll2.Senior <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Senior.rds", p2_dir))
igraphs.Senior.allNodes <- adj2igraph(getRefit(se.mb.sll2.Senior),  vertex.attr=list(name=taxa_names(tagl.Senior)), rmEmptyNodes = F)
igraphs.Senior.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Senior),  vertex.attr=list(name=taxa_names(tagl.Senior)), rmEmptyNodes = T)



# Younger - both Child and Teen samples
samps.Younger <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                       SLL2.meta[, "Age_groups"] %in% c("Child","Teen") &
                                       SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Younger  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Younger)
# se.mb.sll2.Younger <- spiec.easi(tagl.Younger, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Younger, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Younger.rds", p2_dir))
se.mb.sll2.Younger <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Younger.rds", p2_dir))
igraphs.Younger.allNodes <- adj2igraph(getRefit(se.mb.sll2.Younger),  vertex.attr=list(name=taxa_names(tagl.Younger)), rmEmptyNodes = F)
igraphs.Younger.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Younger),  vertex.attr=list(name=taxa_names(tagl.Younger)), rmEmptyNodes = T)



# Older - both Adult and Senior samples
samps.Older <- rownames(SLL2.meta[ ! is.na(SLL2.meta[, "Age_groups"]) &
                                     SLL2.meta[, "Age_groups"] %in% c("Adult","Senior") &
                                     SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Older  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Older)
# se.mb.sll2.Older <- spiec.easi(tagl.Older, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Older, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Older.rds", p2_dir))
se.mb.sll2.Older <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Older.rds", p2_dir))
igraphs.Older.allNodes <- adj2igraph(getRefit(se.mb.sll2.Older),  vertex.attr=list(name=taxa_names(tagl.Older)), rmEmptyNodes = F)
igraphs.Older.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Older),  vertex.attr=list(name=taxa_names(tagl.Older)), rmEmptyNodes = T)



# ************************************************** #
# ************************************************** #

# use package intergraph to convert igraphs to networks, from which I can pull out all the edge connections ####
library(intergraph)

edge.mat.Child   <- as.matrix(asNetwork(igraphs.Child.allNodes))
edge.mat.Teen    <- as.matrix(asNetwork(igraphs.Teen.allNodes))
edge.mat.Adult   <- as.matrix(asNetwork(igraphs.Adult.allNodes))
edge.mat.Senior  <- as.matrix(asNetwork(igraphs.Senior.allNodes))
edge.mat.Younger <- as.matrix(asNetwork(igraphs.Adult.allNodes))
edge.mat.Older   <- as.matrix(asNetwork(igraphs.Senior.allNodes))

# for each genus, get vector of those other genera with which it has an edge
gen_edges.Child   <- list()
gen_edges.Teen    <- list()
gen_edges.Adult   <- list()
gen_edges.Senior  <- list()
gen_edges.Younger <- list()
gen_edges.Older   <- list()
for (gen in rownames(edge.mat.Child)) {
  
  gen_edges.Child[[ gen ]]   <- colnames(edge.mat.Child)[ edge.mat.Child[ gen, ] != 0 ]
  gen_edges.Teen[[ gen ]]    <- colnames(edge.mat.Teen)[ edge.mat.Teen[ gen, ] != 0 ]
  gen_edges.Adult[[ gen ]]   <- colnames(edge.mat.Adult)[ edge.mat.Adult[ gen, ] != 0 ]
  gen_edges.Senior[[ gen ]]  <- colnames(edge.mat.Senior)[ edge.mat.Senior[ gen, ] != 0 ]
  gen_edges.Younger[[ gen ]] <- colnames(edge.mat.Younger)[ edge.mat.Younger[ gen, ] != 0 ]
  gen_edges.Older[[ gen ]]   <- colnames(edge.mat.Older)[ edge.mat.Older[ gen, ] != 0 ]
  
}





# ******************** #
# first test with 1st subsamp 
#   - for each genus, get vector of connected genera that are same in both CF and mC
#   - vector of genera only connected in CF
#   - vector of genera only connected in mC

gen_connect.age_groups <- list()

gen_connect.age_groups[[ "All" ]] <- sapply(names(gen_edges.Child), function(gen) 
  gen_edges.Child[[ gen ]][ gen_edges.Child[[ gen ]] %in% c(gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ] )

gen_connect.age_groups[[ "Child" ]] <- sapply(names(gen_edges.Child), function(gen) 
  gen_edges.Child[[ gen ]][ ! gen_edges.Child[[ gen ]] %in% c(gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Teen" ]] <- sapply(names(gen_edges.Teen), function(gen) 
  gen_edges.Teen[[ gen ]][ ! gen_edges.Teen[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Adult" ]] <- sapply(names(gen_edges.Adult), function(gen) 
  gen_edges.Adult[[ gen ]][ ! gen_edges.Adult[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Teen[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Senior" ]] <- sapply(names(gen_edges.Senior), function(gen) 
  gen_edges.Senior[[ gen ]][ ! gen_edges.Senior[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]]) ])

gen_connect.age_groups[[ "Younger" ]] <- sapply(names(gen_edges.Younger), function(gen) 
  gen_edges.Younger[[ gen ]][ ! gen_edges.Younger[[ gen ]] %in% gen_edges.Older[[ gen ]] ])

gen_connect.age_groups[[ "Older" ]] <- sapply(names(gen_edges.Older), function(gen) 
  gen_edges.Older[[ gen ]][ ! gen_edges.Older[[ gen ]] %in% gen_edges.Younger[[ gen ]] ])


# *************************************** #
gen_con_freqs.age_groups <- function(gc_list, gen) {
  
  cat(sprintf("**** %s connections in All ****\n", gen))
  print(sort(gc_list$All[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Child only ****\n", gen))
  print(sort(gc_list$Child[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Teen only ****\n", gen))
  print(sort(gc_list$Teen[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Adult only ****\n", gen))
  print(sort(gc_list$Adult[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Senior only ****\n", gen))
  print(sort(gc_list$Senior[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Younger only ****\n", gen))
  print(sort(gc_list$Younger[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Older only ****\n", gen))
  print(sort(gc_list$Older[[ gen ]]))
  
}
# *************************************** #


gen_con_freqs.age_groups(gen_connect.age_groups, "Porphyromonas")




# ****************************************************************************************************************** ####
































# ****************************************************************************************************************** ####

# ************************************************** #

# Networks for Community ####

library(nettools)
library(SpiecEasi)

tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))
# keep only those samples that do not have any chronic disorder
tagl_healthy <- subset_samples(tagl, Chronic_disorder=="No")

# first get objects for each age group


# Andalucia
samps.Andalucia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Andalucía" &
                                         SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Andalucia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Andalucia)

# Aragon
samps.Aragon <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Aragón" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Aragon  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Aragon)

# Cantabria
samps.Cantabria <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Cantabria" &
                                         SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Cantabria  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Cantabria)

# Cataluña
samps.Cataluna <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Cataluña" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Cataluna  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Cataluna)

# Madrid
samps.Madrid <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Comunidad de Madrid" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Madrid  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Madrid)

# Valencia
samps.Valencia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Comunidad Valenciana" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Valencia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Valencia)

# Galicia
samps.Galicia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Galicia" &
                                       SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Galicia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Galicia)

# Baleares
samps.Baleares <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Islas Baleares" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Baleares  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Baleares)

# Murcia
samps.Murcia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Murcia" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Murcia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Murcia)

# PV
samps.PV <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="País Vasco" &
                                  SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.PV  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.PV)

# # remove those taxa that dont have at least 100 counts across all CF + mC samples,
# #   then will make tagl.CF.Yes and talg.CF.all_No again below
# tagl <- prune_taxa(taxa_sums(tagl.CF.all_No) + taxa_sums(tagl.CF.Yes) > 100, tagl)

# ************************** #
# from this tutorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
# filter out taxa that dont have at least 15 reads in at least 20 of the CF + mC samples, but with this function,
#    keeps a row of all the counts that were removed, so it can keep the total counts for each sample
library(seqtime)
minCount <- 15
minOcc   <- 20
g10 <- as.data.frame(gloms$Genus[ , sample_names(tagl_healthy) ])
g10[ g10 <= minCount ] <- 0
# filterobj <- filterTaxonMatrix(gloms$Genus[ , sample_names(tagl_healthy) ], minocc=100, keepSum = TRUE, return.filtered.indices = TRUE)
filterobj <- filterTaxonMatrix(g10, minocc=minOcc, keepSum = TRUE, return.filtered.indices = TRUE)
otus.f <- filterobj$mat

# replace the values that were made to be 0s for the filtering
anti.g10 <- as.data.frame(gloms$Genus[ , sample_names(tagl_healthy) ])
anti.g10[ anti.g10 > minCount ] <- 0
# replace them into g10
g10[ anti.g10 != 0 ] <- anti.g10[ anti.g10 != 0 ]
# first update the values that were made 0
upd.otus_f <- otus.f
upd.otus_f[ 1:(nrow(upd.otus_f)-1), ] <- as.data.frame(gloms$Genus[ rownames(upd.otus_f)[ 1:(nrow(upd.otus_f)-1) ], sample_names(tagl_healthy) ])
# then update the summed row from the removed values
upd.otus_f[ nrow(upd.otus_f), ] <- as.numeric(upd.otus_f[ nrow(upd.otus_f), ] + (colSums(upd.otus_f) - colSums(otus.f)))

# finally update all the names appropriately, including pseudo-name for the summed counts, which will be ignored later in the plots
taxa.f <- tax_table(tagl_healthy)[setdiff(1:nrow(tax_table(tagl_healthy)),filterobj$filtered.indices),]
dummyTaxonomy <- c("k__dummy","p__","c__","o__","f__","g__")
taxa.f <- rbind(taxa.f, dummyTaxonomy)
rownames(taxa.f)[nrow(taxa.f)] <- "0"
rownames(upd.otus_f)[nrow(upd.otus_f)] <- "0"

tagl_healthy <- phyloseq(otu_table(upd.otus_f, taxa_are_rows = T),
                 sample_data(SLL2)[ sample_names(tagl_healthy), ],
                 tax_table(taxa.f))
# ************************** #



# first get objects for each community

# Andalucia
samps.Andalucia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Andalucía" &
                                         SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Andalucia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Andalucia)
# se.mb.sll2.Andalucia <- spiec.easi(tagl.Andalucia, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Andalucia, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Andalucia.rds", p2_dir))
se.mb.sll2.Andalucia <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Andalucia.rds", p2_dir))
igraphs.Andalucia.allNodes <- adj2igraph(getRefit(se.mb.sll2.Andalucia),  vertex.attr=list(name=taxa_names(tagl.Andalucia)), rmEmptyNodes = F)
igraphs.Andalucia.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Andalucia),  vertex.attr=list(name=taxa_names(tagl.Andalucia)), rmEmptyNodes = T)


# Aragon
samps.Aragon <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Aragón" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Aragon  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Aragon)
# se.mb.sll2.Aragon <- spiec.easi(tagl.Aragon, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Aragon, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Aragon.rds", p2_dir))
se.mb.sll2.Aragon <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Aragon.rds", p2_dir))
igraphs.Aragon.allNodes <- adj2igraph(getRefit(se.mb.sll2.Aragon),  vertex.attr=list(name=taxa_names(tagl.Aragon)), rmEmptyNodes = F)
igraphs.Aragon.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Aragon),  vertex.attr=list(name=taxa_names(tagl.Aragon)), rmEmptyNodes = T)


# Cantabria
samps.Cantabria <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Cantabria" &
                                         SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Cantabria  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Cantabria)
# se.mb.sll2.Cantabria <- spiec.easi(tagl.Cantabria, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
# saveRDS(se.mb.sll2.Cantabria, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Cantabria.rds", p2_dir))
se.mb.sll2.Cantabria <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Cantabria.rds", p2_dir))
igraphs.Cantabria.allNodes <- adj2igraph(getRefit(se.mb.sll2.Cantabria),  vertex.attr=list(name=taxa_names(tagl.Cantabria)), rmEmptyNodes = F)
igraphs.Cantabria.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Cantabria),  vertex.attr=list(name=taxa_names(tagl.Cantabria)), rmEmptyNodes = T)


# Cataluna
samps.Cataluna <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Cataluña" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Cataluna  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Cataluna)
se.mb.sll2.Cataluna <- spiec.easi(tagl.Cataluna, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Cataluna, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Cataluna.rds", p2_dir))
se.mb.sll2.Cataluna <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Cataluna.rds", p2_dir))
igraphs.Cataluna.allNodes <- adj2igraph(getRefit(se.mb.sll2.Cataluna),  vertex.attr=list(name=taxa_names(tagl.Cataluna)), rmEmptyNodes = F)
igraphs.Cataluna.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Cataluna),  vertex.attr=list(name=taxa_names(tagl.Cataluna)), rmEmptyNodes = T)


# Madrid
samps.Madrid <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Comunidad de Madrid" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Madrid  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Madrid)
se.mb.sll2.Madrid <- spiec.easi(tagl.Madrid, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Madrid, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Madrid.rds", p2_dir))
se.mb.sll2.Madrid <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Madrid.rds", p2_dir))
igraphs.Madrid.allNodes <- adj2igraph(getRefit(se.mb.sll2.Madrid),  vertex.attr=list(name=taxa_names(tagl.Madrid)), rmEmptyNodes = F)
igraphs.Madrid.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Madrid),  vertex.attr=list(name=taxa_names(tagl.Madrid)), rmEmptyNodes = T)


# Valencia
samps.Valencia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Comunidad Valenciana" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Valencia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Valencia)
se.mb.sll2.Valencia <- spiec.easi(tagl.Valencia, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Valencia, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Valencia.rds", p2_dir))
se.mb.sll2.Valencia <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Valencia.rds", p2_dir))
igraphs.Valencia.allNodes <- adj2igraph(getRefit(se.mb.sll2.Valencia),  vertex.attr=list(name=taxa_names(tagl.Valencia)), rmEmptyNodes = F)
igraphs.Valencia.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Valencia),  vertex.attr=list(name=taxa_names(tagl.Valencia)), rmEmptyNodes = T)


# Galicia
samps.Galicia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Galicia" &
                                       SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Galicia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Galicia)
se.mb.sll2.Galicia <- spiec.easi(tagl.Galicia, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Galicia, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Galicia.rds", p2_dir))
se.mb.sll2.Galicia <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Galicia.rds", p2_dir))
igraphs.Galicia.allNodes <- adj2igraph(getRefit(se.mb.sll2.Galicia),  vertex.attr=list(name=taxa_names(tagl.Galicia)), rmEmptyNodes = F)
igraphs.Galicia.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Galicia),  vertex.attr=list(name=taxa_names(tagl.Galicia)), rmEmptyNodes = T)


# Baleares
samps.Baleares <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Islas Baleares" &
                                        SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Baleares  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Baleares)
se.mb.sll2.Baleares <- spiec.easi(tagl.Baleares, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Baleares, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Baleares.rds", p2_dir))
se.mb.sll2.Baleares <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Baleares.rds", p2_dir))
igraphs.Baleares.allNodes <- adj2igraph(getRefit(se.mb.sll2.Baleares),  vertex.attr=list(name=taxa_names(tagl.Baleares)), rmEmptyNodes = F)
igraphs.Baleares.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Baleares),  vertex.attr=list(name=taxa_names(tagl.Baleares)), rmEmptyNodes = T)


# Murcia
samps.Murcia <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="Murcia" &
                                      SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.Murcia  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.Murcia)
se.mb.sll2.Murcia <- spiec.easi(tagl.Murcia, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.Murcia, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Murcia.rds", p2_dir))
se.mb.sll2.Murcia <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.Murcia.rds", p2_dir))
igraphs.Murcia.allNodes <- adj2igraph(getRefit(se.mb.sll2.Murcia),  vertex.attr=list(name=taxa_names(tagl.Murcia)), rmEmptyNodes = F)
igraphs.Murcia.noEmpty <- adj2igraph(getRefit(se.mb.sll2.Murcia),  vertex.attr=list(name=taxa_names(tagl.Murcia)), rmEmptyNodes = T)


# PV
samps.PV <- rownames(SLL2.meta[ SLL2.meta[, "Community"]=="País Vasco" &
                                  SLL2.meta[, "Chronic_disorder"]=="No", ])
tagl.PV  <- subset_samples(tagl_healthy, sample_names(tagl_healthy) %in% samps.PV)
se.mb.sll2.PV <- spiec.easi(tagl.PV, method="mb", lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.mb.sll2.PV, file = sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.PV.rds", p2_dir))
se.mb.sll2.PV <- readRDS(sprintf("%s/R_objects/SpiecEasi/se.mb.sll2.PV.rds", p2_dir))
igraphs.PV.allNodes <- adj2igraph(getRefit(se.mb.sll2.PV),  vertex.attr=list(name=taxa_names(tagl.PV)), rmEmptyNodes = F)
igraphs.PV.noEmpty <- adj2igraph(getRefit(se.mb.sll2.PV),  vertex.attr=list(name=taxa_names(tagl.PV)), rmEmptyNodes = T)






# ************************************************** #
# ************************************************** #

# use package intergraph to convert igraphs to networks, from which I can pull out all the edge connections ####
library(intergraph)

edge.mat.Child   <- as.matrix(asNetwork(igraphs.Child.allNodes))
edge.mat.Teen    <- as.matrix(asNetwork(igraphs.Teen.allNodes))
edge.mat.Adult   <- as.matrix(asNetwork(igraphs.Adult.allNodes))
edge.mat.Senior  <- as.matrix(asNetwork(igraphs.Senior.allNodes))
edge.mat.Younger <- as.matrix(asNetwork(igraphs.Adult.allNodes))
edge.mat.Older   <- as.matrix(asNetwork(igraphs.Senior.allNodes))

# for each genus, get vector of those other genera with which it has an edge
gen_edges.Child   <- list()
gen_edges.Teen    <- list()
gen_edges.Adult   <- list()
gen_edges.Senior  <- list()
gen_edges.Younger <- list()
gen_edges.Older   <- list()
for (gen in rownames(edge.mat.Child)) {
  
  gen_edges.Child[[ gen ]]   <- colnames(edge.mat.Child)[ edge.mat.Child[ gen, ] != 0 ]
  gen_edges.Teen[[ gen ]]    <- colnames(edge.mat.Teen)[ edge.mat.Teen[ gen, ] != 0 ]
  gen_edges.Adult[[ gen ]]   <- colnames(edge.mat.Adult)[ edge.mat.Adult[ gen, ] != 0 ]
  gen_edges.Senior[[ gen ]]  <- colnames(edge.mat.Senior)[ edge.mat.Senior[ gen, ] != 0 ]
  gen_edges.Younger[[ gen ]] <- colnames(edge.mat.Younger)[ edge.mat.Younger[ gen, ] != 0 ]
  gen_edges.Older[[ gen ]]   <- colnames(edge.mat.Older)[ edge.mat.Older[ gen, ] != 0 ]
  
}





# ******************** #
# first test with 1st subsamp 
#   - for each genus, get vector of connected genera that are same in both CF and mC
#   - vector of genera only connected in CF
#   - vector of genera only connected in mC

gen_connect.age_groups <- list()

gen_connect.age_groups[[ "All" ]] <- sapply(names(gen_edges.Child), function(gen) 
  gen_edges.Child[[ gen ]][ gen_edges.Child[[ gen ]] %in% c(gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ] )

gen_connect.age_groups[[ "Child" ]] <- sapply(names(gen_edges.Child), function(gen) 
  gen_edges.Child[[ gen ]][ ! gen_edges.Child[[ gen ]] %in% c(gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Teen" ]] <- sapply(names(gen_edges.Teen), function(gen) 
  gen_edges.Teen[[ gen ]][ ! gen_edges.Teen[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Adult[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Adult" ]] <- sapply(names(gen_edges.Adult), function(gen) 
  gen_edges.Adult[[ gen ]][ ! gen_edges.Adult[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Teen[[ gen ]], gen_edges.Senior[[ gen ]]) ])

gen_connect.age_groups[[ "Senior" ]] <- sapply(names(gen_edges.Senior), function(gen) 
  gen_edges.Senior[[ gen ]][ ! gen_edges.Senior[[ gen ]] %in% c(gen_edges.Child[[ gen ]], gen_edges.Teen[[ gen ]], gen_edges.Adult[[ gen ]]) ])

gen_connect.age_groups[[ "Younger" ]] <- sapply(names(gen_edges.Younger), function(gen) 
  gen_edges.Younger[[ gen ]][ ! gen_edges.Younger[[ gen ]] %in% gen_edges.Older[[ gen ]] ])

gen_connect.age_groups[[ "Older" ]] <- sapply(names(gen_edges.Older), function(gen) 
  gen_edges.Older[[ gen ]][ ! gen_edges.Older[[ gen ]] %in% gen_edges.Younger[[ gen ]] ])


# *************************************** #
gen_con_freqs.age_groups <- function(gc_list, gen) {
  
  cat(sprintf("**** %s connections in All ****\n", gen))
  print(sort(gc_list$All[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Child only ****\n", gen))
  print(sort(gc_list$Child[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Teen only ****\n", gen))
  print(sort(gc_list$Teen[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Adult only ****\n", gen))
  print(sort(gc_list$Adult[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Senior only ****\n", gen))
  print(sort(gc_list$Senior[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Younger only ****\n", gen))
  print(sort(gc_list$Younger[[ gen ]]))
  
  cat(sprintf("\n**** %s connections in Older only ****\n", gen))
  print(sort(gc_list$Older[[ gen ]]))
  
}
# *************************************** #


gen_con_freqs.age_groups(gen_connect.age_groups, "Porphyromonas")

# *************************************** #




Child.net <- get_network_objects("Child", se.mb.sll2.Child, tagl.Child, taxTables.both$Genus, 
                                 addLabels = T, Vprop = T, Eprop = T)


Teen.net <- get_network_objects("Teen", se.mb.sll2.Teen, tagl.Teen, taxTables.both$Genus, 
                                addLabels = T, Vprop = T, Eprop = T)


Adult.net <- get_network_objects("Adult", se.mb.sll2.Adult, tagl.Adult, taxTables.both$Genus, 
                                 addLabels = T, Vprop = T, Eprop = T)


Senior.net <- get_network_objects("Senior", se.mb.sll2.Senior, tagl.Senior, taxTables.both$Genus, 
                                  addLabels = T, Vprop = T, Eprop = T)


Younger.net <- get_network_objects("Younger", se.mb.sll2.Younger, tagl.Younger, taxTables.both$Genus, 
                                   addLabels = T, Vprop = T, Eprop = T)


Older.net <- get_network_objects("Older", se.mb.sll2.Older, tagl.Older, taxTables.both$Genus, 
                                 addLabels = T, Vprop = T, Eprop = T)



# ****************************************************************************************************************** ####


































# ****************************************************************************************************************** #

# ************************************************ #
# Co-occurrence networks a la SLL1 ####
# tutorial here: https://statnet.org/trac/raw-attachment/wiki/Resources/introToSNAinR_sunbelt_2012_tutorial.pdf
# ************************************************ #

*******
  In figure 1 here: https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1004226&type=printable
  their networks have solid edges when rho >= 0.5 is a solid line, rho >= 0.35 is a dotted line ... good way to simply display strength
*************


2
tl <- "Genus"

get_correlation_matrices <- function(samps, tl, glomTabList) {
  
  glomTab <- glomTabList[[ tl ]][ , samps ]
  
  res.matrix <- matrix(NA, nrow = nrow(glomTab), ncol = nrow(glomTab))
  colnames(res.matrix) <- rownames(glomTab)
  rownames(res.matrix) <- rownames(glomTab)
  ps.matrix <- matrix(NA, nrow = nrow(glomTab), ncol = nrow(glomTab))
  colnames(ps.matrix) <- rownames(glomTab)
  rownames(ps.matrix) <- rownames(glomTab)
  
  for (i in rownames(glomTab)) {
    for (j in rownames(glomTab)) {
      correl <- cor.test(glomTab[i,], glomTab[j,], na.rm=T)
      res.matrix[i,j] <- correl$estimate
      ps.matrix[i,j] <- correl$p.value
    }
  }
  
  ps.matrix.adj <- apply(ps.matrix, 2, p.adjust, method='bonferroni', n=nrow(glomTab))
  
  return(list("res"=res.matrix, "ps"=ps.matrix, "ps.adj"=ps.matrix.adj))
}


glomTab.CF <- gloms_clr[[ tl ]][ , CF.samps]
glomTab.mC <- gloms_clr[[ tl ]][ , mC.samps]


res.matrix <- matrix(NA, nrow = nrow(glomTab.CF), ncol = ncol(glomTab.CF))
colnames(res.matrix) <- rownames(glomTab.CF)
rownames(res.matrix) <- rownames(glomTab.CF)
ps.matrix <- matrix(NA, nrow = nrow(glomTab.CF), ncol = ncol(glomTab.CF))
colnames(ps.matrix) <- rownames(glomTab.CF)
rownames(ps.matrix) <- rownames(glomTab.CF)

for (i in rownames(glomTab.CF)) {
  for (j in rownames(glomTab.CF)) {
    correl <- cor.test(glomTab.CF[i,], glomTab.CF[j,], na.rm=T)
    res.matrix[i,j] <- correl$estimate
    ps.matrix[i,j] <- correl$p.value
  }
}

ps.matrix.adj <- apply(ps.matrix, 2, p.adjust, method='bonferroni', n=nrow(glomTab.CF))






# ****************************************************************************************************************** #














# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

# Adonis and ANOVA for Age_groups with matched Teen and Adult subsamplings ####

# ************************************************************************ #
# ************************************************************************ #


# healthySamps <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No", ])
# meta.healthy <- SLL2.meta[ healthySamps, ]
# phy.healthy  <- prune_samples( healthySamps, SLL2)


healthy.nonBottle <- healthySamps[ ! is.na(SLL2.meta[healthySamps, "Water_type_home"]) &
                                     SLL2.meta[healthySamps, "Water_type_home"] != "Embotellada" ]
meta.nonBotHeal <- SLL2.meta[ healthy.nonBottle, ]
phy.nonBotHeal  <- prune_samples( healthy.nonBottle, SLL2)

# 
# ageGroupSubsTests <- run_full_subsampling_calcs.ageGroups("healthy", phy.healthy, meta.healthy, 100)
# ageGroupSubsTests.TAS <- run_full_subsampling_calcs.ageGroups("healthy", phy.healthy, meta.healthy, 100, TAS=T)
# ageGroupSubsTests.TAS_ph_bmi <- run_full_subsampling_calcs.ageGroups("healthy", phy.healthy, meta.healthy, 100, TAS=T, pH_BMI=T, 
#                                                                      chosenControls = lapply(ageGroupSubsTests.TAS, function(x) x$samples))
# ageGroupSubsTests.nonBottle <- run_full_subsampling_calcs.ageGroups("healthy.nonBottle", phy.nonBotHeal, meta.nonBotHeal, 100)
# 
# saveRDS(ageGroupSubsTests, file = sprintf("%s/R_objects/ageGroupSubsTests.rds", p2_dir))
# saveRDS(ageGroupSubsTests.TAS, file = sprintf("%s/R_objects/ageGroupSubsTests.TAS.rds", p2_dir))
# saveRDS(ageGroupSubsTests.TAS_ph_bmi, file = sprintf("%s/R_objects/ageGroupSubsTests.TAS_ph_bmi.rds", p2_dir))
# saveRDS(ageGroupSubsTests.nonBottle, file = sprintf("%s/R_objects/ageGroupSubsTests.nonBottle.rds", p2_dir))
# # ************************************************************************ #

ageGroupSubsTests <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.rds", p2_dir))
ageGroupSubsTests.TAS <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.TAS.rds", p2_dir))
ageGroupSubsTests.TAS_ph_bmi <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.TAS_ph_bmi.rds", p2_dir))
ageGroupSubsTests.nonBottle <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.nonBottle.rds", p2_dir))
# ************************************************************************ #







# ************************************************************************ #

gen.pMeans.ag   <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "Genus", "mean")
phy.pMeans.ag   <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "Phylum", "mean")
conts.pMeans.ag <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "contVar", "mean")

gen.num_sig.ag   <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "Genus", "num_sig")
phy.num_sig.ag   <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "Phylum", "num_sig")
conts.num_sig.ag <- get_subs_anova_per_taxon(ageGroupSubsTests, "Age_groups", "contVar", "num_sig")



# ************************************************************************ #

gen.pMeans.ageTAS   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "Genus", "mean")
phy.pMeans.ageTAS   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "Phylum", "mean")
conts.pMeans.ageTAS <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "contVar", "mean")

gen.num_sig.ageTAS   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "Genus", "num_sig")
phy.num_sig.ageTAS   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "Phylum", "num_sig")
conts.num_sig.ageTAS <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS, "Age_groups", "contVar", "num_sig")
# ************************************************************************ #


gen.pMeans.TAS_ph_bmi   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "Genus", "mean")
phy.pMeans.TAS_ph_bmi   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "Phylum", "mean")
conts.pMeans.TAS_ph_bmi <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "contVar", "mean")

gen.num_sig.TAS_ph_bmi   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "Genus", "num_sig")
phy.num_sig.TAS_ph_bmi   <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "Phylum", "num_sig")
conts.num_sig.TAS_ph_bmi <- get_subs_anova_per_taxon(ageGroupSubsTests.TAS_ph_bmi, "Age_groups", "contVar", "num_sig")
# ************************************************************************ #








# ************************************************************************ #
# mean and num sig for each adonis test ####



ado.Pmeans.ag <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                        function(x) mean(p.adjust(unlist(lapply(ageGroupSubsTests, function(y) 
                          y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.Psds.ag <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                      function(x) sd(p.adjust(unlist(lapply(ageGroupSubsTests, function(y) 
                        y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.num_sig.ag <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                         function(x) sum(p.adjust(unlist(lapply(ageGroupSubsTests, function(y) 
                           y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.ag <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                         function(x) mean(unlist(lapply(ageGroupSubsTests, function(y) 
                           y$Adonis[[ x ]]["Age_groups","R2"]))))

ado.Fmeans.ag <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                        function(x) mean(unlist(lapply(ageGroupSubsTests, function(y) 
                          y$Adonis[[ x ]]["Age_groups","F.Model"]))))


aG.samps <- unique(unlist(lapply(ageGroupSubsTests, function(x) x$samples)))
phy.aG <- prune_samples(aG.samps, SLL2)
mTab.aG <- SLL2.meta[ aG.samps, ]

# ordObj.aG <- subsampling_ordination_objects(aG.samps, phy.aG, distsOnly=T)
# saveRDS(ordObj.aG, file = sprintf("%s/R_objects/ordObj.aG.rds", p2_dir))
ordObj.aG <- readRDS(sprintf("%s/R_objects/ordObj.aG.rds", p2_dir))

ordObj.aG.pcoa <- list()
ordObj.aG.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.aG$Aitchison)
ordObj.aG.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.aG$Weighted_Unifrac)
ordObj.aG.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.aG$Unweighted_Unifrac)
ordObj.aG.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.aG$Bray)
ordObj.aG.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.aG$Jaccard)

plot_adonis_w_covars("Aitchison", phy.aG, mTab.aG, gloms_clr, 
                     shapeBy = "Age_groups", colorBy = "Age_groups", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.aG.pcoa)



# ************************************************************************************** #



ado.Pmeans.ageTAS <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(p.adjust(unlist(lapply(ageGroupSubsTests.TAS, function(y) 
                              y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.Psds.ageTAS <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                          function(x) sd(p.adjust(unlist(lapply(ageGroupSubsTests.TAS, function(y) 
                            y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.num_sig.ageTAS <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) sum(p.adjust(unlist(lapply(ageGroupSubsTests.TAS, function(y) 
                               y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.ageTAS <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) mean(unlist(lapply(ageGroupSubsTests.TAS, function(y) 
                               y$Adonis[[ x ]]["Age_groups","R2"]))))

ado.Fmeans.ageTAS <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(unlist(lapply(ageGroupSubsTests.TAS, function(y) 
                              y$Adonis[[ x ]]["Age_groups","F.Model"]))))



agTAS.samps <- unique(unlist(lapply(ageGroupSubsTests.TAS, function(x) x$samples)))
phy.agTAS <- prune_samples(agTAS.samps, SLL2)
mTab.agTAS <- SLL2.meta[ agTAS.samps, ]

# ordObj.agTAS <- subsampling_ordination_objects(agTAS.samps, phy.agTAS, distsOnly=T)
# saveRDS(ordObj.agTAS, file = sprintf("%s/R_objects/ordObj.agTAS.rds", p2_dir))
ordObj.agTAS <- readRDS(sprintf("%s/R_objects/ordObj.agTAS.rds", p2_dir))

ordObj.agTAS.pcoa <- list()
ordObj.agTAS.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.agTAS$Aitchison)
ordObj.agTAS.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.agTAS$Weighted_Unifrac)
ordObj.agTAS.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.agTAS$Unweighted_Unifrac)
ordObj.agTAS.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.agTAS$Bray)
ordObj.agTAS.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.agTAS$Jaccard)

plot_adonis_w_covars("Aitchison", phy.agTAS, mTab.agTAS, gloms_clr, 
                     shapeBy = "Age_groups", colorBy = "Age_groups", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.agTAS.pcoa)



# ************************************************************************************** #



ado.Pmeans.TAS_ph_bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(p.adjust(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(y) 
                              y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.Psds.TAS_ph_bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                          function(x) sd(p.adjust(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(y) 
                            y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr")))

ado.num_sig.TAS_ph_bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) sum(p.adjust(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(y) 
                               y$Adonis[[ x ]]["Age_groups","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.TAS_ph_bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) mean(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(y) 
                               y$Adonis[[ x ]]["Age_groups","R2"]))))

ado.Fmeans.TAS_ph_bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(y) 
                              y$Adonis[[ x ]]["Age_groups","F.Model"]))))



agTAS_ph_bmi.samps <- unique(unlist(lapply(ageGroupSubsTests.TAS_ph_bmi, function(x) x$samples)))
phy.agTAS_ph_bmi <- prune_samples(agTAS_ph_bmi.samps, SLL2)
mTab.agTAS_ph_bmi <- SLL2.meta[ agTAS_ph_bmi.samps, ]

# ordObj.agTAS_ph_bmi <- subsampling_ordination_objects(agTAS_ph_bmi.samps, phy.agTAS_ph_bmi, distsOnly=T)
# saveRDS(ordObj.agTAS_ph_bmi, file = sprintf("%s/R_objects/ordObj.agTAS_ph_bmi.rds", p2_dir))
ordObj.agTAS_ph_bmi <- readRDS(sprintf("%s/R_objects/ordObj.agTAS_ph_bmi.rds", p2_dir))

ordObj.agTAS_ph_bmi.pcoa <- list()
ordObj.agTAS_ph_bmi.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.agTAS_ph_bmi$Aitchison)
ordObj.agTAS_ph_bmi.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.agTAS_ph_bmi$Weighted_Unifrac)
ordObj.agTAS_ph_bmi.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.agTAS_ph_bmi$Unweighted_Unifrac)
ordObj.agTAS_ph_bmi.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.agTAS_ph_bmi$Bray)
ordObj.agTAS_ph_bmi.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.agTAS_ph_bmi$Jaccard)

plot_adonis_w_covars("Aitchison", phy.agTAS_ph_bmi, mTab.agTAS_ph_bmi, gloms_clr, 
                     shapeBy = "Age_groups", colorBy = "Age_groups", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.agTAS_ph_bmi.pcoa)








# ************************************************************************************** #
# Subsampling plots for age groups ####
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ag <- gen.num_sig.ag[ gen.num_sig.ag >= 75]
sort(gen.pMeans.ag[ names(freqSig.ag) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ag)))

gensToCheck <- unique(names(freqSig.ag))
# as.data.frame(cbind(gen.num_sig.ag[ gensToCheck ], formatC(gen.pMeans.ag[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps <- as.data.frame(cbind(gen.num_sig.ag[ gensToCheck ], gen.pMeans.ag[ gensToCheck ]))
colnames(freq.mean.ps) <- c("num_sig","meanP")
freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig, -freq.mean.ps$meanP)), ]
freq.mean.ps


subsampling.plots.box("Age_groups","Genus","Age_groups", names(rev(sort(freqSig.ag))), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

phy.sig.order <- phy.num_sig.ag[ rev(order(phy.num_sig.ag, -phy.pMeans.ag)) ]
subsampling.plots.box("Age_groups","Phylum","Age_groups", names(phy.sig.order), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

subsampling.plots.box("Age_groups","contVar","Age_groups", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests, dstStruc = "Age_groups",
                      plotType = "box", plot_tukey = F, xAngle = 15)

# ************************************************************************************** #





freqSig.ageTAS <- gen.num_sig.ageTAS[ gen.num_sig.ageTAS >= 75]
sort(gen.pMeans.ageTAS[ names(freqSig.ageTAS) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ageTAS)))

gensToCheck.ageTAS <- unique(names(freqSig.ageTAS))
# as.data.frame(cbind(gen.num_sig.ag[ gensToCheck ], formatC(gen.pMeans.ag[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.ageTAS <- as.data.frame(cbind(gen.num_sig.ageTAS[ gensToCheck.ageTAS ], gen.pMeans.ageTAS[ gensToCheck.ageTAS ]))
colnames(freq.mean.ps.ageTAS) <- c("num_sig","meanP")
freq.mean.ps.ageTAS <- freq.mean.ps.ageTAS[ rev(order(freq.mean.ps.ageTAS$num_sig, -freq.mean.ps.ageTAS$meanP)), ]
freq.mean.ps.ageTAS


subsampling.plots.box("Age_groups","Genus","Age_groups", rownames(freq.mean.ps.ageTAS), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

phy.sig.order <- phy.num_sig.ageTAS[ rev(order(phy.num_sig.ageTAS, -phy.pMeans.ageTAS)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.ageTAS[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age_groups","Phylum","Age_groups", names(phy.sig.order), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

subsampling.plots.box("Age_groups","contVar","Age_groups", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS, dstStruc = "Age_groups",
                      plotType = "box", plot_tukey = F, xAngle = 15)






subsampling.plots.box("Age_groups","Genus","Age", rownames(freq.mean.ps.ageTAS), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 15)
# ************************************************************************************** #






freqSig.TAS_ph_bmi <- gen.num_sig.TAS_ph_bmi[ gen.num_sig.TAS_ph_bmi >= 75]
sort(gen.pMeans.TAS_ph_bmi[ names(freqSig.TAS_ph_bmi) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.TAS_ph_bmi)))

gensToCheck.TAS_ph_bmi <- unique(names(freqSig.TAS_ph_bmi))
# as.data.frame(cbind(gen.num_sig.ag[ gensToCheck ], formatC(gen.pMeans.ag[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.TAS_ph_bmi <- as.data.frame(cbind(gen.num_sig.TAS_ph_bmi[ gensToCheck.TAS_ph_bmi ], gen.pMeans.TAS_ph_bmi[ gensToCheck.TAS_ph_bmi ]))
colnames(freq.mean.ps.TAS_ph_bmi) <- c("num_sig","meanP")
freq.mean.ps.TAS_ph_bmi <- freq.mean.ps.TAS_ph_bmi[ rev(order(freq.mean.ps.TAS_ph_bmi$num_sig, -freq.mean.ps.TAS_ph_bmi$meanP)), ]
freq.mean.ps.TAS_ph_bmi


subsampling.plots.box("Age_groups","Genus","Age_groups", rownames(freq.mean.ps.TAS_ph_bmi), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS_ph_bmi, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

phy.sig.order <- phy.num_sig.TAS_ph_bmi[ rev(order(phy.num_sig.TAS_ph_bmi, -phy.pMeans.TAS_ph_bmi)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.TAS_ph_bmi[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age_groups","Phylum","Age_groups", names(phy.sig.order), groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS_ph_bmi, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 15)

subsampling.plots.box("Age_groups","contVar","Age_groups", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      groupQs, only_cont,
                      gloms_clr, SLL2.meta, SLL2, ageGroupSubsTests.TAS_ph_bmi, dstStruc = "Age_groups",
                      plotType = "box", plot_tukey = F, xAngle = 15)

# ************************************************************************************** #














# ************************************************************************************** #

# ageGroup.pairwise.subsTests <- run_pairwise_subsampling_calcs.ageGroups("healthy", phy.healthy, meta.healthy, 100, "Adonis",
#                                                                         chosenControls = all_subSamps$Age_groups)
# TAS.subSamps <- lapply(ageGroupSubsTests.TAS, function(x) x$samples)
# ageGroup.pairwise.subsTests.TAS <- run_pairwise_subsampling_calcs.ageGroups("healthy", phy.healthy, meta.healthy, 100, "Adonis",
#                                                                             chosenControls = TAS.subSamps, TAS = T)
# 
# saveRDS(ageGroup.pairwise.subsTests, file = sprintf("%s/R_objects/ageGroup.pairwise.subsTests.rds", p2_dir))
# saveRDS(ageGroup.pairwise.subsTests.TAS, file = sprintf("%s/R_objects/ageGroup.pairwise.subsTests.TAS.rds", p2_dir))

ageGroup.pairwise.subsTests <- readRDS(sprintf("%s/R_objects/ageGroup.pairwise.subsTests.rds", p2_dir))
ageGroup.pairwise.subsTests.TAS <- readRDS(sprintf("%s/R_objects/ageGroup.pairwise.subsTests.TAS.rds", p2_dir))




agp.num_sig <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests, function(a) 
    sum(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","Pr(>F)"])), method="fdr") < 0.05 )))

agp.Pmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests, function(a) 
    mean(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","Pr(>F)"])), method="fdr"))))



agp.R2means <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","R2"])))))

agp.Fmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","F.Model"])))))



lapply(list("num_sig"=agp.num_sig, "Pmeans"=agp.Pmeans, "R2means"=agp.R2means, "Fmeans"=agp.Fmeans), function(x) unlist(x[,"Aitchison"]))

# ************************************************************************************** #


agp.TAS.num_sig <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests.TAS, function(a) 
    sum(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","Pr(>F)"])), method="fdr") < 0.05 )))

agp.TAS.Pmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests.TAS, function(a) 
    mean(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","Pr(>F)"])), method="fdr"))))



agp.TAS.R2means <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests.TAS, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","R2"])))))

agp.TAS.Fmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageGroup.pairwise.subsTests.TAS, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_groups","F.Model"])))))



as.data.frame(lapply(list("num_sig"=agp.TAS.num_sig, "Pmeans"=agp.TAS.Pmeans, "R2means"=agp.TAS.R2means, "Fmeans"=agp.TAS.Fmeans), 
                     function(x) unlist(x[,"Aitchison"])))

as.data.frame(lapply(list("num_sig"=ado.num_sig.ageTAS, "Pmeans"=ado.Pmeans.ageTAS, "R2means"=ado.R2means.ageTAS, 
                          "Fmeans"=ado.Fmeans.ageTAS), 
                     function(x) unlist(x["Aitchison"])))

# ************************************************************************************** #










# ****************************************************************************************************************** ####

# Subsample tests for Age_bins ####
meta.healthy$Age_bins <- sapply(rownames(meta.healthy), function(x) 
  ifelse( is.na(meta.healthy[x, "Age"]), NA,
          ifelse(meta.healthy[x, "Age"] >= 0 & meta.healthy[x, "Age"] < 13, "0_13",
                 ifelse(meta.healthy[x, "Age"] >= 13 & meta.healthy[x, "Age"] < 20, "13_20",
                        ifelse(meta.healthy[x, "Age"] >= 20 & meta.healthy[x, "Age"] < 30, "20_30",
                               ifelse(meta.healthy[x, "Age"] >= 30 & meta.healthy[x, "Age"] < 40, "30_40",
                                      ifelse(meta.healthy[x, "Age"] >= 40 & meta.healthy[x, "Age"] < 50, "40_50",
                                             ifelse(meta.healthy[x, "Age"] >= 50 & meta.healthy[x, "Age"] < 60, "50_60", "60+"))))))))

mTab.abTAS <- meta.healthy[ ! is.na(meta.healthy[,"Age_bins"]) & meta.healthy[,"Age_bins"] != "0_13", ]

table(mTab.abTAS$Age_bins, useNA = "always")
table(mTab.abTAS$Age_bins, mTab.abTAS$Gender)
table(mTab.abTAS$Age_bins, mTab.abTAS$Community)

sapply(colnames(table(mTab.abTAS$Age_bins, mTab.abTAS$Gender)), function(x)
  table(mTab.abTAS$Age_bins, mTab.abTAS$Gender)[,x] / rowSums(table(mTab.abTAS$Age_bins, mTab.abTAS$Gender)) * 100)

sapply(colnames(table(mTab.abTAS$Age_bins, mTab.abTAS$Community)), function(x)
  table(mTab.abTAS$Age_bins, mTab.abTAS$Community)[,x] / rowSums(table(mTab.abTAS$Age_bins, mTab.abTAS$Community)) * 100)


phy.abTAS <- phyloseq(otu_table(prune_samples(rownames(mTab.abTAS), SLL2)),
                      sample_data(mTab.abTAS),
                      phy_tree(prune_samples(rownames(mTab.abTAS), SLL2)),
                      tax_table(prune_samples(rownames(mTab.abTAS), SLL2)))

# ******************************************************************* #
# ageBinsSubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.abTAS, mTab.abTAS, 100)
# 
# saveRDS(ageBinsSubsTests, file = sprintf("%s/R_objects/ageBinsSubsTests.rds", p2_dir))
# # ************************************************************************ #

ageBinsSubsTests <- readRDS(sprintf("%s/R_objects/ageBinsSubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.bins   <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "Genus", "mean")
phy.pMeans.bins   <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "Phylum", "mean")
conts.pMeans.bins <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "contVar", "mean")

gen.num_sig.bins   <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "Genus", "num_sig")
phy.num_sig.bins   <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "Phylum", "num_sig")
conts.num_sig.bins <- get_subs_anova_per_taxon(ageBinsSubsTests, "Age_bins", "contVar", "num_sig")




# ************************************************************************************** #



ado.Pmeans.bins <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(p.adjust(unlist(lapply(ageBinsSubsTests, function(y) 
                              y$Adonis[[ x ]]["Age_bins","Pr(>F)"])), method = "fdr")))

ado.Psds.bins <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                          function(x) sd(p.adjust(unlist(lapply(ageBinsSubsTests, function(y) 
                            y$Adonis[[ x ]]["Age_bins","Pr(>F)"])), method = "fdr")))

ado.num_sig.bins <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) sum(p.adjust(unlist(lapply(ageBinsSubsTests, function(y) 
                               y$Adonis[[ x ]]["Age_bins","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.bins <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                             function(x) mean(unlist(lapply(ageBinsSubsTests, function(y) 
                               y$Adonis[[ x ]]["Age_bins","R2"]))))

ado.Fmeans.bins <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) mean(unlist(lapply(ageBinsSubsTests, function(y) 
                              y$Adonis[[ x ]]["Age_bins","F.Model"]))))



bins.samps <- unique(unlist(lapply(ageBinsSubsTests, function(x) x$samples)))

# ordObj.bins <- subsampling_ordination_objects(bins.samps, phy.abTAS, distsOnly=T)
# saveRDS(ordObj.bins, file = sprintf("%s/R_objects/ordObj.bins.rds", p2_dir))
ordObj.bins <- readRDS(sprintf("%s/R_objects/ordObj.bins.rds", p2_dir))

ordObj.bins.pcoa <- list()
ordObj.bins.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.bins$Aitchison)
ordObj.bins.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.bins$Weighted_Unifrac)
ordObj.bins.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.bins$Unweighted_Unifrac)
ordObj.bins.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.bins$Bray)
ordObj.bins.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.bins$Jaccard)

plot_adonis_w_covars("Aitchison", phy.bins, mTab.abTAS, gloms_clr, 
                     shapeBy = "Age_bins", colorBy = "Age_bins", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.bins.pcoa)



# ************************************************************************************** #





# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ab <- gen.num_sig.bins[ gen.num_sig.bins >= 75]
sort(gen.pMeans.bins[ names(freqSig.ab) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ab)))

gensToCheck <- unique(names(freqSig.ab))
# as.data.frame(cbind(gen.num_sig.bins[ gensToCheck ], formatC(gen.pMeans.bins[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps <- as.data.frame(cbind(gen.num_sig.bins[ gensToCheck ], gen.pMeans.bins[ gensToCheck ]))
colnames(freq.mean.ps) <- c("num_sig","meanP")
freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig, -freq.mean.ps$meanP)), ]
freq.mean.ps


subsampling.plots.box("Age_bins","Genus","Age_bins", names(rev(sort(freqSig.ab))), c(groupQs,"Age_bins"), only_cont,
                      gloms_clr, mTab.abTAS, phy.abTAS, ageBinsSubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 30)

phy.sig.order <- phy.num_sig.bins[ rev(order(phy.num_sig.bins, -phy.pMeans.bins)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.bins[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age_bins","Phylum","Age_bins", names(phy.sig.order), c(groupQs,"Age_bins"), only_cont,
                      gloms_clr, mTab.abTAS, phy.abTAS, ageBinsSubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age_bins","contVar","Age_bins", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs,"Age_bins"), only_cont,
                      gloms_clr, mTab.abTAS, phy.abTAS, ageBinsSubsTests, dstStruc = "Age_groups",
                      plotType = "box", plot_tukey = F, xAngle = 30)

# ************************************************************************************** #


subsampling.plots.box("Age_bins","Genus","Age_bins", 
                      c("Anaeroglobus", "Eikenella", "Fretibacterium", "Comamonas", "Olsenella", "Phocaeicola",
                        "Alloprevotella", "Streptobacillus", "Haemophilus", "Prevotella", "Granulicatella", "Bergeyella"),
                      c(groupQs,"Age_bins"), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, ageBinsSubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, xAngle = 30)


# ************************************************************************************** #





# ageBins.subSamps <- lapply(ageBinsSubsTests, function(x) x$samples)
# ageBins.pairwise.subsTests <- run_pairwise_subsampling_calcs.ageBins("healthy", phy.ageBins, mTab.ageBins, 100, "Adonis",
#                                                                      chosenControls = ageBins.subSamps)
# saveRDS(ageBins.pairwise.subsTests, file = sprintf("%s/R_objects/ageBins.pairwise.subsTests.rds", p2_dir))

ageBins.pairwise.subsTests <- readRDS(sprintf("%s/R_objects/ageBins.pairwise.subsTests.rds", p2_dir))
# ************************************************************************************** #


abp.num_sig <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageBins.pairwise.subsTests, function(a) 
    sum(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_bins","Pr(>F)"])), method="fdr") < 0.05 )))

abp.Pmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageBins.pairwise.subsTests, function(a) 
    mean(p.adjust(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_bins","Pr(>F)"])), method="fdr"))))



abp.R2means <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageBins.pairwise.subsTests, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_bins","R2"])))))

abp.Fmeans <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"), function(di) 
  lapply(ageBins.pairwise.subsTests, function(a) 
    mean(unlist(lapply(a, function(i) 
      i$Adonis[[ di ]]["Age_bins","F.Model"])))))



as.data.frame(lapply(list("num_sig"=abp.num_sig, "Pmeans"=abp.Pmeans, "R2means"=abp.R2means, "Fmeans"=abp.Fmeans), 
                     function(x) unlist(x[,"Aitchison"])))

as.data.frame(lapply(list("num_sig"=ado.num_sig.bins, "Pmeans"=ado.Pmeans.bins, "R2means"=ado.R2means.bins, 
                          "Fmeans"=ado.Fmeans.bins), 
                     function(x) unlist(x["Aitchison"])))

# ************************************************************************************** #




# ****************************************************************************************************************** #










# ****************************************************************************************************************** ####

# Subsample tests for Age (continuous) using Age_bins subsamps ####

# ******************************************************************* #
# age_cont.SubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100, 
#                                                          trait = "Age", chosenControls = all_subSamps$Age_bins)
# 
# saveRDS(age_cont.SubsTests, file = sprintf("%s/R_objects/age_cont.SubsTests.rds", p2_dir))
# # ************************************************************************ #

age_cont.SubsTests <- readRDS(sprintf("%s/R_objects/age_cont.SubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.age_cont   <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "Genus", "mean")
phy.pMeans.age_cont   <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "Phylum", "mean")
conts.pMeans.age_cont <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "contVar", "mean")

gen.num_sig.age_cont   <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "Genus", "num_sig")
phy.num_sig.age_cont   <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "Phylum", "num_sig")
conts.num_sig.age_cont <- get_subs_anova_per_taxon(age_cont.SubsTests, "Age", "contVar", "num_sig")




# ************************************************************************************** #



ado.Pmeans.age_cont <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                          function(x) mean(p.adjust(unlist(lapply(age_cont.SubsTests, function(y) 
                            y$Adonis[[ x ]]["Age","Pr(>F)"])), method = "fdr")))

ado.Psds.age_cont <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                        function(x) sd(p.adjust(unlist(lapply(age_cont.SubsTests, function(y) 
                          y$Adonis[[ x ]]["Age","Pr(>F)"])), method = "fdr")))

ado.num_sig.age_cont <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                           function(x) sum(p.adjust(unlist(lapply(age_cont.SubsTests, function(y) 
                             y$Adonis[[ x ]]["Age","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.age_cont <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                           function(x) mean(unlist(lapply(age_cont.SubsTests, function(y) 
                             y$Adonis[[ x ]]["Age","R2"]))))

ado.Fmeans.age_cont <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                          function(x) mean(unlist(lapply(age_cont.SubsTests, function(y) 
                            y$Adonis[[ x ]]["Age","F.Model"]))))



age_cont.samps <- unique(unlist(lapply(age_cont.SubsTests, function(x) x$samples)))
phy.age_cont <- prune_samples( age_cont.samps, phy.healthy)

# ordObj.age_cont <- subsampling_ordination_objects(age_cont.samps, phy.age_cont, print_current_dists=T, dists_and_pcoas_only=T)
# saveRDS(ordObj.age_cont, file = sprintf("%s/R_objects/ordObj.age_cont.rds", p2_dir))
ordObj.age_cont <- readRDS(sprintf("%s/R_objects/ordObj.age_cont.rds", p2_dir))

# ordObj.age_cont.pcoa <- list()
# ordObj.age_cont.pcoa[[ "Aitchison" ]]          <- ape::pcoa(ordObj.age_cont$Aitchison)
# ordObj.age_cont.pcoa[[ "Weighted_Unifrac" ]]   <- ape::pcoa(ordObj.age_cont$Weighted_Unifrac)
# ordObj.age_cont.pcoa[[ "Unweighted_Unifrac" ]] <- ape::pcoa(ordObj.age_cont$Unweighted_Unifrac)
# ordObj.age_cont.pcoa[[ "Bray" ]]               <- ape::pcoa(ordObj.age_cont$Bray)
# ordObj.age_cont.pcoa[[ "Jaccard" ]]            <- ape::pcoa(ordObj.age_cont$Jaccard)

plot_adonis_w_covars("Aitchison", phy.age_cont, meta.healthy[ age_cont.samps, ], gloms_clr, 
                     shapeBy = "Project", colorBy = "Age", 
                     NA,# adonis.pvals, adonis.covs.signifs,
                     pcoaList = ordObj.age_cont$pcoas)



# ************************************************************************************** #

agecont.sigGen <- c("Actinomyces","Alloprevotella","Anaeroglobus","Streptobacillus","Eikenella","Fretibacterium",
                    "Haemophilus","Prevotella","Granulicatella","Comamonas","Olsenella","Bergeyella","Phocaeicola")
as.data.frame(taxa_print("Genus", agecont.sigGen, tl.limit = "Genus")[ with(as.data.frame(taxa_print("Genus", agecont.sigGen, tl.limit = "Genus")),
                                                                            order(Phylum, Class, Order, Family)) ], row.names=F)



# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ab <- gen.num_sig.age_cont[ gen.num_sig.age_cont >= 75]
sort(gen.pMeans.age_cont[ names(freqSig.ab) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ab)))

gensToCheck <- unique(names(freqSig.ab))
# as.data.frame(cbind(gen.num_sig.age_cont[ gensToCheck ], formatC(gen.pMeans.age_cont[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps <- as.data.frame(cbind(gen.num_sig.age_cont[ gensToCheck ], gen.pMeans.age_cont[ gensToCheck ]))
colnames(freq.mean.ps) <- c("num_sig","meanP")
freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig, -freq.mean.ps$meanP)), ]
freq.mean.ps


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.ab))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1)

phy.sig.order <- phy.num_sig.age_cont[ rev(order(phy.num_sig.age_cont, -phy.pMeans.age_cont)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_cont[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1)


# ************************************************************************************** #
# use Loess for geom_smooth function:
subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      useLoess = T)

phy.sig.order <- phy.num_sig.age_cont[ rev(order(phy.num_sig.age_cont, -phy.pMeans.age_cont)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_cont[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      useLoess = T)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      useLoess = T)

# ************************************************************************************** #

divs.age.cont <- subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      useLoess = T, facetRows = 1)

divs.age.bins <- subsampling.plots.box("Age_bins","contVar","Age_bins", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                      c(groupQs,"Age_bins"), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, ageBinsSubsTests, dstStruc = "Age_groups",
                      plotType = "box", plot_tukey = F, xAngle = 30, pts=20, facetRows = 1)


ggarrange(divs.age.cont, divs.age.bins, ncol = 1, nrow = 2, labels = c("(a)", "(b)"))

# *********** #
divs.age.cont <- subsampling.plots.box("Age","contVar","Age", 
                                       c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                                       c(groupQs), only_cont,
                                       gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups",
                                       plotType = "scatter", plot_tukey = F, xAngle = 30,
                                       singleSubSamp = 1, useLoess = T, facetRows = 1)

divs.age.bins <- subsampling.plots.box("Age_bins","contVar","Age_bins", 
                                       c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                                       c(groupQs,"Age_bins"), only_cont,
                                       gloms_clr, meta.healthy, phy.healthy, ageBinsSubsTests, dstStruc = "Age_groups",
                                       plotType = "box", plot_tukey = F, xAngle = 30, 
                                       singleSubSamp = 1, pts=20, facetRows = 1)


ggarrange(divs.age.cont, divs.age.bins, ncol = 1, nrow = 2, labels = c("(a)", "(b)"))

# *********** #
divs.age.cont <- subsampling.plots.box("Age","contVar","Age", 
                                       c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                                       c(groupQs), only_cont,
                                       gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                                       plotType = "scatter", plot_tukey = F, xAngle = 30,
                                       singleSubSamp = 1, useLoess = T, facetRows = 1)

divs.age.bins <- subsampling.plots.box("Age_bins","contVar","Age_bins", 
                                       c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
                                       c(groupQs,"Age_bins"), only_cont,
                                       gloms_clr, meta.healthy, phy.healthy, ageBinsSubsTests, dstStruc = "Age_groups",
                                       plotType = "box", plot_tukey = F, xAngle = 30, 
                                       singleSubSamp = 1, pts=20, facetRows = 1)


ggarrange(divs.age.cont, divs.age.bins, ncol = 1, nrow = 2, labels = c("(a)", "(b)"))

# ************************************************************************************** #

# plot the differing genera with boxes of 13-60 and >60
meta.healthy$Age.Senior_other <- sapply(meta.healthy$Age, function(x) ifelse(x > 60, ">60", "13-60"))
# meta.healthy$Age.Senior_other <- factor(meta.healthy$Age.Senior_other, levels = c("13-60",">60"))

age_diff_genera <- c("Anaeroglobus","Eikenella","Fretibacterium","Comamonas","Olsenella","Phocaeicola",
                     "Alloprevotella","Streptobacillus","Haemophilus","Prevotella","Granulicatella","Bergeyella")
subsampling.plots.box("Age.Senior_other","Genus","Age.Senior_other", age_diff_genera, c(groupQs,"Age.Senior_other"), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "box", plot_tukey = F, include_legend_below=T, sts=15, pts = 20,
                      facetScales = "fixed", flipCoords = T, facetRows = 6, choose_yl="Abundance (centered log ratio)",#xAngle = 30, 
                      singleSubSamp = 1, ignore_pvals = T, adjustAlphas = T)

# ************************************************************************************** #





# ****************************************************************************************************************** #

Anova(lm(Div.Shannon ~ Age + Gender + Population, data = meta.healthy[ ! is.na(meta.healthy$Age),]))
Anova(lm(Div.Shannon ~ splines::bs(Age,  degree = 2) + Gender + Population, data = meta.healthy[ age_cont.SubsTests$`1`$samples,]))

get_lm( c("Age", "Gender","Population"), 
        "contVar", meta.healthy[age_cont.SubsTests$`1`$samples,], NULL, 
        c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts",
          "pH","BMI"), 
        noRemove = T, rerun.nonSig = F,
        polyRegression = "Age", polyKnots = "NULL", polyDegree = 2)

# Subsample tests for Age (continuous) using Age_bins subsamps - polynomial regression ####

# # ******************************************************************* #
# age_quadratic.SubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                               trait = "Age", chosenControls = all_subSamps$Age_bins,
#                                                               polyRegression = "Age", polyKnots = "NULL", polyDegree = 2,
#                                                               anova_only = T)
# 
# saveRDS(age_quadratic.SubsTests, file = sprintf("%s/R_objects/age_quadratic.SubsTests.rds", p2_dir))
# # ************************************************************************ #

age_quadratic.SubsTests <- readRDS(sprintf("%s/R_objects/age_quadratic.SubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.age_quadratic   <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "Genus", "mean")
phy.pMeans.age_quadratic   <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "Phylum", "mean")
conts.pMeans.age_quadratic <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "contVar", "mean")

gen.num_sig.age_quadratic   <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "Genus", "num_sig")
phy.num_sig.age_quadratic   <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "Phylum", "num_sig")
conts.num_sig.age_quadratic <- get_subs_anova_per_taxon(age_quadratic.SubsTests, "Age", "contVar", "num_sig")


# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.aq <- gen.num_sig.age_quadratic[ gen.num_sig.age_quadratic >= 70]
sort(gen.pMeans.age_quadratic[ names(freqSig.aq) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.aq)))

gensToCheck <- unique(names(freqSig.aq))
# as.data.frame(cbind(gen.num_sig.age_quadratic[ gensToCheck ], formatC(gen.pMeans.age_quadratic[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.aq <- as.data.frame(cbind(gen.num_sig.age_quadratic[ gensToCheck ], gen.pMeans.age_quadratic[ gensToCheck ]))
colnames(freq.mean.ps.aq) <- c("num_sig","meanP")
freq.mean.ps.aq <- freq.mean.ps.aq[ rev(order(freq.mean.ps.aq$num_sig, -freq.mean.ps.aq$meanP)), ]
freq.mean.ps.aq


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.aq))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps.aq), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1, useLoess = T)

phy.sig.order <- phy.num_sig.age_quadratic[ rev(order(phy.num_sig.age_quadratic, -phy.pMeans.age_quadratic)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_quadratic[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)


# ************************************************************************************** #
# ****************************************************************************************************************** #





# # ******************************************************************* #
# age_cubic.SubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                           trait = "Age", chosenControls = all_subSamps$Age_bins,
#                                                           polyRegression = "Age", polyKnots = "NULL", polyDegree = 3,
#                                                           anova_only = T)
# 
# saveRDS(age_cubic.SubsTests, file = sprintf("%s/R_objects/age_cubic.SubsTests.rds", p2_dir))
# # ************************************************************************ #

age_cubic.SubsTests <- readRDS(sprintf("%s/R_objects/age_cubic.SubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.age_cubic   <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "Genus", "mean")
phy.pMeans.age_cubic   <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "Phylum", "mean")
conts.pMeans.age_cubic <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "contVar", "mean")

gen.num_sig.age_cubic   <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "Genus", "num_sig")
phy.num_sig.age_cubic   <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "Phylum", "num_sig")
conts.num_sig.age_cubic <- get_subs_anova_per_taxon(age_cubic.SubsTests, "Age", "contVar", "num_sig")

# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ac <- gen.num_sig.age_cubic[ gen.num_sig.age_cubic >= 70]
sort(gen.pMeans.age_cubic[ names(freqSig.ac) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ac)))

gensToCheck <- unique(names(freqSig.ac))
# as.data.frame(cbind(gen.num_sig.age_cubic[ gensToCheck ], formatC(gen.pMeans.age_cubic[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.ac <- as.data.frame(cbind(gen.num_sig.age_cubic[ gensToCheck ], gen.pMeans.age_cubic[ gensToCheck ]))
colnames(freq.mean.ps.ac) <- c("num_sig","meanP")
freq.mean.ps.ac <- freq.mean.ps.ac[ rev(order(freq.mean.ps.ac$num_sig, -freq.mean.ps.ac$meanP)), ]
freq.mean.ps.ac


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.ac))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_cubic.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps.ac), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cubic.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1, useLoess = T)

phy.sig.order <- phy.num_sig.age_cubic[ rev(order(phy.num_sig.age_cubic, -phy.pMeans.age_cubic)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_cubic[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cubic.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cubic.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)


# ************************************************************************************** #
# ****************************************************************************************************************** #




# # ******************************************************************* #
# age_exp.SubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                         trait = "Age", chosenControls = all_subSamps$Age_bins,
#                                                         expModel = "Age",
#                                                         anova_only = T)
# 
# saveRDS(age_exp.SubsTests, file = sprintf("%s/R_objects/age_exp.SubsTests.rds", p2_dir))
# # ************************************************************************ #

age_exp.SubsTests <- readRDS(sprintf("%s/R_objects/age_exp.SubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.age_exp   <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "Genus", "mean")
phy.pMeans.age_exp   <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "Phylum", "mean")
conts.pMeans.age_exp <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "contVar", "mean")

gen.num_sig.age_exp   <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "Genus", "num_sig")
phy.num_sig.age_exp   <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "Phylum", "num_sig")
conts.num_sig.age_exp <- get_subs_anova_per_taxon(age_exp.SubsTests, "Age", "contVar", "num_sig")


# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ae <- gen.num_sig.age_exp[ gen.num_sig.age_exp >= 70]
sort(gen.pMeans.age_exp[ names(freqSig.ae) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ae)))

gensToCheck <- unique(names(freqSig.ae))
# as.data.frame(cbind(gen.num_sig.age_exp[ gensToCheck ], formatC(gen.pMeans.age_exp[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.ae <- as.data.frame(cbind(gen.num_sig.age_exp[ gensToCheck ], gen.pMeans.age_exp[ gensToCheck ]))
colnames(freq.mean.ps.ae) <- c("num_sig","meanP")
freq.mean.ps.ae <- freq.mean.ps.ae[ rev(order(freq.mean.ps.ae$num_sig, -freq.mean.ps.ae$meanP)), ]
freq.mean.ps.ae


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.ae))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_exp.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps.ae), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_exp.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1, useLoess = T)

phy.sig.order <- phy.num_sig.age_exp[ rev(order(phy.num_sig.age_exp, -phy.pMeans.age_exp)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_exp[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_exp.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_exp.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)


# ************************************************************************************** #
# ****************************************************************************************************************** #





# # ******************************************************************* #
# age_log.SubsTests <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                         trait = "Age", chosenControls = all_subSamps$Age_bins,
#                                                         logModel = "Age",
#                                                         anova_only = T)
# # 
# saveRDS(age_log.SubsTests, file = sprintf("%s/R_objects/age_log.SubsTests.rds", p2_dir))
# # ************************************************************************ #

age_log.SubsTests <- readRDS(sprintf("%s/R_objects/age_log.SubsTests.rds", p2_dir))
# ******************************************************************* #


gen.pMeans.age_log   <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "Genus", "mean")
phy.pMeans.age_log   <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "Phylum", "mean")
conts.pMeans.age_log <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "contVar", "mean")

gen.num_sig.age_log   <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "Genus", "num_sig")
phy.num_sig.age_log   <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "Phylum", "num_sig")
conts.num_sig.age_log <- get_subs_anova_per_taxon(age_log.SubsTests, "Age", "contVar", "num_sig")


# ************************************************************************************** #
# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.al <- gen.num_sig.age_log[ gen.num_sig.age_log >= 70]
sort(gen.pMeans.age_log[ names(freqSig.al) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.al)))

gensToCheck <- unique(names(freqSig.al))
# as.data.frame(cbind(gen.num_sig.age_log[ gensToCheck ], formatC(gen.pMeans.age_log[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps.al <- as.data.frame(cbind(gen.num_sig.age_log[ gensToCheck ], gen.pMeans.age_log[ gensToCheck ]))
colnames(freq.mean.ps.al) <- c("num_sig","meanP")
freq.mean.ps.al <- freq.mean.ps.al[ rev(order(freq.mean.ps.al$num_sig, -freq.mean.ps.al$meanP)), ]
freq.mean.ps.al


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.al))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_log.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps.al), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_log.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1, useLoess = T)

phy.sig.order <- phy.num_sig.age_log[ rev(order(phy.num_sig.age_log, -phy.pMeans.age_log)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.age_log[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_log.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_log.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1, useLoess = T)


# ************************************************************************************** #
# ****************************************************************************************************************** #












# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####


# Compare ANOSIM results with distances for various variables ####


# healthySamps <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No", ])
# meta.healthy <- SLL2.meta[ healthySamps, ]
# phy.healthy  <- prune_samples( healthySamps, SLL2)
# 
# # ords.healthy <- subsampling_ordination_objects(healthySamps, phy.healthy, dists_and_ps_Only = T)
# # saveRDS(ords.healthy, sprintf("%s/R_objects/SLL2.ords.healthy.rds", p2_dir))
# ords.healthy <- readRDS(sprintf("%s/R_objects/SLL2.ords.healthy.rds", p2_dir))
# 
# for (dist_meas in c("Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard","Aitchison")) {
#   
#   glomTab <- gloms_clr
#   
#   print(dist_meas)
#   clu <- get_clusters( phy.healthy, "Species", dist_meas, ords.healthy$ps.dists, glomTab, 
#                        subPops = T, subPop_Bdivs = ords.healthy)
#   
#   #add stomatotype to sample_data
#   meta.healthy[ names(clu$cluster_full), sprintf("healthy.%s_%s", clu$clus_var_name, dist_meas) ] <- clu$cluster_full#as.factor(cluster_full)
# }






healthy.nonBottle <- healthySamps[ ! is.na(SLL2.meta[healthySamps, "Water_type_home"]) &
                                     SLL2.meta[healthySamps, "Water_type_home"] != "Embotellada" ]
meta.nonBotHeal <- SLL2.meta[ healthy.nonBottle, ]
phy.nonBotHeal  <- prune_samples( healthy.nonBottle, SLL2)






anosim_disorder <- c("Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes","Circulatory_issues","Hypertension","Headaches",
                     "Lactose_intolerant","Anemia","Lung_issues","Kidney_issues","Hypothyroidism")


# for each binary variable that needs subsampling, will run that 100 times, and do anosim for each

# bin_vars <- c("Gender","Smoker","Braces","Mouth_piercing","Fluoride_toothpaste","Fluoride_supplement","Mouth_wounds",
#               "Chronic_disorder","Celiac","Cystic_fibrosis","Downs_Syndrome","Gingivitis_periodontitis",
#               "Eating_disorder","Medications","Antibiotics","Analgesics","Vitamin_supplements",
#               "Asthma","Wheezing",
#               "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen","Allergy.Animals",
#               "Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
#               "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary",
#               "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
#               "Kissing_partner",
#               "MALDI.Yeast_detected","MALDI.Mold_detected",
#               additional_diseases,"seqGroup")


bin_vars.short <- c("Gender","Smoker",
                    "Braces.binary","Mouth_piercing",
                    "Fluoride_toothpaste","Fluoride_supplement","Mouth_wounds.binary",
                    "Antibiotics","Analgesics","Vitamin_supplements","Asthma","Wheezing","Allergy",
                    "Bite_nails","Hair_in_mouth",
                    "Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
                    "Kissing_partner","MALDI.Yeast_detected","MALDI.Mold_detected","seqGroup")

bin_vars.temp <- c("Smoker","Braces.binary","Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes",
                   "Fluoride_toothpaste","Antibiotics","Wash_hands_after_bathroom","Do_you_feel_well",
                   "Circulatory_issues","Hypertension","Headaches",
                   "Lactose_intolerant","Anemia","Lung_issues","Kidney_issues","Hypothyroidism",
                   "Vitamin_supplements","MALDI.Yeast_detected","MALDI.Mold_detected")

no_sub_samp <- c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner","seqGroup")#,"Allergy","Wash_hands_before_eat"

non_bin_vars <- c("Age_groups","Water_type_home","BMI_group","BMI_official","How_do_you_feel")
Age_groups
Water_type_home
BMI_group - remove "Unknown"
BMI_official - remove "Unknown"
How_do_you_feel - make a binary version as well with Bien and not Bien

# anosim_subSamps <- list()
# # for (binVar in bin_vars.short) {
# for (binVar in bin_vars.temp) {
#   print(binVar)
#   if (binVar %in% no_sub_samp) {
#     # run anosim once, since no need to subsample
#     anosim_subSamps[[ binVar ]] <- ""
# 
#   } else if (binVar %in% anosim_disorder) {
#     # run 100 subsamplings with these variables, get an average pval for anosim later
#     anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, SLL2, SLL2.meta, 100)
#   } else {
#     # run 100 subsamplings with these variables, get an average pval for anosim later
#     anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, phy.healthy, meta.healthy, 100)
#   }
#   saveRDS(anosim_subSamps[[ binVar ]], file = sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, binVar))
# 
# }
# saveRDS(anosim_subSamps, file = sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))
# 
# 
# for (binVar in no_sub_samp) {
#   print(binVar)
#   
#   anosim_subSamps[[ binVar ]] <- list()
#   anosim_subSamps[[ binVar ]][[ "1" ]] <- list()
#   
#   bvs <- rownames(meta.healthy[ ! is.na(meta.healthy[,binVar]), ])
#   anosim_subSamps[[ binVar ]][[ "1" ]][[ "samples" ]] <- bvs
#   # anosim_subSamps[[ binVar ]][[ "1" ]][[ "ordObjects" ]] <- ords.healthy$Aitchison
#   if (binVar %in% c("seqGroup")) { # and will have to add some disorder groups too, since controls only taken from seqGroup "One"
#     anosim_subSamps[[ binVar ]][[ "1" ]][[ "ANOSIM" ]] <- anosim( as.dist(as.matrix(ords.healthy$Aitchison)[bvs, bvs]),
#                                                                  as.character(meta.healthy[ bvs, binVar]))
#   } else {
#     anosim_subSamps[[ binVar ]][[ "1" ]][[ "ANOSIM" ]] <- anosim( as.dist(as.matrix(ords.healthy$Aitchison)[bvs, bvs]),
#                                                                  as.character(meta.healthy[ bvs, binVar]),
#                                                                  strata = meta.healthy[ bvs, "seqGroup"])
#   }
#   
# }



# for (binVar in c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner")) {
#   print(binVar)
#   # run 100 subsamplings with these variables, get an average pval for anosim later
#   anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, phy.healthy, meta.healthy, 100, YesNo.add = T)
#   saveRDS(anosim_subSamps[[ binVar ]], file = sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, binVar))
# }
# saveRDS(anosim_subSamps, file = sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))


# for (binVar in c("MALDI.Yeast_detected","MALDI.Mold_detected")) {
#   print(binVar)
#   # run 100 subsamplings with these variables, get an average pval for anosim later
#   anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, phy.healthy, meta.healthy, 100)
#   saveRDS(anosim_subSamps[[ binVar ]], file = sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, binVar))
# }
# saveRDS(anosim_subSamps, file = sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))


# for (binVar in c("Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats")) {
#   print(binVar)
#   # run 100 subsamplings with these variables, get an average pval for anosim later
#   anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, phy.healthy, meta.healthy, 100, YesNo.add = T)
#   saveRDS(anosim_subSamps[[ binVar ]], file = sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, binVar))
# }
# saveRDS(anosim_subSamps, file = sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))

# for (binVar in c("Water_type_home")) {
#   print(binVar)
#   # run 100 subsamplings with these variables, get an average pval for anosim later
#   anosim_subSamps[[ binVar ]] <- run_subsampling_anosims(binVar, phy.healthy, meta.healthy, 100, waterType = T)
#   saveRDS(anosim_subSamps[[ binVar ]], file = sprintf("%s/R_objects/anosim_subSamps.temp.%s.rds", p2_dir, binVar))
# }
# saveRDS(anosim_subSamps, file = sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))



anosim_subSamps <- readRDS(sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))
anosim_subSamps.sigs <- readRDS(sprintf("%s/R_objects/anosim_subSamps.sigs.rds", p2_dir))



aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(names(anosim_subSamps)), "Binary_variables",
                   plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps)


pag.health <- c("Smoker","Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes",
                "Circulatory_issues","Hypertension","Headaches","Lactose_intolerant",
                "Anemia","Lung_issues","Kidney_issues","Hypothyroidism")
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.health), "Binary_variables",
                   plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps)
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.health), "Binary_variables",
                   plotType="box", whichGetFunc = "subSamps", ast=anosim_subSamps)



pag.misc.1 <- c("Gender","Braces.binary","Fluoride_toothpaste","Wash_hands_after_bathroom",
                "Mouth_wounds.binary","Bite_nails","Chew_pens","Kissing_partner")
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.misc.1), "Binary_variables",
                   plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps)
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.misc.1), "Binary_variables",
                   plotType="box", whichGetFunc = "subSamps", ast=anosim_subSamps)


pag.misc.2 <- c("Antibiotics","Analgesics","Vitamin_supplements","MALDI.Yeast_detected","MALDI.Mold_detected")
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.misc.2), "Binary_variables",
                   plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps)
plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(pag.misc.2), "Binary_variables",
                   plotType="box", whichGetFunc = "subSamps", ast=anosim_subSamps)




plot_anosim_groups(SLL2.meta, as.matrix(aitch), "Aitchison", "ANOSIM",
                   rev(c(pag.misc.1,pag.misc.2)), "Binary_variables",
                   plotType="box", whichGetFunc = "subSamps", ast=anosim_subSamps)


# **************** #


par(mfrow=c(4,7), mar=c(2,2,2,2))

for (col_to_check in names(anosim_subSamps)) {
  lvs <- get_anosim_subSamps_plot_vals.alt(anosim_subSamps, SLL2.meta, aitch, col_to_check, "Binary_variables", plotType="box")
  
  meanP <- mean(p.adjust( unlist(lapply(anosim_subSamps[[ col_to_check ]], function(x) x$ANOSIM$signif)), method = "fdr"))
  meanR <- round(mean( unlist(lapply(anosim_subSamps[[ col_to_check ]], function(x) x$ANOSIM$statistic)) ), 3)
  
  printR <- ifelse(meanP < 0.0005, sprintf("R = %s ***", meanR), 
                   ifelse(meanP < 0.005, sprintf("R = %s **", meanR), 
                          ifelse(meanP < 0.05, sprintf("R = %s *", meanR), 
                                 ifelse(meanP < 0.1, sprintf("R = %s ·", meanR), ""))))
  
  # print(c(col_to_check, printR, meanP))
  boxplot(lvs, main=col_to_check, col=c("red","blue","green"), horizontal=T )
  text(x=55, y=1.45, labels = printR)
}

# **************** #

par(mfrow=c(2,5), mar=c(2, 2.2, 1.2, 1), cex.lab=1.5, cex.axis=1.5)

cols <- c("Downs_Syndrome","Cystic_fibrosis","Celiac","Hypertension","Smoker",
          "MALDI.Yeast_detected","Gender","Antibiotics","Fluoride_toothpaste","Braces.binary")

for (col_to_check in cols) {
  lvs <- get_anosim_subSamps_plot_vals.alt(anosim_subSamps.sigs, SLL2.meta, aitch, col_to_check, "Binary_variables", plotType="box")
  
  meanP <- mean(p.adjust( unlist(lapply(anosim_subSamps.sigs[[ col_to_check ]], function(x) x$ANOSIM$signif)), method = "fdr"))
  meanR <- round(mean( unlist(lapply(anosim_subSamps.sigs[[ col_to_check ]], function(x) x$ANOSIM$statistic)) ), 3)
  
  printR <- ifelse(meanP < 0.0005, sprintf("R = %s ***", meanR), 
                   ifelse(meanP < 0.005, sprintf("R = %s **", meanR), 
                          ifelse(meanP < 0.05, sprintf("R = %s *", meanR), 
                                 ifelse(meanP < 0.1, sprintf("R = %s ·", meanR), ""))))
  printP <- sprintf("p.adj = %s", round(meanP, 4))
  # print(c(col_to_check, printR, meanP))
  # boxplot(lvs, main=col_to_check, col=c("#00CED1","#CD5C5C","#3CB371"), 
  #         horizontal=T )
  # text(x=55, y=1.45, labels = printR)
  
  if (col_to_check %in% cols[1:5]) xl <- c("","","")
  else xl <- c("Yes/F","No/M","Between")
  
  if (col_to_check %in% cols[c(1,6)]) yl <- c("20","30","40","50","60","70")
  else yl <- c("","","","","","")
  
  boxplot(lvs, main=col_to_check, col=c("#00CED1","#CD5C5C","#3CB371"), 
          horizontal=F, ylim=c(15,70), xaxt="n", yaxt="n", 
          cex.main=1.5, frame=F)
  axis(side=1, at=c(1,2,3), labels=xl)
  axis(side=2, at=c(20,30,40,50,60,70), labels=yl)
  text(x=2, y=70, labels = printR, cex=1.75, col="red")
  text(x=2, y=67, labels = printP, cex=1.25, col="black")
  box(bty="l")
}

# **************** #







Smoker
Stomatotype
Gender
Braces.binary
Mouth_piercing
Fluoride_toothpaste - change Yes/No value in subsampling function ***
Fluoride_supplement
Mouth_wounds - make binary
Medications
Antibiotics
Analgesics
Vitamin_supplements
Other_medications_binary
Asthma
Wheezing

Family_unit

DS
CF
Chronic_disorder
Diabetes


Age_groups
Water_type_home
BMI_group - remove "Unknown"
BMI_official - remove "Unknown"
How_do_you_feel - make a binary version as well with Bien and not Bien

# ****************************************************************************************************************** #




# For Age_groups ####
ageGroupSubsTests.TAS <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.TAS.rds", p2_dir))

mTab.ageGroups <- meta.healthy[ unique(unlist(lapply(ageGroupSubsTests.TAS, function(x) x$samples))), ]
mTab.ageGroups$Children <- ifelse(mTab.ageGroups$Age_groups == "Child", "Yes", "No")
mTab.ageGroups$Teens    <- ifelse(mTab.ageGroups$Age_groups == "Teen", "Yes", "No")
mTab.ageGroups$Adults   <- ifelse(mTab.ageGroups$Age_groups == "Adult", "Yes", "No")
mTab.ageGroups$Seniors  <- ifelse(mTab.ageGroups$Age_groups == "Senior", "Yes", "No")

phy.ageGroups <- phyloseq(otu_table(prune_samples(rownames(mTab.ageGroups), SLL2)),
                          sample_data(mTab.ageGroups),
                          phy_tree(prune_samples(rownames(mTab.ageGroups), SLL2)),
                          tax_table(prune_samples(rownames(mTab.ageGroups), SLL2)))

# anosim_subSamps.ageTAS <- list()
# for (ag in c("Teens","Adults","Seniors")) {
#   print(ag)
#   anosim_subSamps.ageTAS[[ ag ]] <- run_subsampling_anosims(ag, phy.ageGroups, mTab.ageGroups, 100, chosenControls = ageGroupSubsTests.TAS)
# }
# saveRDS(anosim_subSamps.ageTAS, sprintf("%s/R_objects/anosim_subSamps.ageTAS.rds", p2_dir))

anosim_subSamps.ageTAS <- readRDS(sprintf("%s/R_objects/anosim_subSamps.ageTAS.rds", p2_dir))

ano_ss.aTAS.means <- t(as.data.frame(list("meanP" = unlist(lapply(anosim_subSamps.ageTAS, function(x) 
  mean(p.adjust(unlist(lapply(x, function(y) y$ANOSIM$signif)), method = "fdr")))),
  "meanR" = unlist(lapply(anosim_subSamps.ageTAS, function(x) 
    mean(unlist(lapply(x, function(y) y$ANOSIM$statistic)))))
)))
ano_ss.aTAS.means

plot_anosim_groups(mTab.ageGroups, as.matrix(aitch), "Aitchison", "ANOSIM",
                   c("Seniors","Adults","Teens"), "Age_groups",
                   plotType="bar", whichGetFunc = "subSamps", ast=anosim_subSamps.ageTAS)

# ********************* #
# ********************* #
# ********************* #


# ********************* #

ag_mean_diffs_prev <- as.data.frame(t(sapply(c("Teens","Adults","Seniors"), function(v)
  sapply(c("Yes","No"), function(gVal) 
    mean(get_dists_per_group(v, gVal, mTab.ageTAS, aitch, allSubs = lapply(anosim_subSamps.ageTAS[[ v ]], function(x) x$samples)))
  ))))
ag_mean_diffs_prev$Diff <- ag_mean_diffs_prev$Yes - ag_mean_diffs_prev$No
# ag_mean_diffs_prev <- ag_mean_diffs_prev[ rev(order(abs(ag_mean_diffs_prev$Diff))), ]

ag_median_diffs_prev <- as.data.frame(t(sapply(c("Teens","Adults","Seniors"), function(v)
  sapply(c("Yes","No"), function(gVal) 
    mean(get_dists_per_group(v, gVal, mTab.ageTAS, aitch, allSubs = lapply(anosim_subSamps.ageTAS[[ v ]], function(x) x$samples)))
  ))))
ag_median_diffs_prev$Diff <- ag_median_diffs_prev$Yes - ag_median_diffs_prev$No
# ********************* #



# ********************* #
ag_mean_diffs <- as.data.frame((sapply(c("Teens","Adults","Seniors"), function(v)
  lapply(anosim_subSamps.ageTAS[[ v ]], function(i) 
    sapply(c("Yes","No"), function(gVal) 
      # lapply(anosim_subSamps.ageTAS[[ v ]], function(i) {
      mean(get_dists_per_group(v, gVal, mTab.ageTAS, aitch, chosenSamps = i$samples))
    )
  )
)))

amd <- as.data.frame(matrix(NA, nrow=3, ncol=3))
rownames(amd) <- c("Teens","Adults","Seniors")
colnames(amd) <- c("Yes","No", "Diff")
amd$Yes  <- unlist(lapply(ag_mean_diffs, function(x) mean(unlist(lapply(x, function(y) y["Yes"])))))
amd$No   <- unlist(lapply(ag_mean_diffs, function(x) mean(unlist(lapply(x, function(y) y["No"])))))
amd$Diff <- unlist(lapply(ag_mean_diffs, function(x) mean(unlist(lapply(x, function(y) y["Yes"] - y["No"])))))

# ********************* #

ag_median_diffs <- as.data.frame((sapply(c("Teens","Adults","Seniors"), function(v)
  lapply(anosim_subSamps.ageTAS[[ v ]], function(i) 
    sapply(c("Yes","No"), function(gVal) 
      # lapply(anosim_subSamps.ageTAS[[ v ]], function(i) {
      median(get_dists_per_group(v, gVal, mTab.ageTAS, aitch, chosenSamps = i$samples))
    )
  )
)))

amdd <- as.data.frame(matrix(NA, nrow=3, ncol=3))
rownames(amdd) <- c("Teens","Adults","Seniors")
colnames(amdd) <- c("Yes","No", "Diff")
amdd$Yes  <- unlist(lapply(ag_median_diffs, function(x) mean(unlist(lapply(x, function(y) y["Yes"])))))
amdd$No   <- unlist(lapply(ag_median_diffs, function(x) mean(unlist(lapply(x, function(y) y["No"])))))
amdd$Diff <- unlist(lapply(ag_median_diffs, function(x) mean(unlist(lapply(x, function(y) y["Yes"] - y["No"])))))

# ********************* #




aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

mTab.ageTAS <- meta.healthy[ unique(unlist(lapply(anosim_subSamps.ageTAS$Teens, function(x) x$samples))), ]
mTab.ageTAS$Teens    <- ifelse(mTab.ageTAS$Age_groups == "Teen", "Yes", "No")
mTab.ageTAS$Adults   <- ifelse(mTab.ageTAS$Age_groups == "Adult", "Yes", "No")
mTab.ageTAS$Seniors  <- ifelse(mTab.ageTAS$Age_groups == "Senior", "Yes", "No")

phy.ageTAS <- phyloseq(otu_table(prune_samples(rownames(mTab.ageTAS), SLL2)),
                       sample_data(mTab.ageTAS),
                       phy_tree(prune_samples(rownames(mTab.ageTAS), SLL2)),
                       tax_table(prune_samples(rownames(mTab.ageTAS), SLL2)))


agap.mean   <- get_Age_anosim_plotVals(c("Teens","Adults","Seniors"), "mean", mTab.ageTAS, anosim_subSamps.ageTAS, aitch)
agap.median <- get_Age_anosim_plotVals(c("Teens","Adults","Seniors"), "median", mTab.ageTAS, anosim_subSamps.ageTAS, aitch)

plot_Age_anosim_plotVals(agap.mean, "mean", "-log(anosim.p.adj)")
plot_Age_anosim_plotVals(agap.mean, "mean", "anosim.R")

plot_Age_anosim_plotVals(agap.median, "median", "-log(anosim.p.adj)")
plot_Age_anosim_plotVals(agap.median, "median", "anosim.R")




# TAS.subSamps <- lapply(anosim_subSamps.ageTAS$Teens, function(x) x$samples)
# ageGroup.pairwise.anosim_subs <- run_pairwise_subsampling_calcs.ageGroups("healthy", phy.ageTAS, mTab.ageTAS, 100, "ANOSIM",
#                                                                           chosenControls = TAS.subSamps, TAS = T)
# saveRDS(ageGroup.pairwise.anosim_subs, file = sprintf("%s/R_objects/ageGroup.pairwise.anosim_subs.rds", p2_dir))

ageGroup.pairwise.anosim_subs <- readRDS(sprintf("%s/R_objects/ageGroup.pairwise.anosim_subs.rds", p2_dir))


ag.pw.mean   <- get_Age_pairwise_anosim_plotVals("Age_groups", "mean", mTab.ageTAS, ageGroup.pairwise.anosim_subs, aitch)
ag.pw.median <- get_Age_pairwise_anosim_plotVals("Age_groups", "median", mTab.ageTAS, ageGroup.pairwise.anosim_subs, aitch)

plot_Age_anosim_plotVals(ag.pw.mean, "mean", "-log(anosim.p.adj)", paired = T)
plot_Age_anosim_plotVals(ag.pw.mean, "mean", "anosim.R", paired = T)

plot_Age_anosim_plotVals(ag.pw.median, "median", "-log(anosim.p.adj)", paired = T)
plot_Age_anosim_plotVals(ag.pw.median, "median", "anosim.R", paired = T)



# ********************* #
# ********************* #
# ********************* #

bd.ageTAS <- lapply(ageGroupSubsTests.TAS, function(x)
  betadisper(as.dist(aitch[x$samples,x$samples]), mTab.ageTAS[x$samples, "Age_groups"]) )

bd.ageTAS.pmean <- mean(p.adjust(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageTAS.pmedian <- median(p.adjust(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageTAS.psd <- sd(p.adjust(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))

bd.ageTAS.Fmean <- mean(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","F value"])))
bd.ageTAS.Fmedian <- median(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","F value"])))
bd.ageTAS.Fsd <- sd(unlist(lapply(bd.ageTAS, function(x) anova(x)["Groups","F value"])))

bd.ageTAS.dists <- list("Teen"   = unlist(lapply(bd.ageTAS, function(x) x$distances[x$group=="Teen"])),
                        "Adult"  = unlist(lapply(bd.ageTAS, function(x) x$distances[x$group=="Adult"])),
                        "Senior" = unlist(lapply(bd.ageTAS, function(x) x$distances[x$group=="Senior"])))

boxplot(bd.ageTAS.dists)

# ****************************************************************************************************************** #




# ****************************************************************************************************************** #


# For Age_bins ####
mTab.ageBins <- meta.healthy[ unique(unlist(lapply(all_subSamps$Age_bins, function(x) x))), ]
mTab.ageBins$`13_20` <- ifelse(mTab.ageBins$Age_bins == "13_20", "Yes", "No")
mTab.ageBins$`20_30` <- ifelse(mTab.ageBins$Age_bins == "20_30", "Yes", "No")
mTab.ageBins$`30_40` <- ifelse(mTab.ageBins$Age_bins == "30_40", "Yes", "No")
mTab.ageBins$`40_50` <- ifelse(mTab.ageBins$Age_bins == "40_50", "Yes", "No")
mTab.ageBins$`50_60` <- ifelse(mTab.ageBins$Age_bins == "50_60", "Yes", "No")
mTab.ageBins$`60+`   <- ifelse(mTab.ageBins$Age_bins == "60+", "Yes", "No")

phy.ageBins <- phyloseq(otu_table(prune_samples(rownames(mTab.ageBins), SLL2)),
                          sample_data(mTab.ageBins),
                          phy_tree(prune_samples(rownames(mTab.ageBins), SLL2)),
                          tax_table(prune_samples(rownames(mTab.ageBins), SLL2)))

# anosim_subSamps.ageBins <- list()
# for (ag in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ag)
#   anosim_subSamps.ageBins[[ ag ]] <- run_subsampling_anosims(ag, phy.ageBins, mTab.ageBins, 100, chosenControls = ageBinsSubsTests)
# }
# saveRDS(anosim_subSamps.ageBins, sprintf("%s/R_objects/anosim_subSamps.ageBins.rds", p2_dir))

anosim_subSamps.ageBins <- readRDS(sprintf("%s/R_objects/anosim_subSamps.ageBins.rds", p2_dir))

ano_ss.ageBins.means <- t(as.data.frame(list("meanP" = unlist(lapply(anosim_subSamps.ageBins, function(x) 
  mean(p.adjust(unlist(lapply(x, function(y) y$ANOSIM$signif)), method = "fdr")))),
  "meanR" = unlist(lapply(anosim_subSamps.ageBins, function(x) 
    mean(unlist(lapply(x, function(y) y$ANOSIM$statistic)))))
)))
ano_ss.ageBins.means

plot_anosim_groups(mTab.ageBins, as.matrix(ordObj.bins$Aitchison), "Aitchison", "ANOSIM",
                   rev(c("13_20","20_30","30_40","40_50","50_60","60+")), "Age_bins",
                   plotType="box", whichGetFunc = "subSamps", ast=anosim_subSamps.ageBins)






# ********************* #
mTab.ageBins <- meta.healthy[ unique(unlist(lapply(anosim_subSamps.ageBins$`13_20`, function(x) x$samples))), ]
mTab.ageBins$`13_20` <- ifelse(mTab.ageBins$Age_bins == "13_20", "Yes", "No")
mTab.ageBins$`20_30` <- ifelse(mTab.ageBins$Age_bins == "20_30", "Yes", "No")
mTab.ageBins$`30_40` <- ifelse(mTab.ageBins$Age_bins == "30_40", "Yes", "No")
mTab.ageBins$`40_50` <- ifelse(mTab.ageBins$Age_bins == "40_50", "Yes", "No")
mTab.ageBins$`50_60` <- ifelse(mTab.ageBins$Age_bins == "50_60", "Yes", "No")
mTab.ageBins$`60+`   <- ifelse(mTab.ageBins$Age_bins == "60+", "Yes", "No")

phy.ageBins <- phyloseq(otu_table(prune_samples(rownames(mTab.ageBins), SLL2)),
                        sample_data(mTab.ageBins),
                        phy_tree(prune_samples(rownames(mTab.ageBins), SLL2)),
                        tax_table(prune_samples(rownames(mTab.ageBins), SLL2)))


abap.mean   <- get_Age_anosim_plotVals(c("13_20","20_30","30_40","40_50","50_60","60+"), 
                                       "mean", mTab.ageBins, anosim_subSamps.ageBins, aitch)
abap.median <- get_Age_anosim_plotVals(c("13_20","20_30","30_40","40_50","50_60","60+"), 
                                       "median", mTab.ageBins, anosim_subSamps.ageBins, aitch)

plot_Age_anosim_plotVals(abap.mean, "mean", "-log(anosim.p.adj)")
plot_Age_anosim_plotVals(abap.mean, "mean", "anosim.R")

plot_Age_anosim_plotVals(abap.median, "median", "-log(anosim.p.adj)")
plot_Age_anosim_plotVals(abap.median, "median", "anosim.R")




# ageBins.subSamps <- lapply(anosim_subSamps.ageBins$`13_20`, function(x) x$samples)
# ageBins.pairwise.anosim_subs <- run_pairwise_subsampling_calcs.ageBins("healthy", phy.ageBins, mTab.ageBins, 100, "ANOSIM",
#                                                                        chosenControls = ageBins.subSamps)
# saveRDS(ageBins.pairwise.anosim_subs, file = sprintf("%s/R_objects/ageBins.pairwise.anosim_subs.rds", p2_dir))

ageBins.pairwise.anosim_subs <- readRDS(sprintf("%s/R_objects/ageBins.pairwise.anosim_subs.rds", p2_dir))


ab.pw.mean   <- get_Age_pairwise_anosim_plotVals("Age_bins", "mean", mTab.ageBins, ageBins.pairwise.anosim_subs, aitch)
ab.pw.median <- get_Age_pairwise_anosim_plotVals("Age_bins", "median", mTab.ageBins, ageBins.pairwise.anosim_subs, aitch)

plot_Age_anosim_plotVals(ab.pw.mean, "mean", "-log(anosim.p.adj)", paired = T)
plot_Age_anosim_plotVals(ab.pw.mean, "mean", "anosim.R", paired = T)

plot_Age_anosim_plotVals(ab.pw.median, "median", "-log(anosim.p.adj)", paired = T)
plot_Age_anosim_plotVals(ab.pw.median, "median", "anosim.R", paired = T)


# ********************* #
# ********************* #
# ********************* #
# ****************************************************************************************************************** #
# homogeneity tests ####

mTab.ageBins <- meta.healthy[ unique(unlist(lapply(all_subSamps$Age_bins, function(x) x))), ]
mTab.ageBins$`13_20` <- ifelse(mTab.ageBins$Age_bins == "13_20", "Yes", "No")
mTab.ageBins$`20_30` <- ifelse(mTab.ageBins$Age_bins == "20_30", "Yes", "No")
mTab.ageBins$`30_40` <- ifelse(mTab.ageBins$Age_bins == "30_40", "Yes", "No")
mTab.ageBins$`40_50` <- ifelse(mTab.ageBins$Age_bins == "40_50", "Yes", "No")
mTab.ageBins$`50_60` <- ifelse(mTab.ageBins$Age_bins == "50_60", "Yes", "No")
mTab.ageBins$`60+`   <- ifelse(mTab.ageBins$Age_bins == "60+", "Yes", "No")

phy.ageBins <- phyloseq(otu_table(prune_samples(rownames(mTab.ageBins), SLL2)),
                        sample_data(mTab.ageBins),
                        phy_tree(prune_samples(rownames(mTab.ageBins), SLL2)),
                        tax_table(prune_samples(rownames(mTab.ageBins), SLL2)))

aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

bd.ageBins <- lapply(all_subSamps$Age_bins, function(x)
  betadisper(as.dist(aitch[x,x]), mTab.ageBins[x, "Age_bins"]) )

bd.ageBins.pmean <- mean(p.adjust(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.pmedian <- median(p.adjust(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.psd <- sd(p.adjust(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))

bd.ageBins.Fmean <- mean(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","F value"])))
bd.ageBins.Fmedian <- median(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","F value"])))
bd.ageBins.Fsd <- sd(unlist(lapply(bd.ageBins, function(x) anova(x)["Groups","F value"])))

bd.ageBins.dists <- list("13_20" = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="13_20"])),
                         "20_30" = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="20_30"])),
                         "30_40" = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="30_40"])),
                         "40_50" = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="40_50"])),
                         "50_60" = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="50_60"])),
                         "60+"   = unlist(lapply(bd.ageBins, function(x) x$distances[x$group=="60+"])))

boxplot(bd.ageBins.dists, notch=T)
# ********************* #
# ********************* #

# get adonis values for each individual age bin

# adonis_subSamps.ageBins_indiv <- list()
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ab)
#   adonis_subSamps.ageBins_indiv[[ ab ]] <- run_subsampling_adonis(ab, phy.ageBins, mTab.ageBins, 100, chosenControls = all_subSamps$Age_bins)
# }
# saveRDS(adonis_subSamps.ageBins_indiv, sprintf("%s/R_objects/adonis_subSamps.ageBins_indiv.rds", p2_dir))

adonis_subSamps.ageBins_indiv <- readRDS(sprintf("%s/R_objects/adonis_subSamps.ageBins_indiv.rds", p2_dir))

asabi.mean_padj <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(p.adjust(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Aitchison[x, "Pr(>F)"])), method = "fdr")))

asabi.mean_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Aitchison[x, "R2"]))))

asabi.mean_F <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Aitchison[x, "F.Model"]))))


asabi.all_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) {
    nSamps <- length(i$samples[ meta.healthy[i$samples, "Age_bins"] == x ])
    c( i$Adonis$Aitchison[x, "R2"], rep(NA, nSamps-1) )
  }))
)

# ********************* #
# ********************* #

# get betadisper results for each individual age bin against others

# betadisp.ageBins_indiv <- list()
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ab)
#   betadisp.ageBins_indiv[[ ab ]] <- lapply(all_subSamps$Age_bins, function(x) {
#     mtab.bd_bins <- meta.healthy[ x, ]
#     mtab.bd_bins$age_bins.binary <- sapply(mtab.bd_bins$Age_bins, function(y) ifelse(y==ab, "Bin", "Other"))
#     betadisper(as.dist(aitch[x,x]), mtab.bd_bins[x, "age_bins.binary"])
#   })
# }
# saveRDS(betadisp.ageBins_indiv, sprintf("%s/R_objects/betadisp.ageBins_indiv.rds", p2_dir))

betadisp.ageBins_indiv <- readRDS(sprintf("%s/R_objects/betadisp.ageBins_indiv.rds", p2_dir))

bd.ageBins_indiv.pmean <- c(
  "13_20" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`13_20`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "20_30" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`20_30`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "30_40" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`30_40`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "40_50" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`40_50`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "50_60" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`50_60`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "60+"   = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv$`60+`,   function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )
bd.ageBins_indiv.Fmean <- c(
  "13_20" = mean(unlist(lapply(betadisp.ageBins_indiv$`13_20`, function(x) anova(x)["Groups","F value"]))),
  "20_30" = mean(unlist(lapply(betadisp.ageBins_indiv$`20_30`, function(x) anova(x)["Groups","F value"]))),
  "30_40" = mean(unlist(lapply(betadisp.ageBins_indiv$`30_40`, function(x) anova(x)["Groups","F value"]))),
  "40_50" = mean(unlist(lapply(betadisp.ageBins_indiv$`40_50`, function(x) anova(x)["Groups","F value"]))),
  "50_60" = mean(unlist(lapply(betadisp.ageBins_indiv$`50_60`, function(x) anova(x)["Groups","F value"]))),
  "60+"   = mean(unlist(lapply(betadisp.ageBins_indiv$`60+`,   function(x) anova(x)["Groups","F value"]))) )

# ********************* #
# ********************* #



bdab.melt <- reshape::melt.list(bd.ageBins.dists)
colnames(bdab.melt) <- c("Distance_to_centroid","Age_bin")

bdab.melt$mean_padj <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_padj[ x ])
bdab.melt$mean_plog <- sapply(bdab.melt$Age_bin, function(x) -log(asabi.mean_padj[ x ]))
bdab.melt$mean_R2   <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_R2[ x ])
bdab.melt$mean_F    <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_F[ x ])

bdab.melt$Age_bin <- sapply(bdab.melt$Age_bin, function(x) {
  bin <- gsub("60\\+",">60",gsub("_","-",x))
  perm.padj <- ifelse(asabi.mean_padj[x] < 0.001, "< 0.001", paste0("= ",round(asabi.mean_padj[ x ], 3)))
  betd.padj <- ifelse(bd.ageBins.pmean[x] < 0.001, "< 0.001", paste0("= ", round(bd.ageBins.pmean[ x ], 3)))
  
  if (asabi.mean_padj[x] > 0.1 & bd.ageBins.pmean[x] > 0.1)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("<span style = 'color:darkgrey;'>permanova p %s</span>", perm.padj),
            sprintf("<span style = 'color:darkgrey;'>betaDisper p %s</span>", betd.padj))
  else if (asabi.mean_padj[x] > 0.05 & asabi.mean_padj[x] <= 0.1 & bd.ageBins.pmean[x] > 0.1)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("<span style = 'color:grey50;'>permanova p %s</span>", perm.padj),
            sprintf("<span style = 'color:darkgrey;'>betaDisper p %s</span>", betd.padj))
  else if (asabi.mean_padj[x] > 0.1 & bd.ageBins.pmean[x] > 0.05 & bd.ageBins.pmean[x] <= 0.1)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("<span style = 'color:darkgrey;'>permanova p %s</span>", perm.padj),
            sprintf("<span style = 'color:grey50;'>betaDisper p %s</span>", betd.padj))
  else if (asabi.mean_padj[x] > 0.05 & bd.ageBins.pmean[x] > 0.05)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("<span style = 'color:grey50;'>permanova p %s</span>", perm.padj),
            sprintf("<span style = 'color:grey50;'>betaDisper p %s</span>", betd.padj))
  else if (asabi.mean_padj[x] > 0.05 & bd.ageBins.pmean[x] < 0.05)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("<span style = 'color:grey50;'>permanova p %s</span>", perm.padj),
            sprintf("betaDisper p %s", betd.padj))
  else if (asabi.mean_padj[x] < 0.05 & bd.ageBins.pmean[x] > 0.05)
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("permanova p %s", perm.padj),
            sprintf("<span style = 'color:grey50;'>betaDisper p %s</span>", betd.padj))
  else
    sprintf("**%s**<br>%s<br>%s", bin,
            sprintf("permanova p %s", perm.padj),
            sprintf("betaDisper p %s", betd.padj))
})
bdab.melt$Age_bin <- factor(bdab.melt$Age_bin, levels = rev(unique(bdab.melt$Age_bin)))

# ggplot(bdab.melt, aes(x=Age_bin, y=Distance_to_centroid, fill=mean_plog)) +
#   geom_boxplot(notch = T) +
#   scale_fill_gradient2(low="white", high="darkblue")

ggplot(bdab.melt, aes(x=Age_bin, y=Distance_to_centroid, fill=mean_R2)) +
  geom_boxplot(notch = T) +
  scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  theme_minimal() +
  coord_flip() +
  labs(y="Distance to spatial median") +
  theme(axis.title.y = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        axis.title.x = element_text(size=17), axis.text.x = element_text(size=15),
        legend.text = element_text(size=13), legend.title = element_text(size=13))


# ****************************************************************************************************************** #

# ORIGINAL VERSION:

bdab.melt <- reshape::melt.list(bd.ageBins.dists)
colnames(bdab.melt) <- c("Distance_to_centroid","Age_bin")

bdab.melt$mean_padj <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_padj[ x ])
bdab.melt$mean_plog <- sapply(bdab.melt$Age_bin, function(x) -log(asabi.mean_padj[ x ]))
bdab.melt$mean_R2   <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_R2[ x ])
bdab.melt$mean_F    <- sapply(bdab.melt$Age_bin, function(x) asabi.mean_F[ x ])

bdab.melt$Age_bin <- sapply(bdab.melt$Age_bin, function(x) {
  bin <- gsub("60\\+",">60",gsub("_","-",x))
  padj <- round(asabi.mean_padj[ x ], 3)
  
  if (padj > 0.1)
    sprintf("**%s**<br><span style = 'color:darkgrey;'>p = %s</span>", bin, padj)
  else if (padj > 0.05)
    sprintf("**%s**<br><span style = 'color:grey40;'>p = %s</span>", bin, padj)
  else
    sprintf("**%s**<br>p = %s", bin, padj)
})
bdab.melt$Age_bin <- factor(bdab.melt$Age_bin, levels = unique(bdab.melt$Age_bin))

# ggplot(bdab.melt, aes(x=Age_bin, y=Distance_to_centroid, fill=mean_plog)) +
#   geom_boxplot(notch = T) +
#   scale_fill_gradient2(low="white", high="darkblue")

ggplot(bdab.melt, aes(x=Age_bin, y=Distance_to_centroid, fill=mean_R2)) +
  geom_boxplot(notch = T) +
  scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17), axis.text.y = element_text(size=15),
        legend.text = element_text(size=13), legend.title = element_text(size=13))

# ********** #

bdab.melt$all_R2 <- unlist(asabi.all_R2) + 0.02
# bdab.melt$all_R2    <- unlist(asabi.all_R2) + (min(bdab.melt$Distance_to_centroid) * (max(unlist(asabi.all_R2), na.rm = T) / max(bdab.melt$Distance_to_centroid)))

scaleFactor <- max(bdab.melt$Distance_to_centroid) / max(bdab.melt$all_R2, na.rm = T)

bdab.melt$Age_bin2 <- factor(bdab.melt$Age_bin, levels = rev(unique(bdab.melt$Age_bin)))

ggplot(bdab.melt, aes(y=Age_bin2)) +
  geom_boxplot(aes(x=Distance_to_centroid), width=0.4, fill="#CB6767", notch = T, position = position_nudge(y=-0.2)) +
  geom_boxplot(aes(x=all_R2 * scaleFactor), width=0.4, fill="#73BDD3", notch = T, position = position_nudge(y=0.2)) +
  # scale_x_continuous(name="Distance to spatial median", sec.axis=sec_axis(~./scaleFactor - (min(bdab.melt$Distance_to_centroid) * (max(unlist(asabi.all_R2), na.rm = T) / max(bdab.melt$Distance_to_centroid))), name="permanova R2")) +
  scale_x_continuous(name="Distance to spatial median", sec.axis=sec_axis(~./scaleFactor - 0.02, name="permanova R2")) +
  # scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        axis.title.x.bottom = element_text(size=17, color="#CB6767"), axis.text.x.bottom = element_text(size=15, color="#CB6767"),
        axis.title.x.top = element_text(size=17, color="#73BDD3"), axis.text.x.top = element_text(size=15, color="#73BDD3"),
        legend.text = element_text(size=13), legend.title = element_text(size=13))

## same but vertical
ggplot(bdab.melt, aes(x=Age_bin)) +
  geom_boxplot(aes(y=Distance_to_centroid), width=0.3, fill="#CB6767", notch = T, position = position_nudge(x=-0.15)) +
  geom_boxplot(aes(y=all_R2 * scaleFactor), width=0.3, fill="#73BDD3", notch = T, position = position_nudge(x=0.15)) +
  # scale_x_continuous(name="Distance to spatial median", sec.axis=sec_axis(~./scaleFactor - (min(bdab.melt$Distance_to_centroid) * (max(unlist(asabi.all_R2), na.rm = T) / max(bdab.melt$Distance_to_centroid))), name="permanova R2")) +
  scale_y_continuous(name="Distance to spatial median", sec.axis=sec_axis(~./scaleFactor - 0.02, name="permanova R2")) +
  # scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  # labs(x="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = ggtext::element_markdown(size=15),
        axis.title.y.left = element_text(size=17, color="#CB6767"), axis.text.y.left = element_text(size=15, color="#CB6767"),
        axis.title.y.right = element_text(size=17, color="#73BDD3"), axis.text.y.right = element_text(size=15, color="#73BDD3"),
        legend.text = element_text(size=13), legend.title = element_text(size=13))


# ************************* #
# facet_wrap of distance and R2 values

bdm1 <- reshape::melt.list(bd.ageBins.dists)
colnames(bdm1) <- c("value","Age_bin")
bdm1$valType <- "(B) Distance to spatial median"

bdm2 <- reshape::melt.list(asabi.all_R2)
colnames(bdm2) <- c("value","Age_bin")
bdm2 <- bdm2[ ! is.na(bdm2$value), ]
bdm2$valType <- "(A) Permanova R<sup>2</sup>"

bdab.melt <- rbind(bdm2, bdm1)


ggplot(bdab.melt, aes(x=Age_bin, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free") +
  scale_fill_manual(values=c("#73BDD3","#CB6767"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17), axis.text.y = element_text(size=15),
        strip.text = element_text(size=13), strip.background = element_rect(fill="white"))



# bdab.melt$Age_bin2 <- factor(bdab.melt$Age_bin, levels = rev(unique(bdab.melt$Age_bin)))

# https://github.com/wilkelab/ggtext/issues/48
bdab.melt$Age_bin2 <- factor(gsub("60\\+","&gt;60", gsub("_","-",bdab.melt$Age_bin)),
                              levels=c("13-20","20-30","30-40","40-50","50-60","&gt;60"))
bdab.melt$bin_and_type <- paste(bdab.melt$Age_bin, bdab.melt$valType, sep = "-")

# ********** #
bin_with_info <- sapply(unique(bdab.melt$bin_and_type), function(x) {
  
  bin  <- strsplit(x, "-")[[1]][1]
  type <- strsplit(x, "-")[[1]][2]
  
  mTab.sub <- meta.healthy[ all_subSamps$Age_bins$`1`, ]
  bin_num  <- nrow(mTab.sub[ ! is.na(mTab.sub$Age_bins) & mTab.sub$Age_bins == bin, ])
  
  ado_p <- asabi.mean_padj[ bin ]
  ado_p.stars <- ifelse(ado_p < 0.0001, "****",
                        ifelse(ado_p < 0.001, "***",
                               ifelse(ado_p < 0.01, "**",
                                      ifelse(ado_p < 0.05, "*", ""))))
  
  alt_x <- unique(bdab.melt[ bdab.melt$Age_bin == bin, "Age_bin2"])
  
  if (type == "(A) Permanova R<sup>2</sup>")
    sprintf("%s (n=%s) %s", alt_x, bin_num, 
            sprintf("<span style = 'color:#be4141;font-size:18pt;'>%s</span>", ado_p.stars) )
  else
    sprintf("%s (n=%s)", alt_x, bin_num)
  
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
# ********** #

bdab.melt$Age_bin.alt <- sapply(bdab.melt$bin_and_type, function(x) bin_with_info[ x ])
bdab.melt$Age_bin.alt <- factor(bdab.melt$Age_bin.alt, levels = unique(bdab.melt$Age_bin.alt))
bdab.melt$Age_bin.alt2 <- factor(bdab.melt$Age_bin.alt, levels = rev(unique(bdab.melt$Age_bin.alt)))

bdmPlot <- ggplot(bdab.melt, aes(x=Age_bin.alt2, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#CB6767","#73BDD3"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text.y = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17),
        axis.text.x = element_text(size=15), 
        # axis.text.x = ggtext::element_markdown(size=15, angle = 45, vjust = 0.75), 
        strip.text = ggtext::element_markdown(size=13), strip.background = element_rect(fill="white")) +
  xlab("Age bins")

# ******************************* #
# SLL2 Paper Figure 1 ####
library(ggpubr)
# divs.age.cont <- subsampling.plots.box("Age","contVar","Age", 
#                                        c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness"), 
#                                        c(groupQs), only_cont,
#                                        gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
#                                        plotType = "scatter", plot_tukey = F, xAngle = 30,
#                                        singleSubSamp = 1, useLoess = T, facetRows = 2, use_ggtitle=F, use_xlab = T,
#                                        specify_facet_labels = c("(C) Shannon Diversity","(D) Simpson Diversity",
#                                                                 "(E) Faiths PD","(F) Species richness"))
# 
# ggarrange(bdmPlot, divs.age.cont, ncol = 2, nrow = 1, widths = c(1,1.25))#, labels = c("(a)", "(b)"))
# # ******************************* #
# 
# divs.age.cont2 <- subsampling.plots.box("Age","contVar","Age", 
#                                        c("Div.Shannon","Div.Simpson","pH","Faiths.PD","Species_Richness","BMI"), 
#                                        c(groupQs), only_cont,
#                                        gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
#                                        plotType = "scatter", plot_tukey = F, xAngle = 30,
#                                        singleSubSamp = 1, useLoess = T, facetRows = 2)
# 
# ggarrange(bdmPlot, divs.age.cont2, ncol = 2, nrow = 1, widths = c(1,1.75), labels = c("(a)", "(b)"))
# # ******************************* #
# # in order to include "Age" as the x-axis label for the scatterplots
# #   and to keep the boxes lined up well, include the label on both rows
# #   so must make separate 2-faceted plots:
# divs.age.cont.sha.sim <- subsampling.plots.box("Age","contVar","Age", 
#                                                c("Div.Shannon","Div.Simpson"), 
#                                                c(groupQs), only_cont,
#                                                gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
#                                                plotType = "scatter", plot_tukey = F, xAngle = 30,
#                                                singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
#                                                specify_facet_labels = c("(C) Shannon Diversity","(D) Simpson Diversity"))
# 
# divs.age.cont.fpd.sr <- subsampling.plots.box("Age","contVar","Age", 
#                                               c("Faiths.PD","Species_Richness"), 
#                                               c(groupQs), only_cont,
#                                               gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
#                                               plotType = "scatter", plot_tukey = F, xAngle = 30,
#                                               singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
#                                               specify_facet_labels = c("(E) Faiths PD","(F) Species richness"))
# 
# ggarrange(bdmPlot, 
#           ggarrange(divs.age.cont.sha.sim, divs.age.cont.fpd.sr, ncol = 1, nrow = 2, widths = c(1,1.25)), 
#           ncol = 2, nrow = 1, widths = c(1,1.25))#, labels = c("(a)", "(b)"))
# ******************************* #
# ******************************* #
# now instead make them 4 separate plots:
divs.age.cont.sha <- subsampling.plots.box("Age","contVar","Age", 
                                               c("Div.Shannon"), 
                                               c(groupQs), only_cont,
                                               gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                                               plotType = "scatter", plot_tukey = F, xAngle = 30, forceFacet=T, choose_yl = "Shannon Diversity",
                                               singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
                                               specify_facet_labels = c("(C) Shannon Diversity"))

divs.age.cont.sim <- subsampling.plots.box("Age","contVar","Age", 
                                               c("Div.Simpson"), 
                                               c(groupQs), only_cont,
                                               gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                                               plotType = "scatter", plot_tukey = F, xAngle = 30, forceFacet=T, choose_yl = "Simpson Diversity",
                                               singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
                                               specify_facet_labels = c("(D) Simpson Diversity"))

divs.age.cont.fpd <- subsampling.plots.box("Age","contVar","Age", 
                                              c("Faiths.PD"), 
                                              c(groupQs), only_cont,
                                              gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                                              plotType = "scatter", plot_tukey = F, xAngle = 30, forceFacet=T, choose_yl = "Faiths PD",
                                              singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
                                              specify_facet_labels = c("(E) Faiths PD"))

divs.age.cont.sr <- subsampling.plots.box("Age","contVar","Age", 
                                              c("Species_Richness"), 
                                              c(groupQs), only_cont,
                                              gloms_clr, meta.healthy, phy.healthy, age_quadratic.SubsTests, dstStruc = "Age_groups",
                                              plotType = "scatter", plot_tukey = F, xAngle = 30, forceFacet=T, choose_yl = "Species richness",
                                              singleSubSamp = 1, useLoess = T, facetRows = 1, use_ggtitle=F, use_xlab = T,
                                              specify_facet_labels = c("(F) Species richness"))

ggarrange(bdmPlot, 
          ggarrange(divs.age.cont.sha, divs.age.cont.sim, divs.age.cont.fpd, divs.age.cont.sr, 
                    ncol = 2, nrow = 2), 
          ncol = 2, nrow = 1, widths = c(1,1.25))#, labels = c("(a)", "(b)"))
# ******************************* #

### Combined Figure 1 and Figure S1:

# plot the differing genera with boxes of 13-60 and >60
meta.healthy$Age.Senior_other <- sapply(meta.healthy$Age, function(x) ifelse(x > 60, ">60", "13-60"))
# meta.healthy$Age.Senior_other <- factor(meta.healthy$Age.Senior_other, levels = c("13-60",">60"))

age_diff_genera <- c("Anaeroglobus","Eikenella","Fretibacterium","Comamonas","Olsenella","Phocaeicola",
                     "Alloprevotella","Streptobacillus","Haemophilus","Prevotella","Granulicatella","Bergeyella")


ggarrange(
  ggarrange(bdmPlot, 
          ggarrange(divs.age.cont.sha, divs.age.cont.sim, divs.age.cont.fpd, divs.age.cont.sr, 
                    ncol = 2, nrow = 2), 
          ncol = 2, nrow = 1, widths = c(1,1.25)), #, labels = c("(a)", "(b)"))
  subsampling.plots.box("Age.Senior_other","Genus","Age.Senior_other", age_diff_genera, c(groupQs,"Age.Senior_other"), only_cont,
                        gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                        plotType = "box", plot_tukey = F, include_legend_below=T, sts=15, pts = 17,
                        facetScales = "fixed", flipCoords = T, facetRows = 6, str.bck = "white",
                        choose_xl = "Senior vs others", choose_yl="Abundance (centered log ratio)",#xAngle = 30, 
                        choose_title = "(G) Age vs indicated Genera",
                        singleSubSamp = 1, ignore_pvals = T, adjustAlphas = T),
  ncol = 2, nrow = 1, widths = c(1.85, 1)
)
# ******************************* #





# ******************************* #
# For supplementary figure of Permanova and Distance to spatial median based on Unifrac (weighted and unweighted):

mTab.ageBins <- meta.healthy[ unique(unlist(lapply(all_subSamps$Age_bins, function(x) x))), ]
mTab.ageBins$`13_20` <- ifelse(mTab.ageBins$Age_bins == "13_20", "Yes", "No")
mTab.ageBins$`20_30` <- ifelse(mTab.ageBins$Age_bins == "20_30", "Yes", "No")
mTab.ageBins$`30_40` <- ifelse(mTab.ageBins$Age_bins == "30_40", "Yes", "No")
mTab.ageBins$`40_50` <- ifelse(mTab.ageBins$Age_bins == "40_50", "Yes", "No")
mTab.ageBins$`50_60` <- ifelse(mTab.ageBins$Age_bins == "50_60", "Yes", "No")
mTab.ageBins$`60+`   <- ifelse(mTab.ageBins$Age_bins == "60+", "Yes", "No")

phy.ageBins <- phyloseq(otu_table(prune_samples(rownames(mTab.ageBins), SLL2)),
                        sample_data(mTab.ageBins),
                        phy_tree(prune_samples(rownames(mTab.ageBins), SLL2)),
                        tax_table(prune_samples(rownames(mTab.ageBins), SLL2)))

# **************************** #####
# Supp figure Age - Weighted_Unifrac ####

ado.sub.wu.pmeans <- mean(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Weighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))
ado.sub.wu.pmedians <- median(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Weighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))
ado.sub.wu.psds <- sd(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Weighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))

ado.sub.wu.R2means <- mean(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Weighted_Unifrac["Age_bins","R2"])))
ado.sub.wu.Fmeans <- mean(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Weighted_Unifrac["Age_bins","F.Model"])))



# aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))
wu <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))


bd.ageBins.wu <- lapply(all_subSamps$Age_bins, function(x)
  betadisper(as.dist(wu[x,x]), mTab.ageBins[x, "Age_bins"]) )

bd.ageBins.wu.pmean <- mean(p.adjust(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.wu.pmedian <- median(p.adjust(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.wu.psd <- sd(p.adjust(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))

bd.ageBins.wu.Fmean <- mean(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","F value"])))
bd.ageBins.wu.Fmedian <- median(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","F value"])))
bd.ageBins.wu.Fsd <- sd(unlist(lapply(bd.ageBins.wu, function(x) anova(x)["Groups","F value"])))

bd.ageBins.wu.dists <- list("13_20" = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="13_20"])),
                            "20_30" = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="20_30"])),
                            "30_40" = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="30_40"])),
                            "40_50" = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="40_50"])),
                            "50_60" = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="50_60"])),
                            "60+"   = unlist(lapply(bd.ageBins.wu, function(x) x$distances[x$group=="60+"])))

boxplot(bd.ageBins.wu.dists, notch=T)


# **************************************** #
# get adonis values for each individual age bin

adonis_subSamps.ageBins_indiv <- readRDS(sprintf("%s/R_objects/adonis_subSamps.ageBins_indiv.rds", p2_dir))

asabi.wu.mean_padj <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(p.adjust(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Weighted_Unifrac[x, "Pr(>F)"])), method = "fdr")))

asabi.wu.mean_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Weighted_Unifrac[x, "R2"]))))

asabi.wu.mean_F <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Weighted_Unifrac[x, "F.Model"]))))


asabi.wu.all_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) {
    nSamps <- length(i$samples[ meta.healthy[i$samples, "Age_bins"] == x ])
    c( i$Adonis$Weighted_Unifrac[x, "R2"], rep(NA, nSamps-1) )
  }))
)



# ********************* #
# ********************* #

# get betadisper results for each individual age bin against others

# betadisp.ageBins_indiv.wu <- list()
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ab)
#   betadisp.ageBins_indiv.wu[[ ab ]] <- lapply(all_subSamps$Age_bins, function(x) {
#     mtab.bd_bins <- meta.healthy[ x, ]
#     mtab.bd_bins$age_bins.binary <- sapply(mtab.bd_bins$Age_bins, function(y) ifelse(y==ab, "Bin", "Other"))
#     betadisper(as.dist(wu[x,x]), mtab.bd_bins[x, "age_bins.binary"])
#   })
# }
# saveRDS(betadisp.ageBins_indiv.wu, sprintf("%s/R_objects/betadisp.ageBins_indiv.wu.rds", p2_dir))

betadisp.ageBins_indiv.wu <- readRDS(sprintf("%s/R_objects/betadisp.ageBins_indiv.wu.rds", p2_dir))

bd.ageBins_indiv.wu.pmean <- c(
  "13_20" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`13_20`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "20_30" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`20_30`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "30_40" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`30_40`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "40_50" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`40_50`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "50_60" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`50_60`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "60+"   = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.wu$`60+`,   function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )
bd.ageBins_indiv.wu.Fmean <- c(
  "13_20" = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`13_20`, function(x) anova(x)["Groups","F value"]))),
  "20_30" = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`20_30`, function(x) anova(x)["Groups","F value"]))),
  "30_40" = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`30_40`, function(x) anova(x)["Groups","F value"]))),
  "40_50" = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`40_50`, function(x) anova(x)["Groups","F value"]))),
  "50_60" = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`50_60`, function(x) anova(x)["Groups","F value"]))),
  "60+"   = mean(unlist(lapply(betadisp.ageBins_indiv.wu$`60+`,   function(x) anova(x)["Groups","F value"]))) )

# ********************* #
# ********************* #



# ************************* #
# facet_wrap of distance and R2 values

bdm1.wu <- reshape::melt.list(bd.ageBins.wu.dists)
colnames(bdm1.wu) <- c("value","Age_bin")
bdm1.wu$valType <- "<span style = 'color:darkgoldenrod;'>(C) Distance to spatial median<br>Weighted Unifrac</span>"

bdm2.wu <- reshape::melt.list(asabi.wu.all_R2)
colnames(bdm2.wu) <- c("value","Age_bin")
bdm2.wu <- bdm2.wu[ ! is.na(bdm2.wu$value), ]
bdm2.wu$valType <- "<span style = 'color:darkgoldenrod;'>(A) Permanova R<sup>2</sup><br>Weighted Unifrac</span>"

bdab.melt.wu <- rbind(bdm2.wu, bdm1.wu)


# https://github.com/wilkelab/ggtext/issues/48
bdab.melt.wu$Age_bin2 <- factor(gsub("60\\+","&gt;60", gsub("_","-",bdab.melt.wu$Age_bin)),
                                levels=c("13-20","20-30","30-40","40-50","50-60","&gt;60"))
bdab.melt.wu$bin_and_type <- paste(bdab.melt.wu$Age_bin, bdab.melt.wu$valType, sep = "-")

# ********** #
bin_with_info.wu <- sapply(unique(bdab.melt.wu$bin_and_type), function(x) {
  
  bin  <- strsplit(x, "-")[[1]][1]
  type <- strsplit(x, "-")[[1]][2]
  
  mTab.sub <- meta.healthy[ all_subSamps$Age_bins$`1`, ]
  bin_num  <- nrow(mTab.sub[ ! is.na(mTab.sub$Age_bins) & mTab.sub$Age_bins == bin, ])
  
  ado_p <- asabi.wu.mean_padj[ bin ]
  ado_p.stars <- ifelse(ado_p < 0.0001, "****",
                        ifelse(ado_p < 0.001, "***",
                               ifelse(ado_p < 0.01, "**",
                                      ifelse(ado_p < 0.05, "*", ""))))
  
  alt_x <- unique(bdab.melt.wu[ bdab.melt.wu$Age_bin == bin, "Age_bin2"])
  
  if (type == "<span style = 'color:darkgoldenrod;'>(A) Permanova R<sup>2</sup><br>Weighted Unifrac</span>")
    sprintf("%s (n=%s) %s", alt_x, bin_num, 
            sprintf("<span style = 'color:#be4141;font-size:18pt;'>%s</span>", ado_p.stars) )
  else
    sprintf("%s (n=%s)", alt_x, bin_num)
  
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
# ********** #

bdab.melt.wu$Age_bin.alt <- sapply(bdab.melt.wu$bin_and_type, function(x) bin_with_info.wu[ x ])
bdab.melt.wu$Age_bin.alt <- factor(bdab.melt.wu$Age_bin.alt, levels = unique(bdab.melt.wu$Age_bin.alt))
bdab.melt.wu$Age_bin.alt2 <- factor(bdab.melt.wu$Age_bin.alt, levels = rev(unique(bdab.melt.wu$Age_bin.alt)))

bdmPlot.wu <- ggplot(bdab.melt.wu, aes(x=Age_bin.alt2, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#CB6767","#73BDD3"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text.y = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17),
        axis.text.x = element_text(size=15), 
        # axis.text.x = ggtext::element_markdown(size=15, angle = 45, vjust = 0.75), 
        strip.text = ggtext::element_markdown(size=15), strip.background = element_rect(fill="white")) +
  xlab("Age bins")

bdmPlot.wu









# **************************** #####
# Supp figure Age - Unweighted_Unifrac ####

ado.sub.uu.pmeans <- mean(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Unweighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))
ado.sub.uu.pmedians <- median(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Unweighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))
ado.sub.uu.psds <- sd(p.adjust(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Unweighted_Unifrac["Age_bins","Pr(>F)"])), method = "fdr"))

ado.sub.uu.R2means <- mean(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Unweighted_Unifrac["Age_bins","R2"])))
ado.sub.uu.Fmeans <- mean(unlist(lapply(adonis_subSamps$Age_bins, function(x) x$Adonis$Unweighted_Unifrac["Age_bins","F.Model"])))



uu <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))

bd.ageBins.uu <- lapply(all_subSamps$Age_bins, function(x)
  betadisper(as.dist(uu[x,x]), mTab.ageBins[x, "Age_bins"]) )

bd.ageBins.uu.pmean <- mean(p.adjust(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.uu.pmedian <- median(p.adjust(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))
bd.ageBins.uu.psd <- sd(p.adjust(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr"))

bd.ageBins.uu.Fmean <- mean(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","F value"])))
bd.ageBins.uu.Fmedian <- median(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","F value"])))
bd.ageBins.uu.Fsd <- sd(unlist(lapply(bd.ageBins.uu, function(x) anova(x)["Groups","F value"])))

bd.ageBins.uu.dists <- list("13_20" = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="13_20"])),
                            "20_30" = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="20_30"])),
                            "30_40" = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="30_40"])),
                            "40_50" = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="40_50"])),
                            "50_60" = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="50_60"])),
                            "60+"   = unlist(lapply(bd.ageBins.uu, function(x) x$distances[x$group=="60+"])))

boxplot(bd.ageBins.uu.dists, notch=T)







# aSubs.uu <- list()
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ab)
#   aSubs.uu[[ ab ]] <- run_subsampling_adonis(ab, phy.ageBins, mTab.ageBins, 100, chosenControls = all_subSamps$Age_bins, which_dists = "Unweighted_Unifrac")
# }
# 
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   for (n in names(adonis_subSamps.ageBins_indiv[[ ab ]])) {
#     adonis_subSamps.ageBins_indiv[[ ab ]][[ n ]]$Adonis$Unweighted_Unifrac <- aSubs.uu[[ ab ]][[ n ]]$Adonis$Unweighted_Unifrac
#   }
# }

# **************************************** #
# get adonis values for each individual age bin

adonis_subSamps.ageBins_indiv <- readRDS(sprintf("%s/R_objects/adonis_subSamps.ageBins_indiv.rds", p2_dir))

asabi.uu.mean_padj <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(p.adjust(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Unweighted_Unifrac[x, "Pr(>F)"])), method = "fdr")))

asabi.uu.mean_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Unweighted_Unifrac[x, "R2"]))))

asabi.uu.mean_F <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  mean(unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) i$Adonis$Unweighted_Unifrac[x, "F.Model"]))))


asabi.uu.all_R2 <- sapply(names(adonis_subSamps.ageBins_indiv), function(x) 
  unlist(lapply(adonis_subSamps.ageBins_indiv[[ x ]], function(i) {
    nSamps <- length(i$samples[ meta.healthy[i$samples, "Age_bins"] == x ])
    c( i$Adonis$Unweighted_Unifrac[x, "R2"], rep(NA, nSamps-1) )
  }))
)



# ********************* #
# ********************* #

# get betadisper results for each individual age bin against others

# betadisp.ageBins_indiv.uu <- list()
# for (ab in c("13_20","20_30","30_40","40_50","50_60","60+")) {
#   print(ab)
#   betadisp.ageBins_indiv.uu[[ ab ]] <- lapply(all_subSamps$Age_bins, function(x) {
#     mtab.bd_bins <- meta.healthy[ x, ]
#     mtab.bd_bins$age_bins.binary <- sapply(mtab.bd_bins$Age_bins, function(y) ifelse(y==ab, "Bin", "Other"))
#     betadisper(as.dist(uu[x,x]), mtab.bd_bins[x, "age_bins.binary"])
#   })
# }
# saveRDS(betadisp.ageBins_indiv.uu, sprintf("%s/R_objects/betadisp.ageBins_indiv.uu.rds", p2_dir))

betadisp.ageBins_indiv.uu <- readRDS(sprintf("%s/R_objects/betadisp.ageBins_indiv.uu.rds", p2_dir))

bd.ageBins_indiv.uu.pmean <- c(
  "13_20" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`13_20`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "20_30" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`20_30`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "30_40" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`30_40`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "40_50" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`40_50`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "50_60" = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`50_60`, function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")),
  "60+"   = mean(p.adjust(unlist(lapply(betadisp.ageBins_indiv.uu$`60+`,   function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )
bd.ageBins_indiv.uu.Fmean <- c(
  "13_20" = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`13_20`, function(x) anova(x)["Groups","F value"]))),
  "20_30" = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`20_30`, function(x) anova(x)["Groups","F value"]))),
  "30_40" = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`30_40`, function(x) anova(x)["Groups","F value"]))),
  "40_50" = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`40_50`, function(x) anova(x)["Groups","F value"]))),
  "50_60" = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`50_60`, function(x) anova(x)["Groups","F value"]))),
  "60+"   = mean(unlist(lapply(betadisp.ageBins_indiv.uu$`60+`,   function(x) anova(x)["Groups","F value"]))) )

# ********************* #
# ********************* #



# ************************* #
# facet_wrap of distance and R2 values

bdm1.uu <- reshape::melt.list(bd.ageBins.uu.dists)
colnames(bdm1.uu) <- c("value","Age_bin")
bdm1.uu$valType <- "<span style = 'color:darkgreen;'>(D) Distance to spatial median<br>Unweighted Unifrac</span>"

bdm2.uu <- reshape::melt.list(asabi.uu.all_R2)
colnames(bdm2.uu) <- c("value","Age_bin")
bdm2.uu <- bdm2.uu[ ! is.na(bdm2.uu$value), ]
bdm2.uu$valType <- "<span style = 'color:darkgreen;'>(B) Permanova R<sup>2</sup><br>Unweighted Unifrac</span>"

bdab.melt.uu <- rbind(bdm2.uu, bdm1.uu)


# https://github.com/wilkelab/ggtext/issues/48
bdab.melt.uu$Age_bin2 <- factor(gsub("60\\+","&gt;60", gsub("_","-",bdab.melt.uu$Age_bin)),
                                levels=c("13-20","20-30","30-40","40-50","50-60","&gt;60"))
bdab.melt.uu$bin_and_type <- paste(bdab.melt.uu$Age_bin, bdab.melt.uu$valType, sep = "-")

# ********** #
bin_with_info.uu <- sapply(unique(bdab.melt.uu$bin_and_type), function(x) {
  
  bin  <- strsplit(x, "-")[[1]][1]
  type <- strsplit(x, "-")[[1]][2]
  
  mTab.sub <- meta.healthy[ all_subSamps$Age_bins$`1`, ]
  bin_num  <- nrow(mTab.sub[ ! is.na(mTab.sub$Age_bins) & mTab.sub$Age_bins == bin, ])
  
  ado_p <- asabi.uu.mean_padj[ bin ]
  ado_p.stars <- ifelse(ado_p < 0.0001, "****",
                        ifelse(ado_p < 0.001, "***",
                               ifelse(ado_p < 0.01, "**",
                                      ifelse(ado_p < 0.05, "*", ""))))
  
  alt_x <- unique(bdab.melt.uu[ bdab.melt.uu$Age_bin == bin, "Age_bin2"])
  
  if (type == "<span style = 'color:darkgreen;'>(B) Permanova R<sup>2</sup><br>Unweighted Unifrac</span>")
    sprintf("%s (n=%s) %s", alt_x, bin_num, 
            sprintf("<span style = 'color:#be4141;font-size:18pt;'>%s</span>", ado_p.stars) )
  else
    sprintf("%s (n=%s)", alt_x, bin_num)
  
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
# ********** #

bdab.melt.uu$Age_bin.alt <- sapply(bdab.melt.uu$bin_and_type, function(x) bin_with_info.uu[ x ])
bdab.melt.uu$Age_bin.alt <- factor(bdab.melt.uu$Age_bin.alt, levels = unique(bdab.melt.uu$Age_bin.alt))
bdab.melt.uu$Age_bin.alt2 <- factor(bdab.melt.uu$Age_bin.alt, levels = rev(unique(bdab.melt.uu$Age_bin.alt)))

bdmPlot.uu <- ggplot(bdab.melt.uu, aes(x=Age_bin.alt2, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#CB6767","#73BDD3"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text.y = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17),
        axis.text.x = element_text(size=15), 
        # axis.text.x = ggtext::element_markdown(size=15, angle = 45, vjust = 0.75), 
        strip.text = ggtext::element_markdown(size=15), strip.background = element_rect(fill="white")) +
  xlab("")

bdmPlot.uu



# ******************************* #
# ******************************* #

library(ggpubr)
ggarrange(bdmPlot.wu, bdmPlot.uu, 
          ncol = 2, nrow = 1)#, #widths = c(1,1.25))#, labels = c("Weighted Unifrac", "Unweighted Unifrac") )

# ******************************* #
# ******************************* #












# ******************************* #
# ******************************* #

bdm3 <- mTab.ageBins[ , c("Div.Shannon","Age_bins")]
colnames(bdm3) <- c("value","Age_bin")
bdm3$valType <- "Shannon diversity"

bdm4 <- mTab.ageBins[ , c("Div.Simpson","Age_bins")]
colnames(bdm4) <- c("value","Age_bin")
bdm4$valType <- "Simpson diversity"

bdab.melt2 <- rbind(bdm1, bdm2, bdm3, bdm4)
# reorder levels for valType
bdab.melt2$valType <- factor(bdab.melt2$valType, levels = c("Distance to spatial median","Shannon diversity",
                                                            "permanova R2","Simpson diversity"))

bdab.melt2$Age_bin2 <- factor(bdab.melt2$Age_bin, levels = rev(unique(bdab.melt2$Age_bin)))

ggplot(bdab.melt2, aes(x=Age_bin2, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#73BDD3","#2E8B57","#E4BA4E","#CB6767"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        # axis.title.x = element_text(size=17), 
        axis.text.x = element_text(size=15), 
        strip.text = element_text(size=13), strip.background = element_rect(fill="white"))



bdab.melt3 <- bdab.melt2
bdab.melt3$valType <- factor(bdab.melt3$valType, levels = c("Distance to spatial median","permanova R2",
                                                            "Shannon diversity","Simpson diversity"))

ggplot(bdab.melt3, aes(x=Age_bin, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free") +
  scale_fill_manual(values=c("#73BDD3","#E4BA4E","#2E8B57","#CB6767"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17), axis.text.y = element_text(size=15),
        strip.text = element_text(size=13), strip.background = element_rect(fill="white"))




bdm5 <- mTab.ageBins[ , c("Faiths.PD","Age_bins")]
colnames(bdm5) <- c("value","Age_bin")
bdm5$valType <- "Faith's PD"

bdm6 <- mTab.ageBins[ , c("Species_Richness","Age_bins")]
colnames(bdm6) <- c("value","Age_bin")
bdm6$valType <- "Species richness"


bdab.melt4 <- rbind(bdm1, bdm2, bdm3, bdm4, bdm5, bdm6)
bdab.melt4$valType <- factor(bdab.melt4$valType, levels = c("Distance to spatial median","Shannon diversity","Faith's PD",
                                                            "permanova R2","Simpson diversity","Species richness"))
# https://github.com/wilkelab/ggtext/issues/48
bdab.melt4$Age_bin2 <- factor(gsub("60\\+","&gt;60", gsub("_","-",bdab.melt4$Age_bin)),
                             levels=c("13-20","20-30","30-40","40-50","50-60","&gt;60"))
bdab.melt4$bin_and_type <- paste(bdab.melt4$Age_bin, bdab.melt4$valType, sep = "-")

# ********** #
bin_with_info2 <- sapply(unique(bdab.melt4$bin_and_type), function(x) {
  
  bin  <- strsplit(x, "-")[[1]][1]
  type <- strsplit(x, "-")[[1]][2]
  
  mTab.sub <- meta.healthy[ all_subSamps$Age_bins$`1`, ]
  bin_num  <- nrow(mTab.sub[ ! is.na(mTab.sub$Age_bins) & mTab.sub$Age_bins == bin, ])
  
  ado_p <- asabi.mean_padj[ bin ]
  ado_p.stars <- ifelse(ado_p < 0.0001, "****",
                        ifelse(ado_p < 0.001, "***",
                               ifelse(ado_p < 0.01, "**",
                                      ifelse(ado_p < 0.05, "*", ""))))
  
  alt_x <- unique(bdab.melt4[ bdab.melt4$Age_bin == bin, "Age_bin2"])
  
  if (type == "permanova R2")
    sprintf("%s (n=%s) %s", alt_x, bin_num, 
            sprintf("<span style = 'color:#be4141;font-size:18pt;'>%s</span>", ado_p.stars) )
  else
    sprintf("%s (n=%s)", alt_x, bin_num)
  
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
# ********** #

bdab.melt4$Age_bin.alt <- sapply(bdab.melt4$bin_and_type, function(x) bin_with_info2[ x ])
bdab.melt4$Age_bin.alt <- factor(bdab.melt4$Age_bin.alt, levels = unique(bdab.melt4$Age_bin.alt))

# ggplot(bdab.melt4, aes(x=Age_bin, y=value, fill=valType)) +
ggplot(bdab.melt4, aes(x=Age_bin.alt, y=value, fill=valType)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  scale_fill_manual(values=c("#008c8c","#E4BA4E","#9933FF","#CB6767","#73BDD3","#3bb16f"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.x = ggtext::element_markdown(size=15, angle = 90),
        # axis.title.y = element_text(size=17),
        axis.text.y = element_text(size=15),
        strip.text = element_text(size=15), strip.background = element_rect(fill="white"))



bdab.melt4$valType2 <- factor(bdab.melt4$valType, levels = c("Distance to spatial median","permanova R2",
                                                             "Shannon diversity","Simpson diversity",
                                                             "Faith's PD","Species richness"))
# bdab.melt4$Age_bin3 <- factor(bdab.melt4$Age_bin, levels = rev(unique(bdab.melt4$Age_bin)))
bdab.melt4$Age_bin.alt2 <- factor(bdab.melt4$Age_bin.alt, levels = rev(unique(bdab.melt4$Age_bin.alt)))

ggplot(bdab.melt4, aes(x=Age_bin.alt2, y=value, fill=valType2)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType2, scales="free", nrow=3) +
  coord_flip() +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E","#3bb16f","#CB6767","#9933FF","#008c8c"), guide=NULL) +
  scale_fill_manual(values=c("#008c8c","#CB6767","#E4BA4E","#73BDD3","#9933FF","#3bb16f"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        # axis.title.x = element_text(size=17), 
        axis.text.x = element_text(size=15), plot.margin = margin(6,10,6,6),
        strip.text = element_text(size=17), strip.background = element_rect(fill="white"))

# ****************************************************************************************************************** #






# ****************************************************************************************************************** ####
# all_subSamps object ####

# # save list with just the subSample groups for each of the 100 iterations for each variable:
# anosim_subSamps <- readRDS(sprintf("%s/R_objects/anosim_subSamps.rds", p2_dir))
# all_subSamps <- sapply(names(anosim_subSamps), function(v) 
#   lapply(anosim_subSamps[[ v ]], function(n) n$samples))
# 
# ageGroupSubsTests <- readRDS(sprintf("%s/R_objects/ageGroupSubsTests.rds", p2_dir))
# all_subSamps[[ "Age_groups" ]] <- lapply(ageGroupSubsTests, function(x) x$samples)
# 
# for (binVar in c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner")) {
#   all_subSamps[[ binVar ]] <- lapply(anosim_subSamps[[ binVar ]], function(x) x$samples)
# }

# saveRDS(all_subSamps, file = sprintf("%s/R_objects/all_subSamps.rds", p2_dir))

all_subSamps <- readRDS(file = sprintf("%s/R_objects/all_subSamps.rds", p2_dir))




# ****************************************************************************************************************** ####
# Compare distances within groups for given variables ####

# ynk <- list()
# for (dist_meas in names(full_ords)) {
#   print(dist_meas)
#   ynk[[ dist_meas ]] <- get_YesVsNo_Kruskals(all_subSamps, SLL2.meta, full_ords[[ dist_meas ]], names(all_subSamps))
# }
# saveRDS(ynk, file = sprintf("%s/R_objects/yes_no_kruskals.rds", p2_dir))

# for (dist_meas in names(full_ords)) {
#   print(dist_meas)
#   ynk.2[[ dist_meas ]] <- get_YesVsNo_Kruskals(all_subSamps, SLL2.meta, full_ords[[ dist_meas ]], c("Age_groups.TAS","Age_bins"))
# }

# ynk.2 <- list()
# for (dist_meas in names(full_ords)) {
#   print(dist_meas)
#   ynk.2[[ dist_meas ]] <- get_YesVsNo_Kruskals(all_subSamps, SLL2.meta, full_ords[[ dist_meas ]], c("Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats"))
#   
#   for (var in c("Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats")) {
#     ynk[[ dist_meas ]][[ var ]] <- ynk.2[[ dist_meas ]][[ var ]]
#   }
# }


ynk <- readRDS(file = sprintf("%s/R_objects/yes_no_kruskals.rds", p2_dir))


plot_YesVsNo_Kruskals(ynk$Aitchison, "Aitchson")
plot_YesVsNo_Kruskals(ynk$Aitchison, "Aitchson", plotVal = "stat")

plot_YesVsNo_Kruskals(ynk$Weighted_Unifrac, "Weighted_Unifrac")
plot_YesVsNo_Kruskals(ynk$Unweighted_Unifrac, "Unweighted_Unifrac")

plot_YesVsNo_Kruskals(ynk$Bray, "Bray")
plot_YesVsNo_Kruskals(ynk$Jaccard, "Jaccard")

round(sort(unlist(lapply(ynk$Aitchison, function(y) mean(p.adjust(unlist(lapply(y, function(x) x$p.value)), method = "fdr"))))), 4)
round(sort(unlist(lapply(ynk$Aitchison, function(y) median(p.adjust(unlist(lapply(y, function(x) x$p.value)), method = "fdr"))))), 4)
sort(unlist(lapply(ynk$Aitchison, function(y) mean(unlist(lapply(y, function(x) x$statistic))))))

# ******************** #
ignoreVars <- c("seqGroup","Age_groups","Age_groups.TAS","Age_bins","Do_you_feel_well","Lactose_intolerant","Anemia",
                "Kidney_issues","Circulatory_issues","Lung_issues","Water_type_home",
                "MALDI.Mold_detected","Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats",
                "Bite_nails","Chew_pens","Mouth_wounds.binary","Headaches","Analgesics",
                "Hypothyroidism","Diabetes","Kissing_partner","Wash_hands_after_bathroom","Vitamin_supplements")


# ******************** #
dont_ignore <- c("Smoker","Cystic_fibrosis","Downs_Syndrome","Celiac","Hypertension","Antibiotics","MALDI.Yeast_detected","Full_MALDI.Candida")

ynmd <- as.data.frame((sapply(dont_ignore, function(v)
  lapply(all_subSamps[[ v ]], function(i) 
    sapply(c("Yes/F","No/M"), function(g) 
      mean(get_dists_per_group(v, ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1]), 
                               SLL2.meta, aitch, chosenSamps = i))
    )
  )
)))

yn_mean_diffs <- as.data.frame(matrix(NA, nrow=length(dont_ignore), ncol=3))
rownames(yn_mean_diffs) <- dont_ignore
colnames(yn_mean_diffs) <- c("Yes/F","No/M", "Diff")
yn_mean_diffs$`Yes/F` <- unlist(lapply(ynmd, function(x) mean(unlist(lapply(x, function(y) y["Yes/F"])))))
yn_mean_diffs$`No/M`  <- unlist(lapply(ynmd, function(x) mean(unlist(lapply(x, function(y) y["No/M"])))))
yn_mean_diffs$Diff    <- unlist(lapply(ynmd, function(x) mean(unlist(lapply(x, function(y) y["Yes/F"] - y["No/M"])))))

yn_mean_diffs <- yn_mean_diffs[ rev(order(abs(yn_mean_diffs$Diff))), ]
yn_mean_diffs$KW.p.adj <- unlist(lapply(ynk$Aitchison, function(y) 
  mean(p.adjust(unlist(lapply(y, function(x) x$p.value)), method = "fdr"))))[ rownames(yn_mean_diffs) ]
yn_mean_diffs$KW.stat <- unlist(lapply(ynk$Aitchison, function(y) 
  mean(unlist(lapply(y, function(x) x$statistic)))))[ rownames(yn_mean_diffs) ]



# ********************* #





ynmd_d <- as.data.frame((sapply(dont_ignore, function(v)
  lapply(all_subSamps[[ v ]], function(i) 
    sapply(c("Yes/F","No/M"), function(g) 
      median(get_dists_per_group(v, ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1]), 
                                 SLL2.meta, aitch, chosenSamps = i))
    )
  )
)))

yn_median_diffs <- as.data.frame(matrix(NA, nrow=length(dont_ignore), ncol=3))
rownames(yn_median_diffs) <- dont_ignore
colnames(yn_median_diffs) <- c("Yes/F","No/M", "Diff")
yn_median_diffs$`Yes/F` <- unlist(lapply(ynmd_d, function(x) mean(unlist(lapply(x, function(y) y["Yes/F"])))))
yn_median_diffs$`No/M`  <- unlist(lapply(ynmd_d, function(x) mean(unlist(lapply(x, function(y) y["No/M"])))))
yn_median_diffs$Diff    <- unlist(lapply(ynmd_d, function(x) mean(unlist(lapply(x, function(y) y["Yes/F"] - y["No/M"])))))

yn_median_diffs <- yn_median_diffs[ rev(order(abs(yn_median_diffs$Diff))), ]
yn_median_diffs$KW.p.adj <- unlist(lapply(ynk$Aitchison, function(y) 
  median(p.adjust(unlist(lapply(y, function(x) x$p.value)), method = "fdr"))))[ rownames(yn_median_diffs) ]
yn_median_diffs$KW.stat <- unlist(lapply(ynk$Aitchison, function(y) 
  median(unlist(lapply(y, function(x) x$statistic)))))[ rownames(yn_median_diffs) ]





ynsd <- as.data.frame((sapply(dont_ignore, function(v)
  lapply(all_subSamps[[ v ]], function(i) 
    sapply(c("Yes/F","No/M"), function(g) 
      sd(get_dists_per_group(v, ifelse(v=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1]), 
                             SLL2.meta, aitch, chosenSamps = i))
    )
  )
)))

yn_sd_diffs <- as.data.frame(matrix(NA, nrow=length(dont_ignore), ncol=2))
rownames(yn_sd_diffs) <- dont_ignore
colnames(yn_sd_diffs) <- c("Yes/F","No/M")
yn_sd_diffs$`Yes/F` <- unlist(lapply(ynsd, function(x) mean(unlist(lapply(x, function(y) y["Yes/F"])))))
yn_sd_diffs$`No/M`  <- unlist(lapply(ynsd, function(x) mean(unlist(lapply(x, function(y) y["No/M"])))))
# ****************************************************************************************************************** #


ads <- list()

ads[["Child"]] <- as.numeric(as.dist(aitch[all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Child"], all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Child"]]))
ads[["Teen"]] <- as.numeric(as.dist(aitch[unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Teen"], unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Teen"]]))
ads[["Adult"]] <- as.numeric(as.dist(aitch[unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Adult"], unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Adult"]]))
ads[["Senior"]] <- as.numeric(as.dist(aitch[all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Senior"], all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Senior"]]))

ads[["ct"]] <- as.numeric(as.matrix(aitch[all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Child"], unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Teen"]]))
ads[["ca"]] <- as.numeric(as.matrix(aitch[all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Child"], unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Adult"]]))
ads[["cs"]] <- as.numeric(as.matrix(aitch[all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Child"], all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Senior"]]))

ads[["ta"]] <- as.numeric(as.matrix(aitch[unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Teen"], unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Adult"]]))
ads[["ts"]] <- as.numeric(as.matrix(aitch[unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Teen"], all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Senior"]]))

ads[["as"]] <- as.numeric(as.matrix(aitch[unique(unlist(all_subSamps$Age_groups))[ meta.healthy[unique(unlist(all_subSamps$Age_groups)),"Age_groups"]=="Adult"], all_subSamps$Age_groups$`1`[ meta.healthy[all_subSamps$Age_groups$`1`,"Age_groups"]=="Senior"]]))
boxplot(ads, notch=T)


# ****************************************************************************************************************** #



# ************* #
# ynd.mean.anosim   <- get_YesVsNo_anosim_plotVals(dont_ignore, 
#                                                  stat = "mean", test="anosim")
# ynd.median.anosim <- get_YesVsNo_anosim_plotVals(dont_ignore, 
#                                                  stat = "median", test="anosim")
# saveRDS(ynd.mean.anosim, sprintf("%s/R_objects/ynd.mean.anosim.rds", p2_dir))
# saveRDS(ynd.median.anosim, sprintf("%s/R_objects/ynd.median.anosim.rds", p2_dir))
# 
# ynd.mean.adonis   <- get_YesVsNo_anosim_plotVals(dont_ignore, 
#                                                  stat = "mean", test="adonis")
# ynd.median.adonis <- get_YesVsNo_anosim_plotVals(dont_ignore, 
#                                                  stat = "median", test="adonis")
# saveRDS(ynd.mean.adonis, sprintf("%s/R_objects/ynd.mean.adonis.rds", p2_dir))
# saveRDS(ynd.median.adonis, sprintf("%s/R_objects/ynd.median.adonis.rds", p2_dir))


ynd.mean.anosim   <- readRDS(sprintf("%s/R_objects/ynd.mean.anosim.rds", p2_dir))
ynd.median.anosim <- readRDS(sprintf("%s/R_objects/ynd.median.anosim.rds", p2_dir))

# plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "KW.p.adj")
# plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "-log(KW.p.adj)")
# plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "KW.stat")
# plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "anosim.p.adj")
plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "-log(anosim.p.adj)")
plot_YesVsNo_anosim_plotVals(ynd.mean.anosim, "mean", "anosim.R")

# plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "KW.p.adj")
# plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "-log(KW.p.adj)")
# plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "KW.stat")
# plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "anosim.p.adj")
plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "-log(anosim.p.adj)")
plot_YesVsNo_anosim_plotVals(ynd.median.anosim, "median", "anosim.R")



ynd.mean.adonis   <- readRDS(sprintf("%s/R_objects/ynd.mean.adonis.rds", p2_dir))
ynd.median.adonis <- readRDS(sprintf("%s/R_objects/ynd.median.adonis.rds", p2_dir))


plot_YesVsNo_anosim_plotVals(ynd.mean.adonis, "mean", "-log(anosim.p.adj)")
plot_YesVsNo_anosim_plotVals(ynd.mean.adonis, "mean", "anosim.R")

plot_YesVsNo_anosim_plotVals(ynd.median.adonis, "median", "-log(anosim.p.adj)")
plot_YesVsNo_anosim_plotVals(ynd.median.adonis, "median", "anosim.R")

# # find and ignore variables that are not significant anywhere
# y1 <- unique(ynd.mean[,2:8])
# y1p <- y1[ , c("KW.p.adj","anosim.p.adj")]
# rownames(y1p) <- y1$Variable
# badVars  <- names(rowSums(y1p < 0.05)[rowSums(y1p < 0.05) == 0])
# okVars   <- names(rowSums(y1p < 0.05)[rowSums(y1p < 0.05) == 1])
# goodVars <- names(rowSums(y1p < 0.05)[rowSums(y1p < 0.05) == 2])
# y1p[ badVars, ]
# y1p[ okVars, ]
# y1p[ goodVars, ]
# 
# 
# y2 <- unique(ynd.median[,2:8])
# y2p <- y2[ , c("KW.p.adj","anosim.p.adj")]
# rownames(y2p) <- y2$Variable
# badVars  <- names(rowSums(y2p < 0.05)[rowSums(y2p < 0.05) == 0])
# okVars   <- names(rowSums(y2p < 0.05)[rowSums(y2p < 0.05) == 1])
# goodVars <- names(rowSums(y2p < 0.05)[rowSums(y2p < 0.05) == 2])
# y2p[ badVars, ]
# y2p[ okVars, ]
# y2p[ goodVars, ]
# ****************************************************************************************************************** #









# ****************************************************************************************************************** ####
# Sub-sampling Stomatotypes ####


# healthySamps <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No", ])
# meta.healthy <- SLL2.meta[ healthySamps, ]
# phy.healthy  <- prune_samples( healthySamps, SLL2)


weighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_w_unifrac.rds", p2_dir))
diag(weighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored
unweighted_Unifrac <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_uw_unifrac.rds", p2_dir))
diag(unweighted_Unifrac) <- NA # change diagonal to NA because the values are already 0s since its each sample against itself, can be ignored

bray <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_bray.rds", p2_dir))
diag(bray) <- NA
jaccard <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_jaccard.rds", p2_dir))
diag(jaccard) <- NA

aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

full_ords <- list("Aitchison"=aitch, "Weighted_Unifrac"=as.matrix(weighted_Unifrac),
                  "Unweighted_Unifrac"=as.matrix(unweighted_Unifrac), "Bray"=bray,
                  "Jaccard"=jaccard)

# ********************************** #
# Stomato_subSamps <- list()
# 
# print(Sys.time())
# for (variable in names(anosim_subSamps)) {
#   
#   print(variable)
#   
#   if (variable %in% no_sub_samp) {
#     
#     Stomato_subSamps[[ variable ]] <- run_subsampling_Stomatotypes(variable, SLL2, SLL2.meta, 1, 
#                                                                    full_ords, glomTab=gloms_clr, 
#                                                                    chosenControls=all_subSamps[[ variable ]], useAntibiotics=F)
#     
#   } else {
#     
#     Stomato_subSamps[[ variable ]] <- run_subsampling_Stomatotypes(variable, SLL2, SLL2.meta, 100, 
#                                                                    full_ords, glomTab=gloms_clr, 
#                                                                    chosenControls=all_subSamps[[ variable ]], useAntibiotics=F)
#     if (variable == "Cystic_fibrosis") {
#       # also run the multinom analysis using Antibiotics as one of the covariates
#       Stomato_subSamps[[ "Cystic_fibrosis.Antibiotics" ]] <- run_subsampling_Stomatotypes(variable, SLL2, SLL2.meta, 100, 
#                                                                                           full_ords, glomTab=gloms_clr, 
#                                                                                           chosenControls=all_subSamps[[ variable ]], 
#                                                                                           useAntibiotics=T)
#     }
#     
#   }
#   
#   
#   saveRDS(Stomato_subSamps, file = sprintf("%s/R_objects/Stomato_subSamps.rds", p2_dir))
#   print(Sys.time())
# }
# saveRDS(Stomato_subSamps, file = sprintf("%s/R_objects/Stomato_subSamps.rds", p2_dir))


# print(Sys.time())
# # for (variable in c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner")) {
# for (variable in c("MALDI.Yeast_detected","MALDI.Mold_detected")) {
# 
#   print(variable)
# 
#   Stomato_subSamps[[ variable ]] <- run_subsampling_Stomatotypes(variable, SLL2, SLL2.meta, 100,
#                                                                  full_ords, glomTab=gloms_clr,
#                                                                  chosenControls=all_subSamps[[ variable ]], useAntibiotics=F)
# 
#   saveRDS(Stomato_subSamps, file = sprintf("%s/R_objects/Stomato_subSamps.rds", p2_dir))
#   print(Sys.time())
# }
# saveRDS(Stomato_subSamps, file = sprintf("%s/R_objects/Stomato_subSamps.rds", p2_dir))

# ********************************** #

Stomato_subSamps <- readRDS(sprintf("%s/R_objects/Stomato_subSamps.rds", p2_dir))

# ********************************** #

# unlist(lapply(Stomato_subSamps, function(x)
#   sum(p.adjust(unlist(lapply(x, function(y) y$Chi_Square$Aitchison$p.value)), method = "fdr")<0.05)))

stomato.num_sig.chi <- sapply(names(Stomato_subSamps$Smoker$`1`$Chi_Square), function(z)
  unlist(lapply(Stomato_subSamps, function(x)
    # sum(p.adjust(unlist(lapply(x, function(y) y$Chi_Square[[ z ]]$p.value)), method = "fdr")<0.05) / length(x) * 100)))
    sum(p.adjust(unlist(lapply(x, function(y) y$Chi_Square[[ z ]]$p.value)), method = "fdr")<0.05) )))
stomato.num_sig.chi <- stomato.num_sig.chi[ rowSums(stomato.num_sig.chi) > 10, ]


stomato.Pmeans.chi <- sapply(names(Stomato_subSamps$Smoker$`1`$Chi_Square), function(z)
  unlist(lapply(Stomato_subSamps, function(x)
    mean(p.adjust(unlist(lapply(x, function(y) y$Chi_Square[[ z ]]$p.value)), method = "fdr")))))
stomato.Pmeans.chi <- stomato.Pmeans.chi[ rownames(stomato.num_sig.chi), ]
stomato.Pmeans.chi[ stomato.Pmeans.chi > 0.1 ] <- ""


stomato.Pmedians.chi <- sapply(names(Stomato_subSamps$Smoker$`1`$Chi_Square), function(z)
  unlist(lapply(Stomato_subSamps, function(x)
    median(p.adjust(unlist(lapply(x, function(y) y$Chi_Square[[ z ]]$p.value)), method = "fdr")))))
stomato.Pmedians.chi <- stomato.Pmedians.chi[ rownames(stomato.num_sig.chi), ]
stomato.Pmedians.chi[ stomato.Pmedians.chi > 0.1 ] <- ""


table(p.adjust(unlist(lapply(Stomato_subSamps$Smoker, function(x) x$Chi_Square$Aitchison$p.value)), method = "fdr")<0.05)
mean(unlist(lapply(Stomato_subSamps$Smoker, function(x) x$Chi_Square$Aitchison$statistic)))

table(unlist(lapply(Stomato_subSamps$Smoker, function(x) length(unique(x$mTab$Stomatotype_Aitchison)))))
table(unlist(lapply(Stomato_subSamps$Smoker, function(x) length(unique(x$mTab$Stomatotype_Weighted_Unifrac)))))
table(unlist(lapply(Stomato_subSamps$Smoker, function(x) length(unique(x$mTab$Stomatotype_Unweighted_Unifrac)))))
table(unlist(lapply(Stomato_subSamps$Smoker, function(x) length(unique(x$mTab$Stomatotype_Bray)))))
table(unlist(lapply(Stomato_subSamps$Smoker, function(x) length(unique(x$mTab$Stomatotype_Jaccard)))))

# ****************************************************************************************************************** #












# ****************************************************************************************************************** ####
# Sub-sampling PERMANOVAS ####

# adonis_subSamps <- list()
# 
# print(Sys.time())
# # for (variable in names(all_subSamps)[ ! names(all_subSamps) %in% c("seqGroup","Age_groups","Smoker")]) {
# for (variable in c("Full_MALDI.Candida")) {
#   print(variable)
#   # adonis_subSamps[[ variable ]] <- run_subsampling_adonis(variable, SLL2, SLL2.meta, 100,
#   #                                                         chosenControls=all_subSamps[[ variable ]],
#   #                                                         chosenOrdObj=full_ords)
#   adonis_subSamps.temp <- run_subsampling_adonis(gsub(".TAS","",variable), SLL2, SLL2.meta, 100,
#                                                  chosenControls=all_subSamps[[ variable ]],
#                                                  chosenOrdObj=full_ords)
#   adonis_subSamps[[ variable ]] <- adonis_subSamps.temp
#   # saveRDS(adonis_subSamps, file = sprintf("%s/R_objects/adonis_subSamps.rds", p2_dir))
#   saveRDS(adonis_subSamps.temp, file = sprintf("%s/R_objects/adonis_subSamps.temp.%s.rds", p2_dir, gsub(".TAS","",variable)))
#   print(Sys.time())
# }
# saveRDS(adonis_subSamps, file = sprintf("%s/R_objects/adonis_subSamps.rds", p2_dir))


                                               
adonis_subSamps <- readRDS(file = sprintf("%s/R_objects/adonis_subSamps.rds", p2_dir))

bd.ignore <- c("Age_groups","Age_groups.TAS","Age_bins")

ado.sub.pmeans <- sapply(names(adonis_subSamps)[ ! names(adonis_subSamps) %in% bd.ignore ], function(var)
  mean(p.adjust(unlist(lapply(adonis_subSamps[[ var ]], function(x) x$Adonis$Aitchison[var,"Pr(>F)"])), method = "fdr")) )

ado.sub.pmedians <- sapply(names(adonis_subSamps)[ ! names(adonis_subSamps) %in% bd.ignore ], function(var)
  median(p.adjust(unlist(lapply(adonis_subSamps[[ var ]], function(x) x$Adonis$Aitchison[var,"Pr(>F)"])), method = "fdr")) )

ado.sub.psds <- sapply(names(adonis_subSamps)[ ! names(adonis_subSamps) %in% bd.ignore ], function(var)
  sd(p.adjust(unlist(lapply(adonis_subSamps[[ var ]], function(x) x$Adonis$Aitchison[var,"Pr(>F)"])), method = "fdr")) )



ado.sub.R2means <- sapply(names(adonis_subSamps)[ ! names(adonis_subSamps) %in% bd.ignore ], function(var)
  mean(unlist(lapply(adonis_subSamps[[ var ]], function(x) x$Adonis$Aitchison[var,"R2"]))) )
ado.sub.Fmeans <- sapply(names(adonis_subSamps)[ ! names(adonis_subSamps) %in% bd.ignore ], function(var)
  mean(unlist(lapply(adonis_subSamps[[ var ]], function(x) x$Adonis$Aitchison[var,"F.Model"]))) )




bd.interest <- c("MALDI.Yeast_detected","Full_MALDI.Candida","Cystic_fibrosis","Downs_Syndrome","Smoker","Celiac",
                 "Hypertension","Antibiotics")#,"Gender","Braces.binary","Fluoride_toothpaste")

ado.sub.all_R2 <- sapply(bd.interest, function(x) 
  unlist(lapply(adonis_subSamps[[ x ]], function(i) {
    nSamps <- length(i$samples)#[ ! is.na(SLL2.meta[i$samples, "Age_bins"]) & SLL2.meta[i$samples, "Age_bins"] == x ])
    c( i$Adonis$Aitchison[x, "R2"], rep(NA, nSamps-1) )
  }))
)
# ado.sub.all_R2 <- sapply(bd.interest, function(x) unlist(lapply(adonis_subSamps[[ x ]], function(i) i$Adonis$Aitchison[x, "R2"])))

# bd.subs <- sapply(bd.interest, function(var)
#   lapply(adonis_subSamps[[ var ]], function(i)
#     betadisper(as.dist(aitch[i$samples, i$samples]), SLL2.meta[i$samples, var]) ))
# 
# bd.subs <- as.list(as.data.frame(bd.subs))
# saveRDS(bd.subs, sprintf("%s/R_objects/bd.subs.rds", p2_dir))

bd.subs <- readRDS(sprintf("%s/R_objects/bd.subs.rds", p2_dir))

bd.subs.pmean <- sapply(names(bd.subs), function(var)
  mean(p.adjust(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )
bd.subs.pmedian <- sapply(names(bd.subs), function(var)
  median(p.adjust(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )
bd.subs.psd <- sapply(names(bd.subs), function(var)
  sd(p.adjust(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","Pr(>F)"])), method = "fdr")) )

bd.subs.Fmean <- sapply(names(bd.subs), function(var)
  mean(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","F value"]))) )
bd.subs.Fmedian <- sapply(names(bd.subs), function(var)
  median(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","F value"]))) )
bd.subs.Fsd <- sapply(names(bd.subs), function(var)
  sd(unlist(lapply(bd.subs[[ var ]], function(x) anova(x)["Groups","F value"]))) )



# **************** #
bd.subs.dists <- sapply(names(bd.subs), function(var) {
  ynlist <- sapply(c("Yes/F","No/M"), function(g)
    unlist(lapply(bd.subs[[ var ]], function(x) 
      x$distances[x$group == ifelse(var=="Gender", strsplit(g, "/")[[1]][2], strsplit(g, "/")[[1]][1]) ]))
  )
  if ( class(ynlist) == "list")
    ynlist
  else
    as.list(as.data.frame(ynlist))
})
bd.subs.dists <- as.list(as.data.frame(bd.subs.dists))

boxplot(bd.subs.dists$Cystic_fibrosis, notch=T)


# **************** #
bds.melt <- reshape::melt.list(bd.subs.dists)
colnames(bds.melt) <- c("Distance_to_centroid","YesNo","Variable")

# bds.melt$YesNo <- factor(bds.melt$YesNo, levels = c("Yes/F","No/M"))

bds.melt$ado_mean_padj <- sapply(bds.melt$Variable, function(x) ado.sub.pmeans[ x ])
bds.melt$ado_mean_plog <- sapply(bds.melt$Variable, function(x) -log(ado.sub.pmeans[ x ]))
bds.melt$ado_mean_R2   <- sapply(bds.melt$Variable, function(x) ado.sub.R2means[ x ])
bds.melt$ado_mean_F    <- sapply(bds.melt$Variable, function(x) ado.sub.Fmeans[ x ])

bds.melt$bd_mean_padj <- sapply(bds.melt$Variable, function(x) bd.subs.pmean[ x ])
bds.melt$bd_mean_plog <- sapply(bds.melt$Variable, function(x) -log(bd.subs.pmean[ x ]))
bds.melt$bd_mean_F    <- sapply(bds.melt$Variable, function(x) bd.subs.Fmean[ x ])


# for ordering boxes
med_diffs <- sapply(unique(bds.melt$Variable), function(x) 
  median(bds.melt[bds.melt$Variable==x & bds.melt$YesNo=="Yes/F","Distance_to_centroid"]) -
    median(bds.melt[bds.melt$Variable==x & bds.melt$YesNo=="No/M","Distance_to_centroid"]))
bds.melt$med_diff <- sapply(bds.melt$Variable, function(x) med_diffs[ x ])

# bdab.melt$Age_bin <- sapply(bdab.melt$Age_bin, function(x) {
#   bin <- gsub("60\\+",">60",gsub("_","-",x))
#   padj <- round(asabi.mean_padj[ x ], 3)
#   sprintf("**%s**<br>p = %s", bin, padj)
# })
# bdab.melt$Age_bin <- factor(bdab.melt$Age_bin, levels = unique(bdab.melt$Age_bin))

var_with_info <- sapply(unique(bds.melt$Variable), function(x) {
  if (x %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well")) {
    mTab.sub <- meta.healthy[ all_subSamps[[ x ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "No", ])
    
  } else if (x %in% c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner",
                      "Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats")) {
    var_num <- 100
    
  } else if (x %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes","Hypertension","Cholesterol",
                      "Depression","Anxiety","Headaches","Lactose_intolerant","Gastritis","Intestinal_issues",
                      "Anemia","Sinusitis","Fibrosis_carrier","Thyroid_issue","Hypothyroidism","Cancer","Transplant",
                      "Immune_issues","Skin_issues","Lung_issues","Circulatory_issues","Kidney_issues",
                      "Tonsil_issues","Samp_collect_issues","Birth_control","Antihistamines")) {
    mTab.sub <- SLL2.meta[ all_subSamps[[ x ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "Yes", ])
    
  } else {
    mTab.sub <- meta.healthy[ all_subSamps[[ x ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , x]) & mTab.sub[ , x] == "Yes", ])
    
  }
  
  ado_p <- unique(bds.melt$ado_mean_padj[ bds.melt$Variable == x ])
  bd_p  <- unique(bds.melt$bd_mean_padj[ bds.melt$Variable == x ])
  
  # sprintf("%s (n=%s)\n%s p.adj = %s\n%s R = %s", x, var_num, stat, round(var_p, 3), stat, round(var_R, 3))
  if (ado_p > 0.1 & bd_p > 0.1)
    sprintf("**%s** (n=%s)<br>%s<br>%s", x, var_num, 
            sprintf("<span style = 'color:darkgrey;'>permanova p.adj = %s</span>",round(ado_p, 3)), 
            sprintf("<span style = 'color:darkgrey;'>betaDisper p.adj = %s</span>",round(bd_p, 3)) )
  else if (ado_p > 0.1)
    sprintf("**%s** (n=%s)<br>%s<br>betaDisper p.adj = %s", x, var_num, 
            sprintf("<span style = 'color:darkgrey;'>permanova p.adj = %s</span>",round(ado_p, 3)), 
            round(bd_p, 3))
  else if (bd_p > 0.1)
    sprintf("**%s** (n=%s)<br>permanova p.adj = %s<br>%s", x, var_num, round(ado_p, 3), 
            sprintf("<span style = 'color:darkgrey;'>betaDisper p.adj = %s</span>",round(bd_p, 3)) )
  else
    sprintf("**%s** (n=%s)<br>permanova p.adj = %s<br>betaDisper p.adj = %s", x, var_num, round(ado_p, 3), round(bd_p, 3))
  
    
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
bds.melt$Variable <- sapply(bds.melt$Variable, function(x) var_with_info[ x ])

bds.melt$YesNo <- gsub("Yes/F","Yes", gsub("No/M","No", bds.melt$YesNo))
saveRDS(bds.melt, sprintf("%s/R_objects/bds.melt.rds", p2_dir))

ggplot(bds.melt, aes(x=reorder(Variable, abs(med_diff), FUN=median), y=Distance_to_centroid, fill=ado_mean_R2)) +
  geom_boxplot(notch = T, aes(color=YesNo)) +
  coord_flip() +
  labs(y="Distance to spatial median") +
  scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  scale_color_manual(values = c("darkgreen","darkred"), name=NULL, guide=guide_legend(reverse = T)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=17), axis.text.x = element_text(size=15),
        axis.title.y = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        legend.text = element_text(size=13), legend.title = element_text(size=13))




# ****************************************************************************************************************** #

bds.melt$all_R2 <- (unlist(ado.sub.all_R2) * 2) + 0.06
# bdab.melt$all_R2    <- unlist(asabi.all_R2) + (min(bdab.melt$Distance_to_centroid) * (max(unlist(asabi.all_R2), na.rm = T) / max(bdab.melt$Distance_to_centroid)))

scaleFactor <- max(bds.melt$Distance_to_centroid) / max(bds.melt$all_R2, na.rm = T)

ggplot(bds.melt, aes(x=reorder(Variable, abs(med_diff), FUN=median))) +
  geom_boxplot(notch = T, aes(color=YesNo, y=Distance_to_centroid)) +
  # see::geom_violinhalf(aes(color=YesNo, y=Distance_to_centroid)) +
  geom_boxplot(notch = T, aes(y=all_R2 * scaleFactor), fill="blue") +
  scale_y_continuous(name="Distance to spatial median", sec.axis=sec_axis(~(./scaleFactor - 0.06)/2, name="permanova R2")) +
  coord_flip() +
  labs(y="Distance to spatial median") +
  scale_fill_gradient2(low="white", high="darkblue", name="mean R2") +
  scale_color_manual(values = c("darkgreen","darkred"), name=NULL, guide=guide_legend(reverse = T)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size=17), axis.text.x = element_text(size=15),
        axis.title.y = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        legend.text = element_text(size=13), legend.title = element_text(size=13))


bds.viohalf.yes <- sapply(unique(bds.melt$Variable), function(x) bds.melt$Distance_to_centroid[bds.melt$Variable == x & bds.melt$YesNo == "Yes"])
bds.viohalf.no  <- sapply(unique(bds.melt$Variable), function(x) bds.melt$Distance_to_centroid[bds.melt$Variable == x & bds.melt$YesNo == "No"])

vioplot::vioplot(bds.viohalf.yes, side="left", plotcenter="line", horizontal=T, add=T)
vioplot::vioplot(bds.viohalf.no, side="right", plotcenter="line", horizontal=T, add=T)


# ****************************************************************************************************************** #
bds1 <- reshape::melt.list(bd.subs.dists)
colnames(bds1) <- c("value","YesNo","Variable")
bds1$YesNo <- gsub("Yes/F","Yes", gsub("No/M","No", bds1$YesNo))
bds1$valType <- "(B) Distance to spatial median"


bds2 <- reshape::melt.list(ado.sub.all_R2)
colnames(bds2) <- c("value","Variable")
bds2 <- bds2[ ! is.na(bds2$value), ]
bds2$YesNo <- "R2"
bds2$valType <- "(A) Permanova R<sup>2</sup>"
bds2 <- bds2[ , c("value","YesNo","Variable","valType")]

bds.combo <- rbind(bds2, bds1)


# for ordering boxes
med_diffs.combo <- sapply(unique(bds.combo$Variable), function(x) 
  median(bds.combo[bds.combo$Variable==x & bds.combo$YesNo=="Yes","value"]) -
    median(bds.combo[bds.combo$Variable==x & bds.combo$YesNo=="No","value"]))
bds.combo$med_diff <- sapply(bds.combo$Variable, function(x) med_diffs.combo[ x ])


bds.combo$Variable2 <- gsub("MALDI.", "", gsub("_", " ", gsub("Full_MALDI.Candida", "Candida detected", bds.combo$Variable)))


ggplot(bds.combo, aes(x=reorder(Variable2, abs(med_diff), FUN=median), y=value, fill=YesNo)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free") +
  scale_fill_manual(values=c("#73BDD3","#CB6767","#E4BA4E"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = ggtext::element_markdown(size=15),
        axis.title.y = element_text(size=17), axis.text.y = element_text(size=15),
        strip.text = element_text(size=13), strip.background = element_rect(fill="white"))



bds.combo$Variable3 <- factor(bds.combo$Variable2, levels = rev(unique(bds.combo$Variable2)))

bds.combo$var_and_type <- paste(bds.combo$Variable, bds.combo$valType, sep = "-")

# ********** #
var_with_info2 <- sapply(unique(bds.combo$var_and_type), function(x) {
  
  var  <- strsplit(x, "-")[[1]][1]
  type <- strsplit(x, "-")[[1]][2]
  
  if (var %in% c("Fluoride_toothpaste","Wash_hands_before_eat","Wash_hands_after_bathroom","Do_you_feel_well")) {
    mTab.sub <- meta.healthy[ all_subSamps[[ var ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , var]) & mTab.sub[ , var] == "No", ])
    
  } else if (var %in% c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner",
                      "Pets","Pets.Dogs","Pets.Cats","Pets.Both_dogs_cats")) {
    var_num <- 100
    
  } else if (var %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Diabetes","Hypertension","Cholesterol",
                      "Depression","Anxiety","Headaches","Lactose_intolerant","Gastritis","Intestinal_issues",
                      "Anemia","Sinusitis","Fibrosis_carrier","Thyroid_issue","Hypothyroidism","Cancer","Transplant",
                      "Immune_issues","Skin_issues","Lung_issues","Circulatory_issues","Kidney_issues",
                      "Tonsil_issues","Samp_collect_issues","Birth_control","Antihistamines")) {
    mTab.sub <- SLL2.meta[ all_subSamps[[ var ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , var]) & mTab.sub[ , var] == "Yes", ])
    
  } else {
    mTab.sub <- meta.healthy[ all_subSamps[[ var ]]$`1`, ]
    var_num  <- nrow(mTab.sub[ ! is.na(mTab.sub[ , var]) & mTab.sub[ , var] == "Yes", ])
    
  }
  
  ado_p <- ado.sub.pmeans[ var ]
  ado_p.stars <- ifelse(ado_p < 0.0001, "****",
                        ifelse(ado_p < 0.001, "***",
                               ifelse(ado_p < 0.01, "**",
                                      ifelse(ado_p < 0.05, "*", ""))))
  bd_p  <- bd.subs.pmean[ var ]
  bd_p.stars <- ifelse(bd_p < 0.0001, "****",
                        ifelse(bd_p < 0.001, "***",
                               ifelse(bd_p < 0.01, "**",
                                      ifelse(bd_p < 0.05, "*", ""))))
  
  alt_x <- unique(bds.combo[ bds.combo$Variable == var, "Variable3"])
  
  
  # if (type == "Distance to spatial median" & bd_p > 0.05)
  #   sprintf("**%s** (n=%s)<br>%s", alt_x, var_num, 
  #           sprintf("<span style = 'color:darkgrey;'>homog P = %s</span>",round(bd_p, 3)) )
  # else if (type == "Distance to spatial median")
  #   sprintf("**%s** (n=%s)<br>homog P = %s", alt_x, var_num, round(bd_p, 3))
  # else if (type == "permanova R2" & ado_p > 0.05)
  #   sprintf("**%s** (n=%s)<br>%s", alt_x, var_num, 
  #           sprintf("<span style = 'color:darkgrey;'>perm P = %s</span>",round(ado_p, 3)) )
  # else if (type == "permanova R2")
  #   sprintf("**%s** (n=%s)<br>perm P = %s", alt_x, var_num, round(ado_p, 3))
  # else 
  #   "fsdf"
  
  if (type == "(B) Distance to spatial median")
    sprintf("**%s** (n=%s) %s", alt_x, var_num, 
            sprintf("<span style = 'color:#2a6f84;font-size:18pt;'>%s</span>", bd_p.stars) )
  else if (type == "(A) Permanova R<sup>2</sup>")
    sprintf("**%s** (n=%s) %s", alt_x, var_num, 
            sprintf("<span style = 'color:#be4141;font-size:18pt;'>%s</span>", ado_p.stars) )
  
  # check these links for properly using the markdown formatting to get bold text:
  #   https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
  #   https://github.com/wilkelab/ggtext
})
# ********** #
bds.combo$Variable4 <- sapply(bds.combo$var_and_type, function(x) var_with_info2[ x ])


ggplot(bds.combo, aes(x=reorder(Variable4, abs(med_diff), FUN=median), y=value, fill=YesNo)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#73BDD3","#CB6767","#E4BA4E"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        # axis.title.x = element_text(size=17), 
        axis.text.x = element_text(size=15), 
        strip.text = ggtext::element_markdown(size=17), strip.background = element_rect(fill="white"))

# ****************************************************************************************************************** #
fig2.a <- ggplot(bds.combo[bds.combo$valType=="(A) Permanova R<sup>2</sup>",], 
       aes(x=reorder(Variable4, abs(med_diff), FUN=median), y=value, fill=YesNo)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#CB6767"), guide=NULL) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        # axis.title.x = element_text(size=17), 
        axis.text.x = element_text(size=15), 
        strip.text = ggtext::element_markdown(size=17), strip.background = element_rect(fill="white"))

fig2.b <- ggplot(bds.combo[bds.combo$valType=="(B) Distance to spatial median",], 
                 aes(x=reorder(Variable4, abs(med_diff), FUN=median), y=value, fill=YesNo)) +
  geom_boxplot(notch = T) +
  facet_wrap(~valType, scales="free", nrow=2) +
  coord_flip() +
  scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=guide_legend(reverse=T)) +
  # scale_fill_manual(values=c("#73BDD3","#E4BA4E"), guide=NULL) +
  # labs(y="Distance to spatial median") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
        # axis.title.x = element_text(size=17), 
        legend.title = element_blank(), legend.text = element_text(size=15), 
        legend.position = "bottom", legend.box.spacing = unit(-0.5, "lines"),
        axis.text.x = element_text(size=15), 
        strip.text = ggtext::element_markdown(size=17), strip.background = element_rect(fill="white"))

ggarrange(fig2.a, fig2.b, ncol = 1, nrow = 2)#, labels = c("(a)", "(b)"))
# ****************************************************************************************************************** #





# ****************************************************************************************************************** ####
# Sub-sampling Driver Genera ####

# mean_drivers_subSamps <- sapply(names(all_subSamps), function(v) {
#   print(v)
#   rank_drivers_for_subSamps(gsub(".TAS","",v), all_subSamps[[ v ]], SLL2.meta, gloms_clr, nrow(gloms_clr$Genus), "Genus")
# })
# saveRDS(mean_drivers_subSamps, file = sprintf("%s/R_objects/mean_drivers_subSamps.rds", p2_dir))

mean_drivers_subSamps <- readRDS(file = sprintf("%s/R_objects/mean_drivers_subSamps.rds", p2_dir))


# rank top 5 genera for each group in each variable
top_drivers_subSamps <- lapply(mean_drivers_subSamps, function(v)
  sapply(colnames(v), function(g)
    sort(sapply(rownames(v), function(tax) v[tax,g]), decreasing = T)[1:10], simplify = F))

# ****************************************************************************************************************** #


# anova_subSamps <- list()
# for (variable in c("Celiac","Smoker","MALDI.Yeast_detected","Full_MALDI.Candida","Gender","Braces.binary","Antibiotics","Hypertension")) {
#   # for (variable in c("Sr","Hardness","Alcalinity")) {
# 
#   print(variable)
# 
#   mTab.anovaSubs <- meta.healthy
#   phy.anovaSubs  <- phy.healthy
# 
#   if (variable %in% c("Celiac","Hypertension")) {
#     mTab.anovaSubs <- SLL2.meta
#     phy.anovaSubs  <- SLL2
#   }
# 
#   covVars <- c("Age","Gender","Population")
#   if (variable == "Gender")
#     covVars <- c("Age","Population")
# 
#   anova_subSamps[[ variable ]] <- run_full_subsampling_calcs.ageBins("healthy", phy.anovaSubs, mTab.anovaSubs, 100,
#                                                                      trait = variable, chosenControls = all_subSamps[[ variable ]],
#                                                                      covVars = covVars, #subSampObjects = waterType.subSamp.objects,
#                                                                      dontReturnMTab = T, dontReturnSamps = T, anova_only = T)
# }
# saveRDS(anova_subSamps, file = sprintf("%s/R_objects/anova_subSamps.rds", p2_dir))

anova_subSamps <- readRDS(file = sprintf("%s/R_objects/anova_subSamps.rds", p2_dir))

# ******************************************************** #

gen.pMeans.anova_subSamps   <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "Genus", "mean"))
phy.pMeans.anova_subSamps   <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "Phylum", "mean"))
conts.pMeans.anova_subSamps <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "contVar", "mean"))

gen.num_sig.anova_subSamps   <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "Genus", "num_sig"))
phy.num_sig.anova_subSamps   <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "Phylum", "num_sig"))
conts.num_sig.anova_subSamps <- lapply(anova_subSamps, function(w) get_subs_anova_per_taxon(w, 1, "contVar", "num_sig"))

# ******************************************************** #

varCheck <- "Celiac";       plot.mTab <- SLL2.meta; plot.phy <- SLL2
varCheck <- "Hypertension"; plot.mTab <- SLL2.meta; plot.phy <- SLL2
varCheck <- "Smoker";               plot.mTab <- meta.healthy; plot.phy <- phy.healthy
varCheck <- "MALDI.Yeast_detected"; plot.mTab <- meta.healthy; plot.phy <- phy.healthy
varCheck <- "Full_MALDI.Candida";   plot.mTab <- meta.healthy; plot.phy <- phy.healthy
varCheck <- "Gender";               plot.mTab <- meta.healthy; plot.phy <- phy.healthy
varCheck <- "Braces.binary";        plot.mTab <- meta.healthy; plot.phy <- phy.healthy
varCheck <- "Antibiotics";          plot.mTab <- meta.healthy; plot.phy <- phy.healthy

plot_freq_sigs(varCheck, gen.num_sig.anova_subSamps[[ varCheck ]], gen.pMeans.anova_subSamps[[ varCheck ]], 50, "Genus", 
               plot.mTab, plot.phy, anova_subSamps[[ varCheck ]], chosenSamps = all_subSamps[[ varCheck ]]$`1`, adjustAlphas = T)

plot_freq_sigs(varCheck, phy.num_sig.anova_subSamps[[ varCheck ]], phy.pMeans.anova_subSamps[[ varCheck ]], 50, "Phylum", 
               plot.mTab, plot.phy, anova_subSamps[[ varCheck ]], chosenSamps = all_subSamps[[ varCheck ]]$`1`, adjustAlphas = T)

plot_freq_sigs(varCheck, conts.num_sig.anova_subSamps[[ varCheck ]], conts.pMeans.anova_subSamps[[ varCheck ]], 50, "contVar", 
               plot.mTab, plot.phy, anova_subSamps[[ varCheck ]], chosenSamps = all_subSamps[[ varCheck ]]$`1`, adjustAlphas = T)

# ****************************************************************************************************************** #








# ****************************************************************************************************************** ####
# ****************************************************************************************************************** #
# BMI ####

aitch <- readRDS(sprintf("%s/R_objects/beta_diversities/SLL2_aitch.rds", p2_dir))

bmi.samps <- rownames(meta.healthy)[ ! is.na(meta.healthy$BMI) ]

# bmi.tests <- list("BMI"=list("Anova"=list()), "BMI_official"=list("Anova"=list()))
# 
# bmi.tests[[ "BMI" ]][[ "Adonis" ]] <- get_adonis("BMI", c("Gender","Population"), "Aitchison", meta.healthy[ bmi.samps, ], 
#                                                  ordObj = as.dist(aitch[ bmi.samps, bmi.samps ]) )
# bmi.tests[[ "BMI_official" ]][[ "Adonis" ]] <- get_adonis("BMI_official", c("Gender","Population"), "Aitchison", meta.healthy[ bmi.samps, ], 
#                                                           ordObj = as.dist(aitch[ bmi.samps, bmi.samps ]) )
# 
# 
# for (tl in c("contVar","Phylum","Class","Order","Family","Genus","Species")) {
#   print(tl)
#   # ******************** #
#   if (tl == "contVar") {
#     dv <- c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","Gene_counts","pH")
#     # ******************** #
#   } else {
#     dv <- NULL
#   }
#   # ******************** #
#   
#   bmi.tests[[ "BMI" ]][[ "Anova" ]][[ tl ]] <- get_lm( c("BMI", "Gender","Population"), 
#                                                        tl, meta.healthy[ bmi.samps, ], 
#                                                        gloms_clr[[ tl ]][ , bmi.samps ], 
#                                                        dv, noRemove = T, rerun.nonSig = F)
#   
#   bmi.tests[[ "BMI_official" ]][[ "Anova" ]][[ tl ]] <- get_lm( c("BMI_official", "Gender","Population"), 
#                                                                 tl, meta.healthy[ bmi.samps, ], 
#                                                                 gloms_clr[[ tl ]][ , bmi.samps ], 
#                                                                 dv, noRemove = T, rerun.nonSig = F)
# }
# saveRDS(bmi.tests, file = sprintf("%s/R_objects/bmi.tests.rds", p2_dir))

bmi.tests <- readRDS(file = sprintf("%s/R_objects/bmi.tests.rds", p2_dir))

# ****************************************************************************************************************** #
# age_bin.subSamp.objects <- get_subsamp_objects(all_subSamps$Age_bins, phy.healthy, meta.healthy, 100)
# saveRDS(age_bin.subSamp.objects, file = sprintf("%s/R_objects/age_bin.subSamp.objects.rds", p2_dir))
# 
# age_bin.subSamp.objects <- readRDS(sprintf("%s/R_objects/age_bin.subSamp.objects.rds", p2_dir))
# 
# # run these tests using the Age_bins subsamples
# bmi_subsTests <- list()
# 
# for (variable in c("BMI","BMI_official","BMI_group","pH")) {
# # for (variable in c("Sr","Hardness","Alcalinity")) {
# 
#   print(variable)
# 
#   covVars <- c("Gender","Population") # BMI and pH were significantly associated with Age, so dont include it as a covVar
#   
#   bmi_subsTests[[ variable ]] <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                                     trait = variable, chosenControls = all_subSamps$Age_bins,
#                                                                     covVars = covVars, subSampObjects = age_bin.subSamp.objects,
#                                                                     dontReturnMTab = T, dontReturnSamps = T)
# }
# saveRDS(bmi_subsTests, file = sprintf("%s/R_objects/bmi_subsTests.rds", p2_dir))

bmi_subsTests <- readRDS(file = sprintf("%s/R_objects/bmi_subsTests.rds", p2_dir))
# ****************************************************************************************************************** #


gen.pMeans.bmi   <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "Genus", "mean")
phy.pMeans.bmi   <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "Phylum", "mean")
conts.pMeans.bmi <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "contVar", "mean")

gen.num_sig.bmi   <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "Genus", "num_sig")
phy.num_sig.bmi   <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "Phylum", "num_sig")
conts.num_sig.bmi <- get_subs_anova_per_taxon(bmi_subsTests$BMI, "BMI", "contVar", "num_sig")




ado.Pmeans.bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                              function(x) mean(p.adjust(unlist(lapply(bmi_subsTests$BMI, function(y) 
                                y$Adonis[[ x ]]["BMI","Pr(>F)"])), method = "fdr")))

ado.Psds.bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                            function(x) sd(p.adjust(unlist(lapply(bmi_subsTests$BMI, function(y) 
                              y$Adonis[[ x ]]["BMI","Pr(>F)"])), method = "fdr")))

ado.num_sig.bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                               function(x) sum(p.adjust(unlist(lapply(bmi_subsTests$BMI, function(y) 
                                 y$Adonis[[ x ]]["BMI","Pr(>F)"])), method = "fdr") < 0.05))



ado.R2means.bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                               function(x) mean(unlist(lapply(bmi_subsTests$BMI, function(y) 
                                 y$Adonis[[ x ]]["BMI","R2"]))))

ado.Fmeans.bmi <- sapply(c("Aitchison","Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard"),
                              function(x) mean(unlist(lapply(bmi_subsTests$BMI, function(y) 
                                y$Adonis[[ x ]]["BMI","F.Model"]))))



# ************************************************************************************** #

# only those for which at least 95/100 subsamplings were signif (ie p < 0.05)
freqSig.ab <- gen.num_sig.bmi[ gen.num_sig.bmi >= 75]
sort(gen.pMeans.bmi[ names(freqSig.ab) ])

# gensToCheck <- unique(c("Actinobacillus","Porphyromonas","Prevotella","Eikenella","Bacteroides",
#                         "Fusobacterium","Treponema","Campylobacter","Brevundimonas",
#                         "Gemella","Alloprevotella","Atopobium","Kingella", names(freqSig.ab)))

gensToCheck <- unique(names(freqSig.ab))
# as.data.frame(cbind(gen.num_sig.bmi[ gensToCheck ], formatC(gen.pMeans.bmi[ gensToCheck ],  format="f", digits=5)))e

freq.mean.ps <- as.data.frame(cbind(gen.num_sig.bmi[ gensToCheck ], gen.pMeans.bmi[ gensToCheck ]))
colnames(freq.mean.ps) <- c("num_sig","meanP")
freq.mean.ps <- freq.mean.ps[ rev(order(freq.mean.ps$num_sig, -freq.mean.ps$meanP)), ]
freq.mean.ps


# subsampling.plots.box("Age","Genus","Age", names(rev(sort(freqSig.ab))), c(groupQs), only_cont,
#                       gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
#                       plotType = "scatter", plot_tukey = F, xAngle = 30)

subsampling.plots.box("Age","Genus","Age", rownames(freq.mean.ps), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30, 
                      singleSubSamp = 1)

phy.sig.order <- phy.num_sig.bmi[ rev(order(phy.num_sig.bmi, -phy.pMeans.bmi)) ]
phy.sig.order <- phy.sig.order[ phy.pMeans.bmi[ names(phy.sig.order) ] < 0.05 ]
subsampling.plots.box("Age","Phylum","Age", names(phy.sig.order), c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups", 
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1)

subsampling.plots.box("Age","contVar","Age", 
                      c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness","pH","BMI"), 
                      c(groupQs), only_cont,
                      gloms_clr, meta.healthy, phy.healthy, age_cont.SubsTests, dstStruc = "Age_groups",
                      plotType = "scatter", plot_tukey = F, xAngle = 30,
                      singleSubSamp = 1)


# ************************************************************************************** #


# ****************************************************************************************************************** #









# ****************************************************************************************************************** ####
# Water values ####

# ****************************************************************************************************************** #

# waterType.subSamp.objects <- get_subsamp_objects(all_subSamps$Water_type_home, phy.healthy, meta.healthy, 100)
# saveRDS(waterType.subSamp.objects, file = sprintf("%s/R_objects/waterType.subSamp.objects.rds", p2_dir))
# 
# waterType.subSamp.objects <- readRDS(sprintf("%s/R_objects/waterType.subSamp.objects.rds", p2_dir))
# 
# # run these tests using the Water_type_home subsamples
# water_tests <- list()
# for (variable in c(cont_water_data)) {
#   # for (variable in c("Sr","Hardness","Alcalinity")) {
#   
#   print(variable)
#   
#   covVars <- c("Age","Gender") # these values are linked to the cities, which are also the sources of the Pop values
#   
#   water_tests[[ variable ]] <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy, meta.healthy, 100,
#                                                                   trait = variable, chosenControls = all_subSamps$Water_type_home,
#                                                                   covVars = covVars, subSampObjects = waterType.subSamp.objects,
#                                                                   dontReturnMTab = T, dontReturnSamps = T, anova_only = T)
# }
# saveRDS(water_tests, file = sprintf("%s/R_objects/water_tests.rds", p2_dir))

water_tests <- readRDS(file = sprintf("%s/R_objects/water_tests.rds", p2_dir))

# ****************************************************************************************************************** #



ignoreConts <- c("Li","Sr")
cont_water_data.good <- cont_water_data[ ! cont_water_data %in% ignoreConts ]

gen.pMeans.water   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "Genus", "mean"))
phy.pMeans.water   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "Phylum", "mean"))
conts.pMeans.water <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "contVar", "mean"))

gen.num_sig.water   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "Genus", "num_sig"))
phy.num_sig.water   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "Phylum", "num_sig"))
conts.num_sig.water <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests[[w]], w, "contVar", "num_sig"))

# ************************************************************************************** #

gen.num_sig.water.igco <- gen.num_sig.water[ , ! colnames(gen.num_sig.water) %in% ignoreConts ]
gen.pMeans.water.igco <- gen.pMeans.water[ , ! colnames(gen.pMeans.water) %in% ignoreConts ]

freqSig.water <- gen.num_sig.water.igco[ rowSums(gen.num_sig.water.igco >= 95) > 0, ]
freqSig.water.ps <- gen.pMeans.water.igco[ rowSums(gen.num_sig.water.igco >= 95) > 0, ]
# freqSig.water.ps[ freqSig.water.ps > 0.05 ] <- 1
freqSig.water.psLog <- -log(freqSig.water.ps)

library(ComplexHeatmap)
library(circlize)
Heatmap(freqSig.water.psLog)

gensToCheck <- rownames(freqSig.water)

# freqSig.water.glms <- list()
# for (gen in gensToCheck) {
#   print(gen)
# 
#   gt <- gloms_clr$Genus
#   rownames(gt) <- gsub("-",".", rownames(gt))
#   dataTab <- cbind(meta.healthy, t(gt[, rownames(meta.healthy) ]))
# 
#   freqSig.water.glms[[ gen ]] <- as.list(as.data.frame( sapply(cont_water_data.good, function(con)
#     lapply(all_subSamps$Water_type_home, function(i)
#       summary(glm(as.formula(sprintf("%s ~ %s + Age + Gender", gen, con)),
#           data=dataTab[ i , ]))
#     )) ))
# }
# saveRDS(freqSig.water.glms, file = sprintf("%s/R_objects/freqSig.water.glms.rds", p2_dir))

freqSig.water.glms <- readRDS(file = sprintf("%s/R_objects/freqSig.water.glms.rds", p2_dir))

get_mean_CIs <- function(gen, con, glms) {
  CIs <- lapply(glms[[ gen ]][[ con ]], function(i) {
    coefs <- i$coefficients
    est <- coefs[ con, "Estimate"]
    se  <- coefs[ con, "Std. Error"]
    
    ci.25  <- est - (1.96*se)
    ci.975 <- est + (1.96*se)
    
    return(list("est"=est, "ci.25"=ci.25, "ci.975"=ci.975))
  })
  
  CIs.est <- mean(unlist(lapply(CIs, function(x) x$est)))
  CIs.25  <- mean(unlist(lapply(CIs, function(x) x$ci.25)))
  CIs.975 <- mean(unlist(lapply(CIs, function(x) x$ci.975)))
  
  return(list("est"=CIs.est, "ci.25"=CIs.25, "ci.975"=CIs.975))
}

freqSig.water.cis <- as.list(as.data.frame( sapply(gensToCheck, function(gen)
  as.list(as.data.frame( sapply(cont_water_data.good, function(con) get_mean_CIs(gen, con, freqSig.water.glms))) )) ))
freqSig.water.cis <- lapply(freqSig.water.cis, function(g) lapply(g, unlist))



freqSig.water.psLog.direction <- freqSig.water.sigs <- matrix(0, nrow = nrow(freqSig.water.psLog), ncol = ncol(freqSig.water.psLog))
rownames(freqSig.water.psLog.direction) <- rownames(freqSig.water.sigs) <- rownames(freqSig.water.psLog)
colnames(freqSig.water.psLog.direction) <- colnames(freqSig.water.sigs) <- colnames(freqSig.water.psLog)

for (gen in rownames(freqSig.water.psLog)) {
  for (con in colnames(freqSig.water.psLog)) {
    
    if (freqSig.water.ps[ gen, con ] > 0.05) {
      freqSig.water.sigs[ gen, con ] <- ""
    } else {
      freqSig.water.sigs[ gen, con ] <- "+"
    }
    
    if (freqSig.water.cis[[ gen ]][[ con ]][ "est" ] < 0) {
      freqSig.water.psLog.direction[ gen, con ] <- (-1) * freqSig.water.psLog[ gen, con ]
    } else {
      freqSig.water.psLog.direction[ gen, con ] <- freqSig.water.psLog[ gen, con ]
    }
  }
}
# freqSig.water.psLog.direction[ freqSig.water.ps > 0.05 ] <- 0.0

# adjust name since its too long
rownames(freqSig.water.psLog.direction) <- rownames(freqSig.water.sigs) <-
  gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium","ANPR",rownames(freqSig.water.sigs))

Heatmap(freqSig.water.psLog.direction)


hm <- Heatmap(freqSig.water.psLog.direction) # get object with default row/col clustering
colminmax <- max(abs(freqSig.water.psLog.direction))

hm2 <- Heatmap(freqSig.water.psLog.direction, 
        row_order = rev(row_order(hm)), # flips rows and cols so that most signif are in topleft
        column_order = rev(column_order(hm)), 
        cell_fun = function(j, i, x, y, w, h, col) { # add pluses to cells that are significant
          grid.text(freqSig.water.sigs[i, j], x, y)
        }, 
        km=3, column_km=2, # break into the primary branches of dendrogram (see default heatmap without reordering rows/cols)
        col = colorRamp2(c(-colminmax, 0, colminmax), c("blue","white","red")), # ensure 0 is midpoint of colors
        heatmap_legend_param = list(title="log(p-value)", direction="horizontal"))
draw(hm2, heatmap_legend_side="bottom")


# ************** #
# Comapring to figure 7 from the SLL1 paper:
# https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click on image to zoom&p=PMC3&id=6284318_40168_2018_592_Fig7_HTML.jpg
same_as_SLL1 <- list(
  # "Acinetobacter"=c(),
  # "Mesorhizobium"=c(),
  "Variovorax"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","HCO3","Na","SO4"),
  "Psuedomonas"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Na","SO4"),
  "Ralstonia"=c("Alcalinity","Ca","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","Na","SO4"),
  # "Curvibacter"=c(),
  "ANPR"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"),
  # "Hyphomicrobium"=c(),
  "Bradyrhizobium"=c("Ca","Hardness","NO3"),
  "Porphyromonas"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","HCO3","Na","SO4"),
  "Catonella"=c()
  # "Filifactor"=c(),
  # "Fretibacterium"=c()
)

diff_from_SLL1 <- list(
  # "Acinetobacter"=list("SLL1"=c(), "SLL2"=c()),
  # "Mesorhizobium"=list("SLL1"=c(), "SLL2"=c()),
  "Variovorax"=list(
    "SLL1"=c("water_pH"), 
    "SLL2"=c("Ca","Cl","F","Hardness","K","Mg","NO3")),
  "Psuedomonas"=list(
    "SLL1"=c(), 
    "SLL2"=c("Ca","Cl","Hardness","K","Mg","NO3","water_pH")),
  "Ralstonia"=list(
    "SLL1"=c("NO3"), 
    "SLL2"=c("K","water_pH")),
  # "Curvibacter"=list("SLL1"=c(), "SLL2"=c()),
  "ANPR"=list(
    "SLL1"=c(), 
    "SLL2"=c("Cl","F","K","Na","NO3")),
  # "Hyphomicrobium"=list("SLL1"=c(), "SLL2"=c()),
  "Bradyrhizobium"=list(
    "SLL1"=c("K"), 
    "SLL2"=c("Alcalinity","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","Hardness","HCO3","Mg","Na","SO4")),
  "Porphyromonas"=list(
    "SLL1"=c("F"),
    "SLL2"=c("Ca","Hardness","Mg","NO3")),
  "Catonella"=list(
    "SLL1"=c("Cl"),
    "SLL2"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"))
  # "Filifactor"=list("SLL1"=c(), "SLL2"=c()),
  # "Fretibacterium"=list("SLL1"=c(), "SLL2"=c())
)

# ************** #
summary(glm(as.formula(sprintf("%s ~ %s","Variovorax", "Conductivity + Age + Gender")), 
            data=cbind(meta.healthy[ all_subSamps$Age_bins$`1`, ], t(gloms_clr$Genus[,all_subSamps$Age_bins$`1`]))))$coefficients

# ************** #
plot_data.cont("Acinetobacter", "Li", "TvsQ", "Genus", rownames(meta.healthy), only_cont, gloms_clr, meta.healthy)

# ************************************************************************************** #






# ****************************************************************************************************************** #
# Water values - no bottled drinkers ####

# run these tests using the Water_type_home subsamples, exclude bottled water drinkers
meta.healthy_nonBottle <- meta.healthy[ ! is.na(meta.healthy$Water_type_home) &
                                          meta.healthy$Water_type_home != "Embotellada", ]
phy.healthy_nonBottle <- phyloseq(otu_table(prune_samples(rownames(meta.healthy_nonBottle), SLL2)),
                                  sample_data(meta.healthy_nonBottle),
                                  phy_tree(prune_samples(rownames(meta.healthy_nonBottle), SLL2)),
                                  tax_table(prune_samples(rownames(meta.healthy_nonBottle), SLL2)))

nonBottle.subSamps <- lapply(all_subSamps$Water_type_home, function(i) i[i %in% rownames(meta.healthy_nonBottle)] )


# water_tests_nonBottle <- list()
# for (variable in c(cont_water_data[ ! cont_water_data %in% c("Li","Sr")])) {
#   # for (variable in c("Sr","Hardness","Alcalinity")) {
#   
#   print(variable)
#   
#   covVars <- c("Age","Gender") # these values are linked to the cities, which are also the sources of the Pop values
#   
#   water_tests_nonBottle[[ variable ]] <- run_full_subsampling_calcs.ageBins("healthy", phy.healthy_nonBottle, meta.healthy_nonBottle, 100,
#                                                                             trait = variable, chosenControls = nonBottle.subSamps,
#                                                                             covVars = covVars, #subSampObjects = waterType.subSamp.objects,
#                                                                             dontReturnMTab = T, dontReturnSamps = T, anova_only = T)
# }
# saveRDS(water_tests_nonBottle, file = sprintf("%s/R_objects/water_tests_nonBottle.rds", p2_dir))

water_tests_nonBottle <- readRDS(file = sprintf("%s/R_objects/water_tests_nonBottle.rds", p2_dir))

# ************************************************************************************** #

ignoreConts <- c("Li","Sr")
cont_water_data.good <- cont_water_data[ ! cont_water_data %in% ignoreConts ]

gen.pMeans.water_nonBottle   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "Genus", "mean"))
phy.pMeans.water_nonBottle   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "Phylum", "mean"))
conts.pMeans.water_nonBottle <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "contVar", "mean"))

gen.num_sig.water_nonBottle   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "Genus", "num_sig"))
phy.num_sig.water_nonBottle   <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "Phylum", "num_sig"))
conts.num_sig.water_nonBottle <- sapply(cont_water_data.good, function(w) get_subs_anova_per_taxon(water_tests_nonBottle[[w]], w, "contVar", "num_sig"))

# ************************************************************************************** #

freqSig.water_nonBottle <- gen.num_sig.water_nonBottle[ rowSums(gen.num_sig.water_nonBottle >= 95) > 0, ]
freqSig.water_nonBottle.ps <- gen.pMeans.water_nonBottle[ rowSums(gen.num_sig.water_nonBottle >= 95) > 0, ]
# freqSig.water_nonBottle.ps[ freqSig.water_nonBottle.ps > 0.05 ] <- 1
freqSig.water_nonBottle.psLog <- -log(freqSig.water_nonBottle.ps)

library(ComplexHeatmap)
library(circlize)
Heatmap(freqSig.water_nonBottle.psLog)

gensToCheck_nonBottle <- rownames(freqSig.water_nonBottle)


# freqSig.water_nonBottle.glms <- list()
# for (gen in gensToCheck_nonBottle) {
#   print(gen)
#   
#   gt <- gloms_clr$Genus
#   rownames(gt) <- gsub("-",".", rownames(gt))
#   dataTab <- cbind(meta.healthy_nonBottle, t(gt[, rownames(meta.healthy_nonBottle) ]))
#   
#   freqSig.water_nonBottle.glms[[ gen ]] <- as.list(as.data.frame( sapply(cont_water_data.good, function(con)
#     lapply(nonBottle.subSamps, function(i)
#       summary(glm(as.formula(sprintf("%s ~ %s + Age + Gender", gen, con)),
#                   data=dataTab[ i , ]))
#     )) ))
# }
# saveRDS(freqSig.water_nonBottle.glms, file = sprintf("%s/R_objects/freqSig.water_nonBottle.glms.rds", p2_dir))

freqSig.water_nonBottle.glms <- readRDS(file = sprintf("%s/R_objects/freqSig.water_nonBottle.glms.rds", p2_dir))



freqSig.water_nonBottle.cis <- as.list(as.data.frame( sapply(gensToCheck_nonBottle, function(gen)
  as.list(as.data.frame( sapply(cont_water_data.good, function(con) get_mean_CIs(gen, con, freqSig.water_nonBottle.glms))) )) ))
freqSig.water_nonBottle.cis <- lapply(freqSig.water_nonBottle.cis, function(g) lapply(g, unlist))



freqSig.water_nonBottle.psLog.direction <- freqSig.water_nonBottle.sigs <- matrix(0, nrow = nrow(freqSig.water_nonBottle.psLog), 
                                                                                  ncol = ncol(freqSig.water_nonBottle.psLog))
rownames(freqSig.water_nonBottle.psLog.direction) <- rownames(freqSig.water_nonBottle.sigs) <- rownames(freqSig.water_nonBottle.psLog)
colnames(freqSig.water_nonBottle.psLog.direction) <- colnames(freqSig.water_nonBottle.sigs) <- colnames(freqSig.water_nonBottle.psLog)

for (gen in rownames(freqSig.water_nonBottle.psLog)) {
  for (con in colnames(freqSig.water_nonBottle.psLog)) {
    
    if (freqSig.water_nonBottle.ps[ gen, con ] > 0.05) {
      freqSig.water_nonBottle.sigs[ gen, con ] <- ""
    } else {
      freqSig.water_nonBottle.sigs[ gen, con ] <- "+"
    }
    
    if (freqSig.water_nonBottle.cis[[ gen ]][[ con ]][ "est" ] < 0) {
      freqSig.water_nonBottle.psLog.direction[ gen, con ] <- (-1) * freqSig.water_nonBottle.psLog[ gen, con ]
    } else {
      freqSig.water_nonBottle.psLog.direction[ gen, con ] <- freqSig.water_nonBottle.psLog[ gen, con ]
    }
  }
}
# freqSig.water_nonBottle.psLog.direction[ freqSig.water_nonBottle.ps > 0.05 ] <- 0.0

# adjust name since its too long
rownames(freqSig.water_nonBottle.psLog.direction) <- rownames(freqSig.water_nonBottle.sigs) <-
  gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium","ANPR",rownames(freqSig.water_nonBottle.sigs))

Heatmap(freqSig.water_nonBottle.psLog.direction)


hm <- Heatmap(freqSig.water_nonBottle.psLog.direction) # get object with default row/col clustering
colminmax <- max(abs(freqSig.water_nonBottle.psLog.direction))

hm2 <- Heatmap(freqSig.water_nonBottle.psLog.direction, 
               row_order = rev(row_order(hm)), # flips rows and cols so that most signif are in topleft
               column_order = rev(column_order(hm)), 
               cell_fun = function(j, i, x, y, w, h, col) { # add pluses to cells that are significant
                 grid.text(freqSig.water_nonBottle.sigs[i, j], x, y)
               }, 
               km=3, column_km=2, # break into the primary branches of dendrogram (see default heatmap without reordering rows/cols)
               col = colorRamp2(c(-colminmax, 0, colminmax), c("blue","white","red")), # ensure 0 is midpoint of colors
               heatmap_legend_param = list(title="log(p-value)", direction="horizontal"))
draw(hm2, heatmap_legend_side="bottom")


# ************** #
# Comapring to figure 7 from the SLL1 paper:
# https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click on image to zoom&p=PMC3&id=6284318_40168_2018_592_Fig7_HTML.jpg
same_as_SLL1 <- list(
  # "Acinetobacter"=c(),
  # "Mesorhizobium"=c(),
  "Variovorax"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","HCO3","Na","SO4"),
  "Psuedomonas"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Na","SO4"),
  "Ralstonia"=c("Alcalinity","Ca","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","Na","NO3","SO4"),
  # "Curvibacter"=c(),
  "ANPR"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"),
  # "Hyphomicrobium"=c(),
  "Bradyrhizobium"=c("Ca","Hardness","NO3"),
  "Porphyromonas"=c(),
  "Catonella"=c()
  # "Filifactor"=c(),
  # "Fretibacterium"=c()
)

diff_from_SLL1 <- list(
  # "Acinetobacter"=list("SLL1"=c(), "SLL2"=c()),
  # "Mesorhizobium"=list("SLL1"=c(), "SLL2"=c()),
  "Variovorax"=list(
    "SLL1"=c("water_pH"), 
    "SLL2"=c("Ca","Cl","F","Hardness","K","Mg","NO3")),
  "Psuedomonas"=list(
    "SLL1"=c(), 
    "SLL2"=c("Ca","Cl","Hardness","K","Mg","NO3","water_pH")),
  "Ralstonia"=list(
    "SLL1"=c(), 
    "SLL2"=c("K","water_pH")),
  # "Curvibacter"=list("SLL1"=c(), "SLL2"=c()),
  "ANPR"=list(
    "SLL1"=c(), 
    "SLL2"=c("Cl","F","Na","NO3")),
  # "Hyphomicrobium"=list("SLL1"=c(), "SLL2"=c()),
  "Bradyrhizobium"=list(
    "SLL1"=c("K"), 
    "SLL2"=c("Alcalinity","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","Hardness","HCO3","Mg","Na","SO4")),
  "Porphyromonas"=list(
    "SLL1"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Na","SO4"),
    "SLL2"=c("NO3")),
  "Catonella"=list(
    "SLL1"=c("Cl"),
    "SLL2"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"))
  # "Filifactor"=list("SLL1"=c(), "SLL2"=c()),
  # "Fretibacterium"=list("SLL1"=c(), "SLL2"=c())
)

# ****************************************************************************************************************** #





# ****************************************************************************************************************** #
# Check water values only in those teens with same age as SLL1 ####
meta.sll1age <- meta.healthy[ ! is.na(meta.healthy$Age) & meta.healthy$Age >= 13 & meta.healthy$Age < 16, ]

# sll1age.anovas <- list()
# for (con in cont_water_data[ ! cont_water_data %in% ignoreConts]) {
#   
#   print( con )
#   sll1age.anovas[[ con ]] <- list()
#   
#   for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
#     
#     # ******************** #
#     sll1age.anovas[[ con ]][[ tl ]] <- get_lm( c(con, "Age", "Gender"), 
#                                                tl, meta.sll1age, gloms_clr[[ tl ]][ , rownames(meta.sll1age) ], 
#                                                NULL, noRemove = T, rerun.nonSig = F)
#     # ******************** #
#   }
# }
# saveRDS(sll1age.anovas, file = sprintf("%s/R_objects/sll1age.anovas.rds", p2_dir))

sll1age.anovas <- readRDS(file = sprintf("%s/R_objects/sll1age.anovas.rds", p2_dir))

# ************************************************************ #
ignoreConts <- c("Li","Sr")

sll1age.sigs <- matrix(0, nrow = nrow(gloms_clr$Genus), ncol = length(cont_water_data[ ! cont_water_data %in% ignoreConts]))
rownames(sll1age.sigs) <- gsub("-",".",rownames(gloms_clr$Genus))
colnames(sll1age.sigs) <- cont_water_data[ ! cont_water_data %in% ignoreConts]

for (gen in rownames(sll1age.sigs)) {
  for (con in colnames(sll1age.sigs)) {
    sll1age.sigs[ gen, con ] <- sll1age.anovas[[ con ]]$Genus[ gen, con ]
  }
}

sll1age.sigs <- sll1age.sigs[ rowSums(sll1age.sigs < 0.005) > 0, ]


sll1age.sigs.Log <- -log(sll1age.sigs)


# ************************************************************ #

# sll1age.glms <- list()
# for (gen in rownames(sll1age.sigs)) {
#   print(gen)
#   sll1age.glms[[ gen ]] <- list()
#   
#   gt <- gloms_clr$Genus[ , rownames(meta.sll1age) ]
#   rownames(gt) <- gsub("-",".", rownames(gt))
#   dataTab <- cbind(meta.sll1age, t(gt[, rownames(meta.sll1age) ]))
#   
#   for (con in colnames(sll1age.sigs)) {
#     sll1age.glms[[ gen ]][[ con ]] <- summary(glm(as.formula(sprintf("%s ~ %s + Age + Gender", gen, con)),
#                                                   data=dataTab))
#   }
# }
# saveRDS(sll1age.glms, file = sprintf("%s/R_objects/sll1age.glms.rds", p2_dir))

sll1age.glms <- readRDS(file = sprintf("%s/R_objects/sll1age.glms.rds", p2_dir))


get_sll1age_CIs <- function(gen, con, glms) {
  coefs <- glms[[ gen ]][[ con ]]$coefficients
  est <- coefs[ con, "Estimate"]
  se  <- coefs[ con, "Std. Error"]
  
  ci.25  <- est - (1.96*se)
  ci.975 <- est + (1.96*se)
  
  return(list("est"=est, "ci.25"=ci.25, "ci.975"=ci.975))
}

sll1age.sigs.cis <- as.list(as.data.frame( sapply(rownames(sll1age.sigs), function(gen)
  as.list(as.data.frame( sapply(colnames(sll1age.sigs), function(con) get_sll1age_CIs(gen, con, sll1age.glms))) )) ))
sll1age.sigs.cis <- lapply(sll1age.sigs.cis, function(g) lapply(g, unlist))


# ************************************************************ #



sll1age.sigs.Log.direction <- sll1age.sigs.signs <- matrix(0, nrow = nrow(sll1age.sigs.Log), ncol = ncol(sll1age.sigs.Log))
rownames(sll1age.sigs.Log.direction) <- rownames(sll1age.sigs.signs) <- rownames(sll1age.sigs.Log)
colnames(sll1age.sigs.Log.direction) <- colnames(sll1age.sigs.signs) <- colnames(sll1age.sigs.Log)

for (gen in rownames(sll1age.sigs.Log)) {
  for (con in colnames(sll1age.sigs.Log)) {
    
    if (sll1age.sigs[ gen, con ] > 0.05) {
      sll1age.sigs.signs[ gen, con ] <- ""
    } else {
      sll1age.sigs.signs[ gen, con ] <- "+"
    }
    
    if (sll1age.sigs.cis[[ gen ]][[ con ]][ "est" ] < 0) {
      sll1age.sigs.Log.direction[ gen, con ] <- (-1) * sll1age.sigs.Log[ gen, con ]
    } else {
      sll1age.sigs.Log.direction[ gen, con ] <- sll1age.sigs.Log[ gen, con ]
    }
  }
}


rownames(sll1age.sigs.Log.direction) <- rownames(sll1age.sigs.signs) <-
  gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium","ANPR",rownames(sll1age.sigs.signs))

library(ComplexHeatmap)
library(circlize)
Heatmap(sll1age.sigs.Log.direction)


hm <- Heatmap(sll1age.sigs.Log.direction) # get object with default row/col clustering
colminmax <- max(abs(sll1age.sigs.Log.direction))

hm2 <- Heatmap(sll1age.sigs.Log.direction, 
               row_order = rev(row_order(hm)), # flips rows and cols so that most signif are in topleft
               column_order = rev(column_order(hm)), 
               cell_fun = function(j, i, x, y, w, h, col) { # add pluses to cells that are significant
                 grid.text(sll1age.sigs.signs[i, j], x, y)
               }, 
               km=3, column_km=3, # break into the primary branches of dendrogram (see default heatmap without reordering rows/cols)
               # col = colorRamp2(c(-colminmax, 0, colminmax), c("blue","white","red")), # ensure 0 is midpoint of colors
               heatmap_legend_param = list(title="log(p-value)", direction="horizontal"))
draw(hm2, heatmap_legend_side="bottom")

# ************** #
# Comapring to figure 7 from the SLL1 paper:
# https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click on image to zoom&p=PMC3&id=6284318_40168_2018_592_Fig7_HTML.jpg
sll1age.same_as_SLL1 <- list(
  "Variovorax"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","HCO3","Na","SO4"),
  "Psuedomonas"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Na","SO4"),
  "Ralstonia"=c("Alcalinity","Ca","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness",
                "HCO3","Mg","Na","NO3","SO4"),
  "ANPR"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"),
  "Bradyrhizobium"=c("Hardness","NO3","K - opposite trends"),
  "Porphyromonas"=c("Alcalinity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","SO4"),
  "Rothia"=c(),
  "Eikenella"=c()
)

sll1age.diff_from_SLL1 <- list(
  "Variovorax"=list(
    "SLL1"=c("water_pH"), 
    "SLL2"=c("Ca","Cl","F","Hardness","K","Mg","NO3")),
  "Psuedomonas"=list(
    "SLL1"=c(), 
    "SLL2"=c("Ca","Cl","Hardness","Mg","NO3","water_pH")),
  "Ralstonia"=list(
    "SLL1"=c(), 
    "SLL2"=c("F","K","water_pH")),
  "ANPR"=list(
    "SLL1"=c(), 
    "SLL2"=c("Cl","F","K","Na","NO3","water_pH")),
  "Bradyrhizobium"=list(
    "SLL1"=c("Ca"), 
    "SLL2"=c("Alcalinity","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Mg","SO4"))#,
  "Porphyromonas"=list(
    "SLL1"=c("Conductivity","F","Na"), 
    "SLL2"=c("NO3")),
  "Rothia"=list(
    "SLL1"=c("K"),
    "SLL2"=c("NO3")),
  "Eikenella"=list(
    "SLL1"=c("water_pH"),
    "SLL2"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"))
)

# ************** #
# save tables of pvalues, estimates, 2.5 and 97.5% CIs for the heatmaps


# ************************************************************************************** #





# ****************************************************************************************************************** #
# Check water values only in those teens with same age as SLL1 - no bottled drinkers ####
ignoreConts <- c("Li","Sr")
meta.sll1age_nonBottle <- meta.healthy_nonBottle[ ! is.na(meta.healthy_nonBottle$Age) & meta.healthy_nonBottle$Age >= 13 & 
                                                    meta.healthy_nonBottle$Age < 16, ]

# sll1age.anovas_nonBottle <- list()
# for (con in cont_water_data[ ! cont_water_data %in% ignoreConts]) {
# 
#   print( con )
#   sll1age.anovas_nonBottle[[ con ]] <- list()
# 
#   for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
#     
#     # ******************** #
#     sll1age.anovas_nonBottle[[ con ]][[ tl ]] <- get_lm( c(con, "Age", "Gender"),
#                                                          tl, meta.sll1age_nonBottle, gloms_clr[[ tl ]][ , rownames(meta.sll1age_nonBottle) ],
#                                                          NULL, noRemove = T, rerun.nonSig = F)
#     # ******************** #
#   }
# }
# saveRDS(sll1age.anovas_nonBottle, file = sprintf("%s/R_objects/sll1age.anovas_nonBottle.rds", p2_dir))

sll1age.anovas_nonBottle <- readRDS(file = sprintf("%s/R_objects/sll1age.anovas_nonBottle.rds", p2_dir))

# ************************************************************ #


sll1age.sigs_nonBottle <- matrix(0, nrow = nrow(gloms_clr$Genus), ncol = length(cont_water_data[ ! cont_water_data %in% ignoreConts]))
rownames(sll1age.sigs_nonBottle) <- gsub("-",".",rownames(gloms_clr$Genus))
colnames(sll1age.sigs_nonBottle) <- cont_water_data[ ! cont_water_data %in% ignoreConts]

for (gen in rownames(sll1age.sigs_nonBottle)) {
  for (con in colnames(sll1age.sigs_nonBottle)) {
    sll1age.sigs_nonBottle[ gen, con ] <- sll1age.anovas_nonBottle[[ con ]]$Genus[ gen, con ]
  }
}

sll1age.sigs_nonBottle <- sll1age.sigs_nonBottle[ rowSums(sll1age.sigs_nonBottle < 0.005) > 0, ]


sll1age.sigs_nonBottle.Log <- -log(sll1age.sigs_nonBottle)


# ************************************************************ #

# sll1age.glms_nonBottle <- list()
# for (gen in rownames(sll1age.sigs_nonBottle)) {
#   print(gen)
#   sll1age.glms_nonBottle[[ gen ]] <- list()
# 
#   gt <- gloms_clr$Genus[ , rownames(meta.sll1age_nonBottle) ]
#   rownames(gt) <- gsub("-",".", rownames(gt))
#   dataTab <- cbind(meta.sll1age_nonBottle, t(gt[, rownames(meta.sll1age_nonBottle) ]))
# 
#   for (con in colnames(sll1age.sigs_nonBottle)) {
#     sll1age.glms_nonBottle[[ gen ]][[ con ]] <- summary(glm(as.formula(sprintf("%s ~ %s + Age + Gender", gen, con)),
#                                                             data=dataTab))
#   }
# }
# saveRDS(sll1age.glms_nonBottle, file = sprintf("%s/R_objects/sll1age.glms_nonBottle.rds", p2_dir))

sll1age.glms_nonBottle <- readRDS(file = sprintf("%s/R_objects/sll1age.glms_nonBottle.rds", p2_dir))


get_sll1age_CIs <- function(gen, con, glms) {
  coefs <- glms[[ gen ]][[ con ]]$coefficients
  est <- coefs[ con, "Estimate"]
  se  <- coefs[ con, "Std. Error"]
  
  ci.25  <- est - (1.96*se)
  ci.975 <- est + (1.96*se)
  
  return(list("est"=est, "ci.25"=ci.25, "ci.975"=ci.975))
}

sll1age.sigs_nonBottle.cis <- as.list(as.data.frame( sapply(rownames(sll1age.sigs_nonBottle), function(gen)
  as.list(as.data.frame( sapply(colnames(sll1age.sigs_nonBottle), function(con) get_sll1age_CIs(gen, con, sll1age.glms_nonBottle))) )) ))
sll1age.sigs_nonBottle.cis <- lapply(sll1age.sigs_nonBottle.cis, function(g) lapply(g, unlist))


# ************************************************************ #



sll1age.sigs_nonBottle.Log.direction <- sll1age.sigs_nonBottle.signs <- matrix(0, nrow = nrow(sll1age.sigs_nonBottle.Log), 
                                                                               ncol = ncol(sll1age.sigs_nonBottle.Log))
rownames(sll1age.sigs_nonBottle.Log.direction) <- rownames(sll1age.sigs_nonBottle.signs) <- rownames(sll1age.sigs_nonBottle.Log)
colnames(sll1age.sigs_nonBottle.Log.direction) <- colnames(sll1age.sigs_nonBottle.signs) <- colnames(sll1age.sigs_nonBottle.Log)

for (gen in rownames(sll1age.sigs_nonBottle.Log)) {
  for (con in colnames(sll1age.sigs_nonBottle.Log)) {
    
    if (sll1age.sigs_nonBottle[ gen, con ] > 0.05) {
      sll1age.sigs_nonBottle.signs[ gen, con ] <- ""
    } else {
      sll1age.sigs_nonBottle.signs[ gen, con ] <- "+"
    }
    
    if (sll1age.sigs_nonBottle.cis[[ gen ]][[ con ]][ "est" ] < 0) {
      sll1age.sigs_nonBottle.Log.direction[ gen, con ] <- (-1) * sll1age.sigs_nonBottle.Log[ gen, con ]
    } else {
      sll1age.sigs_nonBottle.Log.direction[ gen, con ] <- sll1age.sigs_nonBottle.Log[ gen, con ]
    }
  }
}


rownames(sll1age.sigs_nonBottle.Log.direction) <- rownames(sll1age.sigs_nonBottle.signs) <-
  gsub("Allorhizobium.Neorhizobium.Pararhizobium.Rhizobium","ANPR",rownames(sll1age.sigs_nonBottle.signs))

library(ComplexHeatmap)
library(circlize)
Heatmap(sll1age.sigs_nonBottle.Log.direction)


hm <- Heatmap(sll1age.sigs_nonBottle.Log.direction) # get object with default row/col clustering
colminmax <- max(abs(sll1age.sigs_nonBottle.Log.direction))

hm2 <- Heatmap(sll1age.sigs_nonBottle.Log.direction, 
               row_order = rev(row_order(hm)), # flips rows and cols so that most signif are in topleft
               column_order = rev(column_order(hm)), 
               cell_fun = function(j, i, x, y, w, h, col) { # add pluses to cells that are significant
                 grid.text(sll1age.sigs_nonBottle.signs[i, j], x, y)
               }, 
               km=2, column_km=3, # break into the primary branches of dendrogram (see default heatmap without reordering rows/cols)
               # col = colorRamp2(c(-colminmax, 0, colminmax), c("blue","white","red")), # ensure 0 is midpoint of colors
               heatmap_legend_param = list(title="log(p-value)", direction="horizontal"))
draw(hm2, heatmap_legend_side="bottom")

# ************** #
# Comapring to figure 7 from the SLL1 paper:
# https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click on image to zoom&p=PMC3&id=6284318_40168_2018_592_Fig7_HTML.jpg
sll1age.same_as_SLL1 <- list(
  "Variovorax"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","HCO3","Na","SO4"),
  "Psuedomonas"=c("Alcalinity","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Na","SO4"),
  "Ralstonia"=c("Alcalinity","Ca","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness",
                "HCO3","Mg","Na","NO3","SO4"),
  "ANPR"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"),
  "Bradyrhizobium"=c("Hardness","NO3","K - opposite trends"),
  "Catonella"=c()
)

sll1age.diff_from_SLL1 <- list(
  "Variovorax"=list(
    "SLL1"=c("water_pH"), 
    "SLL2"=c("Ca","Cl","F","Hardness","Mg","NO3")),
  "Psuedomonas"=list(
    "SLL1"=c(), 
    "SLL2"=c("Ca","Cl","Hardness","Mg","NO3","water_pH")),
  "Ralstonia"=list(
    "SLL1"=c(), 
    "SLL2"=c()),
  "ANPR"=list(
    "SLL1"=c(), 
    "SLL2"=c("Cl","K","Na","NO3")),
  "Bradyrhizobium"=list(
    "SLL1"=c("Ca"), 
    "SLL2"=c("Alcalinity","Cl","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","F","HCO3","Mg","SO4"))#,
  "Catonella"=list(
    "SLL1"=c("Cl"),
    "SLL2"=c("Alcalinity","Ca","Conductivity","Dry.residue.at.110.C","Dry.residue.at.180.C","Hardness","HCO3","Mg","SO4"))
)

# ************** #
# save tables of pvalues, estimates, 2.5 and 97.5% CIs for the heatmaps


# ************************************************************************************** #














# ****************************************************************************************************************** ####


# ****************************************************************************************************************** #
# Networks for healthy Stomatotypes ####

# first create columns in sample_data for each diversity value
# for (div.estimate in c("Div.Shannon","Div.Simpson","Faiths.PD","Species_Richness")) {
#   
#   divs <- meta.healthy[ , div.estimate]
#   quarts <- quantile(divs, probs=seq(0,1,0.25))
#   
#   # get samples by diversity quartile
#   q1.samples <- rownames(meta.healthy[ meta.healthy[,div.estimate] < quarts['25%'], ] )
#   q2.samples <- rownames(meta.healthy[ (quarts['25%'] <= meta.healthy[,div.estimate]) & (meta.healthy[,div.estimate] < quarts['50%']), ] )
#   q3.samples <- rownames(meta.healthy[ (quarts['50%'] <= meta.healthy[,div.estimate]) & (meta.healthy[,div.estimate] < quarts['75%']), ] )
#   q4.samples <- rownames(meta.healthy[ meta.healthy[,div.estimate] >= quarts['75%'], ] )
#   
#   # add diversity group to sample_data 
#   meta.healthy[q1.samples, sprintf("healthy.Diversity_group_%s", div.estimate)] <- "Low"
#   meta.healthy[c(q2.samples,q3.samples), sprintf("healthy.Diversity_group_%s", div.estimate)] <- "Average"
#   meta.healthy[q4.samples, sprintf("healthy.Diversity_group_%s", div.estimate)] <- "High"
# }



spiec_dir <- sprintf("%s/R_objects/SpiecEasi", p2_dir)
spiec_sub_dir <- sprintf("%s/variable_subSamps", spiec_dir)


tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

# networks.nonSubs <- list()
# 
# for (netVar in c("Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner","seqGroup",
#                  "Water_type_home","BMI_group",
#                  "healthy.Stomatotype_Aitchison","healthy.Stomatotype_Bray","healthy.Stomatotype_Jaccard",
#                  "healthy.Stomatotype_Weighted_Unifrac","healthy.Stomatotype_Unweighted_Unifrac",
#                  "healthy.Diversity_group_Div.Shannon","healthy.Diversity_group_Div.Simpson",
#                  "healthy.Diversity_group_Faiths.PD","healthy.Diversity_group_Species_Richness")) {
#   print(netVar)
#   networks.nonSubs[[ netVar ]] <- get_networks_for_nonSubs(netVar, tagl, meta.healthy, gloms$Genus, minCount=15, minOcc=20)
#   
#   saveRDS(networks.nonSubs, file = sprintf("%s/networks.nonSubs.rds", spiec_dir))
#   print(Sys.time())
# }

networks.nonSubs <- readRDS(sprintf("%s/networks.nonSubs.rds", spiec_dir))
networks.div_groups <- readRDS(sprintf("%s/networks.div_groups.rds", spiec_dir))
networks.stomatotypes <- readRDS(sprintf("%s/networks.stomatotypes.rds", spiec_dir))

# ****************************************************************************************************************** ####




# ************************************************** #

# Networks for all subSamps ####

tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))
# ************************** #


# ************************** #
# for (variable in names(all_subSamps)[ ! names(all_subSamps) %in% c("seqGroup","Age_groups","Smoker","Braces.binary",
#                                                                    "Cystic_fibrosis","Celiac","Downs_Syndrome","Diabetes",
#                                                                    "Fluoride_toothpaste","Antibiotics")]) {
# for (variable in c("Full_MALDI.Candida")) {#"MALDI.Yeast_detected","Gender")) {
#   print(variable)
# 
#   network_objects <- get_networks_for_subSamps(variable, tagl, all_subSamps, SLL2.meta, gloms$Genus, minCount=15, minOcc=20)
# 
#   # ************************** #
#   dir.create(sprintf("%s/%s", spiec_sub_dir, variable), showWarnings = F)
# 
#   saveRDS(network_objects$se.mb.main, file = sprintf("%s/%s/%s.se.mb.main.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.main.allNodes, file = sprintf("%s/%s/%s.igraphs.main.allNodes.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.main.noEmpty, file = sprintf("%s/%s/%s.igraphs.main.noEmpty.rds", spiec_sub_dir, variable, variable))
# 
#   saveRDS(network_objects$se.mb.all_subs, file = sprintf("%s/%s/%s.se.mb.all_subs.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.all_subs.allNodes, file = sprintf("%s/%s/%s.igraphs.all_subs.allNodes.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.all_subs.noEmpty, file = sprintf("%s/%s/%s.igraphs.all_subs.noEmpty.rds", spiec_sub_dir, variable, variable))
# 
#   saveRDS(network_objects$SE_objects, file = sprintf("%s/%s/%s.SE_objects.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.allNodes, file = sprintf("%s/%s/%s.igraphs.allNodes.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$igraphs.noEmpty, file = sprintf("%s/%s/%s.igraphs.noEmpty.rds", spiec_sub_dir, variable, variable))
# 
#   saveRDS(network_objects$main_v_subs.hammings, file = sprintf("%s/%s/%s.main_v_subs.hammings.rds", spiec_sub_dir, variable, variable))
#   saveRDS(network_objects$subs_subs.hammings, file = sprintf("%s/%s/%s.subs_subs.hammings.rds", spiec_sub_dir, variable, variable))
#   print(Sys.time())
# }
# ************************** #



# ****************************************************************************************************************** #

# get connections for each value within each variable

subs_vars <- names(all_subSamps)[ ! names(all_subSamps) %in% c(names(networks.nonSubs))]#,
                                                               # # there are way too few samples for these variables below
                                                               # "Diabetes","Headaches","Lactose_intolerant","Anemia",
                                                               # "Lung_issues","Kidney_issues","Hypothyroidism")]

gen_connects.subs <- list()
gen_con_freqs.subs <- list()

for (v in subs_vars) {#[!subs_vars %in% c("Smoker","Braces.binary","Age_groups")]) {
  print(v)
  gen_connects.subs[[ v ]] <- connect_edges.subs(v, spiec_sub_dir)
  
  # then for each genus, get the frequency of connections to other genera within each value for the given variable
  gen_con_freqs.subs[[ v ]] <- list()
  for (gen in names(gen_connects.subs[[ v ]]$gen_connect$all$`1`)) {
    gen_con_freqs.subs[[ v ]][[ gen ]] <- get_gen_con_freqs(gen_connects.subs[[ v ]]$gen_connect, gen)
  }
}

# ******* #
gen_connects.nonSubs <- list()
gen_con_freqs.nonSubs <- list()

for (v in names(networks.nonSubs)) {
  print(v)
  gen_connects.nonSubs[[ v ]] <- connect_edges.nonSubs(networks.nonSubs, v)
  
  # then for each genus, get the frequency of connections to other genera within each value for the given variable
  gen_con_freqs.nonSubs[[ v ]] <- list()
  for (gen in names(gen_connects.nonSubs[[ v ]]$all)) {
    gen_con_freqs.nonSubs[[ v ]][[ gen ]] <- get_gen_con_freqs.all_No.or.nonSubs(gen_connects.nonSubs[[ v ]], gen)
  }
}


# ******* #
gen_connects.div_groups <- list()
gen_con_freqs.div_groups <- list()

for (v in names(networks.div_groups)) {
  print(v)
  gen_connects.div_groups[[ v ]] <- connect_edges.nonSubs(networks.div_groups, v)
  
  # then for each genus, get the frequency of connections to other genera within each value for the given variable
  gen_con_freqs.div_groups[[ v ]] <- list()
  for (gen in names(gen_connects.div_groups[[ v ]]$all)) {
    gen_con_freqs.div_groups[[ v ]][[ gen ]] <- get_gen_con_freqs.all_No.or.nonSubs(gen_connects.div_groups[[ v ]], gen)
  }
}


# ******* #
gen_connects.stomatotypes <- list()
gen_con_freqs.stomatotypes <- list()

for (v in names(networks.stomatotypes)) {
  print(v)
  gen_connects.stomatotypes[[ v ]] <- connect_edges.nonSubs(networks.stomatotypes, v)
  
  # then for each genus, get the frequency of connections to other genera within each value for the given variable
  gen_con_freqs.stomatotypes[[ v ]] <- list()
  for (gen in names(gen_connects.stomatotypes[[ v ]]$all)) {
    gen_con_freqs.stomatotypes[[ v ]][[ gen ]] <- get_gen_con_freqs.all_No.or.nonSubs(gen_connects.stomatotypes[[ v ]], gen)
  }
}

# ****************************************************************************************************************** #

# then compare variables to each other, first subs to subs, for each genus,
#     to check how many genus-genus connections occur just in subs,
#     and get the absolute value of difference in number of subSamps where same connection occurs as a score (lower is better)
#  if a connection does not occur in both variables subs (0 in at least 1 variable)
#     will give a score of 100 (so a score of 0 means identical, 100 means totally different)

**** should not only give score, but also weight that score based on the values in each variable
because if two variables have Streptococcus and Corynebacterium associated only 2 times, theyd have a "perfect" score
but its actually pretty meaningless
  **** so instead of a perfecto score being 0, do (100-diff)*mean times assoc, 
so that a perfect score would be 100, if both variables had 100 ( so diff=0, and mean=100 => (100-0)*100 )

# assoc_scores <- list()
# assoc_scores[[ "subs_subs" ]] <- get_assoc_scores(subs_vars, gen_con_freqs.subs, "subs_subs")
# assoc_scores[[ "subs_main" ]] <- get_assoc_scores(subs_vars, gen_con_freqs.subs, "subs_main")
# assoc_scores[[ "main_subs" ]] <- get_assoc_scores(subs_vars, gen_con_freqs.subs, "main_subs")
# saveRDS(assoc_scores, file = sprintf("%s/R_objects/SpiecEasi/assoc_scores.rds", p2_dir))

assoc_scores <- readRDS(file = sprintf("%s/R_objects/SpiecEasi/assoc_scores.rds", p2_dir))



sa_names.ss <- names(assoc_scores$subs_subs)[ sapply(names(assoc_scores$subs_subs), function(v1v2) 
  length(unlist(assoc_scores$subs_subs[[ v1v2 ]])[unlist(assoc_scores$subs_subs[[ v1v2 ]]) > 60])) > 0 ]
strong_assocs.ss <- sapply(sa_names.ss, function(x) {
  sa <- sort(unlist(assoc_scores$subs_subs[[ x ]])[ unlist(assoc_scores$subs_subs[[ x]]) > 60 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})

sa_names.sm <- names(assoc_scores$subs_main)[ sapply(names(assoc_scores$subs_main), function(v1v2) 
  length(unlist(assoc_scores$subs_main[[ v1v2 ]])[unlist(assoc_scores$subs_main[[ v1v2 ]]) > 60])) > 0 ]
strong_assocs.sm <- sapply(sa_names.sm, function(x) {
  sa <- sort(unlist(assoc_scores$subs_main[[ x ]])[ unlist(assoc_scores$subs_main[[ x]]) > 60 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})

sa_names.ms <- names(assoc_scores$main_subs)[ sapply(names(assoc_scores$main_subs), function(v1v2) 
  length(unlist(assoc_scores$main_subs[[ v1v2 ]])[unlist(assoc_scores$main_subs[[ v1v2 ]]) > 60])) > 0 ]
strong_assocs.ms <- sapply(sa_names.ms, function(x) {
  sa <- sort(unlist(assoc_scores$main_subs[[ x ]])[ unlist(assoc_scores$main_subs[[ x]]) > 60 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})



# strongest connections in subs_subs
sort(unlist(strong_assocs.ss))
sort(unlist(strong_assocs.sm))
sort(unlist(strong_assocs.ms))
# ****************************************************************************************************************** #

# get sort of anti-assoc_scores, subtracting freqs between variables for particular associations ####

# "Lung_issues","Kidney_issues","Hypothyroidism")]

anti_assoc_vars <- c("Cystic_fibrosis","Downs_Syndrome","Celiac","Smoker","MALDI.Yeast_detected","Full_MALDI.Candida",#"Gender","Braces.binary",
                     "Hypertension","Antibiotics")
gen_con_freqs.subs <- list()

for (v in anti_assoc_vars) {
  print(v)
  gen_connects.subs <- connect_edges.subs(v, spiec_sub_dir)
  
  # then for each genus, get the frequency of connections to other genera within each value for the given variable
  gen_con_freqs.subs[[ v ]] <- list()
  for (gen in names(gen_connects.subs$gen_connect$all$`1`)) {
    gen_con_freqs.subs[[ v ]][[ gen ]] <- get_gen_con_freqs(gen_connects.subs$gen_connect, gen)
  }
}


# anti_assoc_scores <- list()
# anti_assoc_scores[[ "main_main" ]] <- get_assoc_scores(anti_assoc_vars, gen_con_freqs.subs, "main_main", anti_assoc = T)
# anti_assoc_scores[[ "subs_subs" ]] <- get_assoc_scores(anti_assoc_vars, gen_con_freqs.subs, "subs_subs", anti_assoc = T)
# anti_assoc_scores[[ "subs_main" ]] <- get_assoc_scores(anti_assoc_vars, gen_con_freqs.subs, "subs_main", anti_assoc = T)
# anti_assoc_scores[[ "main_subs" ]] <- get_assoc_scores(anti_assoc_vars, gen_con_freqs.subs, "main_subs", anti_assoc = T)
# saveRDS(anti_assoc_scores, file = sprintf("%s/R_objects/SpiecEasi/anti_assoc_scores.rds", p2_dir))

anti_assoc_scores <- readRDS(file = sprintf("%s/R_objects/SpiecEasi/anti_assoc_scores.rds", p2_dir))




sa_names.mm <- names(anti_assoc_scores$main_main)[ sapply(names(anti_assoc_scores$main_main), function(v1v2) 
  length(unlist(anti_assoc_scores$main_main[[ v1v2 ]])[abs(unlist(anti_assoc_scores$main_main[[ v1v2 ]])) > 80])) > 0 ]
min80_anti_assocs.mm <- sapply(sa_names.mm, function(x) {
  sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) > 80 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})



lapply(min80_anti_assocs.mm[ grepl("Smoker", names(min80_anti_assocs.mm)) ], function(x) {
  names(x) <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","ANPR",names(x))
  x
})

uncs <- taxa_print("Genus", c("unclassified.G77","unclassified.G20","unclassified.G40","unclassified.G43",
                              "unclassified.G75","unclassified.G77","unclassified.G29","unclassified.G32",
                              "unclassified.G7","unclassified.G24","unclassified.G14","unclassified.G17",
                              "unclassified.G72","unclassified.G19"), 
                   tl.limit = "Genus")
as.data.frame(uncs[ with(as.data.frame(uncs), order(Phylum, Class, Order, Family))], row.names=F)





lapply(min80_anti_assocs.mm[ grepl("MALDI.Yeast_detected", names(min80_anti_assocs.mm)) ], function(x) {
  names(x) <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","ANPR",names(x))
  x
})



lapply(min80_anti_assocs.mm[ grepl("Full_MALDI.Candida", names(min80_anti_assocs.mm)) ], function(x) {
  names(x) <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","ANPR",names(x))
  x
})


lapply(min80_anti_assocs.mm[ grepl("Hypertension", names(min80_anti_assocs.mm)) ], function(x) {
  names(x) <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","ANPR",names(x))
  x
})



# **************** #


sa_names.mm <- names(anti_assoc_scores$main_main)[ sapply(names(anti_assoc_scores$main_main), function(v1v2) 
  length(unlist(anti_assoc_scores$main_main[[ v1v2 ]])[abs(unlist(anti_assoc_scores$main_main[[ v1v2 ]])) > 95])) > 0 ]
min95_anti_assocs.mm <- sapply(sa_names.mm, function(x) {
  sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) > 95 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})

# **************** #


min100_anti_assocs.mm <- sapply(names(anti_assoc_scores$main_main), function(x) {
  sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) >= 100 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})

min100_anti_assocs.maldi_as_sub <- sapply(names(anti_assoc_scores$main_main), function(x) {
  if (startsWith(x, "MALDI.Yeast_detected"))
    # use "subs" as comparison for MALDI - in this case subs first
    sa <- sort(unlist(anti_assoc_scores$subs_main[[ x ]])[ abs(unlist(anti_assoc_scores$subs_main[[ x]])) >= 100 ])
  else if (endsWith(x, "MALDI.Yeast_detected"))
    # use "subs" as comparison for MALDI - in this case subs last
    sa <- sort(unlist(anti_assoc_scores$main_subs[[ x ]])[ abs(unlist(anti_assoc_scores$main_subs[[ x]])) >= 100 ])
  else
    sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) >= 100 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})


# sa_names.ss <- names(anti_assoc_scores$subs_subs)[ sapply(names(anti_assoc_scores$subs_subs), function(v1v2) 
#   length(unlist(anti_assoc_scores$subs_subs[[ v1v2 ]])[abs(unlist(anti_assoc_scores$subs_subs[[ v1v2 ]])) > 80])) > 0 ]
# strong_anti_assocs.ss <- sapply(sa_names.ss, function(x) {
#   sa <- sort(unlist(anti_assoc_scores$subs_subs[[ x ]])[ abs(unlist(anti_assoc_scores$subs_subs[[ x]])) > 80 ])
#   # remove the repeated values where the genera are simply switched
#   gg_pairs <- sapply(sort(names(sa)), function(gg) 
#     sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
#   gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
#   return(sort(sa[ gg_noDup ]))
# })
# 
# strong_anti_assocs.ss[ grepl("Smoker", names(strong_anti_assocs.ss)) ]
# strong_anti_assocs.ss[ grepl("MALDI.Yeast_detected", names(strong_anti_assocs.ss)) ]

# **************** #

# sa_names.sm <- names(anti_assoc_scores$subs_main)[ sapply(names(anti_assoc_scores$subs_main), function(v1v2) 
#   length(unlist(anti_assoc_scores$subs_main[[ v1v2 ]])[abs(unlist(anti_assoc_scores$subs_main[[ v1v2 ]])) > 80])) > 0 ]
# strong_anti_assocs.sm <- sapply(sa_names.sm, function(x) {
#   sa <- sort(unlist(anti_assoc_scores$subs_main[[ x ]])[ abs(unlist(anti_assoc_scores$subs_main[[ x]])) > 80 ])
#   # remove the repeated values where the genera are simply switched
#   gg_pairs <- sapply(sort(names(sa)), function(gg) 
#     sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
#   gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
#   return(sort(sa[ gg_noDup ]))
# })
# 
# strong_anti_assocs.sm[ grepl("MALDI.Yeast_detected", names(strong_anti_assocs.sm)) ]




all_anti_assocs.mm <- sapply(names(anti_assoc_scores$main_main), function(x) {
  sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) > 0 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})


all_anti_assocs.maldi_as_sub <- sapply(names(anti_assoc_scores$main_main), function(x) {
  if (startsWith(x, "MALDI.Yeast_detected"))
    # use "subs" as comparison for MALDI - in this case subs first
    sa <- sort(unlist(anti_assoc_scores$subs_main[[ x ]])[ abs(unlist(anti_assoc_scores$subs_main[[ x]])) >= 0 ])
  else if (endsWith(x, "MALDI.Yeast_detected"))
    # use "subs" as comparison for MALDI - in this case subs last
    sa <- sort(unlist(anti_assoc_scores$main_subs[[ x ]])[ abs(unlist(anti_assoc_scores$main_subs[[ x]])) >= 0 ])
  else
    sa <- sort(unlist(anti_assoc_scores$main_main[[ x ]])[ abs(unlist(anti_assoc_scores$main_main[[ x]])) >= 0 ])
  # remove the repeated values where the genera are simply switched
  gg_pairs <- sapply(sort(names(sa)), function(gg) 
    sort(strsplit( gsub("unclassified.G","unclassified-G",gg), "\\.")[[1]]))
  gg_noDup <- colnames(gg_pairs)[ ! duplicated(t(gg_pairs))]
  return(sort(sa[ gg_noDup ]))
})


# **************************************************************************** #
# plot network for Smoker to see the positive and negative associations ####

tagl <- phyloseq(otu_table(gloms$Genus), sample_data(SLL2), tax_table(taxTables.both$Genus))

tagl.Smoker.for_network <- get_tagl_for_networks("Smoker", tagl, all_subSamps, meta.healthy, gloms$Genus)

# first get objects for CF samples
samps.Smoker <- rownames(meta.healthy[ meta.healthy[, "Smoker"]=="Yes", ])
tagl.Smoker  <- subset_samples(tagl.Smoker.for_network, sample_names(tagl.Smoker.for_network) %in% samps.Smoker)
se.mb.Smoker <- readRDS(sprintf("%s/Smoker/Smoker.se.mb.main.rds", spiec_sub_dir))


par(mfrow = c(1,2))

Smoker.net <- get_network_objects("Smoker", se.mb.Smoker, tagl.Smoker, taxTables.both$Genus, 
                                  addLabels = T, Vprop = T, Eprop = T,
                                  addLegend = F,
                                  taxa_to_label = c(
                                    "Actinobacillus","Neisseria", 
                                    "Bergeriella","Delftia",
                                    "Delftia","Mogibacterium","Hydrogenophaga",
                                    "Leptotrichia","Stomatobaculum",
                                    "Rikenellaceae_RC9_gut_group","Treponema", 
                                    "Abiotrophia","Porphyromonas","Neisseria","Streptobacillus",
                                    "Actinobacillus","Gemella",
                                    "Anaeroglobus","Dialister","Fretibacterium","Shuttleworthia",
                                    "Aestuariimicrobium","Streptococcus","Centipeda",
                                    "Cardiobacterium","Lautropia","Kingella",
                                    "Mycoplasma","Treponema"
                                  ))


tagl.MALDI.for_network <- get_tagl_for_networks("MALDI.Yeast_detected", tagl, all_subSamps, meta.healthy, gloms$Genus)

samps.MALDI <- unique(unlist(lapply(all_subSamps$MALDI.Yeast_detected, function(i) i)))
samps.MALDI <- samps.MALDI[ meta.healthy[ samps.MALDI, "MALDI.Yeast_detected"]=="Yes" ]
tagl.MALDI  <- subset_samples(tagl.MALDI.for_network, sample_names(tagl.MALDI.for_network) %in% samps.MALDI)
se.mb.MALDI <- readRDS(sprintf("%s/MALDI.Yeast_detected/MALDI.Yeast_detected.se.mb.main.rds", spiec_sub_dir))


MALDI.net <- get_network_objects("MALDI", se.mb.MALDI, tagl.MALDI, taxTables.both$Genus,
                                 addLabels = T, Vprop = T, Eprop = T,
                                 addLegend = F,
                                 taxa_to_label = c(
                                   "Alloscardovia","Veillonella",
                                   "Bacteroides","Fusobacterium","Oceanivirga",
                                   "Bergeriella","Selenomonas", 
                                   "Campylobacter","Megasphaera","Selenomonas", 
                                   "Chryseobacterium","Gemella","Microbacterium",
                                   "Haemophilus","Megasphaera","Streptococcus",
                                   "Alloprevotella","Gemella",
                                   "Butyrivibrio","Oribacterium",
                                   "Capnocytophaga","Corynebacterium","Eikenella","Lautropia",
                                   "Corynebacterium","Porphyromonas",
                                   "Parvimonas","Treponema"
                                 ))

# ********************** #
# get edge weights of particular associations

nc.Smoker <- symBeta(getOptBeta(se.mb.Smoker))
colnames(nc.Smoker) <- rownames(nc.Smoker) <- taxa_names(tagl.Smoker)
inc.Smoker <- graph.adjacency(nc.Smoker, mode="undirected", add.rownames = T, weighted = T)
edge.weights.Smoker <- E(inc.Smoker)$weight
names(edge.weights.Smoker) <- attr(E(inc.Smoker),"vnames")
names(edge.weights.Smoker) <- sapply(names(edge.weights.Smoker), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))

assocs_of_interest.Smoker <- c("Actinobacillus.Neisseria","Abiotrophia.Neisseria","Delftia.Mogibacterium","Fretibacterium.Parvimonas")
sort(edge.weights.Smoker[ assocs_of_interest.Smoker ])
# ********************** #


nc.MALDI <- symBeta(getOptBeta(se.mb.MALDI))
colnames(nc.MALDI) <- rownames(nc.MALDI) <- taxa_names(tagl.MALDI)
inc.MALDI <- graph.adjacency(nc.MALDI, mode="undirected", add.rownames = T, weighted = T)
edge.weights.MALDI <- E(inc.MALDI)$weight
names(edge.weights.MALDI) <- attr(E(inc.MALDI),"vnames")
names(edge.weights.MALDI) <- sapply(names(edge.weights.MALDI), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))

assocs_of_interest.MALDI <- c("Olsenella.Treponema","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium.Fusobacterium","Actinomyces.Prevotella",
                              "Actinomyces.Treponema","Alloprevotella.Granulicatella","Bergeyella.Capnocytophaga","Bradyrhizobium.Brevundimonas",
                              "Curvibacter.Ralstonia","Hydrogenophaga.Variovorax","Hyphomicrobium.Variovorax","Lachnoanaerobaculum.Peptostreptococcus",
                              "Leptotrichia.Selenomonas","Peptostreptococcus.Pseudomonas","Porphyromonas.Pseudomonas","Ralstonia.Ruminococcaceae_UCG-014",
                              "Ruminococcaceae_UCG-014.Stomatobaculum","Selenomonas.Streptococcus")
sort(edge.weights.MALDI[ assocs_of_interest.MALDI ])
# ********************** #

tagl.Candida.for_network <- get_tagl_for_networks("Full_MALDI.Candida", tagl, all_subSamps, meta.healthy, gloms$Genus)

samps.Candida <- unique(unlist(lapply(all_subSamps$Full_MALDI.Candida, function(i) i)))
samps.Candida <- samps.Candida[ meta.healthy[ samps.Candida, "Full_MALDI.Candida"]=="Yes" ]
tagl.Candida  <- subset_samples(tagl.Candida.for_network, sample_names(tagl.Candida.for_network) %in% samps.Candida)
se.mb.Candida <- readRDS(sprintf("%s/Full_MALDI.Candida/Full_MALDI.Candida.se.mb.main.rds", spiec_sub_dir))

nc.Candida <- symBeta(getOptBeta(se.mb.Candida))
colnames(nc.Candida) <- rownames(nc.Candida) <- taxa_names(tagl.Candida)
inc.Candida <- graph.adjacency(nc.Candida, mode="undirected", add.rownames = T, weighted = T)
edge.weights.Candida <- E(inc.Candida)$weight
names(edge.weights.Candida) <- attr(E(inc.Candida),"vnames")
names(edge.weights.Candida) <- sapply(names(edge.weights.Candida), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))

assocs_of_interest.Candida <- c("Olsenella.Treponema","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium.Fusobacterium","Actinomyces.Prevotella",
                              "Actinomyces.Treponema","Alloprevotella.Granulicatella","Bergeyella.Capnocytophaga","Bradyrhizobium.Brevundimonas",
                              "Curvibacter.Ralstonia","Hydrogenophaga.Variovorax","Hyphomicrobium.Variovorax","Lachnoanaerobaculum.Peptostreptococcus",
                              "Leptotrichia.Selenomonas","Peptostreptococcus.Pseudomonas","Porphyromonas.Pseudomonas","Ralstonia.Ruminococcaceae_UCG-014",
                              "Ruminococcaceae_UCG-014.Stomatobaculum","Selenomonas.Streptococcus",
                              "Peptostreptococcus.Prevotella")
sort(edge.weights.Candida[ assocs_of_interest.Candida ])
# ********************** #



tagl.Hypertension.for_network <- get_tagl_for_networks("Hypertension", tagl, all_subSamps, meta.healthy, gloms$Genus)

samps.Hypertension <- unique(unlist(lapply(all_subSamps$Hypertension, function(i) i)))
samps.Hypertension <- samps.Hypertension[ meta.healthy[ samps.Hypertension, "Hypertension"]=="Yes" ]
tagl.Hypertension  <- subset_samples(tagl.Hypertension.for_network, sample_names(tagl.Hypertension.for_network) %in% samps.Hypertension)
se.mb.Hypertension <- readRDS(sprintf("%s/Hypertension/Hypertension.se.mb.main.rds", spiec_sub_dir))

nc.Hypertension <- symBeta(getOptBeta(se.mb.Hypertension))
# this tagl is strange because it leaves out Suttonella, so will have to place that genus into the list of names
colnames(nc.Hypertension) <- rownames(nc.Hypertension) <- colnames(se.mb.Hypertension$est$data)
inc.Hypertension <- graph.adjacency(nc.Hypertension, mode="undirected", add.rownames = T, weighted = T)
edge.weights.Hypertension <- E(inc.Hypertension)$weight
names(edge.weights.Hypertension) <- attr(E(inc.Hypertension),"vnames")
names(edge.weights.Hypertension) <- sapply(names(edge.weights.Hypertension), function(x) paste( sort(strsplit(x, "\\|")[[1]]), collapse="."))

assocs_of_interest.Hypertension <- c("Brevundimonas.Porphyromonas","Butyrivibrio.Olsenella","Dialister.Streptococcus","Hyphomicrobium.Solobacterium",
                                     "Johnsonella.Prevotella","Kingella.Peptococcus","Lachnoanaerobaculum.Selenomonas","Moraxella.Shuttleworthia",
                                     "Mycoplasma.Peptococcus","Mycoplasma.Rothia")
sort(edge.weights.Hypertension[ assocs_of_interest.Hypertension ])
# ********************** #


# **************************************************************************** #

# get betaDisper and adonis values for plot:
bds.melt <- readRDS(sprintf("%s/R_objects/bds.melt.rds", p2_dir))
bds.melt.bdPs <- unique(bds.melt[ , c("Variable","bd_mean_padj","bd_mean_F","ado_mean_padj","ado_mean_R2")])
rownames(bds.melt.bdPs) <- sapply(bds.melt.bdPs$Variable, function(x) gsub("\\*","", strsplit(x, " ")[[1]][1]))

# ****************************************** #
plot_network_uniqueness <- function(score_List) {
  
  score_melt <- reshape::melt.list(lapply(score_List$var_to_var_scores, function(x) sort(unlist(x))))
  colnames(score_melt) <- c("Score","Variable")
   
  score_melt$ado_mean_R2 <- bds.melt.bdPs[ score_melt$Variable, "ado_mean_R2"]
  score_melt$Variable <- bds.melt.bdPs[ score_melt$Variable, "Variable"]
  score_melt$Variable <- sapply(score_melt$Variable, function(x) 
    gsub("detected","absent", 
         gsub("MALDI.", "",
              gsub("_", " ", 
                   gsub("Full_MALDI.Candida", "Candida detected", 
                        strsplit(x, "<br>")[[1]][1])))))
  
  # score_melt$Variable <- sapply(score_melt$Variable, function(x) gsub(" \\(n", "<br>\\(n", x))
  
  ggplot(score_melt, aes(x=reorder(Variable, Score, median), y=Score, fill=ado_mean_R2)) +
    geom_boxplot(notch = F) +
    coord_flip() +
    labs(y="Unique network score") +
    scale_fill_gradient2(low="white", high="darkred", name="mean R2") +
    theme_minimal() +
    theme(axis.title.x = element_text(size=17), axis.text.x = element_text(size=15),
          axis.title.y = element_blank(), axis.text.y = ggtext::element_markdown(size=15),
          legend.text = element_text(size=13), legend.title = element_text(size=13))
}
# ****************************************** #


# final_scores <- final_assoc_scores(all_anti_assocs.mm, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus") 
# saveRDS(final_scores, sprintf("%s/final_scores.rds", spiec_sub_dir))

final_scores <- readRDS(sprintf("%s/final_scores.rds", spiec_sub_dir))

sort(unlist(final_scores$var_scores_total))
sort(unlist(lapply(final_scores$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores$var_to_var_scores, function(x) sort(unlist(x)))

plot_network_uniqueness(final_scores)
# **************************************************************************** #

# final_scores_min80 <- final_assoc_scores(min80_anti_assocs.mm, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus")
# saveRDS(final_scores_min80, sprintf("%s/final_scores_min80.rds", spiec_sub_dir))

final_scores_min80 <- readRDS(sprintf("%s/final_scores_min80.rds", spiec_sub_dir))

sort(unlist(final_scores_min80$var_scores_total))
sort(unlist(lapply(final_scores_min80$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores_min80$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores_min80$var_to_var_scores, function(x) sort(unlist(x)))

plot_network_uniqueness(final_scores_min80)
# **************************************************************************** #


# final_scores_min95 <- final_assoc_scores(min95_anti_assocs.mm, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus")
# saveRDS(final_scores_min95, sprintf("%s/final_scores_min95.rds", spiec_sub_dir))

final_scores_min95 <- readRDS(sprintf("%s/final_scores_min95.rds", spiec_sub_dir))

sort(unlist(final_scores_min95$var_scores_total))
sort(unlist(lapply(final_scores_min95$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores_min95$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores_min95$var_to_var_scores, function(x) sort(unlist(x)))

plot_network_uniqueness(final_scores_min95)
# **************************************************************************** #


# final_scores_min100 <- final_assoc_scores(min100_anti_assocs.mm, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus")
# saveRDS(final_scores_min100, sprintf("%s/final_scores_min100.rds", spiec_sub_dir))

final_scores_min100 <- readRDS(sprintf("%s/final_scores_min100.rds", spiec_sub_dir))

sort(unlist(final_scores_min100$var_scores_total))
sort(unlist(lapply(final_scores_min100$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores_min100$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores_min100$var_to_var_scores, function(x) sort(unlist(x)))

plot_network_uniqueness(final_scores_min100)
# **************************************************************************** #




# **************************************************************************** #
# make Circos plots ####

# *************************************** #
# # first get the co-occurrence matrices for each variable 
# nc.assoc_vars <- list()
# for (v in anti_assoc_vars) {
#   print(v)
#   
#   if (v %in% c("Cystic_fibrosis","Downs_Syndrome","Celiac","Hypertension"))
#     mTab <- SLL2.meta
#   else
#     mTab <- meta.healthy
#   
#   tagl.for_network <- get_tagl_for_networks(v, tagl, all_subSamps, mTab, gloms$Genus)
#   
#   # netSamps <- unique(unlist(lapply(subSamps[[ v ]], function(i) i)))
#   # netSamps <- netSamps[ mTab[ netSamps, v] %in% varVals ]
#   # net_tagl  <- subset_samples(tagl.for_network, sample_names(tagl.for_network) %in% netSamps)
#   se.mb <- readRDS(sprintf("%s/%s/%s.se.mb.%s.rds", spiec_sub_dir, v, v, "main"))
#   # ************************ #
#   
#   nc <- symBeta(getOptBeta(se.mb))
#   colnames(nc) <- rownames(nc) <- taxa_names(tagl.for_network)
#   nc.assoc_vars[[ v ]] <- nc
# }
# saveRDS(nc.assoc_vars, sprintf("%s/nc.assoc_vars.rds", spiec_sub_dir) )

nc.assoc_vars <- readRDS(sprintf("%s/nc.assoc_vars.rds", spiec_sub_dir))

# *********************************** #
# for each variable, get list of assocs that are unique compared to at least 1 other variable

uniq_assocs <- list()
for (v in anti_assoc_vars) {
  print(v)
  
  comps <- names(min100_anti_assocs.mm)[ grepl(v, names(min100_anti_assocs.mm))]
  genPairs <- c()
  for (varPair in comps) {
    pair_vs <- strsplit(varPair, "-")[[1]]
    
    if (v == pair_vs[1])
      genPairs <- c(genPairs, names(min100_anti_assocs.mm[[ varPair ]][ min100_anti_assocs.mm[[ varPair ]] > 0 ]))
    else if (v == pair_vs[2])
      genPairs <- c(genPairs, names(min100_anti_assocs.mm[[ varPair ]][ min100_anti_assocs.mm[[ varPair ]] < 0 ]))
  }
  
  uniq_assocs[[ v ]] <- sort(unique(genPairs))
}
uniq_assocs <- lapply(uniq_assocs, function(x) gsub("unclassified.G","unclassified_G", x))
# *********************************** #

var_colors <- c("#E4BA4E","#CB6767","blue","#4b0082","#78ab78","seagreen","skyblue","#877b5d")#black")
names(var_colors) <- anti_assoc_vars

library(circlize)

# par(mfrow=c(3,3))
par(mfrow=c(2,4))

cor_circos <- list()
for (v in c("Cystic_fibrosis","Hypertension","MALDI.Yeast_detected","Full_MALDI.Candida",
            "Downs_Syndrome","Celiac","Smoker","Antibiotics")) {
# for (v in anti_assoc_vars) {
  
  cc <- as.matrix(nc.assoc_vars[[ v ]])
  cor_circos[[ v ]] <- abs(cc)
  
  for (ua in uniq_assocs[[ v ]]) {
    g1 <- gsub("unclassified_G","unclassified.G", strsplit(ua, "\\.")[[1]][1])
    g2 <- gsub("unclassified_G","unclassified.G", strsplit(ua, "\\.")[[1]][2])
    cor_circos[[ v ]][ g1, g2 ] <- cor_circos[[ v ]][ g1, g2 ] * (-5)
    cor_circos[[ v ]][ g2, g1 ] <- cor_circos[[ v ]][ g2, g1 ] * (-5)
  }
  
  # ********************* #
  # color for links
  col_fun = function(x) ifelse(x < 0, var_colors[v], "grey")
  
  # ********************* #
  # colors for grid (outer circle - go by phylum)
  grid.col <- sort( sapply(rownames(cor_circos[[ v ]]), function(x) 
    ifelse(x=="0", "0", as.character(unique(taxa_print("Genus", x)[ , "Phylum"])))) )
  
  colPal.cols <- c("#e38f87","#9abddf","#90C98A",
                   "#FFBC67","#a854a8","#ab8b67",
                   "darkblue","red","aquamarine","black",
                   "brown","darkorange3","#56423c","gray", "white")
  names(colPal.cols) <- c("Firmicutes","Bacteroidetes","Proteobacteria",
                          "Fusobacteria","Actinobacteria","Epsilonbacteraeota",
                          "Patescibacteria","Spirochaetes","Synergistetes","unclassified.P1",
                          "Verrucomicrobia","Tenericutes","Cyanobacteria","Chloroflexi", "0")
  grid.col.byPhylum <- colPal.cols[ grid.col ]
  print(unique(names(grid.col.byPhylum)))
  names(grid.col.byPhylum) <- names(grid.col)
  # ********************* #
  
  
  chordDiagram(cor_circos[[ v ]], col = col_fun, transparency = 0.25, symmetric = TRUE, 
               order = names(grid.col),#rownames(cor_circos[[ v ]]),#taxa_names(ps_neg),
               annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col.byPhylum)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    # circos.text(mean(xlim), ylim[1] + 0.95, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.25))
    # circos.axis(h = "top", labels.cex = 0.15, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  title(main = gsub("detected","absent",
                    gsub("MALDI.", "", 
                         gsub("_", " ", 
                              gsub("Full_MALDI.Candida", "Candida detected", v)))),
        cex.main = 2.5, 
        line = -1,
        # line = -2,
        col.main = var_colors[ v ])
}




# **************************************************************************** #











# final_scores.maldi_as_sub <- final_assoc_scores(all_anti_assocs.maldi_as_sub, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus")
# saveRDS(final_scores.maldi_as_sub, sprintf("%s/final_scores.maldi_as_sub.rds", spiec_sub_dir))

final_scores.maldi_as_sub <- readRDS(sprintf("%s/final_scores.maldi_as_sub.rds", spiec_sub_dir))

sort(unlist(final_scores.maldi_as_sub$var_scores_total))
sort(unlist(lapply(final_scores.maldi_as_sub$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores.maldi_as_sub$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores.maldi_as_sub$var_to_var_scores, function(x) sort(unlist(x)))
# **************************************************************************** #

# final_scores_min100.maldi_as_sub <- final_assoc_scores(min100_anti_assocs.maldi_as_sub, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus")
# saveRDS(final_scores_min100.maldi_as_sub, sprintf("%s/final_scores_min100.maldi_as_sub.rds", spiec_sub_dir))

final_scores_min100.maldi_as_sub <- readRDS(sprintf("%s/final_scores_min100.maldi_as_sub.rds", spiec_sub_dir))

sort(unlist(final_scores_min100.maldi_as_sub$var_scores_total))
sort(unlist(lapply(final_scores_min100.maldi_as_sub$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores_min100.maldi_as_sub$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores_min100.maldi_as_sub$var_to_var_scores, function(x) sort(unlist(x)))
# **************************************************************************** #





# **************************************************************************** #
# final_scores_withNegs <- final_assoc_scores(all_anti_assocs.mm, "main_main", anti_assoc_vars, tagl, all_subSamps, gloms, "Genus", 
#                                             includeNegs = T)
# saveRDS(final_scores_withNegs, sprintf("%s/final_scores_withNegs.rds", spiec_sub_dir))

final_scores_withNegs <- readRDS(sprintf("%s/final_scores_withNegs.rds", spiec_sub_dir))

sort(unlist(final_scores_withNegs$var_scores_total))
sort(unlist(lapply(final_scores_withNegs$var_to_var_scores, function(x) mean(unlist(x)))))
sort(unlist(lapply(final_scores_withNegs$var_to_var_scores, function(x) median(unlist(x)))))
lapply(final_scores_withNegs$var_to_var_scores, function(x) sort(unlist(x)))
# **************************************************************************** #

# ****************************************************************************************************************** #










# ****************************************************************************************************************** ####

# Functional analyses with t4f ####

library(themetagenomics)

# Predict the taxonomic funcional content of the whole otu table according the general taxonomy abundance:
# Use of tax4Fun instead of Picrust because of the reference database that we have been using during the process (SILVA). Picrust uses Greengenes. 

tmp <- tempdir()
download_ref(tmp, reference='silva_ko', overwrite=FALSE)



# Apply the t4f function: Given a taxonomic abundance table obtained by SILVA reference database, 
# predicts the functional content using a KO precalculated mapping table that maps the taxonomic
# abundance for a given tax_table to functional abudance content across set of functional genes. 

##reference_path=tmp --> Folder path of the silva-to-kegg mapping file 
# FUNCTIONS <- t4f(SLL2@otu_table, rows_are_taxa= TRUE, SLL2@tax_table,
#                  reference_path=tmp, type='uproc', short=TRUE,
#                  cn_normalize=TRUE, sample_normalize=TRUE, drop=TRUE)
# saveRDS(FUNCTIONS, sprintf("%s/R_objects/FUNCTIONS.rds", p2_dir))

FUNCTIONS <- readRDS(sprintf("%s/R_objects/FUNCTIONS.rds", p2_dir))


# FUNCTIONS object contains:

# fxn_table: a matrix of gene counts across samples
# fxn_meta: a list of functional metadata corresponding to fxn_table
# method_meta: a matrix of method specific metadata 


# Example of the output: 

FUNCTIONS$fxn_table[1:5,1:5]
#We can extract the general description for each of the entries: 
#head(FUNCTIONS$fxn_meta$KEGG_Pathways)   

# Each of the rows is a sample and each of the columns a functional gene (ortholog), so that we have
# the abundance content of each of the genes according the overall taxonomic abundance for each of the samples. 

# Kegg orthology (KO) - manually defined, generalised orthology groups that correspond to KEGG pathway nodes and
# BRITE hierarchy nodes in all organisms.   

# We have predicted 3184 orthologs: 
Table_functions_taxa_otu <- FUNCTIONS$fxn_table
dim(FUNCTIONS$fxn_table)
dim(Table_functions_taxa_otu)



##LINEAR MODEL TO OBTAIN SIGNIFICANTLY DIFFERENCTIALLY ABUNDANT ORTHOLOGOUS GROUPS

# Join the FUNCTIONS data with the metadata in order to asses significant differences according the metadata variables. 

#Table_functions_taxa_otu <- Table_functions_taxa_otu[, rownames(Table_functions_taxa_otu) %in% rownames(ps.meta) == TRUE]
ps.meta_path <- merge(SLL2.meta, Table_functions_taxa_otu, by=0, all=T,)
rownames(ps.meta_path) <- ps.meta_path$Row.names; ps.meta_path$Row.names <- NULL




# Now we have obtained a large data frame containing information of the metadata for each of the samples but adding the functional information. 
paths_vars <- colnames(Table_functions_taxa_otu)


# vars to run 
func_vars <- c("Smoker","Cystic_fibrosis","Downs_Syndrome","Celiac","Antibiotics","Hypertension","MALDI.Yeast_detected",
               "MALDI.Mold_detected","Full_MALDI.Candida","Age","Age_quadratic",
               "Analgesics","Anemia","Bite_nails","Circulatory_issues","Diabetes","Fluoride_toothpaste",
               "Gender","Headaches","Hypothyroidism","Kidney_issues","Lactose_intolerant","Lung_issues",
               "Mouth_wounds.binary","Vitamin_supplements")
               # "Diabetes","Fluoride_toothpaste",
               # "Wash_hands_after_bathroom","Circulatory_issues","Headaches","Lactose_intolerant","Anemia",
               # "Lung_issues","Kidney_issues","Hypothyroidism","Vitamin_supplements",
               # "MALDI.Mold_detected","Gender","Mouth_wounds.binary","Analgesics","Bite_nails","Chew_pens","Kissing_partner",
               # "seqGroup","Do_you_feel_well","Age_groups","Age_groups.TAS","Age_bins","Pets","Pets.Dogs","Pets.Cats",
               # "Pets.Both_dogs_cats","Water_type_home")

# Anova.pvals_path <- list()
# 
# for (fVar in func_vars) {
#   print(fVar)
#   Anova.pvals_path[[ fVar ]] <- run_functional_lms(fVar, SLL2, ps.meta_path, 100, paths_vars, 
#                                                    chosenControls=all_subSamps[[ fVar ]],
#                                                    silentAnova = T, useRescale=T)
#   saveRDS(Anova.pvals_path, sprintf("%s/R_objects/Anova.pvals_path.rds", p2_dir))
# }

# for (fVar in func_vars) {
#   Anova.pvals_path[[ fVar ]] <- readRDS(sprintf("%s/R_objects/Anova.pvals_path.%s.rds", p2_dir, fVar))
# }
# saveRDS(Anova.pvals_path, sprintf("%s/R_objects/Anova.pvals_path.rds", p2_dir))
Anova.pvals_path <- readRDS(sprintf("%s/R_objects/Anova.pvals_path.rds", p2_dir))



# # list of adjusted pvals, averaged over the 100 subsamples, for each variable
# Anova.pvals_path.meanAdj <- list()
# 
# for (fVar in names(Anova.pvals_path)) {
#   
#   base::print(fVar)
#   
#   covars <- colnames(Anova.pvals_path[[ fVar ]]$`1`$Anova)
#   paths <- rownames(Anova.pvals_path[[ fVar ]]$`1`$Anova)
#   
#   Anova.pvals_path.meanAdj[[ fVar ]] <- matrix(0, nrow=length(paths), ncol=length(covars))
#   rownames(Anova.pvals_path.meanAdj[[ fVar ]]) <- paths
#   colnames(Anova.pvals_path.meanAdj[[ fVar ]]) <- covars
#   
#   for (pat in paths) {
#     Anova.pvals_path.meanAdj[[ fVar ]][ pat, ] <- colMeans(apply(sapply(covars, function(x) 
#       lapply(Anova.pvals_path[[ fVar ]], function(y) y$Anova[pat, x])), 2, p.adjust, "fdr"), na.rm = T)
#   }
# }
# saveRDS(Anova.pvals_path.meanAdj, sprintf("%s/R_objects/Anova.pvals_path.meanAdj.rds", p2_dir))
Anova.pvals_path.meanAdj <- readRDS(sprintf("%s/R_objects/Anova.pvals_path.meanAdj.rds", p2_dir))


# get table of the pvals from just the variables of interest
path_heatmap <- data.frame(matrix(0, nrow=nrow(Anova.pvals_path.meanAdj$Smoker), ncol=length(names(Anova.pvals_path.meanAdj))))
rownames(path_heatmap) <- rownames(Anova.pvals_path.meanAdj$Smoker)
colnames(path_heatmap) <- names(Anova.pvals_path.meanAdj)

for (fVar in names(Anova.pvals_path.meanAdj)) {
  path_heatmap[ , fVar ] <- Anova.pvals_path.meanAdj[[ fVar ]][ , gsub("_quadratic","",fVar)]
}
# keep only those that were significant in at least 1 variable
path_heatmap.sigs <- path_heatmap[ rowSums(path_heatmap < 0.05) > 0, ]

# heatmap of pvals
d3heatmap::d3heatmap(log(path_heatmap.sigs), key=T, key.location = "tr")
d3heatmap::d3heatmap(log(path_heatmap.sigs[, c("Smoker","Cystic_fibrosis","Downs_Syndrome","MALDI.Yeast_detected",
                                               "Full_MALDI.Candida","Age","Age_quadratic")]),
                     key=T, key.location = "tr", cexCol = 1.1, cexRow = 0.01, srtCol = 25, show_grid = F)

sort(apply(path_heatmap, 2, function(x) length(x[x < 0.05])))




##Create a database of pathways and the orthologous groups involved:

#Extract list of pathways for each orthologous group that was predicted:
ortho_path_complete <- list ()
for (kegg_ortho in rownames(path_heatmap.sigs)) {
  ortho_path_complete[kegg_ortho] <- FUNCTIONS$fxn_meta$KEGG_Pathways[kegg_ortho]
}

length(ortho_path_complete) # 499


list_db_pathways_complete <- list()
for (fVar in colnames(path_heatmap.sigs)) {
  list_db_pathways_complete[[ fVar ]] <- list()
  
  for (list_ortho in rownames(path_heatmap.sigs)) {
    if (path_heatmap.sigs[ list_ortho, fVar ] >= 0.05)
      next
    #print (list_ortho)
    for (i in 1:length(FUNCTIONS$fxn_meta$KEGG_Pathways[list_ortho][[1]])) {
      path_name <- (ortho_path_complete[[list_ortho]][[i]][[3]]) #the third element is the one that contains the specific name of each of the pathways.
      #print (path_name)
      if (path_name %in% list_db_pathways_complete[[ fVar ]]) {
      } else {
        list_db_pathways_complete[[ fVar ]] <- c(list_db_pathways_complete[[ fVar ]], path_name)
      }
    }
  }
}

unlist(lapply(list_db_pathways_complete, length)) #185 different pathways from the 499 ort.
#               Smoker      Cystic_fibrosis       Downs_Syndrome               Celiac          Antibiotics         Hypertension 
#                 168                    0                  112                    0                    0                    0
# MALDI.Yeast_detected MALDI.Mold_detected   Full_MALDI.Candida                  Age        Age_quadratic 
#                   0                   0                    0                   26                   30


path_db_complete <- list()
list_proteins_complete <- list()
for (fVar in colnames(path_heatmap.sigs)) {
  
  path_db_complete[[ fVar ]] <- data.frame(matrix(ncol = 3, nrow = 0))
  #For each of the pathways of our "database"
  for (pathway_description in list_db_pathways_complete[[ fVar ]]) {
    # print (pathway_description)##
    count_times_pathway=0
    list_proteins_complete[[ fVar ]] <- list()
    
    #for each of the orthologs that we have fount that are differential expressed in our dataset 
    for (list_ortho in rownames(path_heatmap.sigs)) {
      if (path_heatmap.sigs[ list_ortho, fVar ] >= 0.05)
        next
      #for each of the pathways attributed to this ortholog 
      #print (list_ortho)
      for (i in 1:length(FUNCTIONS$fxn_meta$KEGG_Pathways[list_ortho][[1]])) {
        #take the specific name 
        path_name_complete <- (ortho_path_complete[[list_ortho]][[i]][[3]])
        
        #count in how many proteins appears the pathway and extract a list containing the identifiers for those pathways 
        if (path_name_complete == pathway_description  ) {
          count_times_pathway=count_times_pathway+1
          list_proteins_complete[[ fVar ]] <- c(list_proteins_complete[[ fVar ]], list_ortho) 
        }
      }
    }
    
    #cat ("The pathway", p, "have appeared in", c, "orthologs.\n\n")
    list_proteins_complete[[ fVar ]] <- paste(list_proteins_complete[[ fVar ]], collapse = ";")
    #print (list_proteins)
    
    path_db_complete[[ fVar ]] = rbind(path_db_complete[[ fVar ]], 
                                       data.frame(pathway_description , count_times_pathway, list_proteins_complete[[ fVar ]]))
    #afegir columna amb la llista de proteines en que apareix dita pathway 
  }
  colnames(path_db_complete[[ fVar ]]) <- c("pathway_description","count_times_pathway","list_proteins_complete")
  #Example: The pathway pentose and glucuronate interonversions have appeared in 7 differentially orthologs". 
  
}
unlist(lapply(path_db_complete, nrow))
#               Smoker      Cystic_fibrosis       Downs_Syndrome               Celiac          Antibiotics         Hypertension 
#                  168                    0                  112                    0                    0                    0 
# MALDI.Yeast_detected  MALDI.Mold_detected   Full_MALDI.Candida                  Age        Age_quadratic 
#                    0                    0                    0                   26                   30









## Stablish relation between DA. orthologous groups and pathways:
# First, for each of the DA ortholous groups, extract its list of pathways:
ortho_path <- list ()
for (kegg_ortho in rownames(path_heatmap.sigs)) {
  ortho_path[kegg_ortho] <- FUNCTIONS$fxn_meta$KEGG_Pathways[kegg_ortho]
}
# ortho_path is a list that contains all the proteins that are DA according some variable and 
# for each of them a sublist of the pathways in which are involved. 

# Then, create a database of the specific names of all the pathways that appear in our DE proteins 

# Extract the pathways in which the proteins (orthologs) are involved: For each protein that is differential 
# abundant according to a variable, extract in which pathways is involved and if the name of the pathway is not 
# included in the pathway's list (db), include it. 

list_db_pathways <- list()
for (fVar in colnames(path_heatmap.sigs)) {
  list_db_pathways[[ fVar ]] <- list()
  
  for (list_ortho in rownames(path_heatmap.sigs)) {
    if (path_heatmap.sigs[ list_ortho, fVar ] >= 0.05)
      next
    #print (list_ortho)
    for (i in 1:length(FUNCTIONS$fxn_meta$KEGG_Pathways[list_ortho][[1]])) {
      path_name <- (ortho_path[[list_ortho]][[i]][[3]]) #the third element is the one that contains the specific name of each of the pathways.
      #print (path_name)
      if (path_name %in% list_db_pathways[[ fVar ]]) {
      } else {
        list_db_pathways[[ fVar ]] <- c(list_db_pathways[[ fVar ]], path_name)
      }
    }
  }
}

unlist(lapply(list_db_pathways_complete, length)) #185 different pathways from the 499 ort.
#               Smoker      Cystic_fibrosis       Downs_Syndrome               Celiac          Antibiotics         Hypertension 
#                 168                    0                  112                    0                    0                    0
# MALDI.Yeast_detected MALDI.Mold_detected   Full_MALDI.Candida                  Age        Age_quadratic 
#                   0                   0                    0                   26                   30







# Considering all the samples: 94 of pathways (All Dx)

path_db <- list()
list_proteins <- list()
for (fVar in colnames(path_heatmap.sigs)) {
  
  ## Priorize pathways
  # For the pathways that we have saved, how many hits have these pathway in my list of DA proteins:
  path_db[[ fVar ]] <- data.frame(matrix(ncol = 3, nrow = 0))
  
  #For each of the pathways of our "database"
  for (pathway_description in list_db_pathways[[ fVar ]]) {
    # print (pathway_description)##
    count_times_pathway=0
    list_proteins[[ fVar ]] <- list()
    
    #for each of the orthologs that we have found that are differential expressed in our dataset 
    for (list_ortho in rownames(path_heatmap.sigs)) {
      if (path_heatmap.sigs[ list_ortho, fVar ] >= 0.05)
        next
      #for each of the pathways attributed to this ortholog 
      #print (list_ortho)
      for (i in 1:length(FUNCTIONS$fxn_meta$KEGG_Pathways[list_ortho][[1]])) {
        #take the specific name 
        path_name <- (ortho_path[[list_ortho]][[i]][[3]])
        
        #count in how many proteins appears the pathway and extract a list containing the identifiers for those pathways 
        if (path_name == pathway_description  ) {
          count_times_pathway=count_times_pathway+1
          list_proteins[[ fVar ]] <- c(list_proteins[[ fVar ]], list_ortho) 
        }
      }
    }
    
    #cat ("The pathway", p, "have appeared in", c, "orthologs.\n\n")
    list_proteins[[ fVar ]] <- paste(list_proteins[[ fVar ]], collapse = ";")
    #print (list_proteins)
    
    path_db[[ fVar ]] = rbind(path_db[[ fVar ]], 
                              data.frame(pathway_description , count_times_pathway, list_proteins[[ fVar ]]))
    #afegir columna amb la llista de proteines en que apareix dita pathway 
  }
  colnames(path_db[[ fVar ]]) <- c("pathway_description","count_times_pathway","list_proteins_complete")
  #Example: The pathway pentose and glucuronate interonversions have appeared in 7 differentially orthologs". 
}
unlist(lapply(path_db, nrow))
#               Smoker      Cystic_fibrosis       Downs_Syndrome               Celiac          Antibiotics         Hypertension 
#                  168                    0                  112                    0                    0                    0 
# MALDI.Yeast_detected  MALDI.Mold_detected   Full_MALDI.Candida                  Age        Age_quadratic 
#                    0                    0                    0                   26                   30



# path_db contains for each of the pathways, in how many DE proteins is a hit and the list of these proteins that are involved in the pathways. 
# Path_db [pathway_description | count_times_pathway (number of DE proteins that have attributed this pathway) | 
#          list_proteins (the name of these proteins)]


# Add a column with the total of orthologous groups that are involved in the pathway (from the 499 that we predicted)
for (fVar in colnames(path_heatmap.sigs)) {
  if (nrow(path_db[[ fVar ]]) == 0)
    next
  
  path_db[[ fVar ]]$total_OG <- 0
  for ( i in 1:nrow(path_db[[ fVar ]])) {
    path <- path_db[[ fVar ]][i,1]
    for (p in 1:nrow(path_db_complete[[ fVar ]])) {
      path_c <- path_db_complete[[ fVar ]][p,1]
      if (path == path_c) {
        path_db[[ fVar ]][i, "total_OG"] <- path_db_complete[[ fVar ]][p,2]
      }
    }
  }
}






##Relate the pathways with colorectal cancer 

# Observe if there are publications relating the pathway with colorectal cancer:
  
# To take just the pathways with more than 10 predicted OG from our taxa and 
# 10 or more % of D.A orthologous groups involved in order to priorize the pathways: 

#path_db[path_db$count_times_pathway >= 10,]
#dim(path_db_more_10)
path_db_more_10 <-  lapply(path_db, function(x)
  x[(as.numeric(x$total_OG) >= 5) & x$count_times_pathway / as.numeric(x$total_OG)*100 >= 5,])


lapply(path_db_more_10, dim)



library(easyPubMed)

# #Iniciate an empty list
# list_path_papers <- list()
# 
# for (fVar in colnames(path_heatmap.sigs)) {
#   print(fVar)
#   list_path_papers[[ fVar ]] <- list()
# 
#   if (fVar == "Smoker") {
#     searchVar <- "smoking"
#   } else if (fVar == "Downs_Syndrome") {
#     searchVar <- "down+syndrome"
#   # } else if (fVar %in% c("Age","Age_quadratic")) {
#   #   searchVar <- "aging"
#   } else {
#     next
#   }
# 
#   #for each of the pathways
#   for (name_pathway in path_db_more_10[[ fVar ]]$pathway_description) {
#     base::print (name_pathway)
#     query <- paste(name_pathway, sprintf("AND+%s[TI]", searchVar), sep="+")
#     #Construct the query to be used for pubmed taking in consideration the name of the pathway and the term colorectal cancer.
#     query <- gsub(" ", "+", query)
# 
#     #first get the api search to the pubmed database
#     query_on_pubmed <- get_pubmed_ids(query)
# 
#     #Retrive PubMed records following a search performed via the get_pubmed_ids()
#     papers <- fetch_pubmed_data(query_on_pubmed)
# 
#     titles <- custom_grep(papers, "ArticleTitle", "char")
# 
#     #For each of the pathways within the list of pathways, insert a list of articles in which are found
#     list_path_papers[[ fVar ]][name_pathway] <- list(titles) #convert the character containing all the article in a list in order to save it.
# 
#   }
# }
# saveRDS(list_path_papers, sprintf("%s/R_objects/list_path_papers.rds", p2_dir))
list_path_papers <- readRDS(sprintf("%s/R_objects/list_path_papers.rds", p2_dir))

sapply(names(list_path_papers), function(x)
  sort(unlist(lapply(list_path_papers[[ x ]], length))[ unlist(lapply(list_path_papers[[ x ]], length)) > 0], 
       decreasing = T))





# List_path_papers is a list of the pathways and the titles of the articles that relate each of the 
#  pathways and a given variable.

# Path_articles: Name of pathway | Number of articles found 

path_articles <- list()
for (fVar in colnames(path_heatmap.sigs)) {
  if ( ! is.null(list_path_papers[[ fVar ]])) {
    path_articles[[ fVar ]] <- data.frame(matrix(ncol=1, nrow = length(list_path_papers[[ fVar ]])))
    rownames(path_articles[[ fVar ]]) <- names(list_path_papers[[ fVar ]])
    colnames(path_articles[[ fVar ]]) <- "Number of articles found"
    
    for (r in rownames(path_articles[[ fVar ]])) {
      path_articles[[ fVar ]][r,"Number of articles found"] <- length(list_path_papers[[ fVar ]][[r]])
    }
  }
  
}

# Construct a dataframe containing for each of the pathways a number of articles that have found. 
# Complete this table with the number of DE orthologs and the list of these DE (to map them in the other tables). 

# To do so, combine path_db and path_articles 
# We want the rownames to be the path descriptions: 

# Bind the dataframes contaning information about 
#   the pathways (path_db [pathway_description | count_times_pathway | list_proteins] ) 
#   and path_articles (pathway description | number of articles found)

otu_functional_pathways <- list()
for (fVar in names(path_articles)) {
  if (nrow(path_articles[[ fVar ]]) == 0)
    next
  
  otu_functional_pathways[[ fVar ]] <- cbind(path_db_more_10[[ fVar ]], path_articles[[ fVar ]])
  colnames(otu_functional_pathways[[ fVar ]]) <- c("pathway_description","Number_DE_orthologs_involved", 
                                                   "List_orthologs","total_OG","Number_articles")
}

lapply(otu_functional_pathways, dim)


# Restructure the information to have:
#          Sig_Dx  P-values Pathway
# Protein1
# protein2


# First create a list containing all the DA proteins that are involved with one of these priorized pathways:
list_de_pathways <- list()

for (fVar in names(otu_functional_pathways)) {
  list_de_pathways[[ fVar ]] <- list()
  
  for (row in 1:nrow(otu_functional_pathways[[ fVar ]])) {
    DE_proteins <- as.matrix(otu_functional_pathways[[ fVar ]])[row,"List_orthologs"]
    list_DE_proteins <- as.list(strsplit(DE_proteins, ";"))
    for (element in list_DE_proteins[[1]]) {
      if (element %in% list_de_pathways[[ fVar ]]) {
        
      } else {
        list_de_pathways[[ fVar ]] <- c(list_de_pathways[[ fVar ]], element)
      }
    }
  }
  
}

lapply(list_de_pathways, length)





# Structure the dataframe:

#Being the number of rows the number of DE orthologs that have been found linked with one of the pathways
orthologs_pathways <- list()


for (fVar in names(list_de_pathways)) {
  base::print(fVar)
  orthologs_pathways[[ fVar ]] <- data.frame(matrix(ncol=3, nrow = length(list_de_pathways[[ fVar ]])))
  colnames(orthologs_pathways[[ fVar ]]) <- c("Sig_Dx", "PValues","Pathways")
  rownames(orthologs_pathways[[ fVar ]]) <- list_de_pathways[[ fVar ]]
  
  
  # Fill the dataframe:
  for (position in 1:nrow(orthologs_pathways[[ fVar ]])) {
    list_sig_dx <- list()
    #Fill the sig_dx according Tukey test
    row <- rownames(orthologs_pathways[[ fVar ]][position,0])
    #print (row)
    Sig_Dx_all <- path_heatmap.sigs[ rownames(path_heatmap.sigs)==row, names(list_de_pathways) ]#df.heatmap[rownames(df.heatmap)==row,]
    Sig_Dx <- Sig_Dx_all[colSums(!is.na(Sig_Dx_all)) > 0]
    
    
    # print (colnames(Sig_Dx))
    list_sig_dx <- list(colnames(Sig_Dx))
    string_sig_dx <- paste(list_sig_dx[[1]], collapse=",")
    orthologs_pathways[[ fVar ]][position,1] <- string_sig_dx
    
    
    #Fill the respective p-values:
    pvalues_list <- list()
    for (combination in list_sig_dx) {
      Sig_pvalue <- path_heatmap.sigs[ row, combination ]#df.heatmap[row,combination]
      
      pvalues_list <- c(pvalues_list,Sig_pvalue)
      #change to string
      string_pvalues <- paste(pvalues_list, collapse=",")
      
    }
    
    orthologs_pathways[[ fVar ]][position,2]  <- string_pvalues
    # print (string_pvalues)
    #Complete the otu_functional_pathways dataframe with the name of the pathways:
    
    pathway_list <- list()
    for (i in 1:nrow(otu_functional_pathways[[ fVar ]])) {
      #
      list_proteins_otu <- as.matrix(otu_functional_pathways[[ fVar ]])[i,3]
      list_de_proteins_otu <- as.list(strsplit(list_proteins_otu,";"))
      list_de_proteins_otu <- list_de_proteins_otu[[1]]
      
      if (row %in% list_de_proteins_otu) {
        pathway_list <- c(pathway_list, as.matrix(otu_functional_pathways[[ fVar ]])[i,1])
        string_pathways <- paste(pathway_list, collapse = ",")
        
        orthologs_pathways[[ fVar ]][position,3] <- string_pathways
      }
    }
    
  }
}


# Eliminate the rows in which we have orthologs that although they were indicated as differential abundant#
#  according Dx when performing anova, in tukey comparisons there isn't any comparison that is significant:

# orthologs_pathways <- orthologs_pathways[!orthologs_pathways$Sig_Dx == "",]





# ********************** #
# # get direction of differences for each
# orth_dirs <- list()
# for (fVar in c("Smoker","Downs_Syndrome")) {
#   print(fVar)
#   orth_dirs[[ fVar ]] <- sapply(unlist(list_de_pathways[[ fVar ]]), function(x) {
#     sapply(1:100, function(i) {
#       mod.aov <- aov(formula = as.formula(sprintf("%s ~ %s", x, fVar)), data = ps.meta_path[ all_subSamps[[ fVar ]][[ i ]], ])
#       thsd <- TukeyHSD(mod.aov)
#       thsd[[ 1 ]][ "Yes-No", "diff"]
#     })
#   })
# }
# 
# # cant use TukeyHSD for continuous value like Age
# orth_dirs[[ "Age" ]] <- sapply(unlist(list_de_pathways[[ "Age" ]]), function(x) {
#   sapply(1:100, function(i) {
#     mod.lm <- lm(formula = as.formula(sprintf("%s ~ %s", x, "Age")), data = ps.meta_path[ all_subSamps[[ "Age_bins" ]][[ i ]], ])
#     coef(mod.lm)[ 2 ]
#   })
# })
# saveRDS(orth_dirs, sprintf("%s/R_objects/orth_dirs.rds", p2_dir))
orth_dirs <- readRDS(sprintf("%s/R_objects/orth_dirs.rds", p2_dir))

# check to see if always same direction for given orths
# values of table() should all be either 0 or 100
lapply(orth_dirs, function(x) table(colSums(x < 0)))
lapply(orth_dirs, function(x) table(colSums(x > 0)))

orth_dirs.vals <- lapply(orth_dirs, colMeans)
# ********************** #

# add "up" or "down" for pathways based on orth_dirs.vals
for (fVar in names(otu_functional_pathways)) {
  counter <- 0
  counter.numart <- 0
  otu_functional_pathways[[ fVar ]]$OG_up   <- 0
  otu_functional_pathways[[ fVar ]]$OG_down <- 0
  
  for (row in rownames(otu_functional_pathways[[ fVar ]])) {
    all_orths <- strsplit(as.character(otu_functional_pathways[[ fVar ]][ row, "List_orthologs"]), ";")[[1]]
    all_orth_dirs <- orth_dirs.vals[[ fVar ]][ all_orths ]
    
    otu_functional_pathways[[ fVar ]][ row, "OG_up"]    <- length(all_orth_dirs[ all_orth_dirs > 0 ])
    otu_functional_pathways[[ fVar ]][ row, "OG_down"]  <- length(all_orth_dirs[ all_orth_dirs < 0 ])
    
    # check if all orths had same direction for a given description
    if ( ! sum(all_orth_dirs < 0) %in% c(0, length(all_orth_dirs))) {
      counter <- counter + 1; 
      # print(c(fVar, 
      #         as.character(otu_functional_pathways[[ fVar ]][ row, "Number_articles"]),
      #         as.character(otu_functional_pathways[[ fVar ]][ row, "Description"])))#row, all_orth_dirs))
      if (otu_functional_pathways[[ fVar ]][ row, "Number_articles"] == 0)
        counter.numart <- counter.numart + 1
    }
  }
  print(c(fVar, counter, nrow(otu_functional_pathways[[ fVar ]]), counter.numart))
}
# ********************** #





for (fVar in names(otu_functional_pathways)) {
  colnames(otu_functional_pathways[[ fVar ]]) <- c("Description","Number_DA_orthologs","List_orthologs",
                                                   "total_OG","Number_articles","OG_up","OG_down")
}

# PLOT

library(dplyr)
otu_functional_pathways$Smoker %>% 
  ggplot(aes(x=Number_DA_orthologs, 
             y=Description, 
             colour=Number_articles, 
             size=Number_DA_orthologs)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Number of differentially abundant orthologous groups", 
       y="Enriched pathway", 
       colour="Number of articles", 
       size="Number_DA_orthologs", 
       title = "Enriched pathways (Smoker vs non-Smoker)") +
  guides(size=FALSE) 

###
otu_functional_pathways.all <- rbind(otu_functional_pathways$Smoker, otu_functional_pathways$Downs_Syndrome)#,
                                     #otu_functional_pathways$Age)#, otu_functional_pathways$Age_quadratic)
otu_functional_pathways.all$Variable <- c(rep("Smoker", nrow(otu_functional_pathways$Smoker)),
                                          rep("Downs_Syndrome", nrow(otu_functional_pathways$Downs_Syndrome)))#,
                                          #rep("Age", nrow(otu_functional_pathways$Age)))#,
                                          # rep("Age_quadratic", nrow(otu_functional_pathways$Age_quadratic)))

otu_functional_pathways.all$Variable <- factor(otu_functional_pathways.all$Variable,
                                               levels=c("Smoker","Downs_Syndrome"))#,"Age"))#,"Age_quadratic"))
desc_by_numDAOrt <- names(sort(sapply(as.character(unique(otu_functional_pathways.all$Description)), 
                                    function(x) 
                                      max(otu_functional_pathways.all$Number_DA_orthologs[otu_functional_pathways.all$Description==x]))))
desc_by_numArt <- names(sort(sapply(as.character(unique(otu_functional_pathways.all$Description)), 
                                    function(x) 
                                      max(otu_functional_pathways.all$Number_articles[otu_functional_pathways.all$Description==x]))))
otu_functional_pathways.all$Description <- factor(otu_functional_pathways.all$Description,
                                                  levels = desc_by_numDAOrt)
                                                  # levels = desc_by_numArt)


otu_functional_pathways.all %>% 
  ggplot(aes(x=Number_DA_orthologs, 
             y=Description, 
             colour=Number_articles, 
             size=Number_DA_orthologs)) +
  geom_point() +
  expand_limits(x=0) +
  facet_wrap(~Variable, nrow=1) +
  labs(x="Number of differentially abundant orthologous groups", 
       y="Enriched pathway", 
       colour="Number of articles", 
       size="Number_DA_orthologs", 
       title = "Enriched pathways") +
  guides(size=FALSE) 




p <- ggplot(otu_functional_pathways.all, 
            aes(x=Description,#reorder(Description, -Number_articles), 
                y=Number_DA_orthologs, fill=Number_articles))  +
  geom_bar(stat="identity") + 
  facet_wrap(~Variable, nrow=1)
# Horizontal bar plot
p + coord_flip()  + 
  ggtitle("Enriched pathways according to the Diagnosis") + 
  ylab("Number of differentially abundant orthologous groups") + 
  xlab("Enriched pathway")


# ************************* #
ofp.all.withDir.up <- otu_functional_pathways.all
ofp.all.withDir.up$OG_dir_val <- ofp.all.withDir.up$OG_up
ofp.all.withDir.up$OG_dir_lab <- "Up"

ofp.all.withDir.down <- otu_functional_pathways.all
ofp.all.withDir.down$OG_dir_val <- -ofp.all.withDir.down$OG_down
ofp.all.withDir.down$OG_dir_lab <- "Down"
# ofp.all.withDir.down$Number_articles <- -ofp.all.withDir.down$Number_articles

ofp.all.withDir <- rbind(ofp.all.withDir.up, ofp.all.withDir.down)
# ************************* #
p <- ggplot(ofp.all.withDir, 
            aes(x=Description,#reorder(Description, -Number_articles), 
                y=OG_dir_val,#Number_DA_orthologs, 
                fill=Number_articles))  +
  geom_bar(stat="identity") + 
  # scale_fill_distiller(palette="Reds", direction = -1) +
  # paletteer::paletteer_c("viridis::plasma") +
  scale_fill_viridis_c(direction = -1) +
  geom_hline(yintercept = 0, color="blue") +
  theme_minimal() +
  theme(strip.text = element_text(size=13), panel.spacing = unit(2, "lines"),
        axis.text.x = element_text(size=13), axis.text.y = element_text(size=11),
        axis.title = element_text(size=15), plot.title = element_text(size=17)) +
  facet_wrap(~Variable, nrow=1 )
# Horizontal bar plot
p + coord_flip()  + 
  ggtitle("Enriched pathways by variable") + 
  ylab("Number and direction of differentially abundant KOs") + 
  xlab("Pathway") + labs(fill="# of articles")

# ****************************************************************************************************************** ####














# ****************************************************************************************************************** ####
# Compare Driver genera associations among diversity groups ####


# which associations are present in all groups for given variable:
# allCons.shannon <- lapply(gen_con_freqs.nonSubs$healthy.Diversity_group_Div.Shannon, function(x) x$all)
allCons.shannon <- lapply(gen_con_freqs.div_groups$healthy.Diversity_group_Div.Shannon, function(x) x$all)
allCons.shannon <- allCons.shannon[ unlist(lapply(allCons.shannon, function(x) length(x) > 0)) ]


# allCons.aitch <- lapply(gen_con_freqs.nonSubs$healthy.Stomatotype_Aitchison, function(x) x$all)
allCons.aitch <- lapply(gen_con_freqs.stomatotypes$healthy.Stomatotype_Aitchison, function(x) x$all)
allCons.aitch <- allCons.aitch[ unlist(lapply(allCons.aitch, function(x) length(x) > 0)) ]

# ****************************************************************************************************************** #

drivers.aitch <- get_drivers("Genus", "Aitchison", gloms_clr, meta.healthy, nTop = 10, all_or_healthy = "healthy")

# drivers.cons.aitch <- lapply(drivers.aitch, function(x) gen_con_freqs.nonSubs$healthy.Stomatotype_Aitchison[ names(x) ])
drivers.cons.aitch <- lapply(drivers.aitch, function(x) gen_con_freqs.stomatotypes$healthy.Stomatotype_Aitchison[ names(x) ])

drivers.cons.aitch$`1`$Acinetobacter
drivers.cons.aitch$`2`$Lachnoanaerobaculum

plot_primary_drivers("Genus", meta.healthy, gloms_clr, "healthy.Stomatotype_Aitchison", 
                     all_or_healthy = "healthy", nTop=3)

plot_primary_drivers("Genus", meta.healthy, gloms_clr, "healthy.Stomatotype_Weighted_Unifrac", 
                     all_or_healthy = "healthy", nTop=3)



drivers.shannon <- get_drivers("Genus", "healthy.Diversity_group_Div.Shannon", gloms_clr, meta.healthy, nTop = 10, notStomatotype = T)

# drivers.cons.shannon <- lapply(drivers.shannon, function(x) gen_con_freqs.nonSubs$healthy.Diversity_group_Div.Shannon[ names(x) ])
drivers.cons.shannon <- lapply(drivers.shannon, function(x) gen_con_freqs.div_groups$healthy.Diversity_group_Div.Shannon[ names(x) ])

drivers.cons.shannon$Low$Streptococcus
drivers.cons.shannon$Average$Prevotella
drivers.cons.shannon$High$Treponema

plot_primary_drivers("Genus", meta.healthy, gloms_clr, "healthy.Diversity_group_Div.Shannon", nTop=3, notStomatotype = T)




drivers.simpson <- get_drivers("Genus", "healthy.Diversity_group_Div.Simpson", gloms_clr, meta.healthy, nTop = 10, notStomatotype = T)

# drivers.cons.simpson <- lapply(drivers.simpson, function(x) gen_con_freqs.nonSubs$healthy.Diversity_group_Div.Simpson[ names(x) ])
drivers.cons.simpson <- lapply(drivers.simpson, function(x) gen_con_freqs.div_groups$healthy.Diversity_group_Div.Simpson[ names(x) ])

drivers.cons.simpson$Low$Streptococcus
drivers.cons.simpson$Average$Lachnoanaerobaculum
drivers.cons.simpson$High$Fusobacterium

plot_primary_drivers("Genus", meta.healthy, gloms_clr, "healthy.Diversity_group_Div.Simpson", nTop=3, notStomatotype = T)



drivers.faiths <- get_drivers("Genus", "healthy.Diversity_group_Faiths.PD", gloms_clr, meta.healthy, nTop = 10, notStomatotype = T)

# drivers.cons.faiths <- lapply(drivers.faiths, function(x) gen_con_freqs.nonSubs$healthy.Diversity_group_Faiths.PD[ names(x) ])
drivers.cons.faiths <- lapply(drivers.faiths, function(x) gen_con_freqs.div_groups$healthy.Diversity_group_Faiths.PD[ names(x) ])

drivers.cons.faiths$Low$Streptococcus
drivers.cons.faiths$Average$Leptotrichia
drivers.cons.faiths$High$Treponema

plot_primary_drivers("Genus", meta.healthy, gloms_clr, "healthy.Diversity_group_Faiths.PD", nTop=3, notStomatotype = T)


# ****************************************************************************************************************** #
# Stomatotypes for Teens only ####

healthyTeens <- rownames(SLL2.meta[ SLL2.meta$Chronic_disorder == "No" & 
                                      ! is.na(SLL2.meta$Age) &
                                      SLL2.meta$Age >= 13 & SLL2.meta$Age < 20, ])
meta.Teens <- SLL2.meta[ healthyTeens, ]
phy.Teens  <- prune_samples( healthyTeens, SLL2)

# ords.Teens <- subsampling_ordination_objects(healthyTeens, phy.Teens, dists_and_ps_Only = T, print_current_dists = T)
# saveRDS(ords.Teens, sprintf("%s/R_objects/SLL2.ords.Teens.rds", p2_dir))
ords.Teens <- readRDS(sprintf("%s/R_objects/SLL2.ords.Teens.rds", p2_dir))


for (dist_meas in c("Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard","Aitchison")) {

  glomTab <- gloms_clr

  print(dist_meas)
  clu <- get_clusters( phy.Teens, "Species", dist_meas, ords.Teens$ps.dists, glomTab,
                       subPops = T, subPop_Bdivs = ords.Teens)

  #add stomatotype to sample_data
  meta.Teens[ names(clu$cluster_full), sprintf("teen.%s_%s", clu$clus_var_name, dist_meas) ] <- clu$cluster_full#as.factor(cluster_full)
}

# *********************** #


# ****************************************************************************************************************** #
# Stomatotypes for nonBottle only ####

healthy_nonBottle <- rownames(meta.healthy[ ! is.na(meta.healthy$Water_type_home) &
                                              meta.healthy$Water_type_home != "Embotellada", ])
meta.healthy_nonBottle <- meta.healthy[ healthy_nonBottle, ]
phy.healthy_nonBottle <- phyloseq(otu_table(prune_samples(rownames(meta.healthy_nonBottle), SLL2)),
                                  sample_data(meta.healthy_nonBottle),
                                  phy_tree(prune_samples(rownames(meta.healthy_nonBottle), SLL2)),
                                  tax_table(prune_samples(rownames(meta.healthy_nonBottle), SLL2)))

nonBottle.subSamps <- lapply(all_subSamps$Water_type_home, function(i) i[i %in% rownames(meta.healthy_nonBottle)] )


# ords.healthy_nonBottle <- subsampling_ordination_objects(healthy_nonBottle, phy.healthy_nonBottle, dists_and_ps_Only = T, print_current_dists = T)
# saveRDS(ords.healthy_nonBottle, sprintf("%s/R_objects/SLL2.ords.healthy_nonBottle.rds", p2_dir))
ords.healthy_nonBottle <- readRDS(sprintf("%s/R_objects/SLL2.ords.healthy_nonBottle.rds", p2_dir))


for (dist_meas in c("Weighted_Unifrac","Unweighted_Unifrac","Bray","Jaccard","Aitchison")) {
  
  glomTab <- gloms_clr
  
  print(dist_meas)
  clu <- get_clusters( phy.healthy_nonBottle, "Species", dist_meas, ords.healthy_nonBottle$ps.dists, glomTab,
                       subPops = T, subPop_Bdivs = ords.healthy_nonBottle)
  
  #add stomatotype to sample_data
  meta.healthy_nonBottle[ names(clu$cluster_full), sprintf("nonBottle.%s_%s", clu$clus_var_name, dist_meas) ] <- clu$cluster_full#as.factor(cluster_full)
}

# *********************** #


get_drivers("Genus", "healthy.Stomatotype_Aitchison", gloms_clr, meta.healthy, nTop = 10, all_or_healthy = "healthy")
get_drivers("Genus", "healthy.Stomatotype_Weighted_Unifrac", gloms_clr, meta.healthy, nTop = 10, all_or_healthy = "healthy")



get_drivers("Genus", "nonBottle.Stomatotype_Aitchison", gloms_clr, meta.healthy_nonBottle, nTop = 10, all_or_healthy = "healthy")
get_drivers("Genus", "nonBottle.Stomatotype_Weighted_Unifrac", gloms_clr, meta.healthy_nonBottle, nTop = 10, all_or_healthy = "healthy")



get_drivers("Genus", "teen.Stomatotype_Aitchison", gloms_clr, meta.Teens, nTop = 10, all_or_healthy = "healthy")
get_drivers("Genus", "teen.Stomatotype_Weighted_Unifrac", gloms_clr, meta.Teens, nTop = 10, all_or_healthy = "healthy")



plot_primary_drivers("Genus", meta.Teens, gloms_clr, "teen.Stomatotype_Aitchison", 
                     all_or_healthy = "healthy", nTop=3)
plot_primary_drivers("Genus", meta.Teens, gloms_clr, "teen.Stomatotype_Weighted_Unifrac",
                     all_or_healthy = "healthy", nTop=3)

# ****************************************************************** #

library(scales)
library(viridis)

pcoa.nonBottle <- ape::pcoa(ords.healthy_nonBottle$Aitchison)

ade4::s.class(pcoa.nonBottle$vectors, 
              fac=as.factor(meta.healthy_nonBottle[ rownames(pcoa.nonBottle$vectors), "nonBottle.Stomatotype_Aitchison" ]),
              grid=F, clabel=1.5, 
              col=hue_pal()(length(unique(meta.healthy_nonBottle[ rownames(pcoa.nonBottle$vectors), "nonBottle.Stomatotype_Aitchison" ]))), 
              sub = sprintf("PCoA separated by %s", "nonBottle.Stomatotype_Aitchison" ))


# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####















# ****************************************************************************************************************** #

# ************************************************ #
# Proportionality of taxa with propr package ####
# tutorial here: https://github.com/tpq/propr
# and here: https://github.com/ggloor/Frontiers_2017/blob/master/Frontiers_supplement.Rmd
# ************************************************ #







# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####
























# ****************************************************************************************************************** #
# Heatmap of correlations for continuous/binary responses ####
# ****************************************************************************************************************** #


# determine which samples to look at
samps.of.int <- "All"

interest <- sample_names(SLL2)#names(cluster_vals)


# ********************************************************************************* #



cor.tables <- list()
for (comparison in c("TvsQ","taxa","questions")) {
  print(comparison)
  
  if (comparison == "questions") {
    # fill tables of p values
    cor.ps <- fill_cor_tables(comparison, NA, interest, only_cont)
    
    # adjust p-values
    cor.ps.adj <- adjust_cor_p_vals(cor.ps[[ "cor" ]], cor.ps[[ "ps" ]], comparison)
    cor.tables[[ "cor.questions" ]] <- cor.ps.adj[[ "cor.adj" ]]
    cor.tables[[ "ps.questions" ]] <- cor.ps.adj[[ "p.adj" ]]
    
    # write cor tables to files -- only print significant correlations
    cor.to.write <- cor.ps.adj[[ "cor.adj" ]]
    cor.to.write[ cor.ps.adj[[ "p.adj" ]] >= 0.05 ] <- ''
    write.csv(cor.to.write, sprintf("%s/figures/Correlations/%s/signif_correlations.%s.csv",
                                    p2_dir, samps.of.int, comparison))
    
  } else {
    for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
      print(tl)
      # fill tables of p values
      cor.ps <- fill_cor_tables(comparison, tl, interest, only_cont)
      
      # adjust p-values
      cor.ps.adj <- adjust_cor_p_vals(cor.ps[[ "cor" ]], cor.ps[[ "ps" ]], comparison)
      cor.tables[[ sprintf("cor.%s.%s",comparison,tl) ]] <- cor.ps.adj[[ "cor.adj" ]]
      cor.tables[[ sprintf("ps.%s.%s",comparison,tl) ]] <- cor.ps.adj[[ "p.adj" ]]
      
      # write cor tables to files -- only print significant correlations
      cor.to.write <- cor.ps.adj[[ "cor.adj" ]]
      cor.to.write[ cor.ps.adj[[ "p.adj" ]] >= 0.05 ] <- ''
      write.csv(cor.to.write, sprintf("%s/figures/Correlations/%s/signif_correlations.%s.%s.csv",
                                      p2_dir, samps.of.int, comparison, tl))
      
    }
  }
}





# ********************************************************************************* #



for (comparison in c("TvsQ","taxa","questions")) {
  print(comparison)
  
  if (comparison == "questions") {
    plot_cor_heatmap(cor.tables[[ "cor.questions" ]],
                     cor.tables[[ "ps.questions" ]],
                     comparison, NA)
    
  } else {
    for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
      plot_cor_heatmap(cor.tables[[ sprintf("cor.%s.%s",comparison,tl) ]],
                       cor.tables[[ sprintf("ps.%s.%s",comparison,tl) ]],
                       comparison, tl)
    }
  }
}


# ********************************************************************************* #





















# ****************************************************************************************************************** #
# Scatterplot to observe correlations ####
# ************************************* #



comparison <- "questions"
# comparison <- "TvsQ"
# comparison <- "taxa"

tl <- "Species"
# tl <- "Genus"

# dc <- get_data.cont(comparison, tl, interest)
dc <- get_data.cont(comparison, tl, sample_names(SLL2), only_cont)

# ********************************** #


n1 <- "Lactobacillus"
n2 <- "MALDI.Num_Yeast_Colonies"

dat <- data.frame(dc[,n1], dc[,n2])
colnames(dat) <- c("x","y")
dat <- dat[!is.na(dat$x) & !is.na(dat$y),]

ggplot(dat, aes(x=x, y=y)) +
  geom_point(shape=1) + 
  # geom_smooth() +
  geom_smooth(method = lm) +
  xlab(n1) + ylab(n2)

# ********************************** #


# ****************************************************************************************************************** #























# ****************************************************************************************************************** #
# Kruskal-Wallis tests ####
# ****************************************************************************************************************** #

#kruskal-wallis to test for any separations within a given survey response


# group_qs <- c( "Family_participants","Grandparent","Parent","Sibling","Grandchild","Child","Partner","Gender",
#                "Country_of_birth","Province_of_birth","City_of_birth","Country_of_birth.mother","Province_of_birt.mother",
#                "City_of_birth.mother","Country_of_birth.father","Province_of_birth.father","City_of_birth.father",
#                "Ethnicity.Caucasian","Ethnicity.Asian","Ethnicity.African","Ethnicity.Arab","Ethnicity.Gypsy",
#                "Ethnicity.Native_American","Ethnicity.No_response","Education.mother","Education.father","Education",
#                "Occupation.mother","Occupation.mother_other","Occupation.father","Occupation.father_other","Municipal_zone",
#                "Moisture_in_home","Pets","Pets.Dogs","Pets.Cats","Pets.Small_furry_animals","Pets.Birds",
#                "Pets.Reptiles_amphibians","Pets.Fish","Pets.Type.Small_furry_animals","Pets.Rabbits","Pets.Rodents",
#                # "Pets.Mammals","Pets.Horses","Smoker","Water_type_home","Braces","Mouth_piercing",
#                "Pets.Mammals","Smoker","Water_type_home","Braces","Mouth_piercing",
#                "Fluoride_toothpaste","Fluoride_supplement","Mouth_wounds","Reason_dental_visit","Reason_dental_visit_other",
#                "Chronic_disorder","Celiac","Cystic_fibrosis","Gingivitis_periodontitis","Downs_Syndrome","Eating_disorder",
#                "Other_disorder_binary","Other_disorder","Medications","Antibiotics","Analgesics","Vitamin_supplements",
#                "Other_medications_binary","Other_medications","Asthma","Wheezing","How_do_you_feel",additional_diseases,
#                "Allergy","Allergy.Mites","Allergy.Humidity","Allergy.Foods","Allergy.Pollen",
#                "Allergy.Animals","Allergy.Sun","Allergy.Medications","Allergy.Nickel","Allergy.Stings","Allergy.Latex",
#                "Allergy.Anisakis","Allergy.Seasonal","Allergy.other_binary","Allergy.other","Other_medical_matters",
#                "Bite_nails","Hair_in_mouth","Chew_pens","Wash_hands_before_eat","Wash_hands_after_bathroom",
#                "Kissing_partner","BMI_group","BMI_official",
#                # "FQ.Mutation","FQ.Lungs.Pseudomonas_aeruginosa","FQ.Lungs.Staphylococcus_aureus",
#                # "FQ.Lungs.Haemophilus_sp","FQ.Lungs.Others","FQ.Lungs.Fungi_binary","FQ.Lungs.Fungi",
#                # "FQ.Infection.Bronchopneumonia","FQ.Infection.Pneumonia","FQ.Infection.Bronchiolitis",
#                # "FQ.Transplant.Lung_transplant","FQ.Transplant.Immunosuppressant","FQ.Transplant.Other_binary",
#                # "FQ.Transplant.Other","FQ.Diet.Specialized","FQ.Diet.Reduced_carbs","FQ.Diet.More_fats","FQ.Diet.More_salt",
#                # "FQ.Diet.Isotonic_drinks","FQ.Diet.Other_binary","FQ.Diet.Other","FQ.Medications.Corticoids",
#                # "FQ.Medications.Pancreatic_enzymes","FQ.Medications.Insulin","FQ.Medications.Bronchodilators",
#                # "FQ.Medications.Immunosuppressants","FQ.Medications.Others_binary","FQ.Medications.Others",
#                # "FQ.Med_format.Oral","FQ.Med_format.Intravenous","FQ.Med_format.Inhaled","FQ.Hygiene.Wash_hands_often",
#                # "FQ.Hygiene.Dont_share_bottles","FQ.Hygiene.Avoid_contact_other_FQ","FQ.Hygiene.Other_binary",
#                # "FQ.Hygiene.Other",
#                "Consumption.Milk.Binary","Consumption.Yogurt.Binary","Consumption.Sweets.Binary",
#                "Consumption.Chewing_gum.Binary","Consumption.Nuts.Binary","Drinks.Decaf_coffee.Binary",
#                "Drinks.Coffee.Binary","Drinks.Tea.Binary","Drinks.Infusion.Binary","Drinks.Soda.Binary",
#                "Drinks.Soda_sugarless.Binary","Drinks.Soda_decaf.Binary","Drinks.RedBull.Binary",
#                "Drinks.Other_sugary_drinks.Binary","Drinks.Alcohol_cold.Binary","Drinks.Alcohol_hot.Binary",
#                "Brushing.Binary","Floss.Binary","Last_dental_visit.Binary","City","Province","Community",
#                "Diversity_group_Div.Shannon","Diversity_group_Div.Simpson","Diversity_group_Weighted_Unifrac",
#                "Diversity_group_Unweighted_Unifrac","Diversity_group_Faiths.PD","Diversity_group_Species_Richness",
#                "Diversity_group_Bray.Curtis","Diversity_group_Canberra","Stomatotype" )
# group_qs <- groupQs


# determine which samples to look at
samps.of.int <- "All"

interest <- sample_names(SLL2)#names(cluster_vals)





# ********************************************************************************* #




kw.p.tables <- list()
for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
  print(cont)
  
  # fill tables of p values
  kw.p.tables[[ sprintf('p_%s',cont) ]] <- fill_kw_p_table(SLL2_rel, cont, interest, groupQs, only_cont)
  
  # adjust p-values
  kw.p.tables[[ sprintf('p.adj_%s',cont) ]] <- adjust_kw_p_vals(kw.p.tables[[ sprintf('p_%s',cont) ]])
  
  # write tables to files
  write.csv(kw.p.tables[[ sprintf('p.adj_%s',cont) ]], 
            sprintf("%s/figures/Kruskal-Wallis/%s/signif_kruskal-wallis.%s.csv", 
                    p2_dir, samps.of.int, cont))
}








# ******************** #

interest <- sample_names(tagl)#names(cluster_vals)
# interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Celiac"]=="Yes")), generate_control_samples("Celiac"))
# interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Cystic_fibrosis"]=="Yes")), generate_control_samples("Cystic_fibrosis"))
# interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Downs_Syndrome"]=="Yes")), generate_control_samples("Downs_Syndrome"))

tlev <- "Genus"; o_rel <- gloms_rel[[ tlev ]]
group_col <- "Downs_Syndrome"
cont_col <- "Pseudomonas"

group_vs_cont_box(SLL2.meta, interest, tlev, o_rel, group_col, cont_col, groupQs, only_cont)

# tlev <- "Species"; o_rel <- SLL2_rel@otu_table
# # tlev <- "Genus"; o_rel <- gloms_rel[[ tlev ]]
# # tlev <- "Phylum"; o_rel <- gloms_rel[[ tlev ]]
# 
# interest <- sample_names(tagl)#names(cluster_vals)
# # interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Celiac"]=="Yes")), generate_control_samples("Celiac"))
# # interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Cystic_fibrosis"]=="Yes")), generate_control_samples("Cystic_fibrosis"))
# # interest <- c(sample_names(subset_samples(SLL2, SLL2@sam_data[,"Downs_Syndrome"]=="Yes")), generate_control_samples("Downs_Syndrome"))
# 
# gqs <- as.matrix(SLL2_rel@sam_data[ interest, groupQs])
# data.mix <- cbind(t(o_rel[ , interest]), SLL2_rel@sam_data[ interest, only_cont])
# data.mix <- as.data.frame(apply(data.mix, 2, factor))
# 
# group_col <- "Celiac"
# cont_col <- "Div.Observed"
# 
# kw.box <- as.data.frame( cbind(as.matrix(data.mix[, cont_col]), gqs[, group_col]) )
# colnames(kw.box) <- c("cont","group")
# kw.box[ kw.box == "No Sabe/No Contesta" ] <- NA
# kw.box <- kw.box[ ! is.na(kw.box$group), ]
# # kw.box <- kw.box[ kw.box[ , "group"] != "No Sabe/No Contesta", ]
# 
# # ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
# #                    fill=reorder(group,-as.numeric(as.character(cont)),median))) +
# ggplot(kw.box, aes(x=group, y=as.numeric(as.character(cont)), fill=group)) +
#   geom_boxplot(notch = T) + 
#   theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab(group_col) + ylab(cont_col) + scale_fill_hue(name=group_col) #+ ylim(10,50)
# # ******************** #



# # ******************** #
# # Check signif difs by anova
# 
# adf <- cbind(t(o_rel[ , interest ]), SLL2_rel@sam_data[ interest, ])
# adf[ adf == "No Sabe/No Contesta" ] <- NA
# # adf <- adf[ ! is.na(adf[, cont_col]), ] # remove NAs from 
# res.aov <- aov(formula = as.numeric(as.matrix(adf[, cont_col])) ~ as.matrix(adf[, group_col]), data = adf)
# # summary(res.aov)
# # TukeyHSD(res.aov)
# TukeyHSD(res.aov)[[1]][ TukeyHSD(res.aov)[[1]][ , "p adj"] < 0.05, ]
# # ******************** #














# ******************** #
# facet box plots for all variables that are significant with a given group_col to show which values stand out ####
table.to.check <- "Species"; o_rel <- gloms_rel[[ table.to.check ]]
# table.to.check <- "Genus"; o_rel <- gloms_rel[[ table.to.check ]]
# table.to.check <- "Phylum"; o_rel <- gloms_rel[[ table.to.check ]]
# table.to.check <- "cont_vars"


group_col <- "Community"
cont_cols <- names(as.data.frame(kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]])[
  kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]][, group_col]!="", group_col])

if (table.to.check != "cont_vars") {
  # order columns by overall abundance if looking at taxa (not cont_vars)
  cont_cols <- names(sort(rowSums(gloms_rel[[ table.to.check ]][cont_cols, ]), decreasing = T))
}

data.mix <- cbind(t(o_rel[ , interest]), SLL2_rel@sam_data[ interest, only_cont])
data.mix <- as.data.frame(apply(data.mix, 2, factor))

kw.box <- as.data.frame( cbind(gqs[, group_col], as.matrix(data.mix[, cont_cols]) ) )
colnames(kw.box)[1] <- "group"
kw.box[ kw.box == "No Sabe/No Contesta" ] <- NA
kw.box <- kw.box[ ! is.na(kw.box$group), ]

ggplot(reshape2::melt(kw.box, id.vars = "group"), aes(x=group, y=as.numeric(as.character(value)), fill=group)) +
  geom_boxplot(notch = T, outlier.alpha = 0.2) + #, outlier.shape = NA) + 
  facet_wrap(~variable, scales = "free") +
  theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  # theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
  #       axis.text.x = element_blank()) +
  xlab(group_col) + ylab(table.to.check) + scale_fill_hue(name=group_col) #+ ylim(10,50)

# ******************** #




adi <- additional_diseases[additional_diseases %in% colnames(kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]])]
adi <- adi[sapply(adi, function(group_col) length(names(as.data.frame(kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]])[
  kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]][, group_col]!="", group_col]))>1)]

for (disease in adi) {
  group_col <- disease
  cont_cols <- names(as.data.frame(kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]])[
    kw.p.tables[[ sprintf("p.adj_%s",table.to.check) ]][, group_col]!="", group_col])
  print(c(disease, length(cont_cols)))
}








# ****************************************************************************************************************** #
# *********************************************************** #
# Chi-squared test ####
# *********************************************************** #

library(graphics)
library(vcd)
library(corrplot)


interest <- sample_names(SLL2)#names(cluster_vals)


# ********************************************************************************* #  





samps.of.int <- "All"

chi.tables <- list()
for (comparison in c("TvsQ","taxa","questions")) {
  print(comparison)
  
  if (comparison == "questions") {
    # fill tables of p values
    chi.tables[[ "chi.questions" ]] <- fill_chi_table(comparison, NA, interest, groupQs)
    
    # write table to file
    write.csv(chi.tables[[ "chi.questions" ]], 
              sprintf("%s/figures/Chi-squared/%s/signif_chi-squared.%s.csv", 
                      p2_dir, samps.of.int, comparison))
  } else {
    for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
      print(sprintf('   %s',tl))
      
      # fill tables of p values
      chi.tables[[ sprintf("chi.%s.%s",comparison,tl) ]] <- fill_chi_table(comparison, tl, interest, groupQs)
      
      # write table to file
      write.csv(chi.tables[[ sprintf("chi.%s.%s",comparison,tl) ]], 
                sprintf("%s/figures/Chi-squared/%s/signif_chi-squared.%s.%s.csv", 
                        p2_dir, samps.of.int, comparison, tl))
    }
  }
}





# ********************************************************************************* #  


comparison <- "questions"
# comparison <- "TvsQ"
# comparison <- "taxa"

tl <- "Phylum"
tl <- "Species"
tl <- "Genus"

rows <- "Community"
cols <- "MALDI.Bacteria_detected"

contingency <- get_contingency_table( rows, cols, comparison, tl, interest, groupQs)
chi <- chisq.test(contingency)
assoc(contingency, shade = T, main = sprintf("Association plot"), labeling_args = list(rot_labels=45))
# ****************************************************************************************************************** #


















# ************************************************ #
# ************************************************ #

### ||  Update Chi-square -- disorder controls -- Yeast and mold detection ####

subs <- list()
subs.samps <- list()

nsubs <- 100
moreControls <- F

for (trait in c("Downs_Syndrome","Cystic_fibrosis")) {
  
  print(trait)
  subs[[ trait ]] <- list()
  subs[[ trait ]][[ "MALDI.Yeast_detected" ]] <- list()
  subs[[ trait ]][[ "MALDI.Mold_detected" ]] <- list()
  
  subs.samps[[ trait ]] <- list()
  
  for (i in 1:nsubs) {
    if (i %in% c(1, 25, 50, 75)) print(sprintf("  %s", i))
    # get appropriate subsampling of control samples for given disorder
    #  - sometimes it may be that the first set of subsampling in generate_control_samples() (by community)
    #    will not leave enough samples for a given age group in the second set of subsampling (by age)
    #    so have to use the error handling below to set control samples
    controls <- NULL
    attempt <- 0
    while( is.null(controls) ) {
      attempt <- attempt + 1
      try(
        if (moreControls) {
          controls <- generate_control_samples.larger(trait)
        } else {
          controls <- generate_control_samples(trait)
        },
        silent = TRUE
      )
    }
    
    dis.samps <- rownames(SLL2.meta[ SLL2.meta[ , trait] == "Yes", ])
    
    mTab <- SLL2.meta[ c(controls, dis.samps), ]
    subs.samps[[ trait ]][[ i ]] <- c(controls, dis.samps)
    
    subs[[ trait ]][[ "MALDI.Yeast_detected" ]][[ i ]] <- chisq.test(mTab[ , "MALDI.Yeast_detected"], mTab[ , trait ])
    subs[[ trait ]][[ "MALDI.Mold_detected" ]][[ i ]] <- chisq.test(mTab[ , "MALDI.Mold_detected"], mTab[ , trait ])
    # subs[[ trait ][[ "MALDI.Mold_detected" ]]][[ i ]] <- get_lm( c(trait, "Gender","Age","Population","seqGroup"), 
    #                                   "contVar", mTab, NULL, c("MALDI.Yeast_detected","MALDI.Mold_detected") )
    # subs[[ trait ]][[ i ]] <- get_lm( c(trait, "Gender","Age","Population","seqGroup"), 
    #                                   "contVar", mTab, NULL, c("MALDI.Yeast_detected","MALDI.Mold_detected") )
    
  }
  
}


table( p.adjust(unlist(lapply(subs$Downs_Syndrome$MALDI.Yeast_detected, function(x) x$p.value)), method = "fdr") < 0.05)
table( p.adjust(unlist(lapply(subs$Downs_Syndrome$MALDI.Mold_detected, function(x) x$p.value)), method = "fdr") < 0.05)

table( p.adjust(unlist(lapply(subs$Cystic_fibrosis$MALDI.Yeast_detected, function(x) x$p.value)), method = "fdr") < 0.05)
table( p.adjust(unlist(lapply(subs$Cystic_fibrosis$MALDI.Mold_detected, function(x) x$p.value)), method = "fdr") < 0.05)

# mean pvalues
mean( p.adjust(unlist(lapply(subs$Cystic_fibrosis$MALDI.Yeast_detected, function(x) x$p.value)), method = "fdr") )

# unadjusted pvalues
table( (unlist(lapply(subs$Downs_Syndrome$MALDI.Yeast_detected, function(x) x$p.value)) ) < 0.05)
table( (unlist(lapply(subs$Downs_Syndrome$MALDI.Mold_detected, function(x) x$p.value)) ) < 0.05)

table( (unlist(lapply(subs$Cystic_fibrosis$MALDI.Yeast_detected, function(x) x$p.value)) ) < 0.05)
table( (unlist(lapply(subs$Cystic_fibrosis$MALDI.Mold_detected, function(x) x$p.value)) ) < 0.05)

# ************************************************ #
# ************************************************ #

















# ************************************************ #
# ************************************************ #

### ||  Correlations -- disorder controls ####

# for the disorder samples only
correl.disorders.pTabs <- list()
correl.disorders.correlTabs <- list()

# for the subsamplings of controls
correl.controls.pTabs <- list()
correl.controls.correlTabs <- list()
correl.controls.sampList <- list()
correl.controls.sigProps <- list()
correl.controls.meanPs <- list()


for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  print(disorder)
  
  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
  
  # for the disorder samples only
  correl.disorders.pTabs[[ disorder ]] <- list()
  correl.disorders.correlTabs[[ disorder ]] <- list()
  
  # for the subsamplings of controls
  correl.controls.pTabs[[ disorder ]] <- list()
  correl.controls.correlTabs[[ disorder ]] <- list()
  correl.controls.sampList[[ disorder ]] <- list()
  correl.controls.sigProps[[ disorder ]] <- list()
  correl.controls.meanPs[[ disorder ]] <- list()
  
  
  for (comparison in c("taxa","questions","TvsQ")) {
    print(comparison)
    
    # if (comparison == "questions") {
    if (comparison %in% c("taxa","questions")) {
      
      # ************** #
      # for the disorder samples only
      disorder.C_P <- fill_cor_tables(comparison, NA, disorder.samples, only_cont)
      
      disorder.C_P[[ "ps.adj" ]] <- apply(disorder.C_P[[ "ps" ]], 2, p.adjust, method='bonferroni')
      #For some questions that have all 0s (occurs when doing cities alone)
      disorder.C_P[[ "ps.adj" ]][ is.na(disorder.C_P[[ "ps.adj" ]]) ] <- 1
      disorder.C_P[[ "cor" ]][ is.na(disorder.C_P[[ "cor" ]]) ] <- 0
      # must make diagonal of the p-vals table equal to 1, in order to ignore question with itself
      diag( disorder.C_P[[ "ps.adj" ]] ) <- 1
      
      correl.disorders.pTabs[[ disorder ]][[ comparison ]] <- disorder.C_P[[ "ps.adj" ]]
      correl.disorders.correlTabs[[ disorder ]][[ comparison ]] <- disorder.C_P[[ "cor" ]]
      # ************** #
      # for the subsamplings of controls
      controls.C_P <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Correlation", comparison, NA)

      correl.controls.pTabs[[ disorder ]][[ comparison ]] <- controls.C_P[[ "pTabsList" ]]
      correl.controls.correlTabs[[ disorder ]][[ comparison ]] <- controls.C_P[[ "correlTab" ]]
      correl.controls.sampList[[ disorder ]][[ comparison ]] <- controls.C_P[[ "sampList" ]]
      correl.controls.sigProps[[ disorder ]][[ comparison ]] <- controls.C_P[[ "sigProps" ]]
      correl.controls.meanPs[[ disorder ]][[ comparison ]] <- controls.C_P[[ "meanPs" ]]
      # ************** #
      
    } else {
      
      for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
        print(sprintf('   %s',tl))
        
        # ************** #
        # for the disorder samples only
        disorder.C_P <- fill_cor_tables(comparison, tl, disorder.samples, only_cont)
        
        disorder.C_P[[ "ps.adj" ]] <- apply(disorder.C_P[[ "ps" ]], 2, p.adjust, method='bonferroni')
        #For some questions that have all 0s (occurs when doing cities alone)
        disorder.C_P[[ "ps.adj" ]][ is.na(disorder.C_P[[ "ps.adj" ]]) ] <- 1
        disorder.C_P[[ "cor" ]][ is.na(disorder.C_P[[ "cor" ]]) ] <- 0
        # # must make diagonal of the p-vals table equal to 1, in order to ignore question with itself
        # if (comparison == "taxa") diag( disorder.C_P[[ "ps.adj" ]] ) <- 1
        
        correl.disorders.pTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- disorder.C_P[[ "ps.adj" ]]
        correl.disorders.correlTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- disorder.C_P[[ "cor" ]]
        # ************** #
        # for the subsamplings of controls
        controls.C_P <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Correlation", comparison, tl)

        correl.controls.pTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- controls.C_P[[ "pTabsList" ]]
        correl.controls.correlTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- controls.C_P[[ "correlTab" ]]
        correl.controls.sampList[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- controls.C_P[[ "sampList" ]]
        correl.controls.sigProps[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- controls.C_P[[ "sigProps" ]]
        correl.controls.meanPs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- controls.C_P[[ "meanPs" ]]
        # ************** #
      }
    }
    
  }

}

# ***************************** #

# # for the disorder samples only
# saveRDS(correl.disorders.pTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.disorders.pTabs.rds")
# saveRDS(correl.disorders.correlTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.disorders.correlTabs.rds")
# 
# # for the subsamplings of controls
# saveRDS(correl.controls.pTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.pTabs.rds")
# saveRDS(correl.controls.correlTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.correlTabs.rds")
# saveRDS(correl.controls.sampList,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.sampList.rds")
# saveRDS(correl.controls.sigProps,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.sigProps.rds")
# saveRDS(correl.controls.meanPs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.meanPs.rds")

# ***************************** #


# for the disorder samples only
correl.disorders.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.disorders.pTabs.rds")
correl.disorders.correlTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.disorders.correlTabs.rds")

# for the subsamplings of controls
correl.controls.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.pTabs.rds")
correl.controls.correlTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.correlTabs.rds")
correl.controls.sampList <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.sampList.rds")
correl.controls.sigProps <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.sigProps.rds")
correl.controls.meanPs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/correlations/correl.controls.meanPs.rds")

# make tables with meanPs and correlation coefficients at cells where sigProp is at least 75, else blank
# then will check where correlations differ between disorder and controls
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  
  for (comp in names(correl.controls.sigProps[[ disorder ]])) {
    
    table.of.signif <- matrix('', nrow = nrow(correl.controls.sigProps[[ disorder ]][[ comp ]]),
                             ncol = ncol(correl.controls.sigProps[[ disorder ]][[ comp ]]))
    rownames(table.of.signif) <- rownames(correl.controls.sigProps[[ disorder ]][[ comp ]])
    colnames(table.of.signif) <- colnames(correl.controls.sigProps[[ disorder ]][[ comp ]])
    
    for (ro in rownames(table.of.signif)) {
      for (co in colnames(table.of.signif)) {
        # check if given correlation is signif in at least 75% of controls and/or significant in the disorder samples
        if ( correl.controls.sigProps[[ disorder ]][[ comp ]][ ro, co ] >= 75 |
             correl.disorders.pTabs[[ disorder ]][[ comp ]][ ro, co ] < 0.05 ) {
          # table.of.signif[ ro, co ] <- correl.controls.meanPs[[ disorder ]][[ comp ]][ ro, co ]
          # table.of.signif[ ro, co ] <- "Significant"
          
          # write correlation coefficients for disorders and controls
          #   format will be:   "<disorder cor>,  <mean control cor>" with "-" in case one is not significant
          if (correl.controls.sigProps[[ disorder ]][[ comp ]][ ro, co ] >= 75 &
              correl.disorders.pTabs[[ disorder ]][[ comp ]][ ro, co ] < 0.05) {
            
            dis.cor <- round(correl.disorders.correlTabs[[ disorder ]][[ comp ]][ ro, co ], 4)
            con.cor <- round(mean(unlist(lapply(correl.controls.correlTabs[[ disorder ]][[ comp ]], 
                                                function(x) x[ro, co]))), 4)
            # check if any of these cases have opposite direction correlations between disorder and controls
            if ((dis.cor > 0 & con.cor < 0) | (dis.cor < 0 & con.cor > 0)) print(c( disorder, comp, ro, co))
            # appears there are no such cases
            
          } else if (correl.disorders.pTabs[[ disorder ]][[ comp ]][ ro, co ] < 0.05) {
            dis.cor <- round(correl.disorders.correlTabs[[ disorder ]][[ comp ]][ ro, co ], 4)
            con.cor <- "INSIGNIFICANT"
            
          } else if (correl.controls.sigProps[[ disorder ]][[ comp ]][ ro, co ] >= 75) {
            dis.cor <- "INSIGNIFICANT"
            con.cor <- round(mean(unlist(lapply(correl.controls.correlTabs[[ disorder ]][[ comp ]], 
                                                function(x) x[ro, co]))), 4)
          }
          
          table.of.signif[ ro, co ] <- sprintf("%s,  %s", dis.cor, con.cor)
          # if both disorders and controls are significant, only print if signs of coeff differ
          if (! "INSIGNIFICANT" %in% c(dis.cor, con.cor) & ( (dis.cor > 0 & con.cor > 0) | (dis.cor < 0 & con.cor < 0)) ) {
            table.of.signif[ ro, co ] <- ""
          }
          
        } 
      }
    }
    
    #at least 1 good p value within rows
    goodrows <- rownames(table.of.signif)[ sapply(rownames(table.of.signif), function(rn)
      sum(table.of.signif[ rn, ] != "") > 0) ]
    #at least 1 good p within cols
    goodcols <- colnames(table.of.signif)[ sapply(colnames(table.of.signif), function(cn)
      sum(table.of.signif[ ,cn ] != "") > 0) ]
    
    # keep only interesting rows in table to write correlations for both disorder and controls
    table.to.write <- as.matrix( table.of.signif[ goodrows, goodcols ] )
    rownames(table.to.write) <- goodrows
    colnames(table.to.write) <- goodcols
    
    dir.create(sprintf("%s/figures/Correlations/disorders/%s", p2_dir, disorder), showWarnings = F)
    write.csv(table.to.write, file = sprintf("%s/figures/Correlations/disorders/%s/%s.%s.disorder_vs_controls_cors.csv", 
                                             p2_dir, disorder, disorder, comp))
    
  }
}




# ************************************************************ #
# scatterplots of given only_cont variables in disorder samples vs controls

# ******************** #
scatterplot.disorder_controls <- function(comparison, tl1, tl2, disorder, contQs, n1, n2) {
  
  # get samples of interest
  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
  
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try(
      controls <- generate_control_samples(disorder),
      silent = TRUE
    )
  }
  if (attempt>1) cat(sprintf("Required %s attempts at subsampling\n\n", attempt))
  # ***************** #
  
  # get values for continuous variable in disoder samples and control samples
  dis.cont1 <- get_data.cont(comparison, tl1, disorder.samples, contQs)
  dis.cont2 <- get_data.cont(comparison, tl2, disorder.samples, contQs)
  con.cont1 <- get_data.cont(comparison, tl1, controls, contQs)
  con.cont2 <- get_data.cont(comparison, tl2, controls, contQs)
  
  # get combined table of Xs and Ys for scatterplot, also column for sample type to separate the plots
  dis.dat  <- data.frame(dis.cont1[,n1], dis.cont2[,n2], rep("Disorder", length(disorder.samples)))
  con.dat  <- data.frame(con.cont1[,n1], con.cont2[,n2], rep("Controls", length(controls)))
  colnames(dis.dat) <- colnames(con.dat) <- c("x","y","sampType")
  dat <- rbind(dis.dat, con.dat)
  dat <- dat[!is.na(dat$x) & !is.na(dat$y),]
  
  # plot
  ggplot(dat, aes(x=x, y=y)) +
    geom_point(shape=1) + 
    # geom_smooth() +
    geom_smooth(method = lm) +
    facet_wrap(~sampType, nrow=2, scales="fixed") +
    theme(strip.text = element_text(size=12), axis.text = element_text(size=12), axis.title = element_text(size=12)) +
    ggtitle(sprintf("%s and controls:\n%s: %s vs %s: %s", disorder, tl1, n1, tl2, n2)) +
    xlab(sprintf("%s: %s", tl1, n1)) + ylab(sprintf("%s: %s", tl2, n2))
  
}
# ******************** #

comparison <- "questions"
# comparison <- "TvsQ"
# comparison <- "taxa"

tl1 <- "Species"
# tl1 <- "Genus"
# tl1 <- "Family"
# tl1 <- "Order"
# tl1 <- "Class"
# tl1 <- "Phylum"

tl2 <- "Species"
# tl2 <- "Genus"

disorder <- "Cystic_fibrosis"
# disorder <- "Celiac"
# disorder <- "Downs_Syndrome"

n1 <- "Kingella"
n2 <- "Porphyromonas"

scatterplot.disorder_controls(comparison, tl1, tl2, disorder, only_cont, n1, n2)
# ************************************************************ #





"Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Porphyromonadaceae" "Porphyromonas"
"Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Prevotellaceae" "Prevotella"
"Bacteria" "Firmicutes" "Bacilli" "Lactobacillales" "Streptococcaceae" "Streptococcus"
"Bacteria" "Firmicutes" "Clostridia" "Clostridiales" "Family_XI" "Parvimonas"
"Bacteria" "Fusobacteria" "Fusobacteriia" "Fusobacteriales" "Fusobacteriaceae" "Fusobacterium"
"Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pasteurellales" "Pasteurellaceae" "Aggregatibacter"
"Bacteria" "Proteobacteria" "Gammaproteobacteria" "Betaproteobacteriales" "Neisseriaceae" "Kingella"
"Bacteria" "Spirochaetes" "Spirochaetia" "Spirochaetales" "Spirochaetaceae" "Treponema"

Biofilms in CF:
"Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales" "Pseudomonadaceae" "Pseudomonas"
"Bacteria" "Proteobacteria" "Gammaproteobacteria" "Xanthomonadales" "Xanthomonadaceae" "Stenotrophomonas"
Burkholderia - None
Achromobacter - None
Mycobacterium - None

Other CF pathogens:
Staphylococcus aureus "Bacteria" "Firmicutes" "Bacilli" "Bacillales" "Staphylococcaceae" "Staphylococcus"
Burkholderia cepacia - None
Haemophilus influenzae "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pasteurellales" "Pasteurellaceae" "Haemophilus"
Klebsiella pneumoniae
Stenotrophomonas maltophilia "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Xanthomonadales" "Xanthomonadaceae" "Stenotrophomonas"
Achromobacter xylosoxidans - None
nontuberculous mycobacteria - None
Enterobacteriaceae - None

Potential treatment of P aeruginosa biofilms:
"Bacteria" "Firmicutes" "Bacilli" "Lactobacillales" "Lactobacillaceae" "Lactobacillus" 



Other periodontitis:
"Bacteria" "Bacteroidetes" "Bacteroidia" "Flavobacteriales" "Weeksellaceae" "Bergeyella"
"Bacteria" "Firmicutes" "Clostridia" "Clostridiales" "Peptostreptococcaceae" "Filifactor"
"Bacteria" "Firmicutes" "Negativicutes" "Selenomonadales" "Veillonellaceae" "Dialister"
"Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pasteurellales" "Pasteurellaceae" "Actinobacillus"




# ************************************************************ #
# . . . Networks of correlations in disorders vs controls ####
# tutorial here: https://statnet.org/trac/raw-attachment/wiki/Resources/introToSNAinR_sunbelt_2012_tutorial.pdf

library(sna)
library(network)

# ****************************** #
network.disorder_controls <- function(tl, disorder, phy, gloms, min.cor=0.25) {
  
  disorder.samples <- rownames(phy@sam_data[ phy@sam_data[,disorder]=="Yes", ])
  
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try(
      controls <- generate_control_samples(disorder),
      silent = TRUE
    )
  }
  if (attempt>1) cat(sprintf("Required %s attempts at subsampling\n\n", attempt))
  
  
  # *********************** #
  net.group <- list()
  for (group in c("disorder","controls")) {
    
    net.group[[ group ]] <- list()
    
    if (group=="disorder") otu_tab <- gloms[[ tl ]][ , disorder.samples]
    if (group=="controls") otu_tab <- gloms[[ tl ]][ , controls]
    
    # must make the otu_table a as.numeric so it can be read by cor.test()
    otu.subset <- apply(otu_tab, 2, as.numeric)
    rownames(otu.subset) <- rownames(otu_tab)
    
    res.matrix <- matrix(NA, nrow=length(rownames(otu_tab)), ncol=length(rownames(otu_tab)))
    colnames(res.matrix) <- rownames(otu_tab)
    rownames(res.matrix) <- rownames(otu_tab)
    ps.matrix <- matrix(NA, nrow=length(rownames(otu_tab)), ncol=length(rownames(otu_tab)))
    colnames(ps.matrix) <- rownames(otu_tab)
    rownames(ps.matrix) <- rownames(otu_tab)
    
    for (i in rownames(otu_tab)) {
      for (j in rownames(otu_tab)) {
        correl <- cor.test(otu.subset[i,], otu.subset[j,], na.rm=T)
        res.matrix[i,j] <- correl$estimate
        ps.matrix[i,j] <- correl$p.value
      }
    }
    
    ps.adj  <- apply(ps.matrix, 2, p.adjust, method='bonferroni', n=length(rownames(otu_tab)))
    
    
    # p < 0.05 and cor_coeff > 0.6 for the graph in sna below****
    # from soil bacteria paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4649028/
    
    # First must prepare a matrix that will serve as the coordinates for vertices in a graph of the network
    # Square matrix where values are 1 if meets requirement of cor_coeff and p-val, 0 if not
    
    top20 <- names(sort(rowSums(gloms[[ tl ]]), decreasing = T))[1:20] 
    sna.gr <- matrix(0, nrow=length(top20), ncol=length(top20))
    rownames(sna.gr) <- top20
    colnames(sna.gr) <- top20
    # Make another matrix 
    
    # min.cor <- 0.25
    for (i in rownames(sna.gr)) {
      for (j in colnames(sna.gr)) {
        if (abs(res.matrix[i,j]) > min.cor & ps.adj[i,j] < 0.05) {
          if (i != j) {
            sna.gr[i,j] <- 1
          }
        }
      }
    }
    
    # Keep only those OTUs that have at least 1 edge (high enough cor coeff and signif)
    sna.gr <- sna.gr[rowSums(sna.gr)>=1,colSums(sna.gr)>=1]
    net.group[[ group ]][[ "sna.net" ]] <- as.network(sna.gr, directed=FALSE)
    
    net.group[[ group ]][[ "deg" ]] <- degree(net.group[[ group ]][[ "sna.net" ]], gmode="graph") # Indegree for MIDs
    phyla <- as.character(taxTables.both[[ tl ]][ colnames(as.sociomatrix(net.group[[ group ]][[ "sna.net" ]])), "Phylum" ])
    
    plist <- list(Actinobacteria="green", Bacteroidetes="orange", Firmicutes="cyan2", Fusobacteria="dodgerblue2",
                  Proteobacteria="red3", TM="black", Spirochaetes="grey", SR="purple")
    net.group[[ group ]][[ "phyla.col" ]] <- as.character(plist[phyla])#as.numeric(as.factor(phyla))
    net.group[[ group ]][[ "phyla" ]] <- phyla
    
    # use as.edgelist(net.group[[ group ]][[ "sna.net" ]]) to plot edge width and color based on correlation coefficient
    edges <- as.edgelist(net.group[[ group ]][[ "sna.net" ]])
    e.cor <- vector(length=nrow(edges))
    for (i in 1:nrow(edges)) {
      snd <- attr(edges, "vnames")[ edges[i, ][1] ]
      rec <- attr(edges, "vnames")[ edges[i, ][2] ]
      e.cor[i] <- res.matrix[snd,rec]
    }
    net.group[[ group ]][[ "e.colors" ]] <- ifelse(e.cor>0,"red","blue")
    net.group[[ group ]][[ "e.sizes" ]] <- ifelse(tl=="Genus", (abs(e.cor*(1/min.cor)))^3, (abs(e.cor*(1/min.cor)))^2)
    net.group[[ group ]][[ "e.cor" ]] <- e.cor
  }
  # *********************** #
  
  return(net.group)
  
  # # plots 
  # # par(mfrow=c(2,1), oma=c(2,2,2,2))
  # gplot(net.group[[ "disorder" ]][[ "sna.net" ]], gmode="graph", 
  #       vertex.col=net.group[[ "disorder" ]][[ "phyla.col" ]], 
  #       #vertex.cex=(net.group[[ "disorder" ]][[ "deg" ]])^1.5/5,
  #       displaylabels=TRUE, label.cex=1.1, label.pos=3, label.border=F, boxed.labels=T,
  #       edge.lwd=net.group[[ "disorder" ]][[ "e.sizes" ]], 
  #       edge.col=net.group[[ "disorder" ]][[ "e.colors" ]]) # mult by (1/min.cor) to ensure each is at least 1 before squaring 
  # 
  # gplot(net.group[[ "controls" ]][[ "sna.net" ]], gmode="graph", 
  #       vertex.col=net.group[[ "controls" ]][[ "phyla.col" ]], 
  #       #vertex.cex=(net.group[[ "controls" ]][[ "deg" ]])^1.5/5,
  #       displaylabels=TRUE, label.cex=1.1, label.pos=3, label.border=F, boxed.labels=T,
  #       edge.lwd=net.group[[ "controls" ]][[ "e.sizes" ]], 
  #       edge.col=net.group[[ "controls" ]][[ "e.colors" ]]) # mult by (1/min.cor) to ensure each is at least 1 before squaring 
  
}
# ****************************** #
plot_networks.disorder_controls <- function(ng, tl, disorder, mode="fruchtermanreingold") {
  # plots 
  par(mfrow=c(2,1), xpd=NA)#, oma=c(0,0,4,0), mar=c(5,4,4,2)+0.1)
  
  dis.coords <- gplot(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "sna.net" ]], gmode="graph", mode=mode, jitter = F,
                      vertex.col=ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "phyla.col" ]],
                      #vertex.cex=(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "deg" ]])^1.5/5,
                      displaylabels=TRUE, label.cex=1.1, label.pos=3, label.border=F, boxed.labels=T, pad=0, label.pad=0,
                      edge.lwd=ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.sizes" ]], 
                      edge.col=ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.colors" ]]) # mult by (1/min.cor) to ensure each is at least 1 before squaring 
  title(disorder)
  
  con.coords <- gplot(ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "sna.net" ]], gmode="graph", mode=mode,
                      vertex.col=ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "phyla.col" ]], 
                      #vertex.cex=(ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "deg" ]])^1.5/5,
                      displaylabels=TRUE, label.cex=1.1, label.pos=3, label.border=F, boxed.labels=T, pad=0, label.pad=0,
                      edge.lwd=ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.sizes" ]], 
                      edge.col=ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.colors" ]]) # mult by (1/min.cor) to ensure each is at least 1 before squaring 
  title("Controls")
  print(dis.coords)
  # print(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "sna.net" ]])
  print(con.coords)
  
  # legends
  xmax <- max(dis.coords[,"x"], con.coords[,"x"])
  xmin <- min(dis.coords[,"x"], con.coords[,"x"])
  ymax <- max(con.coords[,"y"])
  ymin <- min(dis.coords[,"y"])
  e.cor.max <- max(abs(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.cor" ]], 
                         ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.cor" ]])))
  e.cor.min <- min(abs(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.cor" ]], 
                         ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.cor" ]])))
  e.siz.max <- max(abs(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.sizes" ]], 
                         ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.sizes" ]])))
  e.siz.min <- min(abs(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "e.sizes" ]], 
                         ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "e.sizes" ]])))
  phyla <- unique(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "phyla" ]], 
                    ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "phyla" ]]))
  phyla.col <- unique(c(ng[[ disorder ]][[ tl ]][[ "disorder" ]][[ "phyla.col" ]], 
                        ng[[ disorder ]][[ tl ]][[ "controls" ]][[ "phyla.col" ]]))
  
  legend(x=xmax-2,y=ymin+2, phyla,
         fill=phyla.col, bty="n", title="Phylum", cex=1.5)
  legend(x=xmin-1.25,y=ymin+2.5, c("Signif (+) cor", "Signif (-) cor"), 
         fill=c("red","blue"), bty="n", title="Edge color", cex=1.5)
  legend(x=xmin-1.25,y=ymin+1.5, 
         c(sprintf("|cor|=%s",round(e.cor.max,2)), sprintf("|cor|=%s",round(e.cor.min,2))), 
         col=c("black","black"), bty="n", title="Edge width", cex=1.5, lwd=c(e.siz.max,e.siz.min))
}
# ****************************** #


net.groups <- list()
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  cat(disorder, "\n")
  net.groups[[ disorder ]] <- list()
  
  for (tl in c("Genus","Species")) {
    cat("   ", tl, "\n")
    net.groups[[ disorder ]][[ tl ]] <- network.disorder_controls(tl, disorder, SLL2_rel, gloms_rel, min.cor = 0.25)
  }
}
# ****************************** #

plot_networks.disorder_controls(net.groups, "Genus","Celiac", mode="circle")
# ************************************************************ #

























# ************************************************ #
# ************************************************ #

### ||  Kruskal-Wallis tests -- disorder controls ####

kw_disorder_vectors.pTabs <- list()
kw_disorder_vectors.sampList <- list()
kw_disorder_vectors.sigProps <- list()
kw_disorder_vectors.meanPs <- list()

for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  print(disorder)
  
  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
  
  kw_disorder_vectors.pTabs[[ disorder ]] <- list()
  kw_disorder_vectors.sampList[[ disorder ]] <- list()
  kw_disorder_vectors.sigProps[[ disorder ]] <- list()
  kw_disorder_vectors.meanPs[[ disorder ]] <- list()
  
  for (comparison in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
    print(comparison)
    
    kwTabs <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Kruskal-Wallis", comparison, NA)
    
    kw_disorder_vectors.pTabs[[ disorder ]][[ comparison ]] <- kwTabs[[ "pTabsList" ]]
    kw_disorder_vectors.sampList[[ disorder ]][[ comparison ]] <- kwTabs[[ "sampList" ]]
    kw_disorder_vectors.sigProps[[ disorder ]][[ comparison ]] <- kwTabs[[ "sigProps" ]]
    kw_disorder_vectors.meanPs[[ disorder ]][[ comparison ]] <- kwTabs[[ "meanPs" ]]
  }
}

# ***************************** #

# saveRDS(kw_disorder_vectors.pTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.pTabs.rds")
# saveRDS(kw_disorder_vectors.sampList,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.sampList.rds")
# saveRDS(kw_disorder_vectors.sigProps,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.sigProps.rds")
# saveRDS(kw_disorder_vectors.meanPs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.meanPs.rds")

# ***************************** #



kw_disorder_vectors.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.pTabs.rds")
kw_disorder_vectors.sampList <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.sampList.rds")
kw_disorder_vectors.sigProps <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.sigProps.rds")
kw_disorder_vectors.meanPs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorder_vectors.meanPs.rds")

# write csv files with meanPs at cells where sigProp is at least 75, else blank
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  
  for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
    
    table.to.write <- matrix('', nrow = nrow(kw_disorder_vectors.sigProps[[ disorder ]][[ cont ]]),
                             ncol = ncol(kw_disorder_vectors.sigProps[[ disorder ]][[ cont ]]))
    rownames(table.to.write) <- rownames(kw_disorder_vectors.sigProps[[ disorder ]][[ cont ]])
    colnames(table.to.write) <- colnames(kw_disorder_vectors.sigProps[[ disorder ]][[ cont ]])
    
    for (ro in rownames(table.to.write)) {
      for (co in colnames(table.to.write)) {
        if ( kw_disorder_vectors.sigProps[[ disorder ]][[ cont ]][ ro, co ] >= 75) {
          table.to.write[ ro, co ] <- kw_disorder_vectors.meanPs[[ disorder ]][[ cont ]][ ro, co ]
        } 
      }
    }
    
    #at least 1 good p value within rows
    goodrows <- rownames(table.to.write)[ sapply(rownames(table.to.write), function(rn)
      sum(table.to.write[ rn, ] != "") > 0) ]
    # #at least 1 good p within cols
    # goodcols <- colnames(table.to.write)[ sapply(colnames(table.to.write), function(cn)
    #   sum(table.to.write[ ,cn ] != "") > 0) ]
    
    # remove uninteresting rows
    table.to.write <- as.matrix( table.to.write[ goodrows,  ] )
    rownames(table.to.write) <- goodrows
    colnames(table.to.write) <- disorder
    
    dir.create(sprintf("%s/figures/Kruskal-Wallis/disorders/%s", p2_dir, disorder), showWarnings = F)
    write.csv(table.to.write, file = sprintf("%s/figures/Kruskal-Wallis/disorders/%s/%s.%s.signif_75.meanPs.csv", 
                                             p2_dir, disorder, disorder, cont))
  }
}







# ************************************************************ #
# boxplots of particular group_q vs cont data/otu abundance

# ******************** #
group_vs_cont_box.disorder_controls <- function(disorder, tlev, o_rel, cont_col, ggt="", phy=SLL2_rel) {
  
  # get samples of interest
  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
  
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try(
      # controls <- generate_control_samples(disorder),
      controls <- generate_control_samples.larger(disorder),
      silent = TRUE
    )
  }
  if (attempt>1) cat(sprintf("Required %s attempts at subsampling\n\n", attempt))
  # controls <- generate_control_samples(disorder)
  
  # fill tables of p values taking comparable control samples multiple times for each disorder
  interest <- c(disorder.samples, controls)
  # ******************** #
  
  gqs <- as.matrix(phy@sam_data[ interest, groupQs])
  data.mix <- cbind(t(o_rel[ , interest]), phy@sam_data[ interest, only_cont])
  data.mix <- as.data.frame(apply(data.mix, 2, factor))
  
  kw.box <- as.data.frame( cbind(as.matrix(data.mix[, cont_col]), gqs[, disorder]) )
  colnames(kw.box) <- c("cont","group")
  kw.box[ kw.box == "No Sabe/No Contesta" ] <- NA
  kw.box <- kw.box[ ! is.na(kw.box$group), ]
  
  # Check signif difs by anova 
  adf <- cbind(t(o_rel[ , interest ]), phy@sam_data[ interest, ])
  adf[ adf == "No Sabe/No Contesta" ] <- NA
  # adf <- adf[ ! is.na(adf[, cont_col]), ] # remove NAs from 
  res.aov <- aov(formula = as.numeric(as.matrix(adf[, cont_col])) ~ as.matrix(adf[, disorder]), data = adf)
  # summary(res.aov)
  # TukeyHSD(res.aov)
  cat("Check signif difs by anova:\n\n")
  print(TukeyHSD(res.aov)[[1]])
  
  # plot boxes
  # ggplot(kw.box, aes(x=reorder(group,-as.numeric(as.character(cont)),median), y=as.numeric(as.character(cont)), 
  #                    fill=reorder(group,-as.numeric(as.character(cont)),median))) +
  ggplot(kw.box, aes(x=group, y=as.numeric(as.character(cont)), fill=group)) +
    geom_boxplot(notch = T) + 
    ggtitle(ggt) +
    theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab(disorder) + ylab(cont_col) + scale_fill_hue(name=disorder) #+ ylim(10,50)
}
# ******************** #

disorder <- "Celiac"
# disorder <- "Cystic_fibrosis"
# disorder <- "Downs_Syndrome"

tlev <- "Genus"; o_rel <- gloms_rel[[ tlev ]]

cont_col <- "MALDI.Num_Yeast_Colonies"

group_vs_cont_box.disorder_controls(disorder, tlev, o_rel, cont_col)
# ************************************************************ #







# ************************************************************ #
# subsample reads in the Celiac samples and controls to check if low Species_Richness is due to low Gene_counts
library(GUniFrac)
library(vegan)

disorder <- "Celiac"
disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))

Celiac.rarefied.phylo <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, 
                                                      "Kruskal-Wallis", "cont_vars", NA, rarefyPhyloseq = T)
Celiac.rarefied.vegan <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, 
                                                      "Kruskal-Wallis", "cont_vars", NA, rarefyVegan = T)
Celiac.rarefied.GUF <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, 
                                                    "Kruskal-Wallis", "cont_vars", NA, rarefyGUniFrac = T)

# ************************** #
# tail(Celiac.rarefied.vegan$sigProps, n=13)
# tail(Celiac.rarefied.GUF$sigProps, n=13)

Cel.vars <- matrix(c(kw_disorder_vectors.sigProps$Celiac$cont_vars, Celiac.rarefied.phylo$sigProps,
                     Celiac.rarefied.vegan$sigProps, Celiac.rarefied.GUF$sigProps,
                     kw_disorder_vectors.meanPs$Celiac$cont_vars, Celiac.rarefied.phylo$meanPs,
                     Celiac.rarefied.vegan$meanPs, Celiac.rarefied.GUF$meanPs), ncol=8)
rownames(Cel.vars) <- rownames(Celiac.rarefied.vegan$sigProps)
colnames(Cel.vars) <- c("non_rar.sigProp","phylo.sigProp","vegan.sigProp","GUF.sigProp",
                        "non_rar.meanPs","phylo.meanPs","vegan.meanPs","GUF.meanPs")

Cel.vars.sig <- Cel.vars[sapply(rownames(Cel.vars), function(x) sum(Cel.vars[x,1:3]>=70) > 0),]
# in this we should expect the values for Div.Observed, Div.Chao1, Species_Richness and Num_OTUs to
#   all be the same bc the values in theory are the same, but when calculating alpha diversities with phyloseq
#   prunes taxa with 0 counts, which may occur more often with the phyloseq and GUF rarefy tools since they
#   just remove counts randomly ==> leads to more frequent significant diffs for Div.Observed and Div.Chao1
# ************************************************************ #


# ************************** #
# run function for boxplots randomly selecting groups of samples used in the lists above
boxes.rarefied <- function(rarefiedList, cont_col, rarefyTool) {
  # randomly select one of the list of control samples used above
  interest <- rarefiedList$sampList[[ sample(1:length(rarefiedList$sampList), 1) ]]
  # make vector of these controls with Celiac samples
  interest <- sort(c(interest, sample_names(subset_samples(SLL2, SLL2@sam_data[,"Celiac"]=="Yes"))))
  
  # then get updated Gene_counts, Species_Richness, etc
  if (rarefyTool=="Vegan") {
    SLL2.dis <- prune_samples(interest, SLL2)
    minReads <- min(colSums(SLL2.dis@otu_table))
    SLL2.dis.rar <- drarefy(SLL2.dis@otu_table, minReads)
    otu_table(SLL2.dis) <- otu_table(SLL2.dis.rar, taxa_are_rows = T)
    
  } else if (rarefyTool=="GUniFrac") {
    SLL2.dis <- prune_samples(interest, SLL2)
    otr <- Rarefy(t(SLL2.dis@otu_table))
    otu_table(SLL2.dis) <- otu_table(t(otr$otu.tab.rff), taxa_are_rows = T)
  }
  # include updated column for number of 16S gene counts per sample
  SLL2.dis@sam_data[,'Gene_counts'] <- colSums(SLL2.dis@otu_table)
  
  # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
  faith <- pd(t(as.data.frame(SLL2.dis@otu_table)), SLL2.dis@phy_tree, include.root = F)
  faith$PD[ is.na(faith$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
  SLL2.dis@sam_data[,'Faiths.PD'] <- faith$PD
  SLL2.dis@sam_data[,'Species_Richness'] <- faith$SR
  
  # finally, run function to plot boxes
  ggt <- sprintf("Rarefied with %s", rarefyTool)
  group_vs_cont_box.disorder_controls("Celiac", "Species", gloms_rel$Species, cont_col, ggt=ggt, phy=SLL2.dis)
  
}
# ************************** #

cont_col <- "Species_Richness"
# cont_col <- "Faiths.PD"
# cont_col <- "Gene_counts"

boxes.rarefied(Celiac.rarefied.vegan, cont_col, "Vegan")
boxes.rarefied(Celiac.rarefied.GUF,   cont_col, "GUniFrac")



# ************************** #
boxes.rarefied2 <- function(dis, cont_col, phy=SLL2, plotAll=FALSE) {
  
  disSamps <- rownames(phy@sam_data[phy@sam_data[,dis]=="Yes",])
  
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try(
      controls <- generate_control_samples(dis),
      silent = TRUE
    )
  }
  if (attempt>1) print(sprintf("Required %s attempts at subsampling",attempt))
  
  # fill tables of p values taking comparable control samples multiple times for each disorder
  interest <- c(disSamps, controls)
  
  # ************************************************* #
  # then get updated Gene_counts, Species_Richness, etc
  
  ###   for non_rarified samples:
  SLL2.nonRar <- prune_samples(interest, phy)
  # include updated column for number of 16S gene counts per sample
  GC.nonRar <- colSums(SLL2.nonRar@otu_table)
  # include column for number of OTUs identified
  nOTU.nonRar <- apply(SLL2.nonRar@otu_table, 2, function(x) sum( x != 0))
  # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
  faith.nonRar <- pd(t(as.data.frame(SLL2.nonRar@otu_table)), SLL2.nonRar@phy_tree, include.root = F)
  faith.nonRar$PD[ is.na(faith.nonRar$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
  # get other alpha diversity measures
  rich.nonRar <- estimate_richness(prune_taxa(taxa_sums(SLL2.nonRar)>0, SLL2.nonRar),
                                   measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  rich.nonRar[ is.nan(as.matrix(rich.nonRar)) ] <- NA
  
  ###   for phyloseq::rarefy_even_depth
  SLL2.phyRar <- prune_samples(interest, phy)
  SLL2.phyRar <- rarefy_even_depth(SLL2.phyRar, verbose = F)
  # include updated column for number of 16S gene counts per sample
  GC.phyRar <- colSums(SLL2.phyRar@otu_table)
  # include column for number of OTUs identified
  nOTU.phyRar <- apply(SLL2.phyRar@otu_table, 2, function(x) sum( x != 0))
  # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
  faith.phyRar <- pd(t(as.data.frame(SLL2.phyRar@otu_table)), SLL2.phyRar@phy_tree, include.root = F)
  faith.phyRar$PD[ is.na(faith.phyRar$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
  # get other alpha diversity measures
  SLL2.phyRar.round <- SLL2.phyRar
  rich.phyRar <- estimate_richness(prune_taxa(taxa_sums(SLL2.phyRar)>0, SLL2.phyRar),
                                   measures=c("Observed", "Chao1", "ACE","Shannon", "Simpson", "InvSimpson"))
  rich.phyRar[ is.nan(as.matrix(rich.phyRar)) ] <- NA
  
  
  ###   for vegan:
  SLL2.vegan <- prune_samples(interest, phy)
  minReads <- min(colSums(SLL2.vegan@otu_table))
  SLL2.vegan.rar <- drarefy(SLL2.vegan@otu_table, minReads)
  otu_table(SLL2.vegan) <- otu_table(SLL2.vegan.rar, taxa_are_rows = T)
  print(c(mean(colSums(SLL2.vegan@otu_table)), median(colSums(SLL2.vegan@otu_table)),
          sd(colSums(SLL2.vegan@otu_table)), min(colSums(SLL2.vegan@otu_table)),
          max(colSums(SLL2.vegan@otu_table))))
  # include updated column for number of 16S gene counts per sample
  GC.vegan <- colSums(SLL2.vegan@otu_table)
  # include column for number of OTUs identified
  nOTU.vegan <- apply(SLL2.vegan@otu_table, 2, function(x) sum( x != 0))
  # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
  faith.vegan <- pd(t(as.data.frame(SLL2.vegan@otu_table)), SLL2.vegan@phy_tree, include.root = F)
  faith.vegan$PD[ is.na(faith.vegan$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
  # get other alpha diversity measures
  SLL2.vegan.round <- SLL2.vegan
  otu_table(SLL2.vegan.round) <- round(otu_table(SLL2.vegan)) # *** because some few values will be decimals
  rich.vegan <- estimate_richness(prune_taxa(taxa_sums(SLL2.vegan)>0, SLL2.vegan),
                                   measures=c("Shannon", "Simpson", "InvSimpson"))
  rich.vegan2 <- estimate_richness(prune_taxa(taxa_sums(SLL2.vegan.round)>0, SLL2.vegan.round),
                                    measures=c("Observed", "Chao1", "ACE"))
  rich.vegan <- cbind(rich.vegan, rich.vegan2)
  rich.vegan[ is.nan(as.matrix(rich.vegan)) ] <- NA
  
  ###   for GUniFrac:
  SLL2.GUF <- prune_samples(interest, phy)
  otr <- Rarefy(t(SLL2.GUF@otu_table))
  otu_table(SLL2.GUF) <- otu_table(t(otr$otu.tab.rff), taxa_are_rows = T)
  # include updated column for number of 16S gene counts per sample
  GC.GUF <- colSums(SLL2.GUF@otu_table)
  # include column for number of OTUs identified
  nOTU.GUF <- apply(SLL2.GUF@otu_table, 2, function(x) sum( x != 0))
  # include updated column for Faith's Phylogenetic Diversity (alpha diversity that incorporates branch lengths in a tree)
  faith.GUF <- pd(t(as.data.frame(SLL2.GUF@otu_table)), SLL2.GUF@phy_tree, include.root = F)
  faith.GUF$PD[ is.na(faith.GUF$PD) ] <- 0 # for one sample that has only 1 species with non-0 count, which gives NA
  # get other alpha diversity measures
  rich.GUF <- estimate_richness(prune_taxa(taxa_sums(SLL2.GUF)>0, SLL2.GUF),
                                  measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
  rich.GUF[ is.nan(as.matrix(rich.GUF)) ] <- NA
  # ************************************************* #
  
  # Then get table for making boxplots
  cont.vals <- data.frame("Div.Shannon"=c(rich.nonRar$Shannon, rich.phyRar$Shannon, rich.vegan$Shannon, rich.GUF$Shannon),
                          "Div.Simpson"=c(rich.nonRar$Simpson, rich.phyRar$Simpson, rich.vegan$Simpson, rich.GUF$Simpson),
                          "Div.InvSimpson"=c(rich.nonRar$InvSimpson, rich.phyRar$InvSimpson, rich.vegan$InvSimpson, rich.GUF$InvSimpson),
                          "Div.Observed"=c(rich.nonRar$Observed, rich.phyRar$Observed, rich.vegan$Observed, rich.GUF$Observed),
                          "Div.Chao1"=c(rich.nonRar$Chao1, rich.phyRar$Chao1, rich.vegan$Chao1, rich.GUF$Chao1),
                          "Div.ACE"=c(rich.nonRar$ACE, rich.phyRar$ACE, rich.vegan$ACE, rich.GUF$ACE),
                          "Faiths.PD"=c(faith.nonRar$PD, faith.phyRar$PD, faith.vegan$PD, faith.GUF$PD),
                          "Species_Richness"=c(faith.nonRar$SR, faith.phyRar$SR, faith.vegan$SR, faith.GUF$SR), 
                          "Gene_counts"=c(GC.nonRar, GC.phyRar, GC.vegan, GC.GUF), 
                          "Num_OTUs"=c(nOTU.nonRar, nOTU.phyRar, nOTU.vegan, nOTU.GUF), 
                          "disorder"=ifelse(c(names(GC.nonRar), names(GC.phyRar), names(GC.vegan), names(GC.GUF)) %in% 
                                              disSamps, "Yes","No"),
                          "rarefyTool"=c(rep("Non_rarefied",length(interest)), rep("Phyloseq",length(interest)), 
                                         rep("Vegan",length(interest)), rep("GUniFrac",length(interest))))
  cont.vals <- cont.vals[, c(cont_col, "disorder", "rarefyTool")]
  cont.vals$rarefyTool <- factor(cont.vals$rarefyTool, levels = c("Non_rarefied","Phyloseq","Vegan","GUniFrac"))
  
  if (plotAll==TRUE) {
    box.toplot <- reshape2::melt(cont.vals)
    # borders <- as.numeric(as.character(sapply(1:length(unique(box.toplot$variable)), function(x) rep(x,4))))
    
    ggplot(box.toplot, aes(x=disorder, y=value, fill=disorder)) +
      geom_boxplot(notch = T) +
      facet_wrap(variable~rarefyTool, scales = "free", ncol=8) +
      # facet_grid(variable~rarefyTool, scales = "free") +
      theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.y = element_blank()) +#, panel.border = element_rect(color=borders, fill=NA)) +
      xlab(dis) #+ scale_fill_hue(name=disorder) #+ ylim(10,50)
    
  } else {
    cont.vals <- as.data.frame(apply(cont.vals, 2, factor))
    
    box.toplot <- cont.vals[,c(cont_col, "disorder", "rarefyTool")]
    colnames(box.toplot) <- c("cont","group","rarefyTool")
    
    # plot boxes
    facetScales <- ifelse(cont_col=="Gene_counts", "free", "fixed")
    
    ggplot(box.toplot, aes(x=group, y=as.numeric(as.character(cont)), fill=group)) +
      geom_boxplot(notch = T) + 
      # ggtitle(ggt) +
      facet_wrap(~rarefyTool, scales = facetScales) +
      theme(legend.position = "none", axis.text = element_text(size=12), axis.title = element_text(size=15),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      xlab(dis) + ylab(cont_col) #+ scale_fill_hue(name=disorder) #+ ylim(10,50)
  }
  
}
# ************************** #

boxes.rarefied2("Celiac","Div.Shannon")
boxes.rarefied2("Celiac","Div.Simpson")
boxes.rarefied2("Celiac","Div.InvSimpson")

boxes.rarefied2("Celiac","Species_Richness")
boxes.rarefied2("Celiac","Faiths.PD")
boxes.rarefied2("Celiac","Gene_counts")

boxes.rarefied2("Celiac", c("Div.Chao1","Div.ACE","Div.Shannon","Faiths.PD",
                            "Species_Richness","Gene_counts"), plotAll = T)
# ************************************************************ #



sapply(additional_diseases, function(x) 
  table(SLL2@sam_data[sample_names(subset_samples(SLL2, SLL2@sam_data[,"Celiac"]=="Yes")), x]))















# ************************************************ #
# ************************************************ #

### || Chi-square tests -- disorder controls ####

# ***************************** #

chi_disorder_vectors.pTabs <- list()
chi_disorder_vectors.sampList <- list()
chi_disorder_vectors.sigProps <- list()
chi_disorder_vectors.meanPs <- list()

for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  print(disorder)

  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))

  chi_disorder_vectors.pTabs[[ disorder ]] <- list()
  chi_disorder_vectors.sampList[[ disorder ]] <- list()
  chi_disorder_vectors.sigProps[[ disorder ]] <- list()
  chi_disorder_vectors.meanPs[[ disorder ]] <- list()

  for (comparison in c("TvsQ","questions")) {
    print(comparison)

    if (comparison == "questions") {

      chiTabs <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Chi-squared", comparison, NA)

      chi_disorder_vectors.pTabs[[ disorder ]][[ comparison ]] <- chiTabs[[ "pTabsList" ]]
      chi_disorder_vectors.sampList[[ disorder ]][[ comparison ]] <- chiTabs[[ "sampList" ]]
      chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]] <- chiTabs[[ "sigProps" ]]
      chi_disorder_vectors.meanPs[[ disorder ]][[ comparison ]] <- chiTabs[[ "meanPs" ]]

    } else {

      for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
        print(sprintf('   %s',tl))

        chiTabs <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Chi-squared", comparison, tl)

        chi_disorder_vectors.pTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "pTabsList" ]]
        chi_disorder_vectors.sampList[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "sampList" ]]
        chi_disorder_vectors.sigProps[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "sigProps" ]]
        chi_disorder_vectors.meanPs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "meanPs" ]]
      }
    }

  }
}

 # ***************************** #

# saveRDS(chi_disorder_vectors.pTabs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorder_vectors.pTabs.rds")
# saveRDS(chi_disorder_vectors.sampList,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorder_vectors.sampList.rds")
# saveRDS(chi_disorder_vectors.sigProps,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorder_vectors.sigProps.rds")
# saveRDS(chi_disorder_vectors.meanPs,
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorder_vectors.meanPs.rds")

# ***************************** #


chi_disorder_vectors.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorder_vectors.pTabs.rds")
chi_disorder_vectors.sampList <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorder_vectors.sampList.rds")
chi_disorder_vectors.sigProps <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorder_vectors.sigProps.rds")
chi_disorder_vectors.meanPs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorder_vectors.meanPs.rds")

# write csv files with meanPs at cells where sigProp is at least 75, else blank
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  
  for (comparison in names(chi_disorder_vectors.sigProps[[ disorder ]])) {
    
    table.to.write <- matrix('', nrow = nrow(chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]]),
                             ncol = ncol(chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]]))
    rownames(table.to.write) <- rownames(chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]])
    colnames(table.to.write) <- colnames(chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]])
    
    for (ro in rownames(table.to.write)) {
      for (co in colnames(table.to.write)) {
        if ( chi_disorder_vectors.sigProps[[ disorder ]][[ comparison ]][ ro, co ] >= 75) {
          table.to.write[ ro, co ] <- chi_disorder_vectors.meanPs[[ disorder ]][[ comparison ]][ ro, co ]
        } 
      }
    }
    
    #at least 1 good p value within rows
    goodrows <- rownames(table.to.write)[ sapply(rownames(table.to.write), function(rn)
      sum(table.to.write[ rn, ] != "") > 0) ]
    # #at least 1 good p within cols
    # goodcols <- colnames(table.to.write)[ sapply(colnames(table.to.write), function(cn)
    #   sum(table.to.write[ ,cn ] != "") > 0) ]
    
    # remove uninteresting rows
    table.to.write <- as.matrix( table.to.write[ goodrows,  ] )
    rownames(table.to.write) <- goodrows
    colnames(table.to.write) <- disorder
    
    dir.create(sprintf("%s/figures/Chi-squared/disorders/%s", p2_dir, disorder), showWarnings = F)
    write.csv(table.to.write, file = sprintf("%s/figures/Chi-squared/disorders/%s/%s.%s.signif_75.meanPs.csv", 
                                             p2_dir, disorder, disorder, comparison))
  }
}








# ********************************************************************************* #  
get_contingency_table.disorder_controls <- function(disorder, co, comparison, tl, group_qs, phy, glomTab) {
  
  # get samples of interest
  disorder.samples <- rownames(phy@sam_data[ phy@sam_data[,disorder]=="Yes", ])
  
  controls <- NULL
  attempt <- 0
  while( is.null(controls) ) {
    attempt <- attempt + 1
    try(
      controls <- generate_control_samples(disorder),
      silent = TRUE
    )
  }
  if (attempt>1) cat(sprintf("Required %s attempts at subsampling\n\n", attempt))
  # controls <- generate_control_samples(disorder)
  
  # fill tables of p values taking comparable control samples multiple times for each disorder
  interest <- c(disorder.samples, controls)
  # ******************** #
  
  # first get appropriate data
  if (comparison == "TvsQ") {
    if (disorder %in% group_qs) {
      t1 <- as.matrix( phy@sam_data )[ interest, disorder ]
      t2 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ co, ]
    } else {
      t1 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ disorder, ]
      t2 <- as.matrix( phy@sam_data )[ interest, co ]
    }
    
  } else if (comparison == "taxa") {
    t1 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ disorder, ]
    t2 <- apply( glomTab[[ tl ]][ , interest ], 2, function(x) ifelse(x==0, "Absent", "Present") )[ co, ]
    
  } else if (comparison == "questions") {
    t1 <- as.matrix( phy@sam_data )[ interest, disorder ]
    t2 <- as.matrix( phy@sam_data )[ interest, co ]
    
  }
  
  contingency <- table(t1, t2, dnn = c(disorder, co))
  
  # order diversity groups logically
  if (startsWith(disorder, "Diversity_group")) {contingency <- contingency[ c("Low","Average","High"), ]
  } else if (startsWith(co, "Diversity_group")) {contingency <- contingency[ , c("Low","Average","High") ]}
  
  contingency <- contingency[rownames(contingency) != "No Sabe/No Contesta", 
                             colnames(contingency) != "No Sabe/No Contesta"]
  
  # if (nrow(contingency) > ncol(contingency)) contingency <- t(contingency)
  
  return(contingency)
}
# ********************************************************************************* #  


comparison <- "questions"
# comparison <- "TvsQ"
# comparison <- "taxa"

tl <- "Species"
# tl <- "Genus"
# tl <- "Family"
# tl <- "Order"
# tl <- "Class"
# tl <- "Phylum"

disorder <- "Celiac"
# disorder <- "Cystic_fibrosis"
# disorder <- "Downs_Syndrome"

cols <- "MALDI.Yeast.Candida_albicans"

contingency <- get_contingency_table.disorder_controls( disorder, cols, comparison, tl, groupQs)
chi <- chisq.test(contingency)
assoc(contingency, shade = T, main = sprintf("Association plot"), labeling_args = list(rot_labels=45))
# ****************************************************************************************************************** #

















# ************************************************ #
# ************************************************ #

### ||  Ratios of taxa ==> Kruskal-Wallis tests -- disorder controls ####

# *********************************************************** #
get_taxa_ratio <- function(tl, tax1, tax2, interest, countType) {
  
  if (countType == "counts") glom <- gloms
  if (countType == "abunds") glom <- gloms_rel
  
  as.numeric(glom[[ tl ]][ tax1, interest ] / glom[[ tl ]][ tax2, interest ])
}
# *********************************************************** #
get_taxa_ratio.diffs <- function(tl, interest.groups, countType) {
  
  if (countType == "counts") glom <- gloms
  if (countType == "abunds") glom <- gloms_rel
  
  # prepare matrices for values, then fill them
  kw.p <- matrix(NA, nrow=choose(nrow(glom[[ tl ]]),2), ncol=1)
  colnames(kw.p) <- "p-value"
  rownames(kw.p) <- unique(unlist(sapply(1:(nrow(glom[[ tl ]])-1), 
                                         function(t1) sapply( (t1+1):nrow(glom[[ tl ]]), 
                                                              function(t2) sprintf("%s_%s", 
                                                                                   rownames(glom[[ tl ]])[t1], 
                                                                                   rownames(glom[[ tl ]])[t2])))))
  
  # avoid inverse of ratios already calculated (redundant) and ratios to self
  for (tax1 in 1:(nrow(glom[[ tl ]])-1)) {
    for (tax2 in (tax1+1):nrow(glom[[ tl ]])) {
      # interest.groups must be a list, where each entry is a vector of sample names for a given group
      kw.tab <- as.data.frame( unlist(interest.groups) )
      rownames(kw.tab) <- kw.tab[,1]
      colnames(kw.tab) <- "samples"
      
      kw.tab$group <- as.factor(sapply(names(interest.groups), function(x) rep(x, length(interest.groups[[ x ]]))))
      
      kw.tab$ratio <- get_taxa_ratio(tl, rownames(glom[[ tl ]])[tax1], rownames(glom[[ tl ]])[tax2], 
                                     kw.tab$samples, countType)
      
      # ******************************* #
      if (length(table(kw.tab$group, useNA = "no")) == 1) {
        kw.p[sprintf("%s_%s", rownames(glom[[ tl ]])[tax1], rownames(glom[[ tl ]])[tax2]), "p-value"] <- 1
        
      } else if (length(table(kw.tab$group, useNA = "no")) == 2 & min(table(kw.tab$group, useNA = "no")) == 1) {
        kw.p[sprintf("%s_%s", rownames(glom[[ tl ]])[tax1], rownames(glom[[ tl ]])[tax2]), "p-value"] <- 1
        
      } else if ( sum(rowSums(table(kw.tab$group, kw.tab$ratio, useNA = "no")) != 0) < 2) {
        # in the case that "all observations are in the same group"
        # length(table(gq[,j])) == 2 & 
        #          0 %in% rowSums(table(gq[,j], data.cont[,i]))) {
        # print(c(i,j))
        # where there is only a value for the cases, none for cont
        kw.p[sprintf("%s_%s", rownames(glom[[ tl ]])[tax1], rownames(glom[[ tl ]])[tax2]), "p-value"] <- 1
        
      } else {
        kw <- kruskal.test(kw.tab$ratio, kw.tab$group, na.rm=T )#na.action='na.exclude'
        kw.p[sprintf("%s_%s", rownames(glom[[ tl ]])[tax1], rownames(glom[[ tl ]])[tax2]), "p-value"] <- kw$p.value
      }
      # ******************************* #
      
    }
  }
  
  return(kw.p)
  
}
# ************************** #

# get_taxa_ratio("Genus","Genus", "Streptococcus", "Prevotella", sample_names(SLL2), countType = "abunds")
# 
# SLL2@sam_data[,"Streptococcus_Prevotella"] <- as.numeric(SLL2@otu_table["Streptococcus",] / SLL2@sam_data["Prevotella",])


ratios_disorder_vectors.pTabs <- list()
ratios_disorder_vectors.sampList <- list()
ratios_disorder_vectors.sigProps <- list()
ratios_disorder_vectors.meanPs <- list()

for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  print(disorder)
  
  disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
  
  ratios_disorder_vectors.pTabs[[ disorder ]] <- list()
  ratios_disorder_vectors.sampList[[ disorder ]] <- list()
  ratios_disorder_vectors.sigProps[[ disorder ]] <- list()
  ratios_disorder_vectors.meanPs[[ disorder ]] <- list()
  
  for (comparison in c("Phylum","Class","Order","Family","Genus","Species")) {
    print(comparison)
    
    kwTabs <- sub_to_fill_disorder_vectors(disorder, disorder.samples, 100, "Ratios", comparison, NA)
    
    ratios_disorder_vectors.pTabs[[ disorder ]][[ comparison ]] <- kwTabs[[ "pTabsList" ]]
    ratios_disorder_vectors.sampList[[ disorder ]][[ comparison ]] <- kwTabs[[ "sampList" ]]
    ratios_disorder_vectors.sigProps[[ disorder ]][[ comparison ]] <- kwTabs[[ "sigProps" ]]
    ratios_disorder_vectors.meanPs[[ disorder ]][[ comparison ]] <- kwTabs[[ "meanPs" ]]
  }
}

# ***************************** #

# saveRDS(ratios_disorder_vectors.pTabs,
#         file = sprintf("%s/R_objects/disorders_controls/ratios/ratios_disorder_vectors.pTabs.rds", p2_dir))
# saveRDS(ratios_disorder_vectors.sampList,
#         file = sprintf("%s/R_objects/disorders_controls/ratios/ratios_disorder_vectors.sampList.rds", p2_dir))
# saveRDS(ratios_disorder_vectors.sigProps,
#         file = sprintf("%s/R_objects/disorders_controls/ratios/ratios_disorder_vectors.sigProps.rds", p2_dir))
# saveRDS(ratios_disorder_vectors.meanPs,
#         file = sprintf("%s/R_objects/disorders_controls/ratios/ratios_disorder_vectors.meanPs.rds", p2_dir))

# ***************************** #







# ****************************************************************************************************************** #



























# ****************************************************************************************************************** #
# Control samples for comparison to disorders -- making full tables (unnecessary?) ####

### CELIAC
# celiac.ages <- subset_samples(SLL2, Celiac=="Yes")@sam_data$Age
# length(celiac.ages[celiac.ages<19])
# length(celiac.ages[19<=celiac.ages & celiac.ages<56])
# length(celiac.ages[celiac.ages>=56])
# 
# # get counts of samples per community
# table(subset_samples(SLL2, Celiac=="Yes")@sam_data$Community)
# 
# # prepare control samples
# celiac.controls <- sample_names(subset_samples(SLL2, Chronic_disorder=="No"))
# celiac.controls <- celiac.controls[ ! is.na(SLL2.meta[celiac.controls, "Age"])]
# # get proportionate numbers based on comunities
# celiac.controls.And <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Andalucía"]
# celiac.controls.And <- sample(celiac.controls.And, length(celiac.controls.And)*(6/51))
# celiac.controls.Ara <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Aragón"]
# celiac.controls.Ara <- sample(celiac.controls.Ara, length(celiac.controls.Ara)*(1/51))
# celiac.controls.Cat <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Cataluña"]
# celiac.controls.Cat <- sample(celiac.controls.Cat, length(celiac.controls.Cat)*(35/51))
# celiac.controls.Mad <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Comunidad de Madrid"]
# celiac.controls.Mad <- sample(celiac.controls.Mad, length(celiac.controls.Mad)*(3/51))
# celiac.controls.Val <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Comunidad Valenciana"]
# celiac.controls.Val <- sample(celiac.controls.Val, length(celiac.controls.Val)*(1/51))
# celiac.controls.Gal <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="Galicia"]
# celiac.controls.Gal <- sample(celiac.controls.Gal, length(celiac.controls.Gal)*(1/51))
# celiac.controls.PV  <- celiac.controls[SLL2.meta[celiac.controls, "Community"]=="País Vasco"]
# celiac.controls.PV  <- sample(celiac.controls.PV, length(celiac.controls.PV)*(4/51))
# celiac.controls <- sort(c(celiac.controls.And, celiac.controls.Ara, celiac.controls.Cat, 
#                           celiac.controls.Mad, celiac.controls.Val, celiac.controls.Gal, celiac.controls.PV))
# # get proportionate numbers based on ages
# celiac.controls <- sort(c(sample(celiac.controls[SLL2.meta[celiac.controls, "Age"]<19], 17),
#                           sample(celiac.controls[19<=SLL2.meta[celiac.controls, "Age"] &
#                                                    SLL2.meta[celiac.controls, "Age"]<56], 27),
#                           sample(celiac.controls[SLL2.meta[celiac.controls, "Age"]>=56], 7)))
# # verify that distribution of water types are also consistent
# table(subset_samples(SLL2, Celiac=="Yes")@sam_data$Water_type_home)
# table(SLL2.meta[celiac.controls, "Water_type_home"])
# 
# # verify that distributions of communities are also consistent
# table(subset_samples(SLL2, Celiac=="Yes")@sam_data$Community)
# table(SLL2.meta[celiac.controls, "Community"])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### Cystic_fibrosis
# cf.ages <- subset_samples(SLL2, Cystic_fibrosis=="Yes")@sam_data$Age
# length(cf.ages[cf.ages<19])
# length(cf.ages[19<=cf.ages & cf.ages<56])
# #length(cf.ages[cf.ages>=56])
# 
# # get counts of samples per community
# table(subset_samples(SLL2, Cystic_fibrosis=="Yes")@sam_data$Community)
# 
# # prepare controls 
# cf.controls <- sample_names(subset_samples(SLL2, Chronic_disorder=="No"))
# # since only 1 cf sample has filtered tap water, will remove those from controls to avoid bias
# cf.controls <- cf.controls[ ! is.na(SLL2.meta[cf.controls, "Age"]) &
#                               ! SLL2.meta[cf.controls, "Water_type_home"]=="Del Grifo (Filtrada)"]
# # get proportionate numbers based on comunities
# cf.controls.And <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Andalucía"]
# cf.controls.And <- sample(cf.controls.And, length(cf.controls.And)*(4/31))
# cf.controls.Ara <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Aragón"]
# cf.controls.Ara <- sample(cf.controls.Ara, length(cf.controls.Ara)*(1/31))
# cf.controls.Can <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Cantabria"]
# cf.controls.Can <- sample(cf.controls.Can, length(cf.controls.Can)*(11/31))
# cf.controls.Cat <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Cataluña"]
# cf.controls.Cat <- sample(cf.controls.Cat, length(cf.controls.Cat)*(6/31))
# cf.controls.Mad <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Comunidad de Madrid"]
# cf.controls.Mad <- sample(cf.controls.Mad, length(cf.controls.Mad)*(6/31))
# cf.controls.Gal <- cf.controls[SLL2.meta[cf.controls, "Community"]=="Galicia"]
# cf.controls.Gal <- sample(cf.controls.Gal, length(cf.controls.Gal)*(2/31))
# cf.controls.PV  <- cf.controls[SLL2.meta[cf.controls, "Community"]=="País Vasco"]
# cf.controls.PV  <- sample(cf.controls.PV, length(cf.controls.PV)*(1/31))
# cf.controls <- sort(c(cf.controls.And, cf.controls.Ara, cf.controls.Can, cf.controls.Cat,
#                           cf.controls.Mad, cf.controls.Gal, cf.controls.PV))
# # get proportionate numbers based on ages
# cf.controls <- sort(c(sample(cf.controls[SLL2.meta[cf.controls, "Age"]<19], 10),
#                           sample(cf.controls[19<=SLL2.meta[cf.controls, "Age"] &
#                                                    SLL2.meta[cf.controls, "Age"]<56], 21)))
#                           # sample(cf.controls[SLL2.meta[cf.controls, "Age"]>=56], 7)))
# # verify that distribution of water types are also consistent
# table(subset_samples(SLL2, Cystic_fibrosis=="Yes")@sam_data$Water_type_home)
# table(SLL2.meta[cf.controls, "Water_type_home"])
# 
# # verify that distributions of communities are also consistent
# table(subset_samples(SLL2, Cystic_fibrosis=="Yes")@sam_data$Community)
# table(SLL2.meta[cf.controls, "Community"])
# 
# 
# 
# 
# 
# 
# ### Downs_Syndrome
# down.ages <- subset_samples(SLL2, Downs_Syndrome=="Yes")@sam_data$Age
# length(down.ages[down.ages<19])
# length(down.ages[19<=down.ages & down.ages<56])
# #length(down.ages[down.ages>=56])
# 
# # get counts of samples per community
# table(subset_samples(SLL2, Downs_Syndrome=="Yes")@sam_data$Community)
# 
# 
# # prepare control samples
# down.controls <- sample_names(subset_samples(SLL2, Chronic_disorder=="No"))
# down.controls <- down.controls[ ! is.na(SLL2.meta[down.controls, "Age"])]
# # get proportionate numbers based on comunities
# down.controls.And <- down.controls[SLL2.meta[down.controls, "Community"]=="Andalucía"]
# down.controls.And <- sample(down.controls.And, length(down.controls.And)*(3/26))
# down.controls.Cat <- down.controls[SLL2.meta[down.controls, "Community"]=="Cataluña"]
# down.controls.Cat <- sample(down.controls.Cat, length(down.controls.Cat)*(13/26))
# down.controls.Gal <- down.controls[SLL2.meta[down.controls, "Community"]=="Galicia"]
# down.controls.Gal <- sample(down.controls.Gal, length(down.controls.Gal)*(1/26))
# down.controls.PV  <- down.controls[SLL2.meta[down.controls, "Community"]=="País Vasco"]
# down.controls.PV  <- sample(down.controls.PV, length(down.controls.PV)*(9/26))
# down.controls <- sort(c(down.controls.And, down.controls.Cat, down.controls.Gal, down.controls.PV))
# # get proportionate numbers based on ages
# down.controls <- sort(c(sample(down.controls[SLL2.meta[down.controls, "Age"]<19], 11),
#                       sample(down.controls[19<=SLL2.meta[down.controls, "Age"] &
#                                            SLL2.meta[down.controls, "Age"]<56], 15)))
# # sample(down.controls[SLL2.meta[down.controls, "Age"]>=56], 7)))
# # verify that distribution of water types are also consistent
# table(subset_samples(SLL2, Downs_Syndrome=="Yes")@sam_data$Water_type_home)
# table(SLL2.meta[down.controls, "Water_type_home"])
# 
# # verify that distributions of communities are also consistent
# table(subset_samples(SLL2, Downs_Syndrome=="Yes")@sam_data$Community)
# table(SLL2.meta[down.controls, "Community"])

# ************************************************ #







# ************************************************ #
sub_to_fill_disorderTabs <- function(disorder, disorder.samples, nsubs, statTest, comparison, tl) {
  
  dis.pTabsList <- list()
  dis.sampList <- list()
  
  # only need a list here when running correlations
  if (statTest == "Correlation") {
    dis.correlTab <- list()
  } else {
    dis.correlTab <- NULL
  }
  
  # then subsample for controls `nsubs` times
  for (i in 1:nsubs) {
    # get appropriate subsampling of control samples for given disorder
    #  - sometimes it may be that the first set of subsampling in generate_control_samples() (by community)
    #    will not leave enough samples for a given age group in the second set of subsampling (by age)
    #    so have to use the error handling below to set control samples
    controls <- NULL
    attempt <- 0
    while( is.null(controls) ) {
      attempt <- attempt + 1
      try(
        controls <- generate_control_samples(disorder),
        silent = TRUE
      )
    }
    if (attempt>1) print(sprintf("At i = %s, required %s attempts at subsampling", i,attempt))
    # controls <- generate_control_samples(disorder)
    
    # fill tables of p values taking comparable control samples multiple times for each disorder
    interest <- c(disorder.samples, controls)
    # keep track of controls used in each case, so can come back to check them if necessary
    dis.sampList[[ i ]] <- controls
    
    # run given statistical test to get tables of p-values (and correlation coefficients for "Correlation") and
    # multiply table by nsubs as a multiple test correction for each p-value (since running subsamplings nsubs times)
    # ********************************* #
    if (statTest == "Correlation") {
      
      if (comparison=="questions") {
        cors_and_ps <- fill_cor_tables(comparison, NA, interest, only_cont) * nsubs
      } else {
        cors_and_ps <- fill_cor_tables(comparison, tl, interest, only_cont) * nsubs
      }
      
      dis.pTabsList[[ sprintf('p.%s', i) ]] <- cors_and_ps[[ "ps" ]]
      dis.correlTab[[ sprintf('p.%s', i) ]] <- cors_and_ps[[ "cor" ]]
      
      # ********************************* #  
    } else if (statTest == "Kruskal-Wallis") {
      
      dis.pTabsList[[ sprintf('p.%s', i) ]] <- fill_kw_p_table(SLL2, comparison, interest, groupQs, only_cont) * nsubs
      
      # ********************************* #
    } else if (statTest == "Chi-squared") {
      
      if (comparison=="questions") {
        dis.pTabsList[[ sprintf('p.%s', i) ]] <- fill_chi_table(comparison, NA, interest, groupQs, keep_full = TRUE) * nsubs
      } else {
        dis.pTabsList[[ sprintf('p.%s', i) ]] <- fill_chi_table(comparison, tl, interest, groupQs, keep_full = TRUE) * nsubs
      }
      # ********************************* #
    }
    
    # in the case that all cont values were 0 (as may sometimes happen for rare species)
    #   kw test produces a p-value of NaN or NA => will change those to nsubs
    dis.pTabsList[[ sprintf('p.%s', i) ]][ is.nan(dis.pTabsList[[ sprintf('p.%s', i) ]]) ] <- nsubs
    dis.pTabsList[[ sprintf('p.%s', i) ]][ is.na(dis.pTabsList[[ sprintf('p.%s', i) ]]) ] <- nsubs
  }
  
  
  # ************ #
  # get proportion of significant p-values after correction for each KW test
  sigProps <- matrix(0, 
                     nrow = nrow(dis.pTabsList[[ sprintf('p.%s', i) ]]), 
                     ncol = ncol(dis.pTabsList[[ sprintf('p.%s', i) ]]))
  rownames(sigProps) <- rownames(dis.pTabsList[[ sprintf('p.%s', i) ]])
  colnames(sigProps) <- colnames(dis.pTabsList[[ sprintf('p.%s', i) ]])
  
  # also fill table with mean adjusted p-values for each KW test over the nsubs subsamplings
  meanPs <- sigProps
  
  for (ro in rownames(sigProps)) {
    for (co in colnames(sigProps)) {
      # get proportions
      sigProps[ ro, co ] <- 100 * sum(sapply(1:nsubs, function(i)
        dis.pTabsList[[ sprintf('p.%s', i) ]][ ro, co ] < 0.05 )) / nsubs
      # get mean p-values
      meanPs[ ro, co ] <- mean(sapply(1:nsubs, function(i) 
        dis.pTabsList[[ sprintf('p.%s', i) ]][ ro, co ] ))
    }
  }
  # ************ #
  return(list("pTabsList"=dis.pTabsList, "sigProps"=sigProps, "meanPs"=meanPs, 
              "sampList"=dis.sampList, "correlTab"=dis.correlTab))
}
# ************************************************ #






# ************************************************ #
# ************************************************ #

### ||  Kruskal-Wallis tests ###

# kw_disorders <- list()
# kw_disorders.sigProps <- list()
# kw_disorders.meanPs <- list()
# for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
#   print(disorder)
# 
#   disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
# 
#   kw_disorders[[ disorder ]] <- list()
#   kw_disorders.sigProps[[ disorder ]] <- list()
#   kw_disorders.meanPs[[ disorder ]] <- list()
# 
#   for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
#     print(cont)
# 
#     kw_disorders[[ disorder ]][[ cont ]] <- list()
#     kw_disorders.sigProps[[ disorder ]][[ cont ]] <- list()
#     kw_disorders.meanPs[[ disorder ]][[ cont ]] <- list()
# 
#     nsubs <- 100
#     for (i in 1:nsubs) {
#       # get appropriate subsampling of control samples for given disorder
#       #  - sometimes it may be that the first set of subsampling in generate_control_samples() (by community)
#       #    will not leave enough samples for a given age group in the second set of subsampling (by age)
#       #    so have to use the error handling below to set control samples
#       controls <- NULL
#       attempt <- 0
#       while( is.null(controls) ) {
#         attempt <- attempt + 1
#         try(
#           controls <- generate_control_samples(disorder),
#           silent = TRUE
#         )
#       }
#       if (attempt>1) print(sprintf("At i = %s, required %s attempts at subsampling", i,attempt))
#       # controls <- generate_control_samples(disorder)
# 
#       # fill tables of p values taking comparable control samples multiple times for each disorder
#       interest <- c(disorder.samples, controls)
#       # multiply table by nsubs as a multiple test correction for each p-value (since running subsamplings nsubs times)
#       kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]] <- fill_kw_p_table(cont, interest) * nsubs
#       
#       # in the case that all cont values were 0 (as may sometimes happen for rare species)
#       #   kw test produces a p-value of NaN or NA => will change those to nsubs
#       kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]][ 
#         is.nan(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]])] <- nsubs
#       kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]][ 
#         is.na(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]])] <- nsubs
#     }
# 
# 
#     # ************ #
#     # get proportion of significant p-values after correction for each KW test
#     kw_disorders.sigProps[[ disorder ]][[ cont ]] <- matrix(0,
#                                                             nrow = nrow(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]]),
#                                                             ncol = ncol(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]]))
#     rownames(kw_disorders.sigProps[[ disorder ]][[ cont ]]) <- rownames(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]])
#     colnames(kw_disorders.sigProps[[ disorder ]][[ cont ]]) <- colnames(kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]])
# 
#     # also fill table with mean adjusted p-values for each KW test over the nsubs subsamplings
#     kw_disorders.meanPs[[ disorder ]][[ cont ]] <- kw_disorders.sigProps[[ disorder ]][[ cont ]]
# 
#     for (ro in rownames(kw_disorders.sigProps[[ disorder ]][[ cont ]])) {
#       for (co in colnames(kw_disorders.sigProps[[ disorder ]][[ cont ]])) {
#         # get proportions
#         kw_disorders.sigProps[[ disorder ]][[ cont ]][ ro, co ] <- 100 * sum(sapply(1:nsubs, function(i)
#           kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]][ ro, co ] < 0.05 )) / nsubs
#         # get mean p-values
#         kw_disorders.meanPs[[ disorder ]][[ cont ]][ ro, co ] <- mean(sapply(1:nsubs, function(i)
#           kw_disorders[[ disorder ]][[ cont ]][[ sprintf('p.%s', i) ]][ ro, co ] ))
#       }
#     }
#     # ************ #
# 
#   }
# 
# }
#




# 
# kw_disorders.pTabs <- list()
# kw_disorders.sampList <- list()
# kw_disorders.sigProps <- list()
# kw_disorders.meanPs <- list()
# 
# for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
#   print(disorder)
# 
#   disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
# 
#   kw_disorders.pTabs[[ disorder ]] <- list()
#   kw_disorders.sampList[[ disorder ]] <- list()
#   kw_disorders.sigProps[[ disorder ]] <- list()
#   kw_disorders.meanPs[[ disorder ]] <- list()
#   
#   for (comparison in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
#     print(comparison)
#     
#     kwTabs <- sub_to_fill_disorderTabs(disorder, disorder.samples, 100, "Kruskal-Wallis", comparison, NA)
#     
#     kw_disorders.pTabs[[ disorder ]][[ comparison ]] <- kwTabs[[ "pTabsList" ]]
#     kw_disorders.sampList[[ disorder ]][[ comparison ]] <- kwTabs[[ "sampList" ]]
#     kw_disorders.sigProps[[ disorder ]][[ comparison ]] <- kwTabs[[ "sigProps" ]]
#     kw_disorders.meanPs[[ disorder ]][[ comparison ]] <- kwTabs[[ "meanPs" ]]
#   }
# }
# 
# # ***************************** #
# 
# saveRDS(kw_disorders.pTabs, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorders.pTabs.rds")
# saveRDS(kw_disorders.sampList, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorders.sampList.rds")
# saveRDS(kw_disorders.sigProps, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorders.sigProps.rds")
# saveRDS(kw_disorders.meanPs, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/kruskal-wallis/kw_disorders.meanPs.rds")

# ***************************** #

kw_disorders.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/kw_disorders.pTabs.rds")
kw_disorders.sampList <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/kw_disorders.sampList.rds")
kw_disorders.sigProps <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/kw_disorders.sigProps.rds")
kw_disorders.meanPs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/kw_disorders.meanPs.rds")

# write csv files with meanPs at cells where sigProp is at least 75, else blank
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {

  for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
    
    table.to.write <- matrix('', nrow = nrow(kw_disorders.sigProps[[ disorder ]][[ cont ]]),
                             ncol = ncol(kw_disorders.sigProps[[ disorder ]][[ cont ]]))
    rownames(table.to.write) <- rownames(kw_disorders.sigProps[[ disorder ]][[ cont ]])
    colnames(table.to.write) <- colnames(kw_disorders.sigProps[[ disorder ]][[ cont ]])
    
    for (ro in rownames(table.to.write)) {
      for (co in colnames(table.to.write)) {
        if ( kw_disorders.sigProps[[ disorder ]][[ cont ]][ ro, co ] >= 75) {
          table.to.write[ ro, co ] <- kw_disorders.meanPs[[ disorder ]][[ cont ]][ ro, co ]
        } 
      }
    }
    
    #at least 1 good p value within rows
    goodrows <- rownames(table.to.write)[ sapply(rownames(table.to.write), function(ro)
      sum(table.to.write[ ro, ] != "") > 0) ]
    #at least 1 good p within cols
    goodcols <- colnames(table.to.write)[ sapply(colnames(table.to.write), function(co)
      sum(table.to.write[ ,co ] != "") > 0) ]
    
    # remove uninteresting rows and cols
    table.to.write <- table.to.write[ goodrows, goodcols ]
    
    dir.create(sprintf("%s/figures/Kruskal-Wallis/disorders/%s", p2_dir, disorder), showWarnings = F)
    write.csv(table.to.write, file = sprintf("%s/figures/Kruskal-Wallis/disorders/%s/Full_tables.%s.%s.signif_75.meanPs.csv", 
                                             p2_dir, disorder, disorder, cont))
  }
}










# ************************************************ #
# ************************************************ #

### || Chi-square tests ###


# ***************************** #

# chi_disorders.pTabs <- list()
# chi_disorders.sampList <- list()
# chi_disorders.sigProps <- list()
# chi_disorders.meanPs <- list()
# 
# for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
#   print(disorder)
# 
#   disorder.samples <- sample_names(subset_samples(SLL2, SLL2@sam_data[,disorder]=="Yes"))
# 
#   chi_disorders.pTabs[[ disorder ]] <- list()
#   chi_disorders.sampList[[ disorder ]] <- list()
#   chi_disorders.sigProps[[ disorder ]] <- list()
#   chi_disorders.meanPs[[ disorder ]] <- list()
#   
#   for (comparison in c("TvsQ","taxa","questions")) {
#     print(comparison)
#     
#     if (comparison == "questions") {
#       
#       chiTabs <- sub_to_fill_disorderTabs(disorder, disorder.samples, 100, "Chi-squared", comparison, NA)
#       
#       chi_disorders.pTabs[[ disorder ]][[ comparison ]] <- chiTabs[[ "pTabsList" ]]
#       chi_disorders.sampList[[ disorder ]][[ comparison ]] <- chiTabs[[ "sampList" ]]
#       chi_disorders.sigProps[[ disorder ]][[ comparison ]] <- chiTabs[[ "sigProps" ]]
#       chi_disorders.meanPs[[ disorder ]][[ comparison ]] <- chiTabs[[ "meanPs" ]]
#       
#     } else {
#       
#       for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
#         print(sprintf('   %s',tl))
#       
#         chiTabs <- sub_to_fill_disorderTabs(disorder, disorder.samples, 100, "Chi-squared", comparison, tl)
#       
#         chi_disorders.pTabs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "pTabsList" ]]
#         chi_disorders.sampList[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "sampList" ]]
#         chi_disorders.sigProps[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "sigProps" ]]
#         chi_disorders.meanPs[[ disorder ]][[ sprintf("%s.%s",comparison,tl) ]] <- chiTabs[[ "meanPs" ]]
#       }
#     }
#     
#   }
# }
# 
# # ***************************** #
# 
# saveRDS(chi_disorders.pTabs, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorders.pTabs.rds")
# saveRDS(chi_disorders.sampList, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorders.sampList.rds")
# saveRDS(chi_disorders.sigProps, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorders.sigProps.rds")
# saveRDS(chi_disorders.meanPs, 
#         file = "/users/tg/jwillis/SLL/Part_2/R_objects/disorders_controls/chi-squared/chi_disorders.meanPs.rds")

# ***************************** #


chi_disorders.pTabs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorders.pTabs.rds")
chi_disorders.sampList <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorders.sampList.rds")
chi_disorders.sigProps <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorders.sigProps.rds")
chi_disorders.meanPs <- readRDS("/users/tg/jwillis/SLL/Part_2/R_objects/chi_disorders.meanPs.rds")

# write csv files with meanPs at cells where sigProp is at least 75, else blank
for (disorder in c("Celiac","Cystic_fibrosis","Downs_Syndrome")) {
  
  for (comparison in names(chi_disorders.sigProps[[ disorder ]])) {
    
    table.to.write <- matrix('', nrow = nrow(chi_disorders.sigProps[[ disorder ]][[ comparison ]]),
                             ncol = ncol(chi_disorders.sigProps[[ disorder ]][[ comparison ]]))
    rownames(table.to.write) <- rownames(chi_disorders.sigProps[[ disorder ]][[ comparison ]])
    colnames(table.to.write) <- colnames(chi_disorders.sigProps[[ disorder ]][[ comparison ]])
    
    for (ro in rownames(table.to.write)) {
      for (co in colnames(table.to.write)) {
        if ( chi_disorders.sigProps[[ disorder ]][[ comparison ]][ ro, co ] >= 75) {
          table.to.write[ ro, co ] <- chi_disorders.meanPs[[ disorder ]][[ comparison ]][ ro, co ]
        } 
      }
    }
    
    #at least 1 good p value within rows
    goodrows <- rownames(table.to.write)[ sapply(rownames(table.to.write), function(ro)
      sum(table.to.write[ ro, ] != "") > 0) ]
    #at least 1 good p within cols
    goodcols <- colnames(table.to.write)[ sapply(colnames(table.to.write), function(co)
      sum(table.to.write[ ,co ] != "") > 0) ]
    
    # remove uninteresting rows and cols
    table.to.write <- table.to.write[ goodrows, goodcols ]
    
    dir.create(sprintf("%s/figures/Chi-squared/disorders/%s", p2_dir, disorder), showWarnings = F)
    write.csv(table.to.write, file = sprintf("%s/figures/Chi-squared/disorders/%s/Full_tables.%s.%s.signif_75.meanPs.csv", 
                                             p2_dir, disorder, disorder, comparison))
  }
}



# ****************************************************************************************************************** #

















# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####








# ****************************************************************************************************************** #
###### Plotting on map of Spain ######
# ****************************************************************************************************************** #
library(maptools)
library(RColorBrewer)
library(plotrix)
# library(classInt)

# esp.shp <- readShapeSpatial("/users/tg/jwillis/SLL/ESP_adm_shp/ESP_adm0.shp")
esp.shp.com <- readShapeSpatial( sprintf("%s/ESP_adm1.shp", esp_dir) )
esp.shp.prov <- readShapeSpatial( sprintf("%s/ESP_adm2.shp", esp_dir) )
# esp.shp <- readShapeSpatial("/users/tg/jwillis/SLL/ESP_adm_shp/ESP_adm3.shp")
esp.shp.city <- readShapeSpatial( sprintf("%s/ESP_adm4.shp", esp_dir) )
# plot(esp.shp)




plot_map <- function(g, fills, samp_data, otu_data) {
  
  if (fills == "Provinces") { 
    shape <- esp.shp.prov
    shape <- shape[shape$NAME_2 != "Santa Cruz de Tenerife" & shape$NAME_2 != "Las Palmas", ]
    name <- shape$NAME_2
    loc_list <- provinces.sll2
  } else if (fills == "Autonomous Communities") { 
    shape <- esp.shp.com
    shape <- shape[shape$NAME_1 != "Islas Canarias", ]
    name <- shape$NAME_1
    loc_list <- comunidades.sll2
  }
  
  
  cont_vars <- c("pH","Div.Simpson","Div.Shannon","Div.Fisher","Faiths.PD","Species_Richness","Q22.2",
                 "Water_hardness","Weighted_Unifrac","Unweighted_Unifrac","Socioeconomic",cont_water_data)
  
  if (startsWith(g, "Stomatotype")) {
    # % of each autonomous community comprised of each stomatotype
    # By city according to school location (not by postal code of home address bc that was input by the students so may be unreliable):
    ## Apparently they are essentially the same anyway, with only 4 samples in differing locations
    ##    Will use this method for accuracy and because easier to match names to the shp file since I 
    ##    can define the names in this script, instead of relying on names imported from sample data
    g.base <- strsplit(g, '-')[[1]][1]
    shape$Stomatotype1 <- numeric(length=length(name))
    shape$Stomatotype2 <- numeric(length=length(name))
    shape$Stomatotype3 <- numeric(length=length(name))
    names(shape)[startsWith(names(shape), "Stomatotype")] <- c(sprintf("%s-1",g.base), sprintf("%s-2",g.base))
    
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        loc_count <- nrow( samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ] )
        s1_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
                                               rownames(samp_data)[samp_data[ , g.base] == 1 ]) ] )
        s2_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
                                               rownames(samp_data)[samp_data[ , g.base] == 2 ]) ] )
        s3_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
                                               rownames(samp_data)[samp_data[ , g.base] == 3 ]) ] )
        s1_perc <- s1_count/loc_count
        s2_perc <- s2_count/loc_count
        s3_perc <- s3_count/loc_count
        print(c(loc, loc_count, s1_perc, s2_perc, s3_perc))
        
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-1",g.base) ] <- s1_perc
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-2",g.base) ] <- s2_perc
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-3",g.base) ] <- s3_perc
        
      } else {
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-1",g.base) ] <- -1
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-2",g.base) ] <- -1
        attributes(shape)[["data"]][ which(name == loc), sprintf("%s-3",g.base) ] <- -1
      }
    }
    
    # colors <- c('black', rev(brewer.pal(11, "RdBu")))
    colors <- c('black', 'white', brewer.pal(9, "YlOrRd"))
    brks <- c(-1, seq(0, 1, length.out=11))
    legend_labels = c("No data","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100") # leglabs(round(brks, digits=1))
    legend.title <- sprintf('%% of regions with\n%s', g)
    plot_title <- sprintf("Distribution of %s\namongst %s", g, fills)
    leg_x <- 0
    
    # } else if (g %in% c('soft','intermediate','hard','very_hard') ) {
    #   
    # 
    #   shape$hgroup <- numeric(length=length(name))
    #   shape$other_groups <- numeric(length=length(name))
    # 
    #   for (loc in name) {
    #     if (loc %in% names(loc_list)) {
    #       loc_count <- nrow( samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ] )
    #       hg_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
    #                                              rownames(samp_data)[samp_data[ , "Water_hardness_group"] == g ]) ] )
    #       og_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
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
        g_mean <- mean( otu_data[ g, samp_data$School_ID %in% loc_list[loc][[1]] ] )
        g_sd <- sd( otu_data[ g, samp_data$School_ID %in% loc_list[loc][[1]] ] )
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
        cont_vals <- as.numeric(as.matrix(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], g ] ))
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
    
  } else if (g %in% sort(as.character(as.matrix(unique(samp_data[,"Water_type_home"]))))[-4]) {
    
    ## For different types of WATER consumed at home
    for (loc in name) {
      if (loc %in% names(loc_list)) {
        loc_count <- nrow( samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ] )
        agua_count <- tryCatch(nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
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
        # g_mean <- mean( as.numeric(as.matrix(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], g ] )) )
        # g_sd <- sd( as.numeric(as.matrix(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], g ] )) )
        # print(c(loc, g_mean, g_sd))
        # attributes(shape)[["data"]][ which(name == loc), g ] <- g_mean
        
        ###  for g==Q19.1
        loc_count <- nrow( samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ] )
        if (length(intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
                               rownames(samp_data[ samp_data[, g]==1 | samp_data[, g]=="Yes" ]))) == 0) {
          g_count <- 0
        } else {
          g_count <- nrow( samp_data[ intersect(rownames(samp_data[ samp_data$School_ID %in% loc_list[loc][[1]], ]),
                                                rownames(samp_data[ samp_data[, g]==1 | samp_data[, g]=="Yes" ])) ] )
        }
        
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
    leg_x <- 0.45
    
    if (g == "Q19.1") {
      legend.title <- sprintf("%% of region that owns a dog", g)
      plot_title <- sprintf("Distribution of dogs\namongst %s", fills)
    } else if (g == "Q21") {
      legend.title <- sprintf("%% of region with a\nsmoker in the house", g)
      plot_title <- sprintf("Distribution of homes with smokers\namongst %s", fills)
    } else {
      legend.title <- sprintf("%% of region with \n%s", g)
      plot_title <- sprintf("Distribution of homes with smokers\namongst %s", fills)
      # leg_x <- 0
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
  # title(main=plot_title, adj = 0.5, line=-3)
  title(main=g, adj = 0.5, line=-2)
  # boxed.labels(x=coord[,1], y=coord[,2], labels=printNames, bg="white", border=FALSE, cex=1)
  legend(x=leg_x, y=38.8, legend=legend_labels, fill=colors, bty="n", title=legend.title, 
         angle=angles, density=densities, cex = 1.45)
  
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
# g <- cont_water_data[17] # there are 17 values here

plot_map(g, fills, SLL2@sam_data, SLL2@otu_table)

# ****************************************************************************************************************** #




# "Diabetes"                          
# [217] "Hypertension"                       "Cholesterol"                        "Depression"                        
# [220] "Anxiety"                            "Headaches"                          "Lactose_intolerant"                
# [223] "Gastritis"                          "Intestinal_issues"                  "Anemia"                            
# [226] "Sinusitis"                          "Fibrosis_carrier"                   "Thyroid_issue"                     
# [229] "Hypothyroidism"                     "Cancer"                             "Transplant"                        
# [232] "Immune_issues"                      "Skin_issues"                        "Lung_issues"                       
# [235] "Circulatory_issues"                 "Kidney_issues"                      "Tonsil_issues"  


















# ****************************************************************************************************************** ####
# ****************************************************************************************************************** ####

























# ****************************************************************************************************************** #
###### Clustering within disorders / families ######
# ****************************************************************************************************************** #

# First need to get beta diversity values for particular groups of samples
library(foreach)
library(doParallel)

library(fpc)

# *********************************************************** #
B_diversities_subPops <- function(samps, phy) {
  
  phy <- prune_samples(samps, phy)
  
  # remove taxa not seen at least once in at least 25 samples (2.3% of the 1085 samples).
  # This helps protect against an OTU with small mean & trivially large Coefficient of Variation.
  # from /users/tg/jwillis/SLL/Papers/Tools/phyloseq-article-source-files-figs-code-03/phyloseq_plos1_2012-source-doc.html
  phy <- filter_taxa(phy, function(x) sum(x > 0) > 25, prune = TRUE)
  
  # This is required in order to register a parallel "backend" so the calculation can be run in parallel
  registerDoParallel(cores = 10)
  
  # UniFrac diversities
  wu <- UniFrac(phy, weighted = T, parallel = T)
  uu <- UniFrac(phy, weighted = F, parallel = T)
  # zeros seem to cause problems...
  uu[uu==0] <- 0.00001
  
  # Jensen-Shannon Divergence
  jsd <- distance(phy, method = "jsd", parallel = T)
  
  # Bray-Curtis distance combines absolute differences between features (more sensitive to a few large changes)
  bray <- vegdist(t(phy@otu_table), method = "bray")
  # Canberra distance weights all differences equally (more sensitive to many small changes)
  canberra <- vegdist(t(phy@otu_table), method = "canberra")
  
  return(list( "JSD"=jsd, "Weighted_Unifrac"=wu, "Unweighted_Unifrac"=uu, "Bray"=bray, "Canberra"=canberra ))
}
# *********************************************************** #
prediction_strengths_subPops <- function(B_diversities) {
  ps <- list()
  
  for (beta in names(B_diversities)) {
    print(beta)
    ps[[ beta ]] <- prediction.strength(B_diversities[[ beta ]], cutoff = 0.75)
  }

  return(ps)
}
# *********************************************************** #

disorder <- "Cystic_fibrosis"
disorder.samples <- rownames(SLL2@sam_data[ SLL2@sam_data[,disorder]=="Yes", ])

B_divs <- B_diversities_subPops(disorder.samples, SLL2)
ps <- prediction_strengths_subPops(B_divs)


# ***************************************** #

tagl <- prune_samples(disorder.samples, SLL2_rel)
tagl <- filter_taxa(tagl, function(x) sum(x > 0) > 25, prune = TRUE)

dists <- c("JSD","Weighted_Unifrac","Unweighted_Unifrac","Bray","Canberra")

clus.subPops <- list()
for (dist_meas in dists) {
  
  print(dist_meas)
  
  clus.subPops[[ dist_meas ]] <- get_clusters( tagl, "Species", dist_meas, ps, subPops = T, subPop_Bdivs = B_divs )
}
# ***************************************** #






# ***************************************************************************************** #
# Plot 2: principal coordinates analysis (PCoA) ####
library(scales)
library(viridis)

pcoas.subPops <- list()
for (dis in c("JSD","Weighted_Unifrac","Unweighted_Unifrac","Bray","Canberra")) {
  print(dis)
  # pcoas.subPops[[ dis ]] <- dudi.pco(clus.subPops[[ dis ]][[ "dist.matrix" ]], scannf=F, nf=3)
  pcoas.subPops[[ dis ]] <- ape::pcoa(clus.subPops[[ dis ]][[ "dist.matrix" ]])
}


dis <- "Bray"

col.to.plot <- "Transplant"; ctp.factor <- as.factor(SLL2.meta[ rownames(pcoas.subPops[[ dis ]]$li), col.to.plot ])
col.to.plot <- sprintf("Stomatotype %s",dis); ctp.factor <- as.factor(clus.subPops[[ dis ]]$cluster_full)


ade4::s.class(pcoas.subPops[[ dis ]]$li, 
              fac=ctp.factor,
              grid=F, clabel=1.5, 
              col=hue_pal()(length(unique(ctp.factor))), 
              sub = sprintf("PCoA separated by %s", col.to.plot ))

# ***************************************************************************************** #



# ***************************************************************************************** #
tagl <- prune_samples(disorder.samples, SLL2_rel)
tagl <- filter_taxa(tagl, function(x) sum(x > 0) > 25, prune = TRUE)

dis <- "Weighted_Unifrac"; subPop.clus <- clus.subPops[[ dis ]]$cluster_full

kw.subPops <- list()
for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
  kw.subPops[[ sprintf('p_%s',cont) ]] <- fill_kw_p_table(tagl, cont, disorder.samples, disorder, 
                                                          only_cont, externalGroupQ = subPop.clus)
  
  kw.subPops[[ sprintf('p.adj_%s',cont) ]] <- adjust_kw_p_vals(kw.subPops[[ sprintf('p_%s',cont) ]])
}


# ********************************* #
# ********************************* #

kw_full.subPops <- list()
for (cont in c("cont_vars","Phylum","Class","Order","Family","Genus","Species")) {
  kw_full.subPops[[ sprintf('p_%s',cont) ]] <- fill_kw_p_table(tagl, cont, disorder.samples, groupQs, only_cont)
  
  kw_full.subPops[[ sprintf('p.adj_%s',cont) ]] <- adjust_kw_p_vals(kw_full.subPops[[ sprintf('p_%s',cont) ]])
}


# ********************************* #


tlev <- "Species"; o_rel <- gloms_rel[[ tlev ]]
cont_col <- "Shuttleworthia"

group_vs_cont_box(disorder.samples, tlev, o_rel, "CF_Stomatotype-Weighted_Unifrac", cont_col, externalGroupQ = subPop.clus)
# ***************************************************************************************** #






# ***************************************************************************************** #
tagl <- prune_samples(disorder.samples, SLL2_rel)
tagl <- filter_taxa(tagl, function(x) sum(x > 0) > 25, prune = TRUE)

dis <- "Weighted_Unifrac"; subPop.clus <- clus.subPops[[ dis ]]$cluster_full

chi.subPops <- list()
for (comparison in c("TvsQ","taxa","questions")) {
  print(comparison)
  
  if (comparison == "questions") {
    # fill tables of p values
    chi.subPops[[ "questions" ]] <- fill_chi_table(comparison, NA, disorder.samples, groupQs, 
                                                   externalGroupQ = subPop.clus, 
                                                   externalGroupName = sprintf("CF_Stomatotype-%s", dis))
    
  } else {
    for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
      print(sprintf('   %s',tl))
      
      # fill tables of p values
      chi.subPops[[ sprintf("%s.%s",comparison,tl) ]] <- fill_chi_table(comparison, tl, disorder.samples, groupQs, 
                                                                        externalGroupQ = subPop.clus, 
                                                                        externalGroupName = sprintf("CF_Stomatotype-%s", dis))
      
    }
  }
}




# ***************************************************************************************** #









# ***************************************************************************************** #
# ***************************************************************************************** #
tagl <- prune_samples(disorder.samples, SLL2_rel)
tagl <- filter_taxa(tagl, function(x) sum(x > 0) > 25, prune = TRUE)

dis <- "Weighted_Unifrac"; subPop.clus <- clus.subPops[[ dis ]]$cluster_full

cor.subPops <- list()

for (comparison in c("TvsQ","taxa","questions")) {
  print(comparison)
  
  if (comparison %in% c("taxa","questions")) {
    # fill tables of p values
    cor.ps <- fill_cor_tables(comparison, NA, disorder.samples, only_cont)
    
    # adjust p-values
    cor.ps.adj <- adjust_cor_p_vals(cor.ps[[ "cor" ]], cor.ps[[ "ps" ]], comparison)
    cor.subPops[[ sprintf("cor.%s",comparison) ]] <- cor.ps.adj[[ "cor.adj" ]]
    cor.subPops[[ sprintf("ps.%s",comparison) ]] <- cor.ps.adj[[ "p.adj" ]]
    
    # write cor tables to files -- only print significant correlations
    cor.subPops[[ sprintf("cor.to.write.%s",comparison) ]] <- cor.ps.adj[[ "cor.adj" ]]
    cor.subPops[[ sprintf("cor.to.write.%s",comparison) ]][ cor.ps.adj[[ "p.adj" ]] >= 0.05 ] <- ''
    
  } else {
    for (tl in c("Phylum","Class","Order","Family","Genus","Species")) {
      print(sprintf('   %s',tl))
      # fill tables of p values
      cor.ps <- fill_cor_tables(comparison, tl, disorder.samples, only_cont)
      
      # adjust p-values
      cor.ps.adj <- adjust_cor_p_vals(cor.ps[[ "cor" ]], cor.ps[[ "ps" ]], comparison)
      cor.subPops[[ sprintf("cor.%s.%s",comparison,tl) ]] <- cor.ps.adj[[ "cor.adj" ]]
      cor.subPops[[ sprintf("ps.%s.%s",comparison,tl) ]] <- cor.ps.adj[[ "p.adj" ]]
      
      # write cor tables to files -- only print significant correlations
      cor.subPops[[ sprintf("cor.to.write.%s.%s",comparison,tl) ]] <- cor.ps.adj[[ "cor.adj" ]]
      cor.subPops[[ sprintf("cor.to.write.%s.%s",comparison,tl) ]][ cor.ps.adj[[ "p.adj" ]] >= 0.05 ] <- ''
      
    }
  }
}

# ***************************************************************************************** #



# check out medications in CF samples

CF_meds <- sample_names(subset_samples(SLL2, Cystic_fibrosis=="Yes" & ! Other_medications %in% c("0","Illegible or skipped")))
CF_no_meds <- sample_names(subset_samples(SLL2, Cystic_fibrosis=="Yes" & Other_medications %in% c("0","Illegible or skipped")))

table(subPop.clus[ CF_meds ])
table(subPop.clus[ CF_no_meds ])

cbind(SLL2@sam_data[CF_meds, "Other_medications"], subPop.clus[CF_meds])
# ***************************************************************************************** #
































# ******************************************************************** #
# ******************************************************************** #
# ******************************************************************** #

sll2.survey <- read.delim("/users/tg/jwillis/SLL/Part_2/Surveys/SLL_survey_part2_total_repaired.csv", row.names = 1)

s2s <- sll2.survey[colnames(otu.table.rel.2)[ colnames(otu.table.rel.2) %in% rownames(sll2.survey) ], ]



# using the date March 7th 2017 since samples were collected between January 24 and April 26 of 2017, so this is the middle date
# just a simpler method than using the exact dates of collections to get the ages in years. 
# ages <- sapply(as.Date(as.matrix(s2s[, "Q1"])), function(x) as.numeric( difftime("2017-07-10", x, units = "days") / 365) )
ages <- as.numeric(sapply(as.character(s2s[, "Q1"]), function(x) 
  ifelse(nchar(x)==10,
         ifelse( is.na(as.Date(x, "%d-%m-%Y")),
                 ifelse( is.na(as.Date(strsplit(x,'-')[[1]][3], "%Y")),
                         NA,
                         as.numeric( difftime("2017-03-07", as.Date(strsplit(x,'-')[[1]][3], "%Y"), units = "days") / 365)),
                 as.numeric( difftime("2017-03-07", as.Date(x, "%d-%m-%Y"), units = "days") / 365)),
         NA)))
# give those with ages less than 1 year as NA because these were simply typing errors by the participants that put 2014 or 2015 as their birthday
ages[ is.na(ages) ] <- "No Sabe/No Contesta"

s2s[ , "Age"] <- ages



s2s[ as.numeric(s2s$Age)>49 & as.numeric(s2s$Age)<51 & s2s$Q39.1==1 & s2s$Q2==1, "Age"] # 1
m50C <- rownames(s2s)[ as.numeric(s2s$Age)>49 & as.numeric(s2s$Age)<51 & s2s$Q39.1==1 & s2s$Q2==1]

s2s[ as.numeric(s2s$Age)>39 & as.numeric(s2s$Age)<41 & s2s$Q39==0 & s2s$Q2==0, "Age"] # 2
f40 <- rownames(s2s)[ as.numeric(s2s$Age)>39 & as.numeric(s2s$Age)<41 & s2s$Q39==0 & s2s$Q2==0][2]

s2s[ as.numeric(s2s$Age)>54 & as.numeric(s2s$Age)<56 & s2s$Q39==0 & s2s$Q2==1, "Age"] # 2
m55 <- rownames(s2s)[ as.numeric(s2s$Age)>54 & as.numeric(s2s$Age)<56 & s2s$Q39==0 & s2s$Q2==1][2]

s2s[ as.numeric(s2s$Age)>44 & as.numeric(s2s$Age)<46 & s2s$Q39.2==1 & s2s$Q2==0, "Age"] # 1
f45f <- rownames(s2s)[ as.numeric(s2s$Age)>44 & as.numeric(s2s$Age)<46 & s2s$Q39.2==1 & s2s$Q2==0][1]

s2s[ as.numeric(s2s$Age)>13 & as.numeric(s2s$Age)<15 & s2s$Q39.4==1 & s2s$Q2==1, "Age"] # 3
m14d <- rownames(s2s)[ as.numeric(s2s$Age)>13 & as.numeric(s2s$Age)<15 & s2s$Q39.4==1 & s2s$Q2==1][3]

s2s[ as.numeric(s2s$Age)>15 & as.numeric(s2s$Age)<17 & s2s$Q39.2==1 & s2s$Q2==1, "Age"] # 3
m16f <- rownames(s2s)[ as.numeric(s2s$Age)>15 & as.numeric(s2s$Age)<17 & s2s$Q39.2==1 & s2s$Q2==1][1]

s2s[ as.numeric(s2s$Age)>14 & as.numeric(s2s$Age)<16 & s2s$Q39.1==1 & s2s$Q2==0, "Age"] # 2
f15c <- rownames(s2s)[ as.numeric(s2s$Age)>14 & as.numeric(s2s$Age)<16 & s2s$Q39.1==1 & s2s$Q2==0][2]

min(as.numeric(s2s[ s2s$Q39==0 & s2s$Q2==1, "Age"])) # 10.31 years
ym <- rownames(s2s)[as.numeric(s2s$Age) == min(as.numeric(s2s[ s2s$Q39==0 & s2s$Q2==1, "Age"])) & ! is.na(as.numeric(s2s$Age))]

min(as.numeric(s2s[ s2s$Q39.4==1 &s2s$Q2==0, "Age"])) # 7.70 years
yfd <- rownames(s2s)[as.numeric(s2s$Age) == min(as.numeric(s2s[ s2s$Q39.4==1 &s2s$Q2==0, "Age"])) & ! is.na(as.numeric(s2s$Age))]

max(as.numeric(s2s[ s2s$Q39.2==1 &s2s$Q2==1, "Age"], decreasing = T)) # 47.30 years
omf <- rownames(s2s)[as.numeric(s2s$Age) == max(as.numeric(s2s[ s2s$Q39.2==1 &s2s$Q2==1, "Age"])) & ! is.na(as.numeric(s2s$Age))]

s2s[ as.numeric(s2s$Age)>74 & as.numeric(s2s$Age)<76 & s2s$Q39.1==0 & s2s$Q2==0, "Age"] # 4
f75 <- rownames(s2s)[ as.numeric(s2s$Age)>74 & as.numeric(s2s$Age)<76 & s2s$Q39.1==0 & s2s$Q2==0][4]

cards <- c(m50C, f40, m55, f45f, m14d, m16f, f15c, ym, yfd, omf, f75)




ot <- otu.table.rel.2[c("Streptococcus","Prevotella","Haemophilus","Neisseria","Veillonella","Porphyromonas","Gemella"), cards]
ot <- cbind(rowMeans(otu.table.rel.2[c("Streptococcus","Prevotella","Haemophilus","Neisseria","Veillonella","Porphyromonas","Gemella"), ]), ot)
ot.norm <- round(6 * ot/max(ot))

write.csv(ot.norm, "~/Downloads/table.csv")


