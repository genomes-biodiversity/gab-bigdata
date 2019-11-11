##########################################################################################

# Install package and load library:

install.packages(PopGenome)

library(PopGenome)

##########################################################################################

# Read the VCF file into R:

Mala <- readData("VCF", format = "VCF")

# Set the population structure:

Mainland <- c("S53491", "S53492", "S53493", "S53494", "S53496", 
              "S53497", "S53499", "S53500", "S52738")

Isl_I <- c("S52700", "S53339", "S53340", "S53341", "S53342", 
             "S53343", "S53345", "S53346", "S53347", "S53348")

Isl_II <- c("S53349", "S53350", "S53351", "S53352", "S53353", 
           "S53354", "S53355", "S53356", "S53357", "S53358")

OTran_I <- c("S53360", "S53361", "S53362", "S53363", "S53364",
                "S53365", "S53366", "S53367", "S53368", "S53369")

OTran_II <- c("S53389", "S53390", "S53391", "S53392", "S53393",
              "S53394", "S53395", "S53396", "S53397", "S53398")

TTran_I <- c("S53370", "S53371", "S53372", "S53373", "S53374", 
                "S53375", "S53376", "S53377", "S53378")

TTran_II <- c("S53379", "S53380", "S53381", "S53382", "S53383",
            "S53384", "S53385", "S53386", "S53387", "S53388")

Mala <- set.populations(Mala,list(Mainland, Isl_I, Isl_II, OTran_I, OTran_II, TTran_I, 
        TTran_II), diploid = TRUE)

# Check that the population assignment worked:

Mala@populations

##########################################################################################

# DIVERSITY

# Calculate diversity within populations:

Mala <- diversity.stats(Mala)

Mala@nuc.diversity.within

# Calculate per-site diversity (genome length is 872058):

Genome_length <- 872058

Mainland_div <- (493.0327/Genome_length)*100

Mainland_div

# Repeat the above for the remaining diversity values:

Isl_I_div <- (143.1789/Genome_length)*100

Isl_I_div

Isl_II_div <- (163.1105/Genome_length)*100

Isl_II_div

# Plot the diversity:

Nucleotide_diversity <- c(Mainland_div, Isl_I_div, Isl_II_div)

barplot(Nucleotide_diversity, main = "Mala nucleotide diversity; Mainland vs. Islands", 
        ylab = "Nucleotide diversity", ylim = c(0,0.07), 
        names.arg=c("Mainland", "Island_I", "Island_II"), las=1, col="blue")


# Compare the nucleotide diversity between mainland and island populations.  What is the evolutionary outcome of the founder effect?

##########################################################################################

# DIVERSITY

# Compare translocations:

OTran_I_div <- (321.7053/Genome_length)*100

OTran_II_div <- (310.1158/Genome_length)*100

TTran_I_div <- (305.9804/Genome_length)*100

TTran_II_div <- (350.9737/Genome_length)*100

# Plot the diversity:

Nucleotide_diversity <- c(Mainland_div, OTran_I_div, OTran_II_div, TTran_I_div, TTran_II_div)

barplot(Nucleotide_diversity, main = "Mala nucleotide diversity; Mainland vs. Translocated populations", 
        ylab = "Nucleotide diversity", ylim = c(0,0.07), 
        names.arg=c("Mainland", "One_TR_I", "One_TR_I", "Two_TR_I", "Two_TR_II"), las=1, col="darkgreen")

# OTran_I and OTran_II populations underwent 1 translocation, while TTran_I and TTran_II have undergone 2.  
# What is the effect of translocation on nucleotide diversity?

# How do these results compare to the mainland vs. Island results?

##########################################################################################

# DIFFERENTIATION

# Calculate pairwise FST

Mala <- F_ST.stats(Mala)

Mala@nuc.F_ST.pairwise

##########################################################################################
