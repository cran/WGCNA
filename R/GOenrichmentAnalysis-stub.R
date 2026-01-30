# Function that returns GO terms with best enrichment and highest number of genes.

# So that I don't forget: information about GO categories is contained in the GO.db package

GOenrichmentAnalysis = function(...)
{
   stop(
     spaste("This function is deprecated and will be removed in the near future. \n",
            "We suggest using the replacement function enrichmentAnalysis.Entrez \n",
            "in R package anRichment, available from github at\n",
            "https://github.com/plangfelder/anRichment"));
}
