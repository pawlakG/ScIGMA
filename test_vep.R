library(httr)
library(jsonlite)

query <- "17:g.7674220C>T" # TP53 c.700T>C (in hg38 this is chr17:7674220 C>T approx)
# Let's query using hgvs directly using GRCh37 or GRCh38 if known. 
# Better yet, let's just query the VEP API with a known VCF format variant to see the structure.
res <- GET("https://rest.ensembl.org/vep/human/region/17:7578212-7578212:1/C?phenotypes=1",
           add_headers("Accept" = "application/json"))
data <- content(res, "text", encoding = "UTF-8")
print(substr(data, 1, 1000))
