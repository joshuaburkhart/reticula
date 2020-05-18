# R script to download selected samples
# Copy code and run on a local machine to initiate download

# Check for dependencies and install if missing
packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
}
library("rhdf5")
library("tools")

destination_file = "human_matrix_download.h5"
extracted_expression_file = "GSE58135_expression_matrix.tsv"
url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"

# Check if gene expression file was already downloaded and check integrity, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE)
}

# Selected samples to be extracted
samp = c("GSM1401701","GSM1401724","GSM1401678","GSM1401784","GSM1401710","GSM1401743","GSM1401762","GSM1401726","GSM1401759","GSM1401742","GSM1401755","GSM1401794","GSM1401757","GSM1401785","GSM1401730","GSM1401796","GSM1401767","GSM1401708","GSM1401738","GSM1401741","GSM1401653","GSM1401671","GSM1401811","GSM1401721","GSM1401806","GSM1401682","GSM1401771","GSM1401773","GSM1401663","GSM1401686","GSM1401732",
"GSM1401685","GSM1401792","GSM1401795","GSM1401695","GSM1401763","GSM1401665","GSM1401705","GSM1401706","GSM1401694","GSM1401667","GSM1401722","GSM1401718","GSM1401812","GSM1401697","GSM1401674","GSM1401789","GSM1401734","GSM1401699","GSM1401800","GSM1401691","GSM1401657","GSM1401662","GSM1401680","GSM1401793","GSM1401775","GSM1401815","GSM1401740","GSM1401737","GSM1401764","GSM1401729",
"GSM1401693","GSM1401791","GSM1401766","GSM1401783","GSM1401669","GSM1401654","GSM1401703","GSM1401700","GSM1401652","GSM1401670","GSM1401770","GSM1401684","GSM1401704","GSM1401690","GSM1401809","GSM1401761","GSM1401804","GSM1401713","GSM1401696","GSM1401676","GSM1401778","GSM1401672","GSM1401698","GSM1401808","GSM1401782","GSM1401777","GSM1401758","GSM1401745","GSM1401735","GSM1401769",
"GSM1401715","GSM1401649","GSM1401787","GSM1401754","GSM1401739","GSM1401753","GSM1401707","GSM1401683","GSM1401751","GSM1401658","GSM1401749","GSM1401720","GSM1401675","GSM1401716","GSM1401688","GSM1401801","GSM1401781","GSM1401805","GSM1401728","GSM1401650","GSM1401655","GSM1401673","GSM1401717","GSM1401772","GSM1401659","GSM1401768","GSM1401711","GSM1401648","GSM1401733","GSM1401651",
"GSM1401803","GSM1401797","GSM1401709","GSM1401799","GSM1401760","GSM1401756","GSM1401727","GSM1401677","GSM1401802","GSM1401687","GSM1401712","GSM1401774","GSM1401776","GSM1401746","GSM1401664","GSM1401788","GSM1401731","GSM1401798","GSM1401679","GSM1401790","GSM1401668","GSM1401750","GSM1401813","GSM1401714","GSM1401747","GSM1401725","GSM1401666","GSM1401661","GSM1401692","GSM1401689",
"GSM1401810","GSM1401814","GSM1401752","GSM1401779","GSM1401719","GSM1401736","GSM1401744","GSM1401748","GSM1401780","GSM1401786","GSM1401681","GSM1401765","GSM1401702","GSM1401656","GSM1401723","GSM1401660","GSM1401807","")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")

# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

