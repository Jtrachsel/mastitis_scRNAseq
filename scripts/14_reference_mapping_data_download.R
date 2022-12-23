

# dont need this one yet.  no cell type annotations
# BECKER <- LoadH5Seurat('reference_mapping_data/Becker2021_Milkcells.h5Seurat')
# BECKER@meta.data
# https://www.10xgenomics.com/resources/datasets/whole-blood-rbc-lysis-for-pbmcs-neutrophils-granulocytes-3-3-1-standard


# only need blood from this one
# https://data.humancellatlas.org/explore/projects/10201832-7c73-4033-9b65-3ef13d81656a

# tabula sapiens metadata
# https://data.humancellatlas.org/explore/projects/10201832-7c73-4033-9b65-3ef13d81656a

# tabula sapiens 
# https://data.humancellatlas.org/explore/projects/10201832-7c73-4033-9b65-3ef13d81656a
# generated this curl command from the human cell atlas website.
if (!file.exists('reference_mapping_data/TabulaSapiens.h5ad')){
  
  system("curl --location --fail 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp22&format=curl&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2210201832-7c73-4033-9b65-3ef13d81656a%22%5D%7D%2C+%22genusSpecies%22%3A+%7B%22is%22%3A+%5B%22Homo+sapiens%22%5D%7D%2C+%22fileFormat%22%3A+%7B%22is%22%3A+%5B%22h5ad.zip%22%5D%7D%7D&objectKey=manifests%2F021fb812-3225-5aee-a408-5208d5ae1861.b1e0bb04-4675-58d1-be1e-0f2b0eb408fe.curlrc' | curl --config -")
  system('mv 4de7ab62-475c-43f9-a1d4-bfd693e7df07/TabulaSapiens.h5ad.zip reference_mapping_data/')
  system('unzip reference_mapping_data/TabulaSapiens.h5ad.zip')
  
}

# only need whole blood from these data...



####  Human Milk Data 
# These are the urls but you need to login and authenticate, no good way to download programatticalyl 
# nyquist_milk_counts_url <- 'https://singlecell.broadinstitute.org/single_cell/data/public/SCP1671/cellular-and-transcriptional-diversity-over-the-course-of-human-lactation?filename=MIT_Milk_Study_Raw_counts.txt.gz'
# nyquist_milk_metadata_url <- 'https://singlecell.broadinstitute.org/single_cell/data/public/SCP1671/cellular-and-transcriptional-diversity-over-the-course-of-human-lactation?filename=MIT_milk_study_metadata.csv.gz'

# download.file(url=nyquist_milk_counts_url, destfile = 'reference_mapping_data/nyquist_milk_counts.txt.gz')
# download.file(url=nyquist_milk_metadata_url, destfile = 'reference_mapping_data/nyquist_milk_metadata.csv.gz')

# curl()
