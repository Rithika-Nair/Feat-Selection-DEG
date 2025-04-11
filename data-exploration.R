# Load the MAE object from your .rds file
mae_rosmap <- readRDS("rosmap/mae_data/mae.rds")

mae <- mae_GSE71669
# Check the class of the object
class(mae)

# Print the summary of the MAE object
summary(mae)

# List all assays (experiments) contained in the MAE object
experiments(mae)

# Get a quick overview of the metadata associated with samples (colData)
colData(mae)
metadata(mae)

# Loop through each experiment name and extract the dimensions
experiment_names <- names(experiments(mae))

assay_dims <- lapply(experiment_names, function(exp_name) {
  dim(experiments(mae)[[exp_name]])  # Extract each experiment
})

# Assign names to the list for easier reference
names(assay_dims) <- experiment_names

# Access the assay data
assay_data <- experiments(mae)[["mrna"]]
class(assay_data)  

# Check if the data is HDF5-backed
is(assay_data, "HDF5Array")  # Should return TRUE if it's an HDF5-backed object

# Load the HDF5 file containing assay data (if it's separate from the RDS)
h5_file <- HDF5Array::H5File("GSE71669/mae_data/experiments.h5")
h5_file

class(h5_file)
# If the extracted path is another .h5 file, load it
h5_real_data <- HDF5Array::H5ADMatrix("/scratch/st-singha53-1/nair1602/multiomics_project/GSE71669/mae_data/experiments.h5")

# Check the structure
dim(h5_real_data)
