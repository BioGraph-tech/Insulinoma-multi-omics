
library(Seurat)
library(sceasy)
library(reticulate)
# use_condaenv(
#   conda = "C:/ProgramData/miniforge3/condabin/conda.bat",
#   condaenv = "scanpy",
#   required = TRUE
# )
use_condaenv(
  condaenv = "scanpy",
  required = TRUE
)
loompy <- reticulate::import('loompy')
sceasy::convertFormat('h5ad/islet_final_annotated.h5ad', from="anndata", to="seurat",outFile='Rds/islet_final_annotated.Rds')
sceasy::convertFormat('h5ad/scATAC_peak_final_annoated.h5ad', from="anndata", to="seurat",outFile='Rds/scATAC_peak_final_annoated.Rds')
sceasy::convertFormat('integrate_h5ad/islet_beta_PanNET_intergtated.h5ad', from="anndata", to="seurat",outFile='Rds/islet_beta_PanNET_intergtated.Rds')
