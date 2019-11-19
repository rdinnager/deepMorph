get_meshes_from_latent <- function(latent_code, voxel_res = 64L) {
  current_wd <- getwd()
  setwd(file.path(current_wd, "sdf"))
  
  torch <- reticulate::import("torch")
  np <- reticulate::import("numpy")
  trimesh <- reticulate::import("trimesh")
  
  sdf_net <- reticulate::import_from_path("sdf_net", "sdf/model")
  
}