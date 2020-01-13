library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)

beak_files <- list.files("data/oriented_trimmed_meshes_obj_circle/", full.names = TRUE) 

beak_file <- beak_files[16]

beak_mesh <- try(readobj::read.obj(beak_file, convert.rgl = TRUE)[[1]])
rgl::shade3d(beak_mesh)


reorient_beak <- function(beak_file, sdf_folder = "data/sdfs", plot_folder = NULL, save_torch = TRUE) {
  
  if(!dir.exists(sdf_folder)) {
    dir.create(sdf_folder)
  }
  
  bird_file <- basename(beak_file)
  
  sdf_file <- file.path(sdf_folder, bird_file)
  
  if(!is.null(plot_folder)) {
    if(!dir.exists(plot_folder)) {
      dir.create(plot_folder)
    }
  }
  
  
  if(!file.exists(sdf_file)) {
    
    beak_mesh <- try(readobj::read.obj(beak_file, convert.rgl = TRUE)[[1]])
    
    if(inherits(beak_mesh, 'try-error') | inherits(beak_landmarks, 'try-error')) {
      bad_meshes <- unique(c(readr::read_lines("data/bad_mesh.txt"), bird_file))
      readr::write_lines(bad_meshes, "data/bad_mesh.txt")
    } else {
      
        mesh_pnts <- Morpho::vert2points(beak_mesh)
        
        mesh_cent <- apply(mesh_pnts, 2, mean, na.rm = TRUE)
        furthest <- pdist::pdist(mesh_cent, 
                                 mesh_pnts) %>%
          as.matrix() %>%
          max(na.rm = TRUE)
        
        scaled_mesh <- beak_mesh %>%
          rgl::translate3d(-mesh_cent[1], -mesh_cent[2], -mesh_cent[3]) %>%
          Morpho::scalemesh(1/furthest, center = "none")
        
        mesh_pnts <- Morpho::vert2points(scaled_mesh)
        
        max_x_pnt <- mesh_pnts[which.max(mesh_pnts[ , 1]), ]
        furthest <- pdist::pdist(mesh_pnts[which.max(mesh_pnts[ , 1]), , drop = FALSE], 
                                 mesh_pnts[-which.max(mesh_pnts[, 1]), ]) %>%
          as.matrix() %>%
          max(na.rm = TRUE)

        
        sphere <- Rvcg::vcgSphere(subdivision = 5) %>%
          Morpho::scalemesh(furthest) %>%
          rgl::translate3d(max_x_pnt[1], max_x_pnt[2], max_x_pnt[3])
        
        sdf <- Rvcg::vcgClostKD()
        
        rgl::shade3d(scaled_mesh, alpha = 0.5)
        rgl::shade3d(Rvcg::vcgSphere(5), alpha = 0.6, col = "red")
        rgl::shade3d(sphere, alpha = 0.6, col = "red")
          
        
        if(!is.null(plot_folder)) {
          
          rgl::clear3d()
          rgl::shade3d(cut_1, col = "red")
          rgl::rgl.bringtotop()
          rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                                   -0.32, 0.58, 0.74, 0,
                                                   0.69, -0.38, 0.61, 0,
                                                   0, 0, 0, 1), nrow = 4, byrow = TRUE),
                             zoom = 0.8)
          rgl::par3d(windowRect = c(20, 30, 800, 800))
          
          rgl::rgl.snapshot(file.path(plot_folder, gsub(".obj", ".png", bird_file, fixed = TRUE)))
          
        }
        
        
      } else {
        Rvcg::vcgObjWrite(rot$mesh, new_beak_file)
        
        if(!is.null(plot_folder)) {
          rgl::clear3d()
          rgl::shade3d(rot$mesh, col = "red")
          rgl::rgl.bringtotop()
          rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                                   -0.32, 0.58, 0.74, 0,
                                                   0.69, -0.38, 0.61, 0,
                                                   0, 0, 0, 1), nrow = 4, byrow = TRUE),
                             zoom = 0.8)
          rgl::par3d(windowRect = c(20, 30, 800, 800))
          
          rgl::rgl.snapshot(file.path(plot_folder, gsub(".obj", ".png", bird_file, fixed = TRUE)))
        }
      }
    }
  }
  return(invisible(NULL))
}