library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)
library(geometry)
library(proxy)

## Step one is to iterate over each mesh, reorient the beak by the landmarks, and then
## rescale each beak to constant beak length. I will also save the convex hull of 
## each resulting mesh to use in the next step.

beak_files <- list.files("data/objs/", full.names = TRUE) 

beak_file <- beak_files[1]

clip_mesh_by_mesh <- function(mesh_to_clip, mesh_to_clip_with) {
  #rgl::shade3d(mesh_to_clip, alpha = 0.5, col = "red")
  #rgl::shade3d(mesh_to_clip_with, alpha = 0.5, col = "green")
  dists <- Rvcg::vcgClostKD(t(mesh_to_clip$vb), mesh_to_clip_with)
  clipped <- rgl::clipMesh3d(mesh_to_clip, dists$quality)
  clipped
  #rgl::shade3d(clipped)
}

#beak_file <- beak_files[168]
#plot_folder <- "beak_plots_trimmed"
reorient_beak <- function(beak_file, mesh_folder = "data/oriented_scaled_meshes_obj",
                          landmark_folder = "data/landmarks_oriented", hull_folder = "data/conv_hulls") {
  
  if(!dir.exists(mesh_folder)) {
    dir.create(mesh_folder)
  }
  if(!dir.exists(landmark_folder)) {
    dir.create(landmark_folder)
  }
  
  bird_file <- basename(beak_file)

  new_beak_file <- file.path(mesh_folder, bird_file)
  hull_file <- file.path(hull_folder, gsub(".obj", ".csv", bird_file))
  
  if(!dir.exists(hull_folder)) {
    dir.create(hull_folder)
  }
  
  
  
  if(!file.exists(new_beak_file)) {
    
    beak_mesh <- try(readobj::read.obj(beak_file, convert.rgl = TRUE)[[1]])
    suppressMessages(beak_landmarks <- try(readr::read_delim(file.path("data/landmarks/Markup_1", gsub(".obj", ".txt", bird_file, fixed = TRUE)), 
                                                             " ",
                                                             col_names = c("X", "Y", "Z"))))
    
    if(inherits(beak_mesh, 'try-error') | inherits(beak_landmarks, 'try-error')) {
      readr::write_lines(bird_file, "data/bad_models.txt", append = TRUE)
    } else {
      
      main_landmarks <- beak_landmarks[1:4, ] %>%
        as.matrix()
      
      mid_beak <- apply(main_landmarks[3:4, ], 2, mean) %>% matrix(nrow = 1)
      line_1 <- rbind(mid_beak, main_landmarks[1, ])
      
      line_1_length <- dist(line_1)
      line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / line_1_length
      
      top_line <- main_landmarks[2, ] - mid_beak
      
      tt <- top_line %*% t(line_1_vec)
      
      point_on_line_1 <- mid_beak + as.vector(tt) * line_1_vec
      
      landmarks_3 <- rbind(main_landmarks[1, ], point_on_line_1, main_landmarks[2, ])
      landmark_dists <- dist(landmarks_3) %>% as.matrix()
      target_3 <- matrix(c(landmark_dists[1, 2], 0, 0,
                           0, 0, 0,
                           0, 0, landmark_dists[2, 3]),
                         nrow = 3, byrow = TRUE)
      
      L_landmark_folder <- file.path(landmark_folder, "L_landmarks")
      if(!dir.exists(L_landmark_folder)) {
        dir.create(L_landmark_folder)
      }
      
      all_landmark_folder <- file.path(landmark_folder, "all_landmarks")
      if(!dir.exists(all_landmark_folder)) {
        dir.create(all_landmark_folder)
      }
      
      readr::write_csv(as.data.frame(target_3), file.path(L_landmark_folder, gsub(".obj", ".csv", bird_file, fixed = TRUE)),
                       col_names = FALSE)
      
      rot_landmarks <- beak_landmarks %>% as.matrix() %>%
        Morpho::rotonmat(landmarks_3, target_3)
      
      readr::write_csv(as.data.frame(rot_landmarks), file.path(all_landmark_folder, gsub(".obj", ".csv", bird_file, fixed = TRUE)),
                       col_names = FALSE)
      
      rot <- Morpho::rotmesh.onto(beak_mesh, landmarks_3, target_3, adnormals = TRUE)
      
      # rgl::shade3d(rot$mesh, alpha = 0.1)
      # rgl::lines3d(target_3, col = "red")
      # rgl::points3d(rot_landmarks, size = 4, col = "blue")
      # rgl::shade3d(sphere, col = "red", alpha = 0.5)
      
      ## cut by sphere
      dist_scale <- dist(target_3[c(1, 3), ])
      
      sphere <- Rvcg::vcgSphere(5) %>%
        rgl::translate3d(target_3[1, 1], target_3[1, 2], target_3[1, 3]) %>%
        Morpho::scalemesh(dist_scale)
      
      clipped <- clip_mesh_by_mesh(rot$mesh, sphere) %>%
        Rvcg::vcgClean(sel = 0:6, iterate = TRUE, silent = TRUE)
      
      # rgl::shade3d(clipped)
      # rgl::axes3d()
      # 
      # rgl::shade3d(scaled)
      # rgl::axes3d()
      # rgl::points3d(conv_hull, col = "green")
      
      pnts <- Morpho::vert2points(clipped)
      
      centroid <- apply(pnts, 2, mean) %>%
        matrix(nrow = 1)
      dist_scale <- proxy::dist(centroid, pnts) %>%
        as.vector() %>%
        max()
      
      scaled <- clipped %>%
        Morpho::scalemesh(1 / dist_scale, center = "none")
      
      pnts <- Morpho::vert2points(scaled)
      
      conv_hull <- pnts[unique(as.vector(convhulln(pnts))),]
        
      Rvcg::vcgObjWrite(scaled, new_beak_file)
      readr::write_csv(conv_hull %>% as.data.frame(), hull_file)
        
    }
  }
  return(invisible(NULL))
}

library(furrr)

future::plan(future::multisession)

furrr::future_map(beak_files, ~reorient_beak(.x),
                  .progress = TRUE)


#### Find sphere that fits all meshes in it
library(pbapply)
hulls <- list.files("data/conv_hulls", full.names = TRUE)
all_pnts <- pbapply::pblapply(hulls, readr::read_csv) %>%
  dplyr::bind_rows() %>%
  as.matrix()

centroid <- apply(all_pnts, 2, mean)
all_pnts_cent <- t(t(all_pnts) - centroid)

norms <- apply(all_pnts_cent, 1, function(x) sqrt(sum(x^2)))
max_norm <- max(norms)

## sphere containing all meshes has radius 1.5 and centred at 0.616886145, -0.001077347,  0.145010046

##### Generate SDF samples for all meshes
library(rtorch)