library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)

beak_files <- list.files("data/objs/", full.names = TRUE) 

#beak_file <- beak_files[33]
#plot_folder <- "beak_plots_trimmed"
reorient_beak <- function(beak_file, cut_extraneous = TRUE, mesh_folder = "data/oriented_trimmed_meshes_obj",
                          landmark_folder = "data/landmarks_oriented", plot_folder = NULL) {
  
  if(!dir.exists(mesh_folder)) {
    dir.create(mesh_folder)
  }
  if(!dir.exists(landmark_folder)) {
    dir.create(landmark_folder)
  }
  
  bird_file <- basename(beak_file)
  
  if(!dir.exists(mesh_folder)) {
    dir.create(mesh_folder)
  }
  
  new_beak_file <- file.path(mesh_folder, bird_file)
  
  if(!is.null(plot_folder)) {
    if(!dir.exists(plot_folder)) {
      dir.create(plot_folder)
    }
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
      
      line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / dist(line_1)
      
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
      
      if(cut_extraneous) {
        
        beak_width <- abs(rot_landmarks[4, 2] - rot_landmarks[3, 2])
        addit <- 0.01 * beak_width
        
        ## plane made up of hind beak landmarks
        cut_plane_1 <- rot_landmarks[c(2, 4, 3), ] %>% as.matrix()
        cut_1 <- Morpho::cutMeshPlane(rot$mesh, v1 = cut_plane_1[1, ], v2 = cut_plane_1[2, ], v3 = cut_plane_1[3, ])
        
        ## plane cutting up-down aligned with hindmost landmark
        hindmost <- rot_landmarks[which.min(rot_landmarks[ , 1]), , drop = FALSE]
        cut_plane_2 <- rbind(hindmost - c(0, 1, 0),
                                   hindmost - c(0, 0, 1),
                                   hindmost - c(0, -1, 0))
        
        ## approx 30 degree plane cutting towards beak tip and downwards from hind landmarks
        cut_plane_3 <- rbind(rot_landmarks[c(3, 4), ], 
                             c(cos(0.5) + mean(rot_landmarks[c(3, 4), 1]), 0, -1))
        ## horizontal plane aligned with topmost landmark
        topmost <- rot_landmarks[which.max(rot_landmarks[ , 3]), , drop = FALSE]
        cut_plane_4 <- rbind(topmost - c(-1, 0, -addit),
                             topmost - c(0, 1, -addit),
                             topmost - c(1, 0, -addit))
        
        ## vertical plane aligned with leftmost landmark
        leftmost <- rot_landmarks[which.max(rot_landmarks[ , 2]), , drop = FALSE]
        cut_plane_5 <- rbind(leftmost - c(0, -addit, -1),
                             leftmost - c(1, -addit, 0),
                             leftmost - c(0, -addit, 1))
        
        
        ## vertical plane aligned with rightmost landmark
        rightmost <- rot_landmarks[which.min(rot_landmarks[ , 2]), , drop = FALSE]
        cut_plane_6 <- rbind(rightmost - c(0, addit, 1),
                             rightmost - c(1, addit, 0),
                             rightmost - c(0, addit, -1))
        
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_2[1, ], v2 = cut_plane_2[2, ], v3 = cut_plane_2[3, ])
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_3[1, ], v2 = cut_plane_3[2, ], v3 = cut_plane_3[3, ])
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_4[1, ], v2 = cut_plane_4[2, ], v3 = cut_plane_4[3, ])
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_5[1, ], v2 = cut_plane_5[2, ], v3 = cut_plane_5[3, ])
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_6[1, ], v2 = cut_plane_6[2, ], v3 = cut_plane_6[3, ])
        
        cut_1 <- Rvcg::vcgIsolated(cut_1, silent = TRUE, facenum = 0.1*ncol(cut_1$it))
        
        # rgl::shade3d(cut_1, alpha = 0.1)
        # rgl::lines3d(target_3, col = "red")
        # rgl::points3d(rot_landmarks, size = 4, col = "blue")
        
        
        Rvcg::vcgObjWrite(cut_1, new_beak_file)
        
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

for(i in seq_along(beak_files)) {
  reorient_beak(beak_files[i])
  cat("Finished reorienting beak", i, "of", length(beak_files), "->", basename(beak_files[i]), "\n")
}

########### do cuts and rotating ############

for(i in seq_along(beak_files)) {
  reorient_beak(beak_files[i], cut_extraneous = TRUE, mesh_folder = "data/oriented_trimmed_meshes_obj",
                landmark_folder = "data/landmarks_oriented", plot_folder = "beak_plots_trimmed")
  gc()
  cat("Finished reorienting beak", i, "of", length(beak_files), "->", basename(beak_files[i]), "\n")
}
