library(Morpho)
library(rgl)

test_mesh <- readobj::read.obj("data/objs/Abeillia_abeillei_1.obj", convert.rgl = TRUE)[[1]]
test_landmarks <- readr::read_delim("data/landmarks/Markup_1/Abeillia_abeillei_1.txt", " ",
                                    col_names = c("X", "Y", "Z"))
main_landmarks <- test_landmarks[1:4, ] %>%
  as.matrix()

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")

mid_beak <- apply(main_landmarks[3:4, ], 2, mean) %>% matrix(nrow = 1)
line_1 <- rbind(mid_beak, main_landmarks[1, ])

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")

line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / dist(line_1)

top_line <- main_landmarks[2, ] - mid_beak

tt <- top_line %*% t(line_1_vec)

point_on_line_1 <- mid_beak + as.vector(tt) * line_1_vec

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")
rgl::points3d(point_on_line_1, size = 5, col = "blue")

landmarks_3 <- rbind(main_landmarks[1, ], point_on_line_1, main_landmarks[2, ])
landmark_dists <- dist(landmarks_3) %>% as.matrix()
target_3 <- matrix(c(landmark_dists[1, 2], 0, 0,
                     0, 0, 0,
                     0, 0, landmark_dists[2, 3]),
                   nrow = 3, byrow = TRUE)

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")
rgl::points3d(point_on_line_1, size = 5, col = "blue")
rgl::lines3d(target_3, col = "red")

test_rot <- Morpho::rotmesh.onto(test_mesh, landmarks_3, target_3, adnormals = TRUE)

rgl::shade3d(test_rot$mesh, alpha = 0.1)
rgl::lines3d(target_3, col = "red")
