library(reticulate)
library(readr)
library(dplyr)

np <- import("numpy")
umap <- import("umap")

meta <- readr::read_csv("data/resource.csv")
taxo <- readr::read_csv("data/Jetz/2012-03-04206D-master_taxonomy.csv")

latent_codes <- np$load("C:/Users/rdinn/Desktop/latent_codes.npy")

meta_taxo <- meta %>%
  dplyr::left_join(taxo,
                   by = c("Binomal_Jetz" = "TipLabel"))

umapper <- umap$UMAP()

umap_latent <- umapper$fit_transform(latent_codes)
plot(umap_latent)

umapper <- umap$UMAP(min_dist = 0.01)

umap_latent <- umapper$fit_transform(latent_codes)
plot(umap_latent)


names_2020 <- list.files("../data/landmarks_oriented/L_landmarks/")

dat <- meta_taxo %>%
  dplyr::mutate(X = umap_latent[ , 1],
                Y = umap_latent[ , 2])