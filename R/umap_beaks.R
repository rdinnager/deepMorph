library(reticulate)
library(readr)
library(dplyr)
library(ggplot2)
library(directlabels)
library(plotly)
library(purrr)
library(rgl)

np <- import("numpy")
umap <- import("umap")

meta <- readr::read_csv("data/resource.csv")
taxo <- readr::read_csv("data/Jetz/2012-03-04206D-master_taxonomy.csv")
file_names <- readr::read_lines("Marian_results/filenames.txt")

latent_codes <- np$load("C:/Users/rdinn/Desktop/latent_codes.npy")

meta_taxo <- meta %>%
  dplyr::left_join(taxo,
                   by = c("Binomal_Jetz" = "TipLabel"))

fams <- dplyr::tibble(ID = gsub(".obj", "", file_names)) %>%
  dplyr::left_join(meta_taxo %>%
                     dplyr::select(ID, family = BLFamilyEnglish)) %>%
  dplyr::mutate(fam_num = as.numeric(as.factor(family)))

umapper <- umap$UMAP(min_dist = 0.75)

umap_latent <- umapper$fit_transform(latent_codes, fams$fam_num)
plot(umap_latent)

dat <- meta_taxo %>%
  dplyr::right_join(dplyr::tibble(ID = gsub(".obj", "", file_names),
                                  X = umap_latent[ , 1],
                                  Y = umap_latent[ , 2]))

unique(dat$BLFamilyEnglish)

# pdf("figures/test_umap.pdf", height = 30, width = 30)
# 
# ggplot(dat, aes(X, Y)) +
#   geom_point(aes(colour = BLFamilyEnglish), alpha = 0.7, size = 6) +
#   geom_dl(aes(label = BLFamilyEnglish, colour = BLFamilyEnglish), method = "ahull.grid") +
#   scale_color_discrete(guide = "none") +
#   theme_minimal()
# 
# dev.off()

p <- ggplot(dat, aes(X, Y)) +
  geom_point(aes(colour = BLFamilyEnglish), alpha = 0.5) +
  scale_color_discrete(guide = "none") +
  theme_minimal()

pl <- plotly::ggplotly(p)
pl

dat <- dat %>%
  dplyr::mutate(latent_code = purrr::map(1:nrow(latent_codes), ~ latent_codes[.x, ]))


source("R/sdf_tools.R")
setup_SDF()

test <- get_meshes_from_latent(dat$latent_code[[4]], show = TRUE)

dat <- dat %>%
  dplyr::mutate(recon_mesh = pbapply::pblapply(latent_code, get_meshes_from_latent))


