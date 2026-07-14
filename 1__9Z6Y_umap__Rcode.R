library(Rpdb)
library(tidyverse)
library(plotly)
library(umap)
library(igraph)

## Input
data_raw <- read.pdb(".../additional_examples/DRT/9Z6Y.pdb")

## Process the data
data <- tibble(eleid = data_raw$atoms$eleid, recname = data_raw$atoms$recname, elename = data_raw$atoms$elename, resid = data_raw$atoms$resid,
				chainid = data_raw$atoms$chainid, x1 = data_raw$atoms$x1, x2 = data_raw$atoms$x2, x3 = data_raw$atoms$x3) |>
			mutate(x1 = round(x1), x2 = round(x2), x3 = round(x3))
## Prepare the extended data
# Find the center
center_1 <- data |> pull(x1) |> mean()
center_2 <- data |> pull(x2) |> mean()
center_3 <- data |> pull(x3) |> mean()
# Calculate radius
radius <- data |> mutate(radius = sqrt((center_1 - x1)^2) + sqrt((center_2 - x2)^2) + sqrt((center_2 - x2)^2)) |> pull(radius) |> max() |> ceiling() * 1.5
# Generate points, SEE: https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume for example
start_n <- data |> pull(eleid) |> max()
ext_coords <- sample_sphere_volume(dim = 3, n = 50000, radius = radius) |> t() |> as.data.frame()
extension_raw <- tibble(eleid = seq(1 : 50000) + start_n,
						recname = rep("ext", 50000), elename = rep("ext", 50000), resid = rep(NA, 50000),
						chainid = rep("extension", 50000), x1 = ext_coords$V1, x2 = ext_coords$V2, x3 = ext_coords$V3) |>
						rowwise() |>
						mutate(phi = runif(n = 1, min = 0, max = 6.18))  |>
						mutate(cost = runif(n = 1, min = -1, max = 1))   |>
						mutate(u_thing = runif(n = 1, min = 0, max = 1)) |>
						mutate(t = acos(cost), r = radius * u_thing^1/3) |>
						mutate( x1 = (r*sin(t)*cos(phi)), x2 = (r*sin(t)*sin(phi)), x3 = r*cos(t) ) |>
						ungroup() |>
						mutate(x1 = round(x1), x2 = round(x2), x3 = round(x3))
center_ext_1 <- extension_raw |> pull(x1) |> mean()
center_ext_2 <- extension_raw |> pull(x2) |> mean()
center_ext_3 <- extension_raw |> pull(x3) |> mean()
delta_1 <- center_1 - center_ext_1
delta_2 <- center_2 - center_ext_2
delta_3 <- center_3 - center_ext_3
extension <- extension_raw |> mutate(x1 = x1 + delta_1, x2 = x2 + delta_2, x3 = x3 + delta_3)

extension <- extension |> anti_join(data, by = c("x1"="x1", "x2"="x2", "x3"="x3")) |> slice_sample(prop = .5)
data_ext <- bind_rows(data, extension)
# Visualize in 3D
#plot_ly(x=extension$x1, y=extension$x2, z=extension$x3, type="scatter3d", size = .1)
#data_sample <- data |> slice_sample(prop = 0.1)
#plot_ly(x=data_sample$x1, y=data_sample$x2, z=data_sample$x3, type="scatter3d", size = .1, color=data_sample$chainid)
#plot_ly(x=data_ext$x1, y=data_ext$x2, z=data_ext$x3, type="scatter3d", size = .1, color=data_ext$chainid)

## UMAP the data
set.seed(1235)
protein_umap <- umap(data |> select(x1, x2, x3), preserve.seed = TRUE, n_neighbors = 15)
data$umap_1 <- protein_umap$layout[,1]
data$umap_2 <- protein_umap$layout[,2]

set.seed(1235)
protein_umap_ext <- umap(data_ext |> select(x1, x2, x3), preserve.seed = TRUE, n_neighbors = 15)
data_ext$umap_1 <- protein_umap_ext$layout[,1]
data_ext$umap_2 <- protein_umap_ext$layout[,2]

## Plot the data
umap_plot_basic <- ggplot(data, aes(x=umap_1, y=umap_2, color=chainid)) +
	geom_point(size = .2) +
	coord_fixed() +
	theme_classic()
umap_plot_basic
umap_plot_ext <- ggplot(data_ext, aes(x=umap_1, y=umap_2, color=chainid)) +
	geom_point(size = .2) +
	coord_fixed() +
	theme_classic()
umap_plot_ext
umap_plot_ext_no <- ggplot(data_ext |> filter(chainid != "extension"), aes(x=umap_1, y=umap_2, color=chainid)) +
	geom_point(size = .2) +
	coord_fixed() +
	theme_classic()
umap_plot_ext_no

## Save the results
ggsave(".../additional_examples/DRT/9Z6Y_umap_extended.png", plot = umap_plot_basic, scale = .7, width = 7, dpi = 300)
ggsave(".../additional_examples/DRT/9Z6Y_umap_extended_w.png", plot = umap_plot_ext, scale = .7, width = 7, dpi = 300)
ggsave(".../additional_examples/DRT/9Z6Y_umap_extended_wo.png", plot = umap_plot_ext_no, scale = .7, width = 7, dpi = 300)