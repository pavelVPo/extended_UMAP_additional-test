library(Rpdb)
library(tidyverse)
library(plotly)
library(umap)

## Input
data_raw <- read.pdb(".../!_extendedUMAP/additional_examples/DRT/9Z6Y.pdb")

## Process the data
data <- tibble(eleid = data_raw$atoms$eleid, recname = data_raw$atoms$recname, elename = data_raw$atoms$elename, resid = data_raw$atoms$resid,
				chainid = data_raw$atoms$chainid, x1 = data_raw$atoms$x1, x2 = data_raw$atoms$x2, x3 = data_raw$atoms$x3) |>
			mutate(x1 = round(x1), x2 = round(x2), x3 = round(x3))
## Prepare the extended data
# Find the center
center_1 <- data |> pull(x1) |> mean()
center_2 <- data |> pull(x2) |> mean()
center_3 <- data |> pull(x3) |> mean()
# Calculate the radius
radius <- data |> mutate(radius = sqrt((center_1 - x1)^2) + sqrt((center_2 - x2)^2) + sqrt((center_2 - x2)^2)) |> pull(radius) |> max() |> ceiling()
# Generate points, SEE: https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume for example
start_n <- data |> pull(eleid) |> max()
extension <- tibble(eleid = seq(1 : 50000) + start_n,
						recname = rep("ext", 50000), elename = rep("ext", 50000), resid = rep(NA, 50000),
						chainid = rep("extension", 50000), x1 = rep(NA, 50000), x2 = rep(NA, 50000), x3 = rep(NA, 50000)) |>
			 rowwise() |>
			 mutate(phi = runif(n = 1, min = 0, max = 6.18))  |>
			 mutate(cost = runif(n = 1, min = -1, max = 1))   |>
			 mutate(u_thing = runif(n = 1, min = 0, max = 1)) |>
			 mutate(t = acos(cost), r = radius * u_thing^1/3) |>
			 mutate( x1 = center_1 + (r*sin(t)*cos(phi)), x2 = center_2 + (r*sin(t)*sin(phi)), x3 = r*cos(t) ) |>
			 ungroup() |>
			 mutate(x1 = round(x1), x2 = round(x2), x3 = round(x3))

extension <- extension |> anti_join(data, by = c("x1"="x1", "x2"="x2", "x3"="x3")) |> slice_sample(prop = .5)
data_ext <- bind_rows(data, extension)
# Visualize in 3D
#plot_ly(x=extension$x1, y=extension$x2, z=extension$x3, type="scatter3d", size = .1)
#data_sample <- data |> slice_sample(prop = 0.1)
#plot_ly(x=data_sample$x1, y=data_sample$x2, z=data_sample$x3, type="scatter3d", size = .1, color=data_sample$chainid)

# UMAP the data
set.seed(1235)
protein_umap <- umap(data |> select(x1, x2, x3), preserve.seed = TRUE, n_neighbors = 15)
data$umap_1 <- protein_umap$layout[,1]
data$umap_2 <- protein_umap$layout[,2]

set.seed(1235)
protein_umap_ext <- umap(data_ext |> select(x1, x2, x3), preserve.seed = TRUE, n_neighbors = 15)
data_ext$umap_1 <- protein_umap_ext$layout[,1]
data_ext$umap_2 <- protein_umap_ext$layout[,2]

# Plot the data
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
