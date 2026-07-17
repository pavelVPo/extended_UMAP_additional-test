library(tidyverse)
library(magick)
library(patchwork)

## Input
data_9z6y_raw <- read_tsv("C:/.../additional_examples/DRT/9Z6Y_data.tsv") |> filter(chainid != "extension")
data_9z6z_raw <- read_tsv("C:/.../additional_examples/DRT/9Z6Z_data.tsv") |> filter(chainid != "extension")

## Process
# Scale & Center & Flip LR one diagram & Plot to check 
coords_9z6y <- data_9z6y_raw |> select(umap_1, umap_2) |> as.matrix() |> scale(center = c(1,2), scale = c(1,2)) |>
								as.data.frame() |>
								mutate(chainid = data_9z6y_raw$chainid)
coords_9z6z <- data_9z6z_raw |> select(umap_1, umap_2) |> as.matrix() |> scale(center = c(1,2), scale = c(1,2)) |>
								as.data.frame() |>
								mutate(chainid = data_9z6z_raw$chainid) |>
								mutate(umap_1 = -umap_1)
# Plot one diagram
par(pty="s")
png(file = "C:/.../additional_examples/DRT/pics/base_9Z6Y.png", width = 800, height = 800, res = 100)
plot(coords_9z6y$umap_1, coords_9z6y$umap_2, xlim = c(-15, 15), ylim = c(-15, 15), cex = .005, col = as.factor(coords_9z6y$chainid), axes = FALSE, ann = FALSE, bty = "n", asp = 1)
dev.off()
# Read the results, blur the pic and get the data to matrix
dia_one_img <- image_read("C:/.../additional_examples/DRT/pics/base_9Z6Y.png") |> image_blur(radius = 15, sigma = 5)
dia_m <- as.integer(image_data(dia_one_img, channels = "gray"))  |> as.vector()
# Rotate (up to 360 by .1 degree) and plot the other diagram to compare with the first one to find the transformation well aligning them
# Prepare the dataframe to store the results
transforms <- data.frame(pi_theta = rep(NA, 420), sum_delta = rep(NA, 420))
# Set the angle
pi_theta <- 0
for (i in seq(1:100)) {
	pi_theta <- pi_theta + pi/100
	coords_9z6z_try <- data.frame(umap_1 = rep(NA, nrow(coords_9z6z)), umap_2 = rep(NA, nrow(coords_9z6z)), chainid = data_9z6z_raw$chainid)
	coords_9z6z_try$umap_1 <- coords_9z6z |> select(umap_1, umap_2) |> mutate(umap_1 = umap_1*cos(pi_theta) - umap_2*sin(pi_theta)) |> pull(umap_1)
	coords_9z6z_try$umap_2 <- coords_9z6z |> select(umap_1, umap_2) |> mutate(umap_2 = umap_1*sin(pi_theta) + umap_2*cos(pi_theta)) |> pull(umap_2)
	par(pty="s")
	path = str_glue("C:/.../additional_examples/DRT/pics/{i}_9Z6Z_temp.png") 
	png(path, width = 800, height = 800, res = 100)
	plot(coords_9z6z_try$umap_1, coords_9z6z_try$umap_2, xlim = c(-15, 15), ylim = c(-15, 15), cex = .005, col = as.factor(coords_9z6z_try$chainid), axes = FALSE, ann = FALSE, bty = "n", asp = 1)
	dev.off()
	dia_temp_img <- image_read(path) |> image_convert(type = "grayscale") |> image_blur(radius = 15, sigma = 5)
	dia_temp_m <- as.integer(image_data(dia_temp_img, channels = "gray")) |> as.vector()
	delta_m    <- abs(dia_m - dia_temp_m)
	delta 	   <- sum(delta_m)
	transforms[i,1] <- pi_theta
	transforms[i,2] <- delta
}



## Align
transforms <- transforms |> arrange(sum_delta)
transform  <- transforms[1,1]
coords_9z6z_new <- coords_9z6z |> mutate(umap_1 = umap_1*cos(transform) - umap_2*sin(transform), umap_2 = umap_1*sin(transform) + umap_2*cos(transform))
data_9z6y <- data_9z6y_raw
data_9z6z <- data_9z6z_raw
data_9z6y$umap_1 <- coords_9z6y$umap_1
data_9z6y$umap_2 <- coords_9z6y$umap_2
data_9z6z$umap_1 <- coords_9z6z_new$umap_1
data_9z6z$umap_2 <- coords_9z6z_new$umap_2

## Plot
umap_9z6y <- ggplot(data_9z6y, aes(x=umap_1, y=umap_2, color=chainid)) +
	geom_point(size = .025, alpha = .4) +
	coord_fixed(xlim = c(-15,15), ylim = c(-15,15)) +
	theme_void() +
	ggtitle("9Z6Y after extended UMAP procedure", subtitle = "n_neighbors = 15, additional points not shown") +
	theme(legend.position = "none") +
	theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 4))
umap_9z6z <- ggplot(data_9z6z, aes(x=umap_1, y=umap_2, color=chainid)) +
	geom_point(size = .025, alpha = .4) +
	coord_fixed(xlim = c(-15,15), ylim = c(-15,15)) +
	theme_void() +
	ggtitle("9Z6Z after extended UMAP procedure", subtitle = "n_neighbors = 15, additional points not shown,\nflipped, rotated") +
	theme(legend.position = "none") +
	theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 4))
plot <- umap_9z6y + umap_9z6z
plot
ggsave("C:/.../additional_examples/DRT/9Z6Y_9Z6Z_aligned_afterUMAP.png", plot = plot, width = 7, units = "in", dpi = 300)