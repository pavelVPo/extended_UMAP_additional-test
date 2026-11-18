library(Rpdb)
library(tidyverse)
library(plotly)
library(umap)
library(igraph)
library(uniformly)
library(patchwork)
# PDB-files, snapshots of the MD simulation are from supplementaries here:
# https://www.science.org/doi/abs/10.1126/science.1187409?referrer=https%3A%2F%2Fwww.google.com%2F
# files were ordered according to the comments

## Input
input      <- ".../BPTI/input/"
output     <- ".../BPTI/output/"
#rec_1      <- read.pdb(str_glue("{input}1__1187409s_pdb1.pdb"))
rec_2      <- read.pdb(str_glue("{input}2__1187409s_pdb4.pdb"))
rec_3      <- read.pdb(str_glue("{input}3__1187409s_pdb5.pdb"))
rec_4      <- read.pdb(str_glue("{input}4__1187409s_pdb3.pdb"))
rec_5      <- read.pdb(str_glue("{input}5__1187409s_pdb2.pdb"))

## Processing, NB: the following code is OK on condition that all snapshots contains all atoms in the same order
data <- tibble(resname = rec_2$atoms$resname,
					   eleid_2 = rec_2$atoms$eleid, eleid_3 = rec_3$atoms$eleid, eleid_4 = rec_4$atoms$eleid, eleid_5 = rec_5$atoms$eleid,
					   elename_2 = rec_2$atoms$elename, elename_3 = rec_3$atoms$elename, elename_4 = rec_4$atoms$elename, elename_5 = rec_5$atoms$elename,
					   resid_2 = rec_2$atoms$resid, resid_3 = rec_3$atoms$resid, resid_4 = rec_4$atoms$resid, resid_5 = rec_5$atoms$resid,
						V4 = rec_2$atoms$x1, V5 = rec_2$atoms$x2, V6 = rec_2$atoms$x3,
						V7 = rec_3$atoms$x1, V8 = rec_3$atoms$x2, V9 = rec_3$atoms$x3,
						V10 = rec_4$atoms$x1, V11 = rec_4$atoms$x2, V12 = rec_4$atoms$x3,
						V13 = rec_5$atoms$x1, V14 = rec_5$atoms$x2, V15 = rec_5$atoms$x3) |>
			mutate(across(starts_with('V'), ~scale(., center = TRUE, scale = FALSE))) |>
			mutate(across(starts_with('V'), ~round(.)))
# Check for problems
n_problems <- data |> filter(eleid_2 != eleid_3 | eleid_2 != eleid_4 | eleid_2 != eleid_5) |> pull(eleid_2) |> length()
# Find radius of the sphere holding points in its volume
radius <- (( data |> select(starts_with('V')) |> as.matrix() |> max() - data |> select(starts_with('V')) |> as.matrix() |> min() ) / 2) |> ceiling() + 2
# Sample points from p-(15)-cube,
# Uniformly library will be used here, perhaps not really uniform, tbh
by <- join_by(V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15)
sample_raw <- runif_in_cube(n = 50000000, d = 12, O = rep(0, 12), r = radius*3) |>
					ceiling()   |>
					as_tibble() |>
					rename(V4 = V1, V5 = V2, V6 = V3, V7 = V4, V8 = V5, V9 = V6, V10 = V7, V11 = V8, V12 = V9, V13 = V10, V14 = V11, V15 = V12) |>
					distinct()  |>
					anti_join(data, by) |>
					slice_sample(n = 2000)

extension <- sample_raw |> mutate(resname = "space", eleid_2 = 10^10, resid_2 = 10^10) 
data_ext  <- bind_rows(extension |> select(resname, eleid_2, resid_2, starts_with("V")), data |> select(resname, eleid_2, resid_2, starts_with("V")))
# Visualize in 3D to check things
plot_ly(x=data_ext$V13, y=data_ext$V14, z=data_ext$V15, type="scatter3d", size = .025, color=data_ext$resname)

## UMAP the data
neighbors <- 15
set.seed(1235)
#umap_1 <- umap(data_ext |> select(V1, V2, V3), preserve.seed = TRUE, n_neighbors = neighbors)
#data_ext$umap1_1 <- umap_1$layout[,1]
#data_ext$umap1_2 <- umap_1$layout[,2]
umap_2 <- umap(data_ext |> select(V4, V5, V6), preserve.seed = TRUE, n_neighbors = neighbors)
data_ext$umap2_1 <- umap_2$layout[,1]
data_ext$umap2_2 <- umap_2$layout[,2]
umap_3 <- umap(data_ext |> select(V7, V8, V9), preserve.seed = TRUE, n_neighbors = neighbors)
data_ext$umap3_1 <- umap_3$layout[,1]
data_ext$umap3_2 <- umap_3$layout[,2]
umap_4 <- umap(data_ext |> select(V10, V11, V12), preserve.seed = TRUE, n_neighbors = neighbors)
data_ext$umap4_1 <- umap_4$layout[,1]
data_ext$umap4_2 <- umap_4$layout[,2]
umap_5 <- umap(data_ext |> select(V13, V14, V15), preserve.seed = TRUE, n_neighbors = neighbors)
data_ext$umap5_1 <- umap_5$layout[,1]
data_ext$umap5_2 <- umap_5$layout[,2]

umap_dynamic    <- umap(data_ext |> select(starts_with('V')), preserve.seed = TRUE, n_neighbors = neighbors)
data_ext$umap_1 <- umap_dynamic$layout[,1]
data_ext$umap_2 <- umap_dynamic$layout[,2]

umap_dynamic_3d    <- umap(data_ext |> select(starts_with('V')), preserve.seed = TRUE, n_neighbors = neighbors, n_components = 3)
data_ext$umap_x    <- umap_dynamic_3d$layout[,1]
data_ext$umap_y    <- umap_dynamic_3d$layout[,2]
data_ext$umap_z    <- umap_dynamic_3d$layout[,3]
data_ext_no <- data_ext |> filter(resname != "space")
plot_ly(x=data_ext_no$umap_x, y=data_ext_no$umap_y, z=data_ext_no$umap_z,
			type="scatter3d", size = .5, color=data_ext_no$resid_2,
			colors = c("#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8"))

## Plot the data
plot_umap_dynamic <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap_1, y=umap_2)) +
						geom_point(aes(color=factor(resid_2)), size = 1, alpha = 1) +
						scale_colour_manual( values = c(
						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
						coord_fixed() +
						ggtitle("extended UMAP procedure, BPTI folding", subtitle = str_glue("n_neighbors = {neighbors}")) +
						theme_classic() +
						theme(legend.position="none")
plot_umap_dynamic

# plot_umap_1  <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap1_1, y=umap1_2, color=resid_1)) +
# 					geom_point(aes(color=factor(resid_1)), size = 1, alpha = 1) +
# 					scale_colour_manual( values = c(
# 						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
# 						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
# 						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
# 						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
# 						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
# 						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
# 						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
# 					coord_fixed() +
# 					theme_void() +
# 					theme(legend.position="none") +
# 					ggtitle("1", subtitle = waiver())
plot_umap_2  <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap2_1, y=umap2_2, color=resname)) +
					geom_point(aes(color=factor(resid_2)), size = 1, alpha = 1) +
					scale_colour_manual( values = c(
						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
					coord_fixed() +
					theme_void() +
					theme(legend.position="none") +
					ggtitle("2", subtitle = waiver())
plot_umap_3  <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap3_1, y=umap3_2, color=resname)) +
					geom_point(aes(color=factor(resid_2)), size = 1, alpha = 1) +
					scale_colour_manual( values = c(
						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
					coord_fixed() +
					theme_void() +
					theme(legend.position="none") +
					ggtitle("3", subtitle = waiver())
plot_umap_4  <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap4_1, y=umap4_2, color=resname)) +
					geom_point(aes(color=factor(resid_2)), size = 1, alpha = 1) +
					scale_colour_manual( values = c(
						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
					coord_fixed() +
					theme_void() +
					theme(legend.position="none") +
					ggtitle("4", subtitle = waiver())
plot_umap_5  <- ggplot(data_ext |> filter(resname != "space"), aes(x=umap5_1, y=umap5_2, color=resname)) +
					geom_point(aes(color=factor(resid_2)), size = 1, alpha = 1) +
					scale_colour_manual( values = c(
						"#b71c1c", "#c62828", "#d32f2f", "#e53935", "#ef5350", "#e57373", "#ef9a9a", "#ffcdd2", "#ffebee",
						"#e8eaf6", "#c5cae9", "#9fa8da", "#7986cb", "#5c6bc0", "#3f51b5", "#3949ab", "#303f9f", "#1a237e",
						"#0d47a1", "#1565c0", "#1976d2", "#1e88e5", "#2196f3", "#42a5f5", "#64b5f6", "#90caf9", "#bbdefb",
						"#827717", "#9e9d24", "#c0ca33", "#d4e157", "#dce775", "#e6ee9c", "#f0f4c3", "#f9fbe7",
						"grey20","grey30","grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "black",
						"#fff3e0", "#ffe0b2", "#ffcc80", "#ffb74d", "#ffa726", "#ff9800", "#fb8c00", "#f57c00", "#e65100",
						"#795548", "#8d6e63", "#a1887f", "#bcaaa4", "#d7ccc8")) +
					coord_fixed() +
					theme_void() +
					theme(legend.position="none") +
					ggtitle("5", subtitle = waiver())
plots <- plot_umap_2 | plot_umap_3 | plot_umap_4 | plot_umap_5
plots

## Save the results
ggsave(str_glue("{output}BPTI_UMAP_folding_series__wo1.png"), plot = plots, width = 7, height = 3, units = "in", dpi = 300)
ggsave(str_glue("{output}BPTI_UMAP_dynamics__wo1.png"), plot = plot_umap_dynamic, width = 5, height = 7, units = "in", dpi = 300)