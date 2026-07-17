# extended_UMAP_additional-test

Some results obtained using UMAP for R (https://cran.r-project.org/web/packages/umap/) to get the planar depiction of various things extended with the additional data points (https://keldysh.ru/ngpm/2025/proc.pdf).

1. 3D structures of DRT proteins, https://www.rcsb.org/structure/9Z6Y (accessed on 14 of july 2026) and https://www.rcsb.org/structure/9Z6Z (accessed on 16 of july 2026) [DOI: https://doi.org/10.1126/science.aed1656].
  
    Observations:
   
       - Results of the extended UMAP seem to be more refined.

       - Resulting 2D diagrams seem to be mirrored and rotated one against the other. 

Given the second observation, it is relatively easy to align the figures to make their comparative analysis more manageable:

      - Flip one of them LR

      - Rotate until the point clouds will generally match

The results are given in Figure 9Z6Y_9Z6Z_aligned_afterUMAP.png, imperfect alignment here is because of the global changes in structure: 9Z6Z has one longer and one shorter axes.
