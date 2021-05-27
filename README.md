
# overview

<p class="comment">
Raw 10x data were processed with CellRanger pipeline resulting count matrix is analyzed here. [empty Drops method](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) used to distinguish between background (droplets with no RNA) from droplets containing RNA. [Doublet cells](https://ltla.github.io/SingleCellThoughts/software/doublet_detection/bycell.html) were identified through simulation (adding 2 random scRNAseq signature together) - then recluster cells - and count synthetic doublets in neighborhood of each cell. 
We used [Pearson residuals derived from a generalized negative binomial model of UMI counts](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) as input for principal component analysis and differential gene expression analysis. Finally, we corrected for cell cycle influence on gene expression calculating a score based on of [43 genes primarily expressed in G1/S and 55 primarily expressed in G2/M](https://science.sciencemag.org/content/352/6282/189.long).
</p>
