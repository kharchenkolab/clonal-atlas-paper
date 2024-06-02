# Code for the paper "Clonal analysis of murine development reveals novel positional programs directing lineage progression"
A repository with code and analysis for the paper *Erickson, Isaev et al*. AnnData objects and plasmid maps are available on [Zenodo](https://zenodo.org/records/11406618). Python package `scLiTr` that was used for the analysis is available on [GitHub](https://github.com/kharchenkolab/scLiTr). Also, you can find interactive Cellxgene web-explorer on Adameyko lab private [Cellxgene browser](https://adameykolab.hifo.meduniwien.ac.at/cellxgene_public/). Modified TREX algorithm ([original](https://github.com/frisen-lab/TREX)) for clonal barcodes calling is also stored on [GitHub](https://github.com/serjisa/TREX.modified).

## Folders description
- `00_QC` contains notebooks for the quality control of each sample with corresponding pdf-reports (subfolder `QC_reports`),
- `01_Embedding_construction` contains notebooks in which we're building gene expression embeddings and annotate them *before adding the layer with clonal information*,
- `02_Fate_distribution` contains clonal analysis and further hypothesis testing,
- `03_Perturbations` contains a functional and compositional analysis of the results of mosaic perturbations at single-cell and whole clone levels,
- `Supplementary` contains some of the additional analysis for the paper,
- `tools` contains additional functions that are used in Jupyter notebooks.

Please don't hesitate to contact me via `sergey.isaev[at]meduniwien.ac.at` if you have any questions regarding the code or gene expression / clonal objects from the paper.
