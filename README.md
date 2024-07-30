# Code for the paper "Clonal analysis reveals spatial determinants of progenitor fates in embryos"
A repository with code and analysis for the paper *Erickson, Isaev et al*.

## Folders description
- `00_QC` contains notebooks for the quality control of each sample with corresponding pdf-reports (subfolder `QC_reports`),
- `01_Embedding_construction` contains notebooks in which we're building gene expression embeddings and annotate them *before adding the layer with clonal information*,
- `02_Fate_distribution` contains clonal analysis and further hypothesis testing,
- `03_Perturbations` contains a functional and compositional analysis of the results of mosaic perturbations at single-cell and whole clone levels,
- `Supplementary` contains some of the additional analysis for the paper,
- `tools` contains additional functions that are used in Jupyter notebooks.

## Additional data sources
* Raw reads and Cell Ranger outputs are reachable via [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269395),
* Modified TREX algorithm ([original](https://github.com/frisen-lab/TREX)) for clonal barcodes calling is stored on [GitHub](https://github.com/serjisa/TREX.modified),
* AnnData objects (both clonal and gene expression) and plasmid maps (for chimeric references construction) are available on [Zenodo](https://zenodo.org/records/11406618),
* Python package `scLiTr` that was used for the analysis is available on [GitHub](https://github.com/kharchenkolab/scLiTr),
* Interactive Cellxgene web-explorer is located on on Adameyko lab private [Cellxgene browser](https://adameykolab.hifo.meduniwien.ac.at/cellxgene_public/filecrawl/.review/.2024_Erickson_Isaev),
* Interactive web-app *clones2cells* for simultanious clonal and gene embedding exploration is available on [Streamlit Cloud](https://clones2cells.streamlit.app).

Please don't hesitate to contact me via `sergey.isaev[at]meduniwien.ac.at` if you have any questions regarding the code or gene expression / clonal objects from the paper.
