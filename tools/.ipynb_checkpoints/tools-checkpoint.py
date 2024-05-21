import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import scanpy.external as sce
import seaborn as sns
from scipy.sparse import hstack


def SCneighbors(
    adata,
    obs_name,
    GEX_graph_key=None,
    clonal_graph_key_added="shared_clones",
    return_only_adata=False,
    return_clonal_adata=False,
    tqdm_bar=False,
    **kwargs,
):
    from scipy.sparse import csr_matrix
    if tqdm_bar:
        from tqdm import tqdm
    
    adata.obs[obs_name] = adata.obs[obs_name].astype("category")
    var_mapping = dict(zip(
        adata.obs[obs_name].cat.categories[adata.obs[obs_name].cat.categories != "NA"],
        range(len(adata.obs[obs_name].cat.categories) - 1),
    ))
    if GEX_graph_key is None:
        GEX_graph_key = "distances"
    else:
        GEX_graph_key = GEX_graph_key + "_distances"
    
    col_ind = []
    row_ind = []
    data = []
    
    for i in (tqdm(range(len(adata))) if tqdm_bar else range(len(adata))):
        nn = adata.obs[obs_name][(adata.obsp[GEX_graph_key][i] > 0).A[0]].value_counts()
        nn = nn[(nn.index != "NA") & (nn > 0)]
        col_ind += [var_mapping[var] for var in nn.index]
        row_ind += [i] * len(nn)
        data += list(nn.values)
        
    adata_clonal = sc.AnnData(
        X=csr_matrix((data, (row_ind, col_ind))),
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=list(var_mapping.keys())),
    )
    if return_only_adata:
        return adata_clonal
    else:
        if "metric" not in kwargs:
            kwargs["metric"] = "cosine"
        sc.pp.neighbors(
            adata_clonal,
            use_rep="X",
            key_added=clonal_graph_key_added,
            **kwargs,
        )
        if clonal_graph_key_added is None:
            dist_key = "distances"
            conn_key = "connectivities"
            uns_key = "neighbors"
        else:
            dist_key = clonal_graph_key_added + "_distances"
            conn_key = clonal_graph_key_added + "_connectivities"
            uns_key = clonal_graph_key_added

        adata.obsp[dist_key] = adata_clonal.obsp[dist_key].copy()
        adata.obsp[conn_key] = adata_clonal.obsp[conn_key].copy()
        adata.uns[uns_key] = adata_clonal.uns[uns_key].copy()

        if return_clonal_adata:
            return adata_clonal

def seurat_rpca(adata, batch_key):
    """
    This is a wrapper function around Seurat RPCA integration.
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples / batches.
    Returns
    ----------
    AnnData object with corrected expression values.
    """
    import anndata2ri
    import rpy2.robjects as ro
    from scipy.sparse import csr_matrix

    # libraries loading
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    # ro.r("options(future.globals.maxSize = 891289600)")

    # anndata2ri activation
    anndata2ri.activate()

    adata_tmp = adata.copy()
    if "highly_variable" in adata_tmp.var.columns:
        adata_tmp.raw = adata_tmp
        adata_tmp = adata_tmp[:, adata_tmp.var.highly_variable]

    adata_tmp.X = csr_matrix(adata_tmp.X)

    ro.globalenv["adata_tmp"] = adata_tmp

    # converting to Seurat
    ro.r('options(future.globals.maxSize=10737418240)')
    ro.r('seurat.obj <- as.Seurat(adata_tmp, counts = NULL, data = "X")')
    # feature selection is already performed
    ro.r("features <- rownames(seurat.obj)")
    # splitting dataset per batches
    ro.r(f'seurat.objs <- SplitObject(seurat.obj, split.by = "{batch_key}")')
    # centering, scaling and PCA per batch for RPCA-based integration
    ro.r(
        """seurat.objs <- lapply(X = seurat.objs, FUN = function(x) {
                x <- ScaleData(x, features = features, verbose = FALSE)
                x <- RunPCA(x, features = features, verbose = FALSE)
            })"""
    )
    # integration
    ro.r(
        """anchors <- FindIntegrationAnchors(
                object.list = seurat.objs,
                anchor.features = features,
                reduction = "rpca",
                verbose = FALSE
            )"""
    )
    ro.r("seurat.objs.combined <- IntegrateData(anchorset = anchors, verbose = FALSE)")

    # extraction of corrected expressions
    expression = ro.r(
        'GetAssayData(object = seurat.objs.combined, assay = "integrated", slot = "data")'
    )
    obs_names = ro.r("colnames(seurat.objs.combined)")
    var_names = ro.r("rownames(seurat.objs.combined)")

    adata_tmp = adata_tmp[np.array(obs_names), np.array(var_names)]
    adata_tmp.X = expression.T

    return adata_tmp

black_borders_boxplot = {
    "boxprops": {"edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}

def concatenate_mdatas(mdatas):
    keys = list(mdatas[0].mod.keys())
    mdata = {}
    for key in keys:
        adata = mdatas[0][key]
        adata = adata.concatenate(
            [mdata[key] for mdata in mdatas][1:],
            batch_key=None,
            join="outer",
            index_unique=None,
            fill_value=0
        )
        for column in adata.obs.dtypes.to_dict():
            if adata.obs.dtypes.to_dict()[column] == np.dtype("O"):
                adata.obs[column] = adata.obs[column].astype("category")
        mdata[key] = adata
    mdata = mu.MuData(mdata)
    mdata.update()
    return mdata


def assign_gaps(adata, verbose=True):
    partial_clone_ids = [clone_id for clone_id in adata.var_names if "-" in clone_id]
    full_clone_ids = [clone_id for clone_id in adata.var_names if "-" not in clone_id]
    if verbose:
        print(f"Number of clone IDs with gaps before correction: {len(partial_clone_ids)}")
    recovered_clone_ids = {}
    drop = []
    for clone_id in partial_clone_ids:
        if clone_id[0] == "-":
            clone_id_strip = clone_id.lstrip("-")
            clone_ids_matching = pd.Series(full_clone_ids).str.endswith(clone_id_strip)
            clone_ids_matching = np.array(full_clone_ids)[clone_ids_matching]
            if len(clone_ids_matching) == 1:
                recovered_clone_ids[clone_id] = clone_ids_matching[0]
                drop.append(clone_id)
            elif len(clone_ids_matching) > 1:
                drop.append(clone_id)
        else:
            clone_id_strip = clone_id.rstrip("-")
            clone_ids_matching = pd.Series(full_clone_ids).str.startswith(clone_id_strip)
            clone_ids_matching = np.array(full_clone_ids)[clone_ids_matching]
            if len(clone_ids_matching) == 1:
                recovered_clone_ids[clone_id] = clone_ids_matching[0]
                drop.append(clone_id)
            elif len(clone_ids_matching) > 1:
                drop.append(clone_id)
    for clone_id in recovered_clone_ids:
        adata.X[:, adata.var_names == recovered_clone_ids[clone_id]] += adata.X[:, adata.var_names == clone_id]
    adata = adata[:, ~np.isin(adata.var_names, list(recovered_clone_ids.keys()))]
    if verbose:
        print(f"Number of clone IDs with gaps after correction: {len(partial_clone_ids) - len(recovered_clone_ids)}")
    return adata


def get_count_pseudobulk(adata, layer=None, split_by=None, use_raw=False):
    """
    Returns DataFrame with pseudo-bulks generated based on AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with single cell expressions.
    layer : str, optional
        Name of the layer with count data.
    split_by : str, optional
        Name of column in `adata.obs` with variable to split.
    use_raw : bool, optional
        If it's needed to use `adata.raw`
    Returns
    ----------
    DataFrame with pseudo-bulk.
    
    """
    if split_by is None:
        if use_raw:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.raw.X.sum(axis=0)).A[0]})
        elif layer is None:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.X.sum(axis=0)).A[0]})
        else:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.layers[layer].sum(axis=0)).A[0]})
    else:
        df_res = pd.DataFrame()
        for splitter in set(adata.obs[split_by]):
            if use_raw:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].raw.X.sum(axis=0)).A[0]
            elif layer is None:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].X.sum(axis=0)).A[0]
            else:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].layers[layer].sum(axis=0)).A[0]
    if use_raw:
        df_res.index = adata.raw.var_names
    else:
        df_res.index = adata.var_names
    return 

def get_DESeq(adata, group, label="cell_type_l2", split_by="sample_id", reference=False):
    """
    It's a DESeq2 wrapper for DE testing between two groups.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    group : str
        Name of target group.
    label : str, optional
        Name of column in `adata.obs` with target group labeling.
    split_by : str, optional
        Name of column in `adata.obs` with sample labeling.
    reference : str, optional
        Name of column in `adata.obs` with reference group labeling.
        If False use all non-target groups as reference.
        
    Returns
    ----------
    Dictionary with DESeq2 results.
    """
    import pandas as pd
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, Formula
    pandas2ri.activate()
    from rpy2.robjects.packages import importr
    deseq = importr("DESeq2")
    BiocGenerics = importr("BiocGenerics")
    
    to_dataframe = robjects.r('function(x) data.frame(x)')
    
    class py_DESeq2:
        def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
            try:
                assert gene_column in count_matrix.columns, 'Wrong gene id column name'
                gene_id = count_matrix[gene_column]
            except AttributeError:
                sys.exit('Wrong Pandas dataframe?')

            self.dds = None
            self.deseq_result = None
            self.resLFC = None
            self.comparison = None
            self.normalized_count_matrix = None
            self.gene_column = gene_column
            self.gene_id = count_matrix[self.gene_column]
            self.count_matrix = robjects.conversion.py2rpy(count_matrix.drop(gene_column,axis=1))
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
            self.design_formula = Formula(design_formula)

        def run_deseq(self, **kwargs):
            self.dds = deseq.DESeqDataSetFromMatrix(
                countData=self.count_matrix, 
                colData=self.design_matrix,
                design=self.design_formula
            )
            self.dds = deseq.DESeq(self.dds, **kwargs)
            self.normalized_count_matrix = BiocGenerics.counts(self.dds, normalized=True)

        def get_deseq_result(self, **kwargs):
            self.comparison = deseq.resultsNames(self.dds)
            self.deseq_result = deseq.results(self.dds, **kwargs)
            self.deseq_result = to_dataframe(self.deseq_result)
            self.deseq_result = robjects.conversion.rpy2py(self.deseq_result)
            self.deseq_result[self.gene_column] = self.gene_id.values
    
    expr = get_count_pseudobulk(adata[adata.obs[label] == group], layer="counts", split_by=split_by).astype("int32")
    if not reference:
        expr_ref = get_count_pseudobulk(adata[adata.obs[label] != group], layer="counts", split_by=split_by).astype("int32")
    else:
        expr_ref = get_count_pseudobulk(adata[adata.obs[label] == reference], layer="counts", split_by=split_by).astype("int32")
    expr_ref.columns = expr_ref.columns + "_control"
    expression_matrix = expr.join(expr_ref)
    conditions = pd.DataFrame({"treatment" : ["tagret"] * len(expr.columns) + ["control"] * len(expr_ref.columns)}, index=expression_matrix.columns)
    expression_matrix["id"] = expression_matrix.index
    expression_matrix.index = list(range(len(expression_matrix)))
    
    DESeq = py_DESeq2(expression_matrix, conditions, design_formula="~ treatment")
    DESeq.run_deseq()
    DESeq.get_deseq_result()
    
    return DESeq.deseq_result

def get_GSEA(DESeq_results, genesets, outdir):
    """
    Performs GSEA analysis based on DESeq results.
    
    Parameters
    ----------
    DESeq_results : DataFrame
        Output of get_DESeq() function.
    genesets : dict
        Dictionary with gene sets for GSEA analysis.
    outdir : str
        Path to directory with outputs.
        
    Returns
    ----------
    DataFrame with GSEA results.
    """
    import gseapy as gp
    import pandas as pd
    
    if DESeq_results.columns != ["0", "1"]:
        DESeq_df = DESeq_results.dropna().sort_values(by="log2FoldChange", ascending=False)
        rnk = pd.DataFrame({0 : DESeq_df.id, 1 : DESeq_df.log2FoldChange})
    else:
        rnk = DESeq_results

    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=genesets,
        processes=35,
        permutation_num=100,
        outdir=outdir,
        format="pdf",
        seed=0,
        no_plot=True,
        verbose=True,
    )

    results = pre_res.res2d.copy()
    return results

def annotate_symphony(adata, reference_path, layer=None, batch_key=None, reference_key="cell_type", k=5,
                      clean_pred=None, do_preprocessing=True):
    """
    This is a wrapper function around Symphony label transferring.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    reference_path : str
        Path to an .rds-object with pre-built Symphony reference.
    layer : str, optional
        Name of layer with counts. Default: `None` (counts are expected in `adata.X`)
    batch_key : str, optional
        Name of column in `adata.obs` that corresponds to the name of samples / batches.
    reference_key : str, optional
        Name of target column in Symphony reference. Default: `cell_type`.
    k : int, optional
        Number of k nearest neihbors in kNN-classifier.
    clean_pred : float, optional
        Threshold with minimal % of cells within celltype to keep this cell type
        in the final annotation. Default: `None` (don't perform cleaning).
    do_preprocessing : bool, optional
        If counts are stored in layer to work with. Default: `True`.

    Returns
    ----------
    This function adds UMAP projection on the reference and annotated cell types to adata.
    """
    import rpy2.robjects as ro
    import anndata2ri
    
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(symphony))")
    
    anndata2ri.activate()
    
    if (layer is None):
        adata_tmp = sc.AnnData(X=adata.X, obs=adata.obs.copy(), var=pd.DataFrame(index=adata.var_names)).copy()
    else:
        adata_tmp = sc.AnnData(X=adata.layers[layer], obs=adata.obs.copy(), var=pd.DataFrame(index=adata.var_names)).copy()
    ro.globalenv["adata_query"] = adata_tmp
    
    ro.r('seurat.query <- as.Seurat(adata_query, counts = NULL, data = "X")')
    ro.r(f'reference <- readRDS("{reference_path}")')
    prepr = "FALSE"
    if do_preprocessing:
        prepr = "TRUE"
    if type(batch_key) == str:
        ro.r(f'query <- mapQuery(GetAssayData(object = seurat.query), seurat.query[[]], reference, do_normalize = {prepr}, vars = "{batch_key}")')
    elif type(batch_key) == list:
        var = '"' + '", "'.join(batch_key) + '"'
        ro.r(f'query <- mapQuery(GetAssayData(object = seurat.query), seurat.query[[]], reference, do_normalize = {prepr}, vars = c({var}))')
    else:
        ro.r(f'query <- mapQuery(GetAssayData(object = seurat.query), seurat.query[[]], reference, do_normalize = {prepr})')
    ro.r(f'query <- knnPredict(query, reference, reference$meta_data${reference_key}, k = {k})')
    predicted_ct = pd.Series(list(ro.r("as.character(query$meta_data$cell_type_pred_knn)")))
    if not (clean_pred is None):
        predicted_ct[~np.isin(predicted_ct, predicted_ct.value_counts()[predicted_ct.value_counts() / len(predicted_ct) * 100 > clean_pred].index)] = "Other"
    adata.obs[f"predicted_{reference_key}"] = list(predicted_ct)
    adata.obs[f"predicted_{reference_key}"] = adata.obs[f"predicted_{reference_key}"].astype("category")
    adata.obs[f"predicted_{reference_key}_prob"] = ro.r('query$meta_data$cell_type_pred_knn_prob')
    adata.obsm["X_umap_ref"] = ro.r('query$umap')
    
def create_symphony_reference(adata, pred_keys, save_rds_path, save_uwot_path, umap_key="X_umap", batch_key=None,
                              layer=None, n_hvgs=3000, n_pcs=30, do_preprocessing=True):
    """
    This is a wrapper function around Symphony reference building.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    pred_keys : Union[str, list]
        Name or list of names of cell type keys.
    save_rds_path : str
        Path where to save an .rds-object with Symphony reference.
    save_uwot_path : str
        Path where to save UMAP model.
    umap_key : str, optional
        Name of slot where to save Symphony UMAP. Default: `X_umap`
    layer : str, optional
        Name of layer with counts. Default: `None` (counts are expected in `adata.X`)
    batch_key : str, optional
        Name of column in `adata.obs` that corresponds to the name of samples / batches.
    n_pcs : int, optional
        Number of PCs that will be used for UMAP. Default: 30.
    n_hvgs : int, optional
        Number of HVGs to select. Default: 3000.
    do_preprocessing : bool, optional
        If counts are stored in layer to work with. Default: `True`.

    Returns
    ----------
    This function adds UMAP projection on the reference to adata.
    """
    import rpy2.robjects as ro
    import anndata2ri
    import pandas as pd
    
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(symphony))")
    
    anndata2ri.activate()
    
    if type(pred_keys) == list:
        keys = pred_keys + [batch_key]
    else:
        keys = [pred_keys] + [batch_key]
    if (layer is None):
        adata_tmp = sc.AnnData(X=adata.X, obs=pd.DataFrame(adata.obs[keys]), var=pd.DataFrame(index=adata.var_names)).copy()
    else:
        adata_tmp = sc.AnnData(X=adata.layers[layer], obs=pd.DataFrame(adata.obs[keys]), var=pd.DataFrame(index=adata.var_names)).copy()
        
    ro.globalenv["adata_reference"] = adata_tmp
    ro.r('seurat.reference <- as.Seurat(adata_reference, counts = NULL, data = "X")')
    
    add_batch = ""
    if type(batch_key) == str:
        add_batch = f', vars = "{batch_key}", vargenes_groups = "{batch_key}"'
    elif type(batch_key) == list:
        var = '"' + '", "'.join(batch_key) + '"'
        add_batch = f', vars = c({var}), vargenes_groups = "{batch_key[0]}"'
    prepr = "FALSE"
    if do_preprocessing:
        prepr = "TRUE"
    ro.r(f'reference <- buildReference(GetAssayData(object = seurat.reference), seurat.reference[[]], do_normalize = {prepr}, topn = {n_hvgs}, d = {n_pcs}, save_uwot_path = "{save_uwot_path}"{add_batch})')
    ro.r(f'saveRDS(reference, "{save_rds_path}")')
    adata.obsm[umap_key] = ro.r('reference$umap$embedding')

def set_groups_colors(adata, obs_names):
    fig, ax = plt.subplots()
    for obs in obs_names:
        if f"{obs}_colors" in adata.uns:
            del adata.uns[f"{obs}_colors"]
        sc.pl.umap(adata, color=obs, ax=ax, show=False)
    plt.close(fig)
    
def raw_to_save(adata, adata_counts=None, layer_counts="counts"):
    adata_new = adata.raw.to_adata()
    if not (adata_counts is None):
        if layer_counts != "X":
            adata_new.layers["counts"] = adata_counts[adata_new.obs_names, adata_new.var_names].layers[layer_counts].copy()
        else:
            adata_new.layers["counts"] = adata_counts[adata_new.obs_names, adata_new.var_names].X.copy()
    return adata_new

def reverse_dict(d):
    res = {}
    for i in d:
        for j in d[i]:
            res[j] = i
    return res

class Annotation:
    def __init__(self):
        self.order = []
        self.annotations = {}
    def AddLevel(self, annotation, cluster_key, annot_type="ct_to_leiden"):
        self.order.append(cluster_key)
        if annot_type == "leiden_to_ct":
            self.annotations[cluster_key] = annotation
        elif annot_type == "ct_to_leiden":
            self.annotations[cluster_key] = reverse_dict(annotation)
    def Apply(self, adata, celltype_key):
        celltypes = ["Unknown"] * len(adata)
        for cluster_key in self.order:
            counter = 0
            for bc in adata.obs_names:
                cluster = adata.obs.loc[bc][cluster_key]
                if cluster in self.annotations[cluster_key]:
                    celltypes[counter] = self.annotations[cluster_key][cluster]
                counter += 1
        adata.obs[celltype_key] = celltypes
        adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
        self.groups = list(adata.obs[celltype_key].cat.categories)
        self.groups.remove("Unknown")

def compare_groups(adata, cluster_key, first_group, second_group, **kwargs):
    """
    This function performs DEG analysis between two groups.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str
        Name of column in `adata.obs` that corresponds to the name of clusterind.
    first_group : list
        List of the first group cluster names.
    second_group : list
        List of the first group cluster names.
    kwargs : keyword arguments
        Additional arguments for `sc.tl.rank_genes_groups()`.

    Returns
    ----------
    None.
    """
    adata_tmp = adata[np.isin(adata.obs[cluster_key], first_group + second_group)]
    first_name = "(" + ", ".join(first_group) + ")"
    second_name = "(" + ", ".join(second_group) + ")"
    adata_tmp.obs["comp"] = [first_name if i in first_group else second_name for i in adata_tmp.obs[cluster_key]]
    adata_tmp.obs["comp"] = adata_tmp.obs["comp"].astype("category")
    sc.tl.rank_genes_groups(adata_tmp, groupby="comp", groups=[first_name], reference=second_name, **kwargs)
    return sc.get.rank_genes_groups_df(adata_tmp, group=first_name)

def barplot(adata, cluster_key, batch_key):
    """
    This function creates a barplot with the distribution of cells
    from different batches across clusters.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str
        Name of column in `adata.obs` that corresponds to the name of clusterind.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples / batches.

    Returns
    ----------
    None.
    """
    import matplotlib.pyplot as plt
    import seaborn
    
    fig_width = len(adata.obs[cluster_key].cat.categories) * 0.3
    fig, ax = plt.subplots(dpi=150, figsize=(fig_width, 2))
    
    colors = dict(zip(adata.obs[batch_key].cat.categories,
                      adata.uns[f"{batch_key}_colors"]))
    
    sizes = adata.obs.groupby([batch_key, cluster_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=cluster_key, index=batch_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    props.plot.bar(stacked=True, width=1,
                   edgecolor="black", ax=ax, color=colors)
    plt.xticks(rotation=90)
    ax.set_xlabel("")
    ax.legend(loc=(1.01, 0.45), edgecolor="white")

def draw_panel(adata, genes, title="", heigh=4, dpi=150):
    import scvelo as scv
    
    if len(genes) == 1:
        fig, ax = plt.subplots(figsize=(heigh, heigh), dpi=dpi)
        scv.pl.umap(adata, color=genes[0], ax=ax, show=False, cmap=Reds)
        plt.suptitle(title)
    else:
        fig, axes = plt.subplots(ncols=len(genes), figsize=(heigh * len(genes) * 0.8, heigh), dpi=dpi)
        counter = 0
        for gene in genes:
            scv.pl.umap(adata, color=gene, ax=axes[counter], show=False, cmap=Reds)
            counter += 1
        fig.suptitle(title)
        fig.tight_layout()

def draw_counts_qc(adata, title, mt_max=None, umi_min=None, min_genes=None):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    axes[0].hist(
        adata.obs.total_counts,
        bins=10**np.histogram(np.log10(adata.obs.total_counts), bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    ylim = axes[0].get_ylim()[1]
    if not (umi_min is None):
        axes[0].plot([umi_min, umi_min], [0, ylim], linewidth=0.8)
    axes[0].set_ylim(0, ylim)
    axes[0].set_xscale("log"); axes[0].grid()
    axes[0].set_title("Number of UMIs per cell")
    axes[0].set_ylabel("Number of cells")

    axes[1].hist(
        adata.obs.n_genes_by_counts,
        bins=10**np.histogram(np.log10(adata.obs.n_genes_by_counts), bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    ylim = axes[1].get_ylim()[1]
    if not min_genes is None:
        axes[1].plot([min_genes, min_genes], [0, ylim], linewidth=0.8)
    axes[1].set_ylim(0, ylim)
    axes[1].set_xscale("log"); axes[1].grid()
    axes[1].set_title("Number of genes per cell")

    axes[2].hist(
        adata.obs.pct_counts_mt,
        bins=np.histogram(adata.obs.pct_counts_mt, bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    ylim = axes[2].get_ylim()[1]
    if not mt_max is None:
        axes[2].plot([mt_max, mt_max], [0, ylim], linewidth=0.8)
    axes[2].set_ylim(0, ylim)
    axes[2].grid()
    axes[2].set_title("Percent of mitochondrial UMIs")

    fig.suptitle(title, fontsize=15, y=0.93)
    fig.tight_layout()
    
    return fig
    
def draw_counts_umap(adata, mt_max):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 3.5))
    
    sc.pl.umap(adata, color="total_counts", frameon=False, cmap="coolwarm",
               title="Number of UMIs per cell", show=False, ax=axes[0], vmax=30000)

    sc.pl.umap(adata, color="n_genes_by_counts", frameon=False, cmap="coolwarm",
               title="Number of genes per cell", show=False, ax=axes[1], vmax=10000)

    sc.pl.umap(adata, color="pct_counts_mt", frameon=False, cmap="coolwarm",
               title="Percent of mitochondrial counts", show=False, ax=axes[2], vmax=mt_max)
    fig.tight_layout()
    
    return fig
    
def draw_scrublet_statistics(adata, threshold="default"):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    axes[0].hist(
        adata.obs.doublet_score,
        bins=np.histogram(adata.obs.doublet_score, bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    ylim = axes[0].get_ylim()[1]
    if "threshold" in adata.uns["scrublet"] and threshold == "default":
        threshold = adata.uns["scrublet"]["threshold"]
    axes[0].plot([threshold, threshold], [0, ylim], linewidth=0.8)
    axes[0].set_ylim(0, ylim)
    axes[0].grid()
    axes[0].set_title("Scrublet score distribution")

    sc.pl.umap(adata, color="doublet_score", show=False, ax=axes[1], cmap="coolwarm",
                title=f"Doublet score", frameon=False, vmax=threshold)
    
    adata_tmp = adata.copy()
    adata_tmp.obs["predicted_doublet"] = ["True" if ds > threshold else "False" for ds in adata_tmp.obs.doublet_score]
    sc.pl.umap(adata_tmp, color="predicted_doublet", show=False, ax=axes[2],
               title=f"Predicted doublets ({sum(adata_tmp.obs['predicted_doublet'] == 'True')})",
               groups=["True"], legend_loc=None, frameon=False)

    fig.suptitle(f"Scrublet statistics", fontsize=15, y=0.93)
    fig.tight_layout()
    
    return fig
    
def draw_clear_umap(adata, mt_max):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    sc.pl.umap(adata, color="leiden", frameon=False, legend_loc="on data",
               legend_fontoutline=2, legend_fontsize=12,
               title="Clusters (Leiden)", show=False, ax=axes[0])

    sc.pl.umap(adata, color="doublet_score", frameon=False, cmap="coolwarm",
               title="Scrublet score", show=False, ax=axes[1])

    sc.pl.umap(adata, color="pct_counts_mt", frameon=False, cmap="coolwarm",
               title="Percent of mitochondrial counts", show=False, ax=axes[2], vmax=mt_max)

    fig.tight_layout()
    
    return fig


def draw_cloneID_QC(adata_FP, FP):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    for i in range(3):
        cloneIDs_per_cell = (adata_FP.X > (i + 1)).sum(axis=1).A.T[0]
        sns.histplot(cloneIDs_per_cell, ax=axes[i], bins=max(cloneIDs_per_cell), edgecolor="black", linewidth=1, alpha=1, discrete=True)
        axes[i].set_xlim(-1, 10.5)
        axes[i].set_ylabel("")
        axes[i].set_xlabel("Number of cloneIDs per cell")
        axes[i].grid(alpha=0.3)
        axes[i].set_title(f"min UMI of cloneIDs >= {i + 2}")
        axes[i].set_xticks([0, 2, 4, 6, 8, 10])
        axes[i].set_yscale("log")

    axes[0].set_ylabel("Number of cells")
    fig.suptitle(f"Distribution of number of cloneIDs per cell ({FP})", fontsize=15, y=0.93)
    fig.tight_layout()
    return fig
    

def concat_figures(savefig):
    import os
    from PyPDF2 import PdfFileReader, PdfFileWriter, PageObject
    
    os.system(f"mv {savefig} {savefig}_tmp")
    file = open(savefig + "_tmp", "rb")
    pdf = PdfFileReader(file)

    pages = []; width = 0; height = 0
    for pageNum in range(pdf.numPages - 1, -1, -1):
        pageObj = pdf.getPage(pageNum)
        pages.append(pageObj)
        height += pageObj.mediaBox.getHeight()
        if pages[0].mediaBox.getWidth() > width:
            width = pageObj.mediaBox.getWidth()

    merged_page = PageObject.createBlankPage(None, width, height)

    x = 0
    for page in pages:
        merged_page.mergeScaledTranslatedPage(page, 1, 0, x)
        x = float(x) + float(page.mediaBox.getHeight())

    writer = PdfFileWriter()
    writer.addPage(merged_page)

    with open(savefig, "wb") as f:
        writer.write(f)
    os.system(f"rm {savefig}_tmp")
    
colormap_boxplot_1 = {
    "Single" : "#53A346",
    "Multiple in few" : "#5259DB",
    "Multiple in single" : "#E7352C",
    "Unknown" : "#999999",
}
    
def draw_multiple_cloneIDs(adata, FP):
    fig, axes = plt.subplots(figsize=(12, 4), ncols=4)
    colormap = {
        "Single" : "#53A346",
        "Multiple in few" : "#5259DB",
        "Multiple in single" : "#E7352C",
        "Unknown" : "#999999",
    }

    order = ["Unknown", "Single", "Multiple in few", "Multiple in single"]
    sns.boxplot(x=f"cloneID_multiplet_{FP}", y="doublet_score", data=adata.obs, ax=axes[0], order=order,
                showfliers=False, palette=colormap, **black_borders_boxplot, linewidth=1)
    axes[0].set_ylabel("Scrublet score")
    axes[0].set_xlabel("")
    axes[0].grid(alpha=0.3)
    axes[0].tick_params(axis="x", labelrotation=30, labelsize=12)
    
    adata.obs[f"cloneID_multiplet_{FP}"] = adata.obs[f"cloneID_multiplet_{FP}"].astype("category").cat.reorder_categories(order)
    adata.uns[f"cloneID_multiplet_{FP}_colors"] = [colormap[i] for i in adata.obs[f"cloneID_multiplet_{FP}"].cat.categories]
    
    N_single = sum(adata.obs[f"cloneID_multiplet_{FP}"] == "Single")
    N_multiple_few = sum(adata.obs[f"cloneID_multiplet_{FP}"] == "Multiple in few")
    N_miltiple_one = sum(adata.obs[f"cloneID_multiplet_{FP}"] == "Multiple in single")
    N_whole = len(adata)
    
    sc.pl.umap(adata, color=f"cloneID_multiplet_{FP}", groups=["Single"],
               legend_loc=None, ax=axes[1], show=False, frameon=False,
               title=f"Single {FP}\n({N_single} out of {N_whole} cells)")
    sc.pl.umap(adata, color=f"cloneID_multiplet_{FP}", groups=["Multiple in few"],
               legend_loc=None, ax=axes[2], show=False, frameon=False,
               title=f"Multiple {FP} in few\n({N_multiple_few} out of {N_whole} cells)")
    sc.pl.umap(adata, color=f"cloneID_multiplet_{FP}", groups=["Multiple in single"],
               legend_loc=None, ax=axes[3], show=False, frameon=False,
               title=f"Multiple {FP} in single\n({N_miltiple_one} out of {N_whole} cells)")
    
    fig.suptitle(FP, fontsize=15, y=0.93, x=0.55)
    fig.tight_layout()
    return fig
    
def draw_clear_counts_qc(adata):
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    axes[0].hist(
        adata.obs.total_counts,
        bins=10**np.histogram(np.log10(adata.obs.total_counts), bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    axes[0].set_xscale("log"); axes[0].grid()
    axes[0].set_title("Number of UMIs per cell")
    axes[0].set_ylabel("Number of cells")

    axes[1].hist(
        adata.obs.n_genes_by_counts,
        bins=10**np.histogram(np.log10(adata.obs.n_genes_by_counts), bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    axes[1].set_xscale("log"); axes[1].grid()
    axes[1].set_title("Number of genes per cell")

    axes[2].hist(
        adata.obs.pct_counts_mt,
        bins=np.histogram(adata.obs.pct_counts_mt, bins=30)[1],
        edgecolor="black",
        linewidth=1
    )
    axes[2].grid()
    axes[2].set_title("Percent of mitochondrial UMIs")

    fig.suptitle(f"Final QC metrics (total number of cells = {len(adata)})", fontsize=15, y=0.93)
    fig.tight_layout()
    
    return fig

def plot_confusion(cell_difference, mean_logFC, final_assignment, FP1, FP2, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    sns.scatterplot(x=cell_difference, y=mean_logFC, ax=ax, linewidth=0.5, edgecolor="black",
                    s=15, hue=final_assignment, palette={FP1: "green", FP2: "red"})
    ax.set_ylim(*ax.get_ylim())
    ax.set_xlim(-30, 30)
    #ax.plot([0, 0], [-10, 10], color="red")
    #ax.plot([-30, 30], [0, 0], color="red")
    ax.grid(alpha=0.3)
    ax.set_xlabel(f"N({FP1} > {FP2}) — N({FP2} > {FP1})")
    ax.set_ylabel(f"mean logFC({FP1} / {FP2})")
    ax.set_title("Confusion of viral BCs")
    
def plot_confusion_panel(adatas_FP, cell_difference, mean_logFC, final_assignment,
                         FP1, FP2, clear_FP1, clear_FP2, title):
    from matplotlib_venn import venn2
    
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))
    
    intersect = np.isin(adatas_FP[FP1].var_names, adatas_FP[FP2].var_names).sum()
    venn2(
        subsets=(
            len(adatas_FP[FP1].var_names) - intersect,
            len(adatas_FP[FP2].var_names) - intersect,
            intersect,
        ),
        set_labels=(FP1, FP2),
        set_colors=("green", "red"),
        ax=axes[0],
    )
    axes[0].set_title(f"Intersection between\n{FP1} and {FP2}")
    
    plot_confusion(cell_difference, mean_logFC, final_assignment, FP1, FP2, ax=axes[1])
    
    fig.tight_layout()
    
    intersect = np.isin(clear_FP1.var_names, clear_FP2.var_names).sum()
    venn2(
        subsets=(
            len(clear_FP1.var_names) - intersect,
            len(clear_FP2.var_names) - intersect,
            intersect,
        ),
        set_labels=(FP1, FP2),
        set_colors=("green", "red"),
        ax=axes[2],
    )
    axes[2].set_title(f"Intersection between\n{FP1} and {FP2} after correction")
    
    fig.suptitle(title, fontsize=15, y=1.05)
    
    return fig


def predict_doublets(mdata, FP, min_cloneID_umi=2):
    if isinstance(FP, str):
        df = pd.DataFrame(
            mdata[FP].X.A >= min_cloneID_umi,
            index=mdata[FP].obs.index,
            columns=mdata[FP].var.index,
        )
        df_zero = df[df.sum(axis=1) == 0]
        df_double = df[df.sum(axis=1) >= 2]
        df_double_vc = df_double.value_counts()

        cloneID_multiplet = []
        for bc in mdata[FP].obs.index:
            if bc in df_zero.index:
                cloneID_multiplet.append("Unknown")
            elif bc not in df_double.index:
                cloneID_multiplet.append("Single")
            elif df_double_vc[tuple(df_double.loc[bc])] >= 2:
                cloneID_multiplet.append("Multiple in few")
            else:
                cloneID_multiplet.append("Multiple in single")
        mdata[FP].obs[f"clone_id_multiplet"] = cloneID_multiplet
    elif isinstance(FP, list):
        cells = mdata.obs_names
        FP1 = FP[0]; FP2 = FP[1]
        
        df_FP1 = pd.DataFrame(
            mdata[FP1].X.A >= min_cloneID_umi,
            index=mdata[FP1].obs.index,
            columns=mdata[FP1].var.index
        )
        df_FP2 = pd.DataFrame(
            mdata[FP2].X.A >= min_cloneID_umi,
            index=mdata[FP2].obs.index,
            columns=mdata[FP2].var.index
        )
        df_full = pd.concat([df_FP1, df_FP2], axis=1)

        cells_zero = cells[(df_FP1.sum(axis=1) == 0) & (df_FP2.sum(axis=1) == 0)]
        cells_FP1 = cells[(df_FP1.sum(axis=1) >= 1) & (df_FP2.sum(axis=1) == 0)]
        cells_FP2 = cells[(df_FP1.sum(axis=1) == 0) & (df_FP2.sum(axis=1) >= 1)]
        df_double = df_full[(df_FP1.sum(axis=1) >= 1) & (df_FP2.sum(axis=1) >= 1)]
        df_double_vc = df_double.value_counts()

        cloneID_multiplet = []
        for bc in cells:
            if bc in cells_zero:
                cloneID_multiplet.append("Unknown")
            elif bc in cells_FP1:
                cloneID_multiplet.append(FP1)
            elif bc in cells_FP2:
                cloneID_multiplet.append(FP2)
            elif df_double_vc[tuple(df_double.loc[bc])] >= 2:
                cloneID_multiplet.append("Multiple in few")
            else:
                cloneID_multiplet.append("Multiple in single")
        mdata.obs[f"clone_id_multiplet"] = cloneID_multiplet


def clones_QC(adata, clones_obs, title=""):
    fig, axes = plt.subplots(figsize=(12, 4), ncols=3)

    def draw_clone(adata, clones_obs, clone, ax):
        adata_tmp = adata.copy()
        adata_tmp.obs["group"] = ["0" if i != clone else "1" for i in clones_obs]
        sc.pl.umap(adata_tmp, ax=ax, show=False, frameon=False, color="group", groups=["1"],
                   legend_loc=None, title=clone, s=30)

    sns.histplot(
        pd.Series(clones_obs).value_counts()[1:][pd.Series(clones_obs).value_counts()[1:] > 1],
        discrete=True,
        alpha=1,
        edgecolor="black",
        ax=axes[0],
    )
    axes[0].set_yscale("log")
    axes[0].grid(alpha=0.3)
    axes[0].set_ylabel("")
    axes[0].set_title("Clone size distribution")

    draw_clone(adata, clones_obs, pd.Series(clones_obs).value_counts().index[1], ax=axes[1])
    draw_clone(adata, clones_obs, pd.Series(clones_obs).value_counts().index[2], ax=axes[2])
    
    fig.suptitle(title, fontsize=15, y=0.93)
    fig.tight_layout()
    
    return fig
    
    
colormap_boxplot_2 = {
    "GFPbc" : "#53A346",
    "TOMbc" : "#E7352C",
    "Unknown" : "#999999",
    "Multiple in few" : "#5259DB",
    "Multiple in single" : "#815747",
}
    
    
def draw_double_injection_stats_1(adata, FP1, FP2):
    fig, axes = plt.subplots(figsize=(12, 5), ncols=2, dpi=200)

    colormap = {
        FP1 : "#53A346",
        FP2 : "#E7352C",
        "Unknown" : "#999999",
        "Multiple in few" : "#5259DB",
        "Multiple in single" : "#815747",
    }

    order = ["Unknown", FP1, FP2, "Multiple in few", "Multiple in single"]
    sns.boxplot(x=f"cloneID_multiplet", y="doublet_score", data=adata.obs, ax=axes[0], order=order,
                showfliers=False, palette=colormap, **black_borders_boxplot, linewidth=1)
    axes[0].tick_params(axis="x", labelrotation=30, labelsize=12)
    axes[0].set_ylabel("Scrublet score")
    axes[0].set_xlabel("")
    axes[0].grid(alpha=0.3)
    
    axes[1].pie(
        adata.obs.cloneID_multiplet.value_counts().values,
        labels=(
            np.array(adata.obs.cloneID_multiplet.value_counts().index) +
            " (" + np.round(adata.obs.cloneID_multiplet.value_counts().values / len(adata) * 100, 2).astype(str) + "%)"
        ),
        wedgeprops={"edgecolor": "black"},
        textprops={"size": 10},
        colors=[colormap[i] for i in adata.obs.cloneID_multiplet.value_counts().index],
    )
    axes[1].set_title("Composition of the sample")
    
    fig.tight_layout()
    return fig

def draw_double_injection_stats_2(adata, FP1, FP2):
    fig, axes = plt.subplots(figsize=(12, 3), ncols=4)

    colormap = {
        FP1 : "#53A346",
        FP2 : "#E7352C",
        "Unknown" : "#999999",
        "Multiple in few" : "#5259DB",
        "Multiple in single" : "#815747",
    }

    sc.pl.umap(adata, color="cloneID_multiplet", ax=axes[0], show=False, groups=[FP1], frameon=False, legend_loc=None,
               title=f"{FP1}+ {FP2}- cells\n({sum(adata.obs.cloneID_multiplet == FP1)} out of {len(adata)} cells)")
    sc.pl.umap(adata, color="cloneID_multiplet", ax=axes[1], show=False, groups=[FP2], frameon=False, legend_loc=None,
               title=f"{FP1}- {FP2}+ cells\n({sum(adata.obs.cloneID_multiplet == FP2)} out of {len(adata)} cells)")
    sc.pl.umap(adata, color="cloneID_multiplet", ax=axes[2], show=False, groups=["Multiple in few"], frameon=False, legend_loc=None,
               title=f"Multiple in few\n({sum(adata.obs.cloneID_multiplet == 'Multiple in few')} out of {len(adata)} cells)")
    sc.pl.umap(adata, color="cloneID_multiplet", ax=axes[3], show=False, groups=["Multiple in single"], frameon=False, legend_loc=None,
               title=f"Multiple in single\n({sum(adata.obs.cloneID_multiplet == 'Multiple in single')} out of {len(adata)} cells)")
    
    fig.tight_layout()
    return fig
    
    
def draw_expression_distribution(adata, FP1, FP2):
    fig, axes = plt.subplots(ncols=4, figsize=(12, 3))
    
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == FP1, FP1].X.A.T[0], color="green", fill=True,
                ax=axes[0], label=FP1)
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == FP1, FP2].X.A.T[0], color="red", fill=True,
                ax=axes[0], label=FP2)
    axes[0].set_title(f"{FP1}+{FP2}- cells")
    axes[0].set_ylabel("")
    axes[0].grid(alpha=0.3)
    axes[0].set_xlabel("Expression level")
    axes[0].legend()
    
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == FP2, FP1].X.A.T[0], color="green", fill=True,
                ax=axes[1], label=FP1)
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == FP2, FP2].X.A.T[0], color="red", fill=True,
                ax=axes[1], label=FP2)
    axes[1].set_title(f"{FP1}-{FP2}+ cells")
    axes[1].set_ylabel("")
    axes[1].grid(alpha=0.3)
    axes[1].set_xlabel("Expression level")
    axes[1].legend()
    
    sns.kdeplot(adata[np.isin(adata.obs["cloneID_multiplet"], ["Multiple in few", "Multiple in single"]), FP1].X.A.T[0],
                color="green", fill=True, ax=axes[2], label=FP1)
    sns.kdeplot(adata[np.isin(adata.obs["cloneID_multiplet"], ["Multiple in few", "Multiple in single"]), FP2].X.A.T[0],
                color="red", fill=True, ax=axes[2], label=FP2)
    axes[2].set_title(f"{FP1}+{FP2}+ cells")
    axes[2].set_ylabel("")
    axes[2].grid(alpha=0.3)
    axes[2].set_xlabel("Expression level")
    axes[2].legend()
    
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == "Unknown", FP1].X.A.T[0], color="green", fill=True,
                ax=axes[3], label=FP1)
    sns.kdeplot(adata[adata.obs["cloneID_multiplet"] == "Unknown", FP2].X.A.T[0], color="red", fill=True,
                ax=axes[3], label=FP2)
    axes[3].set_title(f"{FP1}-{FP2}- cells")
    axes[3].set_ylabel("")
    axes[3].grid(alpha=0.3)
    axes[3].set_xlabel("Expression level")
    axes[3].legend()
    
    fig.tight_layout()
    return fig


def resolve_confusion(adatas_FP):
    FP1 = list(adatas_FP.keys())[0]
    FP2 = list(adatas_FP.keys())[1]

    ambiguous_bc = adatas_FP[FP1].var_names[np.isin(adatas_FP[FP1].var_names, adatas_FP[FP2].var_names)]

    difference = (
        np.log1p(adatas_FP[FP1][:, ambiguous_bc].X) -
        np.log1p(adatas_FP[FP2][:, ambiguous_bc].X)
    )
    mean_logFC = np.true_divide(
        difference.sum(axis=0).A[0],
        (difference != 0).sum(axis=0).A[0]
    )
    cell_difference = (
        ((adatas_FP[FP1][:, ambiguous_bc].X - adatas_FP[FP2][:, ambiguous_bc].X) > 0).sum(axis=0).A[0] -
        ((adatas_FP[FP2][:, ambiguous_bc].X - adatas_FP[FP1][:, ambiguous_bc].X) > 0).sum(axis=0).A[0]
    )
    final_assignment = []
    for i in range(len(ambiguous_bc)):
        if cell_difference[i] > 0:
            final_assignment.append(FP1)
        elif cell_difference[i] < 0:
            final_assignment.append(FP2)
        elif mean_logFC[i] > 0:
            final_assignment.append(FP1)
        else:
            final_assignment.append(FP2)
    final_assignment = np.array(final_assignment) 
    FP1_bcs = ambiguous_bc[final_assignment == FP1]
    FP2_bcs = ambiguous_bc[final_assignment == FP2]

    clear_FP1 = adatas_FP[FP1][:, ~np.isin(adatas_FP[FP1].var_names, ambiguous_bc)]
    clear_FP2 = adatas_FP[FP2][:, ~np.isin(adatas_FP[FP2].var_names, ambiguous_bc)]
    both_FP = sc.AnnData(
        X=adatas_FP[FP1][:, ambiguous_bc].X + adatas_FP[FP2][:, ambiguous_bc].X,
        obs=adatas_FP[FP1].obs,
        var=adatas_FP[FP1][:, ambiguous_bc].var,
    )
    clear_FP1 = sc.AnnData(
        X=hstack([clear_FP1.X, both_FP[:, FP1_bcs].X]),
        obs=clear_FP1.obs,
        var=pd.DataFrame(index=list(clear_FP1.var_names) + list(FP1_bcs)),
    )
    clear_FP2 = sc.AnnData(
        X=hstack([clear_FP2.X, both_FP[:, FP2_bcs].X]),
        obs=clear_FP2.obs,
        var=pd.DataFrame(index=list(clear_FP2.var_names) + list(FP2_bcs)),
    )

    return {FP1: clear_FP1, FP2: clear_FP2}


def fast_clones(adata, min_cloneID_umi):
    df = pd.DataFrame(
        adata.X.A >= min_cloneID_umi,
        index=adata.obs.index,
        columns=adata.var.index
    )
    clones = {}
    clones_obs = []
    clone_counter = 0
    for bc in adata.obs_names:
        cloneIDs = frozenset(df.columns[df.loc[bc]])
        if len(cloneIDs) == 0:
            clones_obs.append("Unknown")
        elif cloneIDs in clones:
            clones_obs.append(clones[cloneIDs])
        else:
            clones[cloneIDs] = f"clone_{clone_counter}"
            clones_obs.append(clones[cloneIDs])
            clone_counter += 1
    return pd.Series(clones_obs, index=adata.obs_names)


def gex_qc_report(adata_raw, adatas_FP, mt_name="mt-", mt_max=25, savefig="output.pdf",
                  show=False, umi_min=None, min_genes=None, scrublet_threshold=None,
                  min_cloneID_umi=2):
    """
    This function creates a QC report.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with raw counts.
    adatas_FP : dict[str | AnnData]
        Dictionary with Viral IDs UMI matrices.
    clones : DataFrame
        DataFrame with information about clones.
    mt_name : str, optional
        Prefix of mitochondrial genes.
    mt_max : float, optional
        Mitochondrial UMIs content threshold.
    savefig : str, optional
        Where the final figure should be saved.
    show : bool, optional
        If it's needed to show plt figures.
        
    Returns
    ----------
    Annotated data matrix after quality control.
    """
    import matplotlib.backends.backend_pdf
    with sns.axes_style("ticks"):

        adata = adata_raw.copy()
        adata.var_names_make_unique()
        adata.layers["counts"] = adata.X.copy()
        adata.var["mt"] = adata.var_names.str.startswith(mt_name)
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], log1p=False, inplace=True, percent_top=None)

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=5000)
        adata.layers["scaled"] = sc.pp.scale(adata.X, max_value=10, copy=True)
        adata.obsm["X_pca"] = sc.pp.pca(adata[:, adata.var.highly_variable].layers["scaled"])
        sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
        sc.tl.umap(adata)

        figs = []

        figs.append(draw_counts_qc(
            adata, mt_max=mt_max, umi_min=umi_min, min_genes=min_genes,
            title=f"Sample {description['sample_id']}\n\nTotal number of cells before QC = {len(adata)}"
        ))
        figs.append(draw_counts_umap(adata, mt_max))
        
        adata_clear = adata.copy()
        adata_clear.X = adata_clear.layers["counts"].copy()
        adata_clear = adata_clear[adata_clear.obs.pct_counts_mt <= mt_max]
        if not (umi_min is None):
            adata_clear = adata_clear[adata_clear.obs.total_counts >= umi_min]
        if not (min_genes is None):
            adata_clear = adata_clear[adata_clear.obs.n_genes_by_counts >= min_genes]
        
        sc.pp.normalize_total(adata_clear, target_sum=1e4)
        sc.pp.log1p(adata_clear)
        sc.pp.highly_variable_genes(adata_clear, n_top_genes=5000)
        adata_clear.layers["scaled"] = sc.pp.scale(adata_clear.X, max_value=10, copy=True)
        adata_clear.obsm["X_pca"] = sc.pp.pca(adata_clear[:, adata_clear.var.highly_variable].layers["scaled"])
        sc.pp.neighbors(adata_clear, n_pcs=30, n_neighbors=20)
        sc.tl.umap(adata_clear)
        
        from scipy.sparse import vstack, hstack, csr_matrix

        figs.append(draw_clear_counts_qc(adata_clear))
        if len(adatas_FP) == 2:
            FP1 = list(adatas_FP.keys())[0]
            FP2 = list(adatas_FP.keys())[1]

            ambiguous_bc = adatas_FP[FP1].var_names[np.isin(adatas_FP[FP1].var_names, adatas_FP[FP2].var_names)]

            difference = (
                np.log1p(adatas_FP[FP1][:, ambiguous_bc].X) -
                np.log1p(adatas_FP[FP2][:, ambiguous_bc].X)
            )
            mean_logFC = np.true_divide(
                difference.sum(axis=0).A[0],
                (difference != 0).sum(axis=0).A[0]
            )
            cell_difference = (
                ((adatas_FP[FP1][:, ambiguous_bc].X - adatas_FP[FP2][:, ambiguous_bc].X) > 0).sum(axis=0).A[0] -
                ((adatas_FP[FP2][:, ambiguous_bc].X - adatas_FP[FP1][:, ambiguous_bc].X) > 0).sum(axis=0).A[0]
            )
            final_assignment = []
            for i in range(len(ambiguous_bc)):
                if cell_difference[i] > 0:
                    final_assignment.append(FP1)
                elif cell_difference[i] < 0:
                    final_assignment.append(FP2)
                elif mean_logFC[i] > 0:
                    final_assignment.append(FP1)
                else:
                    final_assignment.append(FP2)
            final_assignment = np.array(final_assignment) 
            FP1_bcs = ambiguous_bc[final_assignment == FP1]
            FP2_bcs = ambiguous_bc[final_assignment == FP2]

            clear_FP1 = adatas_FP[FP1][:, ~np.isin(adatas_FP[FP1].var_names, ambiguous_bc)]
            clear_FP2 = adatas_FP[FP2][:, ~np.isin(adatas_FP[FP2].var_names, ambiguous_bc)]
            both_FP = sc.AnnData(
                X=adatas_FP[FP1][:, ambiguous_bc].X + adatas_FP[FP2][:, ambiguous_bc].X,
                obs=adatas_FP[FP1].obs,
                var=adatas_FP[FP1][:, ambiguous_bc].var,
            )
            clear_FP1 = sc.AnnData(
                X=hstack([clear_FP1.X, both_FP[:, FP1_bcs].X]),
                obs=clear_FP1.obs,
                var=pd.DataFrame(index=list(clear_FP1.var_names) + list(FP1_bcs)),
            )
            clear_FP2 = sc.AnnData(
                X=hstack([clear_FP2.X, both_FP[:, FP2_bcs].X]),
                obs=clear_FP2.obs,
                var=pd.DataFrame(index=list(clear_FP2.var_names) + list(FP2_bcs)),
            )

            figs.append(
                plot_confusion_panel(adatas_FP, cell_difference, mean_logFC, final_assignment,
                                     FP1, FP2, clear_FP1, clear_FP2, title="Confusion analysis")
            )
            
        if len(adatas_FP) == 2:
            adatas_FP = {FP1: clear_FP1, FP2: clear_FP2}
        for FP in description["FPs"]:
            figs.append(draw_cloneID_QC(adatas_FP[FP], FP))
            
        if scrublet_threshold is None:            
            sce.pp.scrublet(adata_clear, verbose=False)
        else:
            sce.pp.scrublet(adata_clear, verbose=False, threshold=scrublet_threshold)
        adata_clear.obs["predicted_doublet"] = adata_clear.obs["predicted_doublet"].astype(str)
        figs.append(draw_scrublet_statistics(adata_clear))

        for FP in list(adatas_FP.keys()):
            df = pd.DataFrame(
                adatas_FP[FP].X.A >= min_cloneID_umi,
                index=adatas_FP[FP].obs.index,
                columns=adatas_FP[FP].var.index
            )
            df_zero = df[df.sum(axis=1) == 0]
            df_double = df[df.sum(axis=1) >= 2]
            df_double_vc = df_double.value_counts()

            cloneID_multiplet = []
            for bc in adata_clear.obs.index:
                if bc in df_zero.index:
                    cloneID_multiplet.append("Unknown")
                elif bc not in df_double.index:
                    cloneID_multiplet.append("Single")
                elif df_double_vc[tuple(df_double.loc[bc])] >= 2:
                    cloneID_multiplet.append("Multiple in few")
                else:
                    cloneID_multiplet.append("Multiple in single")

            adata_clear.obs[f"cloneID_multiplet_{FP}"] = cloneID_multiplet
            figs.append(draw_multiple_cloneIDs(adata_clear, FP=FP))
            
        if len(adatas_FP) == 2:
            cells = adatas_FP[FP1].obs_names

            df_FP1 = pd.DataFrame(
                adatas_FP[FP1].X.A >= min_cloneID_umi,
                index=adatas_FP[FP1].obs.index,
                columns=adatas_FP[FP1].var.index
            )
            df_FP2 = pd.DataFrame(
                adatas_FP[FP2].X.A >= min_cloneID_umi,
                index=adatas_FP[FP2].obs.index,
                columns=adatas_FP[FP2].var.index
            )
            df_full = pd.concat([df_FP1, df_FP2], axis=1)

            cells_zero = cells[(df_FP1.sum(axis=1) == 0) & (df_FP2.sum(axis=1) == 0)]
            cells_FP1 = cells[(df_FP1.sum(axis=1) >= 1) & (df_FP2.sum(axis=1) == 0)]
            cells_FP2 = cells[(df_FP1.sum(axis=1) == 0) & (df_FP2.sum(axis=1) >= 1)]
            df_double = df_full[(df_FP1.sum(axis=1) >= 1) & (df_FP2.sum(axis=1) >= 1)]
            df_double_vc = df_double.value_counts()

            cloneID_multiplet = []
            for bc in adata_clear.obs.index:
                if bc in cells_zero:
                    cloneID_multiplet.append("Unknown")
                elif bc in cells_FP1:
                    cloneID_multiplet.append(FP1)
                elif bc in cells_FP2:
                    cloneID_multiplet.append(FP2)
                elif df_double_vc[tuple(df_double.loc[bc])] >= 2:
                    cloneID_multiplet.append("Multiple in few")
                else:
                    cloneID_multiplet.append("Multiple in single")
            adata_clear.obs[f"cloneID_multiplet"] = cloneID_multiplet

            figs.append(draw_double_injection_stats_1(adata_clear, FP1, FP2))
            figs.append(draw_double_injection_stats_2(adata_clear, FP1, FP2))
            figs.append(draw_expression_distribution(adata_clear, FP1, FP2))
            
        for FP in list(adatas_FP.keys()):
            df = pd.DataFrame(
                adatas_FP[FP][adata_clear.obs_names].X.A >= min_cloneID_umi,
                index=adatas_FP[FP][adata_clear.obs_names].obs.index,
                columns=adatas_FP[FP][adata_clear.obs_names].var.index
            )
            clones = {}
            clones_obs = []
            clone_counter = 0
            for bc in adata_clear.obs_names:
                cloneIDs = frozenset(df.columns[df.loc[bc]])
                if len(cloneIDs) == 0:
                    clones_obs.append("Unknown")
                elif cloneIDs in clones:
                    clones_obs.append(clones[cloneIDs])
                else:
                    clones[cloneIDs] = f"clone_{clone_counter}"
                    clones_obs.append(clones[cloneIDs])
                    clone_counter += 1
            figs.append(clones_QC(adata_clear, clones_obs, title=f"{FP} clones statisctics"))

        pdf = matplotlib.backends.backend_pdf.PdfPages(savefig)
        for fig in figs:
            pdf.savefig(fig)
        pdf.close()

        concat_figures(savefig)

        if not show:
            for fig in figs:
                plt.close(fig)
            
        return_adata = adata_raw[adata_clear.obs.index]
        return_adata.obs = adata_clear.obs.copy()
        return (return_adata, adatas_FP)
    
def p2conos(adata, batch_key, n_hvgs=3000, n_pcs=30, space="PCA", layer=None, n_cores=30,
            save_path=None, p2app=None):
    """
    This is a wrapper function around conos integration with pagoda2 preprocessing.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples / batches.
    n_hvgs : int, optional
        Number of highly variable genes.
    n_pcs : int, optional
        Number of PCs that will be used for a graph construction.
    space : str, optional
        Space used to establish putative alignments between samples. May be
        "PCA", "CPCF", "JNMF", "CPCA", "CCA" or "PMA".
    n_cores : int, optional
        Number of CPUs used in analysis.
    save_path : str, optional
        Path of RDS file with conos object to save.
    p2app : str, optional
        Path for pagoda2 app saving.

    Returns
    ----------
    AnnData object with a joint graph representation.
    """
    import rpy2.robjects as ro
    import anndata2ri
    import numpy as np
    
    # libraries loading
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(scater))")
    ro.r("suppressPackageStartupMessages(library(conos))")
    ro.r("suppressPackageStartupMessages(library(pagoda2))")
    
    # anndata2ri activation
    anndata2ri.activate()
    
    adata_tmp = adata.copy()
    if not layer is None:
        adata_tmp.X = adata_tmp.layers[layer].copy()
    
    ro.globalenv["adata_tmp"] = adata_tmp
    
    # converting to Seurat
    ro.r('seurat.obj <- as.Seurat(adata_tmp, counts = NULL, data = "X")')
    # splitting dataset per batches
    ro.r(f'panel <- SplitObject(seurat.obj, split.by = "{batch_key}")')
    # getting expressions
    ro.r("""for (i in 1:length(panel)) {
                panel[[i]] = GetAssayData(panel[[i]], slot = "data")
            }""")
    
    # pagoda2 preprocessing
    ro.r(f"""panel <- lapply(panel, basicP2proc, n.cores = {n_cores}, min.cells.per.gene = 0,
                             n.odgenes = {n_hvgs}, get.largevis = FALSE, make.geneknn = FALSE,
                             get.tsne = FALSE)""")
    
    # conos joint graph building
    ro.r('con <- Conos$new(panel, n.cores = 1)')
    ro.r(f"""con$buildGraph(
        space = "{space}",
        ncomps = {n_pcs},
        n.odgenes = {n_hvgs},
        matching.method = "mNN",
        metric = "angular",
        verbose = TRUE
    )""")
    
    # UMAP
    ro.r(f"con$n.cores = {n_cores}")
    ro.r('con$embedGraph(method = "UMAP")')
    ro.r('con$findCommunities()')
    
    if not save_path is None:
        ro.r(f'saveRDS(con, file = "{save_path}")')
        
    if not p2app is None:
        ro.r(f'p2app = p2app4conos(conos = con, file = "{p2app}", save = TRUE)')
    
    # data extraction
    ro.r(f"""
        matrix.merged <- con$getJointCountMatrix(raw = FALSE)       
        cell.ids <- rownames(matrix.merged)
        gene.ids <- colnames(matrix.merged)
        graph.conn <- igraph::as_adjacency_matrix(con$graph, attr = "weight")[cell.ids, cell.ids]
        graph.dist <- graph.conn
        graph.dist@x <- 1 - graph.dist@x 
    """)
    
    adata_tmp = adata_tmp[np.array(ro.r("cell.ids")), np.array(ro.r("gene.ids"))]
    adata_tmp.obsm["X_umap"] = np.array(ro.r("""con$embedding[cell.ids,] %>% as.data.frame()"""))
    adata_tmp.uns["neighbors"] = dict(
        connectivities=ro.r("graph.conn"),
        distances=ro.r("graph.dist")
    )
    adata_tmp.uns["neighbors"]["params"] = dict(
        n_pcs=n_pcs,
        use_rep="X_conos",
        metric="cosine",
        method="umap",
        n_neighbors=15
    )
    adata_tmp.X = ro.r("matrix.merged")
    
    return adata_tmp

def seurat_cca(adata, batch_key):
    """
    This is a wrapper function around Seurat CCA integration.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples / batches.

    Returns
    ----------
    AnnData object with corrected expression values.
    """
    import rpy2.robjects as ro
    import anndata2ri
    
    # libraries loading
    ro.r("suppressPackageStartupMessages(library(Seurat))")
    ro.r("suppressPackageStartupMessages(library(scater))")
    
    # anndata2ri activation
    anndata2ri.activate()
    
    adata_tmp = adata.copy()
    
    ro.globalenv["adata_tmp"] = adata_tmp
    
    # converting to Seurat
    ro.r('seurat.obj <- as.Seurat(adata_tmp, counts = NULL, data = "X")')
    # feature selection is already performed
    ro.r('features <- rownames(seurat.obj)')
    # splitting dataset per batches
    ro.r(f'seurat.objs <- SplitObject(seurat.obj, split.by = "{batch_key}")')
    # integration
    ro.r('anchors <- FindIntegrationAnchors(object.list = seurat.objs, anchor.features = features, verbose = FALSE)')
    ro.r('seurat.objs.combined <- IntegrateData(anchorset = anchors, verbose = FALSE)')
    
    # extraction of corrected expressions
    expression = ro.r('GetAssayData(object = seurat.objs.combined, assay = "integrated", slot = "data")')
    obs_names = ro.r('colnames(seurat.objs.combined)')
    var_names = ro.r('rownames(seurat.objs.combined)')
    
    adata_tmp = adata_tmp[np.array(obs_names), np.array(var_names)]
    adata_tmp.X = expression.T
    
    return adata_tmp

def beautiful_cmap(initial_cmap="Reds", grey_intensity=0.2, color_intencity=0.1):
    """
    Returns color map for visualization of gene expression on UMAPs. Color
    map will starts from grey, not from white.
    
    Parameters
    ----------
    initial_cmap : str, optional
        What color map will be the base for novel color map.
    grey_intensity : float, optional
        What intensity of grey should be at the start of color map.
    color_intencity : float, optional
        What intensity of color should be after grey at color map

    Returns
    ----------
    ListedColormap object.
    
    """
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    import numpy as np
    
    cm_color = cm.get_cmap(initial_cmap, 128)
    cm_grey = cm.get_cmap("Greys", 128)
    
    c = ListedColormap(
        np.vstack(
            (cm_grey(np.linspace(grey_intensity, grey_intensity, 1)),
             cm_color(np.linspace(color_intencity, 1, 128)))
    ))
    
    return c

