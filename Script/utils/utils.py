import itertools
import math
import matplotlib.pyplot as plt
import scanpy as sc
import snapatac2 as snap
import os
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from tqdm import tqdm

def str_extract(data, pattern):
    import re
    import pandas as pd
    """
    Extracts a part of the string based on a regex pattern.
    :param data: A list or pandas Series containing string elements.
    :param pattern: A regex pattern to search for in each string element.
    :return: A list or pandas Series (depending on the input) of the extracted substrings.
    """
    # Define a helper function to apply regex on each element
    def extract(string):
        match = re.search(pattern, string)
        # Return the matched group if a match is found; otherwise, return None
        return match.group(0) if match else None
    # Apply the helper function to each element in the input
    if isinstance(data, list):
        # If the input is a list, use list comprehension
        return [extract(item) for item in data]
    elif isinstance(data, pd.Series):
        # If the input is a pandas Series, use the .apply() method
        return data.apply(extract).to_list()
    else:
        raise ValueError("Input must be a list or pandas.core.series.Series")



def run_lsi(
    adata,
    n_components=50,
    n_features=None,
    use_snap_feature=True,
    binarize=True,
    drop_first=True,
    random_state=0,
    max_iter=1
):
    """
    TF-IDF + LSI for scATAC

    Parameters
    ----------
    adata : AnnData
        Cell x peak matrix
    n_components : int
        Output LSI dimensions (after drop_first)
    n_features : int or None
        Number of peaks to select
    use_snap_feature : bool
        Use snapatac2.pp.select_features
    binarize : bool
        Whether to binarize matrix
    drop_first : bool
        Drop first LSI component (depth-related)
    random_state : int
    max_iter : int
        Iterations for SnapATAC2 feature selection
    """

    X = adata.X

    # ======================
    # 1. Feature selection
    # ======================
    if n_features is not None:
        if use_snap_feature:
            import snapatac2 as snap

            snap.pp.select_features(
                adata,
                n_features=n_features,
                max_iter=max_iter
            )

            if "selected" in adata.var.columns:
                mask = adata.var["selected"].to_numpy()
            elif "highly_variable" in adata.var.columns:
                mask = adata.var["highly_variable"].to_numpy()
            else:
                raise ValueError(
                    "No 'selected' or 'highly_variable' in adata.var"
                )

            X = X[:, mask]
        else:
            raise ValueError("Only SnapATAC2 feature selection supported")

    # ======================
    # 2. Sparse format
    # ======================
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    # ======================
    # 3. Binarize (重要)
    # ======================
    if binarize:
        X = X.copy()
        X.data[:] = 1.0

    # ======================
    # 4. TF
    # ======================
    cell_sum = np.asarray(X.sum(axis=1)).ravel()
    cell_sum[cell_sum == 0] = 1

    tf = X.multiply(1.0 / cell_sum[:, None])

    # ======================
    # 5. IDF
    # ======================
    peak_sum = np.asarray((X > 0).sum(axis=0)).ravel()
    idf = np.log(1 + X.shape[0] / (peak_sum + 1))

    # ======================
    # 6. TF-IDF
    # ======================
    tfidf = tf.multiply(idf)

    # ======================
    # 7. SVD (LSI)
    # ======================
    svd_dim = n_components + 1 if drop_first else n_components

    svd = TruncatedSVD(
        n_components=svd_dim,
        random_state=random_state
    )

    X_lsi = svd.fit_transform(tfidf)

    # ======================
    # 8. Drop first dim
    # ======================
    if drop_first:
        X_lsi = X_lsi[:, 1:]

    # ======================
    # 9. Save
    # ======================
    adata.obsm["X_lsi"] = X_lsi

    return adata

def plot_lsi_umap_loop(
    adata,
    color="cell_type",
    n_features_list=(50000, 100000, 200000),
    binarize_list=(True, False),
    n_neighbors_list=(10, 30, 50),
    min_dist_list=(0.3,),
    n_components=50,
    drop_first=True,
    random_state=0,
    max_iter=2,
    point_size=8,
    save_dir="lsi_umap_tuning"
):
    """
    Run parameter combinations for LSI + UMAP.
    Plot one figure at a time and save to disk.

    Each plot:
    - shown immediately
    - saved as PNG

    Parameters
    ----------
    save_dir : str
        directory to save figures
    """

    os.makedirs(save_dir, exist_ok=True)

    combos = list(itertools.product(
        n_features_list,
        binarize_list,
        n_neighbors_list,
        min_dist_list
    ))

    print(f"Total combinations: {len(combos)}")

    for i, (n_features, binarize, n_neighbors, min_dist) in enumerate(combos, 1):
        print(f"\n[{i}/{len(combos)}] Running: "
              f"nfeat={n_features}, bin={binarize}, nn={n_neighbors}, md={min_dist}")

        ad = adata.copy()

        try:
            # ---- LSI ----
            ad = run_lsi(
                ad,
                n_components=n_components,
                n_features=n_features,
                use_snap_feature=True,
                binarize=binarize,
                drop_first=drop_first,
                random_state=random_state,
                max_iter=max_iter
            )

            # ---- UMAP ----
            snap.tl.umap(
                ad,
                use_rep="X_lsi",
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                random_state=random_state
            )

            # ---- plot ----
            title = (
                f"nfeat={n_features} | "
                f"bin={binarize} | "
                f"nn={n_neighbors} | "
                f"md={min_dist}"
            )

            fig = sc.pl.umap(
                ad,
                color=color,
                size=point_size,
                title=title,
                return_fig=True,
                show=False
            )

            # ---- save ----
            fname = (
                f"nfeat{n_features}_"
                f"bin{int(binarize)}_"
                f"nn{n_neighbors}_"
                f"md{str(min_dist).replace('.', '')}.png"
            )
            save_path = os.path.join(save_dir, fname)

            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"Saved: {save_path}")

            # ---- show ----
            plt.show()

        except Exception as e:
            print(f"❌ Failed: {e}")

def pseudo_bulk_aggregation(adata, obs_key, agg_func="mean"):
    if obs_key not in adata.obs.columns:
        raise ValueError(f"'{obs_key}' not in adata.obs")

    unique_groups =  adata.obs[obs_key].unique().tolist()
    pseudo_bulk = []
    for group in tqdm(unique_groups):
        subset = adata[adata.obs[obs_key] == group]
        if agg_func == "sum":
            bulk_counts = subset.X.sum(axis=0)
        elif agg_func == "mean":
            bulk_counts = subset.X.mean(axis=0)
        else:
            raise ValueError("agg_func must be 'sum' or 'mean'")

        pseudo_bulk.append(bulk_counts)

    pseudo_bulk_arrary = np.array(pseudo_bulk).squeeze()
    pseudo_bulk_df = pd.DataFrame(pseudo_bulk_arrary.T, 
                                  index=adata.var_names, 
                                  columns=unique_groups)

    return pseudo_bulk_df

def run_umap(adata,nfeatures=2000,batch_key = None,min_dist = 0.3,n_neighbors = 15,n_pcs=20):
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=nfeatures, subset=True)
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    if batch_key is None: 
        sc.pp.neighbors(adata, use_rep='X_pca',n_neighbors=n_neighbors, n_pcs=n_pcs)
    else: 
        sc.external.pp.harmony_integrate(adata, batch_key,basis = 'X_pca',max_iter_harmony = 20)
        sc.pp.neighbors(adata, use_rep='X_pca_harmony',n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata,min_dist=min_dist)
    return(adata)