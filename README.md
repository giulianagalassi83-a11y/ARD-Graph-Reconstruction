# ARD Graph Reconstruction and Clustering

This repository contains the code and data used for graph reconstruction and clustering experiments based on Aggregated Relational Data (ARD).

The repository includes three datasets: **Facebook**, **Cora**, and **Indian Villages**. For each dataset, we provide the data required to reproduce the experiments, along with the Jupyter notebooks implementing each step of the methodology.

## Repository Structure

The repository is organized into three main folders, one for each dataset:

```text
ARD-Graph-Reconstruction/
│
├── Cora/
│   ├── paso0_ARD_matrix_cora.ipynb
│   ├── paso1_synthetic_graphs_cora.ipynb
│   ├── paso2_link_prediction_cora.ipynb
│   ├── paso3_DGI_cora.ipynb
│   ├── paso4_clustering_cora.ipynb
│   └── data/
│
├── Facebook/
│   ├── paso0_ARD_matrix_facebook.ipynb
│   ├── paso1_synthetic_graphs_facebook.ipynb
│   ├── paso2_link_prediction_facebook.ipynb
│   ├── paso3_DGI_facebook.ipynb
│   ├── paso4_clustering_facebook.ipynb
│   └── data/
│
└── Indian Villages/
    ├── paso0_ARD_matrix_indian.ipynb
    ├── paso1_synthetic_graphs_indian.ipynb
    ├── paso2_link_prediction_indian.ipynb
    ├── paso3_DGI_indian.ipynb
    ├── paso4_clustering_indian.ipynb
    └── data/
```

## Methodology

Each dataset follows the same five-step experimental pipeline:

| Notebook | Description |
|---|---|
| `paso0_ARD_matrix` | Construction of the Aggregated Relational Data (ARD) matrix. |
| `paso1_synthetic_graphs` | Generation of synthetic graphs from ARD information. |
| `paso2_link_prediction` | Link prediction and evaluation of graph reconstruction. |
| `paso3_DGI` | Generation of node embeddings using Deep Graph Infomax (DGI). |
| `paso4_clustering` | Clustering analysis based on the learned node embeddings. |

## Data

Each dataset folder contains a `data/` subdirectory with the files required to run the corresponding notebooks.

The data can be downloaded directly from the repository.

## Execution

To reproduce the experiments, download or clone the repository and run the notebooks for the selected dataset in numerical order, starting from `paso0_ARD_matrix.ipynb` and finishing with `paso4_clustering.ipynb`.

Each dataset can be processed independently using its corresponding notebooks and data files.
