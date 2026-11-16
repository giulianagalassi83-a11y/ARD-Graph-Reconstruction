# ARD Graph Reconstruction and Clustering

This repository contains the code and data used for graph reconstruction and clustering experiments based on Aggregated Relational Data (ARD).

The repository includes three datasets: **Facebook**, **Cora**, and **Indian Villages**. For each dataset, we provide the Jupyter notebooks implementing each step of the methodology. The data files are also included for Facebook and Cora.

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
    └── paso4_clustering_indian.ipynb
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

The **Facebook** and **Cora** folders each contain a `data/` subdirectory with the corresponding dataset files, which can be downloaded directly from the repository.

The **Indian Villages** folder contains the experimental notebooks but does not include a `data/` subdirectory.

## Execution

To reproduce the experiments, download or clone the repository and run the notebooks for the selected dataset in numerical order, starting from `paso0_ARD_matrix` and finishing with `paso4_clustering`.

The Facebook and Cora experiments use the data provided in their respective `data/` folders. The Indian Villages experiments require the corresponding data files to be obtained separately.

Each dataset follows an independent experimental pipeline.
