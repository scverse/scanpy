```mermaid

graph LR

    AnnData_Data_Model["AnnData Data Model"]

    Data_Management_Access["Data Management & Access"]

    Core_Analysis_Pipeline["Core Analysis Pipeline"]

    Visualization_Engine["Visualization Engine"]

    Extension_Utilities["Extension & Utilities"]

    Data_Management_Access -- "reads data into and writes data from" --> AnnData_Data_Model

    Core_Analysis_Pipeline -- "transforms and stores processed data, embeddings, and analysis results within" --> AnnData_Data_Model

    Visualization_Engine -- "retrieves data and analysis results from" --> AnnData_Data_Model

    Extension_Utilities -- "orchestrates operations on" --> AnnData_Data_Model

    Extension_Utilities -- "may interact with or modify" --> AnnData_Data_Model

    Data_Management_Access -- "utilizes common services (e.g., file handling, logging, configuration) provided by" --> Extension_Utilities

    Core_Analysis_Pipeline -- "provides processed data and analysis results (e.g., embeddings, cluster assignments) that are visualized by" --> Visualization_Engine

    Core_Analysis_Pipeline -- "relies on general utilities (e.g., numerical operations, logging, configuration) from" --> Extension_Utilities

    Visualization_Engine -- "utilizes services (e.g., plot configuration, saving, logging) from" --> Extension_Utilities

    Extension_Utilities -- "orchestrates workflows involving" --> Data_Management_Access

    Extension_Utilities -- "orchestrates workflows involving" --> Core_Analysis_Pipeline

    Extension_Utilities -- "orchestrates workflows involving" --> Visualization_Engine

    Data_Management_Access -- "accesses common services from" --> Extension_Utilities

    Core_Analysis_Pipeline -- "accesses common services from" --> Extension_Utilities

    Visualization_Engine -- "accesses common services from" --> Extension_Utilities

    click AnnData_Data_Model href "https://github.com/scverse/scanpy/blob/main/.codeboarding//AnnData_Data_Model.md" "Details"

    click Data_Management_Access href "https://github.com/scverse/scanpy/blob/main/.codeboarding//Data_Management_Access.md" "Details"

    click Core_Analysis_Pipeline href "https://github.com/scverse/scanpy/blob/main/.codeboarding//Core_Analysis_Pipeline.md" "Details"

```



[![CodeBoarding](https://img.shields.io/badge/Generated%20by-CodeBoarding-9cf?style=flat-square)](https://github.com/CodeBoarding/GeneratedOnBoardings)[![Demo](https://img.shields.io/badge/Try%20our-Demo-blue?style=flat-square)](https://www.codeboarding.org/demo)[![Contact](https://img.shields.io/badge/Contact%20us%20-%20contact@codeboarding.org-lightgrey?style=flat-square)](mailto:contact@codeboarding.org)



## Details



The `scanpy` architecture is designed as a modular, data-centric scientific computing library, with the `AnnData Data Model` at its core. This model serves as the central hub for all single-cell omics data, around which specialized modules for data handling, analysis, and visualization operate. A dedicated `Extension & Utilities` layer provides cross-cutting concerns and extensibility.



### AnnData Data Model [[Expand]](./AnnData_Data_Model.md)

The foundational in-memory data structure (anndata.AnnData) that serves as the central repository for all single-cell omics data, including expression matrices, cell and gene annotations, dimensionality reduction embeddings, and neighborhood graphs. It is the primary object manipulated by all other modules.





**Related Classes/Methods**:



- `anndata.AnnData`





### Data Management & Access [[Expand]](./Data_Management_Access.md)

Manages all data input/output operations, including reading from and writing to various file formats (e.g., H5AD, 10x Genomics). It also provides access to curated example datasets and utilities for querying and retrieving specific data subsets from the AnnData object.





**Related Classes/Methods**:



- <a href="https://github.com/scverse/scanpy/blob/main/src/scanpy/readwrite.py" target="_blank" rel="noopener noreferrer">`scanpy.readwrite`</a>

- `scanpy.datasets`

- `scanpy.get`

- `scanpy.queries`





### Core Analysis Pipeline [[Expand]](./Core_Analysis_Pipeline.md)

A comprehensive suite of modules responsible for the entire single-cell data analysis workflow. This includes preprocessing (quality control, normalization, scaling, PCA), constructing neighborhood graphs, performing dimensionality reduction (UMAP, t-SNE), clustering, trajectory inference, differential expression analysis, gene scoring, and spatial metrics calculation.





**Related Classes/Methods**:



- `scanpy.preprocessing`

- `scanpy.neighbors`

- `scanpy.tools`

- `scanpy.metrics`





### Visualization Engine

Provides a robust framework and specific functions for generating a wide array of plots to visualize single-cell data and analysis results. This includes core utilities for plot setup, handling common parameters, and high-level functions for scatter plots, heatmaps, dot plots, and specialized trajectory visualizations.





**Related Classes/Methods**:



- `scanpy.plotting`





### Extension & Utilities

A foundational layer providing cross-cutting concerns and extensibility. It includes interfaces for integrating external single-cell analysis tools, a command-line interface for workflow automation, and internal system utilities for configuration management, logging, and general helper functions used across the library.





**Related Classes/Methods**:



- `scanpy.external`

- <a href="https://github.com/scverse/scanpy/blob/main/src/scanpy/cli.py" target="_blank" rel="noopener noreferrer">`scanpy.cli`</a>

- `scanpy._settings`

- <a href="https://github.com/scverse/scanpy/blob/main/src/scanpy/logging.py" target="_blank" rel="noopener noreferrer">`scanpy.logging`</a>

- `scanpy._utils`









### [FAQ](https://github.com/CodeBoarding/GeneratedOnBoardings/tree/main?tab=readme-ov-file#faq)