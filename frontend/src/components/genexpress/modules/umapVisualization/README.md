# UMAP Visualization Module

This module provides UMAP (Uniform Manifold Approximation and Projection) visualization for single-cell data stored in the Django REST framework storage.

## Features

- Fetches UMAP data from storage API endpoints
- Interactive scatter plot visualization using Vega
- Cell selection and highlighting
- Similar styling to differential expression plots
- Responsive design that works with the grid layout system

## Usage

The module is automatically included in the main application grid. It fetches storage data based on the configured storage ID and displays UMAP coordinates as an interactive scatter plot.

## Data Format

The module expects storage data in the following format:

```json
{
    "id": 37,
    "name": "Storage for data id 19",
    "json": {
        "cells": {
            "cell_0": {
                "umap_1": 0.893747524883914,
                "umap_2": 7.595740807670293
            },
            "cell_1": {
                "umap_1": -3.909105958208723,
                "umap_2": 4.754526475663017
            }
        },
        "metadata": {
            "n_cells": 100,
            "umap_key": "X_uce_umap"
        }
    }
}
```

## Components

- `UmapVisualization`: Main component that handles data fetching and UI
- `UmapScatterPlot`: Vega-based scatter plot visualization
- `umapVisualization.styles.ts`: Styled components for layout

## Integration

The module is integrated into the main application through:
- `ModulesKeys.umapVisualization` in constants
- Layout configuration in Redux store
- Grid component inclusion in `GeneExpressGrid`
