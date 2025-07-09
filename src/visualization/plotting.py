"""
Visualization module for protein data.

This module provides plotting functions for visualizing protein data,
including heatmaps and other visualization types using Plotly.
"""

import plotly.graph_objects as go
import plotly.colors as colors
from copy import deepcopy
from pandas import DataFrame


def renderPlot(
    data: DataFrame,
    col_names: list[str],
    row_names: list[str],
    clusters_dict: dict = None,
    cluster_labels: list[str] = None,
):
    """
    Render a heatmap plot for protein data visualization with cluster
    grouping information.

    This function creates an interactive heatmap using Plotly to visualize
    protein expression data across different samples. When cluster information
    is provided, it adds visual indicators to show cluster groupings including
    separator lines, color-coded labels, and cluster annotations.

    Args:
        data (pandas.DataFrame): The data matrix to visualize as a heatmap
        clusters_dict (dict, optional): Dictionary mapping cluster IDs to
            lists of protein names
        cluster_labels (list, optional): List of cluster labels corresponding
            to each row

    Returns:
        None: Displays the plot directly

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame([[1, 2], [3, 4]],
        ...                     columns=['Sample1', 'Sample2'],
        ...                     index=['Protein1', 'Protein2'])
        >>> clusters = {1: ['Protein1'], 2: ['Protein2']}
        >>> renderPlot(data, data.columns, data.index, clusters_dict=clusters)
    """
    # Calculate plot dimensions based on data size
    num_rows = len(data.index)
    num_cols = len(data.columns)
    cell_size = 20  # pixels per cell
    margin = 300  # increased margin for cluster annotations

    width = num_cols * cell_size + margin
    height = num_rows * cell_size + margin

    # Create the heatmap figure
    fig = go.Figure(
        go.Heatmap(
            z=data,
            x=col_names,
            y=row_names,
            colorscale="Viridis",
        )
    )

    # Add cluster visualization if cluster information is provided
    if clusters_dict is not None:
        _add_cluster_annotations(
            fig, row_names, clusters_dict, num_rows, num_cols)

    # Update layout with calculated dimensions and styling
    title = (
        "Protein Expression Heatmap with Cluster Groupings"
        if clusters_dict
        else "Protein Expression Heatmap"
    )

    fig.update_layout(
        width=width,
        height=height,
        xaxis=dict(
            scaleanchor="y",
            scaleratio=1,
            side="bottom",
        ),
        yaxis=dict(
            scaleanchor="x",
            side="left",
            autorange="reversed",
        ),
        title=title,
        title_x=0.5,
    )

    fig.update_yaxes(range=[len(row_names) - 0.5, -0.5])

    # Display the plot
    print("\n---Create deepcopy of fig and displaying!---")
    deepcopy(fig).show()


def _add_cluster_annotations(fig, row_names, clusters_dict, num_rows, num_cols):
    """
    Add cluster visualization elements to the heatmap figure.

    Args:
        fig: Plotly figure object
        row_names: List of row names (protein identifiers)
        clusters_dict: Dictionary mapping cluster IDs to protein lists
        num_rows: Number of rows in the heatmap
        num_cols: Number of columns in the heatmap
    """
    # Create mapping from protein name to cluster ID
    protein_to_cluster = {}
    for cluster_id, proteins in clusters_dict.items():
        for protein in proteins:
            protein_to_cluster[protein] = cluster_id

    # Generate distinct colors for each cluster
    cluster_ids = list(clusters_dict.keys())
    cluster_colors = _generate_cluster_colors(len(cluster_ids))
    cluster_color_map = {
        cluster_id: cluster_colors[i] for i, cluster_id in enumerate(cluster_ids)
    }

    # Add cluster separator lines and annotations
    current_cluster = None
    cluster_boundaries = []

    for i, protein in enumerate(row_names):
        cluster_id = protein_to_cluster.get(protein)

        # Detect cluster boundaries
        if current_cluster is not None and cluster_id != current_cluster:
            # Add separator line between clusters
            boundary_y = num_rows - i - 0.5
            cluster_boundaries.append(
                (boundary_y, current_cluster, cluster_id))

            fig.add_shape(
                type="line",
                x0=-0.5,
                x1=num_cols - 0.5,
                y0=boundary_y,
                y1=boundary_y,
                line=dict(color="white", width=3),
                layer="above",
            )

        current_cluster = cluster_id

    # Add cluster ID annotations on the right side
    _add_cluster_labels(
        fig, row_names, protein_to_cluster, cluster_color_map, num_rows, num_cols
    )


def _generate_cluster_colors(num_clusters):
    """
    Generate distinct colors for each cluster.

    Args:
        num_clusters: Number of clusters to generate colors for

    Returns:
        List of color strings
    """
    if num_clusters <= 10:
        # Use qualitative color palette for small number of clusters
        return colors.qualitative.Set3[:num_clusters]
    else:
        # Generate colors using HSV space for larger number of clusters
        import colorsys

        cluster_colors = []
        for i in range(num_clusters):
            hue = i / num_clusters
            rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.9)
            hex_color = "#%02x%02x%02x" % (
                int(rgb[0] * 255),
                int(rgb[1] * 255),
                int(rgb[2] * 255),
            )
            cluster_colors.append(hex_color)
        return cluster_colors


def _add_cluster_labels(
    fig, row_names, protein_to_cluster, cluster_color_map, num_rows, num_cols
):
    """
    Add cluster ID labels and color coding to the plot.

    Args:
        fig: Plotly figure object
        row_names: List of row names
        protein_to_cluster: Dictionary mapping proteins to cluster IDs
        cluster_color_map: Dictionary mapping cluster IDs to colors
        num_rows: Number of rows in the heatmap
        num_cols: Number of columns in the heatmap
    """
    # Group consecutive proteins by cluster for efficient labeling
    cluster_groups = []
    current_cluster = None
    current_group_start = 0

    for i, protein in enumerate(row_names):
        cluster_id = protein_to_cluster.get(protein)

        if cluster_id != current_cluster:
            if current_cluster is not None:
                # End previous group
                cluster_groups.append(
                    {
                        "cluster_id": current_cluster,
                        "start_idx": current_group_start,
                        "end_idx": i - 1,
                        "color": cluster_color_map.get(current_cluster, "#000000"),
                    }
                )

            # Start new group
            current_cluster = cluster_id
            current_group_start = i

    # Don't forget the last group
    if current_cluster is not None:
        cluster_groups.append(
            {
                "cluster_id": current_cluster,
                "start_idx": current_group_start,
                "end_idx": len(row_names) - 1,
                "color": cluster_color_map.get(current_cluster, "#000000"),
            }
        )

    # Add cluster labels on the right side
    for group in cluster_groups:
        # Calculate middle position of the cluster group
        start_y = num_rows - group["end_idx"] - 1
        end_y = num_rows - group["start_idx"] - 1
        middle_y = (start_y + end_y) / 2

        # Add cluster label annotation
        fig.add_annotation(
            x=num_cols + 0.5,
            y=middle_y,
            text=f"Cluster {group['cluster_id']}",
            showarrow=False,
            font=dict(color=group["color"], size=12, family="Arial Black"),
            xanchor="left",
            yanchor="middle",
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor=group["color"],
            borderwidth=1,
        )

        # Add colored rectangle to highlight cluster region
        fig.add_shape(
            type="rect",
            x0=num_cols - 0.3,
            x1=num_cols - 0.1,
            y0=start_y - 0.4,
            y1=end_y + 0.4,
            fillcolor=group["color"],
            opacity=0.7,
            line=dict(width=0),
            layer="above",
        )
