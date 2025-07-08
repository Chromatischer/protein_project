"""
Visualization module for protein data.

This module provides plotting functions for visualizing protein data,
including heatmaps and other visualization types using Plotly.
"""

import plotly.graph_objects as go


def renderPlot(data, x, y):
    """
    Render a heatmap plot for protein data visualization.

    This function creates an interactive heatmap using Plotly to visualize
    protein expression data across different samples. The plot size is
    automatically calculated based on the data dimensions.

    Args:
        data (pandas.DataFrame): The data matrix to visualize as a heatmap
        x (list): Column labels for the x-axis (sample names)
        y (list): Row labels for the y-axis (protein identifiers)

    Returns:
        None: Displays the plot directly

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame([[1, 2], [3, 4]],
        ...                     columns=['Sample1', 'Sample2'],
        ...                     index=['Protein1', 'Protein2'])
        >>> renderPlot(data, data.columns, data.index)
    """
    # Calculate plot dimensions based on data size
    num_rows = len(data.index)
    num_cols = len(data.columns)
    cell_size = 20  # pixels per cell
    margin = 200  # space for labels, colorbar, etc.

    width = num_cols * cell_size + margin
    height = num_rows * cell_size + margin

    # Create the heatmap figure
    fig = go.Figure(
        go.Heatmap(
            z=data,
            x=x,
            y=y,
            colorscale="Viridis",
        )
    )

    # Update layout with calculated dimensions and aspect ratio
    fig.update_layout(
        width=width,
        height=height,
        xaxis=dict(
            scaleanchor="y",
            scaleratio=1,
        ),
        yaxis=dict(
            scaleanchor="x",
        ),
    )

    # Display the plot
    fig.show()
