"""
Enrichment Plots Module
=======================

The Enrichment Plots module provides functions for visualizing results from gene enrichment analyses. These functions generate various types of plots, allowing users to gain insights into the enriched pathways or concepts.

Functions:
    - ``bar_chart``: Generate a bar chart to visualize enrichment analysis results. The chart displays pathway names on the x-axis and -log10(p-value) on the y-axis.
    
    - ``volcano``: Create a volcano plot to visualize the relationship between odds ratio and -log10(p-value) for enriched pathways. The plot highlights significant upregulated and downregulated pathways.

    - ``dotplot``: Generate a dot plot to visualize odds ratio, proportion of significant genes, and p-values for enriched pathways. The size of each dot represents the proportion of significant genes, and color indicates the p-value.

Example:
    To generate a bar chart for enrichment analysis results:

    ```python
    bar_chart(data, title='Enrichment Analysis Results', xlab='Pathways')
    ```

    This example creates a bar chart displaying pathway names on the x-axis and -log10(p-value) on the y-axis.

For detailed information on each function's parameters and usage examples, refer to the documentation. Tutorials demonstrating the use of Enrichment Plots functions are also available.
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings

def bar_chart(data: pd.DataFrame, 
              title: str, 
              xlab: str, 
              ylab: str = '-log10(p-value)', 
              label_font_size: int = 16,
              max_num_pathways: int = None) -> None:
    """
    Generate a bar chart from enrichment analysis results.

    Parameters:
      - `data (pd.DataFrame)`: DataFrame containing enrichment analysis results, including 'Concept Name' and 'P-value'.
      - `title (str)`: Title of the bar chart.
      - `xlab (str)`: Label for the x-axis.
      - `ylab (str, optional)`: Label for the y-axis. Default is '-log10(p-value)'.
      - `label_font_size (int, optional)`: Font size for axis labels and title. Default is 16.
      - `max_num_pathways (Optional[int], optional)`: Maximum number of pathways to display. Default is None (display all).

    Raises:
      - `ValueError`: If required columns ('Concept Name', 'P-value') are missing in the input DataFrame.

    Returns:
      - None: Displays the bar chart.

    Note:
      - The function uses the '-log10(p-value)' for the y-axis.
    """
    
    plt.rcParams['font.family'] = 'sans-serif'  # Set the default font family
    plt.rcParams['font.sans-serif'] = 'Arial, Helvetica, sans-serif'  # Specify a list of font families
    plt.rcParams['font.size'] = 12  # Set the default font size
    
    # Check if required columns exist
    required_columns = {'Concept Name', 'P-value'}
    missing_columns = required_columns - set(data.columns)
    
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
    
    data_to_plot = data.sort_values(by='P-value').head(max_num_pathways)

    plt.figure(figsize=(10, 6))
    plt.bar(data_to_plot['Concept Name'], -np.log10(data_to_plot['P-value']), color='skyblue')
    plt.xlabel(xlab, fontsize=label_font_size)
    plt.ylabel(ylab, fontsize=label_font_size)
    plt.title(title, fontsize=label_font_size)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

def volcano(data: pd.DataFrame, 
            significance_threshold: float = 0.05, 
            top_n: int = 5, 
            fig_height: int = 9, 
            fig_width: int = 8, 
            xmin: float = None, 
            xmax: float = None, 
            plot_all: bool = False, 
            legend_loc: str = 'best', 
            output_filename: str = None) -> None:
    """
    Generate a volcano plot from enrichment analysis results.

    Parameters:
      - `data (pd.DataFrame)`: DataFrame containing enrichment analysis results, including 'Odds Ratio', 'P-value', and 'Concept Name'.
      - `significance_threshold (float, optional)`: Threshold for considering a point as significant. Default is 0.05.
      - `top_n (int, optional)`: Number of top pathways to label based on p-value. Default is 5.
      - `fig_height (int, optional)`: Height of the figure. Default is 9.
      - `fig_width (int, optional)`: Width of the figure. Default is 8.
      - `xmin (float, optional)`: Minimum x-axis limit. Default is None (auto-calculated).
      - `xmax (float, optional)`: Maximum x-axis limit. Default is None (auto-calculated).
      - `plot_all (bool, optional)`: Whether to plot all points (including non-significant). Default is False.
      - `legend_loc (str, optional)`: Location of the legend. Default is 'best' (no legend).
      - `output_filename (str, optional)`: Filepath to save the plot. Default is None (show the plot).

    Returns:
      - None: Displays the volcano plot.

    Note:
      - The function uses '-log10(P-value)' for the y-axis.
    """


    plt.rcParams['font.family'] = 'sans-serif'  # Set the default font family
    plt.rcParams['font.sans-serif'] = 'Arial, Helvetica, sans-serif'  # Specify a list of font families
    plt.rcParams['font.size'] = 15  # Set the default font size

    # Extract data
    odds_ratio = data['Odds Ratio']
    p_value = data['P-value']
    pathway_names = data['Concept Name']

    # Plotting
    plt.figure(figsize=(fig_width, fig_height))

    # Highlight significant points
    significant_upregulated = (odds_ratio > 1) & (p_value < significance_threshold)
    significant_downregulated = (odds_ratio < 1) & (p_value < significance_threshold)
    notsignificant =  p_value > significance_threshold

    plt.scatter(odds_ratio[significant_upregulated], -np.log10(p_value[significant_upregulated]),
                color='red', s=50, label='Upregulated')
    plt.scatter(odds_ratio[significant_downregulated], -np.log10(p_value[significant_downregulated]),
                color='blue', s=50, label='Downregulated')
    if plot_all:
        plt.scatter(odds_ratio[notsignificant], -np.log10(p_value[notsignificant]),
                color='black', alpha = 0.7, s=50, label='other')

    # Add labels for top N upregulated and downregulated pathways based on p-value
    top_upregulated_indices = np.argsort(p_value[significant_upregulated])[:top_n].index
    top_downregulated_indices = np.argsort(p_value[significant_downregulated])[:top_n].index


    for i in top_upregulated_indices:
        plt.text(odds_ratio[i], -np.log10(p_value[i]), pathway_names[i], fontsize=12, ha='left', va='bottom')

    for i in top_downregulated_indices:
        plt.text(odds_ratio[i], -np.log10(p_value[i]), pathway_names[i], fontsize=12,  ha='right',va='bottom')

    # Add labels and title
    plt.title('Volcano Plot of Pathways')
    plt.xlabel('Odds Ratio')
    plt.ylabel('-log10(P-value)')

    # Set x-axis limits based on data
    if xmin and xmax:
        plt.xlim([xmin, xmax])
    else:
        abs_min = np.abs(1-odds_ratio.min())
        abs_max = np.abs(1-odds_ratio.max())

        max_lim = max(abs_min, abs_max)

        x_axis_limits = [1.0-max_lim , 1.1+max_lim ]
        plt.xlim(x_axis_limits)

    if plot_all:
    # Add significance threshold line
        plt.axhline(-np.log10(significance_threshold), color='gray', linestyle='--', linewidth=1, label='Significance Threshold')

    # Add legend
    if legend_loc:
        plt.legend(loc = legend_loc)

    # Save or show the plot
    if output_filename:
        plt.savefig(output_filename, bbox_inches='tight')
    else:
        plt.show()

def dotplot(data_in: pd.DataFrame, 
            num_to_plot: int, 
            direction: str, 
            output_filename: str = None, 
            title: str = None) -> None:
    """
    Generate a dot plot for visualization of pathway enrichment results.

    Parameters:
      - `data_in (pd.DataFrame)`: DataFrame containing pathway enrichment results.
      - `num_to_plot (int)`: Number of pathways to plot.
      - `direction (str)`: Specifies the direction of pathways to be plotted ('up', 'down', or 'both').
      - `output_filename (str, optional)`: Filepath to save the plot. Default is None (show the plot).
      - `title (str, optional)`: Title of the plot. Default is None.

    Returns:
      - None: Displays the dot plot.

    Note:
      - The dot plot visualizes odds ratio, proportion of significant genes, and p-values for each pathway.
      - For 'both' direction, the plot shows both upregulated and downregulated pathways.
    """


    df = data_in.copy()
    if direction == 'up':
        match = 'upregulated'
        label = 'Up'
    elif direction == 'down':
        match = 'downregulated'
        label = 'Down'
    elif direction == 'both':
        label = 'Up and Down'
    else:
        raise ValueError('Invalid direction: up, down, or both.')

    if direction != 'both':
        df = df[df['Direction'] == match]
        df = df.sort_values(by='P-value').head(num_to_plot).reset_index()
    else:
        num_each = round(num_to_plot/2)
        num_to_plot = 2* num_each
        df1 = df[df['Direction'] == 'upregulated']
        df1 = df1.sort_values(by='P-value').head(num_each).reset_index()
        num_up = len(df1)
        df2 = df[df['Direction'] == 'downregulated']
        df2 = df2.sort_values(by='P-value').head(num_each).reset_index()
        num_down = len(df2)
        df = pd.concat([df1, df2], ignore_index=True)
        df.sort_values(by='P-value').reset_index()


    num_to_plot = len(df)

    # Calculate the maximum length of the pathway names
    max_name_length = df['Concept Name'].apply(len).max()
    # Set a base figure size
    if direction == 'both':
        base_width = 3
    else:
        base_width = 2.5
    base_height = 4

    # Adjust the figure width based on the length of the longest pathway name
    width = max(base_width, base_width + max_name_length * 0.1)
    height = base_height + num_to_plot/10

    plt.rcParams['font.family'] = 'sans-serif'  # Set the default font family
    plt.rcParams['font.sans-serif'] = 'Arial, Helvetica, sans-serif'  # Specify a list of font families
    plt.rcParams['font.size'] = 12  # Set the default font size

    # Plotting with adjusted figure size
    fig, ax = plt.subplots(figsize=(width, height))  # Adjust width and height as needed

    # Create scatter plot with varying marker sizes
    scatter = ax.scatter(df['Odds Ratio'], df.index[::-1], s=df['Proportion of Sig Genes in Set'] * 1000 , c=df['P-value'], cmap='viridis', alpha=0.7, vmin=0, vmax=0.05)
    if direction == 'both':
        ax.axvline(x=1, color='black', linestyle='--', alpha = 0.5)


    # Add labels and title
    ax.set_yticks(df.index[::-1])  # Reverse the index order
    ax.set_yticklabels(df['Concept Name'])
    if direction == 'both':
        ax.set_title(f'Top {num_up} Upregulated and \nTop {num_down} Downregulated Pathways')
    elif title:
        ax.set_title(title)
    else:
        ax.set_title(f'Top {num_to_plot} {label}regulated Pathways')
    ax.set_xlabel('Odds Ratio', fontsize=14)

    cbar_ax = fig.add_axes([.84, 0.16, 0.03, 0.35])  # Adjust the values [left, bottom, width, height]

    # Add colorbar for P-value
    cbar = plt.colorbar(scatter, label='P-value', cax = cbar_ax)
    cbar.set_label('P-value', fontsize = 14)
    # cbar = plt.colorbar(scatter, label='P-value', shrink=0.3, pad= 0.12)


    # Create a legend for scatter plot sizes
    handles, labels = scatter.legend_elements(prop="sizes", num = 4, alpha=0.6)
    new_labs = []
    for label in labels:
        newint = int(label.strip('$\\mathdefault{').strip('}$'))/10
        new_labs.append('$\\mathdefault{' + str(int(newint))+ '\%}$')
        
    legend2 = ax.legend(handles, new_labs, loc="upper left", title="Sig Genes\nin Pathway", bbox_to_anchor=(1.01, 0.95))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, message="This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.")
        plt.tight_layout()

    # Save or show the plot
    if output_filename:
        plt.savefig(output_filename, bbox_inches='tight')
    else:
        plt.show()