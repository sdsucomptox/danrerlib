import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def bar_chart(data, title, xlab, 
              ylab = '-log10(p-value)', 
              label_font_size = 16,
              max_num_pathways = None):
    
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