import pandas as pd

def pretty_print_series(series: pd.Series):
    for element in series:
        print(element)

def save_data(data: pd.Series or pd.DataFrame, location):
    data.to_csv(location, sep = '\t', index=False)
