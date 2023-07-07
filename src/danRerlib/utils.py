import pandas as pd

def pretty_print_series(series: pd.Series):
    for element in series:
        print(element)

def save_series(series: pd.Series, location):
    series.to_csv(location, sep = '\t', index=False)