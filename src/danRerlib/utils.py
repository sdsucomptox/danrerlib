import pandas as pd
import urllib.request

def pretty_print_series(series: pd.Series):
    for element in series:
        print(element)

def save_data(data: pd.Series or pd.DataFrame, location):
    data.to_csv(location, sep = '\t', index=False)

def download_url(url_to_download: str) -> list:
    response = None

    try:
        request = urllib.request.Request(url_to_download)
        response = urllib.request.urlopen(request)

        return response.readlines()

    except urllib.error.HTTPError as e:
        print('Failed to download contents of URL')
        print('Status code: {}'.format(e.code))
        print()

    finally:
        if response != None:
            response.close()

    return []

def download_url_to_file(url_to_download: str, output_file_path: str) -> None:
    response = None
    file_to_save = None

    try:
        request = urllib.request.Request(url_to_download)
        response = urllib.request.urlopen(request)

        file_to_save = open(output_file_path, 'wb')
        file_to_save.write(response.read())

    except urllib.error.HTTPError as e:
        print('Failed to download contents of URL')
        print('Status code: {}'.format(e.code))
        print()

    finally:
        if file_to_save != None:
            file_to_save.close()
        
        if response != None:
            response.close()