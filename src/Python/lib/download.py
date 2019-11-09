import os
import requests


def download_if_not_exists(file_path, file_name, file_url, label):
    if not os.path.exists(file_path + file_name):
        print(f"Downloading {label.upper()}...")
        if len(file_path) != 0:
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
        request_result = requests.get(file_url)
        open(file_path + file_name, 'wb').write(request_result.content)
        print(f"{label.upper()} READY")
