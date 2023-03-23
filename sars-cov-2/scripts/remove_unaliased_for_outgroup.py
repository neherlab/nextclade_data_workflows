import json

def load_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def save_json_file(file_path, data):
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=2)

def update_data(data):
    if isinstance(data, dict):
        if not "Nextclade_pango" in data and "partiallyAliased" in data:
            if "value" in data["partiallyAliased"]: 
                data["partiallyAliased"] = {"value": ""}
        for key in data:
            update_data(data[key])
    elif isinstance(data, list):
        for item in data:
            update_data(item)

def main():
    file_path_in = 'sars-cov-2/auspice/21L/auspice.json'
    file_path_out = 'test2.json'
    data = load_json_file(file_path_in)
    update_data(data)
    save_json_file(file_path_out, data)

if __name__ == "__main__":
    main()