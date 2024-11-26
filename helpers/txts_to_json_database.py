# Combine all text files within the specified folder_path into one JSON to be used as database.

import os
import json

text_files_dict = {}
folder_path = 'your/folder/path'  # TODO: Specify your folder path

for root, dirs, files in os.walk(folder_path):
    for file in files:
        if file.endswith(".txt"):
            file_name = os.path.splitext(file)[0]  # Remove the .txt extension
            file_path = os.path.join(root, file)
            with open(file_path, 'r', encoding='utf-8') as f:
                # Strip whitespaces from each line
                text_files_dict[file_name] = [line.strip() for line in f.readlines()]

# Saving the dictionary as a JSON file in the folder path
json_output_path = os.path.join(folder_path, 'text_files_dict.json')
with open(json_output_path, 'w', encoding='utf-8') as json_file:
    json.dump(text_files_dict, json_file, indent=4)

print(f"JSON saved at: {json_output_path}")