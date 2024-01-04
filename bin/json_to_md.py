import json
import argparse

# def generate_markdown(json_schema, file_name):
#     with open(file_name, 'w') as f:
#         for definition, details in json_schema['definitions'].items():
#             f.write(f"## {details['title']}\n")
#             f.write(f"{details['description']}\n\n")
#             for property, attributes in details['properties'].items():
#                 f.write(f"### {property}\n")
#                 f.write(f"{attributes['description']}\n")
#                 if 'default' in attributes:
#                     f.write(f"Default: {attributes['default']}\n")
#                 f.write("\n")
def generate_markdown(json_schema, file_name):
    with open(file_name, 'w') as f:
        for definition, details in json_schema['definitions'].items():
            f.write(f"## {details['title']}\n")
            f.write(f"{details['description']}\n\n")
            f.write("| Parameter | Description | Default |\n")
            f.write("| --- | --- | --- |\n")
            for property, attributes in details['properties'].items():
                default = attributes.get('default', '')
                f.write(f"| {property} | {attributes['description']} | {default} |\n")
            f.write("\n")
            
parser = argparse.ArgumentParser(description='Generate markdown from JSON schema.')
parser.add_argument('json_file', type=str, help='Input JSON schema file')
parser.add_argument('md_file', type=str, help='Output markdown file')

args = parser.parse_args()

with open(args.json_file) as json_file:
    data = json.load(json_file)
    generate_markdown(data, args.md_file)