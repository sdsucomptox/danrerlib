import os
import toml

def generate_conda_meta(poetry_toml_path, output_dir):
    with open(poetry_toml_path, 'r') as poetry_file:
        poetry_data = toml.load(poetry_file)

    output_path = os.path.join(output_dir, "meta.yaml")

    meta_yaml_content = f"""\
package:
  name: {poetry_data['tool']['poetry']['name']}
  version: {poetry_data['tool']['poetry']['version']}

source:
  path: ..

build:
  number: 0
  noarch: python

requirements:
  build:
    - python
  run:
    - python
"""

    dependencies = poetry_data['tool']['poetry']['dependencies']
    for dep in dependencies:
        meta_yaml_content += f"    - {dep}\n"

    with open(output_path, 'w') as output_file:
        output_file.write(meta_yaml_content)

if __name__ == "__main__":
    os.makedirs("recipe", exist_ok=True)  # Ensure 'recipe' directory exists
    generate_conda_meta("pyproject.toml", "recipe")

