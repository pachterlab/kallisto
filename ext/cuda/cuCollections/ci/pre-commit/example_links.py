#!/usr/bin/env python

import os
import json
import git
import re
import argparse
import zlib
import base64
import subprocess
import sys

def open_file(file_path):
    with open(file_path, 'r') as file:
        data = file.read()
    return data

def write_file(file_path, data):
    with open(file_path, 'w') as file:
        file.write(data)

def create_compiler_explorer_url(file_path, COMPILER_ID, COMPILER_FLAGS):
    source_code = open_file(file_path)
    # Prepare the JSON payload
    payload = {
        "sessions": [
            {
                "id": 1,
                "language": "cuda",
                "source": source_code,
                "compilers": [
                    {
                        "id": COMPILER_ID,
                        "options": COMPILER_FLAGS,
                        "libs": [
                            {
                                "id": "cccl",
                                "version": "trunk"
                            },
                            {
                                "id": "cuco",
                                "version": "dev"
                            }
                        ]
                    }
                ]
            }
        ]
    }

    # Convert the payload to JSON
    config_json = json.dumps(payload, separators=(',', ':'))

    # Compress the JSON using zlib (deflate)
    compressed = zlib.compress(config_json.encode('utf-8'))

    # Base64 encode the compressed data
    encoded = base64.urlsafe_b64encode(compressed).decode('utf-8')

    # Replace problematic unicode characters
    # As per the documentation, map characters in the range \u007F-\uFFFF
    encoded = ''.join(
        c if '\u0000' <= c <= '\u007F' else '\\u{:04x}'.format(ord(c))
        for c in encoded
    )

    # Construct the final URL
    url = f'https://godbolt.org/clientstate/{encoded}'

    return url

def update_readme_with_urls(readme_path, url_table, args):
    # Open and read the README.md file
    readme_content = open_file(readme_path)
    lines = readme_content.splitlines()

    updated_lines = []
    mismatched_files = []
    missing_files = []  # List of files not found in README.md

    # Iterate through each line of the README
    for line in lines:
        updated_line = line
        # Check each file in the url_table
        for file_name, url in url_table.items():
            # Check if the file name exists in the line
            if file_name in line:
                # Look for the specific godbolt link format and replace the URL
                updated_line = re.sub(
                    r"\(see \[live example in godbolt\]\(https?://[^\)]*\)\)",
                    f"(see [live example in godbolt]({url}))",
                    line
                )
                # Also handle empty godbolt links
                if updated_line == line:
                    updated_line = re.sub(
                        r"\(see \[live example in godbolt\]\(\)\)",
                        f"(see [live example in godbolt]({url}))",
                        line
                    )
                if updated_line != line:
                    mismatched_files.append(file_name)

        updated_lines.append(updated_line)

    # After updating lines, check if any example files are missing from README
    for file_name in url_table.keys():
        # If the file name was not found in the README, add it to missing_files
        if not any(file_name in line for line in lines):
            missing_files.append(file_name)

    # Join the lines and write the updated content back to README.md
    updated_content = "\n".join(updated_lines)

    if missing_files:
        print(f"Error: The following example files are missing in README.md: {missing_files}")
        sys.exit("FAILED: README.md is missing entries for the above example files")

    if args.fix:
        write_file(readme_path, updated_content)
    else:
        if mismatched_files:
            write_file(readme_path, updated_content)
            print(mismatched_files)
            sys.exit("FAILED: the above examples have out-of-date Godbolt links")

COMPILER_ID = 'nvcc125u1'  # Desired compiler ID
COMPILER_FLAGS = '-std=c++17 -arch=sm_70 --expt-extended-lambda'  # NVCC compiler flags
EXAMPLES_DIR = "examples"  # Path to the examples directory relative to the repo root

def get_changed_cuda_files():
    # Run the git command to get the changed files
    result = subprocess.run(
        ['git', 'diff', '--cached', '--name-only', '--diff-filter=ACM'],
        capture_output=True,
        text=True
    )

    # Check if the command was successful
    if result.returncode != 0:
        print("Error executing git command:", result.stderr)
        return []

    # Get the output and filter for .cu files in EXAMPLES_DIR
    return result.stdout.splitlines()

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Ensure CUDA example files have up-to-date Godbolt links in the README.")
    parser.add_argument('--fix', action='store_true', 
                        help="Automatically update mismatched or outdated Godbolt links in the README with correct URLs.")
    args = parser.parse_args()

    # Initialize the Git repository
    repo = git.Repo(search_parent_directories=True)
    repo_root = repo.git.rev_parse("--show-toplevel")
    
    # Ensure we are in the ci directory
    ci_dir = os.path.join(repo_root, 'ci')
    os.chdir(ci_dir)

    # Initialize a hash table (dictionary) to store file urls
    url_table = {}

    # Get all CUDA files in the examples directory
    example_files = []
    for root, _, files in os.walk(os.path.join(repo_root, EXAMPLES_DIR)):
        for file_name in files:
            if file_name.endswith('.cu'):
                example_files.append(os.path.relpath(os.path.join(root, file_name), repo_root))

    # Iterate through the example files and create URLs
    for file_path in example_files:
        full_file_path = os.path.join(repo_root, file_path)
        url = create_compiler_explorer_url(full_file_path, COMPILER_ID, COMPILER_FLAGS)
        if url:
            url_table[file_path] = url  # Store the file and its url in the hash table
            # print(f"File: {file_path}: url: {url}")

    # Update README.md with the urls
    readme_path = os.path.join(repo_root, "README.md")
    update_readme_with_urls(readme_path, url_table, args)

if __name__ == "__main__":
    main()
