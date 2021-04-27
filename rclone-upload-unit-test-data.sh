#!/usr/bin/env bash

REMOTE="autometa-test-data"

unit_test_data_files_path="$HOME/Autometa/unit-test-data-files.txt"
test_data_dir="$HOME/Autometa/tests/data"
for filename in `cat $unit_test_data_files_path`;do
    filepath="${test_data_dir}/${filename}"
    # rclone copyto --dry-run --progress $filepath $REMOTE:unit_test_data/${filename}
    echo "Copying $filename"
    rclone copyto --progress $filepath $REMOTE:unit_test_data/${filename}
done