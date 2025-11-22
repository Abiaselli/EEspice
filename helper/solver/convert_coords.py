import os
import sys

# ---------------------------------------------------------
# 1. Setup: Get filenames from arguments or use defaults
# ---------------------------------------------------------
# Usage: python convert_coords.py <input_file> <output_file>

if len(sys.argv) == 3:
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
else:
    # Default filenames if no arguments are provided
    print("No arguments provided. Using default filenames.")
    print("Usage: python convert_coords.py <input_file> <output_file>")
    input_filename = "matrix_data_1based.txt"
    output_filename = "matrix_data_armadillo.txt"

# ---------------------------------------------------------
# 2. The Conversion Logic
# ---------------------------------------------------------
def convert_to_zero_based(in_path, out_path):
    if not os.path.exists(in_path):
        print(f"Error: Input file '{in_path}' not found.")
        return

    print(f"Reading from: {in_path}")
    print(f"Writing to:   {out_path}")

    count = 0
    num_columns = None
    columns_dropped = False
    inconsistent_format_warnings = 0

    try:
        with open(in_path, 'r') as f_in, open(out_path, 'w') as f_out:
            for line in f_in:
                # Strip whitespace and split by whitespace (tabs or spaces)
                parts = line.strip().split()

                # Ensure line has at least 3 parts (row, col, val)
                if len(parts) >= 3:
                    # Detect column count from first valid line
                    if num_columns is None:
                        num_columns = len(parts)
                        print(f"\nFormat detected: {num_columns} column(s) per line")

                        if num_columns == 3:
                            print("✓ Correct format: 3 columns (row, col, value)")
                        elif num_columns > 3:
                            columns_dropped = True
                            print(f"⚠ Warning: {num_columns} columns detected, but only 3 are needed.")
                            print(f"  Column(s) 4-{num_columns} will be ignored and removed from output.")
                            print(f"  Output will contain only: row, col, value\n")

                    # Check for inconsistent column counts
                    if len(parts) != num_columns and inconsistent_format_warnings < 5:
                        print(f"⚠ Warning: Line {count + 1} has {len(parts)} columns (expected {num_columns})")
                        inconsistent_format_warnings += 1
                        if inconsistent_format_warnings == 5:
                            print("  (Further format warnings will be suppressed)")

                    try:
                        # Parse indices as integers
                        row_idx = int(parts[0])
                        col_idx = int(parts[1])

                        # Keep value as string to ensure precision is exactly preserved
                        val_str = parts[2]

                        # Perform the shift (1-based -> 0-based)
                        new_row = row_idx - 1
                        new_col = col_idx - 1

                        # Write to new file (using tabs or spaces)
                        f_out.write(f"{new_row}\t{new_col}\t{val_str}\n")
                        count += 1

                    except ValueError:
                        print(f"Skipping malformed line (non-integer index): {line.strip()}")

        print(f"\n{'='*60}")
        print(f"Success! Processed {count} lines.")
        if columns_dropped:
            print(f"✓ Removed column(s) 4-{num_columns} from all lines")
        print(f"✓ Converted indices from 1-based to 0-based")
        print(f"{'='*60}")

    except Exception as e:
        print(f"An error occurred: {e}")

# ---------------------------------------------------------
# 3. Execution
# ---------------------------------------------------------

if __name__ == "__main__":
    convert_to_zero_based(input_filename, output_filename)