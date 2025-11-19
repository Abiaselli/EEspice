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
    
    try:
        with open(in_path, 'r') as f_in, open(out_path, 'w') as f_out:
            for line in f_in:
                # Strip whitespace and split by whitespace (tabs or spaces)
                parts = line.strip().split()
                
                # Ensure line has at least 3 parts (row, col, val)
                if len(parts) >= 3:
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
        
        print(f"Success! Processed {count} lines.")

    except Exception as e:
        print(f"An error occurred: {e}")

# ---------------------------------------------------------
# 3. Execution
# ---------------------------------------------------------

if __name__ == "__main__":
    convert_to_zero_based(input_filename, output_filename)