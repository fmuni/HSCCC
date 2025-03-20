import sys
import numpy as np
import matplotlib.pyplot as plt

# Check if correct number of arguments is provided
if len(sys.argv) < 3:
    print("Usage: python3 plotFile.py <file_to_plot> <column_to_plot> [<additional_columns>...]")
    sys.exit(1)

# Get file path and columns to plot from arguments
file_path = sys.argv[1]
column_indices = [int(i) - 1 for i in sys.argv[2:]]  # Convert to zero-based index

# Read the file and extract data
with open(file_path, "r") as file:
    lines = file.readlines()

# Find the last commented line before data starts (which contains column names)
header_line = None
for line in reversed(lines):
    if line.startswith("#"):
        header_line = line.strip("#").strip().split()  # Extract column names
        break

# Ensure we have a valid header line
if header_line is None or len(header_line) < 2:
    print("Error: Could not find a valid header line with column names.")
    sys.exit(1)

headers = header_line[1:]  # Skip the first entry (Time)

# Lists to store data
data = {header: [] for header in headers}
time = []

# Parse data
for line in lines:
    if line.startswith("#") or not line.strip():
        continue  # Skip header and empty lines
    values = line.split()
    time.append(float(values[0]))  # First column is time
    for i, header in enumerate(headers):
        data[header].append(float(values[i + 1]))

# Convert lists to numpy arrays
time = np.array(time)
for key in data:
    data[key] = np.array(data[key])

# Validate column indices
if any(i >= len(headers) for i in column_indices):
    print(f"Error: Column index out of range. Available columns are 1 to {len(headers)}.")
    sys.exit(1)

# Plot the selected columns
plt.figure(figsize=(8, 6))
for i in column_indices:
    plt.plot(time, data[headers[i]], label=headers[i])

plt.xlabel("Time")
plt.ylabel("Values")
#plt.title("Time Evolution of Selected Quantities")
plt.legend()
plt.grid(True)
plt.show()

