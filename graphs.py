import os
import matplotlib.pyplot as plt
# TODO: total aa, unique aa, unique bio aa over time
# Read the histogram data from the text file
folder_path = "(successful) output.20230625_055402"  # TODO: Automate filename
histos = os.path.join(folder_path, "histos")
files = [file for file in os.listdir(histos) if file.endswith(".txt")]

for file in files:
    filename = os.path.join(histos, file)
    lines = None
    with open(filename, 'r') as f:
        lines = f.readlines()

    lengths = []
    reads_counts = []
    start = lines.index("Len  Reads  %Reads\n") + 1

    for line in lines[start:]:
        values = line.split()
        length = int(values[0])
        reads_count = int(values[1])
        lengths.append(length)
        reads_counts.append(reads_count)

    # Plotting the histogram without log scale
    plt.subplot(2, 1, 1)
    plt.bar(lengths, reads_counts)
    plt.xlabel('Length')
    plt.ylabel('Reads Count')
    plt.title('Read Length Histogram for ' + file)

    # Plotting the histogram with log scale
    plt.subplot(2, 1, 2)
    plt.bar(lengths, reads_counts)
    plt.xlabel('Length')
    plt.ylabel('Reads Count')
    plt.title('Read Length Histogram (Log Scale)')
    plt.yscale('log')  # Log scale the y-axis

    plt.tight_layout()  # Adjust spacing between subplots
    plt.show()

file_path = folder_path + "/log.txt"

# Lists to store the data
sample_names = []
unique_aa = []
total_aa = []

# Read the text file and extract the data
with open(file_path, 'r') as file:
    lines = file.readlines()
    for line in lines:
        if line.startswith('test'):
            line_data = line.split()
            sample_names.append(line_data[0])
            unique_aa.append(int(line_data[-3]))
            total_aa.append(int(line_data[-2]))

# Plot the data

plt.figure(figsize=(10, 6))
plt.plot(sample_names, unique_aa, marker='o', label='Unique AA')
plt.plot(sample_names, total_aa, marker='o', label='Total AA')
plt.xlabel('Round')
plt.ylabel('Count')
plt.title('Unique and Total Amino Acid Counts Throughout Experiment')
plt.legend()
plt.xticks(range(0, len(sample_names)), [f'#{i}' for i in range(1, len(sample_names) + 1)], rotation=45)
plt.tight_layout()
plt.show()