import os
import matplotlib.pyplot as plt

# Read the histogram data from the text file
folder_path = '(successful) output.20230625_055402/histos'  # Replace with the actual folder path
files = [file for file in os.listdir(folder_path) if file.endswith('.txt')]

for file in files:
    filename = os.path.join(folder_path, file)  # TODO: Automate filename
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