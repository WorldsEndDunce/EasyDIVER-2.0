import matplotlib.pyplot as plt

# Read the histogram data from the text file
filename = '(successful) output.20230625_055402/histos/test1_S1_counts_histo.txt'  # TODO: Automate filename
lines = None
with open(filename, 'r') as file:
    lines = file.readlines()

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
plt.title('Read Length Histogram (Regular Scale)')

# Plotting the histogram with log scale
plt.subplot(2, 1, 2)
plt.bar(lengths, reads_counts)
plt.xlabel('Length')
plt.ylabel('Reads Count')
plt.title('Read Length Histogram (Log Scale)')
plt.yscale('log')  # Log scale the y-axis

plt.tight_layout()  # Adjust spacing between subplots
plt.show()