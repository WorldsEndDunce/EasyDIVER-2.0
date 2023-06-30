import os
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(style="darkgrid")

# TODO: total aa, unique aa, unique bio aa over time
# Read the histogram data from the text file
folder_path = "large_output"  # TODO: Automate filename?
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
    plt.figure(figsize=(12, 8))
    # Plotting the histogram without log scale
    plt.subplot(2, 1, 1)
    plt.bar(lengths, reads_counts, color="steelblue")
    plt.xlabel('Length')
    plt.ylabel('Reads Count')
    plt.title('Read Length Histogram for ' + file)

    # Plotting the histogram with log scale
    plt.subplot(2, 1, 2)
    plt.bar(lengths, reads_counts, color="steelblue")
    plt.xlabel('Length')
    plt.ylabel('Reads Count')
    plt.title('Read Length Histogram (Log Scale)')
    plt.yscale('log')  # Log scale the y-axis

    plt.tight_layout()  # Adjust spacing between subplots
    # plt.show()
    plt.savefig("figures/" + file[:file.rfind(".")] + ".png", dpi=500)
    plt.close()

file_path = folder_path + "/log.txt"

# Lists to store the data
sample_names = []
unique_aa = {}
total_aa = {}
max_round = -1
# Read the text file and extract the data
with open(file_path, 'r') as file:
    lines = file.readlines()
    start = False
    for line in lines:
        if line.startswith("sample"):
            start = True
            continue
        elif start:
            line_data = line.split()
            sample_names.append(line_data[0])
            max_round = max(max_round, int(line_data[0].split("-")[0]))
            unique_aa.setdefault(line_data[0].split("-")[1], []).append(int(line_data[-2])) # TODO: Don't hardcode?
            total_aa.setdefault(line_data[0].split("-")[1], []).append(int(line_data[-1]))

print("Unique amino acid (AA) counts throughout rounds: " + str(unique_aa))
print("Total amino acid (AA) counts throughout rounds: " + str(total_aa))
# Plot the data
plt.figure(figsize=(10, 6))

for key in unique_aa.keys():
    color = "orange"
    if key == "neg":
        color = "green"
    elif key == "in":
        color = "red"
    label_u = "Unique AA " + key
    label_t = "Total AA " + key
    plt.plot(unique_aa[key], marker="o", label=label_u, color=color)
    plt.plot(total_aa[key], marker="o", label=label_t, color=color)

plt.xlabel('Round')
plt.ylabel('Count')
plt.title('Unique and Total Amino Acid Counts Throughout Experiment')
plt.legend()
plt.xticks(range(0, max_round), [f'#{i}' for i in range(1, max_round + 1)], rotation=45)
plt.tight_layout()
# plt.show()
plt.savefig("figures/aa.png", dpi=500)