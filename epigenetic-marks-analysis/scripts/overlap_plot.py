import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv('overlaps.bed', sep="\t")

# Assuming that 'Chromosome' column specifies the chromosome number
# and that you want to visualize the overlaps on a single chromosome
# Filter data for a specific chromosome if needed
df = df[df['Chromosome'] == 'chr1']
# Select the first 10 rows
df = df.iloc[0:9]

# Sort the data by the Start position
df_sorted = df.sort_values(by='Start')

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))

# A list to keep track of the position where we place the text
last_end = 0
for index, row in df_sorted.iterrows():
    # Draw a line for each feature with a rectangle
    ax.plot([row['Start'], row['End']], [index, index], color='brown')
    # Place text in the middle of the rectangle
    mid_point = (row['Start'] + row['End']) / 2
    if mid_point > last_end:  # To avoid text overlap
        ax.text(mid_point, index, row['Tissue'], verticalalignment='center', horizontalalignment='center')
        last_end = mid_point + (last_end - row['Start']) / 2  # Adjust the position for the next text


# Draw a line for functional sequence
ax.plot([df_sorted.iloc[0]['StartMatch'], df_sorted.iloc[0]['EndMatch']], [10, 10], color='blue')
# Place text in the middle of the rectangle
mid_point = (df_sorted.iloc[0]['StartMatch'] + df_sorted.iloc[0]['EndMatch']) / 2

ax.text(mid_point, 10, df_sorted.iloc[0]['ID'], verticalalignment='center', horizontalalignment='center')
last_end = mid_point + (last_end - df_sorted.iloc[0]['StartMatch']) / 2  # Adjust the position for the next text

# Set the limits and labels
ax.set_ylim(-1, len(df_sorted))
ax.set_xlim(df_sorted['Start'].min(), df_sorted['End'].max())
ax.set_xlabel('Genomic Position')
ax.set_yticks([])
ax.set_title('Overlap Graph')

# Show grid
ax.grid(True)

# Show the plot
plt.show()
