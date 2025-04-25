import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# Load the data
data_path = "E:/Docking_Score.csv"
df = pd.read_csv(data_path)
df.columns.values[[8, 40]] = ['1,4-Chrysenequinone', 'Photofrin']
# Assuming the first column contains protein names and the rest are drug names
protein_names = df.iloc[:, 0].values
drug_names = df.columns[1:].values[:50]  # Select the first 50 drugs
binding_scores = df.iloc[:, 1:51].values  # Select scores corresponding to the first 50 drugs

# Get the min and max scores for color bar limits
min_score = np.min(binding_scores)
max_score = np.max(binding_scores)

# Create a figure and axis, set the aspect to 'equal' for square cells
plt.rcParams["font.family"] = "Arial"
plt.figure(figsize=(15, 15))  # Set figure size with equal height and width for square cells

# Define the custom colormap: dark maroon -> light maroon -> light green -> deep green
cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", 
    [
        (0.8, 0.1, 0.3),  # Deep Red
        (1, 0.3, 0.5),           # Red
        (1, 0.5, 0.7),
        (1, 0.7, 0.9),# Light Red
        (0.8, 1, 1),         # Light Cyan
        (0, 1, 1),           # Cyan
        (0, 0.85, 0.85)      # Light Deep Cyan
    ], N=256)

# Plot the matrix with square aspect
im = plt.imshow(binding_scores, cmap=cmap, aspect='equal', vmin=min_score, vmax=max_score)

# Set the title and labels
plt.title('', fontsize=16)
plt.xlabel('', fontsize=14)
plt.ylabel('', fontsize=14)

# Set drug names on the x-axis with rotation
xticks=plt.xticks(np.arange(len(drug_names)), drug_names, rotation=45, ha='right', fontsize=12,fontweight='bold')
plt.setp(xticks[1][:10], color='darkgreen')
plt.yticks(np.arange(len(protein_names)), protein_names, fontsize=12,fontweight='bold')  # Protein names on y-axis


# Add color bar, reduce its size, and adjust height to match plot height
cbar = plt.colorbar(im, fraction=0.03, pad=0.02, shrink=0.17, aspect=9)
cbar.set_label('', fontsize=12,fontweight='bold')

# Set ticks dynamically for the color bar
cbar.set_ticks(np.linspace(min_score, max_score, num=8))
cbar.ax.set_yticklabels([f"{tick:.2f}" for tick in np.linspace(min_score,-5.40, num=8)],fontweight='bold')  # Customize tick labels

plt.tight_layout()  # Adjust layout to prevent clipping of tick labels
# Adjust layout to maximize the height of the plot within the A4 page
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9, hspace=0.2)

# Save the plot to PDF in the Downloads folder
downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
output_pdf = os.path.join(downloads_folder, "binding_affinity_matrix.pdf")
plt.savefig(output_pdf, format='pdf',dpi=600)

plt.show()

print(f"Plot saved as PDF at: {output_pdf}")


####################################### Validation With Independent #################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# Load the data
data_path = "E:/Docking_Score_Independent.csv"
df = pd.read_csv(data_path)
df.columns.values[16] = '1,4-Chrysenequinone'
# Assuming the first column contains protein names and the rest are drug names
protein_names = df.iloc[:, 0].values
drug_names = df.columns[1:].values[:50]  # Select the first 50 drugs
binding_scores = df.iloc[:, 1:51].values  # Select scores corresponding to the first 50 drugs

# Get the min and max scores for color bar limits
min_score = np.min(binding_scores)
max_score = np.max(binding_scores)

# Create a figure and axis, set the aspect to 'equal' for square cells
plt.rcParams["font.family"] = "Arial"
plt.figure(figsize=(15, 15))  # Set figure size with equal height and width for square cells

# Define the custom colormap: dark maroon -> light maroon -> light green -> deep green
cmap = mcolors.LinearSegmentedColormap.from_list("custom_colormap", 
    [
        (0.8, 0.1, 0.3),  # Deep Red
        (1, 0.3, 0.5),           # Red
        (1, 0.5, 0.7),
        (1, 0.7, 0.9),# Light Red
        (0.8, 1, 1),         # Light Cyan
        (0, 1, 1),           # Cyan
        (0, 0.85, 0.85)      # Light Deep Cyan
    ], N=256)

# Plot the matrix with square aspect
im = plt.imshow(binding_scores, cmap=cmap, aspect='equal', vmin=min_score, vmax=max_score)

# Set the title and labels
plt.title('', fontsize=16)
plt.xlabel('', fontsize=14)
plt.ylabel('', fontsize=14)

# Set drug names on the x-axis with rotation
xticks=plt.xticks(np.arange(len(drug_names)), drug_names, rotation=45, ha='right', fontsize=12,fontweight='bold')
plt.setp([xticks[1][i] for i in [0, 1, 2, 5, 8, 9]], color='darkgreen')
#plt.setp([xticks[1][i] for i in [0, 1, 2, 5, 8, 9,11,12,15,19]], color='darkgreen')
plt.yticks(np.arange(len(protein_names)), protein_names, fontsize=12,fontweight='bold')  # Protein names on y-axis


# Add color bar, reduce its size, and adjust height to match plot height
cbar = plt.colorbar(im, fraction=0.011, pad=0.02, shrink=0.17, aspect=9)
cbar.set_label('', fontsize=12,fontweight='bold')

# Set ticks dynamically for the color bar
cbar.set_ticks(np.linspace(min_score, max_score, num=8))
cbar.ax.set_yticklabels([f"{tick:.2f}" for tick in np.linspace(min_score, max_score, num=8)],fontweight='bold')  # Customize tick labels

plt.tight_layout()  # Adjust layout to prevent clipping of tick labels
# Adjust layout to maximize the height of the plot within the A4 page
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9, hspace=0.2)

# Save the plot to PDF in the Downloads folder
downloads_folder = os.path.join(os.path.expanduser("~"), "Downloads")
output_pdf = os.path.join(downloads_folder, "binding_affinity_matrix_independent2.pdf")
plt.savefig(output_pdf, format='pdf',dpi=600)

plt.show()

print(f"Plot saved as PDF at: {output_pdf}")