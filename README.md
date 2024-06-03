# Pymaceuticals

# Pymaceuticals Inc.
---

### Analysis

- From reviewing the box plot data it appears that Capomulin and Ramicane are the most effective; their average final tumor volumes were much lower than Infubinol and Ceftamin.
- There was a strong correlation between mouse weight and tumor volume (around 0.84) observed in the Capomulin data - meaning weight could contribute to tumor growth/size.
- There was an outlier noted in the Infubinol data.  With how far this piece of data was out of range, it was likely data error.
 
# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np
from scipy.stats import linregress

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
data_complete = pd.merge(study_results, mouse_metadata, how="left", on=["Mouse ID", "Mouse ID"])

# Display the data table for preview
data_complete.head()
# Checking the number of mice.
mouse_count = len(data_complete["Mouse ID"].unique())
mouse_count
# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mouse = data_complete.loc[data_complete.duplicated(subset=["Mouse ID", "Timepoint"]), "Mouse ID"].unique()
duplicate_mouse
# Optional: Get all the data for the duplicate mouse ID. 
dupe_mouse_info = data_complete.loc[data_complete["Mouse ID"] == "g989"]
dupe_mouse_info
# Create a clean DataFrame by dropping the duplicate mouse by its ID.

clean_data = data_complete[data_complete["Mouse ID"].isin(duplicate_mouse)==False]
clean_data.head()

# Checking the number of mice in the clean DataFrame.
clean_mouse_count = len(clean_data["Mouse ID"].unique())
clean_mouse_count
## Summary Statistics
clean_data["Drug Regimen"].unique()
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
mean = clean_data["Tumor Volume (mm3)"].groupby(clean_data["Drug Regimen"]).mean()
median = clean_data["Tumor Volume (mm3)"].groupby(clean_data["Drug Regimen"]).median()
variance = clean_data["Tumor Volume (mm3)"].groupby(clean_data["Drug Regimen"]).var()
standard_dev = clean_data["Tumor Volume (mm3)"].groupby(clean_data["Drug Regimen"]).std()
sem = clean_data["Tumor Volume (mm3)"].groupby(clean_data["Drug Regimen"]).sem()
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
stat_summary = pd.DataFrame({"Mean Tumor Volume": mean,
                        "Median Tumor Volume": median,
                        "Tumor Volume Variance": variance,
                        "Tumor Volume Std Dev": standard_dev,
                        "Tumor Volume Std Err": sem})
stat_summary
# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
clean_data.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].aggregate(["mean", "median", "var", "std", "sem"])
## Bar and Pie Charts
mice_per_drug = clean_data["Drug Regimen"].value_counts()
mice_per_drug
# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
mice_per_drug.plot.bar(color="blue")
plt.ylabel("Mice Per Drug")
plt.show()
# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
x_axis = mice_per_drug.index.values
y_axis = mice_per_drug.values
plt.bar(x_axis, y_axis, color="blue", align="edge")
plt.xlabel("Drug Regimen")
plt.ylabel("Mice per drug")
plt.xticks(rotation="vertical")

plt.show()
# Generate a pie plot showing the distribution of female versus male mice using Pandas
mice_gender.plot.pie(autopct="%1.1f%%")
plt.ylabel("Sex")
plt.show()


# Generate a pie plot showing the distribution of female versus male mice using pyplot
labels = ["Male", "Female"]
colors = ["blue", "orange"]
plt.pie(mice_gender, labels=labels, colors=colors, autopct="%1.1f%%")
plt.ylabel("Sex")
plt.show()
## Quartiles, Outliers and Boxplots
# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

drug_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
final_list = clean_data[clean_data["Drug Regimen"].isin(drug_list)]
# Start by getting the last (greatest) timepoint for each mouse

last_timepoint = final_list.groupby('Mouse ID').max()['Timepoint']


# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
last_tumor = pd.merge(data_complete, last_timepoint, on=("Mouse ID", "Timepoint"), how='right')

last_tumor
# Put treatments into a list for for loop (and later for plot labels)
# Create empty list to fill with tumor vol data (for plotting)
vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
# Locate the rows which contain mice on each drug and get the tumor volumes   
# add subset 
for drug in drug_list:
    tumor_info = last_tumor.loc[last_tumor["Drug Regimen"]==drug]["Tumor Volume (mm3)"]
    vol_data.append(tumor_info)
    
    quartiles = tumor_info.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
           
    # Determine outliers using upper and lower bounds
    outliers = tumor_info.loc[(tumor_info > upper_bound) | (tumor_info < lower_bound)]
    print(f"{drug}'s potential outliers: {outliers}")

# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.

plt.boxplot(vol_data, labels=drug_list)
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()
## Line and Scatter Plots
capomulin_data = clean_data.loc[clean_data["Drug Regimen"] == "Capomulin", :]
capomulin_data
# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
cap_mouse = capomulin_data.loc[capomulin_data["Mouse ID"] == "l509"]
cap_x = cap_mouse["Timepoint"]
cap_y = cap_mouse["Tumor Volume (mm3)"]

plt.plot(cap_x, cap_y)
plt.xlabel("Timepoint (Days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Mouse l509 treatment w/ Capomulin")
plt.show()

# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen

avg_cap = capomulin_data["Tumor Volume (mm3)"].groupby(capomulin_data["Mouse ID"]).mean()
mouse_weight = capomulin_data["Weight (g)"].groupby(capomulin_data["Mouse ID"]).mean()

plt.scatter(mouse_weight, avg_cap, color='blue')
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show
## Correlation and Regression
# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
correlation = st.pearsonr(mouse_weight, avg_cap)
print(f"The correlation between mouse weight and the average tumor volume is {round(correlation[0],2)}")

(slope, intercept, rvalue, pvalue, stderr) = linregress(mouse_weight, avg_cap)
regress_values = mouse_weight * slope + intercept
# line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(mouse_weight, avg_cap)
plt.plot(mouse_weight,regress_values,"r-")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()
