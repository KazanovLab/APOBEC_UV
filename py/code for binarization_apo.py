import pandas as pd
import os

input_dir = "/mnt/c/users/bkurt/Desktop/apobec_data/contexts"
output_dir = "/mnt/c/users/bkurt/Desktop/apobec_data/newfiles"

os.makedirs(output_dir, exist_ok=True)

# Function to classify APOBEC mutations
def is_apobec(row):
    ref, alt, context = row["REF"], row["ALT"], row["3_nuc_context"]
    
    if len(context) != 3:
        return 0  # Invalid context, not APOBEC
    
    # APOBEC mutation patterns
    if context[:2] == "TT" and ref == "C" and alt == "T":  
        return 1  # TC → TT

    if context[:2] == "TG" and ref == "C" and alt == "G":  
        return 1  # TC → TG 

    #complementary strand mirors of mutations
    if context[1:] == "AA" and ref == "G" and alt == "A":
        return 1 # GA → AA 
    if context[1:] == "CA" and ref == "G" and alt == "C":
        return 1 # GA → CA
    
    return 0  # Not UV-induced

# Loop through all 12 sample files
for i in range(1, 13):
    file_name = f"processed_Sample{i}.txt"
    input_path = os.path.join(input_dir, file_name)
    
    if os.path.isfile(input_path):
        print(f"Processing {file_name}...")

        data = pd.read_csv(input_path, sep="\t", low_memory=False)

        data["APOBEC"] = data.apply(is_apobec, axis=1)

        clustered_data = data[data["CLUSTER_TYPE"] != "NONE"].copy()

        clustered_data["Group_Number"] = clustered_data["CLUSTER"].str.extract(r"(groupNumber\d+)")

        # Calculate APOBEC fraction per cluster
        cluster_apobec_fraction = (
            clustered_data.groupby("Group_Number")["UV-Induced"]
            .agg(["sum", "count"])
            .reset_index()
        )
        cluster_apobec_fraction["APOBEC_Fraction"] = cluster_apobec_fraction["sum"] / cluster_apobec_fraction["count"]

        data = data.merge(cluster_apobec_fraction[["Group_Number", "APOBEC_Fraction"]],
                          how="left",
                          left_on="CLUSTER",
                          right_on="Group_Number")

        data.drop(columns=["Group_Number"], inplace=True)

        output_path = os.path.join(output_dir, f"sample{i}_with_apobec_and_fraction.txt")
        data.to_csv(output_path, sep="\t", index=False)

        print(f"✅ Sample {i}: Processed and APOBEC fractions calculated and saved.")
    else:
        print(f"⚠️ File {file_name} not found. Skipping...")

print("🎉 All files processed successfully!")
