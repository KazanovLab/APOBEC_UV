import pandas as pd
import os

input_dir = "/mnt/c/users/bkurt/Desktop/Sabancı Marat Kazanov APOBEC/contfiles"
output_dir = "/mnt/c/users/bkurt/Desktop/uv/newfiles"

os.makedirs(output_dir, exist_ok=True)

# Function to classify UV-induced mutations
def is_uv_induced(row):
    ref, alt, context = row["REF"], row["ALT"], row["3_nuc_context"]
    
    if len(context) != 3:
        return 0  # Invalid context, not UV
    
    # UV mutation patterns
    if context[:2] == "TT" and ref == "C" and alt == "T":  
        return 1  # TC → TT

    if context[:2] == "CT" and ref == "C" and alt == "T":  
        return 1  # CC → CT (mutation in second C)

    if context[1:] == "TC" and ref == "C" and alt == "T":  
        return 1  # CC → TC (mutation in first C)

    if context == "TAT" and ref == "T" and alt == "A":  
        return 1  # TTT → TAT (T -> A in TTT)

    if context == "GCT" and ref == "T" and alt == "C":  
        return 1  # GTT → GCT (T -> C in GTT)

    #complementary strand mirors of mutations
    if context[1:] == "AA" and ref == "G" and alt == "A":
        return 1 # GA → AA 
    if context[1:] == "AG" and ref == "G" and alt == "A":
        return 1 # GG → AG
    if context[:2] == "GA" and ref == "G" and alt == "A":
        return 1 # GG → GA
    if context == "ATA" and ref == "A" and alt == "T":
        return 1 # AAA → ATA
    if context == "AGC" and ref == "A" and alt == "G":
        return 1 # AAC → AGC
    
    return 0  # Not UV-induced

# Loop through all 15 sample files
for i in range(1, 16):
    file_name = f"processed_Sample{i}.txt"
    input_path = os.path.join(input_dir, file_name)
    
    if os.path.isfile(input_path):
        print(f"Processing {file_name}...")

        data = pd.read_csv(input_path, sep="\t", low_memory=False)

        data["UV-Induced"] = data.apply(is_uv_induced, axis=1)

        clustered_data = data[data["CLUSTER_TYPE"] != "NONE"].copy()

        clustered_data["Group_Number"] = clustered_data["CLUSTER"].str.extract(r"(groupNumber\d+)")

        # Calculate UV fraction per cluster
        cluster_uv_fraction = (
            clustered_data.groupby("Group_Number")["UV-Induced"]
            .agg(["sum", "count"])
            .reset_index()
        )
        cluster_uv_fraction["UV_Fraction"] = cluster_uv_fraction["sum"] / cluster_uv_fraction["count"]

        data = data.merge(cluster_uv_fraction[["Group_Number", "UV_Fraction"]],
                          how="left",
                          left_on="CLUSTER",
                          right_on="Group_Number")

        data.drop(columns=["Group_Number"], inplace=True)

        output_path = os.path.join(output_dir, f"sample{i}_with_uv_and_fraction.txt")
        data.to_csv(output_path, sep="\t", index=False)

        print(f"✅ Sample {i}: Processed and UV fractions calculated and saved.")
    else:
        print(f"⚠️ File {file_name} not found. Skipping...")

print("🎉 All files processed successfully!")
