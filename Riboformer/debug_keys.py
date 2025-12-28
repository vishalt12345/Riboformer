import os

# CONFIG: Update this path to where your data is
DATA_DIR = "datasets/GSE63858_Validation" 
FASTA_FILE = "sequence.fasta"
WIG_FILE = "biased_profile_f.wig"

def check_keys():
    print("--- DIAGNOSTIC REPORT ---")
    
    # 1. CHECK FASTA KEYS
    fasta_path = os.path.join(DATA_DIR, FASTA_FILE)
    fasta_keys = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Mimic standard parsing: Take everything after > until the first space
                    key = line.strip().split()[0][1:] 
                    fasta_keys.append(f"'{key}'")
                    # Also show what it would look like if we didn't split (common bug source)
                    full_header = line.strip()[1:]
                    print(f"FASTA Raw Header: >{full_header}")
    except FileNotFoundError:
        print("ERROR: FASTA file not found!")
    
    print(f"FASTA Keys found: {fasta_keys}")
    
    # 2. CHECK WIG KEYS
    wig_path = os.path.join(DATA_DIR, WIG_FILE)
    wig_keys = []
    try:
        with open(wig_path, 'r') as f:
            for line in f:
                if "chrom=" in line:
                    # Extract whatever is after chrom=
                    parts = line.split()
                    for p in parts:
                        if p.startswith("chrom="):
                            wig_keys.append(f"'{p.split('=')[1]}'")
                            break
            # Only show unique keys to keep it clean
            wig_keys = list(set(wig_keys))
    except FileNotFoundError:
        print("ERROR: WIG file not found!")
        
    print(f"WIG Keys found:   {wig_keys}")
    
    # 3. VERDICT
    target = "'NC_000913.2'"
    if target in fasta_keys and target in wig_keys:
        print("\nSUCCESS: The keys match! The error might be in how data_processing.py imports them.")
    else:
        print("\nFAILURE: Mismatch detected.")
        if target not in fasta_keys:
            print("-> The FASTA file is missing the key (or it is parsed differently).")
        if target not in wig_keys:
            print("-> The WIG file is missing the key.")

if __name__ == "__main__":
    check_keys()