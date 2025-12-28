import subprocess
import os
import sys
import multiprocessing  # <--- NEW: To detect your CPU cores

# --- CONFIGURATION ---
genome_fasta = "./wigconvert/SRR1693438/sequence.fasta"
raw_fastq = "./wigconvert/SRR1693438/SRR1693438.fastq"       # <--- CONFIRM THIS MATCHES YOUR FILE NAME
trimmed_fastq = "trimmed.fastq"
output_prefix = "biased_profile"
chromosome_name = "NC_000913.2"

def run_command(cmd):
    print(f"RUNNING: {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"ERROR: Command failed: {cmd}")
        sys.exit(1)

def convert_depth_to_wig(depth_file, wig_file, chrom_name):
    # (This function remains the same as before)
    print(f"Formatting {wig_file}...")
    header = f"track type=wiggle_0 name={wig_file} viewLimits=-5:5\n"
    with open(depth_file, 'r') as f_in, open(wig_file, 'w') as f_out:
        f_out.write(header)
        current_step_start = -1
        last_pos = -1
        buffer = []
        for line in f_in:
            parts = line.strip().split()
            if len(parts) < 3: continue
            pos = int(parts[1])
            count = float(parts[2])
            if pos != last_pos + 1:
                if buffer:
                    f_out.write(f"fixedStep chrom={chrom_name} start={current_step_start} step=1\n")
                    for val in buffer: f_out.write(f"{val}\n")
                    buffer = []
                current_step_start = pos
            buffer.append(count)
            last_pos = pos
        if buffer:
            f_out.write(f"fixedStep chrom={chrom_name} start={current_step_start} step=1\n")
            for val in buffer: f_out.write(f"{val}\n")

def main():
    # Detect Cores (Leave 1 free for the OS just in case, minimum 1)
    cores = max(1, multiprocessing.cpu_count() - 1)
    print(f"--- DETECTED {multiprocessing.cpu_count()} CORES. USING {cores} FOR SPEED. ---")

    # 0. Check Tools
    try:
        subprocess.run("cutadapt --version", shell=True, check=True, stdout=subprocess.DEVNULL)
    except:
        print("ERROR: 'cutadapt' is missing! Run: sudo apt install cutadapt")
        sys.exit(1)

    # 1. Build Index
    if not os.path.exists(f"{genome_fasta}.1.ebwt"):
        print("--- Step 1: Building Genome Index ---")
        run_command(f"bowtie-build {genome_fasta} {genome_fasta}")
    
    # 1.5. TRIM ADAPTERS (Accelerated)
    print("--- Step 1.5: Trimming Adapters (Fast Mode) ---")
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    # -j {cores}: This is the speed upgrade
    run_command(f"cutadapt -j {cores} -a {adapter} -m 20 --trim-n -o {trimmed_fastq} {raw_fastq}")

    # 2. Align Reads (Accelerated)
    print("--- Step 2: Aligning Reads (Fast Mode) ---")
    sam_file = "aligned.sam"
    # -p {cores}: This speeds up alignment
    run_command(f"bowtie -p {cores} -v 1 -m 1 --best --strata {genome_fasta} -q {trimmed_fastq} -S {sam_file}")

    # 3. Sort
    print("--- Step 3: Sorting ---")
    bam_file = "aligned.sorted.bam"
    # Samtools sort also supports threads with -@
    run_command(f"samtools view -uS {sam_file} | samtools sort -@ {cores} - -o {bam_file}")
    run_command(f"samtools index {bam_file}")

    # 4. Split Strands & Coverage
    print("--- Step 4 & 5: Processing Strands ---")
    run_command(f"samtools view -b -F 16 {bam_file} > fwd.bam")
    run_command(f"samtools view -b -f 16 {bam_file} > rev.bam")
    run_command("samtools depth -a fwd.bam > fwd.depth")
    run_command("samtools depth -a rev.bam > rev.depth")

    # 6. WIG
    print("--- Step 6: Creating WIG Files ---")
    convert_depth_to_wig("fwd.depth", f"{output_prefix}_f.wig", chromosome_name)
    convert_depth_to_wig("rev.depth", f"{output_prefix}_r.wig", chromosome_name)

    print("\nSUCCESS! Files are ready for Riboformer.")

if __name__ == "__main__":
    main()