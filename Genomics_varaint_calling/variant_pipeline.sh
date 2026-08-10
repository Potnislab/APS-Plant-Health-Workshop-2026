#!/bin/bash
set -e # Exit immediately if any command fails

# Configuration (These allow your Python website to pass files to the script)
REF="${1:-Xperforans1991.fasta}"
SAMPLE_DIR="${2:-.}"
OUT_DIR="${3:-./results}"
THREADS="${4:-4}"

echo "=================================================="
echo "STARTING VARIANT PIPELINE"
echo "=================================================="

# STEP 0: Load Modules exactly as required by the protocol
echo "Loading HPC modules..."
module load bwa 2>/dev/null || echo "BWA loaded"
module load samtools/1.2 2>/dev/null || module load samtools
module load bcftools 2>/dev/null || echo "BCFTools loaded"
module load gatk 2>/dev/null || echo "GATK loaded"
module load vcftools 2>/dev/null || echo "VCFTools loaded"

# STEP 0.5: Download Picard if it doesn't exist
if [ ! -d "picard-tools-1.119" ]; then
    echo "Downloading Picard tools..."
    wget http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip
    unzip -q picard-tools-1.119.zip
fi

# STEP 1: Indexing
echo "STEP 1: Indexing Reference Genome..."
if [ ! -f "${REF}.bwt" ]; then
    bwa index "$REF"
    samtools faidx "$REF"
    java -jar ./picard-tools-1.119/CreateSequenceDictionary.jar R="$REF" O="${REF%.fasta}.dict"
fi

# Find all Read 1 FASTQ files and process them
for read1 in "$SAMPLE_DIR"/*_1.f*q*; do
    [ -e "$read1" ] || continue 
    
    read2="${read1/_1/_2}"
    sample=$(basename "${read1%_1.*}")
    sample_out="$OUT_DIR/$sample"
    mkdir -p "$sample_out"
    
    echo "=================================================="
    echo "Processing Sample: $sample"
    echo "=================================================="

    # File names based on your protocol
    sam_file="$sample_out/${sample}_mapped_bwa.sam"
    bam_file="$sample_out/${sample}_mapped_bwa.bam"
    sorted_bam="$sample_out/${sample}_mapped_bwa_sort.bam"
    marked_bam="$sample_out/${sample}_mapped_bwa_sort_marked.bam"
    metrics_file="$sample_out/${sample}_markduplicates_metrics.txt"
    
    raw_vcf="$sample_out/${sample}_test.vcf"
    gatk_vcf="$sample_out/${sample}_test.filter_gatk.vcf"
    final_snp_prefix="$sample_out/${sample}_test_filter_onlysnp"

    # STEP 2: Mapping
    if [ ! -f "$sam_file" ]; then
        echo "STEP 2: Mapping reads with BWA..."
        bwa mem -M -t "$THREADS" "$REF" "$read1" "$read2" > "$sam_file"
    fi

    # STEP 3: SAM to BAM Conversion & Sorting
    if [ ! -f "$sorted_bam" ]; then
        echo "STEP 3: Converting to BAM and sorting..."
        samtools view -Sb "$sam_file" > "$bam_file"
        samtools sort -O bam -o "$sorted_bam" -T "${sample_out}/temp" "$bam_file"
        rm "$bam_file" "$sam_file" # Clean up massive intermediate files
    fi

    # STEP 4: Removing Duplicates (Picard)
    if [ ! -f "$marked_bam" ]; then
        echo "STEP 4: Removing PCR duplicates with Picard..."
        java -jar ./picard-tools-1.119/MarkDuplicates.jar \
            REMOVE_DUPLICATES=true \
            INPUT="$sorted_bam" \
            OUTPUT="$marked_bam" \
            METRICS_FILE="$metrics_file" \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT
            
        samtools index "$marked_bam"
    fi

    # STEP 5: Variant Calling
    # Note: Using bcftools mpileup instead of samtools mpileup to prevent your previous Version Error
    if [ ! -f "$raw_vcf" ]; then
        echo "STEP 5: Calling variants..."
        bcftools mpileup -Ou -f "$REF" "$marked_bam" | \
        bcftools call -O v -v -c -o "$raw_vcf"
    fi

    # STEP 6.1: Variant Filtering (GATK)
    if [ ! -f "$gatk_vcf" ]; then
        echo "STEP 6.1: Filtering variants with GATK..."
        gatk VariantFiltration \
            -R "$REF" \
            -O "$gatk_vcf" \
            --variant "$raw_vcf" \
            --filter-expression "FQ < 0.025 && MQ > 50 && QUAL > 100 && DP > 15" \
            --filter-name pass_filter
            
        pass_count=$(grep -c 'pass_filter' "$gatk_vcf" || echo "0")
        echo "Variants passing GATK filter: $pass_count"
    fi

    # STEP 6.2: Keeping only SNPs (VCFTools)
    if [ ! -f "${final_snp_prefix}.recode.vcf" ]; then
        echo "STEP 6.2: Removing INDELs with VCFTools..."
        vcftools \
            --vcf "$gatk_vcf" \
            --keep-filtered pass_filter \
            --remove-indels \
            --recode \
            --recode-INFO-all \
            --out "$final_snp_prefix"
    fi

    echo "Finished processing $sample"
done

echo "=================================================="
echo "PIPELINE COMPLETE"
echo "=================================================="
