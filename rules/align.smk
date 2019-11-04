rule minimap2_index:
    input:
        "00_rawdata/reference_sequences/lambda/{reference}.fasta"
    output:
        "00_rawdata/reference_sequences/lambda/{reference}.mmi"
    log:
        "00_rawdata/reference_sequences/lambda/{reference}_MeBaPiNa_minimap2_index.log"
    benchmark:
        "00_rawdata/reference_sequences/lambda/{reference}_MeBaPiNa_minimap2_index.benchmark.tsv"
    params:
        "-x map-ont" ## naopore specific
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -d {output} {input} > {log} 2>&1"

def input_minimap2(wildcards):
    from os import listdir, path
    basecall_dir = "01_processeddata/" + wildcards.run + "/basecall/pass"
    barc_dirs = [path.join(basecall_dir, s) for s in listdir(basecall_dir)]
    return barc_dirs

rule minimap2:
    input:
        input_minimap2,
        dummy_dependency="01_processeddata/{run}/basecall/sequencing_summary.txt", ## used as dummy for the other folders
        target=expand("00_rawdata/reference_sequences/lambda/{reference}.mmi", reference=config["align"]["reference"])
    output:
        "01_processeddata/{run}/align/alignment.sam"
    log:
        "01_processeddata/{run}/align/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align/MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find {input} -type f -name \"*.fastq.gz\") "
        "> {log} 2>&1"

rule minimap2_calibration_strands:
    input:
        dummy_dependency="01_processeddata/{run}/basecall/sequencing_summary.txt", ## used as dummy for the other folders
        target="00_rawdata/reference_sequences/lambda/lambda_3.6kb.mmi"
    output:
        "01_processeddata/{run}/align_calibration_strands/alignment.sam"
    log:
        "01_processeddata/{run}/align_calibration_strands/MeBaPiNa_minimap2.log"
    benchmark:
        "01_processeddata/{run}/align_calibration_strands/MeBaPiNa_minimap2.benchmark.tsv"
    params:
        "-x map-ont", ## naopore specific
        "-a" ## possition accurate CIGAR alignment in SAM output; much slower <- maybe skip?
    conda:
        "../envs/minimap2.yml"
    threads:
        config["machine"]["cpu"]
    shell:
        "minimap2 -t {threads} {params} -o {output} "
        "{input.target} "
        "$(find 01_processeddata/{wildcards.run}/basecall/calibration_strands -type f -name \"*.fastq.gz\") "
        "> {log} 2>&1"
