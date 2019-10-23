rule porechop:
    input:
        "01_basecalling/{run}/pass"
    output:
        directory("01_basecalling/{run}/pass_demultiplexed")
    log:
        "01_basecalling/{run}/MeBaPiNa_demultiplexed.log" ##  > {log} 2>&1 (TLDR both to log) redirects stdout to file and stderr to stdout which is redirected to a file
    benchmark:
        "01_basecalling/{run}/pass_demultiplexed/MeBaPiNa.benchmark.tsv"
    conda:
        "../envs/porechop.yml"
    threads:
        config["machine"]["cpu"]
    params:
        ("--verbosity 2"), ## or nothing to log
    shell:
        "porechop --threads {threads} {params} --input {input} --barcode_dir {output} > {log} 2>&1" 
