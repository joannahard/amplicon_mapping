SAMPLES, = glob_wildcards("data/input/{smp}_r1_paired.fastq.gz")



rule all:
    input:
        expand("data/output/{smp}/mapped/{smp}.mapped.bwa.sam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.mapped.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam.bai", smp=SAMPLES),
        expand("data/output/{smp}/fastqc/", smp=SAMPLES)
        expand("data/output/{smp}/logs/{smp}.qualimap/qualimapReport.html", smp=SAMPLES)        

rule bwa:
    input:
        r1 = "data/input/{smp}_r1_paired.fastq.gz",
        r2 = "data/input/{smp}_r2_paired.fastq.gz",
        ref = config["ref"],
        index = config["ref"] + ".bwt"
    output:
        "data/output/{smp}/mapped/{smp}.mapped.bwa.sam"
    threads: 16
    params:
        "-M"
    log:
        "data/output/{smp}/logs/{smp}.bwa.log",
    shell:
        "bwa mem {params} -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} > {output}"



rule sam_to_bam:
    input:
        "data/output/{smp}/mapped/{smp}.mapped.bwa.sam"
    output:
        "data/output/{smp}/mapped/{smp}.mapped.bwa.bam"
    params:
        java = config["javaopts"]
    log:
        "data/output/{smp}/logs/{smp}.picard.sam2bam.log"
    shell:
        "picard {params.java} SamFormatConverter INPUT={input} OUTPUT={output} > {log} 2>&1;"


rule filter_and_fix:
    input:
        "data/output/{smp}/mapped/{smp}.mapped.bwa.bam"
    output:
        "data/output/{smp}/mapped/{smp}.filtered.bwa.bam"
    params:
        filters = "-b -q 2 -F 8",
        sort = "SORT_ORDER=coordinate",
        read_groups = "CREATE_INDEX=true RGID={smp} RGLB={smp} RGPL=ILLUMINA RGSM={smp} RGCN=\"NA\" RGPU=\"NA\"",
        java = config["javaopts"]
    log:
        filters = "data/output/{smp}/logs/{smp}.samtools.filters.log",
        sort = "data/output/{smp}/logs/{smp}.picard.sortsam.log",
        read_groups = "data/output/{smp}/logs/{smp}.picard.addorreplacereadgroup.log"
    shell:
        "samtools view {params.filters} {input} 2> {log.filters} |"
        "picard {params.java} SortSam {params.sort} INPUT=/dev/stdin OUTPUT=/dev/stdout 2> {log.sort} |"
        "picard {params.java} AddOrReplaceReadGroups {params.read_groups} INPUT=/dev/stdin OUTPUT={output} 2> {log.read_groups}"



rule samtools_sort:
    input:
        "data/output/{smp}/mapped/{smp}.filtered.bwa.bam"
    output:
        "data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index:
    input:
        "data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam"
    output:
        "data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam.bai"
    shell:
        "samtools index {input}"




rule qualimap:
    input: "data/output/{smp}/mapped/{smp}.filtered.bwa.bam"
    output:
      report = "data/output/{smp}/logs/{smp}.qualimap/qualimapReport.html",
      gr = "data/output/{smp}/logs/{smp}.qualimap/genome_results.txt",
      ish = "data/output/{smp}/logs/{smp}.qualimap/raw_data_qualimapReport/insert_size_histogram.txt",
      ch = "data/output/{smp}/logs/{smp}.qualimap/raw_data_qualimapReport/coverage_histogram.txt",
      gc = "data/output/{smp}/logs/{smp}.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"
    log: "data/output/{smp}/logs/{smp}.qualimap/qualimap.log"
    params: "-sd -sdmode 0 --java-mem-size=20G -c -nw 400 -gd hg19"
    threads: 10
    shell: "qualimap bamqc -nt {threads} {params} -bam {input} -outdir $(dirname {output.report}) > {log} 2>&1;"



rule fastqc:
    input:
        r1 = "data/input/{smp}_r1_paired.fastq.gz",
        r2 = "data/input/{smp}_r2_paired.fastq.gz"
    output:
        "data/output/{smp}/fastqc/"
    shell:
        "fastqc --quiet --outdir {output} --extract  -f fastq {input.r1} {input.r2}"


