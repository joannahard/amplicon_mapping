SAMPLES, = glob_wildcards("data/input/{smp}.2.fastq.gz")



rule all:
    input:
        expand("data/output/{smp}/mapped/{smp}.mapped.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam.bai", smp=SAMPLES),
        expand("data/output/{smp}/mapped/{smp}_ac.txt", smp=SAMPLES),

rule bwa:
    input:
        r1 = "data/input/{smp}.1.fastq.gz",
        r2 = "data/input/{smp}.2.fastq.gz",
        ref = config["ref"],
        index = config["ref"] + ".bwt"
    output:
        temp("data/output/{smp}/mapped/{smp}.mapped.bwa.bam")
    threads: 16
    params:
        bwa = "-M",
        java = config["javaopts"]
    log:
        bwa = "data/output/{smp}/logs/{smp}.bwa.log",
        sambam = "data/output/{smp}/logs/{smp}.picard_samtobam.log"
    shell:
        "bwa mem {params.bwa} -t {threads} {input.ref} {input.r1} {input.r2} 2> {log.bwa} > {output} | "
        "picard {params.java} SamFormatConverter INPUT=/dev/stdin OUTPUT={output} 2> {log.sambam}"



rule filter_and_fix:
    input:
        "data/output/{smp}/mapped/{smp}.mapped.bwa.bam"
    output:
        temp("data/output/{smp}/mapped/{smp}.filtered.bwa.bam")
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


rule ac:
    input:
        bam = "data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam",
        index = "data/output/{smp}/mapped/{smp}.filtered.sorted.bwa.bam.bai",
        ref = config["ref"],
        loci = config["loci"]
    output:
        "data/output/{smp}/mapped/{smp}_ac.txt"
    shell:
        "alleleCounter -l {input.loci} -r {input.ref} -b {input.bam} -o {output}"


