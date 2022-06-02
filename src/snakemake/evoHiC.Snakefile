configfile: "config/config.yml"

DATASETS = list(config["evoHiC"]["dataset_to_genome"]) + list(config["evoHiC"]["combined_dataset_to_dataset"])

wildcard_constraints:
  resolution = "\d+|rs"

#
#  Default targets
#

rule all:
  input:
    lambda wildcards:
      rules.all_multiqc.input +
      rules.all_hicexplorer.input + \
      rules.all_mcool.input + \
      rules.all_juicebox.input

rule all_multiqc:
  input:
    "data/evoHiC/qc/multiqc_report.html",
    "data/hicexplorer/qc/multiqc_evoHiC.html"

rule all_hicexplorer:
  input:
    list([
      [
        "data/hicexplorer/h5/" + dataset + "_" + resolution + ".h5",
        "data/hicexplorer/qc/" + dataset + "_" + resolution + "_corrected.png",
        "data/hicexplorer/h5/" + dataset + "_" + resolution + "_corrected.h5",
        "data/hicexplorer/h5/" + dataset + "_" + resolution + "_corrected_domains.bed"
      ] for resolution in ["5000", "rs"]
    ] for dataset in DATASETS)

rule all_mcool:
  input:
    list([
      "data/hicexplorer/cool/" + dataset + "_5000.mcool",
    ] for dataset in DATASETS)

rule all_rda:
  input:
    list([
      [
        "data/hicexplorer/rda/" + dataset + "_" + resolution + "_corrected.rda"
      ] for resolution in ["5000"]
    ] for dataset in DATASETS)

rule all_juicebox:
  input:
    list([
      "data/hic_juicebox/" + dataset + ".hic"
    ] for dataset in DATASETS)

#
#  FASTQ demultiplexing
#

rule je_demultiplex:
  input:
    fastq1 = config["evoHiC"]["fastq"]["read1"],
    fastq2 = config["evoHiC"]["fastq"]["read2"]
  output:
    ["data/evoHiC/fastq/" + library_fastq + ".fastq.gz"
      for dataset in config["evoHiC"]["dataset_to_genome"].keys()
        for library_fastq in [dataset + "_R1", dataset + "_R2"]]
  params:
    barcodes = config["evoHiC"]["barcodes"]
  conda:
    "../../env/je-suite.yaml"
  shell:
    """
    # XT=1 for trimming the T right after the barcode
    je demultiplex F1={input.fastq1} F2={input.fastq2} BF={params.barcodes} XT=1 O=data/evoHiC/fastq
    """

#
#  FASTQ quality control
#

rule fastqc:
  input:
    "data/evoHiC/fastq/{file}.fastq.gz"
  output:
    "data/evoHiC/qc/{file}_fastqc.zip"
  conda:
    "../../env/fastqc.yaml"
  shell:
    """
    fastqc -o data/evoHiC/qc {input}
    """

def library_fastq_qcfiles(wildcards):
  return ["data/evoHiC/qc/" + library_fastq + "_fastqc.zip"
    for dataset in config["evoHiC"]["dataset_to_genome"].keys()
      for library_fastq in [dataset + "_R1", dataset + "_R2"]]

rule fastqc_multiqc:
  input:
    library_fastq_qcfiles
  output:
    "data/evoHiC/qc/multiqc_report.html"
  conda:
    "../../env/multiqc.yaml"
  shell:
    """
    cd data/evoHiC/qc
    rm -rf multiqc_data/
    multiqc --interactive .
    """

#
#  Hi-C read mapping
#

def dataset_to_genome(dataset):
  if dataset in config["evoHiC"]["dataset_to_genome"]:
    return config["evoHiC"]["dataset_to_genome"][dataset]
  else: # for the combined datasets, take the genome defined for Rep1
    return config["evoHiC"]["dataset_to_genome"][dataset + "_Rep1"]

def dataset_to_fasta_genome(dataset):
  genome = dataset_to_genome(dataset)
  return config["genomes"][genome]["fasta"]

rule bwa_mem:
  input:
    "data/evoHiC/fastq/{dataset}_{read}.fastq.gz"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}_{read}.bam")
  wildcard_constraints:
    read = "R1|R2"
  params:
    fasta_genome = lambda wildcards: dataset_to_fasta_genome(wildcards.dataset)
  threads:
    16
  conda:
    "../../env/bwa_samtools.yaml"
  shell:
    """
    bwa mem -E50 -L0 -5 -t {threads} {params.fasta_genome} {input} | samtools view -bT {params.fasta_genome} - > {output}
    """

#
#  Merging Hi-C reads from separate files for read1 and read2 to single files
#

rule samtools_namesort:
  input:
    "data/evoHiC/bam/all_reads/{dataset}.bam"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}.ns.bam")
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    rm -f {output}.tmp.*
    samtools sort -n -@ {threads} -m 4G -O bam {input} -o {output}
    """

rule add_flag1:
  input:
    "data/evoHiC/bam/all_reads/{dataset}_R1.ns.bam"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}_R1.ns.addFlag.bam")
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/sh/addFlag.py {input} 65 > {output}
    """

rule add_flag2:
  input:
    "data/evoHiC/bam/all_reads/{dataset}_R2.ns.bam"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}_R2.ns.addFlag.bam")
  conda:
    "../../env/pysam.yaml"
  shell:
    """
    src/sh/addFlag.py {input} 129 > {output}
    """

rule samtools_merge:
  input:
    "data/evoHiC/bam/all_reads/{dataset}_R1.ns.addFlag.bam",
    "data/evoHiC/bam/all_reads/{dataset}_R2.ns.addFlag.bam"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}.merge.bam")
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools merge -n -@ {threads} {output} {input}
    """

rule samtools_fixmate:
  input:
    "data/evoHiC/bam/all_reads/{dataset}.merge.bam"
  output:
    temp("data/evoHiC/bam/all_reads/{dataset}.merge.fixmate.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools fixmate -p {input} {output}
    """

#
#  Sort reads and index the resulting .BAM files
#

rule samtools_sort:
  input:
    "data/evoHiC/bam/all_reads/{dataset}.bam"
  output:
    "data/evoHiC/bam/all_reads/{dataset}.sort.bam"
  threads:
    4
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    rm -f {output}.tmp.*
    samtools sort -@ {threads} -m 4G -O bam {input} -o {output}
    """

rule samtools_index:
  input:
    "data/evoHiC/bam/all_reads/{dataset}.bam"
  output:
    "data/evoHiC/bam/all_reads/{dataset}.bam.bai"
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools index {input}
    """

#
#  Hi-C read processing using HiCExplorer
#

rule hicFindRestSite:
  output:
    "analysis/rest_site_positions_" + config["evoHiC"]["restrictionEnzyme"] + "_{genome}.bed"
  params:
    fasta_genome = lambda wildcards: config["genomes"][wildcards.genome]["fasta"],
    restrictionSequence = '"' + config["evoHiC"]["restrictionSequence"].replace(" ", "|") + '"'
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicFindRestSite --fasta {params.fasta_genome} --searchPattern {params.restrictionSequence} -o {output}
    """

rule hicBuildMatrix:
  input:
    bam1 = "data/evoHiC/bam/all_reads/{dataset}_R1.bam",
    bam2 = "data/evoHiC/bam/all_reads/{dataset}_R2.bam",
    bed = lambda wildcards: "analysis/rest_site_positions_" + config["evoHiC"]["restrictionEnzyme"] + "_" + dataset_to_genome(wildcards.dataset) + ".bed"
  output:
    bam = "data/hicexplorer/bam/filtered_reads/{dataset}_{resolution}.bam",
    h5 = "data/hicexplorer/h5/{dataset}_{resolution}.h5",
    qc = "data/hicexplorer/qc/{dataset}_{resolution}/hicQC.html"
  wildcard_constraints:
    resolution = "\d+"
  params:
    restrictionSequence = config["evoHiC"]["restrictionSequence"],
    danglingSequence = config["evoHiC"]["danglingSequence"]
  threads:
    8
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicBuildMatrix --samFiles {input.bam1} {input.bam2} --outBam {output.bam} --outFileName {output.h5} \
      --restrictionSequence {params.restrictionSequence} \
      --danglingSequence {params.danglingSequence} \
      --restrictionCutFile {input.bed} \
      --binSize {wildcards.resolution} \
      --QCfolder data/hicexplorer/qc/{wildcards.dataset}_{wildcards.resolution} \
      --threads {threads}
    """

rule hicBuildMatrix_rs:
  input:
    bam1 = "data/evoHiC/bam/all_reads/{dataset}_R1.bam",
    bam2 = "data/evoHiC/bam/all_reads/{dataset}_R2.bam",
    bed = lambda wildcards: "analysis/rest_site_positions_" + config["evoHiC"]["restrictionEnzyme"] + "_" + dataset_to_genome(wildcards.dataset) + ".bed"
  output:
    bam = "data/hicexplorer/bam/filtered_reads/{dataset}_rs.bam",
    h5 = "data/hicexplorer/h5/{dataset}_rs.h5",
    qc = "data/hicexplorer/qc/{dataset}_rs/hicQC.html"
  conda:
    "../../env/hicexplorer.yaml"
  params:
    restrictionSequence = config["evoHiC"]["restrictionSequence"],
    danglingSequence = config["evoHiC"]["danglingSequence"]
  threads:
    8
  shell:
    """
    hicBuildMatrix --samFiles {input.bam1} {input.bam2} --outBam {output.bam} --outFileName {output.h5} \
      --restrictionSequence {params.restrictionSequence} \
      --danglingSequence {params.danglingSequence} \
      --restrictionCutFile {input.bed} \
      --QCfolder data/hicexplorer/qc/{wildcards.dataset}_rs \
      --threads {threads}
    """

def combined_dataset_to_dataset_h5(wildcards):
  if wildcards.dataset in config["evoHiC"]["combined_dataset_to_dataset"]:
    datasets = config["evoHiC"]["combined_dataset_to_dataset"][wildcards.dataset]
    return ["data/hicexplorer/h5/" + dataset + "_" + wildcards.resolution + ".h5" for dataset in datasets]
  else:
    return "/missing_dataset"

rule hicSumMatrices:
  input:
    combined_dataset_to_dataset_h5
  output:
    "data/hicexplorer/h5/{dataset}_{resolution}.h5"
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicSumMatrices -m {input} -o {output}
    """

def dataset_to_chromosomes(dataset):
  genome = dataset_to_genome(dataset)
  return config["genomes"][genome]["chromosomes"]

rule hicCorrectMatrix:
  input:
    "data/hicexplorer/h5/{dataset}_{resolution}.h5"
  output:
    png = "data/hicexplorer/qc/{dataset}_{resolution}_corrected.png",
    h5 = "data/hicexplorer/h5/{dataset}_{resolution}_corrected.h5"
  params:
    chromosomes = lambda wildcards: dataset_to_chromosomes(wildcards.dataset),
    filterThreshold = config["evoHiC"]["filterThreshold"]
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    # note that the first invocation of hicCorrectMatrix modifies the input file(!)
    hicCorrectMatrix diagnostic_plot --chromosomes {params.chromosomes} -m {input} -o {output.png}
    hicCorrectMatrix correct --chromosomes {params.chromosomes} -m {input} --filterThreshold {params.filterThreshold} -o {output.h5}
    """

# when combining replicates: first sum raw count matrices, then normalize the resulting matrix
ruleorder:
  hicCorrectMatrix > hicSumMatrices

rule hicConvertFormat:
  input:
    "data/hicexplorer/h5/{dataset}.h5"
  output:
    "data/hicexplorer/tsv/{dataset}.tsv.gz"
  params:
    output_prefix = "data/hicexplorer/tsv/{dataset}"
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicConvertFormat --matrices {input} --outFileName {params.output_prefix} --inputFormat h5 --outputFormat ginteractions
    gzip {params.output_prefix}.tsv
    """

rule hicFindTADs:
  input:
    "data/hicexplorer/h5/{dataset}_corrected.h5"
  output:
    "data/hicexplorer/h5/{dataset}_corrected_domains.bed"
  threads:
    16
  conda:
    "../../env/hicexplorer.yaml"
  shell:
    """
    hicFindTADs -m {input} --outPrefix data/hicexplorer/h5/{wildcards.dataset}_corrected \
      --correctForMultipleTesting=fdr --numberOfProcessors {threads}
    """

rule hicexplorer_multiqc:
  input:
    ["data/hicexplorer/qc/" + dataset + "_rs/hicQC.html" for dataset in config["evoHiC"]["dataset_to_genome"]]
  output:
    "data/hicexplorer/qc/multiqc_evoHiC.html"
  params:
    directories = [dataset + "_rs" for dataset in config["evoHiC"]["dataset_to_genome"]]
  conda:
    "../../env/multiqc.yaml"
  shell:
    """
    cd data/hicexplorer/qc
    rm -rf multiqc_evoHiC_data/
    multiqc --interactive -n multiqc_evoHiC {params.directories}
    """

#
#  Convert HiCExplorer output to HiGlass .mcool files
#

rule cooler_zoomify_1000:
  input:
    "data/hicexplorer/cool/{dataset}_1000.cool"
  output:
    "data/hicexplorer/cool/{dataset}_1000.mcool"
  conda:
    "../../env/cooler.yaml"
  shell:
    """
    cooler zoomify --resolutions 1000N --balance {input}
    """

rule cooler_zoomify_5000:
  input:
    "data/hicexplorer/cool/{dataset}_5000.cool"
  output:
    "data/hicexplorer/cool/{dataset}_5000.mcool"
  conda:
    "../../env/cooler.yaml"
  shell:
    """
    cooler zoomify --resolutions 5000,10000N --balance {input}
    """

#
#  Convert HiCExplorer output to Juicebox hic files
#

def combined_dataset_to_dataset_bam(wildcards):
  if wildcards.dataset in config["evoHiC"]["combined_dataset_to_dataset"]:
    datasets = config["evoHiC"]["combined_dataset_to_dataset"][wildcards.dataset]
    return ["data/hicexplorer/bam/filtered_reads/" + dataset + "_rs.bam" for dataset in datasets]
  else:
    return "/missing_dataset"

rule samtools_cat_for_juicebox:
  input:
    combined_dataset_to_dataset_bam
  output:
    temp("data/hicexplorer/bam/filtered_reads/{dataset}_rs.bam")
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools cat {input} -o {output}
    """

rule bamtobed_for_juicebox:
  input:
    "data/hicexplorer/bam/filtered_reads/{dataset}_rs.bam"
  output:
    temp("data/hic_juicebox/{dataset}.txt.gz")
  threads:
    8
  conda:
    "../../env/bedtools_pigz.yaml"
  shell:
    """
    bedtools bamtobed -i {input} -bedpe \
      | awk '{{ OFS="\\t"; print $9 != "+", "chr" $1, ($9 == "+" ? $2 + 1 : $3), "0", $10 != "+", "chr" $4, ($10 == "+" ? $5 + 1 : $6), "1" }}' \
      | pigz -c -p {threads} > {output}
    """

rule juicebox_sort:
  input:
    "data/hic_juicebox/{dataset}.txt.gz"
  output:
    temp("data/hic_juicebox/{dataset}.txt.sorted.gz")
  threads:
    8
  conda:
    "../../env/bedtools_pigz.yaml"
  shell:
    """
    pigz -dc -p {threads} {input} \
      | sort -k7,7n \
      | sort --stable -k3,3n \
      | sort --stable -k6,6 \
      | sort --stable -k2,2 \
      | pigz -c -p {threads} > {output}
    """

def dataset_to_chrom_sizes_for_juicebox(dataset):
  genome = dataset_to_genome(dataset)
  return config["genomes"][genome]["chrom_sizes_for_juicebox"]

rule juicebox_pre:
  input:
    "data/hic_juicebox/{dataset}.txt.sorted.gz"
  output:
    "data/hic_juicebox/{dataset}.hic"
  params:
    juicebox = config["software"]["juicebox"],
    chrom_sizes = lambda wildcards: dataset_to_chrom_sizes_for_juicebox(wildcards.dataset)
  shell:
    """
    {params.juicebox} pre {input} {output} {params.chrom_sizes}
    """
