# initialize directory variables
proj_dir = "/home/ethanmach1998/project"
home_dir = "/home/ethanmach1998"
bowtie2_dir = "/home/ethanmach1998/anaconda3/envs/bowtie2_env/bin"
hmmc_dir = proj_dir + "/hmmc_refgenomes"


### workflow for bowtie2 SAM output ###
rule build_bowtie_indices:
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_index.1.bt2"
    input:
        hmmc_dir + "/{taxon_id}/{taxon_id}_combined.fna"
    shell:
        bowtie2_dir + "/bowtie2-build {input} " + hmmc_dir + "/{wildcards.taxon_id}/{wildcards.taxon_id}_index"

rule bowtie_alignment:
    input:
        fastq = proj_dir + "/samples/hmp_mock_454_even/SRR072233.fastq",
        index = hmmc_dir + "/{taxon_id}/{taxon_id}_index.1.bt2"
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bowtie2.sam"
    shell:
        bowtie2_dir + "/bowtie2 -x " + hmmc_dir + "/{wildcards.taxon_id}/{wildcards.taxon_id}_index -U {input.fastq} -S {output}"

rule filter_bt2_sam:
    input:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bowtie2.sam"
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bowtie2_final.sam"
    shell:
        "samtools view -h -F 4 {input} > {output}"


### workflow for BWA SAM output ###
rule build_bwa_indices:
    input:
        hmmc_dir + "/{taxon_id}/{taxon_id}_combined.fna"
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_combined.fna.bwt"
    shell:
        "bwa index {input}"

rule bwa_alignment:
    input:
        ref_genome = hmmc_dir + "/{taxon_id}/{taxon_id}_combined.fna",
        reads = proj_dir + "/samples/hmp_mock_454_even/SRR072233.fastq"
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bwa.sam"
    shell:
        "bwa mem {input.ref_genome} {input.reads} > {output}"

rule filter_bwa_sam_reads:
    input:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bwa.sam"
    output:
        hmmc_dir + "/{taxon_id}/{taxon_id}_bwa_final.sam"
    shell:
        "samtools view -h -F 4 {input} > {output}"


### Reference Genome Manipulation Rules ###
rule combine_fragmented_genomes:
    input: hmmc_dir + "/{taxon_id}/{taxon_id}.fna"
    output: hmmc_dir + "/{taxon_id}/{taxon_id}_combined.fna"
    shell: proj_dir + "/scripts/contig_combiner.py {input}"


### Plotting Rule ###
rule plot_sam_data:
    input:
        bt_sam = hmmc_dir + "/{taxon_id}/{taxon_id}_bowtie2_final.sam",
        bwa_sam = hmmc_dir + "/{taxon_id}/{taxon_id}_bwa_final.sam"
    output:
        bt_frp = hmmc_dir + "/{taxon_id}/{taxon_id}_bowtie2_frp.png",
        bwa_frp = hmmc_dir + "/{taxon_id}/{taxon_id}_bwa_frp.png"
    run:
        shell(proj_dir + "/scripts/sam_to_frp.py {input.bt_sam}")
        shell(proj_dir + "/scripts/sam_to_frp.py {input.bwa_sam}")