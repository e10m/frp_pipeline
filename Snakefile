# initialize directory variables
proj_dir = "/home/ethanmach1998/project/"
home_dir = "/home/ethanmach1998/"
bowtie2_dir = "/home/ethanmach1998/anaconda3/envs/bowtie2_env/bin/"
hmmc_dir = proj_dir + "hmmc_refgenomes/"


### workflow for bowtie2 SAM output ###
rule build_bowtie_indices:
    output: hmmc_dir + "{taxon_id}/{taxon_id}_index.1.bt2"
    input: hmmc_dir + "{taxon_id}/{taxon_id}_combined.fna"
    shell: bowtie2_dir + "bowtie2-build {input} " + hmmc_dir + "{wildcards.taxon_id}/{wildcards.taxon_id}_index"

rule bowtie_alignment:
    input:
        fastq = proj_dir + "samples/hmp_mock_454_even/SRR072233.fastq",
        index = hmmc_dir + "{taxon_id}/{taxon_id}_index.1.bt2"
    output: hmmc_dir + "{taxon_id}/{taxon_id}.sam"
    run:
        index_dir = f"{hmmc_dir}/{wildcards.taxon_id}/bowtie2_indices"

        # align reads and move index files
        shell(f"{bowtie2_dir}bowtie2 -x {hmmc_dir}{wildcards.taxon_id}/{wildcards.taxon_id}_index -U {input.fastq} -S {output}")
        shell(f"mkdir -p {index_dir}")
        shell(f"mv {hmmc_dir}{wildcards.taxon_id}/{wildcards.taxon_id}_index* {index_dir}")


### workflow for fr-hit SAM output ###
rule make_psl:
    input:
        reads = proj_dir + "samples/hmp_mock_454_even/SRR072233.fasta",
        combined_refgenome = hmmc_dir + "{taxon_id}/{taxon_id}_combined.fna"
    output: hmmc_dir + "{taxon_id}/{taxon_id}_frhit.psl"
    shell: "fr-hit -a {input.reads} -d {input.combined_refgenome} -c 80 -f 1 -o {output}"

rule convert_psl2sam:
    input: hmmc_dir + "{taxon_id}/{taxon_id}_frhit.psl"
    output: hmmc_dir + "{taxon_id}/{taxon_id}_frhit.sam"
    shell: "psl2sam.pl 1282_frhit.sam"

rule add_headerPlusMD:
    input:
        combined_refgenome = hmmc_dir + "{taxon_id}/{taxon_id}_combined.fna",
        sam_file = hmmc_dir + "{taxon_id}/{taxon_id}_frhit.sam"
    output: hmmc_dir + "{taxon_id}/{taxon_id}_fr_header.sam"
    shell: "samtools view -hT {input.combined_refgenome} {input.sam_file} > {output}"

rule filter_sam_reads:
    input: hmmc_dir + "{taxon_id}/{taxon_id}_frhit_header.sam"
    output: hmmc_dir + "{taxon_id}/{taxon_id}_frhit_final.sam"
    run:
        output_directory = f"{hmmc_dir}{wildcards.taxon_id}/alignment_files"

        shell(f"mkdir -p {output_directory}")
        shell(f"samtools view -h -F 4 {input} > {output}")
        shell(f"mv {hmmc_dir}{wildcards.taxon_id}/*.sam {output_directory}")
        shell(f"mv {hmmc_dir}{wildcards.taxon_id}/*.psl {output_directory}")


### Misc. rules ###
# rule for linking ref genomes together
rule combine_fragmented_genomes:
    input: hmmc_dir + "{taxon_id}{taxon_id}.fna"
    output: hmmc_dir + "{taxon_id}{taxon_id}_combined.fna"
    shell: proj_dir + "scripts/contig_combiner.py {input}"