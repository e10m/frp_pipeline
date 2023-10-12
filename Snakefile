# initialize working directory variable
working_dir = "/home/ethanmach1998/project/"

# workflow for bowtie2 SAM output
rule build_bowtie_indices:
    output: working_dir + "hmmc_refgenomes/{taxon_id}/{taxon_id}_index.1.bt2"
    input: working_dir + "hmmc_refgenomes/{taxon_id}/{taxon_id}.fna"
    shell: "/home/ethanmach1998/anaconda3/envs/bowtie2_env/bin/bowtie2-build {input} /home/ethanmach1998/project/hmmc_refgenomes/{wildcards.taxon_id}/{wildcards.taxon_id}_index"

rule bowtie_alignment:
    input:
        fastq = working_dir + "samples/hmp_mock_454_even/SRR072233.fastq",
        index = working_dir + "hmmc_refgenomes/{taxon_id}/{taxon_id}_index.1.bt2"
    output: working_dir + "hmmc_refgenomes/{taxon_id}/{taxon_id}.sam"
    shell: "/home/ethanmach1998/anaconda3/envs/bowtie2_env/bin/bowtie2 \
    -x /home/ethanmach1998/project/hmmc_refgenomes/{wildcards.taxon_id}/{wildcards.taxon_id}_index \
    -U {input.fastq} \
    -S {output}"


### workflow for fr-hit SAM output ###
rule make_psl:
    input: working_dir + "samples/hmp_mock_454_even/SRR072233.fasta"
    output: working_dir + "hmmc_refgenomes/{taxon_id}/{taxon_id}_frhit.psl"
    shell: "/usr/local/bin/fr-hit/fr-hit -a {input} \
    -d " + working_dir + "hmmc_refgenomes/{wildcards.taxon_id}/{wildcards.taxon_id}.fna \
    -c 80 \
    -f 1 \
    -o " + working_dir + "hmmc_refgenomes/{wildcards.taxon_id}/{wildcards.taxon_id}_frhit.psl"

# rule psl_to_sam:
