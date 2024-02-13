configfile: "configs/config.yaml"


SEASON, PART, SAMPLE, _, PAIR = glob_wildcards(
    config["input_dir"]
    + "season{season}_{part}_fastp/{sample}_fastp/{sample2}_{pair}_fastp.fq.gz"
)


rule all:
    input:
        # expand(
        # config["output_dir"] + "kraken2/season{season}_{part}/{sample}.kreport",
        # zip,
        # part=PART,
        # sample=SAMPLE,
        # season=SEASON,
        # ),
        # expand(
        # config["output_dir"] + "kraken2/season{season}_{part}/{sample}.kraken",
        # zip,
        # part=PART,
        # sample=SAMPLE,
        # season=SEASON,
        # ),
        # expand(
        # config["output_dir"] + "bracken/season{season}_{part}/{sample}.bracken",
        # zip,
        # part=PART,
        # sample=SAMPLE,
        # season=SEASON,
        # ),
        # expand(
        # config["output_dir"]
        # + "bracken/season{season}_{part}/{sample}_bracken_species.kreport",
        # zip,
        # part=PART,
        # sample=SAMPLE,
        # season=SEASON,
        # ),
        expand(
            config["output_dir"] + "statistics/season{season}_{part}/{sample}.txt",
            zip,
            part=PART,
            sample=SAMPLE,
            season=SEASON,
        ),
        # config["output_dir"] + "kraken2/all_reports.biom",
        # config["output_dir"] + "kraken2/all_reports.txt",
        # config["output_dir"] + "bracken/bracken_species_all",
        # config["output_dir"] + "bracken/all_reports.biom",
        # config["output_dir"] + "bracken/all_reports.txt",
        # config["output_dir"] + "funguild/funguild.taxa.guilds.txt",
        config["output_dir"] + "statistics/beta_diversity.txt",


rule run_kraken:
    input:
        forward_reads=config["input_dir"]
        + "season{season}_{part}_fastp/{sample}_fastp/{sample}_R1_fastp.fq.gz",
        reverse_reads=config["input_dir"]
        + "season{season}_{part}_fastp/{sample}_fastp/{sample}_R2_fastp.fq.gz",
    output:
        kraken_report=config["output_dir"]
        + "kraken2/season{season}_{part}/{sample}.kreport",
        kraken_output=config["output_dir"]
        + "kraken2/season{season}_{part}/{sample}.kraken",
    params:
        db=config["database_dir"],
    shell:
        """
        kraken2 --use-names --db {params.db} --threads 8 --report {output.kraken_report} --paired --gzip-compressed {input.forward_reads} {input.reverse_reads} > {output.kraken_output}
        """


rule kraken2biom:
    input:
        expand(
            config["output_dir"] + "kraken2/season{season}_{part}/{sample}.kraken",
            zip,
            part=PART,
            sample=SAMPLE,
            season=SEASON,
        ),
    output:
        biom_reports=config["output_dir"] + "kraken2/all_reports.biom",
        txt_reports=config["output_dir"] + "kraken2/all_reports.txt",
    conda:
        "configs/biom.yaml"
    shell:
        """
        kraken-biom {input} -o {output.biom_reports}
        biom convert -i {output.biom_reports} -o {output.txt_reports} --to-tsv --header-key taxonomy
        """


rule run_bracken:
    input:
        config["output_dir"] + "kraken2/season{season}_{part}/{sample}.kreport",
    output:
        bracken=config["output_dir"] + "bracken/season{season}_{part}/{sample}.bracken",
        bracken_species=config["output_dir"]
        + "bracken/season{season}_{part}/{sample}_bracken_species.kreport",
    params:
        kmer_distrib=config["kmer_distrib"],
    shell:
        """
        python scripts/handle_bracken.py {input} {params.kmer_distrib} {output.bracken} {output.bracken_species}
        if [ ! -f {output.bracken} ]; then touch {output.bracken}; fi
        if [ ! -f {output.bracken_species} ]; then touch {output.bracken_species}; fi
        """


rule merge_bracken:
    input:
        expand(
            config["output_dir"] + "bracken/season{season}_{part}/{sample}.bracken",
            zip,
            part=PART,
            sample=SAMPLE,
            season=SEASON,
        ),
    output:
        config["output_dir"] + "bracken/bracken_species_all",
    shell:
        """
        python scripts/bracken_combine_outputs.py --files {input} -o {output}
        """


rule bracken2biom:
    input:
        expand(
            config["output_dir"]
            + "bracken/season{season}_{part}/{sample}_bracken_species.kreport",
            zip,
            part=PART,
            sample=SAMPLE,
            season=SEASON,
        ),
    output:
        biom_reports=config["output_dir"] + "bracken/all_reports.biom",
        txt_reports=config["output_dir"] + "bracken/all_reports.txt",
    conda:
        "configs/biom.yaml"
    shell:
        """
        kraken-biom {input} -o {output.biom_reports}
        biom convert -i {output.biom_reports} -o {output.txt_reports} --to-tsv --header-key taxonomy
        """


rule run_funguild:
    input:
        config["output_dir"] + "bracken/all_reports.txt",
    output:
        modified_otu_table=config["output_dir"] + "funguild/funguild.txt",
        taxa=temp(config["output_dir"] + "funguild/funguild.taxa.txt"),
        guilds=config["output_dir"] + "funguild/funguild_guilds.taxa.guilds.txt",
    shell:
        """
            python scripts/change_table_funguild.py {input} {output.modified_otu_table}
            python scripts/FUNGuild/FUNGuild.py taxa -otu {output.modified_otu_table} -format tsv -column taxonomy -classifier unite
            python scripts/FUNGuild/FUNGuild.py guild -taxa {output.taxa} 
        """


rule alpha_diversity:
    input:
        config["output_dir"] + "bracken/season{season}_{part}/{sample}.bracken",
    output:
        config["output_dir"] + "statistics/season{season}_{part}/{sample}.txt",
    shell:
        """
        python scripts/KrakenTools/alpha_diversity.py -f {input} -a Sh > {output}
        """


rule beta_diversity:
    input:
        expand(
            config["output_dir"] + "bracken/season{season}_{part}/{sample}.bracken",
            zip,
            part=PART,
            sample=SAMPLE,
            season=SEASON,
        ),
    output:
        config["output_dir"] + "statistics/beta_diversity.txt",
    shell:
        """
        python scripts/KrakenTools/beta_diversity.py -i {input} --type bracken > {output}
        """
