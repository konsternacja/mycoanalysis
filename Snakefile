configfile: "configs/config.yaml"


(
    SEASON, 
    PART, 
    SAMPLE, 
    _, 
    PAIR
) = glob_wildcards(config['input_dir'] + 
                   'season{season}_{part}_fastp/{sample}_fastp/{sample2}_{pair}_fastp.fq.gz')


onsuccess:
    print("Workflow finished successfully!")


onerror:
    print("Error has occured. Please, check log files for more details.")
    

rule ALL:
    input:
        config['output_dir'] + "kraken2/all_reports.txt",
        config['output_dir'] + "bracken/all_reports.txt",
        config['output_dir'] + "bracken/bracken_species_all",
        config['output_dir'] + "funguild/funguild_guilds.taxa.guilds.txt",
        expand(config['output_dir'] + "statistics/season{season}_{part}/{sample}.txt", zip, part=PART, sample=SAMPLE, season=SEASON),
        config['output_dir'] + "statistics/beta_diversity.txt",


# KRAKEN2 #


rule KRAKEN2:
    input:
        expand(
            config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kraken", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        expand(
            config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kreport", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        expand(
            config['output_dir'] + "kraken2/season{season}_{part}/{sample}.biom", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        config['output_dir'] + "kraken2/all_reports.txt",
    output:
        config['temp_dir'] + "finished_Kraken2"
    threads: 1
    message:
        "[KRAKEN2] Finished Kraken2 analysis"
    shell:
        """
        touch {output}
        """


rule run_kraken2:
    """
    Performes taxonomic classification of reads using Kraken2.
    """
    input:
        forward_reads = config['input_dir'] + "season{season}_{part}_fastp/{sample}_fastp/{sample}_R1_fastp.fq.gz",
        reverse_reads = config['input_dir'] + "season{season}_{part}_fastp/{sample}_fastp/{sample}_R2_fastp.fq.gz",
    output:
        kraken_report = config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kreport",
        kraken_output = config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kraken",
    threads: 8
    params:
        db = config['database_dir']
    shell:
        '''
        kraken2 --use-names --db {params.db} --threads 8 --report {output.kraken_report} --paired --gzip-compressed {input.forward_reads} {input.reverse_reads} > {output.kraken_output}
        '''


rule kraken2biom:
    """
    Converts each Kraken2 report to BIOM format.
    """
    input:
        expand(
            config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kreport", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
    output:
        biom_report = config['output_dir'] + "kraken2/season{season}_{part}/{sample}.biom",
    shell:
        '''
        kraken-biom {input} -o {output.biom_report}
        '''


rule merge_kraken2_reports:
    """
    Merges all Kraken2 reports into one file.
    """
    input:
        expand(
            config['output_dir'] + "kraken2/season{season}_{part}/{sample}.biom", season=SEASON, part=PART, sample=SAMPLE
            ),
    output:
        txt_report = config['output_dir'] + "kraken2/all_reports.txt",
    shell:
        '''
        biom convert -i {input} -o {output.txt_report} --to-tsv --header-key taxonomy
        '''


# BRACKEN #
        

rule BRACKEN:
    input:
        expand(
            config['output_dir'] + "bracken/season{season}_{part}/{sample}.bracken", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        expand(
            config['output_dir'] + "bracken/season{season}_{part}/{sample}_bracken_species.kreport", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        expand(
            config['output_dir'] + "bracken/season{season}_{part}/{sample}.biom", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
        config['output_dir'] + "bracken/all_reports.txt",
        config['output_dir'] + "bracken/bracken_species_all",
    output:
        config['temp_dir'] + "finished_bracken"
    threads: 1
    message:
        "[BRACKEN] Finished Bracken analysis"
    shell:
        """
        touch {output}
        """


rule run_bracken:
    """
    Estimates abundance of species using Bracken.
    """
    input: 
        config['output_dir'] + "kraken2/season{season}_{part}/{sample}.kreport",
    output:
        bracken = config['output_dir'] + "bracken/season{season}_{part}/{sample}.bracken",
        bracken_species = config['output_dir'] + "bracken/season{season}_{part}/{sample}_bracken_species.kreport",
    params:
        kmer_distrib = config['kmer_distrib']
    shell:
        '''
        python scripts/handle_bracken.py {input} {params.kmer_distrib} {output.bracken} {output.bracken_species}
        if [ ! -f {output.bracken} ]; then touch {output.bracken}; fi
        if [ ! -f {output.bracken_species} ]; then touch {output.bracken_species}; fi
        '''


rule merge_bracken:
    """
    Merges all Bracken reports into one file using a bracken script
    """
    input:
        expand(config['output_dir'] + "bracken/season{season}_{part}/{sample}.bracken", zip, part=PART, sample=SAMPLE, season=SEASON),
    output:
        config['output_dir'] + "bracken/bracken_species_all",
    shell:
        '''
        python scripts/bracken_combine_outputs.py --files {input} -o {output}
        '''


rule bracken2biom:
    """
    Converts each Bracken report to BIOM format.
    """
    input:
        expand(
            config['output_dir'] + "bracken/season{season}_{part}/{sample}_bracken_species.kreport", zip, part=PART, sample=SAMPLE, season=SEASON
            ),
    output:
        config['output_dir'] + "bracken/season{season}_{part}/{sample}.biom",
    shell:
        '''
        kraken-biom {input} -o {output}
        '''


rule merge_bracken_reports:
    """
    Merges all Bracken reports into one file.
    """
    input:
        expand(
            config['output_dir'] + "bracken/season{season}_{part}/{sample}.biom", season=SEASON, part=PART, sample=SAMPLE
            ),
    output:
        config['output_dir'] + "bracken/all_reports.txt",
    shell:
        '''
        biom convert -i {input} -o {output} --to-tsv --header-key taxonomy
        '''


# FUNGUILD #


rule FUNGUILD:
    input:
        config['output_dir'] + "funguild/funguild.txt",
        config['output_dir'] + "funguild/funguild.taxa.txt",
        config['output_dir'] + "funguild/funguild_guilds.taxa.guilds.txt",
    output:
        config['temp_dir'] + "finished_funguild"
    threads: 1
    message:
        "[FUNGUILD] Finished FUNGuild analysis"
    shell:
        """
        touch {output}
        """


rule change_table_funguild:
    input:
        config['output_dir'] + "bracken/all_reports.txt",
    output:
        config['output_dir'] + "funguild/funguild.txt",
    shell:
        '''
        python scripts/change_table_funguild.py {input} {output}
        '''


rule funguild_taxa:
    input:
        config['output_dir'] + "funguild/funguild.txt",
    output:
        config['output_dir'] + "funguild/funguild.taxa.txt",
    shell:
        '''
        python scripts/FUNGuild/FUNGuild.py taxa -otu {input} -format tsv -column taxonomy -classifier unite > {output}
        '''


rule funguild_guilds:
    input:
        config['output_dir'] + "funguild/funguild.taxa.txt",
    output:
        config['output_dir'] + "funguild/funguild_guilds.taxa.guilds.txt",
    shell:
        '''
        python scripts/FUNGuild/FUNGuild.py guild -taxa {input} > {output}
        '''


# STATISTICS #


rule alpha_diversity:
    input:
        config['output_dir'] + "bracken/season{season}_{part}/{sample}.bracken"
    output:
        config['output_dir'] + "statistics/season{season}_{part}/{sample}.txt",
    shell:
        '''
        python scripts/KrakenTools/alpha_diversity.py -f {input} -a Sh > {output}
        '''

rule beta_diversity:
    input:
        expand(config['output_dir'] + "bracken/season{season}_{part}/{sample}.bracken", zip, part=PART, sample=SAMPLE, season=SEASON),
    output:
        config['output_dir'] + "statistics/beta_diversity.txt",
    shell:
        '''
        python scripts/KrakenTools/beta_diversity.py -i {input} --type bracken > {output}
