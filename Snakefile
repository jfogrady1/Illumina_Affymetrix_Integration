OUT_DIR = config["output_dir"]


wildcard_constraints:
    breed = "[A-Za-z]+",      # letters only
    chr   = "\d+"             # digits only

rule all:
    input:
        expand(OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected_CHR{chr}.vcf.gz', chr=range(1,30)),
        expand(OUT_DIR + 'work/data_formatting/{breed}.vcf.gz', breed=['Chillingham','WIDDE']),
        OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected.vcf.gz',
        OUT_DIR + 'work/ARS/WIDDE_Chill_Kerry_merged_intersected_ARS1.2.vcf.gz'

# You should have plink on Rodeo
rule convert_chill_bed:
    input:
        raw_chill = config["Chillingham_raw"]["ped_file"]
    output:
        binaries= expand(OUT_DIR + "work/data_formatting/Chillingham.{ext}", ext=['bim','fam','bed'] )
    params:
        out = lambda wc,output:output.binaries[0][:-4],
        chill_prefix = lambda wc, input: input.raw_chill[:-4]
    envmodules:
        'stack/2024-06',
        'gcc/12.2.0',
        'plink/1.9-beta6.27'
    shell:
        """
    plink --cow \
        --file {params.chill_prefix} \
        --keep-allele-order \
        --make-bed \
        --autosome \
        --set-missing-var-ids @_# \
        --double-id \
        --out {params.out} \
        --threads {threads}
        """

rule convert_WIDDE_bed:
    input:
        raw_widde = config["WIDDE_raw"]["ped_file"]
    output:
        binaries= expand(OUT_DIR + "work/data_formatting/WIDDE.{ext}", ext=['bim','fam','bed'] )
    params:
        out = lambda wc,output:output.binaries[0][:-4],
        widde_prefix = lambda wc, input: input.raw_widde[:-4]
    envmodules:
        'stack/2024-06',
        'gcc/12.2.0',
        'plink/1.9-beta6.27'
    shell:
        """
    plink --cow \
        --file {params.widde_prefix} \
        --keep-allele-order \
        --maf 0.05 \
        --geno 0.05 \
        --autosome \
        --make-bed \
        --set-missing-var-ids @_# \
        --double-id \
        --out {params.out} \
        --threads {threads}
        """

rule convert_to_vcf:
    input:
        bed = OUT_DIR + "work/data_formatting/{breed}.bed"
    output:
        vcf = OUT_DIR + 'work/data_formatting/{breed}.vcf.gz',
        index = OUT_DIR + 'work/data_formatting/{breed}.vcf.gz.tbi'
    params:
        in_file = lambda wc, input: input.bed[:-4],
        out = lambda wc,output: output.vcf[:-7]
    conda:
        'bcftools'
    envmodules:
        'stack/2024-06',
        'gcc/12.2.0',
        'plink/1.9-beta6.27',
    shell:
        """
    plink --cow \
        --bfile {params.in_file} \
        --keep-allele-order \
        --recode vcf bgz \
        --out {params.out} \
        --threads {threads}
    bcftools index --tbi --threads {threads} {output.vcf}
    """

rule fix_chill_AB_encoding:
    input:
        vcf = OUT_DIR + 'work/data_formatting/Chillingham.vcf.gz',
        array_annotation = config["illumina_annotation"],
        script = OUT_DIR + 'scripts/Chillingham_AB_coding.R'
    output:
        vcf = OUT_DIR + 'work/data_formatting/Chillingham_updated.vcf.gz'

    envmodules:
        'stack/2024-06',
        'gcc/12.2.0',
        'r/4.5.1',
        'bcftools/1.22'
    threads:
        4
    params:
        vcf = lambda wc, output: output.vcf[:-3]
    shell:
        '''
        Rscript {input.script} {input.array_annotation} {input.vcf} {params.vcf}
        bcftools view {params.vcf} -Oz -o {output.vcf}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''


rule fix_Kerry_AB_encoding:
    input:
        vcf = config['KY_vcf_raw'],
        array_annotation = config["affymetrix_annotation"],
        script = OUT_DIR + 'scripts/Kerry_AB_coding.R'
    output:
        vcf = OUT_DIR + 'work/data_formatting/Kerry_updated.vcf.gz'

    envmodules:
        'stack/2024-06',
        'gcc/12.2.0',
        'r/4.5.1',
        'bcftools/1.22'
    threads:
        4
    params:
        vcf = lambda wc, output: output.vcf[:-3]
    shell:
        '''
        Rscript {input.script} {input.array_annotation} {input.vcf} {params.vcf}
        bcftools view {params.vcf} -Ou | bcftools norm --rm-dup snps -Oz -o {output.vcf}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''

rule allele_check_split_target:
    input:
        vcf = OUT_DIR + 'work/data_formatting/Chillingham_updated.vcf.gz'
    output:
        sep_vcf = OUT_DIR + 'work/data_formatting/Chillingham_updated_CHR{chr}.vcf.gz'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    params:
        chromosome = lambda wc: wc.chr
    resources:
        mem_mb_per_cpu = 2000
    threads:
        10
    shell:
        '''
        bcftools view \
        --regions {params.chromosome} \
        --output-type z \
        --threads {threads} \
        {input.vcf} \
        -o {output.sep_vcf}

        bcftools index --tbi --threads {threads} {output.sep_vcf}
        '''

rule split_Kerry:
    input:
        vcf = OUT_DIR + 'work/data_formatting/Kerry_updated.vcf.gz'
    output:
        sep_vcf = OUT_DIR + 'work/data_formatting/Kerry_updated_CHR{chr}.vcf.gz'
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    params:
        chromosome = lambda wc: wc.chr
    resources:
        mem_mb_per_cpu = 2000
    threads:
        10
    shell:
        '''
        bcftools view \
        --regions {params.chromosome} \
        --output-type z \
        --threads {threads} \
        {input.vcf} \
        -o {output.sep_vcf}

        bcftools index --tbi --threads {threads} {output.sep_vcf}
        '''

rule split_WIDDE:
    input:
        vcf = OUT_DIR + 'work/data_formatting/WIDDE.vcf.gz'
    output:
        sep_vcf = OUT_DIR + 'work/data_formatting/WIDDE_CHR{chr}.vcf.gz'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    params:
        chromosome = lambda wc: wc.chr
    resources:
        mem_mb_per_cpu = 2000
    threads:
        10
    shell:
        '''
        bcftools view \
        --regions {params.chromosome} \
        --output-type z \
        --threads {threads} \
        {input.vcf} \
        -o {output.sep_vcf}
        bcftools index --tbi --threads {threads} {output.sep_vcf}
        '''

# Conform WIDDE as we know the correct orientation of
# Chillingham cattle because we set it above
rule conform_WIDDE:
    input:
        vcf_in = rules.allele_check_split_target.output,
        conform_gt = config['tools']['conform_gt'],
        widde = rules.split_WIDDE.output.sep_vcf,
    output:
        vcf = OUT_DIR + 'work/data_formatting/WIDDE_conformed_CHR{chr}.vcf.gz',
        index = OUT_DIR + 'work/data_formatting/WIDDE_conformed_CHR{chr}.vcf.gz.tbi'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        'openjdk/21.0.3_9',
        'bcftools/1.22'
    params:
        chromosome = lambda wc: wc.chr,
        prefix=lambda wc,output:output.vcf[:-7]
    resources:
        mem_mb_per_cpu = 4000
    threads:
        5
    shell:
        '''
        java -jar {input.conform_gt} \
        ref={input.vcf_in} \
        gt={input.widde} \
        chrom={params.chromosome} \
        match=POS \
        strict=false \
        out={params.prefix}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''

rule conform_Kerry:
    input:
        vcf_in = rules.allele_check_split_target.output,
        conform_gt = config['tools']['conform_gt'],
        widde = rules.split_Kerry.output.sep_vcf,
    output:
        vcf = OUT_DIR + 'work/data_formatting/Kerry_conformed_CHR{chr}.vcf.gz',
        index = OUT_DIR + 'work/data_formatting/Kerry_conformed_CHR{chr}.vcf.gz.tbi'
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        'openjdk/21.0.3_9',
        'bcftools/1.22'
    params:
        chromosome = lambda wc: wc.chr,
        prefix=lambda wc,output:output.vcf[:-7]
    resources:
        mem_mb_per_cpu = 4000
    threads:
        5
    shell:
        '''
        java -jar {input.conform_gt} \
        ref={input.vcf_in} \
        gt={input.widde} \
        chrom={params.chromosome} \
        match=POS \
        strict=false \
        out={params.prefix}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''


rule intersect_arrays:
    input:
        chill = OUT_DIR + 'work/data_formatting/Chillingham_updated_CHR{chr}.vcf.gz',
        widde = rules.conform_WIDDE.output.vcf,
        kerry = rules.conform_Kerry.output.vcf
    output:
        vcf = OUT_DIR + 'work/UMD/intersected_arrays_CHR{chr}.vcf.gz',
        index = OUT_DIR + 'work/UMD/intersected_arrays_CHR{chr}.vcf.gz.tbi'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    resources:
        mem_mb_per_cpu = 2000
    threads:
        10
    shell:
        '''
        bcftools isec \
        --nfiles=3 \
        --write 1 \
        --collapse none \
        --output-type u  \
        --threads {threads} \
        {input.chill} {input.widde} {input.kerry} | bcftools annotate -x FORMAT/GP -Oz --threads {threads} -o {output.vcf}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''

rule WIDDE_Chill_merge:
    input:
        Chill = rules.intersect_arrays.output.vcf,
        Widd = rules.conform_WIDDE.output.vcf,
        kerry = rules.conform_Kerry.output.vcf
    output:
        vcf = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_CHR{chr}.vcf.gz',
        tbi = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_CHR{chr}.vcf.gz.tbi'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    threads:
        8
    shell:
        '''
        bcftools merge \
        {input.Chill} {input.Widd} {input.kerry} \
        --output-type z \
        --output {output.vcf} \
        --threads {threads}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''

rule intersect_merged:
    input:
        merged = rules.WIDDE_Chill_merge.output.vcf,
        sites = rules.intersect_arrays.output.vcf
    output:
	    intersected = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected_CHR{chr}.vcf.gz', 
	    index = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected_CHR{chr}.vcf.gz.tbi'
    envmodules:
        "stack/2024-06",
        "gcc/12.2.0",
        "bcftools/1.22"
    resources:
        mem_mb_per_cpu = 2000
    threads:
        10
    shell:
        '''
        bcftools isec \
        --nfiles=2 \
        --write 1 \
        --collapse none \
        --output-type z  \
        --threads {threads} \
        {input.merged} {input.sites} -o {output.intersected}
        bcftools index --tbi --threads {threads} {output.intersected}
        '''




rule merge_conformed_WIDDE:
    input:
        vcf = expand(OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected_CHR{chr}.vcf.gz', chr=range(1,30)),
        index = expand(OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected_CHR{chr}.vcf.gz.tbi', chr=range(1,30))
    output:
        vcf = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected.vcf.gz',
        index = OUT_DIR + 'work/UMD/WIDDE_Chill_Kerry_merged_intersected.vcf.gz.tbi'
    envmodules:    
        "stack/2024-06",
        "gcc/12.2.0",
        'bcftools/1.22'
    resources:
        mem_mb_per_cpu = 4000
    threads:
        5
    shell:
        '''
        bcftools concat \
         --output-type z \
         --threads {threads} \
         {input.vcf} \
         -o {output.vcf}
        bcftools index --tbi --threads {threads} {output.vcf}
        '''

rule liftover:
    input:
        vcf = rules.merge_conformed_WIDDE.output.vcf,
        chain = config['chain_file'],
        ars_ref = config['ars_ref'],
        umd_ref = config['umd_ref'],
    output:
        vcf = OUT_DIR + 'work/ARS/WIDDE_Chill_Kerry_merged_intersected_ARS1.2.vcf.gz',
        index = OUT_DIR + 'work/ARS/WIDDE_Chill_Kerry_merged_intersected_ARS1.2.vcf.gz.tbi'

    conda:
        'bcftools'
    resources:
        mem_mb_per_cpu = 4000
    threads:
        5
    shell:
        '''
        # Ref alleles need to conform to reference genome (currenlty the conform to chillingham major allele)
        # Remove palindromic because this doesn't convert palindromic sites and throws an error
        bcftools view -e 'REF="A" & ALT="T" | REF="T" & ALT="A" | REF="C" & ALT="G" | REF="G" & ALT="C"' -Ou {input.vcf} | \
        bcftools +fixref -Ou -- -f {input.umd_ref} -m flip | \
        bcftools +liftover -Ou -- -s {input.umd_ref} -f {input.ars_ref} -c {input.chain} | \
        bcftools sort -Ou | bcftools view -m2 -M2 -v snps -Oz -o {output.vcf} -W=tbi
        '''