rule salmon:
	input: 
		r1 = rules.trimgalore.output[0]
		r2 = rules.trimgalore.output[1]
	output: 
		quant = "results/counts/{sample}/quant.sf"
		quant_genes = "results/counts/{sample}/quant.genes.sf"
	log: 
		"logs/salmon/{sample}_quant.log"
	conda:
		"../envs/salmon.yaml"
	threads: 
		threads_mid
	params:
		idx = config["salmon"]["idx"],
		gtf = config["salmon"]["gtf"],
		libtype = config["salmon"]["strand"]
	shell:
		'salmon quant '
		'--libType {params.libtype} '
		'--index {params.idx} '
		'--geneMap {params.gtf} '
		'--output-dir results/counts/{wildcards.sample} '
		'--threads {threads} '
		'--mates1 {input.r1} '
		'--mates2 {input.r2} 2>&1 | tee {log}'

rule star:
	input: 
		rules.trimgalore.output
	output:
		wig = temp("results/bam/{sample}/{sample}.Signal.Unique.str1.out.wig"),
		wig_mult = temp("results/bam/{sample}/{sample}.Signal.UniqueMultiple.str1.out.wig"),
		bam = "results/bam/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	threads: 
		threads_high
	conda:
		"../envs/star.yaml"
	params:
		idx = config['star']['idx'],
		date = time.strftime("%Y-%m-%d")
	shell:
		'STAR '
		'--outFileNamePrefix results/bam/{wildcards.sample}/{wildcards.sample}. '
		'--outSAMtype BAM SortedByCoordinate '
		'--genomeLoad NoSharedMemory '
		'--twopassMode Basic '
		'--outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:{wildcards.sample} DT:{params.date} '
		'--outWigType wiggle '
		'--outWigStrand Unstranded '
		'--outWigNorm RPM '
		'--runThreadN {threads} '
		'--readFilesCommand zcat '
		'--genomeDir {params.idx} '
		'--readFilesIn  {input[0]} {input[1]}'

