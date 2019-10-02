rule kallisto:
	input: 
		get_fastq
	output: 
		"results/counts/{sample}/abundance.tsv"
	log: 
		"logs/kallisto/{sample}.log"
	conda:
		"envs/kallisto.yaml"
	threads: 
		threads_mid
	params:
		idx = config["kallisto"]["idx"],
		strand = config["kallisto"]["strand"]
	shell:
		'kallisto '
		'quant '
		'{params.strand} '
		'--index {params.idx} '
		'--output-dir results/counts/{wildcards.sample} '
		'--threads {threads} '
		'{input[0]} {input[1]} 2>&1 | tee {log}'

rule star:
	input: 
		get_fastq
	output:
		wig = temp("results/bam/{sample}/{sample}.Signal.Unique.str1.out.wig"),
		wig_mult = temp("results/bam/{sample}/{sample}.Signal.UniqueMultiple.str1.out.wig"),
		bam = "results/bam/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	threads: 
		threads_high
	params:
		genome = config['star']['star_genome'],
		date = time.strftime("%Y-%m-%d")
	shell:
		'STAR '
		'--outFileNamePrefix results/bam/{wildcards.sample}/{wildcards.sample}. '
		'--outSAMtype BAM SortedByCoordinate '
		'--genomeLoad NoSharedMemory '
		'--twopassMode Basic '
		'--clip3pNbases 0 '
		'--clip5pNbases 0 '
		'--clip3pAdapterSeq - '
		'--sjdbGTFfile - '
		'--sjdbGTFchrPrefix - '
		'--quantMode GeneCounts '
		'--outFilterMultimapNmax 10 '
		'--outSAMstrandField intronMotif '
		'--outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:{wildcards.sample} DT:{params.date} '
		'--outWigType wiggle '
		'--outWigStrand Unstranded '
		'-outWigNorm RPM '
		'--runThreadN {threads} '
		'--readFilesCommand zcat '
		'--genomeDir {params.genome} '
		'--readFilesIn  {input[0]} {input[1]}'

rule wig2bigwig:
	input:
		rules.star.output.wig_mult
	output:
		"results/bw/{sample}.bw"
	conda:
		"envs/ucsc-wigtobigwig.yaml"
	params:
		chrom_sizes = config['chrom_sizes']
	shell:
		'wigToBigWig '
		'{input} '
		'{params.chrom_sizes} '
		'{output}'
