rule wig2bigwig:
	input:
		rules.star.output.wig_mult
	output:
		"results/bw/{sample}.bw"
	conda:
		"../envs/ucsc-wigtobigwig.yaml"
	params:
		chrom_sizes = config['wig2bigwig']['chrom_sizes']
	shell:
		'wigToBigWig '
		'{input} '
		'{params.chrom_sizes} '
		'{output}'
