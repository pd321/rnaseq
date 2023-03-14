rule trimgalore:
	input:
		get_fastq
	output:
		["results/qc/trimming/{sample}_val_1.fq.gz", "results/qc/trimming/{sample}_val_2.fq.gz"]
	conda:
		"../envs/trimgalore.yaml"
	log:
		"logs/qc/trimgalore/{sample}.log"
	params:
		quality = config['trimgalore']['quality'],
		stringency = config['trimgalore']['stringency'],
		e = config['trimgalore']['e']
	threads: threads_mid if threads_mid < 4 else 4
	shell:
		'trim_galore '
		'--quality {params.quality} '
		'--stringency {params.stringency} '
		'-e {params.e} '
		'--gzip '
		'--output_dir results/bam/ '
		'--cores {threads} '
		'--basename {wildcards.sample} '
		'--paired --no_report_file '
		'{input[0]} {input[1]} 2>&1 | tee {log}'

rule fastqc:
	input:
		lambda wildcards: metadata_df.at[wildcards.sample, wildcards.group]
	output:
		html="results/qc/fastqc/{sample}_{group}_fastqc.html",
		zip="results/qc/fastqc/{sample}_{group}_fastqc.zip"
	log: "logs/qc/fastqc/{sample}_{group}.log"
	threads: config_threads
	wrapper:
		"v1.23.5/bio/fastqc"

rule multiqc:
	input:
		expand(["results/qc/fastqc/{sample}_{group}_fastqc.html", 
			"results/counts/{sample}/quant.sf"], group=["r1", "r2"], sample = samples)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control")
	conda:
		"../envs/multiqc.yaml"
	threads: threads_high
	log:
		"logs/multiqc/multiqc.log"
	shell:
		'multiqc '
		'--force '
		'--outdir results/qc/multiqc '
		'--zip-data-dir results logs 2>&1 | tee {log}'
