rule fastqc:
	input:
		lambda wildcards: metadata_df.at[wildcards.sample, wildcards.group]
	output:
		html="results/qc/fastqc/{sample}_{group}_fastqc.html",
		zip="results/qc/fastqc/{sample}_{group}_fastqc.zip"
	threads: threads_mid
	params:
		"--threads {}".format(threads_mid)
	wrapper:
		"0.36.0/bio/fastqc"

rule multiqc:
	input:
		expand(["results/qc/fastqc/{sample}_{group}_fastqc.html", 
			"results/counts/{sample}/abundance.tsv"], group=["r1", "r2"], sample = samples)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control")
	conda:
		"envs/multiqc.yaml"
	threads: threads_high
	log:
		"logs/multiqc/multiqc.log"
	shell:
		'multiqc '
		'--force '
		'--outdir results/qc/multiqc '
		'--zip-data-dir results logs 2>&1 | tee {log}'
