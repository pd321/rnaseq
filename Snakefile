include: "src/common.smk"
include: "src/qc.smk"
include: "src/align.smk"
include: "src/deg.smk"
include: "src/go.smk"
include: "src/gsea.smk"

deg = expand("results/deg/{contrast}.xls", contrast = contrasts) + ["results/deg/pca.pdf"]
qc = ["results/qc/multiqc/multiqc_report.html"]
out_files = deg + qc

if config['addons']['go']:
	out_files += expand("results/go/{contrast}.xls", contrast = contrasts)

if config['addons']['star']:
	out_files += expand("results/bw/{sample}.bw", contrast = contrasts)

if config['addons']['gsea']:
	out_files += expand("results/gsea/{contrast}_{collection}.done", contrast = contrasts, collection = config['gsea']['collections'].split(","))

rule all:
	input:
		out_files

onsuccess:
	shell("if command -v telegram-notify; then telegram-notify --success --text \'Snakemake:rnaseq:{} Completed\'; fi".format(foldername))

onerror:
	shell("if command -v telegram-notify; then telegram-notify --error --text \'Snakemake:rnaseq:{} Failed\'; fi".format(foldername))
