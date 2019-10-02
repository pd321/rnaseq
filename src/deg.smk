rule deseq2:
	input:
		counts = get_deg_input,
		tx2gene = config['tx2gene'],
		geneinfo = config['geneinfo']
	output:
		deg = expand("results/deg/{contrast}.xls", contrast = contrasts),
		pca_plot = "results/deg/pca.pdf"
	conda:
		"envs/deseq2.yaml"
	threads: 
		threads_high
	params:
		contrasts = config["contrasts"]
	script:
		"R/deseq2.R"
