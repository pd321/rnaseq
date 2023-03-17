rule deseq2:
	input:
		counts = get_deg_input,
		tx2gene = config['deseq2']['tx2gene'],
		geneinfo = config['deseq2']['geneinfo'],
		metadata_file = metadata_file
	output:
		deg = expand("results/deg/{contrast}.xls", contrast = contrasts),
		pca_plot = "results/deg/pca.pdf"
	log: "logs/deseq.log"
	conda:
		"../envs/deseq2.yaml"
	threads: 1
	params:
		contrasts = config['deseq2']["contrasts"]
	script:
		"../scripts/deseq2.R"
