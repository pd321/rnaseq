rule go:
	input: 
		degxls = 'results/deg/{contrast}.xls'
	output: 
		gores = 'results/go/{contrast}.xls'
	threads: threads_low
	log: 'logs/go/{contrast}.log'
	params:
		genome = config['genome'],
		basemean_cutoff = 5,
		padj_cutoff = 0.1
	script:
		"R/go.R"
