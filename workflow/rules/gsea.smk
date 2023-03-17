rule gsea_rnk_create:
	input: 
		degxls = 'results/deg/{contrast}.xls'
	output: 
		gsearnk = temp('results/gsea/{contrast}.rnk')
	threads: 1
	log: 'logs/gsea/{contrast}_rnk_create.log'
	params:
		mean_fpkm_cutoff = 0.5
	script:
		"R/gsea_rnk_create.R"

rule gsea:
	input: rules.gsea_rnk_create.output
	output: 'results/gsea/{contrast}_{collection}.done'
	threads: 1
	params:
		gmt_file = get_gmt_file,
		jar_loc = config['gsea']['jar_loc']
	log: "logs/gsea/{contrast}_{collection}_gsea.log"
	shell: 
		'java '
		'-cp {params.jar_loc} '
		'-Xmx8g xtools.gsea.GseaPreranked '
		'-gmx {params.gmt_file} '
		'-norm meandiv '
		'-nperm 1000 '
		'-rnk {input} '
		'-scoring_scheme weighted '
		'-rpt_label {wildcards.contrast}_{wildcards.collection} '
		'-create_svgs false '
		'-make_sets true '
		'-plot_top_x 100 '
		'-rnd_seed 1024 '
		'-set_max 500 '
		'-set_min 15 '
		'-zip_report false '
		'-out results/gsea '
		'-gui false 2>&1 | tee {log} && touch {output}'
