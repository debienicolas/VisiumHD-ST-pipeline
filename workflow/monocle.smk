log_path = OUTPUT_FOLDER / "logs" 

rule run_monocle:
    input:
        infercnv_adata = rules.run_infercnvpy.output.output_adata,
    output:
        monocle_output_dir = directory(OUTPUT_FOLDER / "monocle"),
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_monocle.txt"
    params:
        config = json.dumps(config["monocle"])
    log:
        log_path / "run_monocle.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_monocle.txt"
    shell:
        """
        Rscript resources/monocle3/run_monocle.R \
            --input_path {input.infercnv_adata} \
            --data_dir {INPUT_FOLDER} \
            --output_path {output.monocle_output_dir} \
            --config '{params.config}' \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """