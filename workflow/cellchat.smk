import json


rule run_cellchat:
    input:
        seurat_obj = rules.annotate_clusters.output.output_file,
        scale_factor_path = Path(INPUT_FOLDER) / "binned_outputs" / config["input_visium_resolution"] / "spatial" / "scalefactors_json.json"
    output:
        output_dir_global_plots = directory(OUTPUT_FOLDER / "cellchat"/ "global_plots"),
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_cellchat.txt"
    params:
        main_output_dir = OUTPUT_FOLDER / "cellchat",
        config = json.dumps(config["cellchat"])
    log:
        LOG_FOLDER / "run_cellchat.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_cellchat.txt"
    shell:
        """
        Rscript resources/cellchat/run_cellchat.R \
            --input_path {input.seurat_obj} \
            --output_path {params.main_output_dir} \
            --data_dir {INPUT_FOLDER} \
            --scale_factor_path {input.scale_factor_path} \
            --species {SPECIES} \
            --config '{params.config}' \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """
