import json

rule run_infercnvpy:
    input:
        visium_data_path=rules.annotate_clusters.output.output_file,
        ficture_postprocessing_complete = OUTPUT_FOLDER / "completed_runs" / "ficture_postprocessing.txt",
    output:
        infercnvpy_output_path=directory(OUTPUT_FOLDER / "infercnvpy"),
        output_adata = OUTPUT_FOLDER / "infercnvpy" / "infercnvpy_output.h5ad",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_infercnvpy.txt"
    params:
        config = json.dumps(config["infercnv"])
    log:
        LOG_FOLDER / "run_infercnvpy.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_infercnvpy.txt"
    shell:
        """
        python resources/infercnvpy/infercnv_py.py \
            --input {input.visium_data_path} \
            --output {output.infercnvpy_output_path} \
            --species {SPECIES} \
            --config '{params.config}' \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """
