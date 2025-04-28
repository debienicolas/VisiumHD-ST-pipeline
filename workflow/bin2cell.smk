import json

rule run_bin2cell:
    input:
        input_dir = Path(INPUT_FOLDER) / "binned_outputs" / config["input_visium_resolution"],
        input_adata = rules.perform_initial_analysis.output.output_adata,
        spatial_dir = Path(INPUT_FOLDER) / "spatial",
        highres_input_image = Path(INPUT_FOLDER) / "highres_input.tif"
    output:
        output_dir = directory(Path(OUTPUT_FOLDER) / "bin2cell"),
        joined_adata = Path(OUTPUT_FOLDER) / "bin2cell" / "bin2cell_unified_labels.h5ad",
        binned_adata = Path(OUTPUT_FOLDER) / "bin2cell" / "bin2cell_output.h5ad",
        completed_run = Path(OUTPUT_FOLDER) / "completed_runs" / "bin2cell.txt"

    log:
        OUTPUT_FOLDER / "logs" / "bin2cell.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "bin2cell.txt"
    params:
        config = json.dumps(config["bin2cell"])
    shell: 
        """
        python resources/binToCell/run_bin2cell.py \
            --input_adata {input.input_adata} \
            --visium_dir {input.input_dir} \
            --spaceranger_spatial_dir {input.spatial_dir} \
            --highres_input_image {input.highres_input_image} \
            --output_dir {output.output_dir} \
            --config '{params.config}' \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """



