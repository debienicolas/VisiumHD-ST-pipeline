import json

rule perform_initial_analysis:
    # performs an initial analysis of the visium data
    # Will output the anndata object ("visium_data.pkl") + plots + stats
    input:
        Path(INPUT_FOLDER) / 'binned_outputs' / config["input_visium_resolution"]
    output:
        completed_output_path = Path(OUTPUT_FOLDER) / "completed_runs" / "initial_analysis.txt",
        analysis_output_path = directory(Path(OUTPUT_FOLDER) / "initial_analysis"),
        output_adata = Path(OUTPUT_FOLDER) / "initial_analysis" / "filtered_adata.h5ad"
    benchmark:
        Path(OUTPUT_FOLDER) / "benchmark" / "initial_analysis.txt"
    log:
        Path(OUTPUT_FOLDER) / "logs" / "initial_analysis.log"
    params:
        config = json.dumps(config["initial_analysis"])
    shell:
        """
            python resources/initial_analysis.py \
            --input_path {input} \
            --output_path {output.analysis_output_path} \
            --species {SPECIES} \
            --config '{params.config}' \
            2>&1 | tee {log} \
            && touch {output.completed_output_path}
        """
    
    