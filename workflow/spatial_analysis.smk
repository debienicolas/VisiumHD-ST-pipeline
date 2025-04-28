import json

rule run_spatial_analysis:
    input:
        input_adata = rules.annotate_clusters.output.output_file,
    output:
        output_dir = directory(Path(OUTPUT_FOLDER) / "spatial_analysis"),
        completed_run = Path(OUTPUT_FOLDER) / "completed_runs" / "spatial_analysis.txt"

    log:
        OUTPUT_FOLDER / "logs" / "spatial_analysis.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "spatial_analysis.txt"
    params:
        config = json.dumps(config["spatial_analysis"])
    shell: 
        """
        """
