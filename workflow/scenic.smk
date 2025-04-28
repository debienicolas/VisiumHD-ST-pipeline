import json

rule run_pyscenic:
    input:
        input_adata = rules.annotate_clusters.output.output_file,
    output:
        output_dir = directory(Path(OUTPUT_FOLDER) / "pyscenic"),
        completed_run = Path(OUTPUT_FOLDER) / "completed_runs" / "pyscenic.txt"

    log:
        OUTPUT_FOLDER / "logs" / "pyscenic.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "pyscenic.txt"
    params:
        config = json.dumps(config["pyscenic"])
    shell: 
        """
        """
