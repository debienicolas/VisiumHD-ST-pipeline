log_path = OUTPUT_FOLDER / "logs" / "celltypist"

rule run_celltypist:
    input:
        anndata_to_annotate = rules.run_bin2cell.output.binned_adata,
    output:
        output_dir = directory( OUTPUT_FOLDER / "celltypist"),
        output_object = OUTPUT_FOLDER / "celltypist" / "celltypist_annotated.h5ad",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_celltypist.txt"
    params:
        model_path = config["celltypist"]["model_path"],
        reference_atlas_path = config["reference_atlas"]["adata_path"],
        cell_type_column = config["reference_atlas"]["cell_type_column"],
    log:
        log_path / "run_celltypist.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_celltypist.txt"
    shell:
        """ python resources/annotation/run_celltypist.py \
            --input_anndata_path {input.anndata_to_annotate} \
            --model_path {params.model_path} \
            --reference_atlas_path {params.reference_atlas_path} \
            --cell_type_column {params.cell_type_column} \
            --output_dir {output.output_dir} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """


