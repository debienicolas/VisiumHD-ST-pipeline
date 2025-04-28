
rule merge_bin2cell_and_ficture:
    input:
        bin2cell_file = rules.run_bin2cell.output.binned_adata,
        ficture_file = rules.ficture_postprocessing.output.output_anndata_path
    output:
        merged_file = OUTPUT_FOLDER / "ficture" / "merged_annotated_anndata.h5ad"
    log:
        OUTPUT_FOLDER / "logs" / "merge_bin2cell_and_ficture.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "merge_bin2cell_and_ficture.txt"
    run:
        import anndata as ad

        # read both files
        bin2cell_adata = ad.read_h5ad(input.bin2cell_file)
        print(bin2cell_adata)
        ficture_adata = ad.read_h5ad(input.ficture_file)
        print(ficture_adata)

        print(bin2cell_adata.obs.index)
        print(bin2cell_adata.obs_names)
        print(bin2cell_adata.obs["object_id"])
        print("--------------------------------")
        print(ficture_adata.obs.index)
        print(ficture_adata.obs_names)

        # merge the two adata objects
        filtered_ficture_adata = ficture_adata[ficture_adata.obs_names.isin(bin2cell_adata.obs_names)]
        # add the factor column to the bin2cell adata
        bin2cell_adata.obs["factor"] = filtered_ficture_adata.obs["factor"]

        # save the merged adata object
        bin2cell_adata.write_h5ad(output.merged_file, compression="gzip")





rule annotate_clusters:
    input:
        input_file = lambda wildcards: (
            rules.merge_bin2cell_and_ficture.output.merged_file if not SEQUENTIAL else rules.ficture_postprocessing.output.output_anndata_path
        ),
        # input_file = rules.ficture_postprocessing.output.output_anndata_path,
        # input_file = Path("results") / OUTPUT_FOLDER / "ficture" / "postprocessed_ficture" / "ficture_anndata.h5ad",
        ficture_output_folder = OUTPUT_FOLDER / "ficture" / "ficture_output",
        chosen_id = rules.ficture_coherence.output.chosen_id
    output:
        completed_run = OUTPUT_FOLDER / "completed_runs" / "annotate_clusters.txt",
        output_file = OUTPUT_FOLDER / "ficture" / "annotated_anndata.h5ad"
    params:
        manual_annot = "--man_annot" if MANUAL_ANNOTATION else ""
    log:
        log_path / "annotate_clusters.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "annotate_clusters.txt"
    shell:
        """
        python resources/annotate_clusters.py \
            --input_file {input.input_file} \
            --ficture_output_folder {input.ficture_output_folder} \
            --chosen_id {input.chosen_id} \
            --output_file {output.output_file} \
            {params.manual_annot}
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """