log_path = Path(OUTPUT_FOLDER) / "logs"

rule preprocess_ficture_prior:
    input:
        reference_adata_path = config["reference_atlas"]["adata_path"]
    output:
        output_folder = directory(OUTPUT_FOLDER / "ficture_prior" / "preprocessed"),
        external_model_path = OUTPUT_FOLDER / "ficture_prior" / "preprocessed" / f"pseudo_bulk_{config['reference_atlas']['cell_type_column']}.tsv.gz",
        external_cmap_path = OUTPUT_FOLDER / "ficture_prior" / "preprocessed" / f"cmap.rgb.tsv",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "preprocess_ficture_prior.txt"
    params:
        config = config["ficture_prior"],
        cell_type_column = config["reference_atlas"]["cell_type_column"],
        reference_adata_path = config["reference_atlas"]["adata_path"]
    log:
        log_path / "preprocess_ficture_prior.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "preprocess_ficture_prior.txt"
    shell:
        """
        python resources/ficture/prior.py \
            --reference_adata_path {params.reference_adata_path} \
            --output_dir {output.output_folder} \
            --cell_type_column {params.cell_type_column}
        2>&1 | tee {log} \
        && touch {output.completed_run}
        """

rule run_ficture_prior:
    input:
        preprocessed_visium_path = lambda wildcards: (
            rules.preprocess_visium.output.output_folder if not BIN2CELL or not SEQUENTIAL else rules.preprocess_bin2cell.output.output_folder
        ),
        scalefactors_path = Path(rules.preprocess_visium.input.input_path) / "spatial" / "scalefactors_json.json",
        min_max_file =  lambda wildcards: (rules.preprocess_visium.output.min_max_file if not BIN2CELL or not SEQUENTIAL else rules.preprocess_bin2cell.output.min_max_file),
        external_model_path = rules.preprocess_ficture_prior.output.external_model_path,
        external_cmap_path = rules.preprocess_ficture_prior.output.external_cmap_path
    params:
        config = config["ficture"],
        mu_scale = 1, ### will always remain 1 as the preprocessing step takes care of scaling the input data to $\mu m$ scale
        n_factors = config["ficture"]["n_factors"],
        anchor_resolution = config["ficture"]["anchor_res"],
        train_width = config["ficture"]["train_width"],
        resolution_int = int(re.findall(r'_(\d+)um', config["input_visium_resolution"])[0]),
        cores = 24,
        min_ct_unit_dge = config["ficture"]["min_ct_unit_dge"],
        min_ct_unit_feature = config["ficture"]["min_ct_unit_feature"],
        min_ct_unit_fit = config["ficture"]["min_ct_unit_fit"]
    output:
        output_dir = directory(OUTPUT_FOLDER / "ficture_prior" / "results"),
        analysis_folder = directory(OUTPUT_FOLDER / "ficture_prior" / "results" / "analysis"),
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_ficture_prior.txt"
    log:
        log_path / "run_ficture_prior.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_ficture_prior.txt"
    shell:
        """
        ficture run_together \
            --in-tsv {input.preprocessed_visium_path}/transcripts.sorted.tsv.gz \
            --out-dir {output.output_dir} \
            --in-minmax {input.min_max_file} \
            --all \
            --mu-scale {params.mu_scale} \
            --n-factor {params.n_factors} \
            --anchor-res {params.anchor_resolution} \
            --threads 1 \
            --n-jobs 10 \
            --plot-each-factor \
            --decode-sub-um-per-pixel {params.resolution_int} \
            --decode-plot-um-per-pixel {params.resolution_int} \
            --minibatch-buffer 100 \
            --min-ct-unit-dge {params.min_ct_unit_dge} \
            --min-ct-feature {params.min_ct_unit_feature} \
            --min-ct-unit-fit {params.min_ct_unit_fit} \
            --decode-from-external-model \
            --external-model {input.external_model_path} \
            --external-cmap {input.external_cmap_path}
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """


rule ficture_prior_postprocessing:
    # This should create a anndata object and merge the factor numbers to the cell type columns 