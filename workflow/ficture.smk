import json
from pathlib import Path
#configfile: "config/ficture.yaml"
import re

log_path = OUTPUT_FOLDER / "logs" / "ficture"

### Wildcard constraints ###
# based on the provided config file values, now supports multiple values for each parameter
wildcard_constraints:
    possible_n_factors = "nF(" + "|".join(config["ficture"]["n_factors"].split(",")) + ")",
    possible_train_width = "d_(" + "|".join(config["ficture"]["train_width"].split(",")) + ")",
    possible_anchor_res = "r_(" + "|".join(config["ficture"]["anchor_res"].split(",")) + ")_\d+\.\d+"



### Preprocess the visium data ###
rule preprocess_visium:
    input:
        input_path=Path(INPUT_FOLDER) / "binned_outputs" / config["input_visium_resolution"],
        input_adata = rules.perform_initial_analysis.output.output_adata
    output:
        output_folder = directory(OUTPUT_FOLDER / "ficture" / "preprocessed_visium"),
        transcripts_sorted = OUTPUT_FOLDER / "ficture" / "preprocessed_visium" / "transcripts.sorted.tsv.gz",
        transcripts_unsorted = OUTPUT_FOLDER / "ficture" / "preprocessed_visium" / "transcripts.unsorted.tsv.gz",
        tissue_positions = OUTPUT_FOLDER / "ficture" / "preprocessed_visium" / "tissue_positions.csv.gz",
        min_max_file = OUTPUT_FOLDER / "ficture" / "preprocessed_visium" / "minmax.tsv",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "preprocess_visium.txt"
    log:
        log_path / "preprocess_visium.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "preprocess_visium.txt"
    shell:
        """
        python resources/ficture/format_visium.py \
            --input {input.input_path} \
            --input_adata {input.input_adata} \
            --output {output.output_folder} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """

rule preprocess_bin2cell:
    input:
        bin2cell_output = Path(OUTPUT_FOLDER) / "bin2cell" / "bin2cell_output.h5ad",
        visium_binned_dir = Path(INPUT_FOLDER) / "binned_outputs" / "square_002um"
    output:
        output_folder = directory(Path(OUTPUT_FOLDER) / "ficture" / "preprocessed_bin2cell"),
        transcripts_sorted = Path(OUTPUT_FOLDER) / "ficture" / "preprocessed_bin2cell" / "transcripts.sorted.tsv.gz",
        transcripts_unsorted = Path(OUTPUT_FOLDER) / "ficture" / "preprocessed_bin2cell" / "transcripts.unsorted.tsv.gz",
        tissue_positions = Path(OUTPUT_FOLDER) / "ficture" / "preprocessed_bin2cell" / "tissue_positions.csv.gz",
        min_max_file = Path(OUTPUT_FOLDER) / "ficture" / "preprocessed_bin2cell" / "minmax.tsv",
        completed_run = Path(OUTPUT_FOLDER) / "completed_runs" / "preprocess_bin2cell.txt"
    log:
        log_path / "preprocess_bin2cell.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "preprocess_bin2cell.txt"
    shell:
        """
        python resources/ficture/format_bin2cell.py \
            --input {input.bin2cell_output} \
            --visium_binned_dir {input.visium_binned_dir} \
            --output {output.output_folder} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """





### Run FICTURE ###
rule run_ficture:
    input:
        preprocessed_visium_path = lambda wildcards: (
            rules.preprocess_visium.output.output_folder if not BIN2CELL or not SEQUENTIAL else rules.preprocess_bin2cell.output.output_folder
        ),
        scalefactors_path = Path(rules.preprocess_visium.input.input_path) / "spatial" / "scalefactors_json.json",
        min_max_file =  lambda wildcards: (rules.preprocess_visium.output.min_max_file if not BIN2CELL or not SEQUENTIAL else rules.preprocess_bin2cell.output.min_max_file),
        tabix_path = "/Users/nicolasdebie/miniconda3/envs/st/bin/tabix" if os.path.exists("/Users/nicolasdebie/miniconda3/envs/st/bin/tabix") else "/hpc/home/nrd20/miniconda3/envs/st/bin/tabix",
        bgzip_path = "/Users/nicolasdebie/miniconda3/envs/st/bin/bgzip" if os.path.exists("/Users/nicolasdebie/miniconda3/envs/st/bin/bgzip") else "/hpc/home/nrd20/miniconda3/envs/st/bin/bgzip",
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
        output_dir = directory(OUTPUT_FOLDER / "ficture" / "ficture_output"),
        analysis_folder = directory(OUTPUT_FOLDER / "ficture" / "ficture_output" / "analysis"),
        completed_run = OUTPUT_FOLDER / "completed_runs" / "run_ficture.txt",
        main_output_files = expand(str(OUTPUT_FOLDER / "ficture" / "ficture_output/analysis/nF{n_factors}.d_{train_width}/nF{n_factors}.d_{train_width}.decode.prj_{train_width}.r_{anchor_res}_{anchor_res_plus_one}.pixel.sorted.tsv.gz"),
            n_factors=config["ficture"]["n_factors"].split(","),
            train_width=config["ficture"]["train_width"].split(","),
            anchor_res=[str(float(x)) for x in config["ficture"]["anchor_res"].split(",")],
            anchor_res_plus_one=[str(float(x) + 1) for x in config["ficture"]["anchor_res"].split(",")])
    log:
        log_path / "run_ficture.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "run_ficture.txt"
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
            --train-width {params.train_width} \
            --threads 1 \
            --n-jobs 10 \
            --plot-each-factor \
            --decode-sub-um-per-pixel {params.resolution_int} \
            --decode-plot-um-per-pixel {params.resolution_int} \
            --minibatch-buffer 100 \
            --min-ct-unit-dge {params.min_ct_unit_dge} \
            --min-ct-feature {params.min_ct_unit_feature} \
            --min-ct-unit-fit {params.min_ct_unit_fit} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """
        #    --bgzip {input.bgzip_path} \
        #    --tabix {input.tabix_path} \



rule ficture_coherence:
    input:
        ficture_output_file = rules.run_ficture.output.completed_run,
        analysis_folder = OUTPUT_FOLDER / "ficture" / "ficture_output" / "analysis"
    output:
        output_file = OUTPUT_FOLDER / "ficture" / "coherence.tsv",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "ficture_coherence.txt",
        chosen_id = OUTPUT_FOLDER / "ficture" / "chosen_id.txt"
    log:
        log_path / "ficture_coherence.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "ficture_coherence.txt"
    shell:
        """
        python -m resources.ficture.coherence \
            --input_path {input.analysis_folder} \
            --output_path {output.output_file} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """



rule ficture_postprocessing:
    input:
        # ficture_output_file = expand(Path("results") / OUTPUT_FOLDER / "ficture" / f"ficture_output/analysis/nF{ficture_config['n_factors']}.d_{ficture_config['train_width']}/nF{ficture_config['n_factors']}.d_{ficture_config['train_width']}.decode.prj_{ficture_config['train_width']}.r_{int(ficture_config['anchor_res'])}_{int(ficture_config['anchor_res'])+1}.pixel.sorted.tsv.gz",
        #             n_factors=ficture_config['n_factors'].split(","),
        #             train_width=ficture_config['train_width'].split(","),
        #             anchor_res=ficture_config['anchor_res']),
        ficture_output_analysis_folder = rules.run_ficture.output.output_dir + "/analysis",
        preprocessed_ficture_file = lambda wildcards: (
            rules.preprocess_visium.output.transcripts_sorted if not BIN2CELL or not SEQUENTIAL else rules.preprocess_bin2cell.output.transcripts_sorted
        ),
        visium_data_path = lambda wildcards: (
            Path(INPUT_FOLDER) / "binned_outputs" / config["input_visium_resolution"] if not BIN2CELL or not SEQUENTIAL else rules.run_bin2cell.output.binned_adata
        ),
        ficture_output_file = rules.run_ficture.output.completed_run,
        cluster_coherence_file = rules.ficture_coherence.output.output_file,
        scale_factors = Path(INPUT_FOLDER) / "binned_outputs" / config["input_visium_resolution"] / "spatial" / "scalefactors_json.json",
        chosen_id = rules.ficture_coherence.output.chosen_id
    output:
        output_file_path = OUTPUT_FOLDER / "ficture" / "postprocessed_ficture" / "transcripts_ficture_joined.tsv.gz",
        output_anndata_path = OUTPUT_FOLDER / "ficture" / "postprocessed_ficture" / "ficture_anndata.h5ad",
        completed_run = OUTPUT_FOLDER / "completed_runs" / "ficture_postprocessing.txt",
        # chosen_id = OUTPUT_FOLDER / "ficture" / "postprocessed_ficture" / "chosen_id.txt"
    log:
        log_path / "ficture_postprocessing.log"
    benchmark:
        OUTPUT_FOLDER / "benchmark" / "ficture_postprocessing.txt"
    params:
        output_folder = OUTPUT_FOLDER / "ficture" / "postprocessed_ficture"
    shell:
        """
        python -m resources.ficture.postproc_ficture \
            --ficture_output_analysis_folder {input.ficture_output_analysis_folder} \
            --preprocessed_ficture_file {input.preprocessed_ficture_file} \
            --input_path_visium {input.visium_data_path} \
            --output_folder {params.output_folder} \
            --cluster_coherence_file {input.cluster_coherence_file} \
            --chosen_id {input.chosen_id} \
            --scale_factors {input.scale_factors} \
            2>&1 | tee {log} \
            && touch {output.completed_run}
        """
            #--output_file_ficture {input.ficture_output_file} \


