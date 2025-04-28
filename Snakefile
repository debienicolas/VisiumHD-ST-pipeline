### Main Snakefile

##### Create a config file for various rules  #####
configfile: "config.yaml"


#INPUT_FOLDER = "input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen"
INPUT_FOLDER = config["input_folder"]

#OUTPUT_FOLDER = Path("results") / "Human_DCIS"
OUTPUT_FOLDER = Path(config["output_folder"])

BIN2CELL = config["bin2cell"]
MANUAL_ANNOTATION = config["manual_annotation"]
SPECIES = config["species"]

SEQUENTIAL = config["sequential"]


### Input and output folder checks ###

# if BIN2CELL and not ficture_config["input_visium_resolution"] == "square_002um":
#     raise ValueError("BIN2CELL is set to True, but the input visium resolution is not set to square_002um")

# if the input folder doesn't exist, raise an error
if not os.path.exists(INPUT_FOLDER):
    raise ValueError(f"Input folder {INPUT_FOLDER} does not exist")

if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)
    os.makedirs(OUTPUT_FOLDER / "completed_runs")

LOG_FOLDER = OUTPUT_FOLDER / "logs"

### include rules ### 
include: "workflow/initial_analysis.smk"
include: "workflow/bin2cell.smk"
include: "workflow/ficture.smk"
include: "workflow/annotate.smk"
include: "workflow/infercnvpy.smk"
include: "workflow/cellchat.smk"
include: "workflow/monocle.smk"
include: "workflow/scenic.smk"
include: "workflow/spatial_analysis.smk"
include: "workflow/cell_typist.smk"

#include: "workflow/spat_gene_assoc.smk"
#include: "workflow/initial_analysis.smk"
#include: "workflow/infercnv.smk"
#include: "workflow/cellphone_db.smk"
#include: "workflow/visualization.smk"


COMPLETED_RUNS = OUTPUT_FOLDER / "completed_runs"

# rule to produce all the downstream analysis results + force the creation of the coherence scores
rule all:
    input:
        COMPLETED_RUNS / "run_infercnvpy.txt",
        COMPLETED_RUNS / "run_cellchat.txt",
        COMPLETED_RUNS / "ficture_coherence.txt",
        COMPLETED_RUNS / "run_monocle.txt",
        COMPLETED_RUNS / "pyscenic.txt",
        COMPLETED_RUNS / "spatial_analysis.txt",
        COMPLETED_RUNS / "run_ficture_prior.txt",
        COMPLETED_RUNS / "run_celltypist.txt"

# rule to go till annotation is finished
rule all_annotation:
    input:
        COMPLETED_RUNS / "run_celltypist.txt",
        COMPLETED_RUNS / "run_ficture_prior.txt",
        COMPLETED_RUNS / "annotate_clusters.txt"

# rule to get the cluster coherence scores
rule get_cluster_coherence:
    input:
        COMPLETED_RUNS / "ficture_coherence.txt"




