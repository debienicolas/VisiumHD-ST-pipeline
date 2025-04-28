import os
import spatialdata_io 

output_path = "input/DCIS_spatial.zarr"
vis = spatialdata_io.visium_hd("input/Visium_HD_Human_Breast_Cancer_Fresh_Frozen")

if not os.path.exists(output_path):
    vis.write(output_path)

print(vis)

