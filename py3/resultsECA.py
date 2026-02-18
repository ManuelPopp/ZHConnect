import os
import argparse
import graphab4py
def parseArguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "-rst", "--result_raster",
        help = "Final output raster of the analyses.",
        type = str
        )
    parser.add_argument(
        "-residx", "--resistance_map_index",
        help = "Index of the resistance map (is then used to load a map from the dat/res folder).",
        type = str
        )
    
    args = parser.parse_args()
    
    return args

args = parseArguments()

result_raster = args.result_raster
resistance_map_index = args.resistance_map_index
resistance_map_list = [
    f for f in os.listdir("/home/poppman/sync/poppman/shared/dami/dat/res") if f.endswith(".tif")
]
resistance_map_list.sort()
resistance_map = resistance_map_list[int(resistance_map_index)]

temp_dir = os.path.join("/home/poppman/sync/poppman/shared/dami/tmp/resultsEC")

graphab4py.set_graphab("/home/poppman/sync/poppman/shared/dami/opt/")
prj = graphab4py.Project()

prj.create_project(
    name = f"{os.path.splitext(result_raster)[0]}_{os.path.splitext(resistance_map)[0]}",
    patches = os.path.join("/home/poppman/sync/poppman/shared/dami/dat/top9", result_raster),
    habitat = 1,
    directory = temp_dir,
    overwrite = True
    )

prj.create_linkset(
    disttype = "cost",
    linkname = "L1",
    threshold = 20.,
    cost_raster = os.path.join(
        "/home/poppman/sync/poppman/shared/dami/dat/res", resistance_map
        )
    )

prj.create_graph(graphname = "G1")
dist_max = prj.convert_distance(500, regression = "linzero")[0]
out = prj.calculate_metric(metric = "EC", d = dist_max, p = 0.05)
ec = out["metric_value"]

with open("/home/poppman/sync/poppman/shared/dami/tmp/resultsEC/outputECA.csv", "a") as f:
    f.write(f"{resistance_map_index},{result_raster},{resistance_map},{ec}\n")