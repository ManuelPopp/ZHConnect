import graphab4py

prj = graphab4py.Project()
prj.create_project(
    name = "high_suitability",
    patches = "/lud11/poppman/shared/dami/dat/spl/high_suitability/sample_1.tif",
    habitat = 1,
    directory = "/lud11/poppman/shared/dami/tmp/", overwrite = True
    )

prj.create_linkset(
    disttype = "cost",
    linkname = "L1",
    threshold = 100,
    cost_raster = "/lud11/poppman/shared/dami/dat/res/resistance_map.tif"
    )

prj.create_graph(graphname = "G1")

prj.enable_distance_conversion(
   save_plot = "/lud11/poppman/shared/dami/fig/Distance_conversion.svg",
   max_euc = 1500,
   regression = "linzero"
   )

prj.convert_distance(500, regression = "linzero")

# Converted distances:
# 250 m = 1.69927692; 500 m = 3.39855384; 750 m = 5.09783076;
# 1000 m = 6.79710767
