import graphab4py

for resistance in ["1_1_1", "1_1_2", "1_1_3", "1_1_4", "1_2_1", "1_2_2", "1_2_3", "1_2_4", "2_1_1", "2_1_2", "2_1_3", "2_1_4", "2_2_1", "2_2_2", "2_2_3", "2_2_4", "3_1_1", "3_1_2", "3_1_3", "3_1_4", "3_2_1", "3_2_2", "3_2_3", "3_2_4"]:
  prj = graphab4py.Project()
  prj.create_project(
    name = resistance,
    patches = "/lud11/poppman/shared/dami/dat/spl/" + resistance + ".shp/sample_1.tif",
    habitat = 1, directory = "/lud11/poppman/shared/dami/tmp"
  )
  
  prj.create_linkset(
    disttype = "cost",
    linkname = "L1",
    threshold = 4.5,
    cost_raster = "/lud11/poppman/shared/dami/dat/res/" + resistance + ".tif"
  )
  
  prj.create_graph(graphname = "G1")
  
  prj.enable_distance_conversion(
    save_plot = "/lud11/poppman/shared/dami/fig/Distance_conversion" + resistance + ".png", max_euc = 1000
  )
  
  d = prj.convert_distance(500, regression = "log")
  
  with open("/lud11/poppman/shared/dami/tst/Distance_conversions.txt", "a") as f:
    f.write(resistance + " " + str(d) + "\n")

