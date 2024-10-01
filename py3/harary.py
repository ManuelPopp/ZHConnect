#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:16:02 2024

[Description]
"""
__author__ = "Manuel"
__date__ = "Tue Feb 27 15:16:02 2024"
__credits__ = ["Manuel R. Popp"]
__license__ = "Unlicense"
__version__ = "2.0.1"
__maintainer__ = "Manuel R. Popp"
__email__ = "requests@cdpopp.de"
__status__ = "Production"

#-----------------------------------------------------------------------------|
# Imports
import os, sys, platform, shutil, glob, argparse, subprocess, datetime, signal

t_start = datetime.datetime.now()

process_ids = []

def sigterm_handler(signum, frame):
    print("Process will be terminated. Cleaning up...")
    for pid in process_ids:
        os.kill(pid, signal.SIGTERM)
    
    sys.exit(0)

signal.signal(signal.SIGTERM, sigterm_handler)

#-----------------------------------------------------------------------------|
# Arguments and settings
def parseArguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-hab", "--habitat", help = "Habitat raster.",
                        type = str)
    parser.add_argument("-d", "--maxdist",
                        help = "Maximum euclidean distance or accumulated" +
                        " cost.", type = str)
    parser.add_argument("-i", "--sample_id", help = "ID of the sample.",
                        type = int)
    parser.add_argument("-r", "--resistance", help = "Resistance raster.",
                        type = str, default = None)
    parser.add_argument("-c", "--cores", help = "Number of cores for Graphab.",
                        type = int, default = 4)
    parser.add_argument("-m", "--memory", help = "RAM provided to Graphab.",
                        type = str, default = "12g")
    parser.add_argument("-mpi", "--mpi", help = "Run in mpi environment.",
                        action = "store_true")
    parser.add_argument("-con", "--con", help = "Pixel connexity for patches.",
                        type = int, default = 8)
    parser.add_argument("-w", "--where", help = "Where the data and output" +
                        "directories are located.",
                        type = str, default = "/storage/poppman/rca")
    
    args = parser.parse_args()
    
    return args

args = parseArguments()

habitat_code = 1
habitat = args.habitat
maxdist = args.maxdist
sampleid = args.sample_id
resist = args.resistance
cores = args.cores
memory = args.memory
mpi = args.mpi
connex = args.con
print(f"MPI was set to {mpi}.")

java = "/usr/bin/java" if platform.system() == "Linux" else "java"
dir_py = os.path.dirname(os.path.realpath(__file__))
dir_main = os.path.dirname(dir_py) if args.where is None else args.where
dir_dat = os.path.join(dir_main, "dat")
dir_spl = os.path.join(dir_dat, "spl", habitat)

round_dist = round(float(maxdist), 2)
dir_out = os.path.join(dir_main, "out", f"{habitat}_{resist}{round_dist}")
dir_tmp = os.path.join(dir_main, "tmp", f"{habitat}_{resist}{round_dist}")

os.umask(0)
os.makedirs(dir_out, exist_ok = True)
os.makedirs(dir_tmp, exist_ok = True)

sample = os.path.join(dir_spl, "sample_{0}.tif".format(sampleid))

resistance = "euclid" if resist == "euclid" else os.path.join(
    dir_dat, "res", resist
    )

name_prj = "_".join(["Graph", habitat, f"{sampleid:0>{5}}"])
name_lnk = "TempLink"
name_gph = "TempGraph"
prj_main = os.path.join(dir_tmp, name_prj)

graphab = os.path.join(dir_main, "opt", "graphab-2.8.jar")

settings = {"java" : java,
            "memory" : memory,
            "cores" : cores,
            "graphab" : graphab
            }

#-----------------------------------------------------------------------------|
# Functions
def base_call(java, memory, cores, graphab, **kwargs):
    '''
    General call to Graphab.
    
    Parameters
    ----------
    java : str
        Location of the java executable.
    memory : str
        Amount of RAM Graphab is allowed to use (number with unit, e.g. "3200m").
    cores : int
        Number of CPU cores Graphab is allowed to use.
    graphab : str
        Path to the graphab .jar file.
    **kwargs : dict, any
        Arguments to append to the Graphab call.
    
    Returns
    -------
    proc_out : bytes
        Process output.
    
    '''
    if "mpi" in kwargs.keys():
        mpi = kwargs["mpi"]
        del kwargs["mpi"]
        
    else:
        mpi = False
    
    if mpi:
        cmd = ["mpirun", java, "-jar", graphab, "-mpi"]
    
    else:
        cmd = [
            java, "-Djava.awt.headless=true", f"-Xmx{memory}",
            "-jar", graphab,
            "-proc", str(cores)
            ]
    
    for key, val in kwargs.items():
        values = val if isinstance(val, list) else [val]
        cmd += ["--{0}".format(key)] + values
    
    print("Running: {}".format(" ".join(cmd)))
    
    #proc_out = subprocess.check_output(cmd)
    process = subprocess.Popen(
        cmd,
        shell = False,
        stdout = subprocess.PIPE,
        stderr=subprocess.PIPE
        )
    pid = process.pid
    print(f"Started subprocess\nProcess ID: {pid}")
    
    global process_ids
    process_ids.append(pid)
    
    proc_out, proc_err = process.communicate()
    process_ids.pop()
    print(proc_err)
    
    return proc_out

def create_project(name,
                   patches,
                   habitat,
                   nomerge = False,
                   nodata = None,
                   minarea = None,
                   maxsize = None,
                   connexity = 8,
                   directory = None,
                   *settings
                   ):
    '''
    Create a Graphab project.
    
    Parameters
    ----------
    name : str
        Name of the project.
    patches : str
        File path to the landscape raster.
    habitat : int
        Integer encoding the habitat type in the landscape raster.
    nomerge : bool, optional
        Whether not tomerge contiguous patches of dierent codes.
        The default is False.
    nodata : int, optional
        Code for NoData values. The default is None.
    minarea : int, optional
        Minimum patch size in ha. The default is None.
    maxsize : int, optional
        Max size in ha. Patches with an area exceeding maxsize will be split.
        The default is None.
    connexity : int in {4, 8}, optional
        Neighbourhood or 4 or 8 pixels in patch definition. The default is 4.
    directory : str, optional
        Directory into which the project is to be created. The default is None.
    *settings : dict
        Dictionary containing Graphab settings.
    
    Returns
    -------
    out : dict
        A dictionary containing process output and project name.
    '''
    directory = os.getcwd() if directory is None else directory
    
    project_settings = [name,
                        patches,
                        f"habitat={habitat}",
                        f"dir={directory}"
                        ]
    
    if nomerge:
        project_settings += ["nomerge"]
    
    if isinstance(nodata, int):
        project_settings += [f"nodata={nodata}"]
    
    if isinstance(minarea, int):
        project_settings += [f"minarea={minarea}"]
    
    if isinstance(maxsize, int):
        project_settings += [f"maxsize={maxsize}"]
    
    if connexity == 8:
        project_settings += ["con8"]
    
    proc_out = base_call(java, memory, cores, graphab,
                         create = project_settings
                         )
    
    out = {"process_output" : proc_out,
           "project_file" : os.path.join(directory, name, name + ".xml")
           }
    
    return out

def create_linkset(project, disttype, linkname, threshold, complete = True,
                   cost_raster = None, *settings):
    '''
    Create a linkset.
    
    Parameters
    ----------
    project : str
        Path to a Graphab project .xml file.
    disttype : str
        Type of distance to use. Either "euclid" or "cost".
    linkname : str
        Name of the linkset.
    threshold : int
        Maximum distance or maximum accumulated cost (depending on the type of
        distance).
    complete : bool, optional
        Whether to create a complete linkset. The default is True.
    cost_raster : str, optional
        Path to an external cost raster file (*.tif). The default is None.
    *settings : dist
        Dictionary containing Graphab settings.
    
    Returns
    -------
    proc_out : bytes
        Process output.
    '''
    link_settings = [f"distance={disttype}", f"name={linkname}"]
    
    if complete:
        link_settings += ["complete"]
    
    link_settings += [f"maxcost={threshold}"]
    
    if cost_raster is not None:
        link_settings += [f"extcost={cost_raster}"]
    
    proc_out = base_call(java, memory, cores, graphab, project = project,
                         linkset = link_settings)
    
    return proc_out

def create_graph(graphname, project, linkset, nointra = True,
                 threshold = None, *settings):
    '''
    Create a graph.
    
    Parameters
    ----------
    graphname : str
        Graph name.
    project : str
        Path to a Graphab project .xml file.
    linkset : str
        Name of the linkset.
    nointra : bool, optional
        Set the "nointra" option. The default is True.
    threshold : int, optional
        Maximum distance or maximum accumulated cost (depending on the type of
        distance). The default is None.
    *settings : dist
        Dictionary containing Graphab settings.
    
    Returns
    -------
    proc_out : bytes
        Process output.
    '''
    graph_settings = [f"name={graphname}"]
    
    if nointra:
        graph_settings += ["nointra"]
    
    if threshold is not None:
        graph_settings += [f"threshold={threshold}"]
    
    proc_out = base_call(java, memory, cores, graphab, project = project,
                         uselinkset = linkset, graph = graph_settings)
    
    return proc_out

def calculate_metric(project, linkset, graph, metric, mtype = "global",
                     **metric_args):
    '''
    Calculate a global metric.
    
    Parameters
    ----------
    project : str
        Path to a Graphab project .xml file.
    linkset : str
        Name of the linkset.
    graph : str
        Graph name.
    metric : str
        Metric name.
    mtype : str {local, global}
        Metric type.
    **metric_args : dict
        Metric paramneters.
    
    Returns
    -------
    out : dict
        A dictionary containing process output and project name.
    '''
    metric_settings = [metric]
    
    for key, val in metric_args.items():
        metric_settings += ["{0}={1}".format(key, val)]
    
    if mtype == "global":
        proc_out = base_call(java, memory, cores, graphab, project = project,
                             uselinkset = linkset, usegraph = graph,
                             gmetric = metric_settings)
        
        out_text = proc_out.decode("utf-8")
        
        metric_value = float(out_text.split(" ")[-1].strip())
        
        out = {"process_output" : proc_out,
               "metric_value" : metric_value}
        
    elif mtype == "local":
        proc_out = base_call(java, memory, cores, graphab, project = project,
                             uselinkset = linkset, usegraph = graph,
                             lmetric = metric_settings)
        
        out = proc_out.decode("utf-8")
        
    elif mtype == "component":
        proc_out = base_call(java, memory, cores, graphab, project = project,
                             uselinkset = linkset, usegraph = graph,
                             cmetric = metric_settings)
        
        out = proc_out.decode("utf-8")
    else:
        raise Exception(f"Invalid metric type '{mtype}'.")
    
    return out

def delta_by_item(project, linkset, graph, metric, select = None,
                select_from_file = None, obj = "patch", mpi = False,
                **metric_args):
    '''
    Calculate a global metric in delta mode on patches or links depending on
    obj parameter for the selected graph.
    
    Parameters
    ----------
    project : str
        Path to a Graphab project .xml file.
    linkset : str
        Name of the linkset.
    graph : str
        Graph name.
    metric : str
        Metric name.
    select : list, optional
        Restrict the calculation to items (patches or links) listed by
        identifier. The default is None.
    select_from_file : str, optional
        Restrict the calculations on items listed in a *.txt file. The file
        must contain one identifier per line. The default is None.
    obj : str {patch, link}, optional
        Type of objects to remove. The default is "patch".
    mpi : bool, optional
        Run in MPI mode (on cluster).
    **metric_args : dict
        Metric paramneters.
    
    Returns
    -------
    out : dict
        A dictionary containing process output and project name.
    '''
    delta_settings = [metric]
    
    for key, val in metric_args.items():
        delta_settings += ["{0}={1}".format(key, val)]
    
    if select is not None:
        delta_settings += ["sel=" + ",".join(select)]
    
    if select_from_file is not None:
        delta_settings += [f"fsel={select_from_file}"]
    
    delta_settings += [f"obj={obj}"]
    
    proc_out = base_call(java, memory, cores, graphab, project = project,
                         uselinkset = linkset, usegraph = graph, mpi = mpi,
                         delta = delta_settings)
    
    return proc_out

#-----------------------------------------------------------------------------|
# Calculate metrics for a variant
## Create info file
info_file = os.path.join(dir_out, "job.info")

if not os.path.isfile(info_file):
    try:
        with open(os.path.join(dir_spl, "patch_meta.info"), "r") as f:
            content = f.readlines()
        
        content.append("Resistance: " + resistance)
        
        with open(info_file, "w") as f:
            f.writelines(content)
    
    except Exception as e:
        print(e)
        print("Unable to write info file.")

## Create project
output = create_project(
    name_prj, sample, habitat_code, directory = dir_tmp, connexity = connex
    )

proj = output["project_file"]

print("Project created.")

## Create linkset
if resist == "euclid":
    create_linkset(
        proj, "euclid", name_lnk, maxdist, complete = True, cost_raster = None
        )
else:
    crst = resistance if resistance.endswith(".tif") else resistance + ".tif"
    create_linkset(
        proj, "cost", name_lnk, maxdist, complete = True,
        cost_raster = crst
        )

print("Linkset created.")

## Create graph
create_graph(name_gph, proj, name_lnk, nointra = True,
                 threshold = maxdist)

t_graph = datetime.datetime.now()
t_delta = t_graph - t_start
print(f"Graph created after {t_delta}.")

## Calculate metrices
### Global metrices
'''
gEC = calculate_metric(proj, name_lnk, name_gph, "EC", d = maxdist, p = 0.05)
gIIC = calculate_metric(proj, name_lnk, name_gph, "IIC")
gH = calculate_metric(proj, name_lnk, name_gph, "H")

### Local metrices
lIF = calculate_metric(proj, name_lnk, name_gph, "IF", mtype = "local",
                       d = maxdist, p = 0.5, beta = 1
                       )

F = calculate_metric(proj, name_lnk, name_gph, "F", mtype = "local",
                       d = maxdist, p = 0.5, beta = 1
                       )

Ec = calculate_metric(proj, name_lnk, name_gph, "Ec", mtype = "local")
'''

### Delta metrices
delta_by_item(proj, name_lnk, name_gph, metric = "H", mpi = mpi)
########################################
# Inserted time calculation
t_H = datetime.datetime.now()
t_delta = t_H - t_start
print(f"Delta H calculated in {t_delta}.")
########################################

'''
output = [
    ["Project", "Metric", "Scope", "Value"],
    [name_prj, "EC", "global", gEC["metric_value"]],
    [name_prj, "IIC", "global", gIIC["metric_value"]],
    [name_prj, "H", "global", gH["metric_value"]]
    ]

with open(os.path.join(dir_out, f"Metrices_{sampleid}.csv"),
          "w", encoding = "UTF8") as f:
    writer = csv.writer(f)
    writer.writerows(output)
'''

if glob.glob(os.path.join(prj_main, "patches.*")) == []:
    raise Exception("No patches files found.")

for src in glob.glob(os.path.join(prj_main, "patches.*")):
    if not src.endswith("tif"):
        src_d, src_f = os.path.split(src)
        src_n, src_e = os.path.splitext(src_f)
        dst = os.path.join(dir_out, f"{src_n}_{sampleid}{src_e}")
        shutil.move(src, dst)

if glob.glob(os.path.join(prj_main, "delta*")) == []:
    raise Exception("No delta files found.")

for src in glob.glob(os.path.join(prj_main, "delta*")):
    src_d, src_f = os.path.split(src)
    src_n, src_e = os.path.splitext(src_f)
    dst = os.path.join(dir_out, f"{src_n}_{sampleid}{src_e}")
    shutil.move(src, dst)
    print(f"Copied files from {src} to {dst}.")

shutil.rmtree(prj_main)
t_end = datetime.datetime.now()
t_delta = t_end - t_start

print(f"Process finished after {t_delta}.")
