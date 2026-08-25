#!/usr/bin/python
'''
This tool is used to test the IP-Glasma+JIMWLK+SubnucleonDiffraction framework for vector meson production. 
It runs a series of simulations with different random seeds, collects the resulting cross section data, 
and compares it to reference data. The reference data has been calculated using a large number of events and is stored in the file "validation_spectra",
using a verified version of the code. 

The computed cross section should also match the published results in https://arxiv.org/pdf/2207.03712 Fig. 1 ("CGC+shape fluct") curve,
except that by default the code uses a smaller lattice for performance reasons, which may lead to small deviations in the cross section. 
The lattice size can be increased by modifying the input file "input_vm_proton".

The script also generates plots  color field (tr V(x)) for selected events.

Example usage:
    python3 test_vector_meson_production.py -maxevents 100 -datadir path_for_cross_sections -subnucleondiffraction_path path_to_subnucleondiffraction_directory

The comparison plot will be saved in the specified directory as "cross_section.pdf". The color field plots will be saved as "output_<event_id>.pdf".

Use -plot_only to skip the simulation and only generate the comparison plots from existing data files in the specified directory.

Before running this, make sure you have downloaded and built the subnucleondiffraction code from https://github.com/hejajama/subnucleondiffraction
'''

import subprocess
import os
import argparse
import numpy as np
import tempfile
import matplotlib.pyplot as plt
subprocess.run("rm Initial_x_*", shell=True)
subprocess.run("rm JIMWLKSnapshot_x_*", shell=True)

subnucleondiffraction_path = "./subnucleondiffraction/"  # can be changed by -subnucleondiffraction_path argument, points to the subnucleondiffraction code
ipglasma_cmd = "../ipglasma"
subprocess.run("cp ../qs2Adj_vs_Tp_vs_Y_200.in .", shell=True)



xpom="0.001705"

reference_cross_section_file = "validation_spectra"



plot_ids = [1,5]

def generate_temp_input(source_input_file="input_vm_proton", seed=0, x_pom="0.001705"):
    with open(source_input_file, "r") as f:
        lines = f.readlines()

    temp_fd, temp_path = tempfile.mkstemp(prefix="input_vm_proton_", suffix=".in", dir=".")
    os.close(temp_fd)

    with open(temp_path, "w") as f:
        for line in lines:
            if line.lstrip().startswith("seed"):
                f.write(f"seed {seed}\n")
            elif line.lstrip().startswith("x_projectile_jimwlk"):
                f.write(f"x_projectile_jimwlk {x_pom}\n")
            elif line.lstrip().startswith("x_target_jimwlk"):
                f.write(f"x_target_jimwlk {x_pom}\n")
            else:
                f.write(line)

    return temp_path


def remove_temp_input(temp_input_path):
    if temp_input_path and os.path.exists(temp_input_path):
        os.remove(temp_input_path)



def PlotComparison(datadir=""):
    print("Reading data from ", datadir, " and comparing to reference cross section data in ", reference_cross_section_file)
    subprocess.run(f"python3 {subnucleondiffraction_path}/python/cross_section.py -maxconf 999 -dir {datadir} -ds 0.005 -steps 36 > {datadir}/cross_section.dat", shell=True)
    events_included = None
    try:
        with open(f"{datadir}/cross_section.dat", "r") as f:
            for line in f:
                if line.startswith("#") and "events included" in line:
                    events_included = int(line.split(":")[-1].strip())
                    break
    except Exception:
        pass

    xs = np.loadtxt(f"{datadir}/cross_section.dat")
    if len(xs) == 0:
        print("No cross section data found in cross_section.dat")
        return
    coh = xs[:, 1]*1e6
    coherr = xs[:, 2]*1e6
    incoh = xs[:, 3]*1e6
    incoherr = xs[:, 4]*1e6
    plt.plot(xs[:, 0], coh, label=f"Test run coh", color="black", linestyle="solid")
    plt.fill_between(xs[:, 0], coh - coherr, coh + coherr, color="gray", alpha=0.5)
    plt.plot(xs[:,0], incoh, label=f"Test run incoh", color="red", linestyle="solid")
    plt.fill_between(xs[:, 0], incoh - incoherr, incoh + incoherr, color="pink", alpha=0.5)

    

    validation = np.loadtxt(reference_cross_section_file)
    plt.errorbar(validation[:, 0], validation[:, 1]*1e6, label=f"Reference coh", yerr=validation[:, 2]*1e6, fmt="o", color="blue", linestyle="dashed")  
    plt.errorbar(validation[:, 0], validation[:, 3]*1e6, label=f"Reference incoh", yerr=validation[:, 4]*1e6, fmt="o", color="orange", linestyle="dashed")
    
    plt.xlabel("t [GeV^2]")
    plt.ylabel(r"$\frac{d\sigma}{dt}$ [$\mu \mathrm{b}/\mathrm{GeV}^2$]")
    plt.ylim(bottom=0.1,top=5e2)
    plt.xlim(left=0, right=2)
    plt.yscale("log")
    plt.title(r"Test run cross section, $N_\mathrm{events} = " + str(events_included) + r", x_\mathbb{P} = " + str(xpom) + "$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{datadir}/cross_section.pdf")
    print("Created cross section plot: ", f"{datadir}/cross_section.pdf")
    plt.close()


parser = argparse.ArgumentParser(description="Test IP-Glasma+JIMWLK+SubnucleonDiffraction for vector meson production.")
parser.add_argument("-plot_only", action="store_true", help="Use existing data files for plotting the comparison.")
parser.add_argument("-maxevents", type=int, default=200, help="Maximum number of events for the runs.")
parser.add_argument("-datadir", type=str, default="./jpsi/", help="Directory to store data files.")
parser.add_argument("-subnucleondiffraction_path", type=str, default=subnucleondiffraction_path, help="Path to the subnucleondiffraction executable.")
args = parser.parse_args()
maxseed = int(args.maxevents/2)
datadir = args.datadir
subnucleondiffraction_path = args.subnucleondiffraction_path

subnucleondiffraction_cmd = f"{subnucleondiffraction_path}/build/bin/subnucleondiffraction -Q2 0 -wavef_file {subnucleondiffraction_path}/gauss-boosted_mzsat.dat -Q2 0 -maxt 2.01 -tstep 0.05 -mcintpoints 2e5"

if args.plot_only:
    PlotComparison(datadir=datadir)
    raise SystemExit(0)

subprocess.run(f"mkdir -p {datadir}", shell=True)
#subprocess.run(f"rm {datadir}/spectra_*", shell=True) # commented out for now in the testing phase
#subprocess.run(f"rm {datadir}/output_*", shell=True)
subprocess.run(f"cp input_vm_proton {datadir}/", shell=True)
subprocess.run(f"rm -f {datadir}/ipglasma_log", shell=True)

for s in range(0, maxseed):
    print("===== RUNNING SEED ", s, " =====")

    #### Update config file
    temp_input_path = generate_temp_input(seed=s, x_pom=xpom)

    with open(temp_input_path, "r") as f:
        lines = f.readlines()

            

    subprocess.run(f"{ipglasma_cmd} {temp_input_path} >> {datadir}/ipglasma_log", shell=True)

    remove_temp_input(temp_input_path)


    id1 = 2 * s +1
    id2 = 2 * s + 2

    

    print("===== RUNNING SUBNUCLEONDIFFRACTION FOR EVENTS  ", id1, " AND ", id2, " =====")
    subprocess.run(f"{subnucleondiffraction_cmd} -dipole 1 ipglasma_binary Final_x_{xpom}_V-{id1} > {datadir}/spectra_{id1}", shell=True)

    subprocess.run(f"{subnucleondiffraction_cmd} -dipole 1 ipglasma_binary Final_x_{xpom}_V-{id2} > {datadir}/spectra_{id2}", shell=True)

    if id1 in plot_ids or id2 in plot_ids:
        plot_id = id1 if id1 in plot_ids else id2
        final_file = f"Final_x_{xpom}_V-{plot_id}"
        subprocess.run(f"{subnucleondiffraction_cmd} -dipole 1 ipglasma_binary {final_file} -print_nucleus -Q2 0 > tmp_output", shell=True)
        try:
            data = np.loadtxt("tmp_output")
            x = data[:, 0]
            y = data[:, 1]
            z = data[:, 5]

            plt.figure()
            plt.scatter(x, y, c=z, cmap="viridis", s=6, marker="s")
            plt.colorbar(label="col6")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title(fr"{final_file} $Tr V(x)/N_c$")
            plt.tight_layout()
            outpufile = f"{datadir}/output_{plot_id}.pdf"
            plt.savefig(outpufile)
            plt.close()
            print("Created plot of V(x) for event ", plot_id, " filename: ", outpufile)
        except Exception as e:
            print(f"Plotting failed for {final_file}: {e}")
        subprocess.run(f"rm -f tmp_output", shell=True)


    # update cross section
    if s % 2 == 0 and s >= 2 :
        PlotComparison(datadir=datadir)
        
    subprocess.run("rm -f Initial_x_*", shell=True)
    subprocess.run("rm -f JIMWLKSnapshot_x_*", shell=True)
    subprocess.run("rm -f Final_x_*", shell=True)