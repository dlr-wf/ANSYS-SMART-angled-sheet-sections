import os
from utility.utils import kill_ansys
from pyMAPDL.pySMART_LSection import SMARTLSection
from ansys.mapdl.core import MapdlPool
from threading import Thread

def worker(mapdl, experiment_name: str, result_path: str,
                    thickness:float):
    """See 01_LSection_MinExample.py for the main function."""

    simu = SMARTLSection(
        mapdl_instance=mapdl,
        result_path=result_path,
        experiment_name=experiment_name,
    )
    simu.setup_experiment()
    simu.define_material_and_elements(70000, 0.33, 4E-11, 2.5)
    simu.build_base_geometry(25, 25, 200, thickness, 90, 10)
    simu.insert_crack_geometry(10, 150)
    simu.mesh_geometry(2*thickness, thickness/8)
    simu.apply_bc_bottom()
    simu.apply_bc_top(F_y=1000)
    simu.visualize_model()
    simu.initialize_crack()
    simu.initialize_crack_growth(R=0.1, cemx=100)
    simu.solve(5)
    simu.export_results()
    simu.export_parameter_file()

    print(f'Simulation {experiment_name} completed.')

def threaded_simulation(pool:MapdlPool, experiment_name: str, result_path: str,
                        thickness: float):
    """Run the worker function in parallel for each thickness."""
    with pool.next() as mapdl:
        worker(mapdl, experiment_name, result_path, thickness)

def main():
    """The minimal example for batch processing a SMART LSection simulation.
    We iterate over a set of parameters (e.g., thickness) and run the simulation for each parameter set.
    We make use of MapdlPool to run multiple simulations in parallel.
    """
    #1. Set up some prerequisites
    root = os.path.dirname(os.path.abspath(__file__))
    result_path = os.path.join(root, 'results')
    ansys_temp_path = os.path.join(result_path, 'ansys_temp')
    if not os.path.exists(result_path):
        os.makedirs(result_path)
    if not os.path.exists(ansys_temp_path):
        os.makedirs(ansys_temp_path)

    #2. Create your parameter sets
    template = lambda x: f'Thickness{x}mm'
    thicknesses = [2, 4, 6, 8, 10]

    #3. Create a MapdlPool instance to run simulations in parallel
    # Note: Adjust the number of instances and processors per according to your system capabilities and available license.
    pool = MapdlPool(4, nproc=4, run_location=ansys_temp_path,remove_temp_dir_on_exit=True)

    #4. Use Threading to run simulations in parallel
    threads = []
    for thickness in thicknesses:
        t = Thread(target=threaded_simulation, args=(pool, template(thickness), result_path, thickness))
        t.start()
        threads.append(t)

    for thread in threads:
        thread.join()

if __name__ == "__main__":
    main()
    kill_ansys()
