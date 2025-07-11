import os
from utility.utils import kill_ansys
from pyMAPDL.pySMART_LSection import SMARTLSection

def main():
    """The minimal example for a SMART LSection simulation.
    Class methods are easily adaptable to your needs.
    """

    # 1. Set up the result path
    root = os.path.dirname(os.path.abspath(__file__))
    result_path = os.path.join(root, 'results')

    # 2. Create the SMARTLSection instance
    simu = SMARTLSection(
        result_path=result_path,
        experiment_name='Test',
    )

    # 3. Create the folder structure and launch an MAPDL instance. You can also provide an existing MAPDL instance.
    simu.setup_experiment()

    # 4. Define the material and elements. Here we use typical properties of a high strength aluminum alloy of the 7000 series.
    simu.define_material_and_elements(70000, 0.33, 4E-11, 2.5)

    # 5. Build the L-Section geometry and insert a wedge crack (single-edge through crack).
    simu.build_base_geometry(50, 100, 300, 4, 90, 10)
    simu.insert_crack_geometry(10, 150)

    # 6. Mesh the geometry, el_size_crack sets the base for the crack growth increment size (in mm).
    # See ANSYS documentation for more details.
    simu.mesh_geometry(8, 0.5)

    # 7. Apply boundary conditions. Bottom is fixed, top has a force applied in the y-direction.
    simu.apply_bc_bottom()
    simu.apply_bc_top(F_y=10000)

    #Optional: visualize the geometry in MAPDL and save it to hard disk.
    simu.visualize_model()

    # 8. Initialize the crack (APDL: CINT) and set the growth parameters for SMART Crack Growth.
    # Let the crack grow until it reaches a length of 100 mm (+ 10 mm for the wedge).
    simu.initialize_crack()
    simu.initialize_crack_growth(R=0.1, cemx=100)

    # 9. Solve the simulation for n load steps
    # Remember that each load step corresponds to a crack growth increment.
    # Set it to a arbitrary high value to ensure the crack grows until your stopping criterion is reached.
    # When the crack reaches the defined stopping criterion (here a length of 100 mm), the simulation stops.
    simu.solve(1000)

    # 10. Export the results to the result path
    simu.export_results()

    # Optional steps: export the parameter file and clean up the solver directory to save disk space.
    simu.export_parameter_file()
    simu.cleanup()

if __name__ == "__main__":
    main()
    kill_ansys()
