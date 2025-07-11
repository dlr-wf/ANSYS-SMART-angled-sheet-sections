import os
import json
import functools
from abc import ABC, abstractmethod
from typing import Optional, Tuple, List

from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.mapdl import MapdlBase as Mapdl
from utility.utils import clear_ansys_directory

PROCESSORS_PER_INSTANCE = 4

class SMARTSimulation(ABC):
    """Base class for ANSYS SMART based simulations using PyMAPDL."""

    def __init__(self, result_path: str, experiment_name: str|None = None,
                 mapdl_instance: Optional[Mapdl] = None,
                 log_to_console: bool = True) -> None:
        """
        Initialize the SmartSimulation class.

        :param result_path: string, path to the result directory where outputs will be stored.
        :param experiment_name: string, name of the experiment for identification. Creates a subdirectory in result_path if provided.
        :param mapdl_instance: Optional[Mapdl], existing MAPDL instance to use; if None, a new instance will be launched.
        :param log_to_console: bool, whether to log messages to the console.
        """

        self.experiment_name = experiment_name
        self.result_path = result_path
        self.mapdl: Optional[Mapdl] = mapdl_instance
        self.log_to_console = log_to_console
        self.local = None

        # helpers
        self.id_storage: dict = {} # stores IDs of entities created in MAPDL
        self.cm_name_storage: dict = {}  # stores APDL component names for the simulation
        self.parameter_storage: dict = {}  # stores parameters for the simulation

        # Simulation Status
        self.setup_finished = False
        self.material_defined = False
        self.geometry_built = False
        self.crack_inserted = False
        self.mesh_generated = False
        self.bc_bottom_applied = False
        self.bc_top_applied = False
        self.crack_initialized = False
        self.crack_growth_initialized = False
        self.SOLU_active = False
        self.simulation_solved = False
        self.data_exported = False

    def setup_experiment(self):
        """Set up the experiment by creating necessary directories and launching MAPDL."""

        self._setup_directories()
        self._launch_mapdl()

        #----------------------------------
        # Status update
        #----------------------------------
        self.setup_finished = True

    def _setup_directories(self) -> None:
        """Create necessary output folders and clean ANSYS directory."""

        if self.experiment_name:
            self.result_path = os.path.join(self.result_path, self.experiment_name)

        self.local = os.path.join(self.result_path, "local")

        for sub in ["local", "input", "nodemaps", "crackdata", "postprocessed", "misc"]:
            os.makedirs(os.path.join(self.result_path, sub), exist_ok=True)

        clear_ansys_directory(self.local)

    def _launch_mapdl(self) -> None:
        """Launch a MAPDL instance if none is provided."""
        if self.mapdl is None:
            print(f"Launching MAPDL instance on {PROCESSORS_PER_INSTANCE} cores...")
            self.mapdl = launch_mapdl(
                nproc=PROCESSORS_PER_INSTANCE,
                run_location=self.local,
                override=True,
                loglevel='INFO' if self.log_to_console else None,
                print_com=True,
                log_apdl=os.path.join(self.result_path, 'input', '00_main.inp'),
            )
        log_path = os.path.join(self.result_path, 'misc', 'debug.log')
        if os.path.isfile(log_path):
            os.remove(log_path)
        self.mapdl._log.log_to_file(log_path, "INFO")
        self.mapdl.clear()

    def cleanup(self) -> None:
        """Clean up the MAPDL instance and any temporary files created during the simulation."""

        if self.mapdl is not None:
            self.mapdl.exit()
            self.mapdl = None
        print("MAPDL instance cleaned up.")

        clear_ansys_directory(self.local)

    def export_parameter_file(self) -> None:
        """Export the parameter set to the input folder."""

        with open(os.path.join(self.result_path, 'input',
                               os.path.basename(self.result_path) + '.json'), 'w') as f:
            json.dump(self.parameter_storage, f, indent=4)

    # ------------------------------------------------------------------
    # model specific methods - must be overwritten by subclasses
    # ------------------------------------------------------------------
    @abstractmethod
    def define_material_and_elements(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def build_base_geometry(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def insert_crack_geometry(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def mesh_geometry(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def apply_bc_top(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def apply_bc_bottom(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def initialize_crack(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def initialize_crack_growth(self, *args, **kwargs) -> None:
        pass

    @abstractmethod
    def solve(self) -> None:
        pass

    @abstractmethod
    def export_results(self) -> None:
        pass
