from typing import Optional
import numpy as np
import os

from ansys.mapdl.core.mapdl import MapdlBase as Mapdl

from pyMAPDL.simulation import SMARTSimulation
from utility.utils import export_nodemaps, export_fracture_parameters

class SMARTLSection(SMARTSimulation):
    """Class for SMART L-Section simulations using PyMAPDL."""

    def __init__(self, result_path: str, experiment_name: str = "SMART_LSection",
                 mapdl_instance: Optional[Mapdl] = None,
                 log_to_console: bool = True) -> None:
        """
        Initialize the SMARTLSection class.
        :param result_path: string, path to the result directory where outputs will be stored.
        :param experiment_name: string, name of the experiment for identification.
        :param mapdl_instance: Optional[Mapdl], existing MAPDL instance to use; if None, a new instance will be launched.
        :param parameters: Optional[ParameterSet], full parameters for the simulation, if any. Can also be exposed via class methods.
        :param log_to_console: bool, whether to log messages to the console.
        """
        super().__init__(result_path, experiment_name, mapdl_instance, log_to_console)

        self.cm_name_storage.update({
            'CrackFront': 'CFR_NODES',  # crack front component name
            'CrackFaceLower': 'CFA_NODES_LOWER',  # lower crack face component name
            'CrackFaceUpper': 'CFA_NODES_UPPER'  # upper crack face component name
        })

    def define_material_and_elements(self, E:float, mu:float, C:float, m:float) -> None:  # to be implemented by subclasses
        """
        Define the material properties for the L-Section simulation.
        The material model is linear elastic and isotropic.
        Fatigue crack growth is defined by the Paris law.

        :param E: Young's modulus of the base material.
        :param mu: Poisson's ratio of the base material.
        :param C: Paris law coefficient for fatigue crack growth.
        :param m: Paris law exponent for fatigue crack growth.

        :return:
        """

        #-------------------------------
        # Preparation
        #-------------------------------
        if not self.setup_finished:
            raise RuntimeError("Setup must be finished before defining material. Call setup_experiment() first.")

        self.parameter_storage.update({
            'E': E,
            'mu': mu,
            'C': C,
            'm': m
        })

        # ------------------------------
        # Base Material Definition
        # ------------------------------

        self.mapdl.prep7()
        self.mapdl.allsel('ALL')

        # Defines the base material as linear elastic isotropic
        self.mapdl.mp(lab='EX',
                 mat=1,
                 c0=E
                 )
        self.mapdl.mp(lab='NUXY',
                 mat=1,
                 c0=mu
                 )

        # Defines the fatigue crack growth material model using the Paris law
        # adjust if necessary (See ANSYS documentation for details)
        self.mapdl.tb("cgcr", 1, 1, "", "PARIS")
        self.mapdl.tbdata(1, C, m)

        # SMART Crack Growth is restricted to SOLID187 elements
        # Define Elements for base
        self.mapdl.et(1,
                 'SOLID187'
                 )
        self.mapdl.type(1)

        # -----------------------------
        # Helper
        # -----------------------------

        # Defines stiff material for beam elements responsible for remote loading
        self.mapdl.mp(lab='EX',
                 mat=10,
                 c0=1e8
                 )
        self.mapdl.mp(lab='PRXY',
                 mat=10,
                 c0=0.33
                 )

        #------------------------------
        # States update
        #------------------------------

        self.material_defined = True

        pass

    def build_base_geometry(self, wl:float, wr:float, h:float, t:float, angle:float, radius:float) -> None:
        """ Build the geometry of the L-Section.
        Constructs a 3D L-Section with specified dimensions and properties.
        The legs are straight and the corner is filleted if a radius is specified.
        The total width of the L-Section is wl + wr + s(t, angle, r), where s is the width of the fillet at the corner.

        :param wl: Width of the left leg of the L-Section.
        :param wr: Width of the right leg of the L-Section.
        :param h: Height of the L-Section.
        :param t: Thickness of the L-Section. Homogeneous (sheet) thickness is assumed.
        :param angle: (Bending) angle of the L-Section in degrees. 0 degrees means a flat section. Only positive angles are allowed.
        :param radius: (Inner) Radius in of the fillet at the corner of the L-Section. 0 means a sharp corner without fillet.
        """

        #-------------------------------
        # Preparation
        #-------------------------------
        if not self.setup_finished:
            raise RuntimeError("Setup must be finished before building geometry. Call setup_experiment() first.")

        # dependent parameters
        if angle == 0:
            angle, radius = 0, 0

        self.parameter_storage.update({
            'wl': wl,
            'wr': wr,
            'h': h,
            't': t,
            'angle': angle,
            'radius': radius
        })

        #------------------------------
        # L-Section Geometry Definition
        #------------------------------

        self.mapdl.prep7()
        self.mapdl.seltol(0.0001)

        self.mapdl.allsel('ALL')

        # suffix defining keypoints relative to the local coordinate system
        # _lu = left upper corner, _ru = right upper corner,
        # _ld = left lower corner, _rd = right lower corner
        # _inner = inner radius, _outer = outer radius

        # Keypoints - left of bending edge
        kp1_ld = self.mapdl.k("", 0, 0, 0)
        kp1_lu = self.mapdl.k("", 0, 0, t)
        kp1_ru = self.mapdl.k("", wl, 0, t)
        kp1_rd = self.mapdl.k("", wl, 0, 0)

        # define offsets for the keypoints
        dx, dz = t, t
        # take care of uniform thickness when bending radius == 0
        if not radius:
            dx = np.tan(np.deg2rad(angle / 2)) * t
            kp1_rd = self.mapdl.k("", wl + dx, 0, 0)
        # when there is no bend_angle, there's no offset in z direction

        # lines -  left of bending angle
        li1_l = self.mapdl.l(kp1_lu, kp1_ld)
        li1_u = self.mapdl.l(kp1_lu, kp1_ru)
        # connecting line - right side
        li_conn = self.mapdl.l(kp1_ru, kp1_rd)
        li1_d = self.mapdl.l(kp1_rd, kp1_ld)

        # create area left of radius
        self.mapdl.al(li1_l, li1_u, li_conn, li1_d)

        # Keypoints - right of bending edge -> ar_rect_02
        # helper coordinate system for the bending edge
        self.mapdl.local(50, 0,
                    xc=wl + radius * np.sin(np.deg2rad(angle)),
                    zc=radius * (1 - np.cos(np.deg2rad(angle))) + dz,
                    thzx=-angle)
        kp2_lu = kp1_ru
        kp2_ld = kp1_rd
        kp2_ru = self.mapdl.k("", wr, 0, 0)
        kp2_rd = self.mapdl.k("", wr, 0, -t)

        # Keypoints + Lines - rounded bending edge
        # when bend_radius == 0, sharp edge
        if angle and radius:
            # new cylindrical helper coordinate system
            self.mapdl.local(51, 1,
                        xc=wl, zc=radius + t,
                        thxy=-90,
                        thzx=90)
            # keypoints for the radii
            kp2_lu = self.mapdl.k("", radius, angle, 0)
            kp2_ld = self.mapdl.k("", radius + t, angle, 0)
            # lines in CYLINDER coordinate system can only sweeped 180° in APDL
            # therefore additional keypoints have to be created when bending angle > 180
            if angle > 180:
                # additional keypoints for the first 180°
                kp_180_inner = self.mapdl.k("", radius, 180, 0)
                kp_180_outer = self.mapdl.k("", radius + t, 180, 0)
                # lines for the first 180°
                l_180_inner = self.mapdl.l(kp1_ru, kp_180_inner)
                l_180_outer = self.mapdl.l(kp1_rd, kp_180_outer)
                l_180_left = li_conn
                l_180_right = self.mapdl.l(kp_180_outer, kp_180_inner)
                li_conn = self.mapdl.l(kp_180_inner, kp_180_outer)
                # create area for the first 180°
                self.mapdl.al(l_180_left, l_180_outer, l_180_right, l_180_inner)
                # redefining so that the lines for (bending angle-180°) can be sweeped
                kp1_ru = kp_180_inner
                kp1_rd = kp_180_outer
            # lines for bending edge
            l_b_inner = self.mapdl.l(kp1_ru, kp2_lu)
            l_b_outer = self.mapdl.l(kp1_rd, kp2_ld)
            # take account for the aditional line for the area right of the bending edge
            self.mapdl.csys(0)
            l_b_left = li_conn
            li_conn = self.mapdl.l(kp2_ld, kp2_lu)
            # area for the bending edge
            self.mapdl.al(l_b_left, l_b_inner, li_conn, l_b_outer)

        # lines - right of bending edge
        li2_u = self.mapdl.l(kp2_lu, kp2_ru)
        li2_r = self.mapdl.l(kp2_ru, kp2_rd)
        li2_d = self.mapdl.l(kp2_rd, kp2_ld)

        # area - right of bending edge
        self.mapdl.al(li_conn, li2_u, li2_r, li2_d)

        # Extrude to create the full bended sheet
        self.mapdl.asel("S", vmin="ALL")
        self.mapdl.vext('ALL', dy=h)

        # glue volumes together, change selection tolerance for more consistent results
        self.mapdl.btol(10E-3)
        self.mapdl.vsel("S", vmin='ALL')
        self.mapdl.vglue("ALL")
        self.mapdl.btol("DEFA")

        # -------------------------------
        # Status update
        # -------------------------------

        self.geometry_built = True

    def insert_crack_geometry(self, a_init:float, cr_pos_h:float) -> None:
        """Insert the crack geometry into the L-Section model.

        The crack is:
        - defined as a wedge-shaped volume that is subtracted from the L-Section.
        - positioned at a specified height at the left edge of the left leg of the L-Section.
        - oriented such that it is normal to the surface of the L-Section.

        The result is a single-edge through crack that extends through the thickness of the L-Section.

        :param a_init: Initial crack length in mm.
        :param cr_pos_h: Position of the crack in the height of the L-Section in mm.
        """

        #-------------------------------
        # Preparation
        #-------------------------------

        # some additional parameters for the crack geometry
        dh = 0.5  # height of the wedge crack opening

        self.parameter_storage.update({
            'a_init': a_init,
            'cr_pos_h': cr_pos_h,
        })

        #------------------------------
        # Crack Geometry Definition
        #------------------------------

        self.mapdl.prep7()
        self.mapdl.seltol(0.0001)
        self.mapdl.csys(0)

        t = self.parameter_storage['t']

        # Define initial Crack growth direction by CS
        self.mapdl.local(30,  # CS ID
                         0,  # cartesian coordinate system
                         0,  # x coordinate of the origin
                         cr_pos_h,  # y coordinate of the origin
                         t / 2,  # z coordinate of the origin
                         0,  # rotation around z-axis
                         0,  # rotation around x-axis
                         0  # rotation around y-axis
                         )

        self.mapdl.allsel('ALL')

        # keypoints crack front
        kp1_cfr_l = self.mapdl.k('', a_init, 0, -t / 2)
        kp1_cfr_r = self.mapdl.k('', a_init, 0, t / 2)

        # additional keypoints for crack opening
        # keypoints for the crack opening at the backface of the wedge crack
        kp1_cfa_lu = self.mapdl.k('', 0, dh, -t / 2)
        kp1_cfa_ld = self.mapdl.k('', 0, -dh, -t / 2)
        kp1_cfa_ru = self.mapdl.k('', 0, dh, t / 2)
        kp1_cfa_rd = self.mapdl.k('', 0, -dh, t / 2)

        # Create wedge
        # linescreate crack front
        li1_cfr = self.mapdl.l(kp1_cfr_l, kp1_cfr_r)
        # select normal crack CS

        # lines create crack opening
        # side left
        li1_cfa_l_u = self.mapdl.l(kp1_cfa_lu, kp1_cfr_l)
        li1_cfa_l_d = self.mapdl.l(kp1_cfa_ld, kp1_cfr_l)
        li1_cfa_l_ud = self.mapdl.l(kp1_cfa_lu, kp1_cfa_ld)

        # side right
        li1_cfa_r_u = self.mapdl.l(kp1_cfa_ru, kp1_cfr_r)
        li1_cfa_r_d = self.mapdl.l(kp1_cfa_rd, kp1_cfr_r)
        li1_cfa_r_ud = self.mapdl.l(kp1_cfa_ru, kp1_cfa_rd)

        # backface additonal lines
        li1_cfa_u_lr = self.mapdl.l(kp1_cfa_lu, kp1_cfa_ru)
        li1_cfa_d_lr = self.mapdl.l(kp1_cfa_ld, kp1_cfa_rd)

        # upper crack face
        ar1_cfa_u = self.mapdl.al(li1_cfa_l_u, li1_cfr, li1_cfa_r_u, li1_cfa_u_lr)

        # lower crack face
        ar1_cfa_d = self.mapdl.al(li1_cfa_l_d, li1_cfr, li1_cfa_r_d, li1_cfa_d_lr)

        # left triangle, right triangle an backface
        ar1_cfa_ltri = self.mapdl.al(li1_cfa_l_u, li1_cfa_l_d, li1_cfa_l_ud)
        ar1_cfa_rtri = self.mapdl.al(li1_cfa_r_u, li1_cfa_r_d, li1_cfa_r_ud)
        ar1_cfa_back = self.mapdl.al(li1_cfa_u_lr, li1_cfa_r_ud, li1_cfa_d_lr, li1_cfa_l_ud)

        # crack volume
        base_vol_id = self.mapdl.geometry.vnum[0]
        crack_vol_id = self.mapdl.va(ar1_cfa_u, ar1_cfa_d, ar1_cfa_ltri, ar1_cfa_rtri, ar1_cfa_back)

        # subtract the volumes to create a crack
        self.mapdl.vsbv(base_vol_id, crack_vol_id)

        # store important crack geometry IDs for later use
        self.id_storage.update({
            'li_cfr': li1_cfr,  # crack front line ID for meshing
            'ar_cfa_lower': ar1_cfa_d,  # lower crack face area ID for meshing
            'ar_cfa_upper': ar1_cfa_u  # upper crack face area ID for meshing
        })

        #-------------------------------
        # Status update
        #-------------------------------

        self.crack_inserted = True

    def mesh_geometry(self, el_size_base:float, el_size_crack:float) -> None:
        """Mesh the geometry of the L-Section.

        This method applies a mesh to the L-Section geometry using SOLID187 elements.

        :param el_size_base: Element size for the base L-Section geometry. Describes edge length of the elements in mm.
        :param el_size_crack: Element size for the crack geometry. Describes edge length of the elements in mm.

        """

        #-------------------------------
        # Preparation
        #-------------------------------
        if not self.geometry_built:
            raise RuntimeError("Geometry must be built before meshing. Call build_base_geometry() first.")

        self.parameter_storage.update({
            'el_size_base': el_size_base,
            'el_size_crack': el_size_crack
        })

        h = self.parameter_storage['h']  # height of the L-Section
        wr = self.parameter_storage['wr']  # width of the right leg of the L-Section

        cm_cfr_nodes = self.cm_name_storage['CrackFront']  # crack front component name
        cm_cfa_lower_nodes = self.cm_name_storage['CrackFaceLower']
        cm_cfa_upper_nodes = self.cm_name_storage['CrackFaceUpper']

        #------------------------------
        # Crack Geometry Meshing
        #------------------------------

        self.mapdl.prep7()
        self.mapdl.csys(0)
        self.mapdl.type(1)  # set element type to SOLID187 of the base material

        if self.crack_inserted:
            # crackfront
            self.mapdl.csys(30)
            self.mapdl.lsel(type_='S', vmin=self.id_storage['li_cfr'])
            self.mapdl.cm('CrackFront', 'LINE')
            self.mapdl.lesize(nl1='ALL',
                         size=el_size_crack,
                         )
            self.mapdl.allsel('ALL')

            # lower crack face
            self.mapdl.asel(type_='S',
                       vmin=self.id_storage['ar_cfa_lower']
                       )
            self.mapdl.cm('CrackFaceLower', 'AREA')
            self.mapdl.aesize(anum='ALL',
                         size=el_size_crack,
                         )
            self.mapdl.allsel('ALL')

            # upper crack face
            self.mapdl.asel(type_='S',
                       vmin=self.id_storage['ar_cfa_upper']
                       )
            self.mapdl.cm('CrackFaceUpper', 'AREA')
            self.mapdl.aesize(anum='ALL',
                         size=el_size_crack,
                         )
            self.mapdl.allsel('ALL')

        #------------------------------
        # Base Geometry Meshing
        #------------------------------

        # mesh the rest of the sample
        self.mapdl.esize(size=el_size_base)
        self.mapdl.mopt('TIMP', 4)
        self.mapdl.vmesh('All')
        self.mapdl.allsel('ALL')

        #-------------------------------
        # Mesh Refinement
        #-------------------------------
        self.mapdl.csys(0)

        # refine the mesh at the bottom and top where the boundary conditions are applied
        self.mapdl.allsel('ALL')
        self.mapdl.nsel('S', 'LOC', 'Y', 0)
        self.mapdl.nsel('A', 'LOC', 'Y', h)
        self.mapdl.nrefine('ALL',
                      level=1,
                      depth=1)

        # refine mesh at the backside of the sample for application of the boundary conditions
        self.mapdl.csys(50)
        self.mapdl.nsel('S', "LOC", "X", wr)
        self.mapdl.nsel('R', "LOC", "Y", h / 2 * 0.95, h / 2 * 1.05)
        self.mapdl.esln('S')
        self.mapdl.erefine('ALL',
                      level=3,
                      depth=1
                      )
        self.mapdl.csys(0)
        self.mapdl.allsel('ALL')

        # refine crack front mesh -> Esize defining parameter for crack growth increment
        self.mapdl.cmsel('S', 'CrackFront')
        self.mapdl.lrefine('All',
                      level=2,
                      depth=3
                      )

        self.mapdl.allsel('ALL')

        self.mapdl.csys(0)
        self.mapdl.allsel('ALL')

        #-------------------------------
        # Crack Front Nodes Components
        #-------------------------------

        if self.crack_inserted:
            # define crackfront nodes for reference
            self.mapdl.cmsel('S',
                             name='CrackFront')
            self.mapdl.nsll('S', 1)
            self.mapdl.cm(cm_cfr_nodes, 'NODE')

            # define crack face nodes for later reference
            self.mapdl.allsel('ALL')
            self.mapdl.csys(30)  # probably not needed
            self.mapdl.cmsel(type_='S',
                             name='CrackFaceLower')
            self.mapdl.nsla('S', 1)
            self.mapdl.cm(cm_cfa_lower_nodes, 'NODE')

            self.mapdl.allsel('ALL')
            self.mapdl.cmsel(type_='S',
                             name='CrackFaceUpper')
            self.mapdl.nsla('S', 1)
            self.mapdl.cm(cm_cfa_upper_nodes, 'NODE')

            # select assistant node for better CG direction initialization
            self.id_storage.update({'cfr_assistant_node': self.mapdl.mesh.nnum[0]})

            self.mapdl.csys(0)
            self.mapdl.allsel('ALL')

        #-------------------------------
        # Status update
        #-------------------------------

        self.mesh_generated = True

        pass

    def apply_bc_bottom(self) -> None:
        """
        Apply the bottom boundary condition to the L-Section model.
        This simplified method applies a fixed support at the bottom of the L-Section.

        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not self.mesh_generated:
            raise RuntimeError("Mesh must be generated before applying boundary conditions. "
                               "Call mesh_geometry() first.")

        #------------------------------
        # Bottom Boundary Condition
        #------------------------------

        self.mapdl.prep7()
        self.mapdl.allsel('ALL')

        # select nodes at the bottom of the L-Section
        self.mapdl.nsel('S', 'LOC', 'Y', 0)

        # apply fixed support at the bottom nodes
        self.mapdl.d('ALL', 'ALL')

        #------------------------------
        # Status update
        #------------------------------

        self.bc_bottom_applied = True

    def apply_bc_top(self,
                     F_x: float = 0,
                        F_y: float = 0,
                        F_z: float = 0,
                        M_x: float = 0,
                        M_y: float = 0,
                        M_z: float = 0,
                        U_x: float|None = 0,
                        U_y: float|None = 0,
                        U_z: float|None = 0,
                        R_x: float|None = 0,
                        R_y: float|None = 0,
                        R_z: float|None = 0,
                     )-> None:
        """ Apply the top boundary condition to the L-Section model.

        This method applies a remote force and moment at the top of the L-Section.
        The force and moment are applied in the local coordinate system defined by the center of gravity of the uncracked section.
        The loads are applied to a remote node that is defined at the top of the L-Section in that center of gravity.
        The loads are transferred via stiff beam elements to the top of the L-Section.

        :param F_x: Force in x-direction (N).
        :param F_y: Force in y-direction (N).
        :param F_z: Force in z-direction (N).
        :param M_x: Moment around x-axis (Nm).
        :param M_y: Moment around y-axis (Nm).
        :param M_z: Moment around z-axis (Nm).
        :param U_x: Displacement in x-direction (m). If None, no displacement is applied.
        :param U_y: Displacement in y-direction (m). If None, no displacement is applied.
        :param U_z: Displacement in z-direction (m). If None, no displacement is applied.
        :param R_x: Rotation around x-axis (rad). If None, no rotation is applied.
        :param R_y: Rotation around y-axis (rad). If None, no rotation is applied.
        :param R_z: Rotation around z-axis (rad). If None, no rotation is applied.

        Attention:
        F and M settings will overwrite settings for U and R.
        -> e.g. if F_y != 0 and U_y != 0, only F_y will be applied.

        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not self.mesh_generated:
            raise RuntimeError("Mesh must be generated before applying boundary conditions. "
                               "Call mesh_geometry() first.")

        h = self.parameter_storage['h']  # height of the L-Section

        #------------------------------
        # Remote Node Definition
        #------------------------------

        self.mapdl.prep7()
        self.mapdl.allsel('ALL')

        self.mapdl.asel('S', 'LOC', 'Y', h)  # select nodes at the top of the L-Section
        self.mapdl.asum()  # necessary, Calculates and prints geometry statistics of the selected areas in APDL
        xc, yc, zc = [self.mapdl.get(entity="AREA", item1="cent", it1num=axis) for axis in ["X", "Y", "Z"]]
        print("Calculating force equivalents for Face Sigma...")
        uncracked_crosssection = self.mapdl.get(entity="AREA", item1="AREA")

        # setup from ANSYS Workbenchs definition of remode loading
        self.mapdl.starset("tid", 10)
        self.mapdl.starset("cid", 11)
        self.mapdl.et("cid", '174')
        self.mapdl.et("tid", '170')
        self.mapdl.keyopt("tid", 2, 1)
        self.mapdl.keyopt("tid", 4, 111111)
        self.mapdl.keyopt("cid", 12, 5)
        self.mapdl.keyopt("cid", 4, 1)  # 1 = Deformable, 2 = Rigid, 3 = Coupled
        self.mapdl.keyopt("cid", 2, 2)

        # create a pilot node at the center of gravity of the uncracked section
        pilot_node = self.mapdl.n("", xc, h, zc)

        # create the beam elements that will transfer the loads to the top of the L-Section
        self.mapdl.type("tid")
        self.mapdl.mat("cid")
        self.mapdl.real("cid")
        self.mapdl.tshap("pilo")
        pilot_element = self.mapdl.mesh.enum[-1] + 1
        self.mapdl.en(pilot_element, pilot_node)
        self.mapdl.tshap("")

        nodes_top = self.mapdl.nsel("S", "LOC", "Y", h)
        nodes_top = [node for node in nodes_top if node != pilot_node]
        for node in nodes_top:
            self.mapdl.e(node, pilot_node)

        #------------------------------
        # Remote Loading Definition
        #------------------------------

        self.mapdl.f(pilot_node, "FX", F_x)
        self.mapdl.f(pilot_node, "FY", F_y)
        self.mapdl.f(pilot_node, "FZ", F_z)
        self.mapdl.f(pilot_node, "MX", M_x)
        self.mapdl.f(pilot_node, "MY", M_y)
        self.mapdl.f(pilot_node, "MZ", M_z)

        if F_x == 0 and U_x is not None:
            self.mapdl.d(pilot_node, "UX", U_x)
        if F_y == 0 and U_y is not None:
            self.mapdl.d(pilot_node, "UY", U_y)
        if F_z == 0 and U_z is not None:
            self.mapdl.d(pilot_node, "UZ", U_z)
        if M_x == 0 and R_x is not None:
            self.mapdl.d(pilot_node, "ROTX", R_x)
        if M_y == 0 and R_y is not None:
            self.mapdl.d(pilot_node, "ROTY", R_y)
        if M_z == 0 and R_z is not None:
            self.mapdl.d(pilot_node, "ROTZ", R_z)

        self.mapdl.type(1)
        self.mapdl.csys(0)
        self.mapdl.allsel('ALL')

        # update the parameters for the remote loading
        self.parameter_storage.update({
            'F_x': F_x,  # force in x-direction
            'F_y': F_y,  # force in y-direction
            'F_z': F_z,  # force in z-direction
            'M_x': M_x,  # moment around x-axis
            'M_y': M_y,  # moment around y-axis
            'M_z': M_z,   # moment around z-axis
            'xc': xc,  # x-coordinate of the center of gravity of the uncracked section
            'yc': yc,  # y-coordinate of the center of gravity of the uncracked section
            'zc': zc,  # z-coordinate of the center of gravity of the uncracked section
            'uncracked_crosssection': uncracked_crosssection,  # area of the uncracked cross-section
        })

        self.id_storage.update({
            'pilot_node': pilot_node,  # ID of the pilot node for remote loading
        })

        #------------------------------
        # Status update
        #------------------------------

        self.bc_top_applied = True
        pass

    def visualize_model(self) -> None:
        """Visualize the L-Section model utilizing pyMAPDL's VTK visualization capabilities.

        The mesh is displayed and saved to the postprocessing folder.
        """

        if not self.mesh_generated:
            raise RuntimeError("Mesh must be generated before visualization. "
                               "Call mesh_geometry() first.")

        # Visualization of the model
        for cpos in ('xy', 'yx', 'xz', 'zx', 'yz', 'zy', 'iso'):
            # continue
            self.mapdl.eplot(plot_bc=True,
                        plot_bc_legend=True,
                        background='w',
                        smooth_shading=True,
                        bc_glyph_size=15,
                        line_width=3,
                        cpos=cpos,
                        title='',
                        window_size=[1024, 768],
                        savefig=os.path.join(self.result_path, 'postprocessed', f'LSection_{cpos}.png'),
                        )
            pass

    def initialize_crack(self, eval_type: str = 'SIFS', contours: int = 6) -> None:
        """
        Initialize the crack in the L-Section model.

        - Defines the crack front components
        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not (self.crack_inserted and self.mesh_generated):
            raise RuntimeError("Crack geometry must be inserted and meshed before initialization. "
                               "Call insert_crack_geometry() and mesh_geometry() first.")

        self.parameter_storage.update({
            'eval_type': eval_type,  # evaluation type for the crack front. Default is 'SIFS'. See ANSYS documentation for details.
            'contours': contours  # number of contours for the crack front evaluation. See ANSYS documentation for details.
        })

        cm_cfr_nodes = self.cm_name_storage['CrackFront']  # crack front component name
        cm_cfa_lower_nodes = self.cm_name_storage['CrackFaceLower']  # lower crack face component name
        cm_cfa_upper_nodes = self.cm_name_storage['CrackFaceUpper']  # upper crack face component name
        node_cfr_assistant = self.id_storage['cfr_assistant_node']  # assistant node for crack front initialization

        #------------------------------
        # Crack Front Initialization
        #------------------------------
        # checks if already in solution routine
        if self.mapdl.parameters.routine != 'SOLUTION':
            self.mapdl.slashsolu()

        # define the basic fracture parameter calculation, see ANSYS Fracture Mechanics documentation
        self.mapdl.cint("new", 1)
        self.mapdl.cint("type", eval_type, 0) # 0 -> default, 1 -> plane stress, 2 -> plane strain
        self.mapdl.cint("ctnc", cm_cfr_nodes)
        self.mapdl.cint('edir', 'CS', 30, 'X', "", node_cfr_assistant)
        self.mapdl.cint("norm", 30, 2)
        self.mapdl.cint("surf", cm_cfa_upper_nodes, cm_cfa_lower_nodes)
        self.mapdl.cint("ncon", contours)

        #------------------------------
        # Status update
        #------------------------------

        self.crack_initialized = True
        pass

    def initialize_crack_growth(self, R: float, cemx:float, da_min: float = 0, da_max:float = 4) -> None:  # to be implemented by subclasses
        """Initialize the crack growth algorithm in the L-Section model.

        This method sets up the crack growth algorithm using the SMART framework.
        For details on the SMART framework, refer to the documentation.

        :param R: Load ratio for the crack growth. This is only used to define the stress intensity factor range inside the SMART Algorithm. See documentation for details.
        :param cemx: Maximum crack length. Simulation will stop when the crack reaches this length. Usually more useful definition than substeps.
        :param da_min: Minimum crack growth increment. This is the smallest increment a crack can grow in one step. Default is 0. Limited by 0.5*the element size of the crack front.
        :param da_max: Maximum crack growth increment. This is the largest increment a crack can grow in one step. Default is 4. Limited by 2*the element size of the crack front.
        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not self.crack_initialized:
            raise RuntimeError("Crack must be initialized before crack growth. "
                               "Call initialize_crack() first.")

        # some predefined parameters for the SMART framework
        method = 'LC'
        smart_eval_contour = 2
        crgrowth_increment = self.parameter_storage['el_size_crack']  # element size for the crack growth increment

        self.parameter_storage.update({
            'R': R,  # stress intensity factor range for the crack growth
            'cemx': cemx,  # maximum crack growth increment
            'da_min': da_min,  # minimum crack growth increment
            'da_max': da_max,  # maximum crack growth increment
            'method': method,  # SMART method for crack growth
            'smart_eval_contour': smart_eval_contour,  # contour for SMART evaluation
            'crgrowth_increment': crgrowth_increment  # crack growth increment size
        })

        #------------------------------
        # Crack Growth Initialization
        #------------------------------
        if self.mapdl.parameters.routine != 'SOLUTION':
            self.mapdl.slashsolu()

        self.mapdl.cgrow("new", 1)
        self.mapdl.cgrow("cid", 1)
        self.mapdl.cgrow("method", "smart")

        self.mapdl.run(f"CGROW,FCOPTION,mtab,1,{smart_eval_contour}")
        self.mapdl.cgrow("fcg", "meth", method)
        self.mapdl.cgrow("fcg", "damn", da_min)
        self.mapdl.cgrow("fcg", "damx", da_max)
        self.mapdl.cgrow("fcg", "srat", R)

        self.mapdl.cgrow("STOP", "cemx", cemx)
        # mapdl.cgrow("STOP", "FBOU") -> can be used to stop the simulation when the crack reaches a free boundary

        #------------------------------
        # Status update
        #------------------------------

        self.crack_growth_initialized = True
        pass

    def solve(self, substeps: int = 1000) -> None:
        """Solve the L-Section simulation.

        This method runs the simulation for the specified number of substeps.
        Each substeps represents on crack growth increment.
        The simulation will stop when the crack reaches stop criteria in the parameters or when the maximum number of substeps is reached.

        :param substeps: Number of substeps to run the simulation. Default is 1000.
        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not self.mesh_generated:
            raise RuntimeError("Mesh must be generated before solving. "
                               "Call mesh_geometry() first.")

        #-------------------------------
        # Setup for solving
        #-------------------------------
        if self.mapdl.parameters.routine != 'SOLUTION':
            self.mapdl.slashsolu()

        self.mapdl.antype("static")
        self.mapdl.eqslv('sparse')
        self.mapdl.rescontrol("", "NONE")
        self.mapdl.nsubst(nsbstp=substeps, nsbmn=substeps)  # number of substeps
        self.mapdl.outres('All', 'All')
        self.mapdl.kbc(1)

        self.mapdl.solve(verbose=True if self.log_to_console else None)

        self.mapdl.finish()

        #-------------------------------
        # Status update
        #-------------------------------

        self.simulation_solved = True

    def export_results(self, export_nodemap_data: bool = True, export_cf_data:bool = True, header:bool = True, cleanup_before_export = True) -> None:
        """Export the results of the L-Section simulation.

        This method exports the results of the simulation to the specified result path.
        The results include crack growth data, stress intensity factors, and other relevant data.

        :param export_nodemap_data: bool, whether to export nodemap data. Default is True.
        :param export_cf_data: bool, whether to export crack front data. Default is True.
        :param header: bool, whether to include a header in the exported data. Default is True.
        :param cleanup_before_export: bool, whether to clean up the result path before exporting. Default is True.
        """

        #-------------------------------
        # Preparation
        #-------------------------------

        if not self.simulation_solved:
            raise RuntimeError("Simulation must be solved before exporting results. "
                               "Call solve() first.")

        nm_folder = self.result_path + '/nodemaps'
        cf_folder = self.result_path + '/crackdata'

        if export_cf_data:
            xc, cr_pos_h, zc = self.parameter_storage['xc'], self.parameter_storage['cr_pos_h'], self.parameter_storage[
                'zc']
            cm_cfr_nodes = self.cm_name_storage['CrackFront']  # crack front component name
            contours = self.parameter_storage['contours']  # number of contours for the crack front evaluation

        #-------------------------------
        # Export Nodemap Data
        #-------------------------------

        if export_nodemap_data:
            try:
                export_nodemaps(self.mapdl, nm_folder, header=True,
                                cleanup_before_export=cleanup_before_export)
            except Exception as e:
                print(f"Error exporting nodemap data: {e}")

        #-------------------------------
        # Export Crack Front Data
        #-------------------------------

        if export_cf_data:
            try:
                export_fracture_parameters(self.mapdl,
                                           cm_cfr_nodes,
                                           cf_folder,
                                           n_contours=contours,
                                           header=True,
                                           cleanup_before_export=cleanup_before_export,
                                           ref_point_cf_ordering=(xc, cr_pos_h, zc))
            except Exception as e:
                print(f"Error exporting crack front data: {e}")