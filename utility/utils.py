# contains some helper functions
import os
import sys
import getpass

import pandas as pd
import numpy as np
import subprocess

from datetime import datetime
from typing import List, Tuple, Optional, Any
from ansys.mapdl.core.mapdl import MapdlBase as Mapdl

from utility.pattern import ExportDataColumns


def distance_point_from_line(p0: np.ndarray, p1: np.ndarray, p_dist: np.ndarray) -> float:
    """Calculates the distance of a point from a line in 3 dimensions.

        :param p0: Point 1 of the line.
        :param p1: Point 2 of the line.
        :param p_dist: Point the distance to should be measured.

        :return: Distance of p_dist from the line defined by p0 and p1.
    """
    return np.linalg.norm(np.cross(p1 - p0, p0 - p_dist)) / np.linalg.norm(p1 - p0)


def clear_ansys_directory(ansys_folder: str):
    """Clears all solver output files from the ansys_folder.
    This is useful for cleaning up the ansys_folder before a new run.

    :param ansys_folder: Path to the folder containing the solver output files.
    """

    files = os.listdir(ansys_folder)
    delete_endings = ('.out', '_.inp', '.db', '.DSP', '.esav', '.mntr', '.iges', '.page', '.sh',
                      '.lock', '.rst', '.err', '.esav', '.full', '.stat', '.log', 'anstmp', '.txt')
    for ansys_file in files:
        if ansys_file.endswith(delete_endings):
            os.remove(os.path.join(ansys_folder, ansys_file))

def kill_ansys():
    """Kill all running ANSYS processes for the current user."""
    print('Killing all ANSYS processes for current user...')

    # exclusion list
    excluded_str = []
    # identifier for python main process
    main_script_str = os.path.join(os.getcwd(), os.path.basename(sys.argv[0]))
    excluded_str.append(main_script_str)
    # exlude open folders
    excluded_str.append('nautilus')


    p = subprocess.Popen(['ps', '-eo', 'pid,user,command'], stdout=subprocess.PIPE)
    out, err = p.communicate()
    for line in out.decode().splitlines():
        # exclud exceptions
        if any([excluded in line for excluded in excluded_str]):
            continue
        if 'ansys' in line.lower():
            pid, user, cmd = line.split(None, 2)
            if user == getpass.getuser():
                try:
                    os.kill(int(pid), 9)
                except ProcessLookupError:
                    pass
                print(f'Killed process {pid} for user {user} with command {cmd}')

def export_nodemaps(mapdl: Mapdl, output_path: str, filename_prefix: str = 'FE_nodemap_Full', component = None,
                    header=True, cleanup_before_export=True) \
        -> List[pd.DataFrame]:
    """Exports nodemaps in a FAT-readable format for all steps in a result set.
       Automatically appends LS_{i}SS_{j} to the filename.

        :param mapdl: MAPDL instance
        :param output_path: Path to export the nodemaps to
        :param filename_prefix: Prefix for the exported nodemap files
        :param component: Component str to export. If None, all nodes will be exported.
        :param header: If True, the first line of the exported file will contain the column names
        :param cleanup_before_export: If True, all files in the output_path will be deleted before export

        :return: List of pandas DataFrames containing the nodemaps
    """
    if cleanup_before_export:
        # deletes old nodemap files in the folder to avoid confusion
        filelist = [f for f in os.listdir(output_path) if f.endswith(".txt")]
        for f in filelist:
            os.remove(os.path.join(output_path, f))

    mapdl.post1()
    n_sets = mapdl.post_processing.nsets
    lst_results = []
    COLUMNS = ExportDataColumns().nodemap

    # extract nodal information for all steps
    for result_nbr in range(n_sets):
        mapdl.set(nset=result_nbr + 1)
        mapdl.allsel('All')

        if component is not None:
            mapdl.cmsel('S',component,"NODE")

        # extract step info
        n_load_step = mapdl.post_processing.load_step
        n_sub_step = mapdl.post_processing.sub_step

        # create new table with metadata
        df = pd.DataFrame(columns=COLUMNS)
        df.attrs['Loadstep'] = n_load_step
        df.attrs['Substep'] = n_sub_step

        # fill table with relevant information
        try:
            df['x'] = mapdl.post_processing.nodal_values('LOC', 'X')
            df['y'] = mapdl.post_processing.nodal_values('LOC', 'Y')
            df['z'] = mapdl.post_processing.nodal_values('LOC', 'Z')
            df['ux'] = mapdl.post_processing.nodal_displacement('X')
            df['uy'] = mapdl.post_processing.nodal_displacement('Y')
            df['uz'] = mapdl.post_processing.nodal_displacement('Z')
            df['epsx'] = mapdl.post_processing.nodal_total_component_strain('X') * 100
            df['epsy'] = mapdl.post_processing.nodal_total_component_strain('Y') * 100
            df['epsxy'] = mapdl.post_processing.nodal_total_component_strain('XY') * 0.5
            df['epseqv'] = mapdl.post_processing.nodal_total_eqv_strain() * 100
            df['sigx'] = mapdl.post_processing.nodal_component_stress('X')
            df['sigy'] = mapdl.post_processing.nodal_component_stress('Y')
            df['sigxy'] = mapdl.post_processing.nodal_component_stress('XY')
            df['sigeqv'] = mapdl.post_processing.nodal_eqv_stress()
        except ValueError as e:
            print(f'Error while exporting nodemap for step {n_load_step} substep {n_sub_step}: {e}')
            continue

        # fill with the node number, sometimes the node number doesn't match the length of the other columns
        # could be because of th pilot node for applying the load
        try:
            df['ID'] = mapdl.post_processing.selected_nodes.nonzero()[0] + 1
        except ValueError as e:
            print(f'Error while exporting nodemap for step {n_load_step} substep {n_sub_step}: {e}')
            continue


        # export
        if os.path.isfile(output_path):
            filename_prefix = os.path.split(output_path)[1].split('.')[0]
        FILEENDSWITH = f'_LS{n_load_step:03d}SS{n_sub_step:03d}.txt'
        # put tabs between columns
        df.to_csv(os.path.join(output_path, filename_prefix + FILEENDSWITH), sep=';', index=False, header=header)
        lst_results.append(df)

    return lst_results

def export_fracture_parameters(mapdl: Mapdl, cm_name: str, output_path: str, crack_id=1,
                               n_contours: int = 6,
                               filename_prefix: str = 'FE_CrackData_FrctPrm',
                               header=True, cleanup_before_export=True,
                               crack_growth_defined=True,
                               ref_point_cf_ordering: tuple = (0, 0, 0),
                               legacy_node_ordering = False) -> List[
    pd.DataFrame]:
    """Exports fracture parameter (CINT) for all steps in a result set.
       Automatically appends LS_{i}SS_{j} to the filename.
        :param mapdl: An instance of MAPDL
        :param cm_name: variable name of the crack front component inside the MAPDL instance
        :param output_path: export path
        :param crack_id: crack id to be evaluated
        :param n_contours: number of contours to be evaluated from 1 to n_contours.
        :param filename_prefix: prefix for the exported files
        :param header: if True, the header is written to the file.
        :param crack_growth_defined: if True, SMART specific values are exported. If False, the corresponding columns are filled with 0.
        :param cleanup_before_export: if True, output_path is cleaned up before export. All txt files are deleted.
        :param ref_point_cf_ordering: reference point (X,Y,Z) for ordering the crack front nodes. The distance between the reference point and the two crack front end nodes is calculated and the crack front is aligned so that the closest end node as first node.
        :return: list of dataframes with the fracture parameters for each step
    """

    def prcint_to_array(inp: str, crack_id: int, n_columns: int) -> np.array:
        """ Transforms the str output from prcint to a numpy array.
        np.array[:,0] is the nodenumber.
        Columns have to be specified but they're usually the number of contours."""
        # strip header

        lst_split = list(map(lambda x: x.rstrip('\n'), inp.split(' ')))
        lst_val = []
        # skip all except nodenumber and value
        for val in lst_split:
            try:
                lst_val.append(float(val))
            except:
                continue

        # debug
        if lst_val[0] == float(crack_id):
            lst_val.pop(0)

        # return sorted by nodenum
        array = np.array(lst_val).reshape(int(len(lst_val) / n_columns), n_columns)
        array_sorted = array[array[:, 0].argsort()]

        # debug
        # remove all zero rows
        array_sorted = array_sorted[array_sorted.any(axis=1)]
        # remove rows with NN == 0
        array_sorted = array_sorted[array_sorted[:, 0] != 0, :]
        return array_sorted

    if cleanup_before_export:
        # deletes crackfront node files in the folder to avoid confusion
        filelist = [f for f in os.listdir(output_path) if f.endswith(".txt")]
        for f in filelist:
            os.remove(os.path.join(output_path, f))

    mapdl.post1()
    n_sets = mapdl.post_processing.nsets
    lst_results = []
    CONTOURS = n_contours
    COLUMNS = ExportDataColumns().crackfront(CONTOURS)

    # extract nodal information for all steps
    for result_nbr in range(n_sets):
        mapdl.set(nset=result_nbr + 1)

        # extract step info
        n_load_step = mapdl.post_processing.load_step
        n_sub_step = mapdl.post_processing.sub_step

        # create new table with metadata
        df = pd.DataFrame(columns=COLUMNS)
        df.attrs['Loadstep'] = n_load_step
        df.attrs['Substep'] = n_sub_step

        # select crackfrontnodes by component
        mapdl.cmsel('S', f'{cm_name}')
        mask = mapdl.post_processing.selected_nodes.nonzero()[0] + 1

        try:
            # fill table with relevant information
            df['ID'] = mapdl.post_processing.selected_nodes.nonzero()[0] + 1
            df['x'] = mapdl.post_processing.nodal_values('LOC', 'X')
            df['y'] = mapdl.post_processing.nodal_values('LOC', 'Y')
            df['z'] = mapdl.post_processing.nodal_values('LOC', 'Z')
            df['sigx'] = mapdl.post_processing.nodal_component_stress('X')
            df['sigy'] = mapdl.post_processing.nodal_component_stress('Y')
            df['sigz'] = mapdl.post_processing.nodal_component_stress('Z')
            df['sigeqv'] = mapdl.post_processing.nodal_eqv_stress()
            df['epelx'] = mapdl.post_processing.nodal_total_component_strain('X')
            df['epely'] = mapdl.post_processing.nodal_total_component_strain('Y')
            df['epelz'] = mapdl.post_processing.nodal_total_component_strain('Z')
            df['epseqv'] = mapdl.post_processing.nodal_total_eqv_strain()
            #TODO: Check if the following approach via prcint is still necessary
            df.loc[:, 'SIFS_K1_Kont1':f'SIFS_K1_Kont{CONTOURS}'] = prcint_to_array(mapdl.prcint(crack_id, dtype='K1'), crack_id,
                                                                                   CONTOURS + 1)[:,
                                                                   1:]
            df.loc[:, 'SIFS_K2_Kont1':f'SIFS_K2_Kont{CONTOURS}'] = prcint_to_array(mapdl.prcint(crack_id, dtype='K2'), crack_id,
                                                                                   CONTOURS + 1)[:,
                                                                   1:]
            df.loc[:, 'SIFS_K3_Kont1':f'SIFS_K3_Kont{CONTOURS}'] = prcint_to_array(mapdl.prcint(crack_id, dtype='K3'), crack_id,
                                                                                   CONTOURS + 1)[:,
                                                                   1:]
            if crack_growth_defined:
                df.loc[:, 'DLTN'] = prcint_to_array(mapdl.prcint(1, dtype='DLTN'), crack_id, 2)[:, 1]
                df.loc[:, 'DLTA'] = prcint_to_array(mapdl.prcint(1, dtype='DLTA'), crack_id, 2)[:, 1]
                df.loc[:, 'DLTK'] = prcint_to_array(mapdl.prcint(1, dtype='DLTK'), crack_id, 2)[:, 1]
                df.loc[:, 'R'] = prcint_to_array(mapdl.prcint(1, dtype='R'), crack_id, 2)[:, 1]
            else:
                df.loc[:, 'DLTN':'R'] = 0
            # print(data)
        except Exception as e:
            print('Error while extracting SIFS data')
            print(e)
            df.loc[:, 'x':'R'] = 0
            continue

        # append postprocessing information -> unit transformation
        # use only the outer 3 contours
        df['K_I_SIFS'] = df.loc[:, f'SIFS_K1_Kont{CONTOURS - 2}':f'SIFS_K1_Kont{CONTOURS}'].mean(axis=1) / np.sqrt(1000)
        df['K_II_SIFS'] = df.loc[:, f'SIFS_K2_Kont{CONTOURS - 2}':f'SIFS_K2_Kont{CONTOURS}'].mean(axis=1) / np.sqrt(
            1000)
        df['K_III_SIFS'] = df.loc[:, f'SIFS_K3_Kont{CONTOURS - 2}':f'SIFS_K3_Kont{CONTOURS}'].mean(axis=1) / np.sqrt(
            1000)

        # extract additional information
        # node ordering and node type
        if legacy_node_ordering:
            ordering_type_tuple = extract_cf_node_order(mapdl, cf_start_ref=ref_point_cf_ordering)
        else:
            ordering_type_tuple = []
            n_nodes = int(mapdl.get_value(entity='CINT', entnum=1, item1='NNOD'))
            node_order = [int(mapdl.get_value(entity='CINT', entnum=1, item1='NODE', it1num=i)) for i in
                          range(1, n_nodes + 1)]

            node1, node2 = node_order[0], node_order[-1]
            # get the absolute position of the nodes
            node1_coords, node2_coords = [], []
            for comp in ['X', 'Y', 'Z']:
                mapdl.nsel('S', 'NODE', vmin=node1)
                node1_coords.append(*mapdl.post_processing.nodal_values(item='LOC', comp=comp))
                mapdl.nsel('S', 'NODE', vmin=node2)
                node2_coords.append(*mapdl.post_processing.nodal_values(item='LOC', comp=comp))
            # calculate the distance to the reference point
            node1_dist = np.linalg.norm(np.array(node1_coords) - np.array(ref_point_cf_ordering))
            node2_dist = np.linalg.norm(np.array(node2_coords) - np.array(ref_point_cf_ordering))
            # set the node with the smaller distance to the first element of the list by inverting the list if necessary
            if node1_dist > node2_dist:
                node_order = node_order[::-1]

            # extract the type of node, corner or midside
            mapdl.cmsel("S", cm_name)
            mapdl.esln()
            mapdl.nsle(nodetype='CORNER')
            n_corner = mapdl.mesh.nnum
            mapdl.cmsel("S", cm_name)
            mapdl.esln()
            mapdl.nsle(nodetype='MID')
            n_midside = mapdl.mesh.nnum
            for node in node_order:
                if node in n_corner:
                    type = 'CORNER'
                elif node in n_midside:
                    type = 'MIDSIDE'
                else:
                    continue
                ordering_type_tuple.append((node, type))

        # put node order rel. to refpoint and type in the dataframe
        for order, (node, type) in enumerate(ordering_type_tuple):
            df.loc[df['ID'] == node, 'rel_order'] = order
            df.loc[df['ID'] == node, 'node_type'] = type

        # export
        if os.path.isfile(output_path):
            filename_prefix = os.path.split(output_path)[1].split('.')[0]
        FILEENDSWITH = f'_LS{n_load_step:03d}SS{n_sub_step:03d}.txt'
        df.to_csv(os.path.join(output_path, filename_prefix + FILEENDSWITH), sep=';', index=False, header=header)
        lst_results.append(df)
        print(f'Exported {filename_prefix+FILEENDSWITH} successfully')
    return lst_results

def extract_cf_node_order(mapdl: Mapdl,
                          cf_cm_name: str = 'CRACKFRONTNODES',
                          cf_start_ref: tuple = (0, 0, 0)) -> Tuple[Tuple[Optional[Any], str], ...]:
    """Extracts the crackfront node order from the nodemaps. Endpoints are specified by first and last element of tuple.
    :param mapdl: mapdl instance
    :param cf_cm_name: name of the crackfront node component
    :param cf_start_ref: reference point for the crackfront node order. The end nodes of the crackfront will be evaluated by their distance to the ref point and the closest end note will be set to the first node of the crackfront.
    :return: ordered tuple of the crackfront node ids and node type (corner,mid).
    """
    print('DEPRECATED: Not necessary anymore. Directly implemented in the export_fracture_parameters function.')
    # extract corner and midside nodes of the crack front and their respective elements
    mapdl.cmsel("S", cf_cm_name)
    nlist = mapdl.mesh.nnum
    mapdl.esln()
    elist = mapdl.mesh.enum
    # find out which nodes belong to which element and if they are corner or midside nodes
    # extract separately so that corner nodes are always the first and last element in the list
    # -> (value[0] and value[-1] respectively)
    n_in_e_dict_corner = {}
    n_in_e_dict_midside = {}
    corner_nodes_list = []
    midside_nodes_list = []
    for e in elist:
        #print(e)
        mapdl.esel("S", "ELEM", vmin=e)
        mapdl.nsle(nodetype='CORNER')
        #print('corner selected')
        mesh_nnum_corner = mapdl.mesh.nnum
        corner_nodes_list.extend(mesh_nnum_corner)
        n_in_e_dict_corner[e] = [n for n in mesh_nnum_corner if n in nlist]
        mapdl.esel("S", "ELEM", vmin=e)
        mapdl.nsle(nodetype='MID')
        mesh_nnum_mid = mapdl.mesh.nnum
        #time.sleep(0.5)
        midside_nodes_list.extend(mesh_nnum_mid)
        n_in_e_dict_midside[e] = [n for n in mesh_nnum_mid if n in nlist]
    # deleting appendage elements with only one node shared with the crack front -> not shape defining
    n_in_e_dict_corner = {k: v for k, v in n_in_e_dict_corner.items() if len(v) > 1}
    # these appendage elements have no midside nodes on the crackfront, therefore ignore n_in_e if empty
    n_in_e_dict_midside = {k: v for k, v in n_in_e_dict_midside.items() if v}

    # now sort the nodes along the crackfront
    # pick starter nodes, here from the first element in the dictionary
    k, sorted_nodes = list(n_in_e_dict_corner)[0], list(n_in_e_dict_corner.values())[0]
    try:
        # if the element has a midside node, inject it back between the corner nodes,
        # which is position [1] for the tuples of length 2
        sorted_nodes.insert(1, *n_in_e_dict_midside[k])
    except TypeError:
        pass

    # iterate till all nodes occurences are sorted or dropped (f.e. two elements share the same nodes)
    while len(n_in_e_dict_corner) > 0:
        element_nodes = list(n_in_e_dict_corner.items())
        #print(f'{len(n_in_e_dict_corner)} elements left')
        for k, v in element_nodes:
            # skip when no nodes are shared
            if not set(sorted_nodes).intersection(v):
                continue
            # skip when two elements share the same nodes
            # pop the second element because it adds no information
            if set(v).issubset(set(sorted_nodes)):
                n_in_e_dict_corner.pop(k)
                continue
            # inject midside nodes back between two corner nodes
            try:
                v.insert(1, *n_in_e_dict_midside[k])
            except TypeError:
                pass
            # extend the sorted_nodes list when two elements share the same corner nodes (first and third element in v)
            # iterated over all possible rotations of the corner nodes (3 Elements -> Midsidenode stays in the middle)
            if sorted_nodes[0] == v[0]:
                sorted_nodes = v[::-1] + sorted_nodes
            elif sorted_nodes[0] == v[-1]:
                sorted_nodes = v + sorted_nodes
            elif sorted_nodes[-1] == v[0]:
                sorted_nodes = sorted_nodes + v
            elif sorted_nodes[-1] == v[-1]:
                sorted_nodes = sorted_nodes + v[::-1]
            else:
                break
            # remove the element for the next loop when its nodes are appened to the sorted_nodes list
            n_in_e_dict_corner.pop(k)

    # remove dublicates in sorted_nodes and preserve order via valueless keys in dict
    # -> python equivalent of ordered sets
    sorted_nodes = list(dict.fromkeys(sorted_nodes))
    sorted_nodes_exp = []

    # order endnodes by absolute position
    # -> this is necessary because the crackfront nodes are not sorted by their absolute position
    # take first and last element of sorted_nodes and find the closest to the reference point
    node1, node2 = sorted_nodes[0], sorted_nodes[-1]
    # get the absolute position of the nodes
    node1_coords, node2_coords = [], []
    for comp in ['X', 'Y', 'Z']:
        mapdl.nsel('S', 'NODE', vmin=node1)
        node1_coords.append(*mapdl.post_processing.nodal_values(item='LOC', comp=comp))
        mapdl.nsel('S', 'NODE', vmin=node2)
        node2_coords.append(*mapdl.post_processing.nodal_values(item='LOC', comp=comp))
    # calculate the distance to the reference point
    node1_dist = np.linalg.norm(np.array(node1_coords) - np.array(cf_start_ref))
    node2_dist = np.linalg.norm(np.array(node2_coords) - np.array(cf_start_ref))
    # set the node with the smaller distance to the first element of the list by inverting the list if necessary
    if node1_dist > node2_dist:
        sorted_nodes = sorted_nodes[::-1]

    # add the corner, midside information back for return
    for n in sorted_nodes:
        if n in midside_nodes_list:
            sorted_nodes_exp.append((n, 'MIDSIDE'))
        elif n in corner_nodes_list:
            sorted_nodes_exp.append((n, 'CORNER'))

    return tuple(sorted_nodes_exp)


