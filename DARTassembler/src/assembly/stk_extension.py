import warnings
import numpy as np
from DARTassembler.src.constants.Periodic_Table import DART_Element as element
import logging
from scipy.spatial.distance import euclidean

# Get stk functions which were refactored at some point
try:
    from stk.molecular.topology_graphs.topology_graph import Edge
    from stk.molecular.topology_graphs.metal_complex.metal_complex import MetalComplex
    from stk.molecular.topology_graphs.metal_complex.vertices import MetalVertex, UnaligningVertex
    from stk.molecular.topology_graphs.topology_graph.vertex import Vertex
    from stk.utilities import get_projection
except (ImportError, ModuleNotFoundError):
    from stk._internal.topology_graphs.edge import Edge
    from stk._internal.topology_graphs.metal_complex.metal_complex import MetalComplex
    from stk._internal.topology_graphs.metal_complex.vertices import MetalVertex, UnaligningVertex
    from stk._internal.topology_graphs.vertex import Vertex
    from stk._internal.utilities.utilities import get_projection



class TridentateLigandVertex(Vertex):
    """
    Some modifications of Vertex elements of stk as a preparation for Tridentate ligands
    No build-in functionality for that in stk
    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class BidentateLigandVertex(Vertex):
    """
    Some modifications of Vertex elements of stk as a preparation for Bidentate ligands
    No build-in functionality for that in stk
    """
    # Todo: This function is a bug which shows for the bidentate ligand consisting of just two oxygens. I assume it's in general for all bi-atomic bidentate ligands. The issue is that in the last call to building_block.with_rotation_between_vectors(), the rotation makes the rotation matrix go NaN.

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        assert (
                building_block.get_num_functional_groups() == 2
        ), (
            f'{building_block} needs to have exactly 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )

        fg0_position, fg1_position = (
            building_block.get_centroid(fg.get_placer_ids())
            for fg in building_block.get_functional_groups()
        )
        edge_position1, edge_position2 = (
            edge.get_position() for edge in edges
        )
        building_block = building_block.with_rotation_between_vectors(
            start=fg1_position - fg0_position,
            target=edge_position2 - edge_position1,
            origin=building_block.get_centroid(),
        )

        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        core_to_placer = placer_centroid - core_centroid

        fg0_position, fg1_position = (
            building_block.get_centroid(fg.get_placer_ids())
            for fg in building_block.get_functional_groups()
        )
        fg_vector = fg1_position - fg0_position

        fg_vector_projection = get_projection(
            start=core_to_placer,
            target=fg_vector,
        )

        edge_centroid = (
                sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=core_to_placer - fg_vector_projection,
            target=edge_centroid - self._position,
            origin=building_block.get_centroid(),
        )

        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class MonodentateLigandVertex(Vertex):
    """
    Places monodentate ligand in a :class:`.MetalComplex`.
    """

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        assert (
                building_block.get_num_functional_groups() == 1
        ), (
            f'{building_block} needs to have exactly 1 functional '
            'group but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        fg, = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(
            atom_ids=fg.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        edge_centroid = (
                sum(edge.get_position() for edge in edges) / len(edges)
        )
        return building_block.with_rotation_between_vectors(
            start=fg_centroid - core_centroid,
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {0: edges[0].get_id()}


class Tridentate(MetalComplex):
    """
    Extension of some basic stk class to model Tridentate Ligands
    """
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)  # _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (TridentateLigandVertex(1, (0, 0, 0)),)  # _ligand_vertex_prototypes = (TridentateLigandVertex(1, (0, 0, 0)),)

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0.1, 0, 0),  # position=(0.1, 0, 0), We know that these rotations worked
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0.1, 0),  # position=(0, 0.1, 0),
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(-0.1, 0, 0)  # position=(-0.1, 0, 0)
        ),
    )


class Bidentate(MetalComplex):
    """
    Extension of some basic stk class to model Tridentate Ligands
    """
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (BidentateLigandVertex(1, (-1.0, 0, - 1.0)),)

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    # This paticular arrangement has the ligand
    # cordination on the Left and Bottom sites
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(-1.0, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0, -1.0),
        ),
    )


class Bidentate_Planar_Right(MetalComplex):
    """
    Extension of some basic stk class to model Tridentate Ligands
    """
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (BidentateLigandVertex(1, (1.6, 0, 0)),)

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    # This paticular arrangement has the ligand
    # cordination on the Left and Bottom sites
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(1.31, 1.31, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(1.31, -1.31, 0),
        ),
    )


class Bidentate_Planar_Left(MetalComplex):
    """
    Extension of some basic stk class to model Tridentate Ligands
    """
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (BidentateLigandVertex(1, (-1.6, 0, 0)),)

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    # This paticular arrangement has the ligand
    # cordination on the Left and Bottom sites
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(-0.5, 0.5, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(-0.5, -0.5, 0),
        ),
    )


class Monodentate(MetalComplex):
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (MonodentateLigandVertex(1, (0, 0, 0)),)
    _edge_prototypes = (Edge(id=0, vertex1=_metal_vertex_prototypes[0], vertex2=_ligand_vertex_prototypes[0], ),)


class complex_topology_three(MetalComplex):
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (
        UnaligningVertex(1, (0, 0, 0)),
        UnaligningVertex(2, (0, 0, 0)),
        UnaligningVertex(3, (0, 0, 0)),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=(0, 0, 0),
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=(0, 0, 0),
        ),

    )


class complex_topology_two(MetalComplex):
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (
        UnaligningVertex(1, (0, 0, 0)),
        UnaligningVertex(2, (0, 0, 0)),
        UnaligningVertex(3, (0, 0, 0)),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=(0, 0, 0),
        ),

    )


########################################################################################################################################################################################################

class Monodentate_Top(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (0, 0, 1.83)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class Monodentate_Bottom(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (0, 0, -1.83)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class Monodentate_Front_Right(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (-1.2, 1.2, 0)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class Monodentate_Front_Left(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (-1.2, -1.2, 0)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class Monodentate_Back_Left(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (-1.2, 1.2, 0)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class Monodentate_Back_Right(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonodentateLigandVertex(1, (1.2, 1.2, 0)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class monodentate_coordinating_distance:
    def __init__(self, metal, ligand, offset: float = 0):
        # todo: it is super important that we cite the source of the covalent radii
        self.metal = metal
        self.ligand = ligand
        self.offset = offset
        self.metal_covalent_radius = float(element(self.metal).covalent_radius) / 100.0
        self.ligand_covalent_radius = float(element(ligand.local_elements[0]).covalent_radius) / 100.0
        logging.debug(f"Monodentate: The metal atomic radius is:{self.metal_covalent_radius}")
        logging.debug(f"Monodentate: The ligand atomic radius is:{self.ligand_covalent_radius}")
        logging.debug("#####")

    def Top(self):
        return [0, 0, abs(self.metal_covalent_radius + self.ligand_covalent_radius + self.offset)]

    def Bottom(self):
        return [0, 0, -1 * abs(self.metal_covalent_radius + self.ligand_covalent_radius + self.offset)]

    def Front_Left(self):
        distance = abs(self.metal_covalent_radius + self.ligand_covalent_radius + self.offset)
        component = np.sqrt(0.5 * (distance ** 2))
        return [-1 * component, -1 * component, 0]

    def Middle_Left(self):
        return [-1 * abs(self.metal_covalent_radius + self.ligand_covalent_radius + self.offset), 0, 0]

    def Back_Left(self):
        distance = abs(self.metal_covalent_radius + self.ligand_covalent_radius + self.offset)
        component = np.sqrt(0.5 * (distance ** 2))
        return [-1 * component, component, 0]


class Bidentate_coordinating_distance:
    def __init__(self, metal, ligand, offset: float = 0):
        self.metal = metal
        self.offset = offset
        self.metal_covalent_radius = float(element(self.metal).covalent_radius) / 100.0
        self.ligand_covalent_radius_0 = float(element(ligand.local_elements[0]).covalent_radius) / 100.0
        self.ligand_covalent_radius_1 = float(element(ligand.local_elements[1]).covalent_radius) / 100.0
        self.ligand_coord_0 = ligand.coordinates[ligand.ligand_to_metal[0]][1]
        self.ligand_coord_1 = ligand.coordinates[ligand.ligand_to_metal[1]][1]
        self.coord_atom_dist = np.linalg.norm(np.array(self.ligand_coord_0) - np.array(self.ligand_coord_1))
        with np.errstate(all='ignore'):     # Ignore all numpy warnings such as taking the sqrt of a negative number
            self.median = 0.5 * np.sqrt((2 * (self.ligand_covalent_radius_0 + self.metal_covalent_radius) ** 2) + (2 * (self.ligand_covalent_radius_1 + self.metal_covalent_radius) ** 2) - self.coord_atom_dist ** 2)

        logging.debug(f"Bidentate: The metal atomic radius is:{self.metal_covalent_radius}")
        logging.debug(f"Bidentate: The ligand atomic radius is:{self.ligand_covalent_radius_0}")
        logging.debug(f"Bidentate: The ligand atomic radius is:{self.ligand_covalent_radius_1}")
        logging.debug("#####")

    def calculate_mean_distance(self):
        if np.isnan(self.median):
            logging.debug("!!!Warning!!! -> nan encountered in Bidentate placer -> Returning default value")
            return 1.6
        else:
            return self.median

    def Left(self):
        return [(-1.0 * self.calculate_mean_distance()), 0, 0]

    def Right(self):
        return [self.calculate_mean_distance(), 0, 0]

    def Bottom(self):
        return [0, 0, -1*self.calculate_mean_distance()]



