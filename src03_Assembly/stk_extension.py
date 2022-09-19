# Module-Cian.1

from stk.molecular.topology_graphs.topology_graph import Edge
from stk.molecular.topology_graphs.metal_complex.metal_complex import MetalComplex
from stk.molecular.topology_graphs.metal_complex.vertices import MetalVertex, UnaligningVertex
from stk.molecular.topology_graphs.topology_graph.vertex import Vertex
from scipy.spatial.distance import euclidean
from stk.utilities import get_projection


class tridentateVertex(Vertex):

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class tridentate(MetalComplex):
    _metal_vertex_prototypes = (MetalVertex(0, (0, 0, 0)),)
    _ligand_vertex_prototypes = (tridentateVertex(1, (0, 0, 0)),)

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0.1, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0.1, 0),
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(-0.1, 0, 0)
        ),
    )


class BiDentateLigandVertex(Vertex):

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


class bidentate(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )
    _ligand_vertex_prototypes = (
        BiDentateLigandVertex(1, (-2.0, 0, - 2.0)),
    )

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
            position=(-2.0, 0, 0),
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=(0, 0, -2.0),
        ),
    )


class MonoDentateLigandVertex(Vertex):
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


class monodentate(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonoDentateLigandVertex(1, (0, 0, 2.5)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class monodentate_flipped(MetalComplex):
    _metal_vertex_prototypes = (
        MetalVertex(0, (0, 0, 0)),
    )

    _ligand_vertex_prototypes = (
        MonoDentateLigandVertex(1, (0, 0, -2.5)),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
    )


class complex_topology(MetalComplex):
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
