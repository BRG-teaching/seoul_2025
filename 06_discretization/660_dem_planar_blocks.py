#! python 3
# r: compas, compas_cgal==0.9.1, compas_libigl==0.7.4, tessagon, compas_dem
# venv: aa

# Standard library imports
from pathlib import Path
from typing import Optional, Union

# COMPAS core imports
import compas
from compas.data import Data
from compas.datastructures import Mesh
from compas.itertools import pairwise
from compas.scene import Scene

# COMPAS geometry imports
from compas.geometry import (
    Box,
    Brep,
    Frame,
    Line,
    Plane,
    Point,
    Polygon,
    Polyline,
    Transformation,
    Translation,
    Vector,
    distance_point_point,
    distance_point_point_sqrd,
    intersection_line_plane,
    intersection_plane_plane,
    is_parallel_vector_vector
)

# COMPAS extensions imports
from compas_cgal.booleans import boolean_difference_mesh_mesh, boolean_union_mesh_mesh
from compas_cgal.meshing import trimesh_dual
from compas_dem.elements import Block
from compas_dem.models import BlockModel
from compas_libigl.mapping import map_mesh
from compas_libigl.parametrisation import trimesh_lsc_mapping
from compas_model.elements import Element
from compas_model.modifiers import Modifier

# Rhino-specific imports
from Rhino.DocObjects import CurveObject
from compas_rhino.conversions import curve_to_compas_polyline, meshobject_to_compas
from compas_rhino.objects import find_object, select_mesh, select_objects
from compas_rhino.objects import select_points
from compas_rhino.conversions import pointobject_to_compas


is_brep = False
tmin = 0.3
tmax = 0.4
compute_interactions = True

# =============================================================================
# RhinoVault mesh
# =============================================================================

guid = select_mesh()
mesh = meshobject_to_compas(guid)

# points form rhino
guids = select_points()
points = []
for guid in guids:
    point = pointobject_to_compas(guid)
    points.append(point)

# =============================================================================
# Remesh and Flatten
# =============================================================================


fixed_vertices = []
for point in points:
    for v in mesh.vertices():
        mesh_point = mesh.vertex_point(v)
        if mesh_point.distance_to_point(point) < 0.01:
            fixed_vertices.append(v)
            break

V, F, DV, DF = trimesh_dual(
    mesh.to_vertices_and_faces(True),
    length_factor=1.0,
    number_of_iterations=10,
    angle_radians=0.9,
    scale_factor=1.0,
    fixed_vertices=fixed_vertices,
)


mesh_remeshed = Mesh.from_vertices_and_faces(V, F)
pattern = Mesh.from_vertices_and_faces(DV, DF)
pattern.unify_cycles()

mv = []
mn = []
for vertex in pattern.vertices():
    mv.append(pattern.vertex_point(vertex))
    mn.append(pattern.vertex_normal(vertex))

# ==============================================================================
# Least-§squares conformal map - to see where the pattern is mapped.
# ==============================================================================

uv = trimesh_lsc_mapping((V, F))
mesh_lscm = mesh_remeshed.copy()
for i in range(mesh_remeshed.number_of_vertices()):
    mesh_lscm.vertex_attributes(i, "xyz", [uv[i][0], uv[i][1], 0])

# ==============================================================================
# Select Pattern and Map
# ==============================================================================



guids = select_objects()
polygons = []

for guid in guids:
    obj = find_object(guid)
    if isinstance(obj, CurveObject):
        polyline = curve_to_compas_polyline(obj.Geometry)
        polygon = Polygon(polyline.points[:-1])
        polygons.append(polygon)

# Session
# polygons = compas.json_load(Path(__file__).parent / "2d_pattern_non_parallel.json")
pattern = Mesh.from_polygons(polygons)
PV, PF = pattern.to_vertices_and_faces()
mv, mf, mn, mb, mg = map_mesh((V, F), (PV, PF), clip_boundaries=False, fixed_vertices=[], tolerance=1e-3)
pattern = Mesh.from_vertices_and_faces(mv, mf)
pattern.unify_cycles()

# ===========================================================================
# Thickness
# =============================================================================

pattern.update_default_vertex_attributes(thickness=0)

zvalues: list[float] = pattern.vertices_attribute(name="z")  # type: ignore
zmin = min(zvalues)
zmax = max(zvalues)

for vertex in pattern.vertices():
    point = Point(*mv[vertex])
    normal = mn[vertex]
    z = (point.z - zmin) / (zmax - zmin)
    thickness = (1 - z) * (tmax - tmin) + tmin
    pattern.vertex_attribute(vertex, name="thickness", value=thickness)

# =============================================================================
# Model
# =============================================================================


class OffsetPlanarBlocks(object):
    """
    Offset a mesh to create blocks with planar sides and chamfered corners.

    Attributes
    ----------
    original_mesh : :class:`compas.datastructures.Mesh`
        The original input mesh before any operations.
    mesh : :class:`compas.datastructures.Mesh`
        A working copy of the mesh with offsets applied.
    offset : float
        The offset distance.
    chamfer : float
        The chamfer distance.
    thickness_scale_bottom : float
        The thickness scale for the bottom.
    thickness_scale_top : float
        The thickness scale for the top.
    project_bottom : bool
        True if bottom face should be projected, False otherwise.
    project_top : bool
        True if top face should be projected, False otherwise.
    tolerance_parallel : float
        The tolerance for parallelism.
    vertex_normals : list[:class:`compas.geometry.Vector`]
        The vertex normals of the mesh.
    edge_frames : dict[tuple[int, int], :class:`compas.geometry.Frame`]
        Frames at each edge of the mesh.
    vertex_frames : dict[tuple[int, int], :class:`compas.geometry.Frame`]
        Frames at each vertex of the mesh.
    parallel_edges : dict[tuple[int, int], bool]
        Dictionary mapping edges to boolean indicating if they are parallel.
    face_vertex_directions_0 : dict[tuple[int, int], :class:`compas.geometry.Line`]
        First vertex direction lines for each face.
    face_vertex_directions_1 : dict[tuple[int, int], :class:`compas.geometry.Line`]
        Second vertex direction lines for each face.
    blocks : list[:class:`compas.datastructures.Mesh`]
        The generated block meshes.
    block_frames : list[:class:`compas.geometry.Frame`]
        The frames of the generated blocks.
    """

    def __init__(
        self,
        mesh,
        offset=0.0,
        chamfer=0.05,
        thickness_scale_bottom=0.0,
        thickness_scale_top=1.0,
        project_bottom=True,
        project_top=False,
        tolerance_parallel=0.5,
        vertex_normals=None,
    ):
        self.original_mesh = mesh
        self.mesh = mesh.copy()
        self.offset = offset
        self.chamfer = chamfer
        self.thickness_scale_bottom = thickness_scale_bottom
        self.thickness_scale_top = thickness_scale_top
        self.project_bottom = project_bottom
        self.project_top = project_top
        self.tolerance_parallel = tolerance_parallel
        self.vertex_normals = vertex_normals

        # Results
        self.edge_frames = {}
        self.vertex_frames = {}
        self.parallel_edges = {}
        self.face_vertex_directions_0 = {}
        self.face_vertex_directions_1 = {}
        self.blocks = []
        self.block_frames = []

        self._check_mesh_thickness()
        self._offset_mesh()
        self._compute_edge_frames()
        self._compute_vertex_frames()
        self._get_vertex_lines()
        self._get_blocks()

    def _check_mesh_thickness(self):
        """Check if mesh has thickness values, set them if needed or raise an error.

        If offset is not zero, assigns that value as thickness to all vertices.
        Otherwise, verifies that thickness attributes are present on all vertices.

        Returns
        -------
        bool
            True if thickness values were set or already existed, False otherwise.

        Raises
        ------
        ValueError
            If mesh has no thickness values and offset is zero.
        """
        if self.offset != 0:
            for vertex in self.mesh.vertices():
                self.mesh.vertex_attribute(vertex, "thickness", self.offset)
        elif not self.mesh.vertices_attribute("thickness"):
            raise ValueError("Mesh must have thickness values.")

    def _offset_mesh(self):
        """Apply offset to mesh vertices based on thickness and scale.

        Displaces each vertex along its normal by its thickness value
        multiplied by the thickness_scale_bottom factor.

        Returns
        -------
        :class:`compas.datastructures.Mesh`
            The mesh with offset applied.
        """
        for vertex in self.mesh.vertices():
            normal = self.vertex_normals[vertex]  # self.mesh.vertex_normal(vertex)
            thickness = self.mesh.vertex_attribute(vertex, "thickness")
            self.mesh.set_vertex_point(vertex, self.mesh.vertex_point(vertex) + normal * thickness * self.thickness_scale_bottom)

    def _compute_edge_frames(self):
        """Compute frames at mesh edges."""
        for edge in self.mesh.edges():
            faces = self.mesh.edge_faces(edge)
            origin = (self.mesh.vertex_point(edge[1]) + self.mesh.vertex_point(edge[0])) / 2
            normal0 = self.mesh.face_normal(faces[0]) if faces[0] is not None else Vector(0, 0, 0)
            normal1 = self.mesh.face_normal(faces[1]) if faces[1] is not None else Vector(0, 0, 0)
            normal = (normal0 + normal1) / 2
            direction = self.mesh.edge_direction(edge)
            self.edge_frames[edge] = Frame(origin, normal, direction)
            self.edge_frames[edge[::-1]] = Frame(origin, normal, -direction)

    def _compute_vertex_frames(self):
        """Compute frames at the corners of the block for chamfering."""
        for face in self.mesh.faces():
            halfedges = self.mesh.face_halfedges(face)

            for i in range(len(halfedges)):
                # indexing
                prev_idx = (i - 1) % len(halfedges)
                curr_idx = i
                prev_edge = halfedges[prev_idx]
                curr_edge = halfedges[curr_idx]

                # orientation
                origin = self.mesh.vertex_point(prev_edge[1])
                zaxis_prev = self.edge_frames[prev_edge].zaxis
                zaxis_curr = self.edge_frames[curr_edge].zaxis
                xaxis_prev = self.edge_frames[prev_edge].xaxis
                xaxis_curr = self.edge_frames[curr_edge].xaxis
                self.parallel_edges[curr_edge] = False

                # when frames are parallel: _ _ , construct 90 degrees rotated frame _|_
                if is_parallel_vector_vector(zaxis_prev, zaxis_curr, tol=self.tolerance_parallel):
                    x_vector = zaxis_prev + zaxis_curr
                    y_vector = xaxis_prev + xaxis_curr
                    self.parallel_edges[curr_edge] = True
                # when frames are not parallel: _ / , construct frame from intersection of planes _\/
                else:
                    plane_intersection = intersection_plane_plane(Plane.from_frame(self.edge_frames[prev_edge]), Plane.from_frame(self.edge_frames[curr_edge]))
                    if plane_intersection:
                        x_vector = zaxis_curr - zaxis_prev
                        y_vector = Point(*plane_intersection[1]) - Point(*plane_intersection[0])
                    else:
                        raise Exception("No Plane-Plane intersection in get_corner_frames method.")

                self.vertex_frames[curr_edge] = Frame(origin, x_vector, y_vector)

                # Calculate angle in radians from dot product of unit vectors
                # Map from dot product to chamfer factor:
                # - When dot = -1 (180°), factor = 0
                # - When dot approaches 1 (0°), factor approaches 1
                dot_product = zaxis_prev.dot(zaxis_curr)  # Range: [-1, 1]
                chamfer_factor = (1 - dot_product) / 2  # Range: [0, 1]

                if not self.mesh.is_vertex_on_boundary(curr_edge[0]) and not is_parallel_vector_vector(zaxis_prev, zaxis_curr, tol=self.tolerance_parallel):
                    self.vertex_frames[curr_edge].translate(self.vertex_frames[curr_edge].zaxis * self.chamfer * chamfer_factor)

    def _get_vertex_lines(self):
        """Calculate vertex normal lines for each face.

        Raises
        ------
        Exception
            If plane-plane or line-plane intersections cannot be computed.
        """
        for face in self.mesh.faces():
            halfedges = self.mesh.face_halfedges(face)

            for i in range(len(halfedges)):
                # ._e0_.e1_.
                prev_edge = halfedges[(i - 1) % len(halfedges)]
                curr_edge = halfedges[i]

                # 1. we intersect vertex and edge frames _.\
                plane_intersection = intersection_plane_plane(Plane.from_frame(self.edge_frames[prev_edge]), Plane.from_frame(self.vertex_frames[curr_edge]))

                if not plane_intersection:
                    raise Exception("No Plane-Plane intersection.")

                # create line from intersection
                line = Line(plane_intersection[0], plane_intersection[1])

                # the intersection line is oriented in the same direction as the vertex normal
                p_o = self.mesh.vertex_point(curr_edge[0]) + self.mesh.vertex_normal(curr_edge[0])
                if distance_point_point_sqrd(line[0] - line.vector, p_o) < distance_point_point_sqrd(line[1] + line.vector, p_o):
                    line = Line(line[1], line[0])

                # We cut the line only for vizualization purposes to show the new face normals from approximate vertex position
                pl = Plane(self.mesh.vertex_point(curr_edge[0]), self.mesh.vertex_normal(curr_edge[0]))
                line_plane_intersection = intersection_line_plane(line, pl)

                if not line_plane_intersection:
                    raise Exception("No Line-Plane intersection.")

                self.face_vertex_directions_0[halfedges[i]] = Line(Point(*line_plane_intersection), Point(*line_plane_intersection) + line.vector)

                # 2. Second we intersect vertex and edge frames \./
                plane_intersection = intersection_plane_plane(Plane.from_frame(self.edge_frames[curr_edge]), Plane.from_frame(self.vertex_frames[curr_edge]))

                if not plane_intersection:
                    raise Exception("No Plane-Plane intersection.")

                line = Line(plane_intersection[0], plane_intersection[1])

                # orient line in the direction of the normal
                p_o = self.mesh.vertex_point(curr_edge[0]) + self.mesh.vertex_normal(curr_edge[0])
                if distance_point_point_sqrd(line[0] - line.vector, p_o) < distance_point_point_sqrd(line[1] + line.vector, p_o):
                    line = Line(line[1], line[0])

                # We cut the line only for vizualization purposes to show the new face normals from approximate vertex position
                pl = Plane(self.mesh.vertex_point(curr_edge[0]), self.mesh.vertex_normal(curr_edge[0]))
                line_plane_intersection = intersection_line_plane(line, pl)

                if not line_plane_intersection:
                    raise Exception("No Line-Plane intersection.")

                self.face_vertex_directions_1[halfedges[i]] = Line(Point(*line_plane_intersection), Point(*line_plane_intersection) + line.vector)

    def _get_blocks(self):
        """Generate blocks with planar sides and chamfered corners."""

        for face in self.mesh.faces():
            vertex_thickness = self.mesh.vertices_attribute("thickness", keys=self.mesh.face_vertices(face))

            # intrados plane from mesh face centroid and its normal
            bottom_points = []
            origin = self.mesh.face_centroid(face)
            plane = Plane(origin, self.mesh.face_normal(face))
            face_vertices = self.mesh.face_vertices(face)

            p0 = plane.projected_point(self.mesh.vertex_point(face_vertices[0]))
            p1 = plane.projected_point(self.mesh.vertex_point(face_vertices[1]))
            x = p1 - p0
            y = x.cross(-plane.normal)
            orientation_frame = Frame(origin, x, y)

            # We iterate over the face vertices and intersect the normal lines with the plane
            face_vertices = self.mesh.face_vertices(face)
            for idx, halfedge in enumerate(self.mesh.face_halfedges(face)):
                # If we want to keep a continuous mesh, construct planes at each vertex using their normals
                if not self.project_bottom:
                    plane = Plane(
                        self.mesh.vertex_point(face_vertices[idx]) + self.mesh.vertex_normal(face_vertices[idx]) * vertex_thickness[idx] * 0,
                        self.mesh.vertex_normal(face_vertices[idx]),
                    )
                    p0 = plane.projected_point(self.mesh.vertex_point(face_vertices[(idx - 1) % len(face_vertices)]))
                    p1 = plane.projected_point(self.mesh.vertex_point(face_vertices[idx]))
                    p2 = plane.projected_point(self.mesh.vertex_point(face_vertices[(idx + 1) % len(face_vertices)]))
                    v0 = p1 - p0
                    v1 = p1 - p2
                    v0.unitize()
                    v1.unitize()

                    if not is_parallel_vector_vector(v0, -v1, tol=self.tolerance_parallel):
                        p0 = self.mesh.vertex_point(face_vertices[(idx - 1) % len(face_vertices)])
                        p1 = self.mesh.vertex_point(face_vertices[idx])
                        p2 = self.mesh.vertex_point(face_vertices[(idx + 1) % len(face_vertices)])
                        v0 = p1 - p0
                        v1 = p1 - p2
                        v0.unitize()
                        v1.unitize()
                        za = v1.cross(v0)
                        za.unitize()
                        plane = Plane(p0 + za * vertex_thickness[idx] * 0, za)

                line = self.face_vertex_directions_0[halfedge]
                result_pt = plane.intersection_with_line(line)
                bottom_points.append(result_pt)

                line = self.face_vertex_directions_1[halfedge]
                result_pt = plane.intersection_with_line(line)
                bottom_points.append(result_pt)

            # extrados plane from mesh face centroid and its normal
            top_points = []
            thickness_average = sum(vertex_thickness) / len(vertex_thickness) * self.thickness_scale_top
            origin = self.mesh.face_centroid(face)
            plane = Plane(origin + self.mesh.face_normal(face) * thickness_average, self.mesh.face_normal(face))

            if self.project_top:
                orientation_frame.translate(orientation_frame.zaxis * thickness_average)
                orientation_frame.flip()

            # We iterate over the face vertices and intersect the normal lines with the plane
            face_vertices = self.mesh.face_vertices(face)
            for idx, halfedge in enumerate(self.mesh.face_halfedges(face)):
                # If we want to keep a continuous mesh, construct planes at each vertex using their normals
                if not self.project_top:
                    plane = Plane(
                        self.mesh.vertex_point(face_vertices[idx]) + self.mesh.vertex_normal(face_vertices[idx]) * vertex_thickness[idx] * self.thickness_scale_top,
                        self.mesh.vertex_normal(face_vertices[idx]),
                    )
                    p0 = plane.projected_point(self.mesh.vertex_point(face_vertices[(idx - 1) % len(face_vertices)]))
                    p1 = plane.projected_point(self.mesh.vertex_point(face_vertices[idx]))
                    p2 = plane.projected_point(self.mesh.vertex_point(face_vertices[(idx + 1) % len(face_vertices)]))
                    v0 = p1 - p0
                    v1 = p1 - p2
                    v0.unitize()
                    v1.unitize()

                    if not is_parallel_vector_vector(v0, -v1, tol=self.tolerance_parallel):
                        p0 = self.mesh.vertex_point(face_vertices[(idx - 1) % len(face_vertices)])
                        p1 = self.mesh.vertex_point(face_vertices[idx])
                        p2 = self.mesh.vertex_point(face_vertices[(idx + 1) % len(face_vertices)])
                        v0 = p1 - p0
                        v1 = p1 - p2
                        v0.unitize()
                        v1.unitize()
                        za = v1.cross(v0)
                        za.unitize()
                        plane = Plane(p0 + za * vertex_thickness[idx] * self.thickness_scale_top, za)

                line = self.face_vertex_directions_0[halfedge]
                result_pt = plane.intersection_with_line(line)
                top_points.append(result_pt)

                line = self.face_vertex_directions_1[halfedge]
                result_pt = plane.intersection_with_line(line)
                top_points.append(result_pt)

            sides = []
            for (a, b), (aa, bb) in zip(pairwise(bottom_points + bottom_points[:1]), pairwise(top_points + top_points[:1])):
                sides.append([a, b, bb, aa])

            polygons = [bottom_points[::-1], top_points] + sides

            polygons_no_duplicates = []

            # remove duplicates
            for points in polygons:
                new_points = [points[0]]
                for i in range(1, len(points)):
                    if distance_point_point(points[i], new_points[-1]) > 1e-3:
                        new_points.append(points[i])
                polygons_no_duplicates.append(new_points)

            block = Mesh.from_polygons(polygons_no_duplicates)
            self.blocks.append(block)
            self.block_frames.append(orientation_frame)


offset_planar_blocks = OffsetPlanarBlocks(
    mesh=pattern,
    offset=0,
    chamfer=0.1,
    thickness_scale_bottom=-0.5,
    thickness_scale_top=1,
    project_bottom=False,
    project_top=True,
    tolerance_parallel=0.25,
    vertex_normals=mn,
)

solid_meshes, block_frames, e_frames = offset_planar_blocks.blocks, offset_planar_blocks.block_frames, offset_planar_blocks.edge_frames

# extract connectivity
connected_elements: list[tuple[int, int, Line, Frame, int, int]] = []

for edge in pattern.edges():
    faces = pattern.edge_faces(edge)
    if faces[0] is None or faces[1] is None:
        continue

    e_frame = e_frames[edge]

    line = pattern.edge_line(edge)
    n = e_frame.zaxis
    o = line.point_at(0.5)
    x = line.vector
    y = n.cross(x)

    connected_elements.append((faces[0], faces[1], line, Frame(o, x, y)))

# add solid blocks
model = BlockModel()
blocks = []
for solid_mesh in solid_meshes:
    if is_brep:
        # polygons = solid_mesh.to_polygons()
        # brep = Brep.from_polygons(polygons)
        brep = Brep.from_mesh(solid_mesh)
        block: Block = Block(geometry=brep)

    else:
        block: Block = Block.from_mesh(solid_mesh)

    model.add_element(block)
    blocks.append(block)

# add interactions
for element in connected_elements:
    model.add_interaction(blocks[element[0]], blocks[element[1]])
    model.graph.edge_attribute((element[0], element[1]), "line", element[2])
    model.graph.edge_attribute((element[0], element[1]), "frame", element[3])

for idx, block in enumerate(blocks):
    model.graph.node_attribute(idx, "frame", block_frames[idx])


# =============================================================================
# Modifiers
# model.add_shear_keys(edges_frames, edges_lines, shape0, shape1)
# =============================================================================


class BooleanDifferenceModifier(Modifier):
    """Boolean difference modifier.

    Parameters
    ----------
    name : str, optional
        The name of the interaction.

    """

    def apply(
        self,
        source,
        targetgeometry: Union[Brep, Mesh],
    ) -> Union[Brep, Mesh]:
        """Boolean difference the source geometry with the shape and boolean union the target geometry with the shape.

        Parameters
        ----------
        source : :class:`compas.geometry.Brep` | :class:`compas.datastructures.Mesh`
            The source of the modification.
        targetgeometry : :class:`compas.geometry.Brep` | :class:`compas.datastructures.Mesh`
            The target of the modification.

        Returns
        -------
        Brep | Mesh
            The modified source geometry.

        """
        if isinstance(source.elementgeometry, Mesh) and isinstance(targetgeometry, Mesh):
            VS, FS = source.elementgeometry.to_vertices_and_faces(True)
            VT, FT = targetgeometry.to_vertices_and_faces(True)
            V, F = boolean_difference_mesh_mesh((VT, FT), (VS, FS))
            mesh = Mesh.from_vertices_and_faces(V, F)
            mesh.attributes = targetgeometry.attributes
            return mesh

        if isinstance(source.elementgeometry, Brep) and isinstance(targetgeometry, Brep):
            result = targetgeometry - source.elementgeometry
            return result

        raise ValueError(f"Source and target geometry must be of the same type. Source: {type(source.elementgeometry)}, Target: {type(targetgeometry)}")


class BooleanUnionModifier(Modifier):
    """Boolean union modifier.

    Parameters
    ----------
    name : str, optional
        The name of the interaction.

    """

    def apply(
        self,
        source,
        targetgeometry: Union[Brep, Mesh],
    ) -> Union[Brep, Mesh]:
        """Boolean difference the source geometry with the shape and boolean union the target geometry with the shape.

        Parameters
        ----------
        source : :class:`compas.geometry.Brep` | :class:`compas.datastructures.Mesh`
            The source of the modification.
        targetgeometry : :class:`compas.geometry.Brep` | :class:`compas.datastructures.Mesh`
            The target of the modification.

        Returns
        -------
        Brep | Mesh
            The modified source geometry.

        """
        if isinstance(source.elementgeometry, Mesh) and isinstance(targetgeometry, Mesh):
            VS, FS = source.elementgeometry.to_vertices_and_faces(True)
            VT, FT = targetgeometry.to_vertices_and_faces(True)
            V, F = boolean_union_mesh_mesh((VT, FT), (VS, FS))
            mesh = Mesh.from_vertices_and_faces(V, F)
            mesh.attributes = targetgeometry.attributes
            return mesh

        if isinstance(source.elementgeometry, Brep) and isinstance(targetgeometry, Brep):
            result = targetgeometry + source.elementgeometry
            return result

        raise ValueError(f"Source and target geometry must be of the same type. Source: {type(source.elementgeometry)}, Target: {type(targetgeometry)}")


class Interface(Element):
    """Class representing interface elements.

    Parameters
    ----------
    shape : :class:`compas.datastructures.Mesh` or :class:`compas.geometry.Brep`
        The base shape of the interface.
    frame : :class:`compas.geometry.Frame`, optional
        The coordinate frame of the interface.
    name : str, optional
        The name of the element.

    Attributes
    ----------
    shape : :class:`compas.datastructures.Mesh` or :class:`compas.geometry.Brep`
        The base shape of the interface.

    """

    _geometry: Union[Mesh, Brep]

    @property
    def __data__(self) -> dict:
        data = super().__data__
        data["geometry"] = self._geometry
        return data

    def __init__(
        self,
        geometry: Union[Mesh, Brep],
        transformation: Optional[Transformation] = None,
        name: Optional[str] = "Interface",
    ) -> None:
        super().__init__(geometry=geometry, transformation=transformation, name=name)

    def compute_elementgeometry(self, include_features: bool = False) -> Mesh:
        return self._geometry

    @classmethod
    def from_box_mesh(cls, frame: Frame, xsize: float, ysize: float, zsize: float) -> "Interface":
        box = Box(xsize, ysize, zsize, frame)
        return cls(box.to_mesh())

    @classmethod
    def from_box_brep(cls, frame: Frame, xsize: float, ysize: float, zsize: float) -> "Interface":
        box = Box(xsize, ysize, zsize, frame)
        return cls(Brep.from_box(box))


if compute_interactions:
    elements = list(model.elements())
    modifier_pairs = []

    for edge in model.graph.edges():
        edge_attrs = model.graph.edge_attributes(edge)
        frame = edge_attrs["frame"]
        line = edge_attrs["line"]
        length = line.length
        interface0 = Interface.from_box_mesh(frame, length * 0.6, 0.07, 0.1) if not is_brep else Interface.from_box_brep(frame, length * 0.6, 0.07, 0.1)
        interface1 = Interface.from_box_mesh(frame, length * 0.65, 0.1, 0.1) if not is_brep else Interface.from_box_brep(frame, length * 0.65, 0.1, 0.1)
        modifier_pairs.append([interface0, interface1, elements[edge[0]], elements[edge[1]]])

    # Add elements and modifiers
    for interface0, interface1, elem0, elem1 in modifier_pairs:
        model.add_element(interface0)
        model.add_element(interface1)
        model.add_modifier(interface0, elem0, BooleanUnionModifier())
        model.add_modifier(interface1, elem1, BooleanDifferenceModifier())

# =============================================================================
# Labels
# =============================================================================

# Path to font data file
HERE = Path(__file__).parent
DATA = HERE
SESSION = DATA / "text.json"


class Label(Data):
    @property
    def __data__(self) -> dict:
        return {
            "name": self.name,
            "frame": self.frame,
            "polylines": self.polylines,
        }

    def __init__(self, name: str = "", frame: Frame = Frame.worldXY(), polylines=None):
        super().__init__(name=name)

        self.frame = frame
        self.polylines: list[Polyline] = polylines or []

        # Store the file path for later lazy loading
        self._session_path = SESSION
        self._session = None
        self._font_data = None
        self._char_map = {}

    def _load_session(self):
        """Lazy load the session data when needed."""
        if self._session is None:
            self._session = compas.json_load(self._session_path)

    def _load_font_data(self):
        """Lazy load the font data when needed."""
        if self._font_data is None:
            self._load_session()
            self._font_data = self._session.get("bbfont", {}).get("letter", [])

    def _get_letter(self, char):
        """Get a letter from the font data, loading it if necessary."""
        # Try to get from cache first
        if char in self._char_map:
            return self._char_map[char]

        # If not in cache, try to load it
        self._load_font_data()

        # Check if it's a code point or character
        code = None
        if isinstance(char, int):
            code = char
        elif isinstance(char, str) and len(char) == 1:
            code = ord(char)

        # Find and cache the letter
        for letter in self._font_data:
            if "_code" in letter:
                try:
                    letter_code = int(letter["_code"])
                    if letter_code == code:
                        self._char_map[char] = letter
                        if 32 <= letter_code <= 126:  # Standard ASCII printable range
                            self._char_map[chr(letter_code)] = letter
                        return letter
                except ValueError:
                    # Handle special cases where _code might be a character
                    if letter["_code"] == char:
                        self._char_map[char] = letter
                        return letter

        # Letter not found
        return None

    @classmethod
    def from_string(cls, string, frame, scale=1.0, line_spacing=1.5, letter_spacing=0.0, space_width=0.5):
        """
        Convert a string to a list of polylines based on the font data in the JSON file.

        Parameters
        ----------
        string : str
            The string to convert to polylines
        frame : Frame
            The frame to place the text in
        scale : float, optional
            Scale factor for the letters, defaults to 1.0
        line_spacing : float, optional
            Factor for vertical spacing between lines, defaults to 1.5
        letter_spacing : float, optional
            Additional spacing factor between letters, defaults to 0.0 for zero spacing
        space_width : float, optional
            Width of space character. Defaults to 0.5.

        Returns
        -------
        list
            List of compas.geometry.Polyline objects.
        """
        label = Label(name=string, frame=frame)
        polylines = []
        current_x = 0
        current_y = 0
        line_height = 1.0

        def get_letter_bounds(letter_data):
            min_x, min_y = float("inf"), float("inf")
            max_x, max_y = float("-inf"), float("-inf")

            path_data = letter_data.get("path", [])
            if not isinstance(path_data, list):
                path_data = [path_data]

            for path in path_data:
                x = float(path.get("_x", 0))
                y = float(path.get("_y", 0))
                min_x, min_y = min(min_x, x), min(min_y, y)
                max_x, max_y = max(max_x, x), max(max_y, y)

                to_points = path.get("to", [])
                if not isinstance(to_points, list):
                    to_points = [to_points]

                for point in to_points:
                    x = float(point.get("_x", 0))
                    y = float(point.get("_y", 0))
                    min_x, min_y = min(min_x, x), min(min_y, y)
                    max_x, max_y = max(max_x, x), max(max_y, y)

            return min_x, min_y, max_x, max_y

        for char in string:
            if char == "\n":
                current_x = 0
                current_y -= line_height * line_spacing * scale
                continue

            if char == " ":
                current_x += space_width * scale
                continue

            letter_data = label._get_letter(char) or label._get_letter(ord(char))

            if not letter_data:
                current_x += 0.3 * scale
                continue

            min_x, min_y, max_x, max_y = get_letter_bounds(letter_data)
            letter_width = max_x - min_x
            letter_height = max_y - min_y

            if letter_height > line_height:
                line_height = letter_height

            start_pos = float(letter_data.get("_start", 0.0))
            letter_x = current_x - start_pos * scale

            path_data = letter_data.get("path", [])
            if not isinstance(path_data, list):
                path_data = [path_data]

            for path in path_data:
                points = []

                start_x = float(path.get("_x", 0)) * scale + letter_x
                start_y = float(path.get("_y", 0)) * scale + current_y
                points.append(Point(start_x, start_y, 0))

                to_points = path.get("to", [])
                if not isinstance(to_points, list):
                    to_points = [to_points]

                for to_point in to_points:
                    x = float(to_point.get("_x", 0)) * scale + letter_x
                    y = float(to_point.get("_y", 0)) * scale + current_y
                    points.append(Point(x, y, 0))

                if len(points) > 1:
                    polylines.append(Polyline(points))

            # Use the _end parameter to determine letter spacing
            # This ensures consistent spacing regardless of the letter's actual width
            end_pos = float(letter_data.get("_end", letter_width))
            current_x += end_pos * scale + letter_spacing * scale

        xform = Transformation.from_frame_to_frame(Frame.worldXY(), label.frame)
        for polyline in polylines:
            polyline.transform(xform)

        label.polylines = polylines

        return label

    def transform(self, xform):
        for polyline in self.polylines:
            polyline.transform(xform)

    def transformed(self, xform):
        new = self.copy()
        new.transform(xform)
        return new

labels = []

for id, block in enumerate(model.elements()):
    if block.name == "Block":
        frame = model.graph.node_attribute(id, "frame")
        labels.append(Label.from_string(str(id), frame.flipped(), 0.1))

# =============================================================================
# Vizualize
# =============================================================================

scene = Scene()

for polygon in polygons:
    scene.add(polygon, name="user_2d_pattern")

# scene.add(mesh)
# scene.add(pattern)
# scene.add(mesh_lscm, name="Mesh", show_points=True)
scene.add(pattern, name="Pattern", show_points=True)

# 3D Blocks
o = Point(5, 5, 0)
for id, block in enumerate(model.elements()):

    if block.name != "Block":
        continue

    geometry = block.modelgeometry

    # Add meshes
    scene.add(geometry, show_lines=True)

    # Add frame
    # scene.add(model.graph.node_attribute(id, "frame"))

    # Add each polyline from the label instead of the label object itself
    transformed_label = labels[id]
    for polyline in transformed_label.polylines:
        scene.add(polyline, color=(255, 0, 0))

# 2D Blocks
x = 0
offset_x = 1
offset_y = -5
for idx, block in enumerate(model.elements()):
    if block.name == "Block":
        geometry = block.modelgeometry.copy()
        frame = model.graph.node_attribute(idx, "frame")

        # Orient mesh to xy frame and move it to the left
        O = Transformation.from_frame_to_frame(frame, Frame.worldXY())
        geometry.transform(O)
        box = geometry.aabb() if isinstance(geometry, Mesh) else geometry.aabb
        xmin, width = abs(box.xmin), box.xsize
        T = Translation.from_vector([x + xmin, offset_y, 0])
        geometry.transform(T)

        # Add geometry
        scene.add(geometry)

        # Add labels
        for polyline in labels[idx].transformed(T * O).polylines:
            scene.add(polyline, color=(255, 0, 0))

        # Update x position
        x += width + offset_x

scene.draw()
