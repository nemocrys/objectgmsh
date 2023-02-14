"""Classes that provide an object oriented way to define your gmsh geometry / mesh."""

from copy import deepcopy
import gmsh
import numpy as np
from objectgmsh.utils import *


factory = gmsh.model.occ
field = gmsh.model.mesh.field


class ObjectgmshError(Exception):
    pass


class GeometryError(Exception):
    pass


class Parameters:
    """Dummy class, used as container to store geometry parameters of
    the shapes.
    """

    def __init__(self):
        self.T_init = 0


class Model:
    """Wrapper for Gmsh model. This class is used to manage the shapes
    and physical groups, and provides high-level access to some major
    functionalities of the gmsh API.
    """

    def __init__(self, name="model"):
        """Create gmsh model.

        Args:
            name (str, optional): Name of the model. Defaults to 'model'.
        """
        self._shapes = []
        self.mesh_restrictions = []
        self.min_field = -1
        self._physical = False
        gmsh.initialize(
            ["-noenv"]
        )  # see https://gitlab.onelab.info/gmsh/gmsh/-/issues/1142 for details about -noenv option
        gmsh.model.add(name)
        gmsh.option.setNumber(
            "Geometry.OCCBoundsUseStl", 1
        )  # for better boundary detection, see https://gitlab.onelab.info/gmsh/gmsh/-/issues/1619

    # TODO allow with / close -> implement __enter__, __exit__
    def close_gmsh(self):
        """Call gmsh.finalize()."""
        gmsh.finalize()

    def __getitem__(self, name):
        """Get shape that is part of this model.

        Args:
            name (str): Name of the shape

        Returns:
            Shape: shape object with shape.name == name.
        """
        for shape in self._shapes:
            if shape.name == name:
                return shape
        raise GeometryError(f"Shape {name} does not exist.")

    def __repr__(self):
        """Get a summary of the model.

        Returns:
            str: A summary of the shapes in the model.
        """
        shapes = [s.name for s in self._shapes]
        return f"Gmsh model created with objectgmsh.\nShapes: {shapes}"

    def show(self):
        """Run gmsh GUI."""
        gmsh.fltk.run()

    def _add_shape(self, shape):
        """Add a shape to the model.

        Args:
            shape (Shape): shape object containing the information about
                           the geometry
        """
        self._shapes.append(shape)

    def _apply_restrictions(self):
        self.min_field = field.add("Min")
        field.setNumbers(
            self.min_field, "FieldsList", [x.field for x in self.mesh_restrictions]
        )
        field.setAsBackgroundMesh(self.min_field)

    def deactivate_characteristic_length(self):
        """Don't take any characteristic length into account. Call this
        when defining the mesh size using MeshControl / the parameter
        mesh_size in the shapes.
        """
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)

    def set_characteristic_length(self, char_length, dimensions=[0]):
        """Set characteristic length for all entities (of a dimension).

        Args:
            char_length (float): characteristic length
            dimensions (list, optional): Dimensions to apply the
                characteristic lenght. Defaults to [0].
        """
        for dim in dimensions:
            gmsh.model.mesh.setSize(gmsh.model.getEntities(dim), char_length)

    def generate_mesh(
        self, dimension=2, order=1, size_factor=1, smoothing=1, optimize=None
    ):
        """Generate the mesh.

        Args:
            dimension (int, optional): Dimension of mesh. Defaults to 2.
            order (int, optional): Element order. Defaults to 1.
            size_factor (int, optional): Increase / decrease mesh size
                by a factor. Defaults to 1.
            smoothing (int, optional): Mesh smoothing. Defaults to 1.
            optimize (str, optional): Mesh optimization algorithm.
                Defaults to None.
        """
        self._apply_restrictions()
        gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", size_factor)
        gmsh.option.setNumber("Mesh.Smoothing", smoothing)
        gmsh.model.mesh.generate(dimension)
        gmsh.model.mesh.setOrder(order)
        if optimize is not None:
            gmsh.model.mesh.optimize(optimize)

    def get_shapes(self, dimension, name=""):
        """Get shapes of certain dimension, optionally filtered by name.

        Args:
            dimension (int): Dimension of shapes
            name (str, optional): Only shapes with this string in their
                                  name.

        Returns:
            list: shape objects with shape.dim = dimension and name in
                  shape.name
        """
        return [s for s in self._shapes if s.dim == dimension and name in s.name]

    def make_physical(self):
        """Convert all shapes into physical groups.

        Raises:
            ObjectgmshError: If this function is called more than once.
        """
        if self._physical:
            raise ObjectgmshError("This model is already physical.")
        for shape in self._shapes:
            shape._make_physical()

    def synchronize(self):
        """Synchronize gmsh geometry kernel."""
        # if self._physical:
        #     raise ObjectgmshError('The model is physical. Synchronizing would break it.')
        factory.synchronize()

    def remove_shape(self, shape, recursive=True):
        """Remove shape from model.

        Args:
            shape (Shape): shape object to be removed
        """
        factory.remove(shape.dimtags, recursive=recursive)
        factory.synchronize()
        self._shapes.remove(shape)

    def write_msh(self, file_name):
        """Write mesh to file.

        Args:
            file_name (str): File name (.msh).
        """
        gmsh.write(file_name)

    def set_const_mesh_sizes(self):
        """Create mesh controls for all shapes where the parameter
        shape.mesh_size is set.
        """
        for shape in self._shapes:
            if shape.mesh_size == 0:
                print(
                    f"Warning: Mesh size = 0 for {shape.name}. Ignoring this shape..."
                )
            else:
                print(shape.name, shape.mesh_size)
                MeshControlConstant(self, shape.mesh_size, shapes=[shape])

    @property
    def symmetry_axis(self):
        sym_ax = []
        for shape in self._shapes:
            if shape.dim == 2:
                sym_ax += shape.get_boundaries_in_box([0, 0], [-1e6, 1e6])
        return sym_ax


class Shape:
    """Wrapper for any kind of shape, that shall be part of the final
    model. Shapes may be 2 or 3D objects, lines or points.
    A shape may consist of multiple gmsh entities of the same dimension.
    """

    def __init__(self, model, dim, name, geo_ids=[]):
        """Create a shape.

        Args:
            model (Model): Gmsh model to which the shape will be added.
            dim (int): Dimension of the shape.
            name (str): Name of the shape.
            geo_ids (list, optional): Geometry IDs of the shape.
                                      They may be added later.
        """
        self.model = model
        self.dim = dim
        self.name = name
        self.geo_ids = deepcopy(geo_ids)
        self.params = Parameters()
        self.ph_id = -1
        self.mesh_size = 0
        self.model._add_shape(self)

    def __iadd__(self, other):
        if type(other) is list:
            self.geo_ids += [x for x in other if x not in self.geo_ids]
        if type(other) is int:
            self.geo_ids.append(other)
        if type(other) is Shape:
            self.geo_ids += [x for x in other.geo_ids if x not in self.geo_ids]
            self.model.remove_shape(other)
        return self

    def __isub__(self, other):
        if type(other) is list:
            self.geo_ids = [x for x in self.geo_ids if x not in other]
        elif type(other) is int and other in self.geo_ids:
            self.geo_ids.remove(other)
        elif type(other) is Shape:
            self.geo_ids = [x for x in self.geo_ids if x not in other.geo_ids]
        else:
            raise ObjectgmshError("substraction doesn't work here")
        return self

    @property
    def geo_id(self):
        """Gmsh geometry id.

        Raises:
            GeometryError: If there are more / less than exactly one id
            assigned to this shape.

        Returns:
            int: Gmsh tag.
        """
        if len(self.geo_ids) == 1:
            return self.geo_ids[0]
        else:
            raise GeometryError(f"This shape has {len(self.geo_ids)} geo ids.")

    @property
    def dimtags(self):
        """Gmsh dimension tags.

        Returns:
            list: Gmsh dim-tags of entities in shape.
        """
        return [(self.dim, x) for x in self.geo_ids]

    @property
    def boundaries(self):
        """Boundaries of the shape (dimension: dimension of shape - 1).

        Returns:
            list: Tags of the external boundaries of the shape.

        Note: Internal boundaries, e.g. a point in between two lines
        that form a shape, are excluded. Only the endpoints of the lines
        are cosidered as boundaries.
        """
        bndry = []
        for geo_id in self.geo_ids:
            bndry += get_boundaries(self.dim, geo_id)
        return [
            x for x in bndry if bndry.count(x) == 1
        ]  # exclude boundaries inside shape

    @property
    def all_boundaries(self):
        """Boundaries of the shape (dimension: dimension of shape - 1).

        Returns:
            list: Tags of the internal and external boundaries of the
                shape.
        """
        bndry = []
        for geo_id in self.geo_ids:
            bndry += get_boundaries(self.dim, geo_id)
        return bndry

    @property
    def bounding_box(self):
        """Get the bounding box of this shape.

        Returns:
            list[float]: [x_min, y_min, z_min, x_max, y_max, z_max]
        """
        boxes = np.array(
            [factory.getBoundingBox(self.dim, tag) for tag in self.geo_ids]
        )
        return [
            boxes[:, 0].min(),
            boxes[:, 1].min(),
            boxes[:, 2].min(),
            boxes[:, 3].max(),
            boxes[:, 4].max(),
            boxes[:, 5].max(),
        ]

    def set_interface(self, shape):
        """Get to know the other shape and remove duplicate boundaries.
        Only required if shapes are in contact with each other (calls
        fragment function). May be problematic if there are three shapes
        sharing one point, see
        https://gitlab.onelab.info/gmsh/gmsh/-/issues/2100. The safe way
        is to set all interfaces at once using the function 
        set_interfaces.

        Args:
            shape (Shape): Other shape that is in contact with this.
        """
        factory.fragment(self.dimtags, shape.dimtags)
        factory.synchronize()

    def set_interfaces(self, shapes):
        """Get to know the other shapes and remove duplicate boundaries.
        Only required if shapes are in contact with each other (calls
        fragment function). See also set_interface function.

        Args:
            shape (list): List of other shapes that are in contact with
            this.
        """
        dimtags = []
        for shape in shapes:
            dimtags += shape.dimtags
        factory.fragment(self.dimtags, dimtags)
        factory.synchronize()

    def get_interface(self, shape):
        """Get boundaries that lay in between two shapes.

        Args:
            shape (Shape): Other shape.

        Returns:
            list: Tags of boundary elements.
        """
        own = self.boundaries
        other = shape.boundaries
        return [x for x in own if x in other]

    def get_boundaries_in_box(self, x, y, z=[0, 0], eps=1e-6, one_only=False):
        """Get boundaries of the shape with a box-select. Only
        boundaries that are completely inside the box are returned.

        Args:
            x (list): [x_min, x_max] limits for x-coordinate
            y (list): [y_min, y_max] limits for y-coordinate
            z (list, optional): [z_min, z_max] limits for z coordinate.
            eps (float, optional): Sensitivity. Defaults to 1e-6.
            one_only (bool, optional): Search for only one boundary;
            raise GeometryError if more / less boundaries are found.

        Returns:
            list: Boundary tags.
            Integer return if only_one was set to true.
        """
        dimtags = gmsh.model.getEntitiesInBoundingBox(
            x[0] - eps,
            y[0] - eps,
            z[0] - eps,
            x[1] + eps,
            y[1] + eps,
            z[1] + eps,
            self.dim - 1,
        )
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.boundaries]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(
                    f"Found {len(tags_filtered)} instead of only one boundary."
                )
            else:
                return tags_filtered[0]
        else:
            return tags_filtered

    def get_part_in_box(self, x, y, z=[0, 0], eps=1e-6, one_only=False):
        """Get part of the shape within a box. Only the parts completely
        inside the box are returned

        Args:
            x (list): [x_min, x_max] limits for x-coordinate
            y (list): [y_min, y_max] limits for y-coordinate
            z (list, optional): [z_min, z_max] limits for z coordinate.
            eps (float, optional): Sensitivity. Defaults to 1e-6.
            one_only (bool, optional): Search for only one part;
            raise GeometryError if more / less parts are found.

        Returns:
            list: Geometry tags.
            Integer return if only_one was set to true.
        """
        dimtags = gmsh.model.getEntitiesInBoundingBox(
            x[0] - eps,
            y[0] - eps,
            z[0] - eps,
            x[1] + eps,
            y[1] + eps,
            z[1] + eps,
            self.dim,
        )
        tags = [x[1] for x in dimtags]
        tags_filtered = [x for x in tags if x in self.geo_ids]
        if one_only:
            if len(tags_filtered) != 1:
                raise GeometryError(
                    f"Found {len(tags_filtered)} instead of only one part."
                )
            else:
                return tags_filtered[0]
        else:
            return tags_filtered

    @property
    def top_boundary(self):
        """Boundary with highest y coordinate.

        Returns:
            int: Tag of the boundary.
        """
        [x_min, _, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_max], [y_max, y_max], [z_min, z_max], one_only=True
        )

    @property
    def bottom_boundary(self):
        """Boundary with lowest y coordinate.

        Returns:
            int: Tag of the boundary.
        """
        [x_min, y_min, z_min, x_max, _, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_max], [y_min, y_min], [z_min, z_max], one_only=True
        )

    @property
    def left_boundary(self):
        """Boundary with lowest y coordinate.

        Returns:
            int: Tag of the boundary.
        """
        [x_min, y_min, z_min, _, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_min, x_min], [y_min, y_max], [z_min, z_max], one_only=True
        )

    @property
    def right_boundary(self):
        """Boundary with highest x coordinate.

        Returns:
            int: Tag of the boundary.
        """
        [_, y_min, z_min, x_max, y_max, z_max] = self.bounding_box
        return self.get_boundaries_in_box(
            [x_max, x_max], [y_min, y_max], [z_min, z_max], one_only=True
        )

    @property
    def left_boundary_2(self):
        """Boundary with lowest x coordinate (alternative implementation).
        
        Returns:
            int: Tag of the boundary
        """
        bnd = self.boundaries
        x_min = np.array([factory.get_bounding_box(self.dim - 1, b)[0] for b in bnd])
        cnt = np.count_nonzero(x_min == x_min.min())
        if cnt == 1:
            return bnd[np.argmin(x_min)]
        elif cnt == 0:
            raise ObjectgmshError("Didn't find a matching boundary.")
        else:
            raise ObjectgmshError("Found more than one boundary.")

    @property
    def right_boundary_2(self):
        """Boundary with highest x coordinate (alternative implementation).
        
        Returns:
            int: Tag of the boundary
        """
        bnd = self.boundaries
        x_max = np.array([factory.get_bounding_box(self.dim - 1, b)[3] for b in bnd])
        cnt = np.count_nonzero(x_max == x_max.max())
        if cnt == 1:
            return bnd[np.argmax(x_max)]
        elif cnt == 0:
            raise ObjectgmshError("Didn't find a matching boundary.")
        else:
            raise ObjectgmshError("Found more than one boundary.")

    def set_characteristic_length(self, char_length):
        """Set caracteristic length recursively on all boundaries and
        their boundaries.

        Args:
            char_length (float): Characteristic length for the mesh
                                 generation.
        """
        boundary = gmsh.model.getBoundary(self.dimtags, False, False, True)
        gmsh.model.mesh.setSize(boundary, char_length)

    def _make_physical(self):
        """Convert shape into physical group."""
        self.ph_id = add_physical_group(self.dim, self.geo_ids, self.name)


class MeshControl:
    """Base class for mesh restrictions."""

    def __init__(self, model):
        self._field = -1
        self._restricted_field = -1
        self.faces_list = []
        self.edges_list = []
        self.volumes_list = []
        model.mesh_restrictions.append(self)

    def restrict_to_shapes(self, shapes):
        """Apply the mesh control to the given shapes only.

        Args:
            shapes (list): List of Shape objects.
        """
        for shape in shapes:
            if shape.dim == 3:
                self.volumes_list += shape.geo_ids
                self.faces_list += shape.all_boundaries
                for face in shape.all_boundaries:
                    dim_tags = gmsh.model.getBoundary([(2, face)], False, False, False)
                    for dim_tag in dim_tags:
                        self.edges_list.append(dim_tag[1])
            if shape.dim == 2:
                self.faces_list += shape.geo_ids
                self.edges_list += shape.all_boundaries
            if shape.dim == 1:
                self.edges_list += shape.geo_ids

    def restrict_to_faces(self, faces):
        """Apply the mesh control to the given faces only.

        Args:
            faces (list): List of faces for restriction (gmsh id).
        """
        self.faces_list += faces

    def restrict_to_edges(self, edges):
        """Apply the mesh control to the given edges only.

        Args:
            edges (list): List of edges for restriction (gmsh id).
        """
        self.edges_list += edges

    def restrict_to_volumes(self, volumes):
        """Apply the mesh control to the given volumes only.

        Args:
            volumes (list): List of volumes for restriction (gmsh id).
        """
        self.volumes_list += volumes

    def restrict(self, shapes, faces, edges, volumes):
        """Apply the mesh control to given parts of the model only

        Args:
            shapes (list): List of Shape objects.
            faces (list): List of faces for restriction (gmsh id).
            edges (list): List of edges for restriction (gmsh id).
            volumes (list): List of volumes for restriction (gmsh id).
        """
        self.restrict_to_shapes(shapes)
        self.restrict_to_faces(faces)
        self.restrict_to_edges(edges)
        self.restrict_to_volumes(volumes)

    @property
    def field(self):
        if self.faces_list == [] and self.edges_list == []:
            if self._field == -1:
                raise ObjectgmshError("Field not set.")
            return self._field
        else:
            self._restricted_field = field.add("Restrict")
            field.setNumber(self._restricted_field, "IField", self._field)
            field.setNumbers(self._restricted_field, "FacesList", self.faces_list)
            field.setNumbers(self._restricted_field, "EdgesList", self.edges_list)
            field.setNumbers(self._restricted_field, "VolumesList", self.volumes_list)
            return self._restricted_field


class MeshControlConstant(MeshControl):
    """Mesh control region with constant mesh size."""

    def __init__(
        self, model, char_length, shapes=[], surfaces=[], edges=[], volumes=[]
    ):
        """Create a mesh control of constant size.

        Args:
            model (Model): Objectgmsh model
            char_length (float): Characteristic length.
            shapes (list, optional): Shapes, to which the mesh control
                is restricted. Defaults to [].
            surfaces (list, optional): Surfaces (gmsh id) to which the
                mesh control is restricted. Defaults to [].
            edges (list, optional): Edges (gmsh id) to which the mesh
                control is restricted. Defaults to [].
            volumes (list, optional): Volumes (gmsh id) to which the
                mesh control is restricted. Defaults to [].
        """
        super().__init__(model)
        self._field = field.add("MathEval")
        field.setString(self._field, "F", str(char_length))
        self.restrict(shapes, surfaces, edges, volumes)


class MeshControlLinear(MeshControl):
    """Mesh control region with linear increasing mesh size."""

    def __init__(
        self,
        model,
        shape,
        min_char_length,
        max_char_length,
        dist_start=0,
        dist_end=None,
        NNodesByEdge=1000,
        shapes=[],
        surfaces=[],
        edges=[],
        volumes=[],
    ):
        """Create a mesh control of linear increasing size.

        Args:
            model (Model): Objectgmsh model.
            shape (Shape): Shape that is the starting point of the
                linearly increasing mesh size.
            min_char_length (float): Minimum mesh size
            max_char_length (float): Maximum mesh size
            dist_start (float, optional): Distance from base shape where
                mesh size starts growing. Defaults to 0.
            dist_end (float, optional): Distance from base shape where
                mesh size stops growing. Defaults to None.
            NNodesByEdge (int, optional): Sampling rate.
                Defaults to 1000.
            shapes (list, optional): Shapes, to which the mesh control
                is restricted. Defaults to [].
            surfaces (list, optional): Surfaces (gmsh id) to which the
                mesh control is restricted. Defaults to [].
            edges (list, optional): Edges (gmsh id) to which the mesh
                control is restricted. Defaults to [].
            volumes (list, optional): Volumes (gmsh id) to which the
                mesh control is restricted. Defaults to [].
        """
        super().__init__(model)

        if dist_end is None:
            dist_end = min_char_length + max_char_length

        if shape.dim == 1:
            edg = shape.geo_ids
        elif shape.dim == 2:
            edg = shape.boundaries
        elif shape.dim == 3:
            srf = shape.boundaries

        dist_field = field.add("Distance")
        field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
        if shape.dim == 3:
            field.setNumbers(dist_field, "FacesList", srf)
        else:
            field.setNumbers(dist_field, "EdgesList", edg)
        self._field = field.add("Threshold")
        field.setNumber(self._field, "IField", dist_field)
        field.setNumber(self._field, "LcMin", min_char_length)
        field.setNumber(self._field, "LcMax", max_char_length)
        field.setNumber(self._field, "DistMin", dist_start)
        field.setNumber(self._field, "DistMax", dist_end)
        self.restrict(shapes, surfaces, edges, volumes)


class MeshControlExponential(MeshControl):
    """Mesh control region with exponentially increasing mesh size."""

    def __init__(
        self,
        model,
        shape,
        char_length,
        exp=1.8,
        fact=1,
        NNodesByEdge=1000,
        shapes=[],
        surfaces=[],
        edges=[],
        volumes=[],
    ):
        """Create a mesh control of exponentially increasing size:
        mesh_size = char_length + fact * distance ^ exp

        Args:
            model (Model): Objectgmsh model.
            shape (Shape): Shape that is the starting point of the
                linearly increasing mesh size.
            char_length (float): Characteristic length at zero distance.
            exp (float, optional): Exponent for distance. Defaults to 1.8.
            fact (float, optional): Factor in front of distance.
                Defaults to 1.
            NNodesByEdge (int, optional): Sampling rate.
                Defaults to 1000.
            shapes (list, optional): Shapes, to which the mesh control
                is restricted. Defaults to [].
            surfaces (list, optional): Surfaces (gmsh id) to which the
                mesh control is restricted. Defaults to [].
            edges (list, optional): Edges (gmsh id) to which the mesh
                control is restricted. Defaults to [].
            volumes (list, optional): Volumes (gmsh id) to which the
                mesh control is restricted. Defaults to [].
        """
        super().__init__(model)

        if shape.dim == 1:
            edg = shape.geo_ids
        elif shape.dim == 2:
            edg = shape.boundaries
        elif shape.dim == 3:
            srf = shape.boundaries

        dist_field = gmsh.model.mesh.field.add("Distance")
        field.setNumber(dist_field, "NNodesByEdge", NNodesByEdge)
        if shape.dim == 3:
            field.setNumbers(dist_field, "FacesList", srf)
        else:
            field.setNumbers(dist_field, "EdgesList", edg)
        self._field = field.add("MathEval")
        field.setString(self._field, "F", f"F{dist_field}^{exp}*{fact}+{char_length}")
        self.restrict(shapes, surfaces, edges, volumes)
