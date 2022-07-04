import gmsh
from objectgmsh import (
    Model,
    Shape,
    MeshControlConstant,
    MeshControlExponential,
    cut,
    rotate,
)
from objectgmsh.objects import MeshControlLinear


# dimensions
cruc_r = 0.06  # crucible radius
cruc_h = 0.03  # crucible height
cruc_hi = 0.015  # crucible height inside
melt_r = 0.025  # melt radius
melt_h = 0.01  # melt height
crys_r = 0.005  # crystal radius
crys_h = 0.1  # crystal height


occ = gmsh.model.occ
model = Model()

# main bodies
crystal = Shape(
    model, 3, "crystal", [rotate(occ.addRectangle(0, 0, 0, crys_r, crys_h))]
)
melt = Shape(
    model, 3, "melt", [rotate(occ.addRectangle(0, -melt_h, 0, melt_r, melt_h))]
)

crucible = occ.addRectangle(0, -melt_h - (cruc_h - cruc_hi), 0, cruc_r, cruc_h)
crucible_hole = occ.addRectangle(0, -melt_h, 0, melt_r, cruc_hi)
cut([(2, crucible)], [(2, crucible_hole)])
crucible = Shape(model, 3, "crucible", [rotate(crucible)])

# create connection between the shapes
crystal.set_interface(melt)
melt.set_interface(crucible)

# set mesh sizes
crystal.mesh_size = 0.0025
melt.mesh_size = 0.005
crucible.mesh_size = 0.01

# detect boundaries
bnd_crystal_out = Shape(
    model,
    2,
    "bnd_crystal_out",
    [x for x in crystal.boundaries if x not in crystal.get_interface(melt)],
)
bnd_melt = Shape(
    model,
    2,
    "bnd_melt_surf",
    [
        x
        for x in melt.boundaries
        if x not in melt.get_interface(crystal) + melt.get_interface(crucible)
    ],
)
bnd_crucible_bottom = Shape(model, 2, "bnd_crucible_bottom", [crucible.bottom_boundary])
bnd_crucible_outside = Shape(
    model,
    2,
    "bnd_crucible_outside",
    [
        x
        for x in crucible.boundaries
        if x not in crucible.get_interface(melt) + bnd_crucible_bottom.geo_ids
    ],
)
if_crystal_melt = Shape(model, 2, "if_crystal_melt", crystal.get_interface(melt))
if_crucible_melt = Shape(model, 2, "if_crucible_melt", crucible.get_interface(melt))

# add physical groups
model.make_physical()

# set mesh constraints
model.deactivate_characteristic_length()
model.set_const_mesh_sizes()

MeshControlExponential(
    model, if_crystal_melt, 0.001, exp=1.7, fact=2
)  # refinement at crystallization front
MeshControlLinear(model, if_crucible_melt, melt.mesh_size, crucible.mesh_size)

# create mesh, show, export
model.generate_mesh(3)
model.show()
model.write_msh("crystal-growth-3D.msh")
