import gmsh
from objectgmsh import Model, Shape, MeshControlConstant, MeshControlExponential, cut, make_wedge, get_wedge_boundaries


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
crystal = make_wedge(Shape(model, 2, "crystal", [occ.addRectangle(0, 0, 0, crys_r, crys_h)]))
melt = make_wedge(Shape(model, 2, "melt", [occ.addRectangle(0, -melt_h, 0, melt_r, melt_h)]))

crucible = occ.addRectangle(0, -melt_h - (cruc_h - cruc_hi), 0, cruc_r, cruc_h)
crucible_hole = occ.addRectangle(0, -melt_h, 0, melt_r, cruc_hi)
cut([(2, crucible)], [(2, crucible_hole)])
crucible = make_wedge(Shape(model, 2, "crucible", [crucible]))

# create connection between the shapes
crystal.set_interface(melt)
melt.set_interface(crucible)

# detect boundaries
wedge_boundaries = get_wedge_boundaries()
bnd_crystal_out = Shape(
    model,
    2,
    "bnd_crystal_out",
    [x for x in crystal.boundaries if x not in crystal.get_interface(melt) + wedge_boundaries],
)
bnd_melt = Shape(
    model,
    2,
    "bnd_melt_surf",
    [
        x
        for x in melt.boundaries
        if x not in melt.get_interface(crystal) + melt.get_interface(crucible) + wedge_boundaries
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
        if x not in crucible.get_interface(melt) + bnd_crucible_bottom.geo_ids + wedge_boundaries
    ],
)
if_crystal_melt = Shape(model, 2, "if_crystal_melt", crystal.get_interface(melt))
if_crucible_melt = Shape(model, 2, "if_crucible_melt", crucible.get_interface(melt))

# add physical groups
model.make_physical()

# set mesh constraints
model.deactivate_characteristic_length()
MeshControlConstant(model, 0.005, [crucible, melt])
MeshControlConstant(model, 0.0025, [crystal])
MeshControlExponential(
    model, if_crystal_melt, 0.001, exp=1.7, shapes=[crystal, melt, crucible]
)

# create mesh, show, export
model.generate_mesh()
model.show()
model.write_msh("crystal-growth-2D-wedge.msh")
