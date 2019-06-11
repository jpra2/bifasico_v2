
from preprocess.load_mesh import Mesh
from bifasico.elems_bifasico import BifasicElems

__all__ = []

mesh = Mesh()
bif_elems = BifasicElems(mesh.data_loaded, mesh.Adjs, mesh.all_centroids, mesh.all_faces_in, mesh.all_kharm, mesh.wirebasket_elems_nv1, mesh.wells_injector, mesh.wells_producer, mesh.tags, mesh.mb, mesh.volumes_d, mesh.volumes_n, mesh.ler_anterior, mesh.mtu)
