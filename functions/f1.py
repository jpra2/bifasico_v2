import scipy.sparse as sp
import numpy as np


def get_Tf_and_b(mobi_in_faces, ids0, ids1, volumes, s_grav_f):
    '''
    mobi_in_faces = mobilidade das faces internas
    ids0 = adjs[:,0]
    ids1 = adjs[:,1]
    volumes = volumes locais
    '''
    lines = []
    cols = []
    data = []
    n = len(volumes)

    lines.append(ids0)
    cols.append(ids1)
    data.append(mobi_in_faces)
    lines.append(ids1)
    cols.append(ids0)
    data.append(mobi_in_faces)
    lines.append(ids0)
    cols.append(ids0)
    data.append(-mobi_in_faces)
    lines.append(ids1)
    cols.append(ids1)
    data.append(-mobi_in_faces)
    lines = np.concatenate(lines)
    cols = np.concatenate(cols)
    data = np.concatenate(data)

    Tf = sp.csc_matrix((data, (lines, cols)), shape=(n, n))

    s_grav_volumes = np.zeros(n)

    lines = []
    cols = []
    data = []

    n2 = len(s_grav_f)

    lines.append(ids0)
    cols.append(np.zeros(n2))
    data.append(s_grav_f)
    lines.append(ids1)
    cols.append(np.zeros(n2))
    data.append(-s_grav_f)
    lines = np.concatenate(lines)
    cols = np.concatenate(cols)
    data = np.concatenate(data)

    b = np.array(sp.csc_matrix((data, (lines, cols)), shape=(n,1)).todense()).flatten()

    return Tf, b

def set_boundary_dirichlet(Tf, b, ids_volsd, values):
    Tf2 = Tf.copy().tolil()
    b2 = b.copy()
    t = Tf2.shape[0]
    n = len(ids_volsd)

    Tf2[ids_volsd] = sp.lil_matrix((n, t))
    Tf2[ids_volsd, ids_volsd] = np.ones(n)
    b2[ids_volsd] = values

    return Tf2, b2

def set_boundary_neuman(b, ids_volsn, values):
    n = len(ids_volsn)
    b[ids_volsn, np.zeros(n, dtype=np.int32)] = values
    return b
