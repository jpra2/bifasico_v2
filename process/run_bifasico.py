from process.loading import mesh, bif_elems

if mesh.ADM:
    pass

elif not mesh.ADM:
    from solucao import sol_direta


import pdb; pdb.set_trace()
print('saiu run bifasico \n')
