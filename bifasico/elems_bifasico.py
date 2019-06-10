import numpy as np
import scipy.sparse as sp
from functions import f1
import time

class BifasicElems:

    def __init__(self, data_loaded, Adjs, all_centroids, all_faces_in, all_kharm, all_volumes, injectors, producers, tags, mb, volumes_d, volumes_n, ler_anterior):
        self.mi_w = float(data_loaded['dados_bifasico']['mi_w'])
        self.mi_o = float(data_loaded['dados_bifasico']['mi_o'])
        self.gama_w = float(data_loaded['dados_bifasico']['gama_w'])
        self.gama_o = float(data_loaded['dados_bifasico']['gama_o'])
        self.Sor = float(data_loaded['dados_bifasico']['Sor'])
        self.Swc = float(data_loaded['dados_bifasico']['Swc'])
        self.nw = float(data_loaded['dados_bifasico']['nwater'])
        self.no = float(data_loaded['dados_bifasico']['noil'])
        self.loops = int(data_loaded['dados_bifasico']['loops'])
        self.total_time = float(data_loaded['dados_bifasico']['total_time'])
        self.gravity = data_loaded['gravity']
        # ler_anterior = data_loaded['ler_anterior']
        self.Adjs = Adjs
        self.tags = tags
        self.all_centroids = all_centroids
        self.all_faces_in = all_faces_in
        self.all_kharm = all_kharm
        self.map_volumes = dict(zip(all_volumes, range(len(all_volumes))))
        self.all_volumes = all_volumes
        self.wells_injector = injectors
        self.wells_producer = producers
        self.ids0 = mb.tag_get_data(tags['ID_reord_tag'], np.array(Adjs[:, 0]), flat=True)
        self.ids1 = mb.tag_get_data(tags['ID_reord_tag'], np.array(Adjs[:, 1]), flat=True)
        self.ids_volsd = mb.tag_get_data(tags['ID_reord_tag'], volumes_d, flat=True)
        self.values_d = mb.tag_get_data(tags['P'], volumes_d, flat=True)
        self.ids_volsn = mb.tag_get_data(tags['ID_reord_tag'], volumes_n, flat=True)
        self.values_n = mb.tag_get_data(tags['Q'], volumes_n, flat=True)

        self.load_sats_ini(mb, tags['SAT'])

        if ler_anterior:
            import pdb; pdb.set_trace()
            pass
        else:
            self.set_lamb()
            self.set_mobi_faces_ini()
        pass
        t0 = time.time()
        self.get_Tf_and_b()
        t1 = time.time()
        print(f'tempo: {t1-t0}')

        import pdb; pdb.set_trace()

    def load_sats_ini(self, mb, sat_tag):
        self.all_sats = mb.tag_get_data(sat_tag, self.all_volumes, flat=True)
        self.all_sats_ant = self.all_sats.copy()

    def pol_interp(self, S):
        # S_temp = (S - self.Swc)/(1 - self.Swc - self.Sor)
        # krw = (S_temp)**(self.nw)
        # kro = (1 - S_temp)**(self.no)
        if S > (1 - self.Sor) and S <= 1:
            krw = 1.0
            kro = 0.0
        elif S < self.Swc and S >= 0:
            krw = 0.0
            kro = 1.0
        else:
            S_temp = (S - self.Swc)/(1 - self.Swc - self.Sor)
            if S_temp < 0 or S_temp > 1:
                print('erro S_temp')
                import pdb; pdb.set_trace()
            krw = (S_temp)**(self.nw)
            kro = (1 - S_temp)**(self.no)

        return krw, kro

    def set_lamb(self):
        """
        seta o lambda
        """
        all_sats = self.all_sats
        all_lamb_w = np.zeros(len(self.all_volumes))
        all_lamb_o = all_lamb_w.copy()
        all_lbt = all_lamb_w.copy()
        all_fw = all_lamb_w.copy()
        all_gamav = all_lamb_w.copy()

        for i, sat in enumerate(all_sats):
            # volume = all_volumes[i]
            krw, kro = self.pol_interp(sat)
            all_lamb_w[i] = krw/self.mi_w
            all_lamb_o[i] = kro/self.mi_o
            all_lbt[i] = all_lamb_o[i] + all_lamb_w[i]
            all_fw[i] = all_lamb_w[i]/float(all_lbt[i])
            gama = (self.gama_w*all_lamb_w[i] + self.gama_o*all_lamb_o[i])/(all_lbt[i])
            all_gamav[i] = gama

        self.all_lamb_w = all_lamb_w
        self.all_lamb_o = all_lamb_o
        self.all_lbt = all_lbt
        self.all_fw = all_fw
        self.all_gamav = all_gamav

        # self.mb.tag_set_data(self.lamb_w_tag, all_volumes, all_lamb_w)
        # self.mb.tag_set_data(self.lamb_o_tag, all_volumes, all_lamb_o)
        # self.mb.tag_set_data(self.lbt_tag, all_volumes, all_lbt)
        # self.mb.tag_set_data(self.fw_tag, all_volumes, all_fw)
        # self.mb.tag_set_data(self.gamav_tag, all_volumes, all_gamav)

    def set_mobi_faces_ini(self):
        lim = 1e-5

        all_lbt = self.all_lbt
        all_centroids = self.all_centroids
        all_fw = self.all_fw
        all_sats = self.all_sats
        all_gamav = self.all_gamav
        all_kharm = self.all_kharm
        all_faces_in = self.all_faces_in
        map_volumes = self.map_volumes
        all_faces_in = self.all_faces_in

        all_mobi_in_faces = np.zeros(len(all_faces_in))
        all_s_gravs = all_mobi_in_faces.copy()
        all_fw_in_face = all_mobi_in_faces.copy()
        all_dfds = all_mobi_in_faces.copy()
        all_gamaf = all_mobi_in_faces.copy()
        Adjs = self.Adjs

        for i, face in enumerate(all_faces_in):
            elems = Adjs[i]
            kharm = all_kharm[i]
            id0 = map_volumes[elems[0]]
            id1 = map_volumes[elems[1]]
            lbt0 = all_lbt[id0]
            lbt1 = all_lbt[id1]
            fw0 = all_fw[id0]
            fw1 = all_fw[id1]
            sat0 = all_sats[id0]
            sat1 = all_sats[id1]
            gama0 = all_gamav[id0]
            gama1 = all_gamav[id1]

            if abs(sat0-sat1) < lim:
                all_dfds[i] = 0.0
            else:
                all_dfds[i] = abs((fw0 - fw1)/(sat0 - sat1))
            if elems[0] in self.wells_injector:
                all_mobi_in_faces[i] = kharm*lbt0
                all_fw_in_face[i] = fw0
                gamaf = gama0
            elif elems[1] in self.wells_injector:
                all_mobi_in_faces[i] = kharm*lbt1
                all_fw_in_face[i] = fw1
                gamaf = gama1
            else:
                all_mobi_in_faces[i] = kharm*(lbt0 + lbt1)/2.0
                all_fw_in_face[i] = (fw0 + fw1)/2.0
                gamaf = (gama0 + gama1)/2.0
            all_s_gravs[i] = gamaf*all_mobi_in_faces[i]*(all_centroids[id1][2] - all_centroids[id0][2])
            all_gamaf[i] = gamaf

        self.all_mobi_in_faces = all_mobi_in_faces
        self.all_s_gravs = all_s_gravs
        self.all_fw_in_face = all_fw_in_face
        self.all_dfds = all_dfds
        self.all_gamaf = all_gamaf

        # self.mb.tag_set_data(self.mobi_in_faces_tag, all_faces_in, all_mobi_in_faces)
        # self.mb.tag_set_data(self.s_grav_tag, all_faces_in, all_s_gravs)
        # self.mb.tag_set_data(self.fw_in_faces_tag, all_faces_in, all_fw_in_face)
        # self.mb.tag_set_data(self.dfds_tag, all_faces_in, all_dfds)
        # self.mb.tag_set_data(self.gamaf_tag, all_faces_in, all_gamaf)

    def get_Tf_and_b(self):
        self.Tf, self.b = f1.get_Tf_and_b(self.all_mobi_in_faces, self.ids0, self.ids1, self.all_volumes, self.all_s_gravs)
        self.Tf2, self.b2 = f1.set_boundary_dirichlet(self.Tf, self.b, self.ids_volsd, self.values_d)
        if len(self.ids_volsn) > 0:
            self.b2 = f1.set_boundary_neuman(b, self.ids_volsn, self.values_n)

    def set_boundary_dirichlet(self, Tf, b, ids_volsd, values):
        Tf2 = Tf.copy().tolil()
        b2 = b.copy()
        t = Tf2.shape[0]
        n = len(ids_volsd)

        Tf2[ids_volsd] = sp.lil_matrix((n, t))
        Tf2[ids_volsd, ids_volsd] = np.ones(n)
        b2[ids_volsd] = values

        return Tf2, b2
