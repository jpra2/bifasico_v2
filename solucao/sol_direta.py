from utils.others_utils import OtherUtils as oth

class SolDireta:

    def __init__(self):
        pass

    def solucao_pressao(self, Tf, b):
        self.P = oth.get_solution(Tf, b)

    def calculate_total_flux(self, volumes, faces):

        p_tag = self.mb.tag_get_handle('P')
        volumes_d = self.mb.get_entities_by_type_and_tag(0, types.MBHEX, np.array([p_tag]), np.array([None]))
        ids_vds = self.mb.tag_get_data(self.ids_volumes_tag, volumes_d, flat=True)
        p_vds = self.mb.tag_get_data(self.pf_tag, volumes_d, flat=True)

        mobi_in_faces = self.mb.tag_get_data(self.mobi_in_faces_tag, faces, flat=True)
        # all_gamaf = self.mb.tag_get_data(self.gamaf_tag, faces, flat=True)
        fws_faces = self.mb.tag_get_data(self.fw_in_faces_tag, faces, flat=True)
        ps = self.mb.tag_get_data(self.pf_tag, volumes, flat=True)
        if self.gravity:
            all_sgravs = self.mb.tag_get_data(self.s_grav_tag, faces, flat=True)
        else:
            all_sgravs = np.zeros(len(faces))

        fluxos = np.zeros(len(volumes))
        fluxos_w = fluxos.copy()
        flux_in_faces = np.zeros(len(faces))
        fluxo_grav_volumes = np.zeros(len(volumes))

        for i, face in enumerate(faces):
            # gamaf = all_gamaf[i]
            elems = self.mb.get_adjacencies(face, 3)
            id0 = self.map_volumes[elems[0]]
            id1 = self.map_volumes[elems[1]]
            mobi = mobi_in_faces[i]
            # s_grav = self.gama*mobi*(self.all_centroids[id1][2] - self.all_centroids[id0][2])
            # s_grav = gamaf*mobi*(self.all_centroids[id1][2] - self.all_centroids[id0][2])
            s_grav = all_sgravs[i]
            fw = fws_faces[i]
            flux = (ps[id1] - ps[id0])*mobi
            if self.gravity == True:
                flux += s_grav
                fluxo_grav_volumes[id0] += s_grav
                fluxo_grav_volumes[id1] -= s_grav

            # flux *= -1

            fluxos[id0] += flux
            fluxos_w[id0] += flux*fw
            fluxos[id1] -= flux
            fluxos_w[id1] -= flux*fw
            flux_in_faces[i] = flux

        self.mb.tag_set_data(self.total_flux_tag, volumes, fluxos)
        self.mb.tag_set_data(self.flux_w_tag, volumes, fluxos_w)
        self.mb.tag_set_data(self.flux_in_faces_tag, faces, flux_in_faces)
        self.mb.tag_set_data(self.s_grav_volume_tag, volumes, fluxo_grav_volumes)
