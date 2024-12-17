from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import pandas as pd
import geopandas as gpd
import shapely
import numpy as np

class QiskitShapelyRenderer(QMplRenderer):
    dfs = []

    def render_poly(self,
                    table: pd.DataFrame,
                    ax,
                    subtracted: bool = False,
                    extra_kw: dict = None):
        self.dfs.append(table)

    def get_net_coordinates(self, resolution=4):
        self.dfs = []
        self.options.resolution = str(resolution)
        self.render_tables(None)

        if len(self.dfs) == 0:
            return gpd.GeoDataFrame(columns=["component", "name", "geometry", "layer", "subtract"])
        gsdf = gpd.GeoDataFrame(pd.concat(self.dfs, ignore_index=True))

        #Discard polygons that are simulation constructs (i.e. "junction" types)
        df_jjs = self.design.qgeometry.tables['junction']
        inds_to_pop = []
        for index, row in df_jjs.iterrows():
            for index_cands, row_cands in gsdf.iterrows():
                if row['component'] == row_cands['component'] and row['name'] == row_cands['name']:
                    inds_to_pop += [index_cands]
        gsdf.drop(inds_to_pop, inplace=True)
        

        return gsdf

    @staticmethod
    def get_rendered_path_poly(design, path_coords, path_width, fillet_radius, resolution=4):
        mplRend = QMplRenderer(None, design, None)
        lePath = mplRend.fillet_path({
            'geometry': shapely.LineString(path_coords.tolist() if isinstance(path_coords, np.ndarray) else path_coords),
            'fillet': fillet_radius,
            'resolution': resolution
        })
        return lePath.buffer(distance=path_width/2,
                          cap_style=shapely.geometry.CAP_STYLE.flat,
                          join_style=shapely.geometry.JOIN_STYLE.mitre,
                          resolution=resolution)
