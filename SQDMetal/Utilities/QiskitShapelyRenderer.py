from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import pandas as pd
import geopandas as gpd

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
        return gpd.GeoDataFrame(pd.concat(self.dfs, ignore_index=True))
