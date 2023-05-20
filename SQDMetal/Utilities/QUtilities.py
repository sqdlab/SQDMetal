# from pint import UnitRegistry
import shapely
import numpy as np

class QUtilities:
    @staticmethod
    def get_units(design):
        unit_conv = design.get_units()
        if unit_conv == 'mm':
            unit_conv = 1e-3
        elif unit_conv == 'um':
            unit_conv = 1e-6
        elif unit_conv == 'nm':
            unit_conv = 1e-9
        else:
            assert False, f"Unrecognised units: {unit_conv}"
        return unit_conv

    @staticmethod
    def parse_value_length(strVal):
        #This is far too slow!?
        # ureg = UnitRegistry().Quantity
        # return ureg(strVal).to('m').magnitude

        #So do this instead...
        strVal = strVal.strip()
        assert len(strVal) > 1, f"Length \'{strVal}\' is invalid."
        if strVal[-2:] == "mm":
            return float(strVal[:-2]+'e-3')
        elif strVal[-2:] == "um":
            return float(strVal[:-2]+'e-6')
        elif strVal[-2:] == "nm":
            return float(strVal[:-2]+'e-9')
        elif strVal[-2:] == "pm":
            return float(strVal[:-2]+'e-12')
        elif strVal[-2:] == "fm":
            return float(strVal[:-2]+'e-15')
        elif strVal[-1:] == "m":
            return float(strVal[:-2])
        else:
            assert len(strVal) > 1, f"Length \'{strVal}\' is invalid."

    @staticmethod
    def chk_within(chk_geom, main_geom, thresh=0.99):
        return shapely.intersection(chk_geom, main_geom).area / chk_geom.area > thresh

    @staticmethod
    def get_comp_bounds(design, objs, units_metres = False):
        if not (isinstance(objs, list) or isinstance(objs, np.ndarray)):
            objs = [objs]
        x_vals = []
        y_vals = []
        for cur_obj in objs:
            paths = design.components[cur_obj].qgeometry_table('path')
            for _, row in paths.iterrows():
                cur_minX, cur_minY, cur_maxX, cur_maxY = row['geometry'].buffer(row['width'] / 2, cap_style=shapely.geometry.CAP_STYLE.flat).bounds
                x_vals += [cur_minX, cur_maxX]
                y_vals += [cur_minY, cur_maxY]
            for cur_poly in design.components[cur_obj].qgeometry_list('poly'):
                cur_minX, cur_minY, cur_maxX, cur_maxY = cur_poly.bounds
                x_vals += [cur_minX, cur_maxX]
                y_vals += [cur_minY, cur_maxY]
        if units_metres:
            unit_conv = QUtilities.get_units(design)
            return unit_conv*min(x_vals), unit_conv*min(y_vals), unit_conv*max(x_vals), unit_conv*max(y_vals)
        else:
            return min(x_vals), min(y_vals), max(x_vals), max(y_vals)
