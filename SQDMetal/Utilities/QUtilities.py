# from pint import UnitRegistry
import shapely

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
