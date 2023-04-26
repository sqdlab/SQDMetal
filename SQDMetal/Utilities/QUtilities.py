from pint import UnitRegistry

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
        ureg = UnitRegistry().Quantity
        return ureg(strVal).to('m').magnitude
