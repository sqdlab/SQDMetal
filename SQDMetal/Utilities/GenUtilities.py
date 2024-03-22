
class GenUtilities:
    @staticmethod
    def add_units(val, sigfigs=-1, noUnits=False):
        if isinstance(val, float) or isinstance(val, int):
            thinspace = u"\u2009"
            def clip_val(value):
                if sigfigs > 0:
                    return '{number:.{width}g}'.format(number=value, width=sigfigs)
                else:
                    return f'{value:.12g}'
            
            if noUnits:
                return f'{clip_val(val)}'

            if abs(val) < 1e-15:
                return f'{clip_val(val*1e18)}{thinspace}a'
            if abs(val) < 1e-12:
                return f'{clip_val(val*1e15)}{thinspace}f'
            if abs(val) < 1e-9:
                return f'{clip_val(val*1e12)}{thinspace}p'
            if abs(val) < 1e-6:
                return f'{clip_val(val*1e9)}{thinspace}n'
            if abs(val) < 1e-3:
                return f'{clip_val(val*1e6)}{thinspace}μ'
            if abs(val) < 1:
                return f'{clip_val(val*1e3)}{thinspace}m'
            if abs(val) < 1000:
                return val
            if abs(val) < 1e6:
                return f'{clip_val(val*1e-3)}{thinspace}k'
            if abs(val) < 1e9:
                return f'{clip_val(val*1e-6)}{thinspace}M'

            return f'{clip_val(val*1e-9)}{thinspace}G'
        else:
            return val
