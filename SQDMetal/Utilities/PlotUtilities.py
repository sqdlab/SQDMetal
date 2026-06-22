import numpy as np
import matplotlib

class HistogramEqualisation:
    '''
    Calculates the colour-map scaling and color-bar ticks for 2D pcolor plots in matplotlib using Histogram Equalisation.
    '''
    def __init__(self, values, centre_zero=True, num_stops=64):
        self._zMin, self._zMax = np.min(values), np.max(values)
        self._quantileStops = np.linspace(0,1,num_stops)
        self._centre_zero = centre_zero
        if centre_zero:
            v_pos = values[values >= 0]
            v_neg = values[values < 0]
            if v_pos.size > 0:
                self._q_pos = np.quantile(v_pos, self._quantileStops)
            if v_neg.size > 0:
                self._q_neg = np.quantile(v_neg, self._quantileStops)

            if v_pos.size > 0:
                self._zinterpInv_pos = lambda x: np.interp(x, self._quantileStops/2+0.5, self._q_pos)
                self._zinterp_pos = lambda x: np.interp(x, self._q_pos, self._quantileStops/2+0.5)
            else:
                self._zinterpInv_pos = lambda x: 0.0
                self._zinterp_pos = lambda x: 0.5
            #
            if v_neg.size > 0:
                self._zinterpInv_neg = lambda x: np.interp(x, self._quantileStops/2, self._q_neg)
                self._zinterp_neg = lambda x: np.interp(x, self._q_neg, self._quantileStops/2)
            else:
                self._zinterpInv_neg = lambda x: 0.0
                self._zinterp_neg = lambda x: 0.5
        else:
            self._q = np.quantile(values, self._quantileStops)
            self._zinterpInv = lambda x: np.interp(x, self._quantileStops, self._q)
            self._zinterp = lambda x: np.interp(x, self._q, self._quantileStops)

    def zinterpInv(self,x): #Converts from [0,1] to data (0.5 mapping to zero)
        init_shape = x.shape
        x = np.ndarray.flatten(x)
        if self._centre_zero:
            pos = x >= 0.5
            neg = ~pos
            out = np.zeros_like(x)
            if np.any(pos):
                out[pos] = self._zinterpInv_pos(x[pos])
            if np.any(neg):
                out[neg] = self._zinterpInv_neg(x[neg])
        else:
            out = self._zinterpInv(x)
        return out.reshape(init_shape)
    def zinterp(self,x):    #Converts from data to [0,1] (zero mapping to 0.5)
        init_shape = x.shape
        x = np.ndarray.flatten(x)
        if self._centre_zero:
            pos = x >= 0
            neg = ~pos
            out = np.zeros_like(x)
            if np.any(pos):
                out[pos] = self._zinterp_pos(x[pos])
            if np.any(neg):
                out[neg] = self._zinterp_neg(x[neg])
        else:
            out = self._zinterp(x)
        return out.reshape(init_shape)
    
    def cmap(self):
        return matplotlib.colors.FuncNorm((self.zinterp,self.zinterpInv), vmin=self._zMin, vmax=self._zMax)
    
    def cTicks(self, num_ticks:int=3):
        assert num_ticks > 0 and num_ticks % 2 == 1, "num_ticks be an odd natural number."
        if self._centre_zero:
            stops = np.concatenate((np.linspace(0, 0.5, int((num_ticks+1)/2)), np.linspace(0.5, 1.0, int((num_ticks+1)/2))[1:]))
            return self.zinterpInv(stops)
        else:
            return self.zinterpInv(np.linspace(0, 1.0, num_ticks))
