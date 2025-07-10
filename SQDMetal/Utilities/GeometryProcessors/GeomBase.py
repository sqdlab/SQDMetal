
class GeomBase:
    @property
    def ParserType(self):
        raise NotImplementedError()

    @property
    def chip_size_x(self):
        raise NotImplementedError()

    @property
    def chip_size_y(self):
        raise NotImplementedError()

    @property
    def chip_size_z(self):
        raise NotImplementedError()

    @property
    def chip_centre(self):
        raise NotImplementedError()

    def process_layers(self, metallic_layers, ground_plane, **kwargs):
        raise NotImplementedError()
