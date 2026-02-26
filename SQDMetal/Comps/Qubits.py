# -*- coding: utf-8 -*-

# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Author: Divita Gautam
# Creation Date: 21/05/2024
# Description: Collection of classes to draw qubits.

import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit
from shapely.geometry.base import CAP_STYLE
from shapely.geometry import Polygon
from SQDMetal.Utilities.QUtilities import QUtilities 
from shapely.geometry import LineString

class TransmonTapered(BaseQubit):
    """Transmon pocket with 'Teeth' connection pads.

    Inherits `BaseQubit` class

    Description:
        Create a standard pocket transmon qubit for a ground plane with teeth
        Here we use the 'Teeth' shape which ones connected to the top pad and one connection pad.

    Options:
        Convention: Values (unless noted) are strings with units included,
        (e.g., '30um')

    Pocket:
        * pad_gap            - the distance between the two charge islands, which is also the
          resulting 'length' of the pseudo junction
        * inductor_width     - width of the pseudo junction between the two charge islands (should be same as josephson junction width)
        * inductor_height    - height of the pseudo junction between the two charge islands 
        * pad_width          - the width (x-axis) of the charge island pads, except the circle radius from both sides
        * pad_height         - the size (y-axis) of the charge island pads
        * pocket_width       - size of the pocket (cut out in ground) along x-axis
        * pocket_height      - size of the pocket (cut out in ground) along y-axis
        * fillet_radius_gap  - The radius of the curved edges of the pocket
        * fillet_resolution  - number of points used to calculate the fillet/curve
        * chrgln_pin_x_offset =     - horizontal distance of the pin from the qubit center
        * chrgln_pin_y_offset =     - vertical distance of the pin from the pocket edge
        * coupled_pad_gap    - the distance between the two teeth shape
        * coupled_pad_width  - the width (x-axis) of the teeth shape on the island pads
        * coupled_pad_height - the size (y-axis) of the teeth shape on the island pads
                            

    Connector lines:
        * pad_gap        - space between the connector pad and the charge island it is
          nearest to
        * pad_width      - width (x-axis) of the connector pad
        * pad_height     - height (y-axis) of the connector pad
        * pad_cpw_shift  - shift the connector pad cpw line by this much away from qubit
        * pad_cpw_extent - how long should the pad be - edge that is parallel to pocket
        * cpw_width      - center trace width of the CPW line
        * cpw_gap        - dielectric gap width of the CPW line
        * cpw_extend     - depth the connector line extends into ground (past the pocket edge)
        * pocket_extent  - How deep into the pocket should we penetrate with the cpw connector
        * fillet_radius - The radius of the curved edges of the connector pad
          (into the ground plane)
        * pocket_rise    - How far up or down relative to the center of the transmon should we
          elevate the cpw connection point on the ground plane
        * loc_W / H      - which 'quadrant' of the pocket the connector is set to, +/- 1 (check
          if diagram is correct)


    Sketch:
        Below is a sketch of the qubit
        ::

                 +1              0             +1
                _________________________________
            -1  |                |               |  +1      Y
                |           | | |_| | |          |          ^
                |        ___| |_____| |____      |          |
                |       /     island       \     |          |----->  X
                |       \__________________/     |
                |             |       |          |
                |               |   |            |
                |                |_|             |
                |                |               |
                |                |               |
                |  pocket        x               |                                              
                |                |               |  
                |                ||              |    
                |               |  |             |                      
                |        _____|______|_____      |
                |       /                  \     |
                |       \__________________/     |
                |                                |
                |                                |
            -1  |________________________________|                          
                            

    .. image::
        transmon_pocket_teeth.png

    .. meta::
        Transmon Pocket Teeth

    """

    #_img = 'transmon_pocket1.png'

    # Default drawing options
    default_options = Dict(
        pad_gap='100um',
        inductor_width='20um',
        inductor_height='5um',
        pad_width='400um',
        pad_height='90um',
        pocket_width='650um',
        pocket_height='650um',
        fillet_radius_gap='30um',
        fillet_resolution=4,
        chrgln_pin_x_offset = '30um',  # User-defined horizontal distance from the qubit center
        chrgln_pin_y_offset = '50um', # User-defined vertical distance from the pocket edge
        #Taperred part of qubit
        taper_width_top='100um',
        taper_width_base='150um',
        taper_height='30um',
        fillet_radius='1um',     
        # coupled_pad belongs to the teeth part. Teeth will have same height/width and are symmetric.
        coupled_pad_height='150um',
        coupled_pad_width='20um',
        coupled_pad_gap='50um',  # One can arrange the gap between the teeth.
        # orientation = 90 has dipole aligned along the +X axis, while 0 aligns to the +Y axis
        _default_connection_pads=Dict(
            pad_gap='15um',
            pad_width='20um',
            pad_height='150um',
            pad_cpw_shift='0um',
            pad_cpw_extent='25um',
            cpw_width='10um',
            cpw_gap='6um',
            # edge_curve=0.1,
            # : cpw_extend: how far into the ground to extend the CPW line from the coupling pads
            cpw_extend='100um',
            pocket_extent='5um',
            pocket_rise='0um',
            fillet_radius='1um',
            fillet_resolution=4,
            pin_y_distance='50um',
            loc_W='+1',  # width location  only +-1 or 0,
            loc_H='+1',  # height location  only +-1 or 0
        ))
    # """Default drawing options"""

    component_metadata = Dict(short_name='Pocket',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')
    # """Component metadata"""

    # TOOLTIP = """Transmon pocket with teeth pads."""

    # def make(self):
    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    **kwargs):
        super().__init__(design, name, options, **kwargs)
        
        # Compute and validate slope
        taper_width_top = float(options.get('taper_width_top', '100um')[:-2])
        taper_width_base = float(options.get('taper_width_base', '150um')[:-2])
        taper_height = float(options.get('taper_height', '30um')[:-2])

        computed_slope = (2 * taper_height) / (taper_width_base - taper_width_top)
        assert computed_slope < 0.4, f"Slope {computed_slope:.2f} is too steep! Adjust taper dimensions for the slope to be less than 0.4."
    def make(self):   
        self.make_pocket()
        self.make_connection_pads()

    #Function to create trapered trapezoid
    def create_trapezoid(self, center_x, top_width, base_width, height, y_offset, rfillet, startx):
        half_top = float(top_width) / 2
        half_base = float(base_width) / 2
        coordinates = [
            (center_x - half_base, y_offset),
            (center_x - half_top, y_offset + height),
            (center_x + half_top, y_offset + height),
            (center_x + half_base, y_offset)
        ]
        
        # Create a Polygon object from the coordinates
        trapezoid_polygon = Polygon(coordinates)
        
        # Convert the coordinates to numpy array
        raw_data = np.array(trapezoid_polygon.exterior.coords)

        start_x = center_x - startx
        start_y = y_offset

        end_x = center_x + startx
        end_y = y_offset

        start = [[start_x, start_y]]
        end = [[end_x, end_y]]

        # Concatenate arrays after converting start and end points to 2D arrays
        raw_data_2 = np.vstack((start, raw_data, end))

        raw_data_1=np.delete(raw_data_2,5,0)

        # Curve the edges of the polygon
        data = QUtilities.calc_filleted_path(raw_data_1, rfillet, 9)
        
        # Create a new Polygon object from the curved coordinates
        trapezoid_with_curved_edges = Polygon(data)

        return trapezoid_with_curved_edges
    
        # """Makes standard transmon in a pocket."""

    def make_pocket(self):
            
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        #  pcop = self.p.coupled_pads[name]  # parser on connector options

        # since we will reuse these options, parse them once and define them as variables
        pad_width = p.pad_width
        pad_height = p.pad_height
        pad_gap = p.pad_gap
        coupled_pad_height = p.coupled_pad_height
        coupled_pad_width = p.coupled_pad_width
        coupled_pad_gap = p.coupled_pad_gap

        # make the pads as rectangles (shapely polygons)
        pad = draw.rectangle(pad_width, pad_height)

        pad_top = draw.translate(pad, 0, +(pad_height + pad_gap) / 2.)
        # Here, you make your pads round. Not sharp shape on the left and right sides and also this should be the same for the bottom pad as the top pad.
        circ_left_top = draw.Point(-pad_width / 2., +(pad_height + pad_gap) /
                                   2.).buffer(pad_height / 2,
                                              resolution=16,
                                              cap_style=CAP_STYLE.round)
        circ_right_top = draw.Point(pad_width / 2., +(pad_height + pad_gap) /
                                    2.).buffer(pad_height / 2,
                                               resolution=16,
                                               cap_style=CAP_STYLE.round)
        # In here you create the teeth part and then you union them as one with the pad. Teeth only belong to top pad.
        coupled_pad = draw.rectangle(coupled_pad_width,
                                     coupled_pad_height + pad_height)
        coupler_pad_round = draw.Point(0., (coupled_pad_height + pad_height) /
                                       2).buffer(coupled_pad_width / 2,
                                                 resolution=16,
                                                 cap_style=CAP_STYLE.round)
        coupled_pad = draw.union(coupled_pad, coupler_pad_round)
        coupled_pad_left = draw.translate(
            coupled_pad, -(coupled_pad_width / 2. + coupled_pad_gap / 2.),
            +coupled_pad_height / 2. + pad_height + pad_gap / 2. -
            pad_height / 2)
        coupled_pad_right = draw.translate(
            coupled_pad, (coupled_pad_width / 2. + coupled_pad_gap / 2.),
            +coupled_pad_height / 2. + pad_height + pad_gap / 2. -
            pad_height / 2)
        pad_top_tmp1 = draw.union([circ_left_top, pad_top, circ_right_top])

        # Add trapezoid to the top pad
        trapezoid_top = self.create_trapezoid(0, p.taper_width_top, p.taper_width_base, p.taper_height, (pad_gap / 2 )- float(p.taper_height), p.fillet_radius,p.pad_width/2)
        # pad_top_tmp = draw.union(pad_top_tmp1, trapezoid_top)
        trapezoid_top_rotated = draw.rotate(trapezoid_top, 180, origin=(0,pad_gap/4+(pad_gap/2- float(p.taper_height))/2)) # Rotate trapezoid by 180 degrees

        # simple union can throw error if geometries do not overlap
        # So we check if the geometries overlap before union
        # pad_top_tmp = draw.union(pad_top_tmp1, trapezoid_top_rotated)
        pad_top_tmp = draw.union(
            pad_top_tmp1.buffer(0), 
            trapezoid_top_rotated.buffer(0)
        )



        # The coupler pads are only created if low_W=0 and low_H=+1
        for name in self.options.connection_pads:
            if self.options.connection_pads[name][
                    'loc_W'] == 0 and self.options.connection_pads[name][
                        'loc_H'] == +1:
                pad_top_tmp = draw.union([
                    circ_left_top, coupled_pad_left, pad_top_tmp, coupled_pad_right,
                    circ_right_top
                ])
        pad_top = pad_top_tmp

        #Curving the edges of the teeth where it joins te pad
        pad_top = pad_top.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(-p.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)
        pad_top = pad_top.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(p.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)

        # Round part for the bottom pad. And again you should unite all of them.
        pad_bot = draw.translate(pad, 0, -(pad_height + pad_gap) / 2.)
        circ_left_bot = draw.Point(-pad_width / 2, -(pad_height + pad_gap) /
                                   2.).buffer(pad_height / 2,
                                              resolution=16,
                                              cap_style=CAP_STYLE.round)
        circ_right_bot = draw.Point(pad_width / 2, -(pad_height + pad_gap) /
                                    2.).buffer(pad_height / 2,
                                               resolution=16,
                                               cap_style=CAP_STYLE.round)
        pad_bot = draw.union([pad_bot, circ_left_bot, circ_right_bot])

        # Add trapezoid to the bottom pad
        trapezoid_bot = self.create_trapezoid(0, p.taper_width_top, p.taper_width_base, p.taper_height, -(pad_gap) / 2, p.fillet_radius,p.pad_width/2)
        pad_bot = draw.union(pad_bot, trapezoid_bot)

###############################################################################################
        # Pins from thr center of the qubit pads 
        #Cordinates for the center of the pin
        taper_width_top = p.taper_width_top 

        top_pin = LineString([
            [0, pad_gap / 2],
            [0, pad_gap / 2 - p.taper_height/2 ]  
              # Extend outward
        ])

        bottom_pin = LineString([
            [0, -pad_gap / 2],
            [0, -pad_gap / 2 + p.taper_height/2 ] # Start from pad bottom edge
              # Extend outward
        ])

###############################################################################################
        ##Pins outside the qubit pocket 
        # Define the pin positions relative to (0,0)
        top_pocket_pin = LineString([
            [ p.pocket_width / 2 + p.chrgln_pin_y_offset, p.chrgln_pin_x_offset],
            [ p.pocket_width / 2 + p.chrgln_pin_y_offset,0]  # Start point (near pocket)
             # Extend outward
        ])

        bottom_pocket_pin = LineString([
            [-p.pocket_width / 2 - p.chrgln_pin_y_offset, p.chrgln_pin_x_offset],
            [-p.pocket_width / 2 - p.chrgln_pin_y_offset,0]  # Start point (near pocket)
              # Extend outward
        ])

###############################################################################################
        # rect_jj = draw.LineString([(0, -pad_gap / 2+p.taper_height), (0, +pad_gap / 2-p.taper_height)])
        # the draw.rectangle representing the josephson junction
        # rect_jj = draw.rectangle(p.inductor_width, p.inductor_height)
        rect_jj = draw.LineString([(0, -p.inductor_height / 2), (0, p.inductor_height / 2)])
        
        rect_pk = draw.rectangle(p.pocket_width, p.pocket_height)

        # to curve the edges of the qubit pocket
        rect_pk = rect_pk.buffer(p.fillet_radius_gap, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)
        

        # Rotate and translate all qgeometry as needed.
        # Done with utility functions in Metal 'draw_utility' for easy rotation/translation
        # NOTE: Should modify so rotate/translate accepts qgeometry, would allow for
        # smoother implementation.
        polys = [rect_jj, pad_top, pad_bot, rect_pk,top_pin, bottom_pin,top_pocket_pin, bottom_pocket_pin]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk,top_pin, bottom_pin,top_pocket_pin, bottom_pocket_pin] = polys

        # Use the geometry to create Metal qgeometry
        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        # self.add_qgeometry('poly', dict(
        #     rect_jj=rect_jj), helper=True)
        self.add_qgeometry('junction',
                           dict(rect_jj=rect_jj),
                           width=p.inductor_width)
        
        #Pins from the center of the qubit pads
        self.add_pin("pin_island", points=list(top_pin.coords), width=taper_width_top, input_as_norm=True)
        self.add_pin("pin_reservior", points=list(bottom_pin.coords), width=taper_width_top, input_as_norm=True)

         # Add pins to the component
        self.add_pin("bottom_pin", points=list(top_pocket_pin.coords), width=taper_width_top)
        self.add_pin("top_pin", points=list(bottom_pocket_pin.coords), width=taper_width_top, input_as_norm=True)

        ##############################################################################################################
        
        # # Get the x-coordinates of the qubit pads (same as existing qubit pads)
        # pad_x_position = p.pad_height / 2  # Center of the qubit pad

        # # Define the new pins' start and end points aligned with the qubit pads
        # extra_top_pin_start = np.array([pad_x_position, (p.pocket_width / 2) + p.pin_y_distance])
        # extra_top_pin_end = np.array([pad_x_position, (p.pocket_width  / 2) + p.pin_y_distance + pin_length])

        # extra_bottom_pin_start = np.array([pad_x_position, (-p.pocket_width  / 2) - p.pin_y_distance])
        # extra_bottom_pin_end = np.array([pad_x_position, (-p.pocket_width  / 2) - p.pin_y_distance - pin_length])

        # print("Pocket Height:", p.pocket_height)
        # print("Pin Y Distance:", p.pin_y_distance)
        # print("Expected Pin Locations:", extra_top_pin_start, extra_top_pin_end, extra_bottom_pin_start, extra_bottom_pin_end)

        # # Add the new pins outside the pocket, aligned with the qubit pads
        # self.add_pin("extra_pin_top", points=[extra_top_pin_start, extra_top_pin_end], width=taper_width_top, input_as_norm=True)
        # self.add_pin("extra_pin_bottom", points=[extra_bottom_pin_start, extra_bottom_pin_end], width=taper_width_top, input_as_norm=True)

        #########################################################################################################################
    def make_connection_pads(self):
        # """Makes standard transmon in a pocket."""
        for name in self.options.connection_pads:
            self.make_connection_pad(name)

    def make_connection_pad(self, name: str):
        # """Makes n individual connector.

        # Args:
        #     name (str) : Name of the connector
        # """

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.connection_pads[name]  # parser on connector options

        # define commonly used variables once
        cpw_width = pc.cpw_width
        cpw_extend = pc.cpw_extend
        pad_width = pc.pad_width
        pad_height = pc.pad_height
        pad_cpw_shift = pc.pad_cpw_shift
        pocket_rise = pc.pocket_rise
        pocket_extent = pc.pocket_extent

        loc_W = float(pc.loc_W)
        loc_W, loc_H = float(pc.loc_W), float(pc.loc_H)
        if float(loc_W) not in [-1., +1., 0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a transmon qubit with loc_W and'
                ' loc_H that are not +1, -1, or 0? Are you sure you want to do this?'
            )

        # Define the geometry
        # Connector pad

        if float(loc_W) != 0:
            connector_pad = draw.rectangle(pad_width, pad_height,
                                           -pad_width / 2, pad_height / 2)
            connector_pad = connector_pad.buffer(pc.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=pc.fillet_resolution) 
            # Connector CPW wire
            connector_wire_path = draw.wkt.loads(f"""LINESTRING (\
                0 {pad_cpw_shift+cpw_width/2}, \
                {pc.pad_cpw_extent}                           {pad_cpw_shift+cpw_width/2}, \
                {(p.pocket_width-p.pad_width)/2-pocket_extent} {pad_cpw_shift+cpw_width/2+pocket_rise}, \
                {(p.pocket_width-p.pad_width)/2+cpw_extend}    {pad_cpw_shift+cpw_width/2+pocket_rise}\
                                            )""")
        else:
            connector_pad = draw.rectangle(pad_width, pad_height, 0,
                                           pad_height / 2)
            connector_pad = connector_pad.buffer(pc.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=pc.fillet_resolution)
            
        ## Made chnages here to add buffer to curve edges 
            # connector_pad=connector_pad_1.buffer(p.edge_curve)
            connector_wire_path = draw.LineString(
                [[0, pad_height],
                 [
                     0,
                     (p.pocket_width / 2 - p.pad_height - p.pad_gap / 2 -
                      pc.pad_gap) + cpw_extend
                 ]])

        # Position the connector, rotate and translate
        objects = [connector_pad, connector_wire_path]

        if loc_W == 0:
            loc_Woff = 1
        else:
            loc_Woff = loc_W

        objects = draw.scale(objects, loc_Woff, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            loc_W * (p.pad_width) / 2.,
            loc_H * (p.pad_height + p.pad_gap / 2 + pc.pad_gap))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [connector_pad, connector_wire_path] = objects

        self.add_qgeometry('poly', {f'{name}_connector_pad': connector_pad})
        self.add_qgeometry('path', {f'{name}_wire': connector_wire_path},
                           width=cpw_width)
        self.add_qgeometry('path', {f'{name}_wire_sub': connector_wire_path},
                           width=cpw_width + 2 * pc.cpw_gap,
                           subtract=True)

        ############################################################

        # add pins
        points = np.array(connector_wire_path.coords)
        self.add_pin(name,
                     points=points[-2:],
                     width=cpw_width,
                     input_as_norm=True)


###################################################################################################################################


class TransmonTaperedInsets(BaseQubit):
    """Transmon pocket with tapered connection pads, with insets for improved coupling. 
    The qubit-resonator coupling has also been modified. 

    Inherits `BaseQubit` class

    Description:
        Create a standard pocket (floating) transmon qubit with tapered (optional) pads
        The qubit-resonator coupling has been modified to minimise participation ratios
        Insets (optional) on the pads for additional couplings

    Options:
        Convention: Values (unless noted) are strings with units included,
        (e.g., '30um')

    Pocket:
        * pad_gap            - the distance between the two charge islands, which is also the
          resulting 'length' of the pseudo junction
        * inductor_width     - width of the pseudo junction between the two charge islands (should be same as josephson junction width)
        * inductor_height    - height of the pseudo junction between the two charge islands 
        * pad_width          - the width (x-axis) of the charge island pads, except the circle radius from both sides
        * pad_height         - the size (y-axis) of the charge island pads
        * pocket_width       - size of the pocket (cut out in ground) along x-axis
        * pocket_height      - size of the pocket (cut out in ground) along y-axis
        * fillet_radius_gap  - The radius of the curved edges of the pocket
        * fillet_resolution  - number of points used to calculate the fillet/curve
        * chrgln_pin_x_offset =     - horizontal distance of the pin from the qubit center
        * chrgln_pin_y_offset =     - vertical distance of the pin from the pocket edge
        * coupled_pad_gap    - the distance between the two teeth shape
        * coupled_pad_width  - the width (x-axis) of the teeth shape on the island pads
        * coupled_pad_height - the size (y-axis) of the teeth shape on the island pads
                            
    Connector lines:
        * pad_gap        - space between the connector pad and the charge island it is
          nearest to
        * pad_width      - width (x-axis) of the connector pad
        * pad_height     - height (y-axis) of the connector pad
        * pad_cpw_shift  - shift the connector pad cpw line by this much away from qubit
        * pad_cpw_extent - how long should the pad be - edge that is parallel to pocket
        * cpw_width      - center trace width of the CPW line
        * cpw_gap        - dielectric gap width of the CPW line
        * cpw_extend     - depth the connector line extends into ground (past the pocket edge)
        * pocket_extent  - How deep into the pocket should we penetrate with the cpw connector
        * fillet_radius - The radius of the curved edges of the connector pad
          (into the ground plane)
        * fillet_radius_inner - The radius of curvature on inner corners
        * fillet_radius_outer - The radius of curvature on outer corners
        * pocket_rise    - How far up or down relative to the center of the transmon should we
          elevate the cpw connection point on the ground plane
        * loc_W / H      - which 'quadrant' of the pocket the connector is set to, +/- 1 (check
          if diagram is correct)


    Sketch:
        Below is a sketch of the qubit
        ::

                 +1              0             +1
                _________________________________
            -1  |                |    coupler    |  +1      Y
                |               |_|              |          ^
                |        __________________      |          |
                |       /     island       \     |          |----->  X
                |       >                  <     |
                |       \__________________/     |
                |             |       |          |
                |               |   |   taper    |
                |                |_|             |
                |                |               |
                |                |               |
                |  pocket        x               |                                              
                |                |               |  
                |                ||              |    
                |               |  |             |                      
                |        _____|______|_____      |
                |       /                  \     |
                |       > insets           <     |
                |       \__________________/     |
                |                                |
                |                                |
            -1  |________________________________|                          
                                    ^^^^ pocket lower tighten
    .. meta::
        Transmon Pocket Tapered With Insets

    """

    # Default drawing options
    default_options = Dict(
        # JJ box
        inductor_width="20um",
        inductor_height="22um",
        # Pins
        chrgln_pin_x_offset="30um",  # User-defined horizontal distance from the qubit center
        chrgln_pin_y_offset="50um",  # User-defined vertical distance from the pocket edge
        junction_pin_inset="0um",    # amount junction pins are inset from pad edge (where patches are placed)
        # Pads
        pad_gap="100um",
        pad_width="800um",
        pad_height="110um",
        # Pocket
        pocket_width="1000um",
        pocket_height="800um",
        pocket_lower_tighten="120um",
        # Tapered part of qubit
        taper_width_top="40um",
        taper_width_base="200um",
        taper_height="40um",
        taper_fillet_radius="3um",
        fillet_resolution_tapered=64,
        # Insets
        inset_width="0um",
        inset_depth="100um",
        inset_fillet_radius="10um",
        # Fillet settings
        fillet_radius="50um",
        fillet_resolution=16,
        fillet_radius_gap="50um",
        fillet_resolution_gap=16,
        # Coupler - resonator
        coupled_pad_height="0um",
        coupled_pad_width="0um",
        coupled_pad_gap="0um",
        # orientation = 90 has dipole aligned along the +X axis, while 0 aligns to the +Y axis
        orientation=0,
        _default_connection_pads=Dict(
            pad_gap="100um",
            pad_height="30um",
            pad_width="200um",
            pad_cpw_shift="0um",
            pad_cpw_extent="25um",
            cpw_width="10um",
            cpw_gap="6um",
            cpw_extend="200um",
            pocket_extent="5um",
            pocket_rise="0um",
            fillet_radius_inner="15um",
            fillet_radius_outer="75um",
            fillet_resolution=16,
            pin_y_distance="0um",
            loc_W="0",  # width location  only +-1 or 0,
            loc_H="1",  # height location  only +-1 or 0
        ),
    )
    # """Default drawing options"""
    component_metadata = Dict(
        short_name="Pocket",
        _qgeometry_table_path="True",
        _qgeometry_table_poly="True",
        _qgeometry_table_junction="True",
    )
    # """Component metadata"""

    # TOOLTIP = """Transmon pocket with tapering."""
    def __init__(self, design, name: str = None, options: Dict = None, **kwargs):
        super().__init__(design, name, options, **kwargs)

    def make(self):
        self.make_pocket()
        self.make_connection_pads()

    # Function to create trapezoid
    def create_trapezoid(
        self, center_x, top_width, base_width, height, y_offset, rfillet, startx
    ):
        half_top = float(top_width) / 2
        half_base = float(base_width) / 2
        coordinates = [
            (center_x - half_base, y_offset),
            (center_x - half_top, y_offset + height),
            (center_x + half_top, y_offset + height),
            (center_x + half_base, y_offset),
        ]

        # Create a Polygon object from the coordinates
        trapezoid_polygon = Polygon(coordinates)

        # Convert the coordinates to numpy array
        raw_data = np.array(trapezoid_polygon.exterior.coords)

        start_x = center_x - startx
        start_y = y_offset

        end_x = center_x + startx
        end_y = y_offset

        start = [[start_x, start_y]]
        end = [[end_x, end_y]]

        # Concatenate arrays after converting start and end points to 2D arrays
        raw_data_2 = np.vstack((start, raw_data, end))

        raw_data_1 = np.delete(raw_data_2, 5, 0)

        # Curve the edges of the polygon
        data = QUtilities.calc_filleted_path(raw_data_1, rfillet, 9)

        # Create a new Polygon object from the curved coordinates
        trapezoid_with_curved_edges = Polygon(data)

        return trapezoid_with_curved_edges

    # """Makes standard transmon in a pocket."""
    def make_pocket(self):

        # self.p allows us to directly access parsed values (string -> numbers) from the user option
        p = self.p
        #  pcop = self.p.coupled_pads[name]  # parser on connector options

        # since we will reuse these options, parse them once and define them as variables
        pad_width = p.pad_width
        pad_height = p.pad_height
        pad_gap = p.pad_gap
        coupled_pad_height = p.coupled_pad_height
        coupled_pad_width = p.coupled_pad_width
        coupled_pad_gap = p.coupled_pad_gap

        # TOP PAD
        # make the pads as rectangles (shapely polygons)
        pad = draw.rectangle(pad_width, pad_height)
        pad_top = draw.translate(pad, 0, +(pad_height + pad_gap) / 2.0)
        # Here, you make your pads round. Not sharp shape on the left and right sides and also this should be the same for the bottom pad as the top pad.
        circ_left_top = draw.Point(
            -pad_width / 2.0, +(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        circ_right_top = draw.Point(
            pad_width / 2.0, +(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        pad_top_tmp1 = draw.union([circ_left_top, pad_top, circ_right_top])
        # TAPER - Add trapezoid to the top pad
        if p.taper_height > 0:
            trapezoid_top = self.create_trapezoid(
                0,
                p.taper_width_top,
                p.taper_width_base,
                p.taper_height,
                (pad_gap / 2) - float(p.taper_height),
                p.taper_fillet_radius,
                p.pad_width / 2,
            )
            # pad_top_tmp = draw.union(pad_top_tmp1, trapezoid_top)
            trapezoid_top_rotated = draw.rotate(
                trapezoid_top,
                180,
                origin=(0, pad_gap / 4 + (pad_gap / 2 - float(p.taper_height)) / 2),
            )  # Rotate trapezoid by 180 degrees
            # create union
            pad_top_tmp = draw.union(
                pad_top_tmp1.buffer(0), trapezoid_top_rotated.buffer(0)
            )
        else:
            pad_top_tmp = pad_top_tmp1

        # COUPLER REGION - TOP PAD
        # In here you create the teeth part and then you union them as one with the pad. Teeth only belong to top pad.
        coupled_pad = draw.rectangle(coupled_pad_width, coupled_pad_height + pad_height)
        coupler_pad_round = draw.Point(
            0.0, (coupled_pad_height + pad_height) / 2
        ).buffer(coupled_pad_width / 2, resolution=16, cap_style=CAP_STYLE.round)
        coupled_pad = draw.union(coupled_pad, coupler_pad_round)
        coupled_pad_left = draw.translate(
            coupled_pad,
            -(coupled_pad_width / 2.0 + coupled_pad_gap / 2.0),
            +coupled_pad_height / 2.0 + pad_height + pad_gap / 2.0 - pad_height / 2,
        )
        coupled_pad_right = draw.translate(
            coupled_pad,
            (coupled_pad_width / 2.0 + coupled_pad_gap / 2.0),
            +coupled_pad_height / 2.0 + pad_height + pad_gap / 2.0 - pad_height / 2,
        )
        # The coupler pads are only created if low_W=0 and low_H=+1
        for name in self.options.connection_pads:
            if (
                self.options.connection_pads[name]["loc_W"] == 0
                and self.options.connection_pads[name]["loc_H"] == +1
            ):
                pad_top_tmp = draw.union(
                    [
                        circ_left_top,
                        coupled_pad_left,
                        pad_top_tmp,
                        coupled_pad_right,
                        circ_right_top,
                    ]
                )
        pad_top = pad_top_tmp
        # Curving the edges of the teeth where it joins the pad
        # outer corners
        pad_top = pad_top.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(
            -p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # inner corners
        pad_top = pad_top.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(
            p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )

        # BOTTOM PAD
        # Round part for the bottom pad. And again you should unite all of them.
        pad_bot = draw.translate(pad, 0, -(pad_height + pad_gap) / 2.0)
        circ_left_bot = draw.Point(
            -pad_width / 2, -(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        circ_right_bot = draw.Point(
            pad_width / 2, -(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        pad_bot = draw.union([pad_bot, circ_left_bot, circ_right_bot])
        # TAPER - Add trapezoid to the bottom pad
        if p.taper_height > 0:
            trapezoid_bot = self.create_trapezoid(
                center_x=0,
                top_width=p.taper_width_top,
                base_width=p.taper_width_base,
                height=p.taper_height,
                y_offset=-(pad_gap / 2),
                rfillet=p.taper_fillet_radius,
                startx=p.pad_width / 2,
            )
            pad_bot = draw.union(pad_bot, trapezoid_bot)

        # outer corners
        pad_bot = pad_bot.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(
            -p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # inner corners
        pad_bot = pad_bot.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(
            p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # pad_bot = pad_bot.buffer(p.fillet_radius, join_style=1).buffer(-p.fillet_radius, join_style=1)

        ###############################################################################################
        # INSETS - cut into the pads for additional coupling
        # define trapezoid
        if p.inset_depth != 0 and p.inset_width != 0:
            inset = draw.rectangle(
                p.inset_depth - (2 * p.inset_fillet_radius),
                p.inset_width - (2 * p.inset_fillet_radius),
            )
            inset_rounded = inset.buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=9,
            )
            inset = inset_rounded
            # define x and y locations for insets
            inset_y_top = (pad_gap / 2) + (p.pad_height / 2)
            inset_y_bot = -inset_y_top
            inset_x_left = (
                -(pad_width / 2)
                - (pad_height / 2)
                + (p.inset_depth / 2)
                - p.inset_fillet_radius
            )
            inset_x_right = -inset_x_left
            # position trapezoids
            inset_left_top = draw.translate(inset, inset_x_left, inset_y_top)
            inset_right_top = draw.translate(inset, inset_x_right, inset_y_top)
            inset_left_bot = draw.translate(inset, inset_x_left, inset_y_bot)
            inset_right_bot = draw.translate(inset, inset_x_right, inset_y_bot)
            # Subtract the insets from the pads
            pad_top_tmp1 = draw.subtract(pad_top, inset_left_top)
            pad_top_tmp2 = draw.subtract(pad_top_tmp1, inset_right_top)
            pad_top = pad_top_tmp2
            pad_bot_tmp1 = draw.subtract(pad_bot, inset_left_bot)
            pad_bot_tmp2 = draw.subtract(pad_bot_tmp1, inset_right_bot)
            pad_bot = pad_bot_tmp2
            # re-round the new pads (inner corners only)
            pad_bot = pad_bot.buffer(
                -p.inset_fillet_radius, join_style=2, cap_style=3
            ).buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=p.fillet_resolution * 3,
            )
            pad_top = pad_top.buffer(
                -p.inset_fillet_radius, join_style=2, cap_style=3
            ).buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=p.fillet_resolution * 3,
            )
            pad_bot = pad_bot.buffer(0)
            pad_top = pad_top.buffer(0)

        ###############################################################################################
        # Pins from the center of the qubit pads
        # Coordinates for the center of the pin
        taper_width_top = p.taper_width_top
        if p.taper_height == 0:
            pin_route_l = 1e-6
        else:
            pin_route_l = p.taper_height
        top_pin = LineString(
            [
                [0, pad_gap / 2 + p.junction_pin_inset],
                [0, pad_gap / 2 + p.junction_pin_inset - pin_route_l / 2],
                # Extend outward
            ]
        )
        bottom_pin = LineString(
            [
                [0, -pad_gap / 2 - p.junction_pin_inset],
                [0, -pad_gap / 2 - p.junction_pin_inset + pin_route_l / 2],  # Start from pad bottom edge
                # Extend outward
            ]
        )

        ###############################################################################################
        # Pins outside the qubit pocket
        # Define the pin positions relative to (0,0)
        # 4 pins, 2 on each side, aligned with the qubit pads and taper, and starting from the edge of the pocket.
        top_right_pocket_pin = LineString(
            [
                [p.pocket_width / 2 + p.chrgln_pin_y_offset, p.pad_gap/2 + p.pad_height/2 + pin_route_l /2],  # Start point (near pocket)
                [
                    p.pocket_width / 2 + p.chrgln_pin_y_offset,
                    p.pad_gap/2 + p.pad_height/2 - pin_route_l / 2,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        bottom_right_pocket_pin = LineString(
            [
                [p.pocket_width / 2 + p.chrgln_pin_y_offset, -p.pad_gap/2 - p.pad_height/2 - pin_route_l],  # Start point (near pocket)
                [
                    p.pocket_width / 2 + p.chrgln_pin_y_offset,
                    -p.pad_gap/2 - p.pad_height/2  ,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        top_left_pocket_pin = LineString(
            [
                [-p.pocket_width / 2 - p.chrgln_pin_y_offset, p.pad_gap/2 + p.pad_height/2 + pin_route_l /2],  # Start point (near pocket)
                [
                    -p.pocket_width / 2 - p.chrgln_pin_y_offset,
                    p.pad_gap/2 + p.pad_height/2 - pin_route_l / 2,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        bottom_left_pocket_pin = LineString(
            [
                [-p.pocket_width / 2 - p.chrgln_pin_y_offset, -p.pad_gap/2 - p.pad_height/2 - pin_route_l ],  # Start point (near pocket)
                [
                    -p.pocket_width / 2 - p.chrgln_pin_y_offset,
                    -p.pad_gap/2 - p.pad_height/2,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        ###############################################################################################
        # Josephson Junction
        rect_jj = draw.LineString(
            [(0, -p.inductor_height / 2), (0, p.inductor_height / 2)]
        )

        # Pocket
        rect_pk = draw.rectangle(
            p.pocket_width - 2 * p.fillet_radius_gap,
            p.pocket_height - p.pocket_lower_tighten - 2 * p.fillet_radius_gap,
            yoff=p.pocket_lower_tighten / 2,
        )

        # to curve the edges of the qubit pocket
        rect_pk = rect_pk.buffer(
            p.fillet_radius_gap,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )

        # Rotate and translate all qgeometry as needed.
        polys = [
            rect_jj,
            pad_top,
            pad_bot,
            rect_pk,
            top_pin,
            bottom_pin,
            top_right_pocket_pin,
            bottom_right_pocket_pin,
            top_left_pocket_pin,
            bottom_left_pocket_pin,
        ]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [
            rect_jj,
            pad_top,
            pad_bot,
            rect_pk,
            top_pin,
            bottom_pin,
            top_right_pocket_pin,
            bottom_right_pocket_pin,
            top_left_pocket_pin,
            bottom_left_pocket_pin,
        ] = polys

        # Add shapes as qiskit geometries
        self.add_qgeometry("poly", dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry("poly", dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry("junction", dict(rect_jj=rect_jj), width=p.inductor_width)

        # Pins from the center of the qubit pads
        self.add_pin(
            "pin_island",
            points=list(top_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )
        self.add_pin(
            "pin_reservior",
            points=list(bottom_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )

        # Add pins to the component
        self.add_pin(
            "top_right_pin", points=list(top_right_pocket_pin.coords), width=taper_width_top
        )
        self.add_pin(
            "bottom_right_pin",
            points=list(bottom_right_pocket_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )
        self.add_pin(
            "top_left_pin", points=list(top_left_pocket_pin.coords), width=taper_width_top
        )
        self.add_pin(
            "bottom_left_pin",
            points=list(bottom_left_pocket_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )

    def make_connection_pads(self):
        # """Makes standard transmon in a pocket."""
        for name in self.options.connection_pads:
            self.make_connection_pad(name)

    def make_connection_pad(self, name: str):
        # """Makes n individual connector.

        # Args:
        #     name (str) : Name of the connector
        # """

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.connection_pads[name]  # parser on connector options

        # define commonly used variables once
        r = pc.fillet_radius_inner
        r_outer = pc.fillet_radius_outer
        cpw_width = pc.cpw_width
        cpw_extend = pc.cpw_extend
        pad_width_fillet = pc.pad_width - 2 * r
        pad_height_fillet = pc.pad_height - 2 * r
        pad_width = pc.pad_width
        pad_height = pc.pad_height
        pad_cpw_shift = pc.pad_cpw_shift
        pocket_rise = pc.pocket_rise
        pocket_extent = pc.pocket_extent

        assert (
            pad_width_fillet >= 0
        ), f"Error: pad_width {pad_width} is too small for the fillet radius {r}. Either increase the width or decrease the fillet radius of the connection pads."
        assert (
            pad_height_fillet >= 0
        ), f"Error: pad_height {pad_height} is too small for the fillet radius {r}. Either increase the height or decrease the fillet radius of the connection pads."
        assert (
            pad_width / 2 - cpw_width / 2 - r >= r_outer
        ), f"Error: fillet_radius_outer {r_outer} is too large for the pad_width {pad_width} and cpw_width {cpw_width}. Either decrease the fillet_radius_outer or increase the pad_width."

        loc_W = float(pc.loc_W)
        loc_W, loc_H = float(pc.loc_W), float(pc.loc_H)
        if float(loc_W) not in [-1.0, +1.0, 0] or float(loc_H) not in [-1.0, +1.0]:
            self.logger.info(
                "Warning: Did you mean to define a transmon qubit with loc_W and"
                " loc_H that are not +1, -1, or 0? Are you sure you want to do this?"
            )

        # Define the geometry
        # Connector pad
        if float(loc_W) != 0:
            connector_pad = draw.rectangle(
                pad_width_fillet, pad_height_fillet, -pad_width / 2, pad_height / 2
            )
            connector_pad = connector_pad.buffer(
                r,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=pc.fillet_resolution,
            )
            # add corners
            if r_outer > 0:
                print("Skipping: connector pad rounding not yet implemented for W != 0")

            # Connector CPW wire
            connector_wire_path = draw.wkt.loads(
                f"""LINESTRING (\
                0 {pad_cpw_shift+cpw_width/2}, \
                {pc.pad_cpw_extent}                           {pad_cpw_shift+cpw_width/2}, \
                {(p.pocket_width-p.pad_width)/2-pocket_extent} {pad_cpw_shift+cpw_width/2+pocket_rise}, \
                {(p.pocket_width-p.pad_width)/2+cpw_extend}    {pad_cpw_shift+cpw_width/2+pocket_rise}\
                                            )"""
            )
        else:
            connector_pad = draw.rectangle(
                pad_width_fillet, pad_height_fillet, 0, pad_height / 2
            )
            connector_pad = connector_pad.buffer(
                r,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=pc.fillet_resolution,
            )
            # add corners
            if r_outer > 0:
                # --- Add the small corner box first ---
                connector_pad_corners_sq = draw.rectangle(
                    r_outer,
                    r_outer,
                    -(cpw_width / 2) - (r_outer / 2),
                    pad_height + (r_outer / 2),
                )
                # --- Create a quarter circle to cut the outer corner ---
                corner_circle = draw.Point(
                    -(cpw_width / 2) - (r_outer), pad_height + r_outer
                ).buffer(r_outer, resolution=32)
                corner_subtract_l = draw.subtract(
                    connector_pad_corners_sq, corner_circle
                )
                # repeat on right side
                connector_pad_corners_sq = draw.rectangle(
                    r_outer,
                    r_outer,
                    (cpw_width / 2) + (r_outer / 2),
                    pad_height + (r_outer / 2),
                )
                corner_circle = draw.Point(
                    (cpw_width / 2) + (r_outer), pad_height + r_outer
                ).buffer(r_outer, resolution=32)
                corner_subtract_r = draw.subtract(
                    connector_pad_corners_sq, corner_circle
                )
                # Union the corners
                connector_pad = draw.union(
                    [connector_pad, corner_subtract_l, corner_subtract_r]
                )
            # CPW path
            connector_wire_path = draw.LineString(
                [
                    [0, pad_height],
                    [
                        0,
                        (p.pocket_width / 2 - p.pad_height - p.pad_gap / 2 - pc.pad_gap)
                        + cpw_extend,
                    ],
                ]
            )

        # Position the connector, rotate and translate
        objects = [connector_pad, connector_wire_path]

        if loc_W == 0:
            loc_Woff = 1
        else:
            loc_Woff = loc_W

        objects = draw.scale(objects, loc_Woff, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            loc_W * (p.pad_width) / 2.0,
            loc_H * (p.pad_height + p.pad_gap / 2 + pc.pad_gap),
        )
        objects = draw.rotate_position(objects, p.orientation, [p.pos_x, p.pos_y])
        [connector_pad, connector_wire_path] = objects

        self.add_qgeometry("poly", {f"{name}_connector_pad": connector_pad})
        self.add_qgeometry(
            "path", {f"{name}_wire": connector_wire_path}, width=cpw_width
        )
        self.add_qgeometry(
            "path",
            {f"{name}_wire_sub": connector_wire_path},
            width=cpw_width + 2 * pc.cpw_gap,
            subtract=True,
        )

        ############################################################

        # add pins
        points = np.array(connector_wire_path.coords)
        self.add_pin(name, points=points[-2:], width=cpw_width, input_as_norm=True)

    def get_resonator_length_mm(self, connection_pad="readout"):
        """
        Gives total length of resonator segment.
        """
        return QUtilities.calc_points_on_path([0], self.design, component_name=self.name, trace_name=f"{connection_pad}_wire")[-1]


class _FluxoniumPocket(BaseQubit):
    # -*- coding: utf-8 -*-

    # This code is part of Qiskit.
    #
    # (C) Copyright IBM 2017, 2021.
    #
    # This code is licensed under the Apache License, Version 2.0. You may
    # obtain a copy of this license in the LICENSE.txt file in the root directory
    # of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
    #
    # Any modifications or derivative works of this code must retain this
    # copyright notice, and modified files need to carry a notice indicating
    # that they have been altered from the originals.

    # This class was created by Figen YILMAZ, Christian Kraglund Andersen
    """The base `_FluxoniumPocket` class.

    Inherits `BaseQubit` class.

    Description:
        Create a standard pocket fluxonium qubit for a ground plane,
        with two pads connected by a junction and a kinetic inductor 
        between two pads (see drawing below), creates a magnetic loop
        in between the main JJ and the JJ array (can be a nanowire too).

    Connector lines were added using the `flux-bias line` & `charge line`
    & `readout line` dictionaries, the 'fake_flux_bias_line' is in use 
    for Maxwell capacitance matrix, only. Each connector pad has a name and 
    a list of default properties (except for the fake_flux_bias_line).

    Sketch:
        Below is a sketch of the qubit
        ::
                               0
        
                             | | |  charge_line
                 +1          \___/           +1
                _______________________________
            -1 |              ___              | +1        Y
               |             /   \             |           ^   
               |             \   /             |           |
               |              | |___           |           |----->  X
               |              |_|   |    ______|
               |               |    |   |  ____|-- fake_flux_bias_line
               |               x    |   | |____|__
               |               |    |   |_________-- flux_bias_line
               |              | |___|          |
               |              | |              |
               |             /   \             |
               |             \___/             |
            -1 |_______________________________|  +1
                 -1        ___________      -1
                          /  _______  \
                         /  /       \  \            
                         \  \_______/  /
                          \___  |  ___/
                              | | | 
                          readout_line

    .. image::
        _FluxoniumPocket.png
            
    .. meta::
        Fluxonium Pocket

    Default Options:
        * pad_gap: '30um' -- The distance between the two charge islands, which is also the resulting 'length' of the pseudo junction
        * jj_width: '20um' -- Width of the pseudo junction between the two charge islands (if in doubt, make the same as pad_gap). Really just for simulating in HFSS / other EM software
        * jj_orientation: '1' -- Direction of the main JJ, fab related but all the JJs has to look at the same direction
        * pad_width: '15um' -- The width (x-axis) of the charge island pads
        * pad_height: '100um' -- The size (y-axis) of the charge island pads
        * pad_radius: '60um' -- Radius of the circles at the end of the pads
        * l_width: '1um' -- Width of the kinetic inductor along the x-axis
        * array_length: '130um' -- Length of the array along the y-axis
        * l_arm_width: '2um' -- Width of the arm of the kinetic inductor along the y-axis
        * l_arm_length: '20um' -- Length of the arm of the kinetic inductor along the x-axis
        * l_inductance: '200nH' -- The kinetic energy for the inductor
        * L_j: '34.38nH' -- Josephson junction inductor value
        * C_j: '1.4fF' -- Josephson junction capacitance value
        * inductor_orientation: '-1' -- Direction of the array, fab related but all the JJs has to look at the same direction
        * pocket_width: '800um' -- Size of the pocket (cut out in ground) along x-axis
        * pocket_height: '800um' -- Size of the pocket (cut out in ground) along y-axis
        * nanowire_inductor: 'True' -- Boolean for nanowire instead of array
        * gds_cell_inductor: 'gds_cell_inductor' -- GDS cell name
        * orientation: '0' -- Degree of qubit rotation
        * flux_bias_line_options=Dict
            * make_fbl = True -- Boolean to make the flux bias line 
            * fbl_sep: '100um' -- The separation between the flux bias line and the inductor along the x-axis
            * fbl_height: '50um' -- The height of the flux bias line along the y-axis
            * cpw_width: 'cpw_width' -- The width of the flux bias line
            * cpw_gap: 'cpw_gap' -- The dielectric gap width of the flux bias line
        * charge_line_options=Dict
            * make_cl = True -- Boolean to make the charge line
            * cl_length: '80um' -- Length of the charge line along the y-axis
            * cl_sep: '15um' -- The separation between the connection pad and the pocket the y-axis
            * cpw_width: 'cpw_width' -- The width of the charge line
            * cpw_gap: 'cpw_gap' -- The dielectric gap width of the charge line
            * loc_W / H -- which 'quadrant' of the pocket the connector is set to, +/- 1 (check
              if diagram is correct)
        * readout_line_options=Dict(
            * make_rol = True -- Boolean to make the readout line
            * pad_sep: '50um' -- The separation between the connection pad and the capacitor pad the y-axis
            * pad_width: '150um' -- Width of the connection pad along the x-axis  
            * pad_height:'50um', -- Height of the connection pad along the y-axis
            * cpw_width: 'cpw_width', -- The width of the charge line
            * cpw_gap: 'cpw_gap' -- The dielectric gap width of the readout line
            * loc_W / H -- which 'quadrant' of the pocket the connector is set to, +/- 1 (check
              if diagram is correct)
    """

    component_metadata = Dict(short_name='_FluxoniumPocket',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True',
                              )
    """Component metadata"""

    # Default drawing options
    default_options = Dict(
        pad_gap='30um',
        jj_width='8um',
        jj_orientation='1',
        pad_width='15um',
        pad_height='100um',
        pad_radius='60um',
        l_width='1um',
        array_length='130um',
        l_arm_width='2um',
        l_arm_length='20um',
        l_inductance='218.823nH',
        L_j='40.7885nH',
        C_j='0.72fF',
        inductor_orientation='-1',
        pocket_width='900um',
        pocket_height='550um',
        nanowire_inductor='True',
        gds_cell_inductor='gds_cell_inductor',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl=False,
            fbl_sep='85um',
            fbl_height='50um',
            cpw_width='10um',
            cpw_gap='11.233um',
        ),
        charge_line_options=Dict(
            make_cl=False,
            cl_length='100um',
            cl_sep='-15um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='+1',  # height location only +1 or -1
        ),
        readout_line_options=Dict(
            make_rol=False,
            pad_sep='85um',
            pad_width='400um',
            pad_height='120um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='-1',  # height location only -1 or +1
        ))
    """Default drawing options"""

    TOOLTIP = """The base `_FluxoniumPocket` class."""

    def make(self):
        """Define the way the options are turned into QGeometry.

        The make function implements the logic that creates the geometry
        (poly, path, etc.) from the qcomponent.options dictionary of
        parameters, and the adds them to the design, using
        qcomponent.add_qgeometry(...), adding in extra needed
        information, such as layer, subtract, etc.
        """
        self.make_pocket()

        if self.p.flux_bias_line_options.make_fbl == True:
            self.make_flux_bias_line()
        if self.p.charge_line_options.make_cl == True:
            self.make_charge_line()
        if self.p.readout_line_options.make_rol == True:
            self.make_readout_line()

    def make_pocket(self):
        """Makes standard fluxonium in a pocket."""

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        # since we will reuse these options, parse them once and define them as variables
        pocket_width = p.pocket_width
        pocket_height = p.pocket_height
        pad_height = p.pad_height
        pad_width = p.pad_width
        pad_radius = p.pad_radius
        pad_gap = p.pad_gap
        l_arm_length = p.l_arm_length
        l_arm_width = p.l_arm_width
        l_width = p.l_width
        nanowire_inductor = p.nanowire

        # Drawing the kinectic inductor
        if nanowire_inductor == True:
            # MY COMMENT: I don't know how l_length is defined here since it is only defined
            inductor = draw.LineString(
                [(l_arm_length, l_length/2), (l_arm_length, -l_length/2)])
            # in the else below. regardless, this just defines the inductor as a straight simple line.
        else:
            l_length = p.array_length
            # This one is for JJ array
            # MY COMMENT: This just makes the inductor longer/shorter, inductor_orientation is just a scaling factor
            io = float(p.inductor_orientation)
            inductor = draw.LineString(
                [(l_arm_length-l_arm_width, io*l_length/2), (l_arm_length-l_arm_width, -io*l_length/2)])

        # Draw 'the arms' and make them curvy, first top arm and then same goes for the bottom
        l_arm_up = draw.Polygon([
            (pad_width/2, l_length/2+l_arm_width),  # point a
            (l_arm_length, l_length/2+l_arm_width),  # point b
            (l_arm_length, l_length/2),  # point c
            (pad_width/2, l_length/2),  # point d
        ])
        # Having semicircle with subtracting the geometries
        l_arm_up_fillet = draw.Point(
            l_arm_length, l_length/2).buffer(l_arm_width)
        cut_ply_up = draw.Polygon([
            (l_arm_length-l_arm_width*2, l_length/2+l_arm_width*2),   # point o
            (l_arm_length, l_length/2+l_arm_width*2),    # point p
            (l_arm_length, l_length/2),   # point r same with point c
            (l_arm_length+l_arm_width*2, l_length/2),   # point s
            (l_arm_length+l_arm_width*2, l_length/2-l_arm_width*2),  # point t
            (l_arm_length-l_arm_width*2, l_length/2-l_arm_width*2),  # point u
        ])
        # Having semicircle with subtracting the geometries
        l_arm_up_fillet = draw.subtract(l_arm_up_fillet, cut_ply_up)
        # Here the bottom arm
        l_arm_bot = draw.Polygon([
            (pad_width/2, -l_length/2),  # point e
            (l_arm_length, -l_length/2),  # point f
            (l_arm_length, -(l_length/2+l_arm_width)),  # point g
            (pad_width/2, -(l_length/2+l_arm_width)),  # point h
        ])
        l_arm_bot_fillet = draw.Point(
            l_arm_length, -l_length/2).buffer(l_arm_width)
        cut_ply_bot = draw.Polygon([
            ((l_arm_length-l_arm_width*2), -(l_length/2+l_arm_width*2)),   # point o
            (l_arm_length, -(l_length/2+l_arm_width*2)),    # point p
            (l_arm_length, -l_length/2),   # point r same with point c
            (l_arm_length+l_arm_width*2, -l_length/2),   # point s
            (l_arm_length+l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point t
            (l_arm_length-l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point u
        ])
        # Having semicircle with subtracting the geometries
        l_arm_bot_fillet = draw.subtract(l_arm_bot_fillet, cut_ply_bot)

        # Draw the pads (shapely polygons)
        pad_rect_top = draw.rectangle(
            pad_width, pad_height, 0, (pad_gap+pad_height)/2)
        pad_circle_top = draw.Point(
            0, (pad_radius+pad_height)).buffer(pad_radius)
        pad_top = draw.union(pad_rect_top, pad_circle_top,
                             l_arm_up, l_arm_up_fillet)
        pad_rect_bot = draw.rectangle(
            pad_width, pad_height, 0, -(pad_gap+pad_height)/2)
        pad_circle_bot = draw.Point(
            0, -(pad_radius+pad_height)).buffer(pad_radius)
        pad_bot = draw.union(pad_rect_bot, pad_circle_bot,
                             l_arm_bot, l_arm_bot_fillet)

        # Draw the junction
        # one can change the JJ orientation. Fab related detail.
        jj_o = float(p.jj_orientation)
        rect_jj = draw.LineString([(0, -pad_gap/2*jj_o), (0, +pad_gap/2*jj_o)])
        # the draw.rectangle representing the josephson junction
        # rect_jj = draw.rectangle(p.jj_width, pad_gap)

        # Draw the pocket
        rect_pk = draw.rectangle(pocket_width, pocket_height)

        # Rotate and translate all qgeometry as needed.
        polys = [rect_jj, pad_top, pad_bot, rect_pk, inductor]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk, inductor] = polys

        # Use the geometry to create Metal qgeometry
        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        # self.add_qgeometry('poly', dict(
        #     rect_jj=rect_jj), helper=True)
        self.add_qgeometry('junction',  # kinetic inductor
                           dict(inductor=inductor),
                           width=l_width,
                           hfss_inductance=p.l_inductance,
                           gds_cell_name=p.gds_cell_inductor)
        self.add_qgeometry('junction',  # the main JJ
                           dict(rect_jj=rect_jj),
                           width=p.jj_width,
                           hfss_inductance=p.L_j,
                           hfss_capacitance=p.C_j)

    def make_flux_bias_line(self):
        """ Adds flux bias line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option

        p = self.p
        pfb = self.p.flux_bias_line_options  # parser on connector options

        # define commonly used variables once
        fbl_sep = pfb.fbl_sep
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        # Define the geometry
        # Flux Bias Line
        # The position of flux bias line on the x-axis, starting point inside the pocket
        d = p.pocket_width/2
        # Draw the top line of the flux-bias line
        flux_bias_lineup = draw.Polygon([
            (d, fbl_height/2),   # point a
            (d, fbl_height/2+cpw_width),    # point b
            (fbl_sep, fbl_height/2+cpw_width),   # point c
            (fbl_sep, fbl_height/2),   # point d
        ])
        # Draw the middle line of the flux-bias line
        flux_bias_linemid = draw.Polygon([
            (fbl_sep-cpw_width, fbl_height/2),   # point e
            (fbl_sep, fbl_height/2),    # point f
            (fbl_sep, -fbl_height/2),   # point g
            (fbl_sep-cpw_width, -fbl_height/2),   # point h
        ])

        # Here we make flux-bias line curvy for the top side and also bottom side and the union all of them
        circle_top = draw.Point(fbl_sep, fbl_height/2).buffer(cpw_width)
        cut_ply = draw.Polygon([
            (fbl_sep*2, fbl_height+cpw_width),   # point o
            (fbl_sep, fbl_height+cpw_width),    # point p
            (fbl_sep, fbl_height/2),   # point r same with point d
            (fbl_sep/2, fbl_height/2),   # point s
            (fbl_sep/2, -fbl_height+cpw_width),  # point t
            (fbl_sep*2, -fbl_height+cpw_width),  # point u
        ])
        circle_top = draw.subtract(circle_top, cut_ply)
        # same goes for bottom edge
        circle_bot = draw.Point(fbl_sep, -fbl_height/2).buffer(cpw_width)
        cut_ply2 = draw.Polygon([
            (fbl_sep-cpw_width, -fbl_height/2),   # point v same with h or i
            (fbl_sep-cpw_width, fbl_height*2),    # point y
            (fbl_sep, fbl_height),   # point z
            (fbl_sep, fbl_height),   # point w
            (fbl_sep*2, fbl_height),  # point x
            (fbl_sep*2, -fbl_height/2),  # point k
        ])
        circle_bot = draw.subtract(circle_bot, cut_ply2)
        flux_bias_linebot = draw.Polygon([
            (d+fbl_sep/2, -fbl_height/2),   # point i
            (d+fbl_sep/2, -fbl_height/2-cpw_width),    # point k
            (fbl_sep, -fbl_height/2-cpw_width),   # point l
            (fbl_sep, -fbl_height/2),   # point m
        ])
        flux_bias_line = draw.union(
            flux_bias_lineup, flux_bias_linemid,  flux_bias_linebot, circle_top, circle_bot)

        # Flux Bias line's gap part, inside the GND
        flux_bias_line_gap = draw.rectangle(
            fbl_sep/2, cpw_width+cpw_gap*2, d+fbl_sep/4, -fbl_height/2-cpw_width/2)

        # Flux-Bias Line CPW wire
        port_line = draw.LineString([((d+fbl_sep/2), 0),
                                    ((d+fbl_sep/2), -(fbl_height+cpw_width))])

        # This port line is a fake port line, it is only in use during LOM analyses because we need to have an ungrounded line for the flux-bias
        fake_port_line = draw.LineString([(d, (fbl_height*2+cpw_width*2)),
                                          (d, -(fbl_height+cpw_width))])

        objects = [flux_bias_line, flux_bias_line_gap,
                   port_line, fake_port_line]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line, fake_port_line] = objects

        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})
        self.add_qgeometry(
            'poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)

        ####################################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line',
                     port_line_cords, cpw_width)

        fake_port_line_cords = list(
            draw.shapely.geometry.shape(fake_port_line).coords)
        self.add_pin('fake_flux_bias_line',
                     fake_port_line_cords, cpw_width)

    def make_charge_line(self):
        """ Adds charge line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.charge_line_options  # parser on charge line options

        # define commonly used variables once
        cl_length = pc.cl_length
        cl_sep = pc.cl_sep
        cpw_width = pc.cpw_width
        cpw_gap = pc.cpw_gap

        # For this design the loc_W has to be in 0 but loc_H can be -1 or +1,
        # For all other directions One can change the orientation of the qubit.
        loc_W = float(pc.loc_W)
        loc_W, loc_H = float(pc.loc_W), float(pc.loc_H)
        if float(loc_W) not in [0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a fluxonium qubit with loc_W is not 0 and'
                ' loc_H is not +1 or -1 ? Are you sure you want to do this?'
            )

        # Define the geometry
        # Charge Line
        charge_line = draw.rectangle(cpw_width, cl_length, 0, 0)
        charge_line_gap = draw.rectangle(
            cpw_width+2*cpw_gap, cl_length+cpw_gap/2, 0, -cpw_gap/4)

        # Making the charge_line and charge_line_gap circle and union them to charge_line and it's gap
        charge_line_round = draw.Point(0., (-cl_length)/2).buffer(cpw_width)
        charge_line_gap_round = draw.Point(
            0., -(cl_length+cpw_gap)/2.2).buffer((cpw_width+2*cpw_gap))
        charge_line = draw.union(charge_line, charge_line_round)
        charge_line_gap = draw.union(charge_line_gap, charge_line_gap_round)

        # Charge Line CPW wire
        port_line = draw.LineString([(-cpw_width/2, cl_length/2),
                                     (cpw_width/2, cl_length/2)])

        # Position the charge_line, rotate and translate
        objects = [charge_line, charge_line_gap, port_line]
        objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            0,
            loc_H * (p.pocket_height/2 + cl_sep + cl_length))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [charge_line, charge_line_gap, port_line] = objects

        self.add_qgeometry('poly', {'charge_line': charge_line})
        self.add_qgeometry(
            'poly', {'charge_line_gap': charge_line_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H == 1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('charge_line',
                     points, cpw_width)

    def make_readout_line(self):
        """ Adds readout line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pr = self.p.readout_line_options  # parser on readout line options

        # define commonly used variables once
        pad_sep = pr.pad_sep
        pad_width = pr.pad_width
        pad_height = pr.pad_height
        cpw_width = pr.cpw_width
        cpw_gap = pr.cpw_gap

        # For this design the loc_W has to be in 0 but loc_H can be -1 or +1,
        # For all other directions One can change the orientation of the qubit.
        loc_W = float(pr.loc_W)
        loc_W, loc_H = float(pr.loc_W), float(pr.loc_H)
        if float(loc_W) not in [0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a fluxonium qubit with loc_W is not 0 and'
                ' loc_H is not +1 or -1 ? Are you sure you want to do this?'
            )

        # Define the geometry
        # Readout pad
        readout_pad = draw.rectangle(pad_width, pad_height, 0, 0)
        # making the pad circle for left side
        readout_pad_circle_left = draw.Point(-pad_width/2,
                                             0).buffer(pad_height/2)
        # making the pad circle for right side
        readout_pad_circle_right = draw.Point(
            pad_width/2, 0).buffer(pad_height/2)

        # Readout pad's gap
        readout_pad_gap = draw.rectangle(
            pad_width+2*cpw_gap, pad_height+2*cpw_gap, 0, 0)
        readout_pad_gap_circle_left = draw.Point(-(pad_width+2*cpw_gap)/2, 0).buffer(
            (pad_height+2*cpw_gap)/2)  # making the pad's gap circle for left side
        readout_pad_gap_circle_right = draw.Point((pad_width+2*cpw_gap)/2, 0).buffer(
            (pad_height+2*cpw_gap)/2)  # making the pad's gap circle for right side

        # Defining the geometry for the readout pad line and it's gap
        readout_line = draw.rectangle(cpw_width, pad_height, 0, pad_height)
        readout_line_gap = draw.rectangle(
            cpw_width+2*cpw_gap, pad_height, 0, pad_height)

        # Here, we union the readout pad and readout line and second line exactly the same for the gap
        readout_padNline = draw.union(
            readout_pad, readout_pad_circle_left, readout_pad_circle_right, readout_line)
        readout_padNline_gap = draw.union(
            readout_pad_gap, readout_pad_gap_circle_left, readout_pad_gap_circle_right, readout_line_gap)

        # Readout Line CPW wire
        port_line = draw.LineString([(cpw_width/2, pad_height*1.5),
                                     (-cpw_width/2, pad_height*1.5)])

        # Position the readout, rotate and translate
        objects = [readout_padNline, readout_padNline_gap, port_line]
        objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            0,
            loc_H * (p.pocket_height/2 + pad_sep))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [readout_padNline, readout_padNline_gap, port_line] = objects

        self.add_qgeometry('poly', {'readout_padNline': readout_padNline})
        self.add_qgeometry(
            'poly', {'readout_padNline_gap': readout_padNline_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H == - \
            1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('readout_line',
                     points, cpw_width)


class FluxoniumPocket(_FluxoniumPocket):

    """The connectors to the inductor were shaped like the dorsal fin of a shark/dolphin.
    made them rectangular.

    Adding these things to Default Options to accomadate a connecting piece on the capacitor to the JJ.
    Has the option to make the pads purely rectangular with rounded sides

    Can accomadate a second flux line on the left side. Currently the two flux lines are mirror images of 
    each other but that can be changed in the future.

    This class inherits _FluxoniumPocket by Figen Yilmaz and Christian Kragland Andersen.
    This class was created by Alexander Nguyen and David Sommers.
    """
    default_options = Dict(_FluxoniumPocket.default_options,
                           top_wire_connector=False,  # by default, the connector won't exist
                           top_wire_center_x='0.0039mm',
                           top_wire_center_y='0.0097mm',
                           top_wire_height='0.012mm',
                           top_wire_width='0.0034mm',
                           bot_wire_connector=False,  # by default, the connector won't exist
                           bot_wire_center_x='-0.0039mm',
                           bot_wire_center_y='-0.0097mm',
                           bot_wire_height='0.015mm',
                           bot_wire_width='0.0039mm',
                           make_rol_left=False,  # make the mirror image of the flux line on the other side
                           # make the edges of the rectangular part round.
                           round_edge=False,
                           teeth_options=Dict(
                               make_teeth=False,  # make teeth in readout resonator pad. Will be FUNKY unless pad_radius=0
                               # using same defaults and variable names as transmon_pocket_teeth whenever possible
                               coupled_pad_width='10um',
                               coupled_pad_height='40um',  # changed this default to be a lot shorter
                               # distance between centers of the teeth. edge to edge distance depends on coupled_pad_width
                               coupled_pad_gap='100um',
                               pad_gap='15um'  # this is the gap between the CPW and the capacitor. NotE: because of rounded end the top of the rounded end is half a cpw_width closer than stated gap
                           ))

    def make(self):
        """Define the way the options are turned into QGeometry.

        The make function implements the logic that creates the geometry
        (poly, path, etc.) from the qcomponent.options dictionary of
        parameters, and the adds them to the design, using
        qcomponent.add_qgeometry(...), adding in extra needed
        information, such as layer, subtract, etc.
        """
        self.make_pocket()

        if self.p.flux_bias_line_options.make_fbl == True:
            self.make_flux_bias_line()
        if self.p.charge_line_options.make_cl == True:
            self.make_charge_line()
        if self.p.readout_line_options.make_rol == True:
            self.make_readout_line()

        if self.p.make_rol_left:
            self.make_flux_bias_line2()

    def make_pocket(self):
        """Makes standard fluxonium in a pocket."""

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p

        # since we will reuse these options, parse them once and define them as variables
        pocket_width = p.pocket_width
        pocket_height = p.pocket_height
        pad_height = p.pad_height
        pad_width = p.pad_width
        pad_radius = p.pad_radius
        pad_gap = p.pad_gap
        # length of the horizontal arms on the right hand side where JJ array connects
        l_arm_length = p.l_arm_length
        l_arm_width = p.l_arm_width
        l_width = p.l_width
        nanowire_inductor = p.nanowire
        teeth_options = p.teeth_options

        # Drawing the kinectic inductor
        if nanowire_inductor == True:
            # MY COMMENT: I don't know how l_length is defined here since it is only defined
            inductor = draw.LineString(
                [(l_arm_length, l_length/2), (l_arm_length, -l_length/2)])
            # in the else below. regardless, this just defines the inductor as a straight simple line.
            if p.round_edge:  # note that the inductor is backwards when using the default value of inductor_orientation of -1
                inductor = draw.LineString([(((pad_width+pad_height)/2)+(l_arm_length/4), (pad_gap+pad_height-l_arm_width)/2), ((
                    (pad_width+pad_height)/2)+(l_arm_length/4), -((pad_gap+pad_height-l_arm_width)/2))])  # ((pad_gap+pad_height)/2)-l_arm_width
        else:
            l_length = p.array_length
            # This one is for JJ array
            # MY COMMENT: This just makes the inductor longer/shorter, inductor_orientation is just a scaling factor
            io = float(p.inductor_orientation)
            inductor = draw.LineString(
                [(l_arm_length-l_arm_width, io*l_length/2), (l_arm_length-l_arm_width, -io*l_length/2)])
            if p.round_edge:
                inductor = draw.LineString([(((pad_width+pad_height)/2)+(l_arm_length/4), io*(pad_gap+pad_height-l_arm_width)/2), ((
                    (pad_width+pad_height)/2)+(l_arm_length/4), -io*((pad_gap+pad_height-l_arm_width)/2))])

        # Draw 'the arms' and make them curvy, first top arm and then same goes for the bottom
        l_arm_up = draw.Polygon([
            (pad_width/2, l_length/2+l_arm_width),  # point a
            (l_arm_length, l_length/2+l_arm_width),  # point b
            (l_arm_length, l_length/2),  # point c
            (pad_width/2, l_length/2),  # point d
        ])
        if p.round_edge:
            # l_arm_up=draw.translate(l_arm_up, (pad_height/2)+(l_arm_length/4), pad_height/2)
            l_arm_up = draw.rectangle(
                l_arm_length, l_arm_width, pad_width, (pad_gap+pad_height)/2)
            # self.add_qgeometry('poly', dict(l_arm_up=l_arm_up))

        """l_arm_up_fillet = draw.Point(l_arm_length, l_length/2).buffer(l_arm_width) # Having semicircle with subtracting the geometries
        cut_ply_up = draw.Polygon([
             (l_arm_length-l_arm_width*2, l_length/2+l_arm_width*2),   # point o
             (l_arm_length, l_length/2+l_arm_width*2),    # point p
             (l_arm_length, l_length/2),   # point r same with point c
             (l_arm_length+l_arm_width*2, l_length/2),   # point s
             (l_arm_length+l_arm_width*2, l_length/2-l_arm_width*2),  # point t
             (l_arm_length-l_arm_width*2, l_length/2-l_arm_width*2),  # point u
        ])
        l_arm_up_fillet = draw.subtract(l_arm_up_fillet, cut_ply_up) # Having semicircle with subtracting the geometries
        """
        # Here the bottom arm
        l_arm_bot = draw.Polygon([
            (pad_width/2, -l_length/2),  # point e
            (l_arm_length, -l_length/2),  # point f
            (l_arm_length, -(l_length/2+l_arm_width)),  # point g
            (pad_width/2, -(l_length/2+l_arm_width)),  # point h
        ])
        if p.round_edge:
            l_arm_bot = draw.rectangle(
                l_arm_length, l_arm_width, pad_width, -(pad_gap+pad_height)/2)
        """l_arm_bot_fillet = draw.Point(l_arm_length, -l_length/2).buffer(l_arm_width)
        cut_ply_bot = draw.Polygon([
             ((l_arm_length-l_arm_width*2), -(l_length/2+l_arm_width*2)),   # point o
             (l_arm_length, -(l_length/2+l_arm_width*2)),    # point p
             (l_arm_length, -l_length/2),   # point r same with point c
             (l_arm_length+l_arm_width*2, -l_length/2),   # point s
             (l_arm_length+l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point t
             (l_arm_length-l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point u
        ])
        l_arm_bot_fillet = draw.subtract(l_arm_bot_fillet, cut_ply_bot) # Having semicircle with subtracting the geometries"""

        # Draw the pads (shapely polygons)
        pad_rect_top = draw.rectangle(
            pad_width, pad_height, 0, (pad_gap+pad_height)/2)
        pad_circle_top = draw.Point(
            0, (pad_radius+pad_height)).buffer(pad_radius)
        pad_top = draw.union(pad_rect_top, pad_circle_top,
                             l_arm_up)  # , l_arm_up_fillet)
        if p.round_edge:
            # origin is middle of gap so go up half gap and half height
            pad_right_circle_top = draw.Point(
                pad_width/2, (pad_gap+pad_height)/2).buffer(pad_height/2)
            pad_left_circle_top = draw.Point(-(pad_width)/2,
                                             (pad_gap+pad_height)/2).buffer(pad_height/2)
            pad_top = draw.union(
                pad_top, pad_right_circle_top, pad_left_circle_top)
        pad_rect_bot = draw.rectangle(
            pad_width, pad_height, 0, -(pad_gap+pad_height)/2)
        pad_circle_bot = draw.Point(
            0, -(pad_radius+pad_height)).buffer(pad_radius)
        pad_bot = draw.union(pad_rect_bot, pad_circle_bot,
                             l_arm_bot)  # , l_arm_bot_fillet)
        if p.round_edge:
            pad_right_circle_bot = draw.Point(
                pad_width/2, -(pad_gap+pad_height)/2).buffer(pad_height/2)
            pad_left_circle_bot = draw.Point(-(pad_width) /
                                             2, -(pad_gap+pad_height)/2).buffer(pad_height/2)
            pad_bot = draw.union(
                pad_bot, pad_right_circle_bot, pad_left_circle_bot)

        if teeth_options.make_teeth:  # makes teeth to insert readout resonator. best if pad_radius=0
            tooth_left = draw.rectangle(teeth_options.coupled_pad_width,
                                        teeth_options.coupled_pad_height, -teeth_options.coupled_pad_gap/2, -(pad_height+((pad_gap+teeth_options.coupled_pad_height)/2)))
            coupler_pad_round_left = draw.Point(-teeth_options.coupled_pad_gap/2, -(teeth_options.coupled_pad_height + pad_height+(pad_gap/2))).buffer(
                teeth_options.coupled_pad_width / 2,
                resolution=16,
                cap_style=CAP_STYLE.round)
            tooth_right = draw.translate(
                tooth_left, teeth_options.coupled_pad_gap, 0)
            coupler_pad_round_right = draw.translate(
                coupler_pad_round_left, teeth_options.coupled_pad_gap, 0)
            pad_bot = draw.union(
                pad_bot, tooth_left, coupler_pad_round_left, tooth_right, coupler_pad_round_right)

        if teeth_options.make_teeth:#makes teeth to insert readout resonator. best if pad_radius=0
            tooth_left=draw.rectangle(teeth_options.coupled_pad_width,
                                     teeth_options.coupled_pad_height, -teeth_options.coupled_pad_gap/2, -(pad_height+((pad_gap+teeth_options.coupled_pad_height)/2)))
            coupler_pad_round_left = draw.Point(-teeth_options.coupled_pad_gap/2, -(teeth_options.coupled_pad_height + pad_height+(pad_gap/2))).buffer(
                teeth_options.coupled_pad_width / 2,
                resolution=16,
                cap_style=CAP_STYLE.round)
            tooth_right=draw.translate(tooth_left, teeth_options.coupled_pad_gap , 0)
            coupler_pad_round_right=draw.translate(coupler_pad_round_left, teeth_options.coupled_pad_gap , 0)
            pad_bot=draw.union(pad_bot, tooth_left, coupler_pad_round_left, tooth_right, coupler_pad_round_right)

        if p.top_wire_connector:
            connector_top = draw.rectangle(
                p.top_wire_width,
                p.top_wire_height,
                p.top_wire_center_x,
                p.top_wire_center_y
            )  # make "finger" for top capactior pad
            pad_top = draw.union(pad_top, connector_top)

        if p.bot_wire_connector:
            connector_bot = draw.rectangle(
                p.bot_wire_width,
                p.bot_wire_height,
                p.bot_wire_center_x,
                p.bot_wire_center_y
            )  # make "finger" for bot capactior pad
            pad_bot = draw.union(pad_bot, connector_bot)

        # Draw the junction
        # one can change the JJ orientation. Fab related detail.
        jj_o = float(p.jj_orientation)
        rect_jj = draw.LineString([(0, -pad_gap/2*jj_o), (0, +pad_gap/2*jj_o)])

        if p.top_wire_connector & p.bot_wire_connector:
            center_x = (p.bot_wire_center_x+p.top_wire_center_x) / \
                2  # trying to find the center between
            # "fingers" between the capacitors
            # bottom edge for transparent box
            bot_edge = p.bot_wire_center_y+(p.bot_wire_height/2)
            top_edge = p.top_wire_center_y-(p.top_wire_height/2)
            rect_jj = draw.LineString(
                [(center_x, bot_edge), (center_x, top_edge)])

            right_edge = p.top_wire_center_x+(p.top_wire_width/2)
            left_edge = p.bot_wire_center_x-(p.bot_wire_width/2)
            p.jj_width = right_edge-left_edge
            # if necessary, consider the case where there is only one connector on top or bottom

        # the draw.rectangle representing the josephson junction
        # rect_jj = draw.rectangle(p.jj_width, pad_gap)

        # Draw the pocket
        rect_pk = draw.rectangle(pocket_width, pocket_height)

        # Rotate and translate all qgeometry as needed.
        polys = [rect_jj, pad_top, pad_bot, rect_pk, inductor]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [rect_jj, pad_top, pad_bot, rect_pk, inductor] = polys

        # Use the geometry to create Metal qgeometry
        self.add_qgeometry('poly', dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry('poly', dict(rect_pk=rect_pk), subtract=True)
        # self.add_qgeometry('poly', dict(
        #     rect_jj=rect_jj), helper=True)
        self.add_qgeometry('junction',  # kinetic inductor
                           dict(inductor=inductor),
                           width=l_width,
                           hfss_inductance=p.l_inductance,
                           gds_cell_name=p.gds_cell_inductor)
        self.add_qgeometry('junction',  # the main JJ
                           dict(rect_jj=rect_jj),
                           width=p.jj_width,
                           hfss_inductance=p.L_j,
                           hfss_capacitance=p.C_j)

    def make_flux_bias_line2(self):
        """ Adds flux bias line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pfb = self.p.flux_bias_line_options  # parser on connector options

        # define commonly used variables once
        fbl_sep = pfb.fbl_sep
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        # Define the geometry
        # Flux Bias Line
        # The position of flux bias line on the x-axis, starting point inside the pocket
        d = p.pocket_width/2
        # Draw the top line of the flux-bias line
        flux_bias_lineup = draw.Polygon([
            (-d, fbl_height/2),   # point a
            (-d, fbl_height/2+cpw_width),    # point b
            (-fbl_sep, fbl_height/2+cpw_width),   # point c
            (-fbl_sep, fbl_height/2),   # point d
        ])
        # Draw the middle line of the flux-bias line
        flux_bias_linemid = draw.Polygon([
            (-(fbl_sep-cpw_width), fbl_height/2),   # point e
            (-fbl_sep, fbl_height/2),    # point f
            (-fbl_sep, -fbl_height/2),  # point g
            (-(fbl_sep-cpw_width), -fbl_height/2),   # point h

        ])

        # Here we make flux-bias line curvy for the top side and also bottom side and the union all of them
        circle_top = draw.Point(-fbl_sep, fbl_height/2).buffer(cpw_width)
        cut_ply = draw.Polygon([
            (-fbl_sep*2, fbl_height+cpw_width),   # point o
            (-fbl_sep, fbl_height+cpw_width),    # point p
            (-fbl_sep, fbl_height/2),   # point r same with point d
            (-fbl_sep/2, fbl_height/2),   # point s
            (-fbl_sep/2, -fbl_height+cpw_width),  # point t
            (-fbl_sep*2, -fbl_height+cpw_width),  # point u
        ])
        circle_top = draw.subtract(circle_top, cut_ply)
        # same goes for bottom edge
        circle_bot = draw.Point(-fbl_sep, -fbl_height/2).buffer(cpw_width)
        cut_ply2 = draw.Polygon([
            (-(fbl_sep-cpw_width), -fbl_height/2),   # point v same with h or i
            (-(fbl_sep-cpw_width), fbl_height*2),    # point y
            (-fbl_sep, fbl_height),   # point z
            (-fbl_sep, fbl_height),   # point w
            (-fbl_sep*2, fbl_height),  # point x
            (-fbl_sep*2, -fbl_height/2),  # point k
        ])
        circle_bot = draw.subtract(circle_bot, cut_ply2)
        flux_bias_linebot = draw.Polygon([
            (-(d+fbl_sep/2), -fbl_height/2),   # point i
            (-(d+fbl_sep/2), -fbl_height/2-cpw_width),    # point k
            (-fbl_sep, -fbl_height/2-cpw_width),   # point l
            (-fbl_sep, -fbl_height/2),   # point m
        ])
        flux_bias_line = draw.union(
            flux_bias_lineup, flux_bias_linemid,  flux_bias_linebot, circle_top, circle_bot)

        # Flux Bias line's gap part, inside the GND
        flux_bias_line_gap = draw.rectangle(
            fbl_sep/2, cpw_width+cpw_gap*2, -(d+fbl_sep/4), -fbl_height/2-cpw_width/2)

        # Flux-Bias Line CPW wire
        port_line = draw.LineString([(-(d+fbl_sep/2), 0),
                                    (-(d+fbl_sep/2), -(fbl_height+cpw_width))])

        # This port line is a fake port line, it is only in use during LOM analyses because we need to have an ungrounded line for the flux-bias
        fake_port_line = draw.LineString([(-d, (fbl_height*2+cpw_width*2)),
                                          (-d, -(fbl_height+cpw_width))])

        objects = [flux_bias_line, flux_bias_line_gap,
                   port_line, fake_port_line]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line,
            fake_port_line] = objects  # flux_bias_line,

        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})

        self.add_qgeometry(
            'poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)

        ####################################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line2',
                     port_line_cords, cpw_width)

        fake_port_line_cords = list(
            draw.shapely.geometry.shape(fake_port_line).coords)
        self.add_pin('fake_flux_bias_line2',
                     fake_port_line_cords, cpw_width)

    def make_readout_line(self):
        """ Adds readout line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pr = self.p.readout_line_options  # parser on readout line options
        teeth_options = p.teeth_options  # to make resonator go inside qubit pocket
        pocket_height = p.pocket_height
        pad_height_qubit = p.pad_height
        pad_gap = p.pad_gap  # height of the gap between two capacitor pads in qubit

        # define commonly used variables once
        pad_sep = pr.pad_sep
        pad_width = pr.pad_width
        pad_height = pr.pad_height
        cpw_width = pr.cpw_width
        cpw_gap = pr.cpw_gap

        # For this design the loc_W has to be in 0 but loc_H can be -1 or +1,
        # For all other directions One can change the orientation of the qubit.
        loc_W = float(pr.loc_W)
        loc_W, loc_H = float(pr.loc_W), float(pr.loc_H)
        if float(loc_W) not in [0] or float(loc_H) not in [-1., +1.]:
            self.logger.info(
                'Warning: Did you mean to define a fluxonium qubit with loc_W is not 0 and'
                ' loc_H is not +1 or -1 ? Are you sure you want to do this?'
            )

        # Define the geometry
        # Readout pad
        readout_pad = draw.rectangle(pad_width, pad_height, 0, 0)
        # making the pad circle for left side
        readout_pad_circle_left = draw.Point(-pad_width/2,
                                             0).buffer(pad_height/2)
        # making the pad circle for right side
        readout_pad_circle_right = draw.Point(
            pad_width/2, 0).buffer(pad_height/2)

        # Readout pad's gap
        readout_pad_gap = draw.rectangle(
            pad_width+2*cpw_gap, pad_height+2*cpw_gap, 0, 0)
        readout_pad_gap_circle_left = draw.Point(-(pad_width+2*cpw_gap)/2, 0).buffer(
            (pad_height+2*cpw_gap)/2)  # making the pad's gap circle for left side
        readout_pad_gap_circle_right = draw.Point((pad_width+2*cpw_gap)/2, 0).buffer(
            (pad_height+2*cpw_gap)/2)  # making the pad's gap circle for right side

        # Defining the geometry for the readout pad line and it's gap
        readout_line = draw.rectangle(cpw_width, pad_height, 0, pad_height)
        readout_line_gap = draw.rectangle(
            cpw_width+2*cpw_gap, pad_height, 0, pad_height)

        # Here, we union the readout pad and readout line and second line exactly the same for the gap
        readout_padNline = draw.union(
            readout_pad, readout_pad_circle_left, readout_pad_circle_right, readout_line)
        readout_padNline_gap = draw.union(
            readout_pad_gap, readout_pad_gap_circle_left, readout_pad_gap_circle_right, readout_line_gap)

        # Readout Line CPW wire
        port_line = draw.LineString([(cpw_width/2, pad_height*1.5),
                                     (-cpw_width/2, pad_height*1.5)])

        # Position the readout, rotate and translate
        objects = [readout_padNline, readout_padNline_gap, port_line]
        objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            0,
            loc_H * (p.pocket_height/2 + pad_sep))
        objects = draw.rotate_position(objects, p.orientation,
                                       [p.pos_x, p.pos_y])
        [readout_padNline, readout_padNline_gap, port_line] = objects

        # if we are using teeth and putting the resonator inside the qubit pocket get rid of readout pad
        if not teeth_options.make_teeth:
            self.add_qgeometry('poly', {'readout_padNline': readout_padNline})
            self.add_qgeometry(
                'poly', {'readout_padNline_gap': readout_padNline_gap}, subtract=True)
        else:  # the readout line has to be in the middle i.e. colinear with y axis when rotation=0
            readout_line = draw.rectangle(cpw_width,
                                          ((pocket_height-pad_gap)/2) -
                                          pad_height_qubit-teeth_options.pad_gap,
                                          0,
                                          ((pocket_height+pad_gap)/4)+((pad_height_qubit+teeth_options.pad_gap)/2))  # the center of the rectangle moves up half a pocket_height-half height of rectangle
            # readout_line_gap = draw.rectangle(cpw_width+2*cpw_gap, pad_height, 0, pad_height)
            port_line = draw.LineString([(cpw_width/2, pocket_height/2),
                                         (-cpw_width/2, pocket_height/2)])
            readout_line_circle = draw.Point(
                0, (pad_gap/2)+pad_height_qubit+teeth_options.pad_gap).buffer(cpw_width/2)
            readout_line = draw.union(readout_line, readout_line_circle)

            objects = [readout_line, port_line]
            # reflects the shape across the x-axis if loc_H is negative
            objects = draw.scale(objects, 1, loc_H, origin=(0, 0))
            """objects = draw.translate(
                objects,
                0,
                loc_H * (p.pocket_height/2 + pad_sep))"""
            objects = draw.rotate_position(objects, p.orientation,
                                           [p.pos_x, p.pos_y])
            [readout_line, port_line] = objects
            self.add_qgeometry('poly', {'readout_line': readout_line})

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H == - \
            1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('readout_line',
                     points, cpw_width)
