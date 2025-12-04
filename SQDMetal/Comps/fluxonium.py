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
"""Fluxonium Pocket"""

from operator import length_hint
import numpy as np
from qiskit_metal import draw, Dict
from math import *
from qiskit_metal.draw.basic import buffer
from qiskit_metal.qlibrary.core import BaseQubit



        
        
class FluxoniumPocket(BaseQubit):
    """The base `FluxoniumPocket` class.

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
        FluxoniumPocket.png
            
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

    component_metadata = Dict(short_name='FluxoniumPocket',
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
        l_arm_width = '2um',
        l_arm_length='20um',
        l_inductance='218.823nH',
        L_j = '40.7885nH',
        C_j = '0.72fF',
        inductor_orientation='-1',
        pocket_width='900um',
        pocket_height='550um',
        nanowire_inductor='True',
        gds_cell_inductor='gds_cell_inductor',
        # 90 has dipole aligned along the +X axis,
        # while 0 has dipole aligned along the +Y axis
        orientation='0',
        flux_bias_line_options=Dict(
            make_fbl = False,
            fbl_sep='85um',
            fbl_height ='50um',
            cpw_width ='10um',
            cpw_gap = '11.233um',
        ),
        charge_line_options=Dict(
            make_cl = False,
            cl_length = '100um',
            cl_sep ='-15um',
            cpw_width='cpw_width',
            cpw_gap= 'cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='+1',  # height location only +1 or -1
        ),
        readout_line_options=Dict(
            make_rol = False,
            pad_sep='85um',
            pad_width = '400um',
            pad_height = '120um',
            cpw_width='cpw_width',
            cpw_gap='cpw_gap',
            loc_W='0',  # width location only 0,
            loc_H='-1',  # height location only -1 or +1
        ))
    """Default drawing options"""

    TOOLTIP = """The base `FluxoniumPocket` class."""

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
            inductor = draw.LineString([(l_arm_length, l_length/2), (l_arm_length, -l_length/2)])#MY COMMENT: I don't know how l_length is defined here since it is only defined
            #in the else below. regardless, this just defines the inductor as a straight simple line.
        else:
            l_length = p.array_length
            # This one is for JJ array
            io = float(p.inductor_orientation)#MY COMMENT: This just makes the inductor longer/shorter, inductor_orientation is just a scaling factor
            inductor = draw.LineString([(l_arm_length-l_arm_width, io*l_length/2), (l_arm_length-l_arm_width, -io*l_length/2)])

        # Draw 'the arms' and make them curvy, first top arm and then same goes for the bottom
        l_arm_up = draw.Polygon([
            (pad_width/2, l_length/2+l_arm_width), # point a
            (l_arm_length, l_length/2+l_arm_width), # point b
            (l_arm_length, l_length/2), # point c
            (pad_width/2, l_length/2), # point d
            ])
        l_arm_up_fillet = draw.Point(l_arm_length, l_length/2).buffer(l_arm_width) # Having semicircle with subtracting the geometries
        cut_ply_up = draw.Polygon([
             (l_arm_length-l_arm_width*2, l_length/2+l_arm_width*2),   # point o
             (l_arm_length, l_length/2+l_arm_width*2),    # point p
             (l_arm_length, l_length/2),   # point r same with point c
             (l_arm_length+l_arm_width*2, l_length/2),   # point s
             (l_arm_length+l_arm_width*2, l_length/2-l_arm_width*2),  # point t
             (l_arm_length-l_arm_width*2, l_length/2-l_arm_width*2),  # point u
        ])
        l_arm_up_fillet = draw.subtract(l_arm_up_fillet, cut_ply_up) # Having semicircle with subtracting the geometries
        # Here the bottom arm
        l_arm_bot = draw.Polygon([
            (pad_width/2, -l_length/2), # point e
            (l_arm_length, -l_length/2), # point f
            (l_arm_length, -(l_length/2+l_arm_width)), # point g
            (pad_width/2, -(l_length/2+l_arm_width)), # point h
            ])
        l_arm_bot_fillet = draw.Point(l_arm_length, -l_length/2).buffer(l_arm_width)
        cut_ply_bot = draw.Polygon([
             ((l_arm_length-l_arm_width*2), -(l_length/2+l_arm_width*2)),   # point o
             (l_arm_length, -(l_length/2+l_arm_width*2)),    # point p
             (l_arm_length, -l_length/2),   # point r same with point c
             (l_arm_length+l_arm_width*2, -l_length/2),   # point s
             (l_arm_length+l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point t
             (l_arm_length-l_arm_width*2, -(l_length/2-l_arm_width*2)),  # point u
        ])
        l_arm_bot_fillet = draw.subtract(l_arm_bot_fillet, cut_ply_bot) # Having semicircle with subtracting the geometries

        # Draw the pads (shapely polygons)
        pad_rect_top = draw.rectangle(pad_width, pad_height, 0, (pad_gap+pad_height)/2)
        pad_circle_top = draw.Point(0, (pad_radius+pad_height)).buffer(pad_radius)
        pad_top = draw.union(pad_rect_top, pad_circle_top, l_arm_up, l_arm_up_fillet)
        pad_rect_bot = draw.rectangle(pad_width, pad_height, 0, -(pad_gap+pad_height)/2)
        pad_circle_bot = draw.Point(0, -(pad_radius+pad_height)).buffer(pad_radius)
        pad_bot = draw.union(pad_rect_bot, pad_circle_bot, l_arm_bot, l_arm_bot_fillet)

        # Draw the junction
        jj_o = float(p.jj_orientation) # one can change the JJ orientation. Fab related detail.
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
        self.add_qgeometry('junction', # kinetic inductor
                           dict(inductor=inductor),
                           width = l_width,
                           hfss_inductance = p.l_inductance,
                           gds_cell_name=p.gds_cell_inductor)
        self.add_qgeometry('junction', # the main JJ
                           dict(rect_jj=rect_jj),
                           width=p.jj_width,
                           hfss_inductance = p.L_j,
                           hfss_capacitance = p.C_j)

    def make_flux_bias_line(self):
        """ Adds flux bias line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
 
        p = self.p
        pfb = self.p.flux_bias_line_options # parser on connector options

        # define commonly used variables once
        fbl_sep = pfb.fbl_sep
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        # Define the geometry
        # Flux Bias Line
        d = p.pocket_width/2 # The position of flux bias line on the x-axis, starting point inside the pocket
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
             (fbl_sep, fbl_height+cpw_width ),    # point p
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
        flux_bias_line = draw.union(flux_bias_lineup, flux_bias_linemid,  flux_bias_linebot, circle_top, circle_bot)

        # Flux Bias line's gap part, inside the GND
        flux_bias_line_gap = draw.rectangle(fbl_sep/2, cpw_width+cpw_gap*2, d+fbl_sep/4 ,-fbl_height/2-cpw_width/2)

        # Flux-Bias Line CPW wire
        port_line = draw.LineString([((d+fbl_sep/2), 0), 
                                    ((d+fbl_sep/2), -(fbl_height+cpw_width))])
        
        # This port line is a fake port line, it is only in use during LOM analyses because we need to have an ungrounded line for the flux-bias 
        fake_port_line = draw.LineString([(d, (fbl_height*2+cpw_width*2)), 
                                    (d, -(fbl_height+cpw_width))])

        objects = [flux_bias_line, flux_bias_line_gap, port_line, fake_port_line]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line, fake_port_line] = objects

        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})
        self.add_qgeometry('poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)      

        ####################################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line', 
                    port_line_cords, cpw_width)
        
        fake_port_line_cords = list(draw.shapely.geometry.shape(fake_port_line).coords)
        self.add_pin('fake_flux_bias_line', 
                    fake_port_line_cords, cpw_width)

    def make_charge_line(self):
        """ Adds charge line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.charge_line_options # parser on charge line options
 
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
        charge_line_gap = draw.rectangle(cpw_width+2*cpw_gap, cl_length+cpw_gap/2, 0, -cpw_gap/4)

        # Making the charge_line and charge_line_gap circle and union them to charge_line and it's gap
        charge_line_round = draw.Point(0., (-cl_length)/2).buffer(cpw_width)
        charge_line_gap_round = draw.Point(0., -(cl_length+cpw_gap)/2.2).buffer((cpw_width+2*cpw_gap))
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
        self.add_qgeometry('poly', {'charge_line_gap': charge_line_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H==1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('charge_line', 
                    points, cpw_width)

    
    def make_readout_line(self):
        """ Adds readout line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pr = self.p.readout_line_options # parser on readout line options

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
        readout_pad_circle_left = draw.Point(-pad_width/2, 0).buffer(pad_height/2) # making the pad circle for left side
        readout_pad_circle_right = draw.Point(pad_width/2, 0).buffer(pad_height/2) # making the pad circle for right side

        # Readout pad's gap
        readout_pad_gap = draw.rectangle(pad_width+2*cpw_gap, pad_height+2*cpw_gap, 0, 0)
        readout_pad_gap_circle_left = draw.Point(-(pad_width+2*cpw_gap)/2, 0).buffer((pad_height+2*cpw_gap)/2) # making the pad's gap circle for left side
        readout_pad_gap_circle_right = draw.Point((pad_width+2*cpw_gap)/2, 0).buffer((pad_height+2*cpw_gap)/2) # making the pad's gap circle for right side
        
        # Defining the geometry for the readout pad line and it's gap
        readout_line = draw.rectangle(cpw_width, pad_height, 0, pad_height)
        readout_line_gap = draw.rectangle(cpw_width+2*cpw_gap, pad_height, 0, pad_height)

        # Here, we union the readout pad and readout line and second line exactly the same for the gap
        readout_padNline = draw.union(readout_pad, readout_pad_circle_left, readout_pad_circle_right, readout_line)
        readout_padNline_gap = draw.union(readout_pad_gap, readout_pad_gap_circle_left, readout_pad_gap_circle_right, readout_line_gap)
    
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
        self.add_qgeometry('poly', {'readout_padNline_gap': readout_padNline_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords if loc_H==-1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('readout_line', 
                    points, cpw_width)


class FluxoniumPocketSQD(FluxoniumPocket):

    """The connectors to the inductor were shaped like the dorsal fin of a shark/dolphin.
    made them rectangular.

    Adding these things to Default Options to accomadate a connecting piece on the capacitor to the JJ
    """
    default_options = Dict(FluxoniumPocket.default_options,
                           top_wire_connector=False,#by default, the connector won't exist
                           top_wire_center_x='0.0039mm',
                           top_wire_center_y='0.0097mm',
                           top_wire_height='0.012mm',
                           top_wire_width='0.0034mm',
                           bot_wire_connector=False,#by default, the connector won't exist
                           bot_wire_center_x='-0.0039mm',
                           bot_wire_center_y='-0.0097mm',
                           bot_wire_height='0.015mm',
                           bot_wire_width='0.0039mm',
                           make_rol_left=False#make the mirror image of the flux line on the other side
                           )
    
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
        l_arm_length = p.l_arm_length
        l_arm_width = p.l_arm_width
        l_width = p.l_width
        nanowire_inductor = p.nanowire

        # Drawing the kinectic inductor
        if nanowire_inductor == True:
            inductor = draw.LineString([(l_arm_length, l_length/2), (l_arm_length, -l_length/2)])#MY COMMENT: I don't know how l_length is defined here since it is only defined
            #in the else below. regardless, this just defines the inductor as a straight simple line.
        else:
            l_length = p.array_length
            # This one is for JJ array
            io = float(p.inductor_orientation)#MY COMMENT: This just makes the inductor longer/shorter, inductor_orientation is just a scaling factor
            inductor = draw.LineString([(l_arm_length-l_arm_width, io*l_length/2), (l_arm_length-l_arm_width, -io*l_length/2)])

        # Draw 'the arms' and make them curvy, first top arm and then same goes for the bottom
        l_arm_up = draw.Polygon([
            (pad_width/2, l_length/2+l_arm_width), # point a
            (l_arm_length, l_length/2+l_arm_width), # point b
            (l_arm_length, l_length/2), # point c
            (pad_width/2, l_length/2), # point d
            ])
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
            (pad_width/2, -l_length/2), # point e
            (l_arm_length, -l_length/2), # point f
            (l_arm_length, -(l_length/2+l_arm_width)), # point g
            (pad_width/2, -(l_length/2+l_arm_width)), # point h
            ])
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
        pad_rect_top = draw.rectangle(pad_width, pad_height, 0, (pad_gap+pad_height)/2)
        pad_circle_top = draw.Point(0, (pad_radius+pad_height)).buffer(pad_radius)
        pad_top = draw.union(pad_rect_top, pad_circle_top, l_arm_up)#, l_arm_up_fillet)
        pad_rect_bot = draw.rectangle(pad_width, pad_height, 0, -(pad_gap+pad_height)/2)
        pad_circle_bot = draw.Point(0, -(pad_radius+pad_height)).buffer(pad_radius)
        pad_bot = draw.union(pad_rect_bot, pad_circle_bot, l_arm_bot)#, l_arm_bot_fillet)

        if p.top_wire_connector:
            connector_top = draw.rectangle(
                  p.top_wire_width,
                  p.top_wire_height,
                  p.top_wire_center_x,
                  p.top_wire_center_y
                  )#make "finger" for top capactior pad
            pad_top=draw.union(pad_top, connector_top)

        if p.bot_wire_connector:
            connector_bot = draw.rectangle(
                  p.bot_wire_width,
                  p.bot_wire_height,
                  p.bot_wire_center_x,
                  p.bot_wire_center_y
                  )#make "finger" for bot capactior pad
            pad_bot=draw.union(pad_bot, connector_bot)


        # Draw the junction
        jj_o = float(p.jj_orientation) # one can change the JJ orientation. Fab related detail.
        rect_jj = draw.LineString([(0, -pad_gap/2*jj_o), (0, +pad_gap/2*jj_o)])

        if p.top_wire_connector & p.bot_wire_connector:
            center_x=(p.bot_wire_center_x+p.top_wire_center_x)/2#trying to find the center between
            #"fingers" between the capacitors
            bot_edge = p.bot_wire_center_y+(p.bot_wire_height/2)#bottom edge for transparent box
            top_edge = p.top_wire_center_y-(p.top_wire_height/2)
            rect_jj = draw.LineString([(center_x, bot_edge), (center_x, top_edge)])

            right_edge = p.top_wire_center_x+(p.top_wire_width/2)
            left_edge = p.bot_wire_center_x-(p.bot_wire_width/2)
            p.jj_width = right_edge-left_edge
            #if necessary, consider the case where there is only one connector on top or bottom

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
        self.add_qgeometry('junction', # kinetic inductor
                           dict(inductor=inductor),
                           width = l_width,
                           hfss_inductance = p.l_inductance,
                           gds_cell_name=p.gds_cell_inductor)
        self.add_qgeometry('junction', # the main JJ
                           dict(rect_jj=rect_jj),
                           width=p.jj_width,
                           hfss_inductance = p.L_j,
                           hfss_capacitance = p.C_j)
    
    def make_flux_bias_line2(self):
        """ Adds flux bias line to fluxonium pocket."""
        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pfb = self.p.flux_bias_line_options # parser on connector options

        # define commonly used variables once
        fbl_sep = pfb.fbl_sep
        fbl_height = pfb.fbl_height
        cpw_width = pfb.cpw_width
        cpw_gap = pfb.cpw_gap

        # Define the geometry
        # Flux Bias Line
        d = p.pocket_width/2 # The position of flux bias line on the x-axis, starting point inside the pocket
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
             (-fbl_sep, -fbl_height/2),# point g
             (-(fbl_sep-cpw_width), -fbl_height/2),   # point h
             
        ])

        # Here we make flux-bias line curvy for the top side and also bottom side and the union all of them
        circle_top = draw.Point(-fbl_sep, fbl_height/2).buffer(cpw_width)
        cut_ply = draw.Polygon([
             (-fbl_sep*2, fbl_height+cpw_width),   # point o
             (-fbl_sep, fbl_height+cpw_width ),    # point p
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
        flux_bias_line = draw.union(flux_bias_lineup, flux_bias_linemid,  flux_bias_linebot, circle_top, circle_bot)

        # Flux Bias line's gap part, inside the GND
        flux_bias_line_gap = draw.rectangle(fbl_sep/2, cpw_width+cpw_gap*2, -(d+fbl_sep/4) ,-fbl_height/2-cpw_width/2)

        # Flux-Bias Line CPW wire
        port_line = draw.LineString([(-(d+fbl_sep/2), 0), 
                                    (-(d+fbl_sep/2), -(fbl_height+cpw_width))])
        
        # This port line is a fake port line, it is only in use during LOM analyses because we need to have an ungrounded line for the flux-bias 
        fake_port_line = draw.LineString([(-d, (fbl_height*2+cpw_width*2)), 
                                    (-d, -(fbl_height+cpw_width))])

        objects = [flux_bias_line, flux_bias_line_gap, port_line, fake_port_line]#
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(objects, p.pos_x, p.pos_y)
        [flux_bias_line, flux_bias_line_gap, port_line, fake_port_line] = objects#flux_bias_line,

        self.add_qgeometry('poly', {'flux_bias_line': flux_bias_line})



        self.add_qgeometry('poly', {'flux_bias_line_gap': flux_bias_line_gap}, subtract=True)      

        ####################################################################

        # add pins
        """"port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        self.add_pin('flux_bias_line', 
                    port_line_cords, cpw_width)
        
        fake_port_line_cords = list(draw.shapely.geometry.shape(fake_port_line).coords)
        self.add_pin('fake_flux_bias_line', 
                    fake_port_line_cords, cpw_width)"""