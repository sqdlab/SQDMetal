# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 01/05/2023
# Description: Collection of classes to create joints/pins to be used in routing/positioning.

from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
from qiskit_metal.toolbox_python.attr_dict import Dict
from SQDMetal.Utilities.QUtilities import QUtilities


class Joint(QComponent):
    """Creates a point pin at a specified location.

    Inherits QComponent class.

    The positioning is specified via: (pos_x, pos_y) with the pin facing the direction given by 'orientation'.
    The 'width' parameter specifies the pin width.

    Pins:
        The solitary pin is called 'a'. Its width matches the supplied 'width' parameter.

    .. image::
        Cap3Interdigital.png

    .. meta::
        A pin at a point.

    Default Options:
        * pos_x='0um', pos_y='0um', orientation=0
        * width='1um'
    """

    default_options = Dict(pos_x='0um', pos_y='0um', orientation=0, width='1um')

    def make(self):
        p = self.p

        ptJoint = np.array([p.pos_x, p.pos_y])
        norm_vec = np.array([-np.cos(p.orientation/180*np.pi), -np.sin(p.orientation/180*np.pi)])

        self.add_pin('a', [(ptJoint+norm_vec).tolist(), ptJoint.tolist()], width=p.width, input_as_norm=True)

class JointExtend(QComponent):
    """Creates a point pin some distance away from another pin. Useful when needing to connect a component to at some offset
    from a point.

    Inherits QComponent class.

    The pin on the target component is specified via:
        * jointObj - Name of target component in the design to which the pin is to attach.
        * jointPin - Name of the pin on the target component. If left as blank (i.e. ''), the component's centre (pos_x, pos_y)
                     will be taken.
        * pin_width - Width associated with pin. Ignored if jointPin is specified
    
    The positioning is specified via:
        * dist_extend - Distance in which to place the pin away from the target component pin
        * orientation - Direction in which to place the pin away from the target component pin
        * extend_off_pin_dir - If True, orientation is ignored and the direction matches the target pin's normal direction.
        * pin_orientation - If None, this option is ignored and the resulting pin direction is given by the direction vector
                            (specified by either 'orientation' or 'extend_off_pin_dir'). Otherwise, the resulting pin direction
                            is given by 'pin_orientation'.
    This class ignores pos_x and pos_y...
        
    Pins:
        The solitary pin is called 'a'. Its width matches the path width of the target component's pin.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Pin placed some distance away from some other pin.

    Default Options:
        * jointObj=''
        * jointPin='a'
        * orientation=0
        * dist_extend='10um'
        * pin_orientation=None
        * extend_off_pin_dir=False
        * pin_width='10um'
    """

    default_options = Dict(jointObj='', jointPin='a', orientation=0, dist_extend='10um', pin_orientation=None, extend_off_pin_dir=False, pin_width='10um')

    def make(self):
        p = self.p

        if self.options.jointPin == '':
            ptJoint = np.array([self._design.components[self.options.jointObj].p.pos_x, self._design.components[self.options.jointObj].p.pos_y])
            angl = self._design.components[self.options.jointObj].p.orientation
            norm = np.array([np.cos(angl/180*np.pi), np.sin(angl/180*np.pi)])
            width = p.pin_width
        else:
            jointPin = self._design.components[self.options.jointObj].pins[self.options.jointPin]
            ptJoint = jointPin['middle']
            norm = jointPin['normal']
            width = jointPin['width']
        if p.extend_off_pin_dir:
            pt = ptJoint + p.dist_extend * norm
        else:
            pt = ptJoint + p.dist_extend * np.array([np.cos(p.orientation/180*np.pi), np.sin(p.orientation/180*np.pi)])

        self.options.pos_x, self.options.pos_y = pt
        if self.options.pin_orientation != None:
            ptJoint = pt - np.array([np.cos(self.options.pin_orientation/180*np.pi), np.sin(self.options.pin_orientation/180*np.pi)])

        self.add_pin('a', [ptJoint.tolist(), pt.tolist()], width=width, input_as_norm=True)

class RouteJoint(QComponent):
    """Creates a point pin positioned on a routing path object. Useful when needing to connect a component to some point
    along a path.

    Inherits QComponent class.

    The path object is specified via:
        * pathObj  - Name of component in the design to which the pin is to attach.
        * pathObjTraceName - Name of the trace (can be omitted if it's a QRoute object) that forms the subgeometry of the
                             object. For example, if a component uses QRoute paths to draw sections of capacitor legs, then
                             the name of this subgeometry must be given here.

    The positioning along the path is specified via:
        * frac_line - The fraction along the path (from 0 to 1) from the start to which to attach the pin
        * dist_line - Distance along the path from the start to which to attach the pin. Note that this is only valid if frac_line
                      is given as a number outside the interval [0,1].
        * is_right_hand  - If True, the orientation of the pin is normal to the path and faces the right hand side; left hand
                           side if given as False.
        * attach_on_side - If True, the pin is offset to match the width of the wire/path so that it's touching the side of
                           the path. If False, the pin position is right on the centre of the path.
    This class ignores pos_x, pos_y and orientation...

    On a side note while using attach_on_side; it works great on straight edges as an adjoining wire can sit flush on its side.
    Evidently, this does not work well on curved edges - in which case, attach_on_side should be set to False so that the wire
    connecting to the path will reach its centre and leave a flush connection when rendered/merged.
        
    Pins:
        The solitary pin is called 'a'. Its width matches the path width.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Pin placed along some point on a path.

    Default Options:
        * pathObj=''
        * pathObjTraceName=''
        * frac_line=0.5
        * dist_line=0
        * is_right_hand=True
        * attach_on_side=False
    """

    default_options = Dict(pathObj='', pathObjTraceName='', frac_line=0.5, dist_line=0, is_right_hand=True, attach_on_side=False)
    
    def make(self):
        p = self.p

        if self.options.pathObj == '':
            ptLine = np.array([0,0])
            norm_vec = np.array([1,0])
            width = 1e-6
        else:
            if p.frac_line < 0 or p.frac_line > 1:
                leDist = p.dist_line
                leFrac = False
            else:
                leDist = p.frac_line
                leFrac = True

            final_pts, normals, width, gap, total_dist = QUtilities.calc_points_on_path([leDist], self._design, self.options.pathObj, trace_name=self.options.pathObjTraceName, dists_are_fractional=leFrac)
            ptLine = final_pts[0]
            norm_vec = normals[0]
            if not p.is_right_hand:
                norm_vec = -norm_vec        
        if p.attach_on_side:
            if p.is_right_hand:
                ptLine += normals[0]*width*0.5
            else:
                ptLine -= normals[0]*width*0.5
        
        self.options.pos_x, self.options.pos_y = ptLine

        self.add_pin('a', [(ptLine-norm_vec).tolist(), ptLine.tolist()], width=width, input_as_norm=True)
