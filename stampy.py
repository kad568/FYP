# abaqus imports
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
from odbAccess import openOdb
from abaqusConstants import *
from symbolicConstants import *
from odbAccess import *
# other imports
from dataclasses import dataclass, asdict
from json import dump
import time
import os
import copy
import shutil
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from enum import Enum
import meshio
import pyvista as pv
# from bayes_opt import BayesianOptimization, JSONLogger, Events
import sys
import json

SCRIPT_PARENT_PATH = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build"

class SolverType(Enum):

    STANDARD = 0
    EXPLICIT = 1


@dataclass
class InputeDex:

    simulation_output_path: str = None
    solver_type: SolverType = SolverType.STANDARD # options are "explicit" or "standard"
    all_part_rotation: int = 90 # degrees

    cup_height: float = None # mm 

    # blank
    blank_radius: float = None # mm
    blank_thickness: float = None # mm
    integration_points: int = None

    trim_depth: float = 2 # mm

    # die
    die_height: float = None # mm
    die_profile_radius: float = None # mm
    die_min_radius: float = None # mm
    die_max_radius: float = None # mm

    # blank holder
    blank_holder_height: float = None # mm
    blank_holder_profile_radius: float = None # mm
    blank_holder_min_radius: float = None # mm
    blank_holder_max_radius: float = None # mm
    blank_holder_die_gap: float = None # mm, minimum blank thickness to prevent overlap

    # punch
    cup_design_height: float = None # mm
    punch_depth: float = None # mm
    punch_profile_radius: float = None # mm
    punch_min_radius: float = None # mm

    # blank material
    blank_material_name: str = None
    density: float = None # kg/mm^3
    youngs_modulus: float = None # MPa
    posissons_ratio: float = None 
    plastic_material_data: tuple = None # (stress [MPa], strain)

    # BCs
    punch_velocity: float = None # mm/s
    punch_depth: float = None # mm
    mass_scalling: float = None
    BHF: float = None # N

    # interaction properties
    friction_coefficient: float = None

    # mesh
    mesh_size: float = None # mm
    local_mesh_size: float = None # mm

    # ideal part
    ideal_cup_radius: float = None # mm
    ideal_cup_height: float = None # mm
    ideal_cup_profile_radius: float = None # mm


def create_blank_part(blank_radius, blank_thickness, part_rotation):

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(0.0, blank_thickness/2), point2=(blank_radius, blank_thickness/2))
    s1.HorizontalConstraint(entity=g[3], addUndoState=False)
    p = mdb.models['Model-1'].Part(name="blank", dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts["blank"]
    p.BaseShellRevolve(sketch=s1, angle=part_rotation, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts["blank"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    
    p = mdb.models['Model-1'].parts['blank']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=119.354, 
        farPlane=210.572, width=49.2805, height=20.9298, viewOffsetX=0.194489, 
        viewOffsetY=4.84156)
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    p.deleteMesh(regions=pickedRegions)
    p = mdb.models['Model-1'].parts['blank']
    f, e1, d2 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=e1[1], 
        sketchPlaneSide=SIDE1, origin=(0, blank_thickness/2, 0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=113.14, gridSpacing=2.82, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['blank']
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.CircleByCenterPerimeter(center=(0, 0), point1=(10, 
        0)) # change inner circle radius
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    e, d1 = p.edges, p.datums
    p.PartitionFaceBySketch(sketchUpEdge=e[1], faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    if part_rotation == 360:

        strip_start_radius = 20
        strip_end_radius = 25

        p = mdb.models['Model-1'].parts['blank']
        f1, e1, d2 = p.faces, p.edges, p.datums
        t = p.MakeSketchTransform(sketchPlane=f1[0], sketchUpEdge=e1[3], 
            sketchPlaneSide=SIDE1, origin=(0.0, blank_thickness/2, 0.0))
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
            sheetSize=310.7, gridSpacing=7.76, transform=t)
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        s1.setPrimaryObject(option=SUPERIMPOSE)
        p = mdb.models['Model-1'].parts['blank']
        p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
        s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(strip_start_radius, 0))
        s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(strip_end_radius, 0))
        p = mdb.models['Model-1'].parts['blank']
        f = p.faces
        pickedFaces = f.getSequenceFromMask(mask=('[#3 ]', ), )
        e, d1 = p.edges, p.datums
        p.PartitionFaceBySketch(sketchUpEdge=e[3], faces=pickedFaces, sketch=s1)
        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

        p = mdb.models['Model-1'].parts['blank']
        f = p.faces
        pickedRegions = f.getSequenceFromMask(mask=('[#8 ]', ), )
        p = mdb.models['Model-1'].parts['blank']
        f, e, d1 = p.faces, p.edges, p.datums
        t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=e[7], 
            sketchPlaneSide=SIDE1, origin=(0.0, blank_thickness/2, 0.0))
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
            sheetSize=225.65, gridSpacing=5.64, transform=t)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=SUPERIMPOSE)
        p = mdb.models['Model-1'].parts['blank']
        p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

        s.Line(point1=(0.0, 0.0), point2=(-1 * blank_radius, 0))
        s.VerticalConstraint(entity=g[21], addUndoState=False)
        p = mdb.models['Model-1'].parts['blank']
        f = p.faces
        pickedFaces = f.getSequenceFromMask(mask=('[#f ]', ), )
        e1, d2 = p.edges, p.datums
        p.PartitionFaceBySketch(sketchUpEdge=e1[9], faces=pickedFaces, sketch=s)
        s.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']



def create_blank_surfaces(part_rotation):

    if part_rotation == 360:

        p = mdb.models['Model-1'].parts['blank']
        s = p.faces
        side2Faces = s.getSequenceFromMask(mask=('[#ff ]', ), )
        p.Surface(side2Faces=side2Faces, name='blank_top_surface')
        p = mdb.models['Model-1'].parts['blank']
        s = p.faces
        side1Faces = s.getSequenceFromMask(mask=('[#ff ]', ), )
        p.Surface(side1Faces=side1Faces, name='blank_bottom_surface')

    else:

        p = mdb.models['Model-1'].parts['blank']
        s = p.faces
        side2Faces = s.getSequenceFromMask(mask=('[#f ]', ), )
        p.Surface(side2Faces=side2Faces, name='blank_top_surface')
        p = mdb.models['Model-1'].parts['blank']
        s = p.faces
        side1Faces = s.getSequenceFromMask(mask=('[#f ]', ), )
        p.Surface(side1Faces=side1Faces, name='blank_bottom_surface')


def create_die_part(die_height, die_profile_radius, die_min_radius, die_max_radius, part_rotation):

    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s.FixedConstraint(entity=g[2])
    s.Line(point1=(die_min_radius, -1 * die_height), point2=(die_min_radius, 0.0))
    s.VerticalConstraint(entity=g[3], addUndoState=False)
    s.Line(point1=(die_min_radius, 0.0), point2=(die_max_radius, 0.0))
    s.HorizontalConstraint(entity=g[4], addUndoState=False)
    s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s.FilletByRadius(radius=die_profile_radius, curve1=g[3], nearPoint1=(die_min_radius, 
        (-1 * die_height)/2), curve2=g[4], nearPoint2=((die_max_radius+ die_min_radius)/2, 
        0.0))
    p = mdb.models['Model-1'].Part(name="die", dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts["die"]
    p.BaseShellRevolve(sketch=s, angle=part_rotation, flipRevolveDirection=OFF)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts["die"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


def create_die_surface(all_part_rotation):

    p1 = mdb.models['Model-1'].parts["die"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts["die"]
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#7 ]', ), )
    if all_part_rotation == 360:
        p.Surface(side1Faces=side1Faces, name="die_surface")
    else:
        p.Surface(side2Faces=side1Faces, name="die_surface")



def die_ref_point(all_part_rotation):

    p = mdb.models['Model-1'].parts["die"]
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes

    if all_part_rotation == 360:
        p.ReferencePoint(point=v1[0])
    else:
        p.ReferencePoint(point=v1[6])


def create_blank_holder(blank_holder_height, blank_holder_profile_radius, blank_holder_min_radius, blank_holder_max_radius, blank_holder_die_gap, part_rotation):

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(blank_holder_min_radius, blank_holder_height+blank_holder_die_gap), point2=(blank_holder_min_radius, blank_holder_die_gap))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    s1.Line(point1=(blank_holder_min_radius, blank_holder_die_gap), point2=(blank_holder_max_radius, blank_holder_die_gap))
    s1.HorizontalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.FilletByRadius(radius=blank_holder_profile_radius, curve1=g[3], nearPoint1=(blank_holder_min_radius, 
        (blank_holder_height+ 2* blank_holder_die_gap)/2), curve2=g[4], nearPoint2=((blank_holder_max_radius + blank_holder_min_radius)/2, 
        blank_holder_die_gap))
    p = mdb.models['Model-1'].Part(name="blank_holder", dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts["blank_holder"]
    p.BaseShellRevolve(sketch=s1, angle=part_rotation, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts["blank_holder"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


def create_blank_holder_surface():

    p1 = mdb.models['Model-1'].parts["blank_holder"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts["blank_holder"]
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#7 ]', ), )
    p.Surface(side1Faces=side1Faces, name="blank_holder_surface")


def blank_holder_ref_point():

    p = mdb.models['Model-1'].parts["blank_holder"]
    v2, e1, d2, n1 = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=p.InterestingPoint(edge=e1[1], rule=MIDDLE))


def create_punch(punch_depth, punch_profile_radius, punch_min_radius, blank_thickness, part_rotation):

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(punch_min_radius, punch_depth + blank_thickness), point2=(punch_min_radius, blank_thickness))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    s1.Line(point1=(punch_min_radius, blank_thickness), point2=(0.0, blank_thickness))
    s1.HorizontalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.FilletByRadius(radius=punch_profile_radius, curve1=g[3], nearPoint1=(punch_min_radius, 
        (punch_depth + 2 * blank_thickness)/2), curve2=g[4], nearPoint2=(punch_min_radius/2, 
        blank_thickness))
    p = mdb.models['Model-1'].Part(name="punch", dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts["punch"]
    p.BaseShellRevolve(sketch=s1, angle=part_rotation, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts["punch"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

def create_ideal_part(cup_height, ideal_cup_profile_radius, ideal_cup_radius, part_rotation):

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(ideal_cup_radius, cup_height), point2=(ideal_cup_radius, 0))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    s1.Line(point1=(ideal_cup_radius, 0), point2=(0.0, 0))
    s1.HorizontalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.FilletByRadius(radius=ideal_cup_profile_radius, curve1=g[3], nearPoint1=(ideal_cup_radius, 
        cup_height/2), curve2=g[4], nearPoint2=(ideal_cup_radius/2, 
        0))
    p = mdb.models['Model-1'].Part(name="idealpart", dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts["idealpart"]
    p.BaseShellRevolve(sketch=s1, angle=part_rotation, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts["idealpart"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

    # save mesh
    p = mdb.models['Model-1'].parts["idealpart"]
    p.seedPart(size=0.2, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['Model-1'].parts['idealpart']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    print(f'{os.getcwd()}/ideal_part.obj')
    session.writeOBJFile(
        fileName=f'{os.getcwd()}/ideal_part.obj', 
        canvasObjects= (session.viewports['Viewport: 1'], ))



def punch_ref_point():

    p = mdb.models['Model-1'].parts["punch"]
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=v1[0])


def create_punch_surface(all_part_rotation):

    p = mdb.models['Model-1'].parts["punch"]
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#7 ]', ), )

    if all_part_rotation == 360:
        p.Surface(side1Faces=side1Faces, name="punch_surface")
    else:
        p.Surface(side2Faces=side1Faces, name="punch_surface")

def create_stopper(all_part_rotation, max_blank_holder_radius, blank_thickness):

    stopper_size = 4

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.rectangle(point1=(max_blank_holder_radius - stopper_size, blank_thickness * 0.5), point2=(max_blank_holder_radius, 0.0))
    p = mdb.models['Model-1'].Part(name='stopper', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['stopper']
    p.BaseSolidRevolve(sketch=s1, angle=all_part_rotation, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['stopper']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

    # # ref point
    # a = mdb.models['Model-1'].rootAssembly
    # e11 = a.instances['stopper-1'].edges
    # a.ReferencePoint(point=a.instances['stopper-1'].InterestingPoint(edge=e11[1], 
    #     rule=MIDDLE))


def create_stopper_surface():

    p = mdb.models['Model-1'].parts['stopper']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#8 ]', ), )
    p.Surface(side1Faces=side1Faces, name='stopper_top_surf')


    p = mdb.models['Model-1'].parts['stopper']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#2 ]', ), )
    p.Surface(side1Faces=side1Faces, name='stopper_bottom_surf')


def create_part_assembly():

    a1 = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts["blank"]
    a1.Instance(name=f'{"blank"}-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts["blank_holder"]
    a1.Instance(name=f'{"blank_holder"}-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts["die"]
    a1.Instance(name=f'{"die"}-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts["punch"]
    a1.Instance(name=f'{"punch"}-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts["stopper"]
    a1.Instance(name=f'{"stopper"}-1', part=p, dependent=ON)


def create_blank_material(blank_material_name, density, youngs_modulus, posissons_ratio, plastic_material_data):

    mdb.models['Model-1'].Material(name=blank_material_name)
    mdb.models['Model-1'].materials[blank_material_name].Density(table=((density, ), 
        ))
    mdb.models['Model-1'].materials[blank_material_name].Elastic(table=((youngs_modulus, 
        posissons_ratio), ))
    mdb.models['Model-1'].materials[blank_material_name].Plastic(scaleStress=None, 
        table=plastic_material_data)
    
    # stopper material
    mdb.models['Model-1'].Material(name='stopper_material', 
        objectToCopy=mdb.models['Model-1'].materials[blank_material_name])
    mdb.models['Model-1'].materials['stopper_material'].elastic.setValues(table=((
        youngs_modulus * 10, posissons_ratio), ))


def create_surface_interactions(solver_type: SolverType, friction_coefficient):



    # create interaction properties
    mdb.models['Model-1'].ContactProperty('die_blank')
    mdb.models['Model-1'].ContactProperty('punch_blank')
    mdb.models['Model-1'].ContactProperty('blank_holder_blank')
    mdb.models['Model-1'].ContactProperty('die_blank_holder_stopper')

    mdb.models['Model-1'].interactionProperties['die_blank'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((friction_coefficient, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    
    mdb.models['Model-1'].interactionProperties['punch_blank'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((friction_coefficient, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)

    mdb.models['Model-1'].interactionProperties['blank_holder_blank'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((friction_coefficient, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)

    # die blank interaction
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['die-1'].surfaces['die_surface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['blank-1'].surfaces['blank_bottom_surface']
    if solver_type == SolverType.EXPLICIT:
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='die_blank', 
            createStepName='Initial', main = region1, secondary = region2, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='die_blank', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
    elif solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name ='die_blank', 
            createStepName='Initial', main = region1, secondary = region2, sliding=FINITE, 
            interactionProperty='die_blank', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
    
    # punch blank interaction
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['punch-1'].surfaces['punch_surface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['blank-1'].surfaces['blank_top_surface']
    if solver_type == SolverType.EXPLICIT:
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='punch_blank', 
            createStepName='Initial', main = region1, secondary = region2, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='punch_blank', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)
    elif solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name ='punch_blank', 
            createStepName='Initial', main = region1, secondary = region2, sliding=FINITE, 
            interactionProperty='punch_blank', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)
    # blank holder blank interaction
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['blank_holder-1'].surfaces['blank_holder_surface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['blank-1'].surfaces['blank_top_surface']
    if solver_type == SolverType.EXPLICIT:
            
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='blank_holder_blank', 
            createStepName='Initial', main = region1, secondary = region2, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='blank_holder_blank', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)
    elif solver_type == SolverType.STANDARD:

        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name ='blank_holder_blank', 
                createStepName='Initial', main = region1, secondary = region2, sliding=FINITE, 
                interactionProperty='blank_holder_blank', initialClearance=OMIT, 
                datumAxis=None, clearanceRegion=None)


    # blank_die interaction
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['die-1'].surfaces['die_surface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['blank-1'].surfaces['blank_bottom_surface']
    if solver_type == SolverType.EXPLICIT:
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='die_blank', 
            createStepName='Initial', main = region1, secondary = region2, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='die_blank', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
    elif solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name ='die_blank', 
            createStepName='Initial', main = region1, secondary = region2, sliding=FINITE, 
            interactionProperty='die_blank', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
        
    # # blank_holder_die contact
    # a = mdb.models['Model-1'].rootAssembly
    # region1=a.instances['blank_holder-1'].surfaces['blank_holder_surface']
    # a = mdb.models['Model-1'].rootAssembly
    # region2=a.instances['die-1'].surfaces['die_surface']
    # if solver_type == SolverType.STANDARD:
    #     mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='blank_holder_die', 
    #         createStepName='Initial', main=region1, secondary=region2, 
    #         sliding=FINITE, thickness=ON, interactionProperty='blank_holder_die', 
    #         adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
    #         clearanceRegion=None)
    # elif solver_type == SolverType.EXPLICIT:   
    #     mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='blank_holder_die', 
    #         createStepName='Initial', main = region2, secondary = region1, 
    #         mechanicalConstraint=PENALTY, sliding=FINITE, 
    #         interactionProperty='blank_holder_die', initialClearance=OMIT, datumAxis=None, 
    #         clearanceRegion=None)


    # blank_holder_stopper contact
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['blank_holder-1'].surfaces['blank_holder_surface']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['stopper-1'].surfaces['stopper_top_surf']
    if solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='blank_holder_stopper', 
            createStepName='Initial', main=region1, secondary=region2, 
            sliding=FINITE, thickness=ON, interactionProperty='die_blank_holder_stopper', 
            adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
    elif solver_type == SolverType.EXPLICIT:   
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='blank_holder_stopper', 
            createStepName='Initial', main = region2, secondary = region1, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='die_blank_holder_stopper', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)

    # die_stopper contact
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['stopper-1'].surfaces['stopper_bottom_surf']
    a = mdb.models['Model-1'].rootAssembly
    region2=a.instances['die-1'].surfaces['die_surface']
    if solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='die_stopper', 
            createStepName='Initial', main=region2, secondary=region1, 
            sliding=FINITE, thickness=ON, interactionProperty='die_blank_holder_stopper', 
            adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)
    elif solver_type == SolverType.EXPLICIT:   
        mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='die_stopper', 
            createStepName='Initial', main = region2, secondary = region1, 
            mechanicalConstraint=PENALTY, sliding=FINITE, 
            interactionProperty='die_blank_holder_stopper', initialClearance=OMIT, datumAxis=None, 
            clearanceRegion=None)

        
def create_boundary_conditions(solver_type, punch_speed, punch_depth, mass_scaling, trim_buffer, all_parts_rotation, BHF):

    load_time = (punch_depth + trim_buffer) / punch_speed

    punch_release_speed = punch_speed
    blank_holder_release_speed = punch_release_speed * 0.5

    release_height = (punch_depth + trim_buffer) / 2
    unload_time = release_height / punch_release_speed

    if solver_type == SolverType.EXPLICIT:

        mdb.models['Model-1'].ExplicitDynamicsStep(name='load', previous='Initial', 
            timePeriod=load_time, improvedDtMethod=ON)
        mdb.models['Model-1'].steps['load'].setValues(massScaling=((SEMI_AUTOMATIC, 
            MODEL, AT_BEGINNING, 0.0, mass_scaling, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None), 
            ), improvedDtMethod=ON)

    elif solver_type == SolverType.STANDARD:
        mdb.models['Model-1'].StaticStep(
            name='load', 
            previous='Initial',
            timePeriod=load_time,
            maxNumInc=1, # 1000
            initialInc=0.1, 
            maxInc=1.0, 
            nlgeom=ON
        )

        mdb.models['Model-1'].steps['load'].Restart(
            frequency=1,
            overlay=OFF
        )

    if all_parts_rotation == 90:

        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['blank-1'].edges
        edges1 = e1.getSequenceFromMask(mask=('[#18 ]', ), )
        region = a.Set(edges=edges1, name='Set-1')
        mdb.models['Model-1'].XsymmBC(name='BC-1', createStepName='load', 
            region=region, localCsys=None)
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['blank-1'].edges

        edges1 = e1.getSequenceFromMask(mask=('[#22 ]', ), )
        region = a.Set(edges=edges1, name='Set-2')
        mdb.models['Model-1'].ZsymmBC(name='BC-2', createStepName='load', 
            region=region, localCsys=None)
        a = mdb.models['Model-1'].rootAssembly
        r1 = a.instances['die-1'].referencePoints
        refPoints1=(r1[3], )
    
    elif all_parts_rotation == 180:
        
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['blank-1'].edges
        edges1 = e1.getSequenceFromMask(mask=('[#4c8 ]', ), )
        region = a.Set(edges=edges1, name='Set-1')
        mdb.models['Model-1'].ZsymmBC(name='BC-1', createStepName='load', 
            region=region, localCsys=None)
        a = mdb.models['Model-1'].rootAssembly
        e1 = a.instances['blank-1'].edges

        edges1 = e1.getSequenceFromMask(mask=('[#a22 ]', ), )
        region = a.Set(edges=edges1, name='Set-2')
        mdb.models['Model-1'].ZsymmBC(name='BC-2', createStepName='load', 
            region=region, localCsys=None)
        a = mdb.models['Model-1'].rootAssembly
        r1 = a.instances['die-1'].referencePoints
        refPoints1=(r1[3], )

    elif all_parts_rotation == 360:

        pass

        # a = mdb.models['Model-1'].rootAssembly
        # v1 = a.instances['blank-1'].vertices
        # verts1 = v1.getSequenceFromMask(mask=('[#40 ]', ), )
        # region = a.Set(vertices=verts1, name='Set-1')
        # mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='load', 
        #     region=region, u1=0.0, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        #     amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        #     localCsys=None)


    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['die-1'].referencePoints
    refPoints1=(r1[3], )

    region = a.Set(referencePoints=refPoints1, name='Set-3')
    mdb.models['Model-1'].DisplacementBC(name='BC-3', createStepName='load', 
        region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['blank_holder-1'].referencePoints
    refPoints1=(r1[3], )

    region = a.Set(referencePoints=refPoints1, name='Set-4')
    mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='load', 
        region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['punch-1'].referencePoints
    refPoints1=(r1[2], )

    region = a.Set(referencePoints=refPoints1, name='Set-5')
    mdb.models['Model-1'].DisplacementBC(name='BC-5', createStepName='load', 
        region=region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['punch-1'].referencePoints
    refPoints1=(r1[2], )

    region = a.Set(referencePoints=refPoints1, name='Set-6')
    mdb.models['Model-1'].VelocityBC(name='BC-6', createStepName='load', 
        region=region, v1=UNSET, v2= -1 * punch_speed, v3=UNSET, vr1=UNSET, vr2=UNSET, 
        vr3=UNSET, amplitude=UNSET, localCsys=None, distributionType=UNIFORM, 
        fieldName='')
    
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['blank_holder-1'].referencePoints
    refPoints1=(r1[3], )
    region = regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].ConcentratedForce(name='BHF', createStepName='load', 
        region=region, cf2=-1 * BHF, distributionType=UNIFORM, field='', 
        localCsys=None)
    mdb.models['Model-1'].boundaryConditions['BC-4'].setValues(u2=UNSET)

    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['stopper-1'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#10 ]', ), )
    region = regionToolset.Region(faces=faces1)
    mdb.models['Model-1'].XsymmBC(name='BC-77', createStepName='Initial', 
        region=region, localCsys=None)

    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['stopper-1'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#20 ]', ), )
    region = regionToolset.Region(faces=faces1)
    mdb.models['Model-1'].ZsymmBC(name='BC-88', createStepName='Initial', 
        region=region, localCsys=None)

def create_mesh(solver_type: SolverType, mesh_size, local_seeding_size):

    if solver_type == SolverType.EXPLICIT:
        
        elemType1_rb = mesh.ElemType(elemCode=R3D4, elemLibrary=EXPLICIT)
        elemType2_rb = mesh.ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)
        elemType1_shell = mesh.ElemType(elemCode=S4R, elemLibrary=EXPLICIT, secondOrderAccuracy=ON)
        elemType2_shell = mesh.ElemType(elemCode=S3R, elemLibrary=EXPLICIT)
    elif solver_type == SolverType.STANDARD:

        elemType1_rb = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elemType2_rb = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        elemType1_shell = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
        elemType2_shell = mesh.ElemType(elemCode=S3R, elemLibrary=STANDARD)

    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p1 = mdb.models['Model-1'].parts['die']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)

    p = mdb.models['Model-1'].parts['die']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#7 ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1_rb, elemType2_rb))
    p = mdb.models['Model-1'].parts['die']
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models['Model-1'].parts['die']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#7 ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP)
    p = mdb.models['Model-1'].parts['die']
    p.generateMesh()
    p = mdb.models['Model-1'].parts['blank_holder']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p = mdb.models['Model-1'].parts['blank_holder']
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)


    p = mdb.models['Model-1'].parts['blank_holder']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#7 ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1_rb, elemType2_rb))
    p = mdb.models['Model-1'].parts['blank_holder']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#7 ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP)
    p = mdb.models['Model-1'].parts['blank_holder']
    p.generateMesh()
    p = mdb.models['Model-1'].parts['blank']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p = mdb.models['Model-1'].parts['blank']
    p.seedPart(size=local_seeding_size, deviationFactor=0.1, minSizeFactor=0.1)
    edges = p.edges.getByBoundingBox(-5, -5, -10, 5, 5, 10) 
    p.seedEdgeBySize(edges=edges, size=local_seeding_size, constraint=FIXED) # removed the / 3 for mesh size



    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1_shell,))
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP)
    p = mdb.models['Model-1'].parts['blank']
    p.generateMesh()

    # # local seeding blank
    # p = mdb.models['Model-1'].parts['blank']
    # f = p.faces
    # faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    # pickedRegions =(faces, )
    # p.setElementType(regions=pickedRegions, elemTypes=(elemType1_shell, elemType2_shell))
    # p = mdb.models['Model-1'].parts['blank']
    # f = p.faces
    # pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    # p.setMeshControls(regions=pickedRegions, technique=SWEEP)
    # p = mdb.models['Model-1'].parts['blank']
    # p.seedPart(size=local_seeding_size, deviationFactor=0.1, minSizeFactor=0.1)
    # p.generateMesh()

    p1 = mdb.models['Model-1'].parts['punch']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts['punch']
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)

    p = mdb.models['Model-1'].parts['punch']
    p.generateMesh()
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#7 ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1_rb, elemType2_rb))
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON, mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
        engineeringFeatures=OFF, mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)

    # stopper
    p = mdb.models['Model-1'].parts['stopper']
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, 
        kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD) # only standard available for stopper
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['stopper']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(cells, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))
    p = mdb.models['Model-1'].parts['stopper']
    c = p.cells
    pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP, 
        algorithm=ADVANCING_FRONT)
    p = mdb.models['Model-1'].parts['stopper']
    p.generateMesh()

    ## local seeding
    p = mdb.models['Model-1'].parts['punch']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
    p.deleteMesh(regions=pickedRegions)
    p = mdb.models['Model-1'].parts['punch']
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#8 ]', ), )
    p.seedEdgeBySize(edges=pickedEdges, size=local_seeding_size, deviationFactor=0.1, 
        minSizeFactor=0.1, constraint=FINER)
    p = mdb.models['Model-1'].parts['punch']
    p.generateMesh()

    # die local seeding
    p = mdb.models['Model-1'].parts['die']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
    p.deleteMesh(regions=pickedRegions)
    p = mdb.models['Model-1'].parts['die']
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#10 ]', ), )
    p.seedEdgeBySize(edges=pickedEdges, size=local_seeding_size, deviationFactor=0.1, 
        minSizeFactor=0.1, constraint=FINER)
    p = mdb.models['Model-1'].parts['die']
    p.generateMesh()

    # blankholder local seeding
    p = mdb.models['Model-1'].parts['blank_holder']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
    p.deleteMesh(regions=pickedRegions)
    p = mdb.models['Model-1'].parts['blank_holder']
    e = p.edges
    pickedEdges = e.getSequenceFromMask(mask=('[#10 ]', ), )
    p.seedEdgeBySize(edges=pickedEdges, size=local_seeding_size, deviationFactor=0.1, 
        minSizeFactor=0.1, constraint=FINER)
    p = mdb.models['Model-1'].parts['blank_holder']
    p.generateMesh()

    
def apply_material_properties(material_name, integration_points, blank_thickness, part_rotation):

    p = mdb.models['Model-1'].parts['blank']
    mdb.models['Model-1'].HomogeneousShellSection(name='blank_section', 
        preIntegrate=OFF, material=material_name, thicknessType=UNIFORM, 
        thickness=blank_thickness, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=integration_points)
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces

    # set stopper section
    p = mdb.models['Model-1'].parts['stopper']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    p.Set(cells=cells, name='stopper_set')
    mdb.models['Model-1'].HomogeneousSolidSection(name='stopper_section', 
        material='stopper_material', thickness=None)

    # set stopper material
    p = mdb.models['Model-1'].parts['stopper']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(cells=cells)
    p = mdb.models['Model-1'].parts['stopper']
    p.SectionAssignment(region=region, sectionName='stopper_section', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    a = mdb.models['Model-1'].rootAssembly

    if part_rotation == 360:
        faces = f.getSequenceFromMask(mask=('[#ff ]', ), )
    else:
        faces = f.getSequenceFromMask(mask=('[#f ]', ), )

    region = regionToolset.Region(faces=faces)
    p = mdb.models['Model-1'].parts['blank']
    p.SectionAssignment(region=region, sectionName='blank_section', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    regionDef=mdb.models['Model-1'].rootAssembly.sets['Set-5']
    mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-2', 
        createStepName='load', variables=('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3','RF1', 'RF2', 'RF3', 'RM1', 'RM2', 
        'RM3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
    
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    p.Set(faces=faces, name='blank_set')
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    mdb.models['Model-1'].FieldOutputRequest(name='F-Output-2', 
        createStepName='load', variables=('S','MISES', 'PEEQ', 'STH'), region=MODEL, 
        exteriorOnly=OFF, sectionPoints=DEFAULT, rebar=EXCLUDE)

    
def run_sim_job(solver_type: SolverType, sim_out_path, nCPU):

    job_name = "stamping_sim"

    if solver_type == SolverType.EXPLICIT:
        mdb.Job(
            name=job_name, model='Model-1', description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
            memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE, 
            nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch=sim_out_path,
            resultsFormat=ODB, numCpus=nCPU, numDomains=nCPU, parallelizationMethodExplicit=DOMAIN, numGPUs=1
        )
    elif solver_type == SolverType.STANDARD:
        mdb.Job(
        #     name=job_name, model='Model-1', 
        #     description='', type=ANALYSIS, 
        #     scratch=sim_out_path,
        #     parallelizationMethodExplicit=THREADS,
        #     numCpus=nCPU, numDomains=1, 
        #     numGPUs=0
        name=job_name, model='Model-1', description='', type=ANALYSIS, 
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                scratch=sim_out_path, resultsFormat=ODB, numThreadsPerMpiProcess=1, 
                multiprocessingMode=DEFAULT, numCpus=nCPU, numDomains=nCPU, numGPUs=1
        )

    mdb.jobs[job_name].submit()

    # Wait for the Abaqus job to complete before proceeding
    mdb.jobs[job_name].waitForCompletion()

    with open("stamping_sim.sta", "r") as f:
        lines = f.readlines()

    last_line = lines[-1]

    if "SUCCESSFULLY" in last_line:
        output = 1
    else:
        output = 0
    
    with open("sim_status.txt", "w") as f:
        f.write(str(output))


def spring_back_analysis(sim_out_path, nCPU):

    mdb.Model(name='Model-1-spring_back', objectToCopy=mdb.models['Model-1'])
    p = mdb.models['Model-1-spring_back'].parts['blank']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1-spring_back'].parts['blank_holder']
    del mdb.models['Model-1-spring_back'].parts['die']
    del mdb.models['Model-1-spring_back'].parts['punch']
    del mdb.models['Model-1-spring_back'].parts['stopper']
    a = mdb.models['Model-1-spring_back'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a = mdb.models['Model-1-spring_back'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF, optimizationTasks=OFF, 
        geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON)
    a = mdb.models['Model-1-spring_back'].rootAssembly
    a.deleteFeatures(('blank_holder-1', 'die-1', 'punch-1', "stopper-1"))
    del mdb.models['Model-1-spring_back'].interactions['blank_holder_blank']
    del mdb.models['Model-1-spring_back'].interactions['die_blank']
    del mdb.models['Model-1-spring_back'].interactions['punch_blank']
    del mdb.models['Model-1-spring_back'].interactions['die_stopper']
    del mdb.models['Model-1-spring_back'].interactions['blank_holder_stopper']

    del mdb.models['Model-1-spring_back'].boundaryConditions['BC-77']
    del mdb.models['Model-1-spring_back'].boundaryConditions['BC-88']

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
        constraints=OFF, connectors=OFF, engineeringFeatures=OFF, 
        adaptiveMeshConstraints=ON)
    del mdb.models['Model-1-spring_back'].steps['load']
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    mdb.models['Model-1-spring_back'].StaticStep(name='Step-1', previous='Initial', 
        maxNumInc=10000, initialInc=0.001, maxInc=0.1, nlgeom=ON)
    
    a = mdb.models['Model-1-spring_back'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#18 ]', ), )
    region = a.Set(edges=edges1, name='Set-9')
    mdb.models['Model-1-spring_back'].XsymmBC(name='BC-1', createStepName='Step-1', 
        region=region, localCsys=None)
    a = mdb.models['Model-1-spring_back'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#22 ]', ), )
    region = a.Set(edges=edges1, name='Set-10')
    mdb.models['Model-1-spring_back'].ZsymmBC(name='BC-2', createStepName='Step-1', 
        region=region, localCsys=None)

    a = mdb.models['Model-1-spring_back'].rootAssembly
    v1 = a.instances['blank-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#10 ]', ), )
    region = a.Set(vertices=verts1, name='Set-11')
    mdb.models['Model-1-spring_back'].DisplacementBC(name='BC-3', 
        createStepName='Step-1', region=region, u1=0.0, u2=0.0, u3=0.0, 
        ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)
    instances=(mdb.models['Model-1-spring_back'].rootAssembly.instances['blank-1'], 
        )

    mdb.models['Model-1'].FieldOutputRequest(name='F-Output-1', 
        createStepName='load', variables=('S','MISES', 'PEEQ', 'U'), region=MODEL, 
        exteriorOnly=OFF, sectionPoints=DEFAULT, rebar=EXCLUDE)


    mdb.models['Model-1-spring_back'].InitialState(
        updateReferenceConfiguration=ON, fileName='stamping_sim', 
        endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
        createStepName='Initial', instances=instances)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF)
    mdb.Job(name='spring_back', model='Model-1-spring_back', description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch=sim_out_path, 
        resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=nCPU, numDomains=nCPU, numGPUs=1)
    mdb.jobs['spring_back'].submit(consistencyChecking=OFF)

    mdb.jobs['spring_back'].waitForCompletion()


def demeri_test(sim_out_path, nCPU):

    mdb.Model(name='Model-demeri', objectToCopy=mdb.models['Model-1'])
    p = mdb.models['Model-demeri'].parts['blank']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-demeri'].parts['blank_holder']
    del mdb.models['Model-demeri'].parts['die']
    del mdb.models['Model-demeri'].parts['punch']
    a = mdb.models['Model-demeri'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    a = mdb.models['Model-demeri'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF, optimizationTasks=OFF, 
        geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
        constraints=ON, connectors=ON, engineeringFeatures=ON)
    a = mdb.models['Model-demeri'].rootAssembly
    a.deleteFeatures(('blank_holder-1', 'die-1', 'punch-1', ))
    del mdb.models['Model-demeri'].interactions['blank_holder_blank']
    del mdb.models['Model-demeri'].interactions['die_blank']
    del mdb.models['Model-demeri'].interactions['punch_blank']
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
        constraints=OFF, connectors=OFF, engineeringFeatures=OFF, 
        adaptiveMeshConstraints=ON)
    del mdb.models['Model-demeri'].steps['load']
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    mdb.models['Model-demeri'].StaticStep(name='Step-1', previous='Initial', 
        maxNumInc=10000, initialInc=0.001, maxInc=0.1, nlgeom=ON)


    # a = mdb.models['Model-1-spring_back'].rootAssembly
    # e1 = a.instances['blank-1'].edges
    # edges1 = e1.getSequenceFromMask(mask=('[#30 ]', ), )
    # region = a.Set(edges=edges1, name='Set-7')
    # a = mdb.models['Model-1-spring_back'].rootAssembly
    # e1 = a.instances['blank-1'].edges
    # edges1 = e1.getSequenceFromMask(mask=('[#30 ]', ), )
    # region = a.Set(edges=edges1, name='Set-9')
    # mdb.models['Model-1-spring_back'].XsymmBC(name='BC-1', createStepName='Step-1', 
    #     region=region, localCsys=None)
    # a = mdb.models['Model-1-spring_back'].rootAssembly
    # e1 = a.instances['blank-1'].edges
    # edges1 = e1.getSequenceFromMask(mask=('[#44 ]', ), )
    # region = a.Set(edges=edges1, name='Set-10')
    # mdb.models['Model-1-spring_back'].ZsymmBC(name='BC-2', createStepName='Step-1', 
    #     region=region, localCsys=None)

    # a = mdb.models['Model-1-spring_back'].rootAssembly
    # v1 = a.instances['blank-1'].vertices
    # verts1 = v1.getSequenceFromMask(mask=('[#20 ]', ), )
    # region = a.Set(vertices=verts1, name='Set-11')
    # mdb.models['Model-1-spring_back'].DisplacementBC(name='BC-3', 
    #     createStepName='Step-1', region=region, u1=0.0, u2=0.0, u3=0.0, 
    #     ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
    #     distributionType=UNIFORM, fieldName='', localCsys=None)

    instances=(mdb.models['Model-demeri'].rootAssembly.instances['blank-1'], 
        )


    mdb.models['Model-demeri'].InitialState(
        updateReferenceConfiguration=ON, fileName='stamping_sim', 
        endStep=LAST_STEP, endIncrement=STEP_END, name='Predefined Field-1', 
        createStepName='Initial', instances=instances)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF)


def run_sim(input_dex: InputeDex):
    
    # reset abaqus

    Mdb()
    session.viewports['Viewport: 1'].setValues(displayedObject=None)

    # set work directory
    os.chdir(input_dex.simulation_output_path)

    if input_dex.solver_type == SolverType.EXPLICIT:
        input_dex.solver_type = "explicit"
    elif input_dex.solver_type == SolverType.STANDARD:
        input_dex.solver_type = "standard"

    # save simulation param inputs
    input_dex_output_path = f"{input_dex.simulation_output_path}/inputs.json"
    with open(input_dex_output_path, "w") as file:
        dump(asdict(input_dex), file)

    if input_dex.solver_type == "explicit":
        input_dex.solver_type = SolverType.EXPLICIT
    elif input_dex.solver_type == "standard":
        input_dex.solver_type = SolverType.STANDARD

    # create blank
    create_blank_part(input_dex.blank_radius, input_dex.blank_thickness, input_dex.all_part_rotation)
    create_blank_surfaces(input_dex.all_part_rotation)

    # create die
    create_die_part(input_dex.die_height, input_dex.die_profile_radius, input_dex.die_min_radius, input_dex.die_max_radius, input_dex.all_part_rotation)
    create_die_surface(input_dex.all_part_rotation)
    die_ref_point(input_dex.all_part_rotation)  

    # create bkank holder
    create_blank_holder(input_dex.blank_holder_height, input_dex.blank_holder_profile_radius, input_dex.blank_holder_min_radius, input_dex.blank_holder_max_radius, input_dex.blank_holder_die_gap, input_dex.all_part_rotation)
    create_blank_holder_surface()
    blank_holder_ref_point()

    # create punch
    create_punch(input_dex.punch_depth, input_dex.punch_profile_radius, input_dex.punch_min_radius, input_dex.blank_thickness, input_dex.all_part_rotation)
    punch_ref_point()
    create_punch_surface(input_dex.all_part_rotation)

    create_stopper(input_dex.all_part_rotation, input_dex.blank_holder_max_radius, input_dex.blank_thickness)
    create_stopper_surface()


    # create ideal part
    create_ideal_part(input_dex.ideal_cup_height, input_dex.ideal_cup_profile_radius, input_dex.ideal_cup_radius, input_dex.all_part_rotation)

    # assemble parts
    create_part_assembly()

    # define material
    create_blank_material(input_dex.blank_material_name, input_dex.density, input_dex.youngs_modulus, input_dex.posissons_ratio, input_dex.plastic_material_data)

    # die_ref_point(input_dex.all_part_rotation)
    # blank_holder_ref_point()
    # punch_ref_point()

    # # meshing
    # create_mesh(input_dex.solver_type, input_dex.mesh_size, input_dex.local_mesh_size)

    # create_blank_surfaces(input_dex.all_part_rotation)
    # create_die_surface(input_dex.all_part_rotation)
    # die_ref_point(input_dex.all_part_rotation)
    # create_blank_holder_surface()
    # blank_holder_ref_point()
    # punch_ref_point()
    # create_punch_surface(input_dex.all_part_rotation)
    # create_stopper_surface()


    # define surface interactions
    create_surface_interactions(input_dex.solver_type, input_dex.friction_coefficient)

    # loading and unloading phases
    create_boundary_conditions(input_dex.solver_type, input_dex.punch_velocity, input_dex.punch_depth, input_dex.mass_scalling, input_dex.trim_depth, input_dex.all_part_rotation, input_dex.BHF)

    # meshing
    create_mesh(input_dex.solver_type, input_dex.mesh_size, input_dex.local_mesh_size)

    apply_material_properties(input_dex.blank_material_name, input_dex.integration_points, input_dex.blank_thickness, input_dex.all_part_rotation)

    nCPU = 14 # vis suite

    # # remote desktop
    # nCPU = 4

    # # mech labs
    # nCPU = 8

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)
  
    run_sim_job(input_dex.solver_type,input_dex.simulation_output_path, nCPU)

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)

    spring_back_analysis(input_dex.simulation_output_path, nCPU)

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)

    # # demeri_test(input_dex.simulation_output_path, nCPU)

    # # cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    # # mdb.saveAs(pathName=cae_path)


def run_batch_sim(batch: list[InputeDex], main_sim_path: str):

    os.makedirs(main_sim_path, exist_ok=False) 

    results = []

    for i, sim in enumerate(batch, start=1):
        sim.simulation_output_path = os.path.join(main_sim_path, f"sim_{i}")

        os.makedirs(sim.simulation_output_path, exist_ok=False)

        run_sim(sim)
        plt.close('all')
        min_thickness = post_pro_thickness()
        plot_thickness_variation()
        post_pro_punch_reaction()
        post_pro_energy()
        plot_energy_data()
        plot_rf_data()
        energy_check()
        post_pro_strain()
        plot_strains()
        post_pro_strain_eq()
        post_pro_misses_stress()
        post_pro_spring_back_dev()

        # compare final stamped shape 
        export_final_node_positions(sim.punch_depth, sim.trim_depth, sim.die_profile_radius)
        height_diplacement = get_stamped_cup_height(sim.ideal_cup_height)
        rmse = compare_meshes()

        # ADD BACK FOR THE BO
        return min_thickness, rmse, height_diplacement

def compare_meshes():

    ideal_mesh = pv.read("ideal_part.obj")

    data = np.loadtxt("trim_results.csv", delimiter=",", skiprows=1)
    points = data[:, 1:4]
    stamped_mesh = pv.PolyData(points)

    ideal_centroid_y = np.mean(ideal_mesh.points[:, 1])
    stamped_centroid_y = np.mean(stamped_mesh.points[:, 1])
    y_translation = ideal_centroid_y - stamped_centroid_y
    stamped_mesh.translate((0, y_translation, 0), inplace=True)
    stamped_mesh_aligned = stamped_mesh

    # stamped_mesh_aligned = stamped_mesh.align(ideal_mesh)

    _, closest_points_to_ideal = ideal_mesh.find_closest_cell(stamped_mesh_aligned.points, return_closest_point=True)

    distances_to_ideal = np.linalg.norm(stamped_mesh_aligned.points - closest_points_to_ideal, axis=1)

    rmse = np.sqrt(np.mean(distances_to_ideal**2))

    print(f"RMSE: {rmse:.2f} mm")

    with open("compare_stamped_shape.txt", "w") as f:

        f.write("RMSE\n")
        f.write(f"{rmse}")

    # Plot the aligned meshes
    # plotter = pv.Plotter()
    # plotter.add_mesh(ideal_mesh, color="blue", opacity=0.5, label="Ideal Mesh")
    # plotter.add_mesh(stamped_mesh_aligned, color="red", opacity=0.5, label="Stamped Mesh (Aligned)")
    # plotter.add_legend()
    # plotter.show()

    # # Plot histogram of deviation distances
    # plt.figure(figsize=(8, 5))
    # plt.hist(distances, bins=50, color="blue", alpha=0.7)
    # plt.axvline(rmse, color="red", linestyle="dashed", linewidth=2, label=f"RMSE: {rmse:.2f} mm")
    # plt.axvline(max_deviation, color="black", linestyle="dashed", linewidth=2, label=f"Max Deviation: {max_deviation:.2f} mm")
    # plt.xlabel("Deviation Distance (mm)")
    # plt.ylabel("Frequency")
    # plt.title("Deviation of Stamped Mesh to Ideal Mesh")
    # plt.legend()
    # plt.grid()

    # # Show the plot
    # plt.show()
        
    # # Create PyVista plot with color-coded deviations
    # stamped_mesh_aligned["Deviation"] = distances  # Assign deviation as a scalar field

    # plotter = pv.Plotter()
    # plotter.add_mesh(stamped_mesh_aligned, scalars="Deviation", cmap="jet", point_size=5, render_points_as_spheres=True)
    # # plotter.add_mesh(ideal_mesh, color="white", opacity=0.3, label="Ideal Mesh")
    # plotter.add_scalar_bar(title="Deviation Distance (mm)")
    # plotter.show()

    return rmse

def post_pro_thickness():

    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_frame = odb.steps["load"].frames[-1]  # Last frame

    nodes_blank = odb.rootAssembly.instances["BLANK-1"].nodes

    sth_field = last_frame.fieldOutputs["STH"]
    displacement_field = last_frame.fieldOutputs["U"]

    min_thickness = 10

    output_file = os.path.join(cwd, "thickness_results.csv")
    with open(output_file, "w") as f:
        f.write("ElementID,SectionThickness,RadialDistance\n")  # CSV header
        
        for (value_sth, value_nodes) in zip(sth_field.values, nodes_blank):
            element_id = value_sth.elementLabel
            thickness = value_sth.data
            coord = value_nodes.coordinates

            if thickness < min_thickness:
                min_thickness = thickness

            X = coord[0]
            Y = coord[1]
            radial_distance =  sqrt(X**2 + Y**2)
            
            f.write(f"{element_id},{thickness},{radial_distance}\n")
    
    print("saved thickness data")

    odb.close()

    return min_thickness

def post_pro_punch_reaction():

    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    key = odb.steps["load"].historyRegions.keys()[-1]
    rf_data = odb.steps["load"].historyRegions[key].historyOutputs["RF2"].data
    rf_displacement = odb.steps["load"].historyRegions[key].historyOutputs["U2"].data

    output_file = os.path.join(cwd, "reaction_force_results.csv")
    with open(output_file, "w") as f:
        f.write("time,displacement,rf\n")  # CSV header
        
        for (value_rf, value_displacement) in zip(rf_data, rf_displacement):
            
            f.write(f"{value_rf[0]},{value_displacement[1]},{value_rf[1]}\n")
    
    print("saved reaction force data")

    odb.close()

def post_pro_energy():

    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    KE = odb.steps["load"].historyRegions["Assembly ASSEMBLY"].historyOutputs["ALLKE"]
    IE = odb.steps["load"].historyRegions["Assembly ASSEMBLY"].historyOutputs["ALLIE"]

    output_file = os.path.join(cwd, "energy_results.csv")
    with open(output_file, "w") as f:
        f.write("time,IE,KE\n")  # CSV header
        
        for (ke_value, ie_value) in zip(KE.data, IE.data):
            
            f.write(f"{ke_value[0]},{ie_value[1]},{ke_value[1]}\n")
    
    print("saved energy data")

    # Close ODB
    odb.close()

def plot_thickness_variation():

    cwd = os.getcwd()

    thickness_file_path = f"{cwd}/thickness_results.csv"

    print(thickness_file_path)

    thickness_df = pd.read_csv(thickness_file_path)

    thickness_df =thickness_df.sort_values(by="RadialDistance")

    plt.figure()

    plt.scatter(thickness_df["RadialDistance"], thickness_df["SectionThickness"], s=2)
    plt.xlabel("Radial Distance (mm)")
    plt.ylabel("Thickness (mm)")
    plt.grid(True)

    plt.savefig("thickness_plot.png")
 

def plot_rf_data():

    cwd = os.getcwd()

    rf_file_path = f"{cwd}/reaction_force_results.csv"

    rf_df = pd.read_csv(rf_file_path)

    rf_df["displacement"] = rf_df["displacement"] * -1
    rf_df["rf"] = rf_df["rf"] * -1e-3

    plt.figure()

    plt.plot(rf_df["displacement"], rf_df["rf"] * 4) # for quarter model
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Reaction Force (KN)")
    plt.grid(True)

    plt.savefig("rf_plot.png")

def plot_energy_data():

    cwd = os.getcwd()

    energy_file_path = f"{cwd}/energy_results.csv"

    energy_df = pd.read_csv(energy_file_path)

    energy_df["IE"] = energy_df["IE"] / 1e3
    energy_df["KE"] = energy_df["KE"] / 1e3


    plt.figure()
    plt.plot(energy_df["time"], energy_df["IE"])
    plt.plot(energy_df["time"], energy_df["KE"])
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (KJ)")
    plt.legend(["IE", "KE"])
    plt.grid(True)
    
    plt.savefig("energy_plot.png")

def energy_check():

    cwd = os.getcwd()

    energy_file_path = f"{cwd}/energy_results.csv"

    energy_df = pd.read_csv(energy_file_path)

    max_ke = max(energy_df["KE"]) 
    max_ie = max(energy_df["IE"])

    ke_ie_pct = max_ke * 100 / max_ie

    with open("energy_check.txt", "w") as f:

        f.write(f"ke/ie ratio = {round(ke_ie_pct, 4)}%")

def post_pro_strain():

    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_frame = odb.steps["load"].frames[-1]  # Last frame

    true_strain_field = last_frame.fieldOutputs["LE"]

    output_file = os.path.join(cwd, "strain_results.csv")
    with open(output_file, "w") as f:
        f.write("ElementID,LE11,LE22\n")  # CSV header
        
        for value_strain in true_strain_field.values:
            element_id = value_strain.elementLabel
            strain_data = value_strain.data

            f.write(f"{element_id},{strain_data[0]},{strain_data[1]}\n")
    
    print("saved strain data")

    # Close ODB
    odb.close()

def plot_strains():

    cwd = os.getcwd()

    strain_file_path = f"{cwd}/strain_results.csv"

    strain_df = pd.read_csv(strain_file_path)

    plt.figure()
    plt.scatter(strain_df["LE22"], strain_df["LE11"], s=2)
    plt.xlabel("minor strain (LE22)")
    plt.ylabel("major strain (LE11)")
    plt.grid(True)
    
    plt.savefig("strain_plot.png")

def post_pro_strain_eq():

    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_frame = odb.steps["load"].frames[-1]  # Last frame

    strain_field = last_frame.fieldOutputs["PEEQ"]

    output_file = os.path.join(cwd, "strain_eq_results.csv")
    with open(output_file, "w") as f:
        f.write("ElementID,PEEQ\n")  # CSV header
        
        for strain_value in strain_field.values:
            element_id = strain_value.elementLabel
            peeq = strain_value.data

            
            f.write(f"{element_id},{peeq}\n")
    
    print("saved peeq data")

    # Close ODB
    odb.close()

def post_pro_misses_stress():


    cwd = os.getcwd()

    odb_filename = "stamping_sim.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_frame = odb.steps["load"].frames[-1]  # Last frame

    stress_field = last_frame.fieldOutputs['S']

    mises_field = stress_field.getScalarField(invariant=MISES)

    output_file = os.path.join(cwd, "mises_stress_results.csv")
    with open(output_file, "w") as f:
        f.write("ElementID,MISES\n")
        
        for mises_value in mises_field.values:
            element_id = mises_value.elementLabel
            mises = mises_value.data

            
            f.write(f"{element_id},{mises}\n")
    
    print("saved mises stress data")

    odb.close()

def post_pro_spring_back_dev():

    cwd = os.getcwd()

    odb_filename = "spring_back.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_frame = odb.steps["Step-1"].frames[-1]  # Last frame

    displacement_field = last_frame.fieldOutputs['U']

    displacement_mag_field = displacement_field.getScalarField(invariant=MAGNITUDE)


    output_file = os.path.join(cwd, "springback_displacement_results.csv")
    with open(output_file, "w") as f:
        f.write("ElementID,U_mag,X,Y,Z\n")
        
        for (dis_value, mag_value) in zip(displacement_field.values, displacement_mag_field.values ):
            element_id = dis_value.elementLabel
            X = dis_value.dataDouble[0]
            Y = dis_value.dataDouble[1]
            Z = dis_value.dataDouble[2]
            mag = mag_value.data

            
            f.write(f"{element_id},{mag},{X},{Y},{Z}\n")
    
    print("saved springback displacement data")

    odb.close()


def export_final_node_positions(punch_depth, trim_depth, die_profile_radius):

    cwd = os.getcwd()

    odb_filename = "spring_back.odb"
    odb_path = os.path.join(cwd, odb_filename)

    odb = openOdb(odb_path)

    last_step = odb.steps["Step-1"]

    last_frame = last_step.frames[-1]

    disp_field = last_frame.fieldOutputs['U']

    nodes_blank = odb.rootAssembly.instances["BLANK-1"].nodes

    output_file = os.path.join(cwd, "all_stamped_nodes.csv")
    with open(output_file, "w") as f:
        f.write("NodeLabel,X,Y,Z\n")

        for (value, node_value) in zip(disp_field.values, nodes_blank):
            node_label = value.nodeLabel
            x_displacement = value.dataDouble[0]  # U1
            y_displacement = value.dataDouble[1]  # U2
            z_displacement = value.dataDouble[2]  # U3

            x_initial, y_initial, z_initial = node_value.coordinates

            x_new = x_initial + x_displacement
            y_new = y_initial + y_displacement
            z_new = z_initial + z_displacement


            f.write(f"{node_label},{x_new},{y_new},{z_new}\n")

    trim_length = punch_depth - trim_depth - die_profile_radius
    
    nodes = pd.read_csv("all_stamped_nodes.csv")

    max_height  = nodes["Y"].max()

    filtered_nodes = nodes[nodes['Y'] <=  max_height - trim_length]

    nodes.to_csv("trim_results.csv", index=False) # change back to filtered / check filtered ???????????????????????
    odb.close()

    print("saved trimmed mesh")

def get_stamped_cup_height(ideal_cup_height):

    cwd = os.getcwd()

    nodes = pd.read_csv("trim_results.csv")

    cup_height = nodes["Y"].max() - nodes["Y"].min() 

    stamped_ideal_cup_displacement = cup_height - ideal_cup_height


    output_file_path = os.path.join(cwd, "stamped_ideal_cup_displacement.txt")

    with open(output_file_path, "w") as f:

        f.write(str(stamped_ideal_cup_displacement))

    return stamped_ideal_cup_displacement



def mesh_conv_final_max():

    import pandas as pd

    sim_numbers = [1, 2, 3, 4, 5, 6, 7, 8]
    data_list = ["thickness_results.csv","strain_eq_results.csv", "springback_displacement_results.csv"]
    var_list = ["SectionThickness","PEEQ", "U_mag"]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp263"

    for sim in sim_numbers:

        print(sim)

        sim_path = rf"{main_sim_output}\sim_{sim}"
        # print(sim_path)

        for (data, var) in zip(data_list, var_list):

            data_path = rf"{sim_path}\{data}"

            # print(data_path)

            data_temp = pd.read_csv(data_path)

            max_val = data_temp[var].max()

            print(var, max_val)

def mesh_conv_final_min():

    import pandas as pd

    sim_numbers = [8]
    data_list = ["thickness_results.csv", "reaction_force_results.csv"]
    var_list = ["SectionThickness", "rf"]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp263"

    for sim in sim_numbers:

        print(sim)

        sim_path = rf"{main_sim_output}\sim_{sim}"
        # print(sim_path)

        for (data, var) in zip(data_list, var_list):

            data_path = rf"{sim_path}\{data}"

            # print(data_path)

            data_temp = pd.read_csv(data_path)

            max_val = data_temp[var].min()

            print(var, max_val)


def plot_exp():

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    # Load data from file (adjust the filename as needed)
    data = np.loadtxt('punch_rf_experimental.txt', delimiter=',')

    # Split data into x and y arrays
    x = data[:, 0]
    y = data[:, 1]

    # Create the plot
    plt.plot(x, y, marker='o', linestyle='-', color='b', label='Experimental Data')


    cwd = r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp263\sim_5"

    rf_file_path = f"{cwd}/reaction_force_results.csv"

    rf_df = pd.read_csv(rf_file_path)
    
    rf_df["displacement"] = rf_df["displacement"] * -1
    rf_df["rf"] = rf_df["rf"] * -1e-3

    plt.plot(rf_df["displacement"], rf_df["rf"] * 4, color="r", label="Simulation Data") # for quarter model
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Reaction Force (KN)")
    plt.legend()
    plt.grid(True)

    plt.savefig("rf_plot_exp_sim.png")

def cylindrical_cup_function(ideal_cup_height, ideal_cup_radius, ideal_cup_profile_radius, punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap, blank_radius):


    input_dex = InputeDex()

    input_dex.solver_type = SolverType.STANDARD
    
    input_dex.all_part_rotation = 90
    input_dex.integration_points = 5

    # ideal part dimensions
    input_dex.ideal_cup_radius = ideal_cup_radius
    input_dex.ideal_cup_height = ideal_cup_height
    input_dex.ideal_cup_profile_radius = ideal_cup_profile_radius

    # blank inputs
    input_dex.blank_radius = blank_radius
    input_dex.blank_thickness = 1

    # punch inputs
    input_dex.punch_profile_radius = punch_profile_radius
    input_dex.punch_min_radius = punch_min_radius

    # die inputs
    input_dex.die_profile_radius = die_profile_radius
    input_dex.die_min_radius = punch_min_radius + die_punch_gap
    input_dex.die_max_radius = input_dex.blank_radius + 15

    # blank holder inputs
    input_dex.blank_holder_height = 10
    input_dex.blank_holder_profile_radius = 0.1
    input_dex.blank_holder_min_radius = input_dex.die_min_radius
    input_dex.blank_holder_max_radius = input_dex.blank_radius + 10
    input_dex.blank_holder_die_gap = input_dex.blank_thickness

    # punch
    input_dex.punch_depth = punch_depth
    input_dex.die_height = punch_depth + 5

    # material inputs
    input_dex.blank_material_name = "AA5754-O"
    input_dex.density = 2.7e-06 
    input_dex.youngs_modulus = 68000.0
    input_dex.posissons_ratio = 0.33
    input_dex.plastic_material_data = (
        (91.74, 0.0), (149.8284957, 0.02020202), (184.0636139, 0.04040404), (
        208.8813246, 0.060606061), (227.7563027, 0.080808081), (242.4816409, 
        0.101010101), (254.160565, 0.121212121), (263.5336322, 0.141414141), (
        271.1246064, 0.161616162), (277.3170581, 0.181818182), (282.3989825, 
        0.202020202), (286.590712, 0.222222222), (290.0633061, 0.242424242), (
        292.9511489, 0.262626263), (295.3608473, 0.282828283), (297.3776686, 
        0.303030303), (299.0702914, 0.323232323), (300.4943701, 0.343434343), (
        301.6952452, 0.363636364), (302.7100293, 0.383838384), (303.5692279, 
        0.404040404), (304.2980115, 0.424242424), (304.9172192, 0.444444444), (
        305.4441588, 0.464646465), (305.8932459, 0.484848485), (306.2765196, 
        0.505050505), (306.604058, 0.525252525), (306.8843178, 0.545454545), (
        307.1244092, 0.565656566), (307.3303228, 0.585858586), (307.5071151, 
        0.606060606), (307.6590614, 0.626262626), (307.7897827, 0.646464646), (
        307.9023503, 0.666666667), (307.9993734, 0.686868687), (308.0830718, 
        0.707070707), (308.1553359, 0.727272727), (308.2177785, 0.747474747), (
        308.2717764, 0.767676768), (308.318507, 0.787878788), (308.3589778, 
        0.808080808), (308.394052, 0.828282828), (308.4244699, 0.848484848), (
        308.4508671, 0.868686869), (308.4737898, 0.888888889), (308.4937077, 
        0.909090909), (308.5110251, 0.929292929), (308.5260903, 0.949494949), (
        308.5392038, 0.96969697), (308.5506247, 0.98989899), (308.560577, 
        1.01010101), (308.5692539, 1.03030303), (308.5768228, 1.050505051), (
        308.5834286, 1.070707071), (308.5891965, 1.090909091), (308.5942353, 
        1.111111111), (308.5986392, 1.131313131), (308.6024899, 1.151515152), (
        308.6058585, 1.171717172), (308.6088065, 1.191919192), (308.6113876, 
        1.212121212), (308.6136483, 1.232323232), (308.6156294, 1.252525253), (
        308.6173659, 1.272727273), (308.6188888, 1.292929293), (308.6202248, 
        1.313131313), (308.6213974, 1.333333333), (308.6224267, 1.353535354), (
        308.6233308, 1.373737374), (308.6241251, 1.393939394), (308.6248232, 
        1.414141414), (308.625437, 1.434343434), (308.6259767, 1.454545455), (
        308.6264516, 1.474747475), (308.6268695, 1.494949495), (308.6272374, 
        1.515151515), (308.6275614, 1.535353535), (308.6278468, 1.555555556), (
        308.6280983, 1.575757576), (308.6283199, 1.595959596), (308.6285154, 
        1.616161616), (308.6286877, 1.636363636), (308.6288397, 1.656565657), (
        308.6289739, 1.676767677), (308.6290923, 1.696969697), (308.6291969, 
        1.717171717), (308.6292892, 1.737373737), (308.6293708, 1.757575758), (
        308.6294429, 1.777777778), (308.6295066, 1.797979798), (308.6295629, 
        1.818181818), (308.6296127, 1.838383838), (308.6296568, 1.858585859), (
        308.6296958, 1.878787879), (308.6297302, 1.898989899), (308.6297608, 
        1.919191919), (308.6297878, 1.939393939), (308.6298117, 1.95959596), (
        308.6298329, 1.97979798), (308.6298517, 2.0)
        )
    
    # BCs
    input_dex.punch_velocity = 1.1
    input_dex.friction_coefficient = 0.09
    input_dex.BHF = 6000 / 4

    # mesh
    input_dex.mesh_size = 5
    input_dex.local_mesh_size = 5 * 0.75 # CHNAGE BACK PLEASE

    input_dex_set = [input_dex]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3\sim_runs/1"

    def get_sim_out_path(main_sim_output):
        if os.path.isdir(main_sim_output):
            sim_number = int(main_sim_output[-1])
            sim_number += 1
            main_sim_output = main_sim_output[:-1] + str(sim_number)
            return get_sim_out_path(main_sim_output)
        
        return main_sim_output
    
    main_sim_output = get_sim_out_path(main_sim_output)

    min_thickness, rmse, height_difference = run_batch_sim(input_dex_set, main_sim_output)

    with open(fr"{main_sim_output}\sim_1\sim_status.txt", "r") as f:
        sim_status = f.read().strip()

    if sim_status == "0":
        height_difference = None

    return min_thickness, rmse, height_difference


def run_cup():

    # Read command-line arguments passed from subprocess
    args = sys.argv[8:]
    
    # Convert input parameters
    ideal_height = float(args[0])
    ideal_radius = float(args[1])
    ideal_profile_radius = float(args[2])
    punch_depth = float(args[3])
    punch_profile_radius = float(args[4])
    punch_min_radius = float(args[5])
    die_profile_radius = float(args[6])
    die_punch_gap = float(args[7])
    blank_radius = float(args[8])

    # Run simulation
    _, _, height_difference = cylindrical_cup_function(
        ideal_height, ideal_radius, ideal_profile_radius, 
        punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap, blank_radius
    )

    if height_difference == None:
        data_out = {"height difference": None}
    else:
        data_out = {"height difference": -abs(height_difference)}
    
    with open(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\stamping_op_build3\sim_output.json", "w") as f:
        json.dump(data_out, f)

def main():

    input_dex = InputeDex()

    input_dex.solver_type = SolverType.STANDARD
    
    input_dex.all_part_rotation = 90
    input_dex.integration_points = 5

    input_dex.trim_depth = 0

    # ideal part dimensions
    input_dex.ideal_cup_radius = 33/2
    input_dex.ideal_cup_height = 20
    input_dex.ideal_cup_profile_radius = 5

    # blank inputs
    input_dex.blank_radius = 30
    input_dex.blank_thickness = 1
    input_dex.integration_points = 15

    # punch inputs
    input_dex.punch_profile_radius = 5 # 5
    input_dex.punch_min_radius = 33/2

    # die inputs
    input_dex.die_profile_radius = 5
    input_dex.die_min_radius = 35.3/2
    input_dex.die_max_radius = 40

    # blank holder inputs
    input_dex.blank_holder_height = 10
    input_dex.blank_holder_profile_radius = 0.1
    input_dex.blank_holder_min_radius = 33.6/2
    input_dex.blank_holder_max_radius = 70/2
    input_dex.blank_holder_die_gap = input_dex.blank_thickness

    # punch
    input_dex.cup_design_height = 10
    input_dex.punch_depth = 35
    input_dex.die_height = 40

    # material inputs
    input_dex.blank_material_name = "AA5754-O"
    input_dex.density = 2.7e-06 
    input_dex.youngs_modulus = 68000.0
    input_dex.posissons_ratio = 0.33
    input_dex.plastic_material_data = (
        (91.74, 0.0), (149.8284957, 0.02020202), (184.0636139, 0.04040404), (
        208.8813246, 0.060606061), (227.7563027, 0.080808081), (242.4816409, 
        0.101010101), (254.160565, 0.121212121), (263.5336322, 0.141414141), (
        271.1246064, 0.161616162), (277.3170581, 0.181818182), (282.3989825, 
        0.202020202), (286.590712, 0.222222222), (290.0633061, 0.242424242), (
        292.9511489, 0.262626263), (295.3608473, 0.282828283), (297.3776686, 
        0.303030303), (299.0702914, 0.323232323), (300.4943701, 0.343434343), (
        301.6952452, 0.363636364), (302.7100293, 0.383838384), (303.5692279, 
        0.404040404), (304.2980115, 0.424242424), (304.9172192, 0.444444444), (
        305.4441588, 0.464646465), (305.8932459, 0.484848485), (306.2765196, 
        0.505050505), (306.604058, 0.525252525), (306.8843178, 0.545454545), (
        307.1244092, 0.565656566), (307.3303228, 0.585858586), (307.5071151, 
        0.606060606), (307.6590614, 0.626262626), (307.7897827, 0.646464646), (
        307.9023503, 0.666666667), (307.9993734, 0.686868687), (308.0830718, 
        0.707070707), (308.1553359, 0.727272727), (308.2177785, 0.747474747), (
        308.2717764, 0.767676768), (308.318507, 0.787878788), (308.3589778, 
        0.808080808), (308.394052, 0.828282828), (308.4244699, 0.848484848), (
        308.4508671, 0.868686869), (308.4737898, 0.888888889), (308.4937077, 
        0.909090909), (308.5110251, 0.929292929), (308.5260903, 0.949494949), (
        308.5392038, 0.96969697), (308.5506247, 0.98989899), (308.560577, 
        1.01010101), (308.5692539, 1.03030303), (308.5768228, 1.050505051), (
        308.5834286, 1.070707071), (308.5891965, 1.090909091), (308.5942353, 
        1.111111111), (308.5986392, 1.131313131), (308.6024899, 1.151515152), (
        308.6058585, 1.171717172), (308.6088065, 1.191919192), (308.6113876, 
        1.212121212), (308.6136483, 1.232323232), (308.6156294, 1.252525253), (
        308.6173659, 1.272727273), (308.6188888, 1.292929293), (308.6202248, 
        1.313131313), (308.6213974, 1.333333333), (308.6224267, 1.353535354), (
        308.6233308, 1.373737374), (308.6241251, 1.393939394), (308.6248232, 
        1.414141414), (308.625437, 1.434343434), (308.6259767, 1.454545455), (
        308.6264516, 1.474747475), (308.6268695, 1.494949495), (308.6272374, 
        1.515151515), (308.6275614, 1.535353535), (308.6278468, 1.555555556), (
        308.6280983, 1.575757576), (308.6283199, 1.595959596), (308.6285154, 
        1.616161616), (308.6286877, 1.636363636), (308.6288397, 1.656565657), (
        308.6289739, 1.676767677), (308.6290923, 1.696969697), (308.6291969, 
        1.717171717), (308.6292892, 1.737373737), (308.6293708, 1.757575758), (
        308.6294429, 1.777777778), (308.6295066, 1.797979798), (308.6295629, 
        1.818181818), (308.6296127, 1.838383838), (308.6296568, 1.858585859), (
        308.6296958, 1.878787879), (308.6297302, 1.898989899), (308.6297608, 
        1.919191919), (308.6297878, 1.939393939), (308.6298117, 1.95959596), (
        308.6298329, 1.97979798), (308.6298517, 2.0)
        )
    
    # BCs
    input_dex.punch_velocity = 1.1
    input_dex.mass_scalling = 5e-6

    input_dex.friction_coefficient = 0.09
    input_dex.BHF = 6000 / 4

    # mesh
    input_dex.mesh_size = 5
    input_dex.local_mesh_size = 5 * 0.75

    # input dex 2
    input_dex_2: InputeDex = copy.deepcopy(input_dex)
    input_dex_2.mesh_size = 4
    input_dex_2.local_mesh_size = 4 * 0.75

    # input dex 3
    input_dex_3: InputeDex = copy.deepcopy(input_dex)
    input_dex_3.mesh_size = 3
    input_dex_3.local_mesh_size = 3 * 0.75

    # input dex 4
    input_dex_4: InputeDex = copy.deepcopy(input_dex)
    input_dex_4.mesh_size = 2
    input_dex_4.local_mesh_size = 2 * 0.75

    # input dex 5
    input_dex_5: InputeDex = copy.deepcopy(input_dex)
    input_dex_5.mesh_size = 1
    input_dex_5.local_mesh_size = 1 * 0.75

    # input dex 6
    input_dex_6: InputeDex = copy.deepcopy(input_dex)
    input_dex_6.mesh_size = 0.75
    input_dex_6.local_mesh_size = 0.75 * 0.75

    # input dex 7
    input_dex_7: InputeDex = copy.deepcopy(input_dex)
    input_dex_7.mesh_size = 0.625
    input_dex_7.local_mesh_size = 0.625 * 0.75

    # input dex 8
    input_dex_8: InputeDex = copy.deepcopy(input_dex)
    input_dex_8.mesh_size = 0.5
    input_dex_8.local_mesh_size = 0.5 * 0.75

    # input dex 9
    input_dex_9: InputeDex = copy.deepcopy(input_dex)
    input_dex_9.mesh_size = 5
    input_dex_9.local_mesh_size = 5 * 0.5

    # input dex 10
    input_dex_10: InputeDex = copy.deepcopy(input_dex)
    input_dex_10.mesh_size = 4
    input_dex_10.local_mesh_size = 4 * 0.5

    # input dex 11
    input_dex_11: InputeDex = copy.deepcopy(input_dex)
    input_dex_11.mesh_size = 3
    input_dex_11.local_mesh_size = 3 * 0.5

    # input dex 12
    input_dex_12: InputeDex = copy.deepcopy(input_dex)
    input_dex_12.mesh_size = 2
    input_dex_12.local_mesh_size = 2 * 0.5

    # input dex 13
    input_dex_13: InputeDex = copy.deepcopy(input_dex)
    input_dex_13.mesh_size = 5
    input_dex_13.local_mesh_size = 5 * 0.333

    # input dex 14
    input_dex_14: InputeDex = copy.deepcopy(input_dex)
    input_dex_14.mesh_size = 4
    input_dex_14.local_mesh_size = 4 * 0.333

    # input dex 15
    input_dex_15: InputeDex = copy.deepcopy(input_dex)
    input_dex_15.mesh_size = 3
    input_dex_15.local_mesh_size = 3 * 0.333

    # input dex 16
    input_dex_16: InputeDex = copy.deepcopy(input_dex)
    input_dex_16.mesh_size = 2
    input_dex_16.local_mesh_size = 2 * 0.333

    # input dex 12
    input_dex_17: InputeDex = copy.deepcopy(input_dex)
    input_dex_17.mesh_size = 1
    input_dex_17.local_mesh_size = 1 * 0.5

    # input dex 12
    input_dex_18: InputeDex = copy.deepcopy(input_dex)
    input_dex_18.mesh_size = 1
    input_dex_18.local_mesh_size = 1 * 0.33

    input_dex_set = [input_dex]
    # input_dex_set = [input_dex, input_dex_2, input_dex_3]
    # input_dex_set = [input_dex, input_dex_2, input_dex_3, input_dex_4, input_dex_5, input_dex_6, input_dex_7, input_dex_8, input_dex_9, input_dex_10, input_dex_11, input_dex_12,
    #                  input_dex_13, input_dex_14, input_dex_15, input_dex_16, input_dex_17, input_dex_18]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp317"
    # main_sim_output =  r"C:\Users\Kadmiel McForrester\OneDrive - University of Bath\Documents\build\batch_temp92"

    run_batch_sim(input_dex_set, main_sim_output)


if __name__ == "__main__":
    
    run_cup()

    # main()


    # cup height post pro
    # sims = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    # sims = [2] # 294
    # sims = [1] # 295
    # for sim in sims:
    #     os.chdir(rf"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp295\sim_{sim}")
    #     compare_meshes()


    # notes
# add retun back to run batch for optimisation loop
# fix the non completed simulations