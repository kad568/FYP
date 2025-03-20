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
    mass_scalling:float = None

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
        sketchPlaneSide=SIDE1, origin=(16.976527, 0.55, 16.976527))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=113.14, gridSpacing=2.82, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['Model-1'].parts['blank']
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.CircleByCenterPerimeter(center=(16.976527, -16.976527), point1=(8.46, 
        -7.755)) # change inner circle radius
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    e, d1 = p.edges, p.datums
    p.PartitionFaceBySketch(sketchUpEdge=e[1], faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']


def create_blank_surfaces():

    p = mdb.models['Model-1'].parts['blank']
    s = p.faces
    side2Faces = s.getSequenceFromMask(mask=('[#3 ]', ), )
    p.Surface(side2Faces=side2Faces, name='blank_top_surface')
    p = mdb.models['Model-1'].parts['blank']
    s = p.faces
    side1Faces = s.getSequenceFromMask(mask=('[#3 ]', ), )
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


def create_die_surface():

    p1 = mdb.models['Model-1'].parts["die"]
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts["die"]
    s = p.faces
    side2Faces = s.getSequenceFromMask(mask=('[#7 ]', ), )
    p.Surface(side2Faces=side2Faces, name="die_surface")


def die_ref_point():

    p = mdb.models['Model-1'].parts["die"]
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
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


def create_punch_surface():

    p = mdb.models['Model-1'].parts["punch"]
    s = p.faces
    side2Faces = s.getSequenceFromMask(mask=('[#7 ]', ), )
    p.Surface(side2Faces=side2Faces, name="punch_surface")


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


def create_blank_material(blank_material_name, density, youngs_modulus, posissons_ratio, plastic_material_data):

    mdb.models['Model-1'].Material(name=blank_material_name)
    mdb.models['Model-1'].materials[blank_material_name].Density(table=((density, ), 
        ))
    mdb.models['Model-1'].materials[blank_material_name].Elastic(table=((youngs_modulus, 
        posissons_ratio), ))
    mdb.models['Model-1'].materials[blank_material_name].Plastic(scaleStress=None, 
        table=plastic_material_data)


def create_surface_interactions(solver_type: SolverType, friction_coefficient):



    # create interaction properties
    mdb.models['Model-1'].ContactProperty('die_blank')
    mdb.models['Model-1'].ContactProperty('punch_blank')
    mdb.models['Model-1'].ContactProperty('blank_holder_blank')

    mdb.models['Model-1'].interactionProperties['blank_holder_blank'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((friction_coefficient, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    
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
        
def create_boundary_conditions(solver_type, punch_speed, punch_depth, mass_scaling, trim_buffer):

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
            maxNumInc=1000,
            initialInc=0.1, 
            maxInc=1.0, 
            nlgeom=ON
        )

        mdb.models['Model-1'].steps['load'].Restart(
            frequency=1,
            overlay=OFF
        )


    # session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='load')
    # session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        # predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#30 ]', ), )
    region = a.Set(edges=edges1, name='Set-1')
    mdb.models['Model-1'].XsymmBC(name='BC-1', createStepName='load', 
        region=region, localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['blank-1'].edges

    edges1 = e1.getSequenceFromMask(mask=('[#44 ]', ), )
    region = a.Set(edges=edges1, name='Set-2')
    mdb.models['Model-1'].ZsymmBC(name='BC-2', createStepName='load', 
        region=region, localCsys=None)
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
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    edges = p.edges.getByBoundingBox(-5, -5, -10, 5, 5, 10) 
    p.seedEdgeBySize(edges=edges, size=mesh_size / 3, constraint=FIXED)



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

    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1_shell, elemType2_shell))
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, technique=SWEEP)
    p = mdb.models['Model-1'].parts['blank']
    p.seedPart(size=local_seeding_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

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

    
def apply_material_properties(material_name, integration_points, blank_thickness):

    p = mdb.models['Model-1'].parts['blank']
    mdb.models['Model-1'].HomogeneousShellSection(name='blank_section', 
        preIntegrate=OFF, material=material_name, thicknessType=UNIFORM, 
        thickness=blank_thickness, thicknessField='', nodalThicknessField='', 
        idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT, 
        thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
        integrationRule=SIMPSON, numIntPts=integration_points)
    p = mdb.models['Model-1'].parts['blank']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#3 ]', ), )
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
    faces = f.getSequenceFromMask(mask=('[#3 ]', ), )
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

def spring_back_analysis(sim_out_path, nCPU):

    mdb.Model(name='Model-1-spring_back', objectToCopy=mdb.models['Model-1'])
    p = mdb.models['Model-1-spring_back'].parts['blank']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1-spring_back'].parts['blank_holder']
    del mdb.models['Model-1-spring_back'].parts['die']
    del mdb.models['Model-1-spring_back'].parts['punch']
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
    a.deleteFeatures(('blank_holder-1', 'die-1', 'punch-1', ))
    del mdb.models['Model-1-spring_back'].interactions['blank_holder_blank']
    del mdb.models['Model-1-spring_back'].interactions['die_blank']
    del mdb.models['Model-1-spring_back'].interactions['punch_blank']
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=OFF, 
        constraints=OFF, connectors=OFF, engineeringFeatures=OFF, 
        adaptiveMeshConstraints=ON)
    del mdb.models['Model-1-spring_back'].steps['load']
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
    mdb.models['Model-1-spring_back'].StaticStep(name='Step-1', previous='Initial', 
        maxNumInc=10000, initialInc=0.001, maxInc=0.1, nlgeom=ON)
    a = mdb.models['Model-1-spring_back'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#30 ]', ), )
    region = a.Set(edges=edges1, name='Set-7')
    a = mdb.models['Model-1-spring_back'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#30 ]', ), )
    region = a.Set(edges=edges1, name='Set-9')
    mdb.models['Model-1-spring_back'].XsymmBC(name='BC-1', createStepName='Step-1', 
        region=region, localCsys=None)
    a = mdb.models['Model-1-spring_back'].rootAssembly
    e1 = a.instances['blank-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#44 ]', ), )
    region = a.Set(edges=edges1, name='Set-10')
    mdb.models['Model-1-spring_back'].ZsymmBC(name='BC-2', createStepName='Step-1', 
        region=region, localCsys=None)

    a = mdb.models['Model-1-spring_back'].rootAssembly
    v1 = a.instances['blank-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#20 ]', ), )
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
    create_blank_surfaces()

    # create die
    create_die_part(input_dex.die_height, input_dex.die_profile_radius, input_dex.die_min_radius, input_dex.die_max_radius, input_dex.all_part_rotation)
    create_die_surface()
    die_ref_point()  

    # create bkank holder
    create_blank_holder(input_dex.blank_holder_height, input_dex.blank_holder_profile_radius, input_dex.blank_holder_min_radius, input_dex.blank_holder_max_radius, input_dex.blank_holder_die_gap, input_dex.all_part_rotation)
    create_blank_holder_surface()
    blank_holder_ref_point()

    # create punch
    create_punch(input_dex.punch_depth, input_dex.punch_profile_radius, input_dex.punch_min_radius, input_dex.blank_thickness, input_dex.all_part_rotation)
    punch_ref_point()
    create_punch_surface()

    # create ideal part
    create_ideal_part(input_dex.ideal_cup_height, input_dex.ideal_cup_profile_radius, input_dex.ideal_cup_radius, input_dex.all_part_rotation)

    # assemble parts
    create_part_assembly()

    # define material
    create_blank_material(input_dex.blank_material_name, input_dex.density, input_dex.youngs_modulus, input_dex.posissons_ratio, input_dex.plastic_material_data)

    # define surface interactions
    create_surface_interactions(input_dex.solver_type, input_dex.friction_coefficient)

    # loading and unloading phases
    create_boundary_conditions(input_dex.solver_type, input_dex.punch_velocity, input_dex.punch_depth, input_dex.mass_scalling, input_dex.trim_depth)

    # meshing
    create_mesh(input_dex.solver_type, input_dex.mesh_size, input_dex.local_mesh_size)

    apply_material_properties(input_dex.blank_material_name, input_dex.integration_points, input_dex.blank_thickness)

    nCPU = 14 # vis suite
    # nCPU = 4

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)
  
    run_sim_job(input_dex.solver_type,input_dex.simulation_output_path, nCPU)

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)

    spring_back_analysis(input_dex.simulation_output_path, nCPU)

    cae_path = f"{input_dex.simulation_output_path}/stamping_sim.cae"
    mdb.saveAs(pathName=cae_path)


def run_batch_sim(batch: list[InputeDex], main_sim_path: str):

    os.makedirs(main_sim_path, exist_ok=False) 

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
        rmse, max_deviation = compare_meshes()

        return min_thickness, rmse, max_deviation

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

    _, closest_points = ideal_mesh.find_closest_cell(stamped_mesh_aligned.points, return_closest_point=True)

    distances = np.linalg.norm(stamped_mesh_aligned.points - closest_points, axis=1)

    rmse = np.sqrt(np.mean(distances**2))
    max_deviation = np.max(distances)

    print(f"RMSE: {rmse:.2f} mm")
    print(f"Max Deviation: {max_deviation:.2f} mm")

    with open("compare_stamped_shape.txt", "w") as f:

        f.write("RMSE,max deviation\n")
        f.write(f"{rmse},{max_deviation}")

    # # Plot the aligned meshes
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

    return rmse, max_deviation

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

    plt.plot(rf_df["displacement"], rf_df["rf"])
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

    filtered_nodes.to_csv("trim_results.csv", index=False)
    odb.close()

    print("saved trimmed mesh")

def cylindrical_cup_function(ideal_cup_height, ideal_cup_radius, ideal_cup_profile_radius, punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap):


    input_dex = InputeDex()

    input_dex.solver_type = SolverType.STANDARD
    
    input_dex.all_part_rotation = 90

    # ideal part dimensions
    input_dex.ideal_cup_radius = ideal_cup_radius
    input_dex.ideal_cup_height = ideal_cup_height
    input_dex.ideal_cup_profile_radius = ideal_cup_profile_radius

    # blank inputs
    input_dex.blank_radius = 40
    input_dex.blank_thickness = 1.1
    input_dex.integration_points = 15

    # punch inputs
    input_dex.punch_profile_radius = punch_profile_radius
    input_dex.punch_min_radius = punch_min_radius

    # die inputs
    input_dex.die_profile_radius = die_profile_radius
    input_dex.die_min_radius = punch_min_radius + die_punch_gap
    input_dex.die_max_radius = input_dex.blank_radius + 5

    # blank holder inputs
    input_dex.blank_holder_height = 10
    input_dex.blank_holder_profile_radius = 6.5
    input_dex.blank_holder_min_radius = input_dex.die_min_radius
    input_dex.blank_holder_max_radius = input_dex.blank_radius + 5
    input_dex.blank_holder_die_gap = input_dex.blank_thickness + 0.5

    # punch
    input_dex.punch_depth = punch_depth
    input_dex.die_height = punch_depth + 5

    # material inputs
    input_dex.blank_material_name = "AA1050"
    input_dex.density = 2.7e-06 
    input_dex.youngs_modulus = 70000.0
    input_dex.posissons_ratio = 0.33
    input_dex.plastic_material_data = ((35.0, 0.0), (66.56588883, 0.01), (73.19453891, 0.02), (
        77.6998518, 0.03), (81.21516632, 0.04), (84.1399572, 0.05), (
        86.66656829, 0.06), (88.90387465, 0.07), (90.92007807, 0.08), (
        92.76100124, 0.09), (94.45905775, 0.1), (96.03810066, 0.11), (
        97.51624224, 0.12), (98.90759062, 0.13), (100.2233697, 0.14), (
        101.4726693, 0.15), (102.6629639, 0.16), (103.8004812, 0.17), (
        104.8904701, 0.18), (105.9373994, 0.19), (106.9451084, 0.2), (
        107.9169222, 0.21), (108.855741, 0.22), (109.7641117, 0.23), (
        110.6442836, 0.24), (111.4982544, 0.25), (112.3278069, 0.26), (
        113.1345397, 0.27), (113.9198918, 0.28), (114.6851642, 0.29), (
        115.431537, 0.3), (116.1600841, 0.31), (116.8717864, 0.32), (
        117.567542, 0.33), (118.2481756, 0.34), (118.9144466, 0.35), (
        119.5670557, 0.36), (120.2066512, 0.37), (120.833834, 0.38), (
        121.4491624, 0.39), (122.053156, 0.4), (122.6462992, 0.41), (
        123.2290448, 0.42), (123.8018162, 0.43), (124.3650101, 0.44), (
        124.9189991, 0.45), (125.4641332, 0.46), (126.000742, 0.47), (
        126.529136, 0.48), (127.0496083, 0.49), (127.5624356, 0.5), (128.06788, 
        0.51), (128.5661893, 0.52), (129.0575987, 0.53), (129.5423312, 0.54), (
        130.0205987, 0.55), (130.4926026, 0.56), (130.9585348, 0.57), (
        131.4185777, 0.58), (131.8729055, 0.59), (132.3216843, 0.6), (
        132.7650726, 0.61), (133.2032221, 0.62), (133.6362775, 0.63), (
        134.0643776, 0.64), (134.4876553, 0.65), (134.9062378, 0.66), (
        135.3202471, 0.67), (135.7298005, 0.68), (136.1350102, 0.69), (
        136.5359844, 0.7), (136.9328269, 0.71), (137.3256375, 0.72), (
        137.7145125, 0.73), (138.0995443, 0.74), (138.4808221, 0.75), (
        138.8584318, 0.76), (139.2324563, 0.77), (139.6029755, 0.78), (
        139.9700666, 0.79), (140.3338041, 0.8), (140.6942601, 0.81), (
        141.051504, 0.82), (141.4056032, 0.83), (141.7566227, 0.84), (
        142.1046256, 0.85), (142.4496727, 0.86), (142.7918231, 0.87), (
        143.131134, 0.88), (143.4676608, 0.89), (143.8014574, 0.9), (
        144.1325758, 0.91), (144.4610665, 0.92), (144.7869787, 0.93), (
        145.11036, 0.94), (145.4312566, 0.95), (145.7497136, 0.96), (
        146.0657745, 0.97), (146.3794819, 0.98), (146.690877, 0.99), (147.0, 
        1.0), (147.3068899, 1.01), (147.6115847, 1.02), (147.9141214, 1.03), (
        148.214536, 1.04), (148.5128636, 1.05), (148.8091384, 1.06), (
        149.1033937, 1.07), (149.3956619, 1.08), (149.6859746, 1.09), (
        149.9743627, 1.1), (150.2608563, 1.11), (150.5454847, 1.12), (
        150.8282766, 1.13), (151.1092599, 1.14), (151.3884619, 1.15), (
        151.6659091, 1.16), (151.9416278, 1.17), (152.2156431, 1.18), (
        152.48798, 1.19), (152.7586628, 1.2), (153.027715, 1.21), (153.29516, 
        1.22), (153.5610204, 1.23), (153.8253182, 1.24), (154.0880753, 1.25), (
        154.3493128, 1.26), (154.6090514, 1.27), (154.8673114, 1.28), (
        155.1241128, 1.29), (155.3794749, 1.3), (155.6334168, 1.31), (
        155.8859572, 1.32), (156.1371143, 1.33), (156.3869061, 1.34), (
        156.63535, 1.35), (156.8824632, 1.36), (157.1282625, 1.37), (
        157.3727646, 1.38), (157.6159854, 1.39), (157.857941, 1.4), (
        158.0986467, 1.41), (158.338118, 1.42), (158.5763697, 1.43), (
        158.8134166, 1.44), (159.049273, 1.45), (159.283953, 1.46), (
        159.5174705, 1.47), (159.7498392, 1.48), (159.9810723, 1.49), (
        160.211183, 1.5), (160.4401842, 1.51), (160.6680885, 1.52), (
        160.8949084, 1.53), (161.120656, 1.54), (161.3453433, 1.55), (
        161.568982, 1.56), (161.7915839, 1.57), (162.0131601, 1.58), (
        162.233722, 1.59), (162.4532804, 1.6), (162.6718462, 1.61), (162.88943, 
        1.62), (163.1060422, 1.63), (163.321693, 1.64), (163.5363926, 1.65), (
        163.7501509, 1.66), (163.9629777, 1.67), (164.1748825, 1.68), (
        164.3858748, 1.69), (164.5959639, 1.7), (164.8051589, 1.71), (
        165.0134688, 1.72), (165.2209025, 1.73), (165.4274688, 1.74), (
        165.6331761, 1.75), (165.838033, 1.76), (166.0420477, 1.77), (
        166.2452285, 1.78), (166.4475834, 1.79), (166.6491203, 1.8), (
        166.8498472, 1.81), (167.0497716, 1.82), (167.2489012, 1.83), (
        167.4472434, 1.84), (167.6448057, 1.85), (167.8415953, 1.86), (
        168.0376193, 1.87), (168.2328847, 1.88), (168.4273986, 1.89), (
        168.6211678, 1.9), (168.814199, 1.91), (169.0064988, 1.92), (
        169.1980739, 1.93), (169.3889307, 1.94), (169.5790756, 1.95), (
        169.7685148, 1.96), (169.9572546, 1.97), (170.145301, 1.98), (
        170.3326602, 1.99), (170.519338, 2.0), (170.7053403, 2.01), (
        170.8906729, 2.02), (171.0753416, 2.03), (171.2593518, 2.04), (
        171.4427093, 2.05), (171.6254194, 2.06), (171.8074877, 2.07), (
        171.9889194, 2.08), (172.1697198, 2.09), (172.349894, 2.1), (
        172.5294474, 2.11), (172.7083848, 2.12), (172.8867113, 2.13), (
        173.0644319, 2.14), (173.2415515, 2.15), (173.4180747, 2.16), (
        173.5940065, 2.17), (173.7693514, 2.18), (173.9441142, 2.19), (
        174.1182993, 2.2), (174.2919114, 2.21), (174.4649549, 2.22), (
        174.6374341, 2.23), (174.8093536, 2.24), (174.9807175, 2.25), (
        175.1515301, 2.26), (175.3217956, 2.27), (175.4915182, 2.28), (
        175.6607019, 2.29), (175.8293509, 2.3), (175.9974691, 2.31), (
        176.1650605, 2.32), (176.3321289, 2.33), (176.4986784, 2.34), (
        176.6647126, 2.35), (176.8302353, 2.36), (176.9952504, 2.37), (
        177.1597614, 2.38), (177.323772, 2.39), (177.4872859, 2.4), (
        177.6503066, 2.41), (177.8128376, 2.42), (177.9748824, 2.43), (
        178.1364444, 2.44), (178.2975271, 2.45), (178.4581339, 2.46), (
        178.618268, 2.47), (178.7779327, 2.48), (178.9371314, 2.49), (
        179.0958672, 2.5), (179.2541434, 2.51), (179.411963, 2.52), (
        179.5693292, 2.53), (179.7262452, 2.54), (179.8827138, 2.55), (
        180.0387383, 2.56), (180.1943215, 2.57), (180.3494664, 2.58), (
        180.504176, 2.59), (180.6584531, 2.6), (180.8123006, 2.61), (
        180.9657213, 2.62), (181.1187181, 2.63), (181.2712937, 2.64), (
        181.4234509, 2.65), (181.5751923, 2.66), (181.7265208, 2.67), (
        181.8774389, 2.68), (182.0279492, 2.69), (182.1780545, 2.7), (
        182.3277572, 2.71), (182.47706, 2.72), (182.6259654, 2.73), (
        182.7744758, 2.74), (182.9225938, 2.75), (183.0703218, 2.76), (
        183.2176623, 2.77), (183.3646176, 2.78), (183.5111902, 2.79), (
        183.6573824, 2.8), (183.8031965, 2.81), (183.948635, 2.82), (
        184.0936999, 2.83), (184.2383938, 2.84), (184.3827187, 2.85), (
        184.526677, 2.86), (184.6702707, 2.87), (184.8135022, 2.88), (
        184.9563736, 2.89), (185.0988871, 2.9), (185.2410446, 2.91), (
        185.3828485, 2.92), (185.5243007, 2.93), (185.6654033, 2.94), (
        185.8061584, 2.95), (185.946568, 2.96), (186.0866341, 2.97))
    
    # BCs
    input_dex.punch_velocity = 10
    input_dex.mass_scalling = 5e-6

    input_dex.friction_coefficient = 0.15

    # mesh
    input_dex.mesh_size = 3
    input_dex.local_mesh_size = 2

    input_dex_set = [input_dex]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\sim_op_2/1"
    # main_sim_output =  r"C:\Users\Kadmiel McForrester\OneDrive - University of Bath\Documents\build\batch_temp92"

    def get_sim_out_path(main_sim_output):
        if os.path.isdir(main_sim_output):
            sim_number = int(main_sim_output[-1])
            sim_number += 1
            main_sim_output = main_sim_output[:-1] + str(sim_number)
            return get_sim_out_path(main_sim_output)
        
        return main_sim_output
    
    main_sim_output = get_sim_out_path(main_sim_output)

    min_thickness, rmse, max_deviation = run_batch_sim(input_dex_set, main_sim_output)

    return min_thickness, rmse, max_deviation


def cup_baysien_optimisation():

    ideal_height = 10
    ideal_radius = 15
    ideal_profile_radius = 6.5

    def cup_function(punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap):

        # use suprocess here to run the sciprt, so no loop
        _, rmse, _ = cylindrical_cup_function(ideal_height, ideal_radius, ideal_profile_radius, punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap)
        print("sim complete")

        return -rmse

    # Bounded region of parameter space
    pbounds = {"punch_depth": (15, 25), "punch_profile_radius": (4, 7), "punch_min_radius": (10, 20),
                "die_profile_radius": (6, 9), "die_punch_gap": (1,5)}

    optimizer = BayesianOptimization(
        f=cup_function,
        pbounds=pbounds,
        random_state=1,
    )


    logger = JSONLogger(path=r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\op_logs.log")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    optimizer.maximize(
    init_points=2,
    n_iter=3,
    )

    print(optimizer.max)


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

    # Run simulation
    _, rmse, _ = cylindrical_cup_function(
        ideal_height, ideal_radius, ideal_profile_radius, 
        punch_depth, punch_profile_radius, punch_min_radius,
        die_profile_radius, die_punch_gap
    )

    data_out = {"rmse": -rmse}
    with open(r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\sim_output.json", "w") as f:
        json.dump(data_out, f)

# def main():

#     input_dex = InputeDex()

#     input_dex.solver_type = SolverType.STANDARD
    
#     input_dex.all_part_rotation = 90

#     # ideal part dimensions
#     input_dex.ideal_cup_radius = 20
#     input_dex.ideal_cup_height = 15
#     input_dex.ideal_cup_profile_radius = 5

#     # blank inputs
#     input_dex.blank_radius = 40
#     input_dex.blank_thickness = 1.1
#     input_dex.integration_points = 15

#     # punch inputs
#     input_dex.punch_profile_radius = 5 # 5
#     input_dex.punch_min_radius = 20

#     # die inputs
#     input_dex.die_profile_radius = 6.5
#     input_dex.die_min_radius = 22
#     input_dex.die_max_radius = input_dex.blank_radius + 5

#     # blank holder inputs
#     input_dex.blank_holder_height = 5
#     input_dex.blank_holder_profile_radius = 6.5
#     input_dex.blank_holder_min_radius = 22
#     input_dex.blank_holder_max_radius = input_dex.blank_radius + 5
#     input_dex.blank_holder_die_gap = input_dex.blank_thickness + 1

#     # punch
#     input_dex.cup_design_height = 20
#     input_dex.punch_depth = input_dex.cup_design_height + input_dex.trim_depth + input_dex.die_profile_radius
#     input_dex.die_height = input_dex.punch_depth + 10

#     # material inputs
#     input_dex.blank_material_name = "AA1050"
#     input_dex.density = 2.7e-06 
#     input_dex.youngs_modulus = 70000.0
#     input_dex.posissons_ratio = 0.33
#     input_dex.plastic_material_data = ((35.0, 0.0), (66.56588883, 0.01), (73.19453891, 0.02), (
#         77.6998518, 0.03), (81.21516632, 0.04), (84.1399572, 0.05), (
#         86.66656829, 0.06), (88.90387465, 0.07), (90.92007807, 0.08), (
#         92.76100124, 0.09), (94.45905775, 0.1), (96.03810066, 0.11), (
#         97.51624224, 0.12), (98.90759062, 0.13), (100.2233697, 0.14), (
#         101.4726693, 0.15), (102.6629639, 0.16), (103.8004812, 0.17), (
#         104.8904701, 0.18), (105.9373994, 0.19), (106.9451084, 0.2), (
#         107.9169222, 0.21), (108.855741, 0.22), (109.7641117, 0.23), (
#         110.6442836, 0.24), (111.4982544, 0.25), (112.3278069, 0.26), (
#         113.1345397, 0.27), (113.9198918, 0.28), (114.6851642, 0.29), (
#         115.431537, 0.3), (116.1600841, 0.31), (116.8717864, 0.32), (
#         117.567542, 0.33), (118.2481756, 0.34), (118.9144466, 0.35), (
#         119.5670557, 0.36), (120.2066512, 0.37), (120.833834, 0.38), (
#         121.4491624, 0.39), (122.053156, 0.4), (122.6462992, 0.41), (
#         123.2290448, 0.42), (123.8018162, 0.43), (124.3650101, 0.44), (
#         124.9189991, 0.45), (125.4641332, 0.46), (126.000742, 0.47), (
#         126.529136, 0.48), (127.0496083, 0.49), (127.5624356, 0.5), (128.06788, 
#         0.51), (128.5661893, 0.52), (129.0575987, 0.53), (129.5423312, 0.54), (
#         130.0205987, 0.55), (130.4926026, 0.56), (130.9585348, 0.57), (
#         131.4185777, 0.58), (131.8729055, 0.59), (132.3216843, 0.6), (
#         132.7650726, 0.61), (133.2032221, 0.62), (133.6362775, 0.63), (
#         134.0643776, 0.64), (134.4876553, 0.65), (134.9062378, 0.66), (
#         135.3202471, 0.67), (135.7298005, 0.68), (136.1350102, 0.69), (
#         136.5359844, 0.7), (136.9328269, 0.71), (137.3256375, 0.72), (
#         137.7145125, 0.73), (138.0995443, 0.74), (138.4808221, 0.75), (
#         138.8584318, 0.76), (139.2324563, 0.77), (139.6029755, 0.78), (
#         139.9700666, 0.79), (140.3338041, 0.8), (140.6942601, 0.81), (
#         141.051504, 0.82), (141.4056032, 0.83), (141.7566227, 0.84), (
#         142.1046256, 0.85), (142.4496727, 0.86), (142.7918231, 0.87), (
#         143.131134, 0.88), (143.4676608, 0.89), (143.8014574, 0.9), (
#         144.1325758, 0.91), (144.4610665, 0.92), (144.7869787, 0.93), (
#         145.11036, 0.94), (145.4312566, 0.95), (145.7497136, 0.96), (
#         146.0657745, 0.97), (146.3794819, 0.98), (146.690877, 0.99), (147.0, 
#         1.0), (147.3068899, 1.01), (147.6115847, 1.02), (147.9141214, 1.03), (
#         148.214536, 1.04), (148.5128636, 1.05), (148.8091384, 1.06), (
#         149.1033937, 1.07), (149.3956619, 1.08), (149.6859746, 1.09), (
#         149.9743627, 1.1), (150.2608563, 1.11), (150.5454847, 1.12), (
#         150.8282766, 1.13), (151.1092599, 1.14), (151.3884619, 1.15), (
#         151.6659091, 1.16), (151.9416278, 1.17), (152.2156431, 1.18), (
#         152.48798, 1.19), (152.7586628, 1.2), (153.027715, 1.21), (153.29516, 
#         1.22), (153.5610204, 1.23), (153.8253182, 1.24), (154.0880753, 1.25), (
#         154.3493128, 1.26), (154.6090514, 1.27), (154.8673114, 1.28), (
#         155.1241128, 1.29), (155.3794749, 1.3), (155.6334168, 1.31), (
#         155.8859572, 1.32), (156.1371143, 1.33), (156.3869061, 1.34), (
#         156.63535, 1.35), (156.8824632, 1.36), (157.1282625, 1.37), (
#         157.3727646, 1.38), (157.6159854, 1.39), (157.857941, 1.4), (
#         158.0986467, 1.41), (158.338118, 1.42), (158.5763697, 1.43), (
#         158.8134166, 1.44), (159.049273, 1.45), (159.283953, 1.46), (
#         159.5174705, 1.47), (159.7498392, 1.48), (159.9810723, 1.49), (
#         160.211183, 1.5), (160.4401842, 1.51), (160.6680885, 1.52), (
#         160.8949084, 1.53), (161.120656, 1.54), (161.3453433, 1.55), (
#         161.568982, 1.56), (161.7915839, 1.57), (162.0131601, 1.58), (
#         162.233722, 1.59), (162.4532804, 1.6), (162.6718462, 1.61), (162.88943, 
#         1.62), (163.1060422, 1.63), (163.321693, 1.64), (163.5363926, 1.65), (
#         163.7501509, 1.66), (163.9629777, 1.67), (164.1748825, 1.68), (
#         164.3858748, 1.69), (164.5959639, 1.7), (164.8051589, 1.71), (
#         165.0134688, 1.72), (165.2209025, 1.73), (165.4274688, 1.74), (
#         165.6331761, 1.75), (165.838033, 1.76), (166.0420477, 1.77), (
#         166.2452285, 1.78), (166.4475834, 1.79), (166.6491203, 1.8), (
#         166.8498472, 1.81), (167.0497716, 1.82), (167.2489012, 1.83), (
#         167.4472434, 1.84), (167.6448057, 1.85), (167.8415953, 1.86), (
#         168.0376193, 1.87), (168.2328847, 1.88), (168.4273986, 1.89), (
#         168.6211678, 1.9), (168.814199, 1.91), (169.0064988, 1.92), (
#         169.1980739, 1.93), (169.3889307, 1.94), (169.5790756, 1.95), (
#         169.7685148, 1.96), (169.9572546, 1.97), (170.145301, 1.98), (
#         170.3326602, 1.99), (170.519338, 2.0), (170.7053403, 2.01), (
#         170.8906729, 2.02), (171.0753416, 2.03), (171.2593518, 2.04), (
#         171.4427093, 2.05), (171.6254194, 2.06), (171.8074877, 2.07), (
#         171.9889194, 2.08), (172.1697198, 2.09), (172.349894, 2.1), (
#         172.5294474, 2.11), (172.7083848, 2.12), (172.8867113, 2.13), (
#         173.0644319, 2.14), (173.2415515, 2.15), (173.4180747, 2.16), (
#         173.5940065, 2.17), (173.7693514, 2.18), (173.9441142, 2.19), (
#         174.1182993, 2.2), (174.2919114, 2.21), (174.4649549, 2.22), (
#         174.6374341, 2.23), (174.8093536, 2.24), (174.9807175, 2.25), (
#         175.1515301, 2.26), (175.3217956, 2.27), (175.4915182, 2.28), (
#         175.6607019, 2.29), (175.8293509, 2.3), (175.9974691, 2.31), (
#         176.1650605, 2.32), (176.3321289, 2.33), (176.4986784, 2.34), (
#         176.6647126, 2.35), (176.8302353, 2.36), (176.9952504, 2.37), (
#         177.1597614, 2.38), (177.323772, 2.39), (177.4872859, 2.4), (
#         177.6503066, 2.41), (177.8128376, 2.42), (177.9748824, 2.43), (
#         178.1364444, 2.44), (178.2975271, 2.45), (178.4581339, 2.46), (
#         178.618268, 2.47), (178.7779327, 2.48), (178.9371314, 2.49), (
#         179.0958672, 2.5), (179.2541434, 2.51), (179.411963, 2.52), (
#         179.5693292, 2.53), (179.7262452, 2.54), (179.8827138, 2.55), (
#         180.0387383, 2.56), (180.1943215, 2.57), (180.3494664, 2.58), (
#         180.504176, 2.59), (180.6584531, 2.6), (180.8123006, 2.61), (
#         180.9657213, 2.62), (181.1187181, 2.63), (181.2712937, 2.64), (
#         181.4234509, 2.65), (181.5751923, 2.66), (181.7265208, 2.67), (
#         181.8774389, 2.68), (182.0279492, 2.69), (182.1780545, 2.7), (
#         182.3277572, 2.71), (182.47706, 2.72), (182.6259654, 2.73), (
#         182.7744758, 2.74), (182.9225938, 2.75), (183.0703218, 2.76), (
#         183.2176623, 2.77), (183.3646176, 2.78), (183.5111902, 2.79), (
#         183.6573824, 2.8), (183.8031965, 2.81), (183.948635, 2.82), (
#         184.0936999, 2.83), (184.2383938, 2.84), (184.3827187, 2.85), (
#         184.526677, 2.86), (184.6702707, 2.87), (184.8135022, 2.88), (
#         184.9563736, 2.89), (185.0988871, 2.9), (185.2410446, 2.91), (
#         185.3828485, 2.92), (185.5243007, 2.93), (185.6654033, 2.94), (
#         185.8061584, 2.95), (185.946568, 2.96), (186.0866341, 2.97))
    
#     # BCs
#     input_dex.punch_velocity = 10
#     input_dex.mass_scalling = 5e-6

#     # mesh
#     input_dex.mesh_size = 5

#     # input dex 2
#     input_dex_2: InputeDex = copy.deepcopy(input_dex)
#     input_dex_2.mesh_size = 3

#     # input dex 3
#     input_dex_3: InputeDex = copy.deepcopy(input_dex)
#     input_dex_3.mesh_size = 2

#     # input dex 4
#     input_dex_4: InputeDex = copy.deepcopy(input_dex)
#     input_dex_4.mesh_size = 0.5

#     # input dex 5
#     input_dex_5: InputeDex = copy.deepcopy(input_dex)
#     input_dex_5.mesh_size = 4

#     input_dex_set = [input_dex_3]

#     main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\batch_temp103"
#     # main_sim_output =  r"C:\Users\Kadmiel McForrester\OneDrive - University of Bath\Documents\build\batch_temp92"

#     min_thickness, rmse, max_deviation = run_batch_sim(input_dex_set, main_sim_output)



def main():

    input_dex = InputeDex()

    input_dex.solver_type = SolverType.STANDARD
    
    input_dex.all_part_rotation = 90

    # ideal part dimensions
    input_dex.ideal_cup_radius = 20
    input_dex.ideal_cup_height = 10
    input_dex.ideal_cup_profile_radius = 5

    # blank inputs
    input_dex.blank_radius = 55
    input_dex.blank_thickness = 1.1
    input_dex.integration_points = 15

    # punch inputs
    input_dex.punch_profile_radius = 5 # 5
    input_dex.punch_min_radius = 20

    # die inputs
    input_dex.die_profile_radius = 6.5
    input_dex.die_min_radius = 22
    input_dex.die_max_radius = input_dex.blank_radius + 10

    # blank holder inputs
    input_dex.blank_holder_height = 10
    input_dex.blank_holder_profile_radius = 6.5
    input_dex.blank_holder_min_radius = 22
    input_dex.blank_holder_max_radius = input_dex.blank_radius + 10
    input_dex.blank_holder_die_gap = input_dex.blank_thickness + 0.5

    # punch
    input_dex.cup_design_height = 10
    input_dex.punch_depth = 20
    input_dex.die_height = input_dex.punch_depth + 5

    # material inputs
    input_dex.blank_material_name = "AA1050"
    input_dex.density = 2.7e-06 
    input_dex.youngs_modulus = 70000.0
    input_dex.posissons_ratio = 0.33
    input_dex.plastic_material_data = ((35.0, 0.0), (66.56588883, 0.01), (73.19453891, 0.02), (
        77.6998518, 0.03), (81.21516632, 0.04), (84.1399572, 0.05), (
        86.66656829, 0.06), (88.90387465, 0.07), (90.92007807, 0.08), (
        92.76100124, 0.09), (94.45905775, 0.1), (96.03810066, 0.11), (
        97.51624224, 0.12), (98.90759062, 0.13), (100.2233697, 0.14), (
        101.4726693, 0.15), (102.6629639, 0.16), (103.8004812, 0.17), (
        104.8904701, 0.18), (105.9373994, 0.19), (106.9451084, 0.2), (
        107.9169222, 0.21), (108.855741, 0.22), (109.7641117, 0.23), (
        110.6442836, 0.24), (111.4982544, 0.25), (112.3278069, 0.26), (
        113.1345397, 0.27), (113.9198918, 0.28), (114.6851642, 0.29), (
        115.431537, 0.3), (116.1600841, 0.31), (116.8717864, 0.32), (
        117.567542, 0.33), (118.2481756, 0.34), (118.9144466, 0.35), (
        119.5670557, 0.36), (120.2066512, 0.37), (120.833834, 0.38), (
        121.4491624, 0.39), (122.053156, 0.4), (122.6462992, 0.41), (
        123.2290448, 0.42), (123.8018162, 0.43), (124.3650101, 0.44), (
        124.9189991, 0.45), (125.4641332, 0.46), (126.000742, 0.47), (
        126.529136, 0.48), (127.0496083, 0.49), (127.5624356, 0.5), (128.06788, 
        0.51), (128.5661893, 0.52), (129.0575987, 0.53), (129.5423312, 0.54), (
        130.0205987, 0.55), (130.4926026, 0.56), (130.9585348, 0.57), (
        131.4185777, 0.58), (131.8729055, 0.59), (132.3216843, 0.6), (
        132.7650726, 0.61), (133.2032221, 0.62), (133.6362775, 0.63), (
        134.0643776, 0.64), (134.4876553, 0.65), (134.9062378, 0.66), (
        135.3202471, 0.67), (135.7298005, 0.68), (136.1350102, 0.69), (
        136.5359844, 0.7), (136.9328269, 0.71), (137.3256375, 0.72), (
        137.7145125, 0.73), (138.0995443, 0.74), (138.4808221, 0.75), (
        138.8584318, 0.76), (139.2324563, 0.77), (139.6029755, 0.78), (
        139.9700666, 0.79), (140.3338041, 0.8), (140.6942601, 0.81), (
        141.051504, 0.82), (141.4056032, 0.83), (141.7566227, 0.84), (
        142.1046256, 0.85), (142.4496727, 0.86), (142.7918231, 0.87), (
        143.131134, 0.88), (143.4676608, 0.89), (143.8014574, 0.9), (
        144.1325758, 0.91), (144.4610665, 0.92), (144.7869787, 0.93), (
        145.11036, 0.94), (145.4312566, 0.95), (145.7497136, 0.96), (
        146.0657745, 0.97), (146.3794819, 0.98), (146.690877, 0.99), (147.0, 
        1.0), (147.3068899, 1.01), (147.6115847, 1.02), (147.9141214, 1.03), (
        148.214536, 1.04), (148.5128636, 1.05), (148.8091384, 1.06), (
        149.1033937, 1.07), (149.3956619, 1.08), (149.6859746, 1.09), (
        149.9743627, 1.1), (150.2608563, 1.11), (150.5454847, 1.12), (
        150.8282766, 1.13), (151.1092599, 1.14), (151.3884619, 1.15), (
        151.6659091, 1.16), (151.9416278, 1.17), (152.2156431, 1.18), (
        152.48798, 1.19), (152.7586628, 1.2), (153.027715, 1.21), (153.29516, 
        1.22), (153.5610204, 1.23), (153.8253182, 1.24), (154.0880753, 1.25), (
        154.3493128, 1.26), (154.6090514, 1.27), (154.8673114, 1.28), (
        155.1241128, 1.29), (155.3794749, 1.3), (155.6334168, 1.31), (
        155.8859572, 1.32), (156.1371143, 1.33), (156.3869061, 1.34), (
        156.63535, 1.35), (156.8824632, 1.36), (157.1282625, 1.37), (
        157.3727646, 1.38), (157.6159854, 1.39), (157.857941, 1.4), (
        158.0986467, 1.41), (158.338118, 1.42), (158.5763697, 1.43), (
        158.8134166, 1.44), (159.049273, 1.45), (159.283953, 1.46), (
        159.5174705, 1.47), (159.7498392, 1.48), (159.9810723, 1.49), (
        160.211183, 1.5), (160.4401842, 1.51), (160.6680885, 1.52), (
        160.8949084, 1.53), (161.120656, 1.54), (161.3453433, 1.55), (
        161.568982, 1.56), (161.7915839, 1.57), (162.0131601, 1.58), (
        162.233722, 1.59), (162.4532804, 1.6), (162.6718462, 1.61), (162.88943, 
        1.62), (163.1060422, 1.63), (163.321693, 1.64), (163.5363926, 1.65), (
        163.7501509, 1.66), (163.9629777, 1.67), (164.1748825, 1.68), (
        164.3858748, 1.69), (164.5959639, 1.7), (164.8051589, 1.71), (
        165.0134688, 1.72), (165.2209025, 1.73), (165.4274688, 1.74), (
        165.6331761, 1.75), (165.838033, 1.76), (166.0420477, 1.77), (
        166.2452285, 1.78), (166.4475834, 1.79), (166.6491203, 1.8), (
        166.8498472, 1.81), (167.0497716, 1.82), (167.2489012, 1.83), (
        167.4472434, 1.84), (167.6448057, 1.85), (167.8415953, 1.86), (
        168.0376193, 1.87), (168.2328847, 1.88), (168.4273986, 1.89), (
        168.6211678, 1.9), (168.814199, 1.91), (169.0064988, 1.92), (
        169.1980739, 1.93), (169.3889307, 1.94), (169.5790756, 1.95), (
        169.7685148, 1.96), (169.9572546, 1.97), (170.145301, 1.98), (
        170.3326602, 1.99), (170.519338, 2.0), (170.7053403, 2.01), (
        170.8906729, 2.02), (171.0753416, 2.03), (171.2593518, 2.04), (
        171.4427093, 2.05), (171.6254194, 2.06), (171.8074877, 2.07), (
        171.9889194, 2.08), (172.1697198, 2.09), (172.349894, 2.1), (
        172.5294474, 2.11), (172.7083848, 2.12), (172.8867113, 2.13), (
        173.0644319, 2.14), (173.2415515, 2.15), (173.4180747, 2.16), (
        173.5940065, 2.17), (173.7693514, 2.18), (173.9441142, 2.19), (
        174.1182993, 2.2), (174.2919114, 2.21), (174.4649549, 2.22), (
        174.6374341, 2.23), (174.8093536, 2.24), (174.9807175, 2.25), (
        175.1515301, 2.26), (175.3217956, 2.27), (175.4915182, 2.28), (
        175.6607019, 2.29), (175.8293509, 2.3), (175.9974691, 2.31), (
        176.1650605, 2.32), (176.3321289, 2.33), (176.4986784, 2.34), (
        176.6647126, 2.35), (176.8302353, 2.36), (176.9952504, 2.37), (
        177.1597614, 2.38), (177.323772, 2.39), (177.4872859, 2.4), (
        177.6503066, 2.41), (177.8128376, 2.42), (177.9748824, 2.43), (
        178.1364444, 2.44), (178.2975271, 2.45), (178.4581339, 2.46), (
        178.618268, 2.47), (178.7779327, 2.48), (178.9371314, 2.49), (
        179.0958672, 2.5), (179.2541434, 2.51), (179.411963, 2.52), (
        179.5693292, 2.53), (179.7262452, 2.54), (179.8827138, 2.55), (
        180.0387383, 2.56), (180.1943215, 2.57), (180.3494664, 2.58), (
        180.504176, 2.59), (180.6584531, 2.6), (180.8123006, 2.61), (
        180.9657213, 2.62), (181.1187181, 2.63), (181.2712937, 2.64), (
        181.4234509, 2.65), (181.5751923, 2.66), (181.7265208, 2.67), (
        181.8774389, 2.68), (182.0279492, 2.69), (182.1780545, 2.7), (
        182.3277572, 2.71), (182.47706, 2.72), (182.6259654, 2.73), (
        182.7744758, 2.74), (182.9225938, 2.75), (183.0703218, 2.76), (
        183.2176623, 2.77), (183.3646176, 2.78), (183.5111902, 2.79), (
        183.6573824, 2.8), (183.8031965, 2.81), (183.948635, 2.82), (
        184.0936999, 2.83), (184.2383938, 2.84), (184.3827187, 2.85), (
        184.526677, 2.86), (184.6702707, 2.87), (184.8135022, 2.88), (
        184.9563736, 2.89), (185.0988871, 2.9), (185.2410446, 2.91), (
        185.3828485, 2.92), (185.5243007, 2.93), (185.6654033, 2.94), (
        185.8061584, 2.95), (185.946568, 2.96), (186.0866341, 2.97))
    
    # BCs
    input_dex.punch_velocity = 15
    input_dex.mass_scalling = 5e-6

    input_dex.friction_coefficient = 0.15

    # mesh
    input_dex.mesh_size = 2
    input_dex.local_mesh_size = 0.75

    # input dex 2
    input_dex_2: InputeDex = copy.deepcopy(input_dex)
    input_dex_2.mesh_size = 3

    # input dex 3
    input_dex_3: InputeDex = copy.deepcopy(input_dex)
    input_dex_3.mesh_size = 2

    # input dex 4
    input_dex_4: InputeDex = copy.deepcopy(input_dex)
    input_dex_4.mesh_size = 1

    # input dex 5
    input_dex_5: InputeDex = copy.deepcopy(input_dex)
    input_dex_5.mesh_size = 0.75

    # input dex 5
    input_dex_6: InputeDex = copy.deepcopy(input_dex)
    input_dex_6.mesh_size = 0.3

    input_dex_set = [input_dex]

    main_sim_output =  r"C:\Users\kam97\OneDrive - University of Bath\Documents\build\build_temp140"
    # main_sim_output =  r"C:\Users\Kadmiel McForrester\OneDrive - University of Bath\Documents\build\batch_temp92"

    run_batch_sim(input_dex_set, main_sim_output)

if __name__ == "__main__":
    
    run_cup()

    # cup_baysien_optimisation()

    # main()

    # create_mesh(SolverType.STANDARD, 2)


    # export_final_node_positions(41.5, 30)
    # os.chdir(r"C:\Users\Kadmiel McForrester\OneDrive - University of Bath\Documents\build\batch_temp103\sim_1")
    # compare_meshes()
    
    # post_pro_strain_eq()
    # post_pro_misses_stress()
    # post_pro_spring_back_dev()
    # plt.close('all')
    # post_pro_thickness()
    # plot_thickness_variation()
    # post_pro_punch_reaction()
    # post_pro_energy()
    # plot_energy_data()
    # plot_rf_data()
    # energy_check()
    # post_pro_strain()
    # plot_strains()