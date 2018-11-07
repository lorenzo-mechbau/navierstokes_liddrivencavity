#================================================================================================================================
#  Start Program
#================================================================================================================================

useHermite = False

numberOfElements = 5

cavitySize = 1.0
elementSize = cavitySize/numberOfElements

Re = 1000.0

#RBS = False
RBS = True

outputFrequency = 1 # Result output frequency

setupOutput = True
progressDiagnostics = True
debugLevel = 3

startTime = 0.0
stopTime  = 50.001
timeStep  = 0.1

# Lid velocity parameters
A = 0.5
B = 0.5
C = 10.0

# Material properties
fluidDensity  = 1.0
fluidDynamicViscosity = fluidDensity*(A+B)/Re

fluidPInit = 0.0

# Set solver parameters
fluidDynamicSolverTheta    = [0.5]
nonlinearMaximumIterations      = 100000000 #default: 100000
nonlinearRelativeTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearSolutionTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearAbsoluteTolerance      = 1.0E-8    #default: 1.0E-10
nonlinearMaxFunctionEvaluations = 10000
nonlinearLinesearchAlpha        = 1.0
linearMaximumIterations      = 100000000 #default: 100000
linearRelativeTolerance      = 1.0E-6    #default: 1.0E-05
linearAbsoluteTolerance      = 1.0E-6    #default: 1.0E-10
linearDivergenceTolerance    = 1.0E5     #default: 1.0E5
linearRestartValue           = 30        #default: 30

#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================


if (useHermite):
    numberOfNodesXi = 2
else:
    numberOfNodesXi = 3

numberOfFluidNodes = (numberOfElements*(numberOfNodesXi-1)+1)*(numberOfElements*(numberOfNodesXi-1)+1)
numberOfFluidElements = numberOfElements* numberOfElements

fluidCoordinateSystemUserNumber = 1

fluidRegionUserNumber = 1

linearBasisUserNumber = 1
quadraticBasisUserNumber = 2
hermiteBasisUserNumber = 3
interfaceQuadraticBasisUserNumber = 4
interfaceHermiteBasisUserNumber = 5

fluidMeshUserNumber = 1

fluidDecompositionUserNumber = 1

fluidGeometricFieldUserNumber     = 1
fluidEquationsSetFieldUserNumber = 2
fluidDependentFieldUserNumber = 3
fluidMaterialsFieldUserNumber = 4
fluidIndependentFieldUserNumber = 5
bcCellMLModelsFieldUserNumber = 6
bcCellMLStateFieldUserNumber = 7
bcCellMLParametersFieldUserNumber = 8
bcCellMLIntermediateFieldUserNumber = 9

fluidEquationsSetUserNumber  = 1

bcCellMLUserNumber = 1

fluidProblemUserNumber = 1

#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,csv,time,sys,os,pdb
from opencmiss.iron import iron

# Path from command line argument or cd
if len(sys.argv) > 1:
    file_root_directory = sys.argv[1]
else:
    file_root_directory = os.path.dirname(__file__)

# Filename from command line argument or default
if len(sys.argv) > 2:
    cellml_fileName = str(sys.argv[2])
else:
    cellml_fileName = "input/fixedlidvelocity.cellml"

cellml_file = os.path.join(file_root_directory, cellml_fileName)

# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
iron.OutputSetOn("Testing")

# Get the number of computational nodes and this computational node number
#computationEnvironment = iron.ComputationEnvironment()
#numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
#computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = iron.EquationsOutputTypes.NONE
#fluidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#fluidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#fluidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
fluidDynamicSolverOutputType = iron.SolverOutputTypes.NONE
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
fluidNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidLinearSolverOutputType = iron.SolverOutputTypes.NONE
#fluidLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidLinearSolverOutputType = iron.SolverOutputTypes.MATRIX

if (setupOutput):
    print('SUMMARY')
    print('=======')
    print(' ')
    print('  Temporal parameters')
    print('  -------------------')
    print(' ')
    print('  Start time:     %.3f s' % (startTime))
    print('  Stop time:      %.3f s' % (stopTime))
    print('  Time increment: %.5f s' % (timeStep))
    print(' ')
    print('  Material parameters')
    print('  -------------------')
    print(' ')
    print('    Fluid:')
    print('      Dynamic viscosity: {0:.3f} kg.m^-1.s^-1'.format(fluidDynamicViscosity))
    print('      Density: {0:.3f} kg.m^-3'.format(fluidDensity))
    print(' ')
    print('  Mesh parameters')
    print('  -------------------')
    print(' ')
    print('    Use Hermite: {}'.format(useHermite))
    print('    Fluid:')
    print('      Number of elements: {0:d}'.format(numberOfElements))
    print('      Number of nodes: {0:d}'.format(numberOfFluidNodes))
    print('      Number of elements: {0:d}'.format(numberOfFluidElements))

#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

# Create a RC coordinate system for the fluid region
fluidCoordinateSystem = iron.CoordinateSystem()
fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber)
fluidCoordinateSystem.DimensionSet(2)
fluidCoordinateSystem.CreateFinish()
if (progressDiagnostics):
    print('Coordinate systems ... Done')

#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

# Create a fluid region
fluidRegion = iron.Region()
fluidRegion.CreateStart(fluidRegionUserNumber,iron.WorldRegion)
fluidRegion.label = 'FluidRegion'
fluidRegion.coordinateSystem = fluidCoordinateSystem
fluidRegion.CreateFinish()

if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')

linearBasis = iron.Basis()
linearBasis.CreateStart(linearBasisUserNumber)
linearBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
linearBasis.numberOfXi = 2
linearBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
linearBasis.quadratureNumberOfGaussXi = [3]*2
linearBasis.CreateFinish()
if (useHermite):
    numberOfNodesXi = 2
    hermiteBasis = iron.Basis()
    hermiteBasis.CreateStart(hermiteBasisUserNumber)
    hermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    hermiteBasis.numberOfXi = 2
    hermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*2
    hermiteBasis.quadratureNumberOfGaussXi = [4]*2
    hermiteBasis.CreateFinish()
else:
    numberOfNodesXi = 3
    quadraticBasis = iron.Basis()
    quadraticBasis.CreateStart(quadraticBasisUserNumber)
    quadraticBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    quadraticBasis.numberOfXi = 2
    quadraticBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*2
    quadraticBasis.quadratureNumberOfGaussXi = [3]*2
    quadraticBasis.CreateFinish()

if (progressDiagnostics):
    print('Basis functions ... Done')

#================================================================================================================================
#  Mesh
#================================================================================================================================

if (progressDiagnostics):
    print('Meshes ...')

fluidNodes = iron.Nodes()
fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
fluidNodes.CreateFinish()

fluidMesh = iron.Mesh()
fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,2)
fluidMesh.NumberOfElementsSet(numberOfFluidElements)
fluidMesh.NumberOfComponentsSet(2)

fluidLinearElements = iron.MeshElements()
fluidLinearElements.CreateStart(fluidMesh,1,linearBasis)
if (useHermite):
    fluidHermiteElements = iron.MeshElements()
    fluidHermiteElements.CreateStart(fluidMesh,2,hermiteBasis)
else:
    fluidQuadraticElements = iron.MeshElements()
    fluidQuadraticElements.CreateStart(fluidMesh,2,quadraticBasis)

# Fluid mesh elements
if (debugLevel > 2):
    print('  Fluid Elements:')
for yElementIdx in range(1,numberOfElements+1):
    for xElementIdx in range(1,numberOfElements+1):
            elementNumber = xElementIdx+(yElementIdx-1)*numberOfElements
            localNodes1= (xElementIdx-1)*(numberOfNodesXi-1)+1+\
                         (yElementIdx-1)*(numberOfElements*(numberOfNodesXi-1)+1)*(numberOfNodesXi-1)
            localNodes3=localNodes1+(numberOfNodesXi-1)
            localNodes7=localNodes1+(numberOfElements*(numberOfNodesXi-1)+1)*(numberOfNodesXi-1)
            localNodes9=localNodes7+(numberOfNodesXi-1)
            fluidLinearElements.NodesSet(elementNumber,[localNodes1,localNodes3,localNodes7,localNodes9])
            if (useHermite):
                fluidHermiteElements.NodesSet(elementNumber,[localNodes1,localNodes3,localNodes7,localNodes9])
                if (debugLevel > 2):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d' % \
                          (elementNumber,localNodes1,localNodes3,localNodes7,localNodes9))
            else:
                localNodes2=localNodes1+1
                localNodes4=localNodes1+(numberOfElements*(numberOfNodesXi-1)+1)
                localNodes5=localNodes4+1
                localNodes6=localNodes5+1
                localNodes8=localNodes7+1
                fluidQuadraticElements.NodesSet(elementNumber,[localNodes1,localNodes2,localNodes3,localNodes4, \
                                                               localNodes5,localNodes6,localNodes7,localNodes8,localNodes9])
                if (debugLevel > 2):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (elementNumber,localNodes1,localNodes2,localNodes3,localNodes4,localNodes5,\
                           localNodes6,localNodes7,localNodes8,localNodes9))
else:
    print('Not implemented.')
    exit

fluidLinearElements.CreateFinish()
if (useHermite):
    fluidHermiteElements.CreateFinish()
else:
    fluidQuadraticElements.CreateFinish()

fluidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')

# Create a decomposition for the fluid mesh
fluidDecomposition = iron.Decomposition()
fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
fluidDecomposition.CalculateFacesSet(True)
fluidDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')

#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')

# Start to create a default (geometric) field on the fluid region
fluidGeometricField = iron.Field()
fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
# Set the decomposition to use
if (useHermite):
    fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
# Set the decomposition
fluidGeometricField.MeshDecompositionSet(fluidDecomposition)
# Set the scaling to use
fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGeometry')
# Set the domain to be used by the field components.
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
# Finish creating the second field
fluidGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')

if (progressDiagnostics):
    print('Geometric Parameters ...')

if (debugLevel > 2):
    print('  Fluid Nodes:')
for yNodeIdx in range(1,numberOfElements*(numberOfNodesXi-1)+2):
    for xNodeIdx in range(1,numberOfElements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx+(yNodeIdx-1)*(numberOfElements*(numberOfNodesXi-1)+1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            xPosition = float(xNodeIdx-1)/float(numberOfElements*(numberOfNodesXi-1))*cavitySize
            yPosition = float(yNodeIdx-1)/float(numberOfElements*(numberOfNodesXi-1))*cavitySize
            fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                         1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
            fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                         1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))
            if (useHermite):
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,1.0)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,0.0)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,1.0)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,0.0)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,0.0)
                if (debugLevel > 2):
                    print('        S1 derivative    = [ %.2f, %.2f ]' % (1.0,0.0))
                    print('        S2 derivative    = [ %.2f, %.2f ]' % (0.0,1.0))
                    print('        S1xS2 derivative = [ %.2f, %.2f ]' % (0.0,0.0))

# Update fields
fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

#================================================================================================================================
#  Equations Set
#================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

# Create the equations set for the fluid region - Navier-Stokes
fluidEquationsSetField = iron.Field()
fluidEquationsSet = iron.EquationsSet()
if RBS:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
else:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.TRANSIENT_NAVIER_STOKES]
fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                              fluidEquationsSetSpecification,fluidEquationsSetFieldUserNumber,
                              fluidEquationsSetField)
fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
fluidEquationsSet.CreateFinish()

if RBS:
    # Set max CFL number (default 1.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0E20)
    # Set time increment (default 0.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,3,timeStep)
    # Set stabilisation type (default 1.0 = RBS)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,4,1.0)

if (progressDiagnostics):
    print('Equations Sets ... Done')


#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

# Create the equations set dependent field variables for dynamic Navier-Stokes
fluidDependentField = iron.Field()
fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber,fluidDependentField)
fluidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,3):
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,2)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,2)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,1)
# Finish the equations set dependent field variables
fluidEquationsSet.DependentCreateFinish()

# Initialise the fluid dependent field
for componentIdx in range(1,3):
    fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)
# Initialise pressure component
fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,fluidPInit)

# Update dependent field
fluidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')

#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

# Create the equations set materials field variables for dynamic Navier-Stokes
fluidMaterialsField = iron.Field()
fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber,fluidMaterialsField)
# Finish the equations set materials field variables
fluidEquationsSet.MaterialsCreateFinish()
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,fluidDynamicViscosity)
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,fluidDensity)

if (progressDiagnostics):
    print('Materials Fields ... Done')

#================================================================================================================================
#  Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

# Fluid equations
fluidEquations = iron.Equations()
fluidEquationsSet.EquationsCreateStart(fluidEquations)
fluidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
fluidEquations.outputType = fluidEquationsOutputType
fluidEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

#================================================================================================================================
#  CellML
#================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

# Create CellML equations for the temporal boundary conditions
bcCellML = iron.CellML()
bcCellML.CreateStart(bcCellMLUserNumber,fluidRegion)
bcCellMLIdx = bcCellML.ModelImport(cellml_file)
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/A")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/B")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/C")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/x")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/y")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/lidx")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/lidy")
bcCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
bcCellML.FieldMapsCreateStart()
# Map geometric field to x, y and z
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/x",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/y",iron.FieldParameterSetTypes.VALUES)
# Map fluid velocity to lidx and lidy to ensure depndent field isn't cleared when the velocities are copied back
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/lidx",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/lidy",iron.FieldParameterSetTypes.VALUES)
# Map inletx, inlety to dependent field
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/lidx",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/lidy",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES)
bcCellML.FieldMapsCreateFinish()


# Create the CellML models field
bcCellMLModelsField = iron.Field()
bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"BCModelMap")
bcCellML.ModelsFieldCreateFinish()

# Only evaluate BC on inlet nodes
bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0)
if (debugLevel > 2):
    print('  CellML Boundary Conditions:')
    print('    Lid Model Set:')
for xNodeIdx in range(2,numberOfElements*(numberOfNodesXi-1)+1):
    nodeNumber = xNodeIdx+(numberOfElements*(numberOfNodesXi-1)+1)*numberOfElements*(numberOfNodesXi-1)
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
    if (nodeDomain == computationalNodeNumber):
        bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                       1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
        if (debugLevel > 2):
            print('      Node        %d:' % (nodeNumber))

# Create the CellML state field
bcCellMLStateField = iron.Field()
bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCState")
bcCellML.StateFieldCreateFinish()

# Create the CellML parameters field
bcCellMLParametersField = iron.Field()
bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"BCParameters")
bcCellML.ParametersFieldCreateFinish()

# Get the component numbers
AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/A")
BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/B")
CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/C")
# Set up the parameters field
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,AComponentNumber,A)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,BComponentNumber,B)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,CComponentNumber,C)

# Create the CELL intermediate field
bcCellMLIntermediateField = iron.Field()
bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCIntermediate")
bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

#================================================================================================================================
#  Problem
#================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a fluid problem
fluidProblem = iron.Problem()
if RBS:
    fluidProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                 iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                 iron.ProblemSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
else:
    fluidProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                 iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                 iron.ProblemSubtypes.TRANSIENT_NAVIER_STOKES]
fluidProblem.CreateStart(fluidProblemUserNumber,fluidProblemSpecification)
fluidProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

#================================================================================================================================
#  Control Loop
#================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fluid problem control loop
fluidControlLoop = iron.ControlLoop()
fluidProblem.ControlLoopCreateStart()
fluidProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],fluidControlLoop)
fluidControlLoop.LabelSet('TimeLoop')
fluidControlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
fluidControlLoop.TimesSet(startTime,stopTime,timeStep)
fluidControlLoop.TimeInputSet(2)
fluidControlLoop.TimeOutputSet(outputFrequency)
fluidProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = iron.Solver()
fluidDynamicSolver = iron.Solver()
fluidNonlinearSolver = iron.Solver()
fluidLinearSolver = iron.Solver()

fluidProblem.SolversCreateStart()
# Solvers for a Navier Stokes problem
# Get the BC CellML solver
fluidProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
# Get the dynamic solver
fluidProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fluidDynamicSolver)
fluidDynamicSolver.OutputTypeSet(fluidDynamicSolverOutputType)
fluidDynamicSolver.DynamicThetaSet(fluidDynamicSolverTheta)
# Get the dynamic nonlinear solver
fluidDynamicSolver.DynamicNonlinearSolverGet(fluidNonlinearSolver)
fluidNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
fluidNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
#fluidNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
fluidNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
fluidNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.QUADRATIC)
fluidNonlinearSolver.OutputTypeSet(fluidNonlinearSolverOutputType)
fluidNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
fluidNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
fluidNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
fluidNonlinearSolver.NewtonSolutionToleranceSet(nonlinearSolutionTolerance)
#fluidNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
# Get the dynamic nonlinear linear solver
fluidNonlinearSolver.NewtonLinearSolverGet(fluidLinearSolver)
#fluidLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
#fluidLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
#fluidLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
#fluidLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
#fluidLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
fluidLinearSolver.OutputTypeSet(fluidLinearSolverOutputType)
fluidLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
fluidLinearSolver.LibraryTypeSet(iron.SolverLibraries.MUMPS)
# Finish the creation of the problem solver
fluidProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

#================================================================================================================================
#  CellML Equations
#================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

# Create CellML equations and add BC equations to the solver
bcEquations = iron.CellMLEquations()
fluidProblem.CellMLEquationsCreateStart()
bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
fluidProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fluid problem solver equations
fluidProblem.SolverEquationsCreateStart()
# Get the fluid dynamic solver equations
fluidSolverEquations = iron.SolverEquations()
fluidDynamicSolver.SolverEquationsGet(fluidSolverEquations)
fluidSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
fluidEquationsSetIndex = fluidSolverEquations.EquationsSetAdd(fluidEquationsSet)
# Finish the creation of the fluid problem solver equations
fluidProblem.SolverEquationsCreateFinish()

if (progressDiagnostics):
    print('Solver Equations ...')

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fluid boundary conditions
fluidBoundaryConditions = iron.BoundaryConditions()
fluidSolverEquations.BoundaryConditionsCreateStart(fluidBoundaryConditions)
if (debugLevel > 2):
    print('  Fluid Boundary Conditions:')
    print('    Wall Boundary conditions:')
# Set boundary conditions on the bottom of the cavity
for xNodeIdx in range(1,numberOfElements*(numberOfNodesXi-1)+2):
    nodeNumber = xNodeIdx
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
    if (nodeDomain == computationalNodeNumber):
        fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                        iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                        nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                        iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                        nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debugLevel > 2):
            print('      Node        %d:' % (nodeNumber))
            print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))
        if (useHermite):
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
# Set boundary conditions on the left and right edges of the cavity
nodeNumbers = [0,0]
nodeDomains = [0,0]
for yNodeIdx in range(2,numberOfElements*(numberOfNodesXi-1)+2):
    nodeNumbers[0] = (yNodeIdx-1)*(numberOfElements*(numberOfNodesXi-1)+1)+1
    nodeDomains[0] = fluidDecomposition.NodeDomainGet(nodeNumbers[0],2)
    nodeNumbers[1] = yNodeIdx*(numberOfElements*(numberOfNodesXi-1)+1)
    nodeDomains[1] = fluidDecomposition.NodeDomainGet(nodeNumbers[1],2)
    for sideIdx in [0,1]:
        if (nodeDomains[sideIdx] == computationalNodeNumber):
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumbers[sideIdx],1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumbers[sideIdx],2,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumbers[sideIdx]))
                print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))
            if (useHermite):
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                                nodeNumbers[sideIdx],1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                                nodeNumbers[sideIdx],2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                                nodeNumbers[sideIdx],1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                                nodeNumbers[sideIdx],2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                                nodeNumbers[sideIdx],1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                                nodeNumbers[sideIdx],2,iron.BoundaryConditionsTypes.FIXED,0.0)
if (debugLevel > 2):
    print('    Lid Boundary conditions:')
# Set boundary conditions on the lid of the cavity
for xNodeIdx in range(2,numberOfElements*(numberOfNodesXi-1)+1):
    nodeNumber = xNodeIdx+(numberOfElements*(numberOfNodesXi-1)+1)*numberOfElements*(numberOfNodesXi-1)
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
    if (nodeDomain == computationalNodeNumber):
        fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                        iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                        nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
        fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                        iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                        nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
        if (debugLevel > 2):
            print('      Node        %d:' % (nodeNumber))
            print('         Velocity         = [ %.2f, %.2f ]' % (A+B,0.0))
        if (useHermite):
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
# Set pressure boundary condition on the bottom left node
if (debugLevel > 2):
    print('    Pressure Boundary conditions:')
nodeDomain = fluidDecomposition.NodeDomainGet(1,1)
if (nodeDomain == computationalNodeNumber):
    fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                    1,3,iron.BoundaryConditionsTypes.FIXED,0.0)
if (debugLevel > 2):
    print('      Node        %d:' % (1))

# Finish fluid boundary conditions
fluidSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#quit()

# Solve the problem
print('Solving problem...')
start = time.time()
fluidProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' %elapsed)
print('Problem solved!')
print('#')

print("Exporting CMGUI data")
# Export results
fields = iron.Fields()
fields.CreateRegion(fluidRegion)
fields.NodesExport("LidDrivenEnd","FORTRAN")
fields.ElementsExport("LidDrivenEnd","FORTRAN")
fields.Finalise()

#================================================================================================================================
#  Finish Program
#================================================================================================================================
# Finalise OpenCMISS-Iron
iron.Finalise()
