# -*- coding: utf-8 -*-
# env/geo3D_sim04
#########################

# author: arkriger - 2023 - 2026
# github: https://github.com/AdrianKriger/geo3D

# Consolidated helper functions for geo3D. OpenFOAM configurations.
# Supports both RANS (Steady-state) and URANS (Transient) strategies.
#########################

import os
import math
import numpy as np
import tempfile
import trimesh

# Mappings for user inputs
NU_MAP = {"sea-level": 1.5e-05, "elevated": 1.8e-05}
#Z0_MAP = {"coastal": 0.001, "open": 0.03, "suburb": 0.1, "city": 0.3}
# Updated Z0_MAP based on European Wind Atlas standards
Z0_MAP = {
    "water":                     0.0002, # Class 0: Open water
    "airport_runway":            0.0024, # Class 0.5: Concrete/mowed grass
    "open":                      0.03,   # Class 1: Open fields, no fences
    "agri_sheltered":            0.1,    # Class 2: Agri with 500m hedgerow spacing
    "rural":                     0.2,    # Class 2.5: Many houses/shrubs
    "village, town and suburb":  0.4,    # Class 3: Villages, small towns, or forests
    "urban suburb and city":     0.8,    # Class 3.5: Larger cities with tall buildings
    "metro":                     1.2     # Class 4: Skyscrapers/Metropolitan areas
}

def write_openfoam_case(case_path, extent, x_off, y_off, max_h, buildings, nu, z0, wind_speed, wind_deg, mode='RANS'):
    """
    Consolidated master case writer.
    mode: 'RANS' (Steady simpleFoam) or 'URANS' (Transient incompressibleFluid)
    """
    nu_val = NU_MAP.get(nu, 1.5e-05)
    z0_val = Z0_MAP.get(z0, 0.3)
    mode = mode.upper()

    # 1. Create Folder Tree
    for folder in ["0", "constant/geometry", "system"]:
        os.makedirs(os.path.join(case_path, folder), exist_ok=True)

    # STL Rotation & Export (Shared logic)
    rot = wind_deg - 270
    rot_rad = np.radians(rot)
    rotation_matrix = trimesh.transformations.rotation_matrix(rot_rad, [0, 0, 1])
    # 3. Apply the rotation to your final_model
    buildings.apply_transform(rotation_matrix)
    
    buildings.update_faces(buildings.unique_faces())
    buildings.fix_normals(multibody=True)

    buildings.export(os.path.join(case_path, "constant/geometry/buildings.obj"))
    
    # 2. 0/ Directory (Physics)
    write_u_file(case_path, wind_speed, z0_val)
    write_p_file(case_path)
    write_k_file(case_path, wind_speed, z0_val)
    write_nut_file(case_path, z0_val)
    
    if mode == 'URANS':
        write_omega_file(case_path, wind_speed, z0_val)
    else:
        write_epsilon_file(case_path, wind_speed, max_h, z0_val)
    
    # 3. constant/ Directory
    write_physical_properties(case_path, nu_val)
    write_momentum_transport(case_path, mode)
    
    # 4. system/ Directory
    write_block_mesh(case_path, extent, x_off, y_off, max_h)
    write_control_dict(case_path, mode)
    write_fv_schemes(case_path, mode)
    write_fv_solution(case_path, mode)
    write_mesh_quality_dict(case_path)
    write_surface_features(case_path)
    write_snappy_hex_mesh(case_path, extent, x_off, y_off, max_h)

def write_header(file_class, file_object):
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {file_class};
    object      {file_object};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"""

# --- PHYSICS FILES (0/) ---

def write_u_file(case_path, wind_speed, z0_value):
    header = write_header("volVectorField", "U")
    content = f"""{header}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform ({wind_speed} 0 0);
boundaryField
{{
    inlet {{ type atmBoundaryLayerInletVelocity; flowDir (1 0 0); zDir (0 0 1); Uref {wind_speed}; Zref 10.0; z0 uniform {z0_value}; zGround uniform 0.0; value $internalField; }}
    outlet {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
    ground {{ type noSlip; }}
    buildings {{ type noSlip; }}
    frontAndBack {{ type symmetry; }}
    sky {{ type slip; }}
}}"""
    with open(os.path.join(case_path, "0/U"), "w") as f: f.write(content)

def write_p_file(case_path):
    header = write_header("volScalarField", "p")
    content = f"""{header}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    inlet {{ type zeroGradient; }}
    outlet {{ type totalPressure; p0 uniform 0; value uniform 0; }}
    ground {{ type zeroGradient; }}
    buildings {{ type zeroGradient; }}
    frontAndBack {{ type symmetry; }}
    sky {{ type zeroGradient; }}
}}"""
    with open(os.path.join(case_path, "0/p"), "w") as f: f.write(content)

def write_k_file(case_path, wind_speed, z0_value):
    header = write_header("volScalarField", "k")
    
    kappa = 0.41
    Cmu = 0.09
    u_star = (wind_speed * kappa) / math.log(10.0 / z0_value)
    
    # k is theoretically constant in the ABL: k = uStar^2 / sqrt(Cmu)
    k_value = (u_star**2) / math.sqrt(Cmu)

    content = f"""{header}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k_value};

boundaryField
{{
    inlet
    {{
        type            atmBoundaryLayerInletK;
        z0              uniform {z0_value};
        flowDir         (1 0 0);  
        zDir            (0 0 1);
        zGround uniform 0.0;
        Zref 10.0;
        Uref            {wind_speed};
        uStar           {u_star};
        kappa           {kappa};
        value           $internalField;
    }}
    outlet          {{ type zeroGradient; }}
    ground          {{ type kqRWallFunction; value $internalField; }}
    buildings       {{ type kqRWallFunction; value $internalField; }}
    frontAndBack    {{ type symmetry; }}
    sky             {{ type slip; }}
}}
"""
    with open(os.path.join(case_path, "0/k"), "w") as f:
        f.write(content)

def write_epsilon_file(case_path, wind_speed, max_h, z0_value):
    header = write_header("volScalarField", "epsilon")
    eps_init = (0.09**0.75 * (1.5 * (wind_speed * 0.05)**2)**1.5) / (0.41 * 10)
    content = f"""{header}
dimensions      [0 2 -3 0 0 0 0];
internalField   uniform {eps_init};
boundaryField
{{
    inlet {{ type atmBoundaryLayerInletEpsilon; zDir (0 0 1); flowDir (1 0 0); Uref {wind_speed}; Zref 10.0; z0 uniform {z0_value}; zGround uniform 0.0; value $internalField; }}
    outlet {{ type zeroGradient; }}
    ground {{ type epsilonWallFunction; value $internalField; }}
    buildings {{ type epsilonWallFunction; value $internalField; }}
    frontAndBack {{ type symmetry; }}
    sky {{ type slip; }}
}}"""
    with open(os.path.join(case_path, "0/epsilon"), "w") as f: f.write(content)

def write_omega_file(case_path, wind_speed, z0_value):
    header = write_header("volScalarField", "omega")
    #omega_init = (math.sqrt(1.5 * (wind_speed * 0.05)**2)) / (0.09 * 10)
    # Standard constants for ABL
    kappa = 0.41
    Cmu = 0.09
    
    # Calculate uStar (friction velocity) based on the log law: 
    # U = (uStar/kappa) * ln(z_ref / z0) -> assuming z_ref is 10m
    u_star = (wind_speed * kappa) / math.log(10.0 / z0_value)
    
    # Internal field initialization remains a sensible average
    omega_init = u_star / (math.sqrt(Cmu) * kappa * (10.0 + z0_value))

    content = f"""{header}
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {omega_init};
boundaryField
{{
    inlet {{ type            atmBoundaryLayerInletK;
            z0              uniform {z0_value};
            flowDir         (1 0 0); 
            zDir            (0 0 1);
            zGround uniform 0.0;
            Zref 10.0;
            Uref            {wind_speed};
            uStar           {u_star};
            kappa           {kappa};
            value           $internalField; 
            }}
    outlet {{ type zeroGradient; }}
    ground {{ type omegaWallFunction; value $internalField; }}
    buildings {{ type omegaWallFunction; value $internalField; }}
    frontAndBack {{ type symmetry; }}
    sky {{ type slip; }}
}}"""
    with open(os.path.join(case_path, "0/omega"), "w") as f: f.write(content)

def write_nut_file(case_path, z0_value):
    header = write_header("volScalarField", "nut")
    content = f"""{header}
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;
boundaryField
{{
    inlet {{ type calculated; value $internalField; }}
    outlet {{ type calculated; value $internalField; }}
    ground {{ type nutkAtmRoughWallFunction; z0 uniform {z0_value}; value $internalField; }}
    buildings {{ type nutkWallFunction; value $internalField; }}
    frontAndBack {{ type symmetry; }}
    sky {{ type slip; }}
}}"""
    with open(os.path.join(case_path, "0/nut"), "w") as f: f.write(content)

# --- CONSTANT DIRECTORY ---

def write_momentum_transport(case_path, mode):
    header = write_header("dictionary", "momentumTransport")
    model = "kOmegaSST" if mode == 'URANS' else "kEpsilon"
    content = f"""{header}
simulationType RAS;
RAS
{{
    model           {model};
    turbulence      on;
    printCoeffs     on;
}}"""
    with open(os.path.join(case_path, "constant/momentumTransport"), "w") as f: f.write(content)

def write_physical_properties(case_path, nu_value):
    header = write_header("dictionary", "physicalProperties")
    content = f"""{header}
viscosityModel   constant;
nu               [0 2 -1 0 0 0 0] {nu_value};"""
    with open(os.path.join(case_path, "constant/physicalProperties"), "w") as f: f.write(content)

# --- SYSTEM DIRECTORY ---

def write_control_dict(case_path, mode):
    header = write_header("dictionary", "controlDict")
    
    if mode == 'URANS':
        # For URANS, we usually write based on time (e.g., every 5 or 10 seconds)
        app_logic = (
            "application foamRun;\n"
            "solver incompressibleFluid;\n"
            "stopAt endTime;\n"
            "endTime 200;\n"
            "deltaT 0.05;\n"
            "writeControl adjustableRunTime;\n"
            "writeInterval 5;\n"  # Writes every 5 seconds of simulation time
            "purgeWrite 3;\n"      # Keeps only last 3 sets of results
            "adjustTimeStep yes;\n"
            "maxCo 5.0;"
        )
    else:
        # For RANS (simpleFoam), we write based on iterations
        app_logic = (
            "application simpleFoam;\n"
            "stopAt endTime;\n"
            "endTime 1000;\n"
            "deltaT 1;\n"
            "writeControl timeStep;\n"
            "writeInterval 100;\n" # Writes every 100 iterations
            "purgeWrite 0;"        # Keep all iterations for convergence checking
        )
    
    content = f"""{header}
{app_logic}

startTime       0;
writeFormat     ascii;
writePrecision  6;
runTimeModifiable true;
libs ("libatmosphericModels.so");"""

    with open(os.path.join(case_path, "system/controlDict"), "w") as f: 
        f.write(content)

def write_fv_schemes(case_path, mode):
    header = write_header("dictionary", "fvSchemes")
    ddt = "backward" if mode == 'URANS' else "steadyState"
    v_scheme = "limitedLinearV 1" if mode == 'URANS' else "linearUpwind limited"
    turb_extra = "div(phi,omega)  Gauss limitedLinear 1;" if mode == 'URANS' else ""

    content = f"""{header}
ddtSchemes       {{ default {ddt}; }}
gradSchemes      {{ default Gauss linear; limited cellLimited Gauss linear 1; }}
divSchemes
{{
    default         none;
    div(phi,U)      Gauss {v_scheme};
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,epsilon) Gauss limitedLinear 1;
    {turb_extra}
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes    {{ default Gauss linear corrected; }}
interpolationSchemes {{ default linear; }}
snGradSchemes       {{ default corrected; }}
wallDist
{{
    method meshWave;
}}
"""
    with open(os.path.join(case_path, "system/fvSchemes"), "w") as f: f.write(content)

def write_fv_solution(case_path, mode):
    header = write_header("dictionary", "fvSolution")
    if mode == 'URANS':
        algo_block = "PIMPLE { nOuterCorrectors 2; nCorrectors 2; nNonOrthogonalCorrectors 2; momentumPredictor yes; }"
        turb_target = "(U|k|omega)"
    else:
        algo_block = "SIMPLE { residualControl { p 1e-5; U 1e-6; k 1e-6; epsilon 1e-6; } nNonOrthogonalCorrectors 2; pRefCell 0; pRefValue 0; }"
        turb_target = "(U|k|epsilon)"

    content = f"""{header}
solvers
{{
    p {{ solver GAMG; tolerance 1e-6; relTol 0.05; smoother GaussSeidel; }}
    pFinal {{ $p; relTol 0; }}
    turb {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }}
    "{turb_target}" {{ $turb; }}
    "{turb_target}Final" {{ $turb; relTol 0; }}
}}
{algo_block}
relaxationFactors
{{
    fields {{ p 0.3; }}
    equations {{ U 0.7; \"{turb_target}.*\" 0.7; }}
}}"""
    with open(os.path.join(case_path, "system/fvSolution"), "w") as f: f.write(content)

def write_block_mesh(case_path, extent, x_off, y_off, max_h):
    rad_x = (extent[2] - extent[0]) / 2.0
    rad_y = (extent[3] - extent[1]) / 2.0
    x_min, x_max = -rad_x * 2, rad_x * 5
    y_min, y_max = -rad_y * 2, rad_y * 2
    z_max = max_h * 5
    header = write_header("dictionary", "blockMeshDict")
    content = f"""{header}
convertToMeters 1;
vertices ( ({x_min} {y_min} 0) ({x_max} {y_min} 0) ({x_max} {y_max} 0) ({x_min} {y_max} 0) ({x_min} {y_min} {z_max}) ({x_max} {y_min} {z_max}) ({x_max} {y_max} {z_max}) ({x_min} {y_max} {z_max}) );
blocks ( hex (0 1 2 3 4 5 6 7) (50 40 20) simpleGrading (1 1 1) );
boundary ( 
    inlet {{ type patch; faces ((0 3 7 4)); }} 
    outlet {{ type patch; faces ((1 5 6 2)); }} 
    ground {{ type wall; faces ((0 1 2 3)); }} 
    frontAndBack {{ type symmetry; faces ((0 4 5 1) (3 2 6 7)); }} 
    sky {{ type patch; faces ((4 5 6 7)); }} 
);"""
    with open(os.path.join(case_path, "system/blockMeshDict"), "w") as f: f.write(content)

def write_snappy_hex_mesh(case_path, extent, x_off, y_off, max_h):
    rad_x, rad_y = (extent[2] - extent[0]) / 2.0, (extent[3] - extent[1]) / 2.0
    safe_radius = max(rad_x, rad_y) * 1.2
    x_min, y_min = -safe_radius, -safe_radius
    x_max, y_max = safe_radius, safe_radius
    macro_z, loc_z = max_h * 1.5, max_h * 3.0
    header = write_header("dictionary", "snappyHexMeshDict")
    content = f"""{header}
#includeEtc "caseDicts/mesh/generation/snappyHexMeshDict.cfg"
castellatedMesh on; snap on; addLayers off;
geometry {{
    buildings {{ type triSurfaceMesh; file "buildings.obj"; }}
    refinementBox {{ type box; min ({x_min} {y_min} 0); max ({x_max} {y_max} {macro_z}); }}
    pedestrianZone {{ type box; min ({x_min} {y_min} 0); max ({x_max} {y_max} 2.0); }}
}};
castellatedMeshControls {{
    features ( {{ file "buildings.eMesh"; level 1; }} );
    refinementSurfaces {{ buildings {{ level (3 3); patchInfo {{ type wall; }} }} }}
    refinementRegions {{ refinementBox {{ mode inside; level 2; }} pedestrianZone {{ mode inside; level 3; }} }}
    locationInMesh (0.1 0.1 {loc_z});
    nCellsBetweenLevels 3;
}}
snapControls {{ explicitFeatureSnap true; }}
addLayersControls {{ relativeSizes true; layers {{ "buildings" {{ nSurfaceLayers 2; }} }} expansionRatio 1.2; 
                                                                                            finalLayerThickness 0.5; 
                                                                                            minThickness 0.1;}}
meshQualityControls {{ #includeEtc "caseDicts/mesh/generation/meshQualityDict.cfg" }}
mergeTolerance 1e-6;"""
    with open(os.path.join(case_path, "system/snappyHexMeshDict"), "w") as f: f.write(content)

def write_mesh_quality_dict(case_path):
    header = write_header("dictionary", "meshQualityDict")
    with open(os.path.join(case_path, "system/meshQualityDict"), "w") as f: 
        f.write(header + "\n#includeEtc \"caseDicts/mesh/generation/meshQualityDict.cfg\"")

def write_surface_features(case_path):
    header = write_header("dictionary", "surfaceFeaturesDict")
    with open(os.path.join(case_path, "system/surfaceFeaturesDict"), "w") as f:
        f.write(header + "\nsurfaces (\"buildings.obj\");\n#includeEtc \"caseDicts/surface/surfaceFeaturesDict.cfg\"")