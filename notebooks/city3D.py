# -*- coding: utf-8 -*-
# env/geo3D_distV2
#########################
# helper functions to create LoD1 3D City Model from volunteered public data (OpenStreetMap) with elevation via a raster DEM.

# author: arkriger - 2023 - 2025
# github: https://github.com/AdrianKriger/geo3D

# script credit:
#    - building height from osm building:level: https://github.com/ualsg/hdb3d-code/blob/master/hdb2d.py - Filip Biljecki <filip@nus.edu.sg>
#    - extruder: https://github.com/cityjson/misc-example-code/blob/master/extruder/extruder.py - Hugo Ledoux <h.ledoux@tudelft.nl>

# additional thanks:
#    - cityjson community: https://github.com/cityjson
#########################

import json
import fiona
import copy

import numpy as np
import pandas as pd

import shapely.geometry as sg
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, LinearRing, shape, mapping
from shapely.ops import snap
from shapely.ops import transform

import pyproj

from openlocationcode import openlocationcode as olc

from cjio import cityjson, geom_help

dps = 3

def calc_Bldheight(data, is_geojson=True, output_file='./data/fp_j.geojson'):
    """Calculate building height and write to GeoJSON from either a GeoJSON dictionary or a GeoDataFrame."""
    
    storeyheight = 2.8  # Default storey height assumption
    footprints = {"type": "FeatureCollection", "features": []}

    # Decide whether data comes from GeoJSON or a GeoDataFrame
    iterable = data["features"] if is_geojson else data.iterrows()

    for item in iterable:
        f = {"type": "Feature", "properties": {}}
        
        # Extract row based on input type
        row = item if is_geojson else item[1]
        properties = row["properties"] if is_geojson else row.to_dict()

        # Handle 'tags' dictionary (GeoJSON or within DataFrame)
        tags = properties.get("tags", {}) if is_geojson else row.get("tags", {})
        if not isinstance(tags, dict):
            tags = {}  # Ensure tags is always a dictionary

        # Skip if 'building:levels' is missing OR its value is None/empty
        if is_geojson:
            has_building_levels = "building:levels" in tags and tags.get("building:levels") not in [None, "", "null"]
        else:
            has_building_levels = (
                ("building:levels" in row and pd.notna(row["building:levels"])) or
                ("tags" in row and isinstance(row["tags"], dict) and "building:levels" in row["tags"] and row["tags"]["building:levels"] not in [None, "", "null"])
            )
        
        if not has_building_levels:
            continue

        # Store OSM ID - Extract from properties (GeoJSON) or columns/tags (GeoDataFrame)
        f["properties"]["osm_id"] = properties.get("id") if is_geojson else row.get("osm_way_id") or row.get("osm_id") or row.get("id")

        # Harvest attributes (ensure only non-null values are stored)
        attributes = [
            'building', 'building:use', 'building:levels', 'building:flats', 'building:units',
            'beds', 'rooms', 'residential', 'amenity', 'social_facility'
        ]
        
        for key in attributes:
            value = tags.get(key) if is_geojson else tags.get(key) if key in tags else row.get(key)
            if value not in [None, "", {}, []]:  # Store only if valid
                f["properties"][key] = value

        # Extract address components (first check columns, then 'tags' dictionary)
        address_keys = [
            'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
            'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
        ]
        
        properties = row["properties"] if is_geojson else row.to_dict()
        tags = properties.get('tags', {}) if is_geojson else (
            row["tags"] if isinstance(row.get("tags"), dict) else {}  # Extract from 'tags' column if DataFrame
        )
        
        # Collect address parts and filter out None/empty values
        address_parts = [tags.get(key) for key in address_keys if tags.get(key) is not None and tags.get(key) != ""]
        
        # Only add address if valid parts exist
        if address_parts:
            f["properties"]["address"] = " ".join(address_parts)
        else:
            f["properties"]["address"] = None  # Set None only if no valid parts exist

        # Convert geometry to a valid polygon
        osm_shape = shape(row["geometry"]) if is_geojson else row["geometry"]

        if osm_shape.geom_type == 'LineString': 
            osm_shape = Polygon(osm_shape)
        elif osm_shape.geom_type == 'MultiPolygon': 
            osm_shape = osm_shape.geoms[0]  # Use first polygon

        f["geometry"] = mapping(osm_shape)
        f["properties"]["footprint"] = mapping(osm_shape)

        # Compute Plus Code
        p = osm_shape.representative_point()
        f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)

        # Compute building height
        levels = float(tags.get('building:levels', 1)) if is_geojson else \
                 float(tags.get('building:levels', 1)) if 'building:levels' in tags else \
                 float(row.get('building:levels', 1))  # Default to 1 if missing

        if tags.get('building') == 'cabin':
            f["properties"]['building_height'] = round(levels * storeyheight, 2)
        else:
            f["properties"]['building_height'] = round(levels * storeyheight + 1.3, 2)

        footprints["features"].append(f)

    # Store data as GeoJSON
    with open(output_file, 'w') as outfile:
        json.dump(footprints, outfile, indent=2)
        

def process_geometry(geometry):
    """Ensure valid polygon geometry for buildings."""
    if geometry.geom_type == 'LineString':
        return Polygon(geometry)
    elif geometry.geom_type == 'MultiPolygon': 
        return Polygon(geometry.geoms[0])
    else:
        return geometry

def extract_address(row):
    """Extract and format address components from the 'tags' dictionary in OSM data."""
    address_keys = [
        'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
        'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
    ]
    tags = row.get('tags', {}) if isinstance(row.get('tags'), dict) else {}
    address_parts = [tags.get(key) for key in address_keys if tags.get(key) is not None]
    return " ".join(address_parts) if address_parts else None


def calculate_building_heights(row, storeyheight=2.8):
    """Compute ground, building, and roof heights based on building type."""
    ground_height = round(row.get("mean", 0), 2)
    tags = row.get('tags', {}) if isinstance(row.get('tags'), dict) else {}
    
    levels = tags.get('building:levels', row.get('building:levels', 1))
    levels = float(levels) if str(levels).replace('.', '').isdigit() else 1.0
    
    building_type = tags.get('building', row.get('building'))
    
    if building_type == 'cabin':
        return {
            'ground_height': ground_height,
            'building_height': round(levels * storeyheight, 2),
            'roof_height': round(levels * storeyheight + ground_height, 2)
        }

    if building_type == 'bridge':
        min_height = tags.get('min_height', row.get('min_height'))
        min_height = float(min_height) if str(min_height).replace('.', '').isdigit() else (
            float(tags.get('building:min_level', row.get('building:min_level', 0))) * storeyheight
        )
        return {
            'ground_height': ground_height,
            'bottom_bridge_height': round(min_height + ground_height, 2),
            'building_height': round(levels * storeyheight, 2),
            'roof_height': round(levels * storeyheight + ground_height, 2)
        }
    
    elif building_type == 'roof':
        return {
            'ground_height': ground_height,
            'bottom_roof_height': round(levels * storeyheight + ground_height, 2),
            'roof_height': round(levels * storeyheight + ground_height + 1.3, 2)
        }
    
    else:
        return {
            'ground_height': ground_height,
            'building_height': round(levels * storeyheight + 1.3, 2),
            'roof_height': round(levels * storeyheight + 1.3 + ground_height, 2)
        }

def write_geojson(ts, jparams):
    """Process buildings and write results to GeoJSON."""
    storeyheight = 2.8
    footprints = {"type": "FeatureCollection", "features": []}
    
    for _, row in ts.iterrows():
        if row.geometry.geom_type == 'LineString' and len(row.geometry.coords) < 3:
            continue  # Skip invalid geometries
        
        if row.get('type') == 'node' or row.get('tags') is None or 'building:levels' not in row.get('tags'):
            continue  # Skip rows that don't meet the condition
        
        f = {"type": "Feature", "properties": {}}
        
        #- harvest OSM 'id'
        f["properties"]["osm_id"] = (
            row.get("osm_way_id") if row.get("osm_way_id") is not None and pd.notna(row.get("osm_way_id")) 
            else row.get("osm_id") if row.get("osm_id") is not None and pd.notna(row.get("osm_id")) 
            else row.get("id")
        )
        
        # Extract address
        f["properties"]["address"] = extract_address(row)
        
        # Extract building type
        f["properties"]["building"] = row.get("building")
        
        # Harvest attributes only if they exist
        for key in [
            'building:use', 'building:levels', 'building:flats', 'building:units',
            'beds', 'rooms', 'residential', 'amenity', 'social_facility'
        ]:
            value = row.get('tags', {}).get(key)
            if value is not None:
                f["properties"][key] = value
        
        # Process geometry
        osm_shape = process_geometry(row["geometry"])
        f["geometry"] = mapping(osm_shape)
        f["properties"]["footprint"] = mapping(osm_shape)
        
        # Compute plus_code
        p = osm_shape.representative_point()
        f["properties"]["plus_code"] = olc.encode(p.y, p.x, 11)
        
        # Compute height attributes
        height_attributes = calculate_building_heights(row, storeyheight)
        for key, value in height_attributes.items():
            if key not in f["properties"]:
                f["properties"][key] = value
        
        footprints['features'].append(f)
    
    # Store the data as GeoJSON
    with open(jparams['osm_bldings'], 'w') as outfile:
        json.dump(footprints, outfile, indent=2)
        

def getBldVertices(dis, gt_forward, rb):
    """
    retrieve vertices from building footprints ~ without duplicates 
    - these vertices already have a z attribute
    """  
    all_coords = []
    min_zbld = []
    dps = 3
    segs = set()
    
    for ids, row in dis.iterrows():
        oring = list(row.geometry.exterior.coords)
        
        if row.geometry.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
        all_coords.extend(coords_rounded)
        zbld = [z for x, y, z in coords_rounded]
        min_zbld.append(min(zbld))
        
        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
        
        for interior in row.geometry.interiors:
            iring = list(interior.coords)
            
            if interior.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
            all_coords.extend(coords_rounded)
            
            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
    
    c = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    
    ac = pd.DataFrame(all_coords, 
                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
        
    return ac, c, min_zbld 

def rasterQuery2(mx, my, gt_forward, rb):
    
    px = int((mx - gt_forward[0]) / gt_forward[1])
    py = int((my - gt_forward[3]) / gt_forward[5])

    intval = rb.ReadAsArray(px, py, 1, 1)

    return intval[0][0]

##- 
def getAOIVertices(aoi, gt_forward, rb): 
    """
    retrieve vertices from aoi ~ without duplicates 
    - these vertices are assigned a z attribute
    """   
    aoi_coords = []
    dps = 3
    segs = set()
    
    for ids, row in aoi.iterrows():
        oring = list(row.geometry.exterior.coords)
        
        if row.geometry.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
        aoi_coords.extend(coords_rounded)
        
        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
        
        for interior in row.geometry.interiors:
            iring = list(interior.coords)
            
            if interior.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
            aoi_coords.extend(coords_rounded)
            
            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
    
    ca = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    
    acoi = pd.DataFrame(aoi_coords, 
                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
    
    return acoi, ca



def concatCoords(gdf, ac):
    df2 = pd.concat([gdf, ac])
    
    return df2

def createSgmts(ac, c, gdf, idx):
    """
    create a segment list for Triangle
    - indices of vertices [from, to]
    """
    
    l = len(gdf) #- 1
    lr = 0
    idx01 = []
    
    for i, row in c.iterrows():
        frx, fry = row.coords[0], row.coords[1]
        tox, toy = row.coords[2], row.coords[3]

        [index_f] = (ac[(ac['x'] == frx) & (ac['y'] == fry)].index.values)
        [index_t] = (ac[(ac['x'] == tox) & (ac['y'] == toy)].index.values)
        idx.append([l + index_f, l + index_t])
        idx01.append([lr + index_f, lr + index_t])
    
    return idx, idx01


# # -- create CityJSON
def doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result): 
    
    #-- create the JSON data structure for the City Model
    cm = {}
    cm["type"] = "CityJSON"
    cm["version"] = "1.1"
    #cm["transform"] = {
        #"scale": [0.0, 0.0, 0.0],
        #"translate": [1.0, 1.0, 1.0]
    #},
    cm["CityObjects"] = {}
    cm["vertices"] = []
    #-- Metadata is added manually
    cm["metadata"] = {
    "title": jparams['cjsn_title'],
    "referenceDate": jparams['cjsn_referenceDate'],
    #"dataSource": jparams['cjsn_source'],
    #"geographicLocation": jparams['cjsn_Locatn'],
    "referenceSystem": jparams['cjsn_referenceSystem'],
    "geographicalExtent": [
        extent[0],
        extent[1],
        minz ,
        extent[1],
        extent[1],
        maxz
      ],
    "datasetPointOfContact": {
        "contactName": jparams['cjsn_contactName'],
        "emailAddress": jparams['cjsn_emailAddress'],
        "contactType": jparams['cjsn_contactType'],
        "website": jparams['cjsn_website']
        },
    "+metadata-extended": {
        "lineage":
            [{"featureIDs": ["TINRelief"],
              "source": [
                  {
                      "description": jparams['cjsn_+meta-description'],
                      "sourceSpatialResolution": jparams['cjsn_+meta-sourceSpatialResolution'],
                      "sourceReferenceSystem": jparams['cjsn_+meta-sourceReferenceSystem'],
                      "sourceCitation":jparams['cjsn_+meta-sourceCitation'],
                      }],
              "processStep": {
                  "description" : "Processing of raster DEM using osm_LoD1_3DCityModel workflow",
                  "processor": {
                      "contactName": jparams['cjsn_contactName'],
                      "contactType": jparams['cjsn_contactType'],
                      "website": jparams['cjsn_website']
                      }
                  }
            },
            {"featureIDs": ["Building"],
              "source": [
                  {
                      "description": "OpenStreetMap contributors",
                      "sourceReferenceSystem": "urn:ogc:def:crs:EPSG:4326",
                      "sourceCitation": "https://www.openstreetmap.org",
                  }],
              "processStep": {
                  "description" : "Processing of building vector contributions using osm_LoD1_3DCityModel workflow",
                  "processor": {
                      "contactName": jparams['cjsn_contactName'],
                      "contactType": jparams['cjsn_contactType'],
                      "website": "https://github.com/AdrianKriger/osm_LoD1_3DCityModel"
                      }
                  }
            }]
        }
    #"metadataStandard": jparams['metaStan'],
    #"metadataStandardVersion": jparams['metaStanV']
    }
      ##-- do terrain
    add_terrain_v(pts, cm)
    grd = {}
    grd['type'] = 'TINRelief'
    grd['geometry'] = [] #-- a cityobject can have >1 
      #-- the geometry
    g = {} 
    g['type'] = 'CompositeSurface'
    g['lod'] = 1
    g['boundaries'] = []
    allsurfaces = [] #-- list of surfaces
    add_terrain_b(TerrainT, allsurfaces)
    g['boundaries'] = allsurfaces
      #-- add the geom 
    grd['geometry'].append(g)
      #-- insert the terrain as one new city object
    cm['CityObjects']['terrain01'] = grd
    

    count = 0
      #-- then buildings
    for (i, geom) in enumerate(lsgeom):

        poly = list(result[lsattributes[i]['osm_id']].values())
        footprint = geom
        footprint = sg.polygon.orient(footprint, 1)

        #-- one building
        oneb = {}
        oneb['type'] = 'Building'
        oneb['attributes'] = {}
        for (k, v)in list(lsattributes[i].items()):
            if v is None:
                del lsattributes[i][k]
        for a in lsattributes[i]:
            oneb['attributes'][a] = lsattributes[i][a]                   
        oneb['geometry'] = [] #-- a cityobject can have > 1
        
        #-- the geometry
        g = {} 
        g['type'] = 'Solid'
        g['lod'] = 1
        allsurfaces = [] #-- list of surfaces forming the shell of the solid
        #-- exterior ring of each footprint
        oring = list(footprint.exterior.coords)
        oring.pop() #-- remove last point since first==last
        
        if footprint.exterior.is_ccw == False:
            #-- to get proper orientation of the normals
            oring.reverse()
        
        if lsattributes[i]['building'] == 'bridge':
            edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_bridge_height'], 
                          allsurfaces, cm, edges)
            count = count + 1

        if lsattributes[i]['building'] == 'roof':
            edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_roof_height'], 
                          allsurfaces, cm, edges)
            count = count + 1

        if lsattributes[i]['building'] != 'bridge' and lsattributes[i]['building'] != 'roof':
            new_edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            new_edges = [[min_zbld[i-count]] + sub_list for sub_list in new_edges]
            extrude_walls(oring, lsattributes[i]['roof_height'], min_zbld[i-count], 
                          allsurfaces, cm, new_edges)
       
        #-- interior rings of each footprint
        irings = []
        interiors = list(footprint.interiors)
        for each in interiors:
            iring = list(each.coords)
            iring.pop() #-- remove last point since first==last
            
            if each.is_ccw == True:
                #-- to get proper orientation of the normals
                iring.reverse() 
            
            irings.append(iring)
            extrude_int_walls(iring, lsattributes[i]['roof_height'], min_zbld[i-count], allsurfaces, cm)
            
        #-- top-bottom surfaces
        if lsattributes[i]['building'] == 'bridge':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                                False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, lsattributes[i]['bottom_bridge_height'], 
                                True, allsurfaces, cm)
        if lsattributes[i]['building'] == 'roof':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                                False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, lsattributes[i]['bottom_roof_height'], 
                                True, allsurfaces, cm)
        if lsattributes[i]['building'] != 'bridge' and lsattributes[i]['building'] != 'roof':
            extrude_roof_ground(oring, irings, lsattributes[i]['roof_height'], 
                            False, allsurfaces, cm)
            extrude_roof_ground(oring, irings, min_zbld[i-count], True, allsurfaces, cm)

        #-- add the extruded geometry to the geometry
        g['boundaries'] = []
        g['boundaries'].append(allsurfaces)
        
        #-- add the geom to the building 
        oneb['geometry'].append(g)
        #-- insert the building as one new city object
        cm['CityObjects'][lsattributes[i]['osm_id']] = oneb

    return cm

def add_terrain_v(pts, cm):
    for p in pts:
        cm['vertices'].append([p[0], p[1], p[2]])
    
def add_terrain_b(Terr, allsurfaces):
    for i in Terr:
        allsurfaces.append([[i[0], i[1], i[2]]]) 
        
def extrude_roof_ground(orng, irngs, height, reverse, allsurfaces, cm):
    oring = copy.deepcopy(orng)
    irings = copy.deepcopy(irngs)
    #irings2 = []
    if reverse == True:
        oring.reverse()
        for each in irings:
            each.reverse()
    for (i, pt) in enumerate(oring):
        cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
        oring[i] = (len(cm['vertices']) - 1)
    for (i, iring) in enumerate(irings):
        for (j, pt) in enumerate(iring):
            cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
            irings[i][j] = (len(cm['vertices']) - 1)
    output = []
    output.append(oring)
    for each in irings:
        output.append(each)
    allsurfaces.append(output)
    
def extrude_walls(ring, height, ground, allsurfaces, cm, edges):  
    #-- each edge become a wall, ie a rectangle
    for (j, v) in enumerate(ring[:-1]):
        #- if iether the left or right vertex has more than 2 heights [grnd and roof] incident:
        if len(edges[j]) > 2 or len(edges[j+1]) > 2:
            cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), edges[j][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
            c = 0
            #- traverse up [grnd-roof]:
            for i, o in enumerate(edges[j+1][1:]):
                cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), o])
                c = c + 1
            #- traverse down [roof-grnd]:
            for i in edges[j][::-1][:-1]:
                cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), i])
                c = c + 1
            t = len(cm['vertices'])
            c = c + 2
            b = c
            l = []
            for i in range(c):
                l.append(t-b)
                b = b - 1 
            allsurfaces.append([l])

        #- if iether the left and right vertex has only 2 heights [grnd and roof] incident: 
        if len(edges[j]) == 2 and len(edges[j+1]) == 2:
            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][1]])
            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][1]])
            t = len(cm['vertices'])
            allsurfaces.append([[t-4, t-3, t-2, t-1]])
    
    #- last edge polygon
    if len(edges[-1]) == 2 and len(edges[0]) == 2:
        cm['vertices'].append([round(ring[-1][0], dps),  round(ring[-1][1], dps), edges[-1][0]]) 
        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
        cm['vertices'].append([round(ring[0][0], dps),  round(ring[0][1], dps),  edges[0][1]])
        cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), edges[-1][1]])
        t = len(cm['vertices'])
        allsurfaces.append([[t-4, t-3, t-2, t-1]])
        
    #- last edge polygon   
    if len(edges[-1]) > 2 or len(edges[0]) > 2:
        c = 0
        cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   edges[-1][0]])
        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
        for i, o in enumerate(edges[0][1:]):
            cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), o])
            c = c + 1
        for i in edges[-1][::-1][:-1]:
            cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   i])
            c = c + 1
        t = len(cm['vertices'])
        c = c + 2
        b = c
        l = []
        for i in range(c): 
            l.append(t-b)
            b = b - 1 
        allsurfaces.append([l])
               
def extrude_int_walls(ring, height, ground, allsurfaces, cm):
    #-- each edge become a wall, ie a rectangle
    for (j, v) in enumerate(ring[:-1]):
        #l = []
        cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), height])
        cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), height])
        t = len(cm['vertices'])
        allsurfaces.append([[t-4, t-3, t-2, t-1]])    
    #-- last-first edge
    #l = []
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), ground])
    #values.append(0)
    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), ground])
    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), height])
    #values.append(0)
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), height])
    t = len(cm['vertices'])
    allsurfaces.append([[t-4, t-3, t-2, t-1]])
    
def output_cityjson(extent, minz, maxz, TerrainT, pts, jparams, min_zbld, acoi, result):
    """
    basic function to produce LoD1 City Model
    - buildings and terrain
    """
      ##- open buildings ---fiona object
    c = fiona.open(jparams['osm_bldings'])
    lsgeom = [] #-- list of the geometries
    lsattributes = [] #-- list of the attributes
    for each in c:
        lsgeom.append(shape(each['geometry'])) #-- geom are casted to Fiona's 
        lsattributes.append(each['properties'])
               
    #- 3D Model
    cm = doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result)    
    
    json_str = json.dumps(cm)#, indent=2)
    fout = open(jparams['cjsn_out'], "w")                 
    fout.write(json_str)  
    ##- close fiona object
    c.close() 
    #clean cityjson
    cm = cityjson.load(jparams['cjsn_out'])               
    cityjson.save(cm, jparams['cjsn_solid']) 