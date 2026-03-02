# -*- coding: utf-8 -*-
# env/geo3D_gthbRepo02
#########################
# helper functions to create LoD1 3D City Model from volunteered public data (OpenStreetMap) with elevation via a raster DEM.

# author: arkriger - 2023 - 2026
# github: https://github.com/AdrianKriger/geo3D

# script credit:
#    - building height from osm building:level: https://github.com/ualsg/hdb3d-code/blob/master/hdb2d.py - Filip Biljecki <filip@nus.edu.sg>
#    - extruder: https://github.com/cityjson/misc-example-code/blob/master/extruder/extruder.py - Hugo Ledoux <h.ledoux@tudelft.nl>

# additional thanks:
#    - cityjson community: https://github.com/cityjson
#########################
import os
import math
import json
import fiona
import copy
import requests
from typing import Optional, Any, Union

import numpy as np
import pandas as pd

import shapely.wkt
import shapely.geometry as sg
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, LinearRing, MultiPolygon, MultiLineString, MultiPoint, shape, mapping
from shapely.ops import snap, transform
from shapely.validation import make_valid

import pyproj 
from pyproj import CRS, Transformer 
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

from osgeo import gdal

from collections import Counter

from openlocationcode import openlocationcode as olc

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")
from cjio import cityjson, geom_help

from IPython.display import IFrame, display

dps = 3
WGS84 = "EPSG:4326"

class GeoDataFrameLite(pd.DataFrame):
    """A lightweight GeoDataFrame-like wrapper with .crs support."""

    _metadata = ["_crs"]

    @property
    def _constructor(self):
        return GeoDataFrameLite

    @property
    def crs(self) -> Optional[CRS]:
        """Return the CRS object, or None if unset."""
        return getattr(self, "_crs", None)

    @crs.setter
    def crs(self, crs_input: Any):
        """Set CRS from user input (EPSG, WKT, PROJ string, CRS object)."""
        self._crs = CRS.from_user_input(crs_input)

    def to_json(self, indent: int = None) -> str:
        """Serialize to GeoJSON FeatureCollection."""
        features = []
        for _, row in self.iterrows():
            geom = row.get("geometry")
            props = {k: v for k, v in row.items() if k != "geometry"}
            features.append({
                "type": "Feature",
                "properties": props,
                "geometry": mapping(geom) if geom is not None else None
            })
        fc = {"type": "FeatureCollection", "features": features}
        return json.dumps(fc, indent=indent)

    @classmethod
    def from_json(cls, json_input: Union[str, dict]) -> "GeoDataFrameLite":
        """Read GeoJSON string or dict into GeoDataFrameLite."""
        if isinstance(json_input, str):
            data = json.loads(json_input)
        else:
            data = json_input

        if data.get("type") != "FeatureCollection":
            raise ValueError("Expected GeoJSON FeatureCollection")

        rows = []
        for feat in data["features"]:
            geom = shape(feat["geometry"]) if feat.get("geometry") else None
            props = feat.get("properties", {})
            props["geometry"] = geom
            rows.append(props)

        df = cls(rows)
        return df

    def estimate_utm_crs(self, datum_name: str = "WGS 84") -> CRS:
        """
        Estimate the best UTM CRS for the current geometries using pyproj.database.query_utm_crs_info.
        Works like GeoPandas. Returns a pyproj CRS object.
        """
        if "geometry" not in self.columns or self.empty:
            return None

        # Compute combined bounding box
        bounds = [g.bounds for g in self["geometry"] if g is not None]
        if not bounds:
            return None
        minx = min(b[0] for b in bounds)
        miny = min(b[1] for b in bounds)
        maxx = max(b[2] for b in bounds)
        maxy = max(b[3] for b in bounds)

        # Build AreaOfInterest for pyproj query
        aoi = AreaOfInterest(
            west_lon_degree=minx,
            south_lat_degree=miny,
            east_lon_degree=maxx,
            north_lat_degree=maxy,
        )

        # Query UTM CRS info
        utm_crs_list = query_utm_crs_info(datum_name=datum_name, area_of_interest=aoi)
        if not utm_crs_list:
            raise ValueError("No suitable UTM CRS found for the bounding box.")

        # Return pyproj CRS object of the first recommended UTM
        return CRS.from_epsg(utm_crs_list[0].code)

    def to_crs(self, crs_input: Any) -> "GeoDataFrameLite":
        """
        Reproject all geometries to a new CRS.
        Returns a new GeoDataFrameLite with transformed geometries.
        """
        if "geometry" not in self.columns or self.empty:
            return self.copy()

        if self.crs is None:
            raise ValueError("Current CRS is not set. Set df.crs before calling to_crs().")

        new_crs = CRS.from_user_input(crs_input)
        transformer = Transformer.from_crs(self.crs, new_crs, always_xy=True)

        def _reproject(geom):
            if geom is None:
                return None
            return transform(transformer.transform, geom)

        df = self.copy()
        df["geometry"] = df["geometry"].apply(_reproject)
        df.crs = new_crs
        return df

def osm2gdf(data):
    """Convert Overpass JSON to list of dicts with geometry + properties."""
    nodes = {}
    ways = {}
    
    # First pass: collect nodes and ways
    for el in data:
        if el["type"] == "node":
            if "lat" in el and "lon" in el:
                nodes[el["id"]] = (el["lon"], el["lat"])
            elif "geometry" in el and len(el["geometry"]) > 0:
                pt = el["geometry"][0]
                nodes[el["id"]] = (pt["lon"], pt["lat"])
        elif el["type"] == "way":
            if "geometry" in el:
                coords = [(pt["lon"], pt["lat"]) for pt in el["geometry"]]
            else:
                coords = [nodes.get(n) for n in el.get("nodes", [])]
                if None in coords:
                    continue
            ways[el["id"]] = coords
    
    def assemble_rings(way_segments):
        """Assemble multiple way segments into one or more closed rings."""
        if not way_segments:
            return []
        
        rings = []
        remaining = [list(seg) for seg in way_segments]
        
        while remaining:
            ring = remaining.pop(0)
            made_progress = True
            while made_progress and remaining:
                made_progress = False
                for i, segment in enumerate(remaining):
                    if ring[-1] == segment[0]:
                        ring.extend(segment[1:])
                        remaining.pop(i)
                        made_progress = True
                        break
                    elif ring[-1] == segment[-1]:
                        ring.extend(reversed(segment[:-1]))
                        remaining.pop(i)
                        made_progress = True
                        break
                    elif ring[0] == segment[-1]:
                        ring = segment[:-1] + ring
                        remaining.pop(i)
                        made_progress = True
                        break
                    elif ring[0] == segment[0]:
                        ring = list(reversed(segment))[:-1] + ring
                        remaining.pop(i)
                        made_progress = True
                        break
                if ring[0] == ring[-1]:
                    break
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            if len(ring) >= 4:
                rings.append(ring)
        return rings
    
    results = []
    for el in data:
        tags = el.get("tags", {})
        geom = None
        
        # Node
        if el["type"] == "node":
            if "geometry" in el and len(el["geometry"]) > 0:
                pt = el["geometry"][0]
                geom = Point(pt["lon"], pt["lat"])
            elif "lat" in el and "lon" in el:
                geom = Point(el["lon"], el["lat"])
        
        # Way
        elif el["type"] == "way":
            coords = ways.get(el["id"])
            if not coords:
                continue
            if coords[0] == coords[-1] and len(coords) >= 4:
                geom = Polygon(coords)
            else:
                geom = LineString(coords)
        
        # Relation
        elif el["type"] == "relation":
            # Multipolygon
            if tags.get("type") == "multipolygon":
                outer_ways = []
                inner_ways = []
                for m in el.get("members", []):
                    if m["type"] == "way":
                        coords = None
                        if "geometry" in m:
                            coords = [(pt["lon"], pt["lat"]) for pt in m["geometry"]]
                        elif m.get("ref") in ways:
                            coords = ways[m.get("ref")]
                        if not coords:
                            continue
                        if m.get("role") == "outer":
                            outer_ways.append(coords)
                        elif m.get("role") == "inner":
                            inner_ways.append(coords)
                
                outer_rings = assemble_rings(outer_ways)
                outers = [Polygon(ring) for ring in outer_rings]
                inner_rings = assemble_rings(inner_ways)
                inners = [Polygon(ring) for ring in inner_rings]
                
                if outers:
                    polys = []
                    for outer in outers:
                        holes = [inner.exterior.coords for inner in inners if inner.within(outer)]
                        polys.append(Polygon(outer.exterior.coords, holes))
                    geom = MultiPolygon(polys) if len(polys) > 1 else polys[0]
            
            # Boundary / other polygon-like
            elif tags.get("type") in ["boundary", None]:
                outer_ways = []
                for m in el.get("members", []):
                    if m["type"] == "way":
                        coords = None
                        if "geometry" in m:
                            coords = [(pt["lon"], pt["lat"]) for pt in m["geometry"]]
                        elif m.get("ref") in ways:
                            coords = ways[m.get("ref")]
                        if coords:
                            outer_ways.append(coords)
                if outer_ways:
                    rings = assemble_rings(outer_ways)
                    if rings:
                        if len(rings) == 1:
                            geom = Polygon(rings[0])
                        else:
                            geom = MultiPolygon([Polygon(ring) for ring in rings])
            
            # Linear relation (e.g., bus route, train line)
            elif tags.get("type") == "route":
                line_segments = []
                for m in el.get("members", []):
                    if m["type"] == "way":
                        coords = None
                        if "geometry" in m:
                            coords = [(pt["lon"], pt["lat"]) for pt in m["geometry"]]
                        elif m.get("ref") in ways:
                            coords = ways[m.get("ref")]
                        if coords:
                            line_segments.append(coords)
                
                # Assemble lines
                assembled_lines = []
                remaining = [list(seg) for seg in line_segments]
                while remaining:
                    line = remaining.pop(0)
                    made_progress = True
                    while made_progress and remaining:
                        made_progress = False
                        for i, seg in enumerate(remaining):
                            if line[-1] == seg[0]:
                                line.extend(seg[1:])
                                remaining.pop(i)
                                made_progress = True
                                break
                            elif line[-1] == seg[-1]:
                                line.extend(reversed(seg[:-1]))
                                remaining.pop(i)
                                made_progress = True
                                break
                            elif line[0] == seg[-1]:
                                line = seg[:-1] + line
                                remaining.pop(i)
                                made_progress = True
                                break
                            elif line[0] == seg[0]:
                                line = list(reversed(seg))[:-1] + line
                                remaining.pop(i)
                                made_progress = True
                                break
                    assembled_lines.append(LineString(line))
                if assembled_lines:
                    geom = assembled_lines[0] if len(assembled_lines) == 1 else MultiLineString(assembled_lines)
        
        if geom:
            results.append({
                "geometry": geom,
                "properties": tags,
                "id": el["id"],
                "osm_type": el["type"]
            })
    
    return results

def overpass_to_gdf(query, url="https://overpass-api.de/api/interpreter", geojson=False):
    """Run an Overpass query and return GeoDataFrame or PyDeck-ready GeoJSON."""
    r = requests.get(url, params={"data": query})
    data = r.json()["elements"]

    shapes = osm2gdf(data)

    geoms = [s["geometry"] for s in shapes]
    props = [s["properties"] for s in shapes]
    osm_ids = [s["id"] for s in shapes]
    osm_types = [r["osm_type"] for r in shapes]

    #df = pd.DataFrame(props)
    df = GeoDataFrameLite(props)
    df['geometry'] = geoms
    df['osm_id'] = osm_ids
    df['osm_type'] = osm_types

    df.crs = "EPSG:4326"

    if geojson:
        features = []
        for _, row in df.iterrows():
            geom = row.get("geometry")
            if geom is None or geom.is_empty:
                continue
        
            # Guarantee all geometries are polygons or multipolygons
            if geom.geom_type not in ("Polygon", "MultiPolygon"):
                continue

            # Build feature
            feature = {
                "type": "Feature",
                "geometry": mapping(geom),
                "properties": {
                    k: (v.item() if hasattr(v, "item") else v)
                    for k, v in row.items()
                    if k != "geometry"
                }
            }
            features.append(feature)
    
        geojson = {
            "type": "FeatureCollection",
            "features": features
        }

        return geojson
    else:
        return df

def to_wgs84_point(point, src_crs=None):
    """
    Convert a Shapely point to WGS84 coordinates.
    """
    _WGS84_CRS = CRS.from_epsg(4326)
    
    if src_crs and src_crs != "EPSG:4326":

        #- create a transformer from the source CRS to WGS 84
        transformer = Transformer.from_crs(src_crs, _WGS84_CRS, always_xy=True)
        #- perform the transformation on the point's coordinates
        lon_wgs84, lat_wgs84 = transformer.transform(point.x, point.y)

        return lat_wgs84, lon_wgs84
    else: 
        return point.y, point.x

def safe_json_value(val):
    """Convert pandas/NumPy values safely for JSON serialization."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    if isinstance(val, (np.integer, np.floating)):
        return float(val)
    return val

def bldHeights(gdf): 
    """
    Consolidates building height calculations, plus code and column naming and selection.
    Accepts a single GeoDataFrame as input.
    
    Args:
        gdf (gpd.GeoDataFrameLite): The input GeoDataFrameLite with building data.
        output_file (str): The path to the output GeoJSON file.
    """
    
    def round_coords(coords, ndigits=3):
    ##- recursively round coordinate tuples inside nested lists.
    ##- works for Polygon (with holes).
        if isinstance(coords, (list, tuple)):
        # If this is a coordinate pair (x, y)
            if len(coords) == 2 and all(isinstance(c, (int, float)) for c in coords):
                return (round(coords[0], ndigits), round(coords[1], ndigits))
        # Otherwise recurse deeper
            return [round_coords(c, ndigits) for c in coords]
        return coords
    
    crs = gdf.crs

    # 1. Filter out rows with missing 'building:levels'
    if 'building:levels' not in gdf.columns:
        print("Warning: 'building:levels' column not found. Skipping height calculations.")
        return
        
    filtered_gdf = gdf[gdf['building:levels'].notna() & (gdf['building:levels'] != '')].copy()
    # Apply the mapping and then cull the precision
    if filtered_gdf.empty:
        print("No buildings with valid 'building:levels' found. No GeoJSON will be created.")
        return

    # Process geometry to ensure all are polygons
    filtered_gdf['geometry'] = filtered_gdf['geometry'].apply(process_geometry)
    #filtered_gdf['footprint'] = filtered_gdf['geometry'].apply(lambda g: mapping(g)["coordinates"])
    filtered_gdf['footprint'] = filtered_gdf['geometry'].apply(lambda g: round_coords(mapping(g)["coordinates"], 3))
    filtered_gdf = filtered_gdf[filtered_gdf['geometry'].notna()]

    # 2. Add new columns using vectorized operations (faster than iterrows)
    #print("Calculating building heights and processing data...")
    height_df = filtered_gdf.apply(process_levels, axis=1, result_type='expand')
    height_cols_to_drop = [col for col in height_df.columns if col in filtered_gdf.columns]
    if height_cols_to_drop:
        filtered_gdf = filtered_gdf.drop(columns=height_cols_to_drop)
    filtered_gdf = filtered_gdf.assign(**height_df)

    # 3. Add address and plus_code columns
    filtered_gdf['address'] = filtered_gdf.apply(extract_address, axis=1)
    
    def get_plus_code(row, crs):
        point = row.geometry.representative_point()
        lat, lon = to_wgs84_point(point, crs)
        return olc.encode(lat, lon, 11)
    
    filtered_gdf['plus_code'] = filtered_gdf.apply(get_plus_code, axis=1, args=(crs,))

    # 4. Prepare for GeoJSON export by selecting and renaming columns
    output_cols = [
        #'id', 
        'osm_id', 'address', 'building', 'building:levels', 'building:use',
        'building:flats', 'building:units', 'beds', 'rooms', 'residential',
        'amenity', 'social_facility', 'operator', 'building_height', #'roof_height',
        #'ground_height', 'bottom_bridge_height', 'bottom_roof_height',
        'min_height', 'plus_code', 'footprint', 'geometry'
    ]
    
    # Ensure the output columns are unique before reindexing
    final_output_cols = [c for c in output_cols if c in filtered_gdf.columns]
    
    final_gdf = filtered_gdf[final_output_cols].copy()
    
    if crs == None:
        final_gdf.crs = "EPSG:4326"
        
    return final_gdf

def process_geometry(geometry):
    """Ensure valid polygon geometry for buildings."""
    if geometry.geom_type == 'LineString':
        return Polygon(geometry)
    elif geometry.geom_type == 'MultiPolygon': 
        return Polygon(geometry.geoms[0])
    else:
        return geometry

def extract_address(row):
    """
    Extract and format address components from a DataFrame row.
    """
    address_keys = [
        'name', 'addr:housename', 'addr:flats', 'addr:housenumber', 'addr:street',
        'addr:suburb', 'addr:postcode', 'addr:city', 'addr:province'
    ]
    # Filter for valid, non-null values from the row
    address_parts = [
        str(row.get(key)) for key in address_keys 
        if row.get(key) not in [None, ""] and pd.notna(row.get(key))
    ]
    return " ".join(address_parts) if address_parts else None

def parse_levels(row, default=1.0):
    """
    Safely parse the 'building:levels' value from a DataFrame row.
    Handles numeric strings, floats, and missing values.
    """
    val = row.get('building:levels', default)
    try:
        if isinstance(val, str):
            # Try to convert string to float, stripping whitespace
            return float(val.strip())
        return float(val)
    except (ValueError, TypeError):
        return float(default)

def process_levels(row, storeyheight=2.8):
    """
    Compute building and roof heights based on building type and levels.
    """
    # ensure numeric conversion with defaults
    levels = pd.to_numeric(parse_levels(row), errors='coerce') or 0
    building_type = row.get('building', None)
    ground_height = pd.to_numeric(row.get('mean', 0), errors='coerce') or 0
    min_height = pd.to_numeric(row.get('min_height', 0), errors='coerce') or 0
    
    # Set default values
    building_height = round(levels * storeyheight + 1.3, 2)
    roof_height = round(building_height + ground_height, 2)
    bottom_bridge_height = None
    bottom_roof_height = None
    
    if building_type == 'cabin':
        building_height = round(levels * storeyheight, 2)
    
    elif building_type == 'bridge':
        min_height = row.get('min_height', 0)
        building_height = round(levels * 2.0, 2)
        
    elif building_type == 'roof':
        building_height = round(1.3, 3)
        min_height = round(levels * storeyheight, 2)
        
    return {
        'building_height': building_height,
        'min_height': min_height
    }
    
def gdf_to_pathlayer(df, color_col='colour'):
    paths = []
    for _, row in df.iterrows():
        geom = row['geometry']
        if geom is None:
            continue
        if geom.geom_type == 'LineString':
            paths.append({'path': list(geom.coords), 'colour': row[color_col]})
        elif geom.geom_type == 'MultiLineString':
            for line in geom.geoms:
                paths.append({'path': list(line.coords), 'colour': row[color_col]})
    return paths

def rasterQuery2(mx, my, gt_forward, rb):
    
    px = int((mx - gt_forward[0]) / gt_forward[1])
    py = int((my - gt_forward[3]) / gt_forward[5])

    intval = rb.ReadAsArray(px, py, 1, 1)

    return intval[0][0]

def prepareTri(gdf, blds, aoi):
    """
    Triangle prep
    """
    # 1. Combine all line constraints into one set
    # We include building exteriors and road boundaries
    all_polys = list(blds.geometry) + list(aoi.geometry)

    # 2. Snap to a 1cm grid to prevent floating point "near-misses"
    # This is the secret to stopping the kernel crashes.
    # We iterate through the list and use Shapely's .simplify method directly.
    # 0.01 (1cm) is great for UTM, but ensure you are in a metric CRS.
    #clean_polys = [p.simplify(0.2, preserve_topology=True) for p in all_polys]
    clean_polys = [p for p in all_polys]

    # 3. If you want to snap strictly to a grid (rounding coordinates):
    #from shapely.wkt import loads, dumps

    #def snap_geom(geom, precision=0.01):
    # This rounds the coordinates in the WKT to force a grid snap
    #    return loads(dumps(geom, rounding_precision=2))
    #snapped_polys = [snap_geom(p) for p in clean_polys]

    # 3. Use unary_union to resolve all overlaps and intersections 
    # This ensures that where a road meets a building, they share a vertex.
    merged_lines = shapely.unary_union([p.boundary for p in clean_polys])

    coords = []
    segments = []
    point_map = {}

    def get_pt_idx(pt):
        #rounded = (round(pt[0], 3), round(pt[1], 3))
        rounded = (pt[0], pt[1])
        if rounded not in point_map:
            point_map[rounded] = len(coords)
            coords.append(rounded)
        return point_map[rounded]

    # Iterate through the cleaned, merged lines
    for line in merged_lines.geoms:
        line_coords = list(line.coords)
        for i in range(len(line_coords) - 1):
            p1 = get_pt_idx(line_coords[i])
            p2 = get_pt_idx(line_coords[i+1])
            if p1 != p2:
                segments.append([p1, p2])

    # 2. NOW: Add the high-fidelity DEM points (from your GDF)
    # We add them to 'coords' so Triangle includes them in the mesh.
    # Because these points have no segments attached, they are treated as 
    # "loose points" that the triangulation must honor.
    for pt in gdf.geometry:
        # We use get_pt_idx to avoid adding a DEM point that 
        # might already exist on a road or building edge.
        get_pt_idx((pt.x, pt.y))

    # Define Regions (0: Terrain and Building)
    # We use representative points of your original 'one' (roads) and 'aoi' (terrain)
    regions = []

    # 2. ADD buildings to 'regions' with high IDs
    building_start_id = 100000
    for idx, row in blds.iterrows():
        rp = row.geometry.representative_point()
        # Region ID = 100000 + idx
        regions.append([rp.x, rp.y, building_start_id + idx, 0])

    # Add a point for the general terrain (outside roads/buildings)
    # Find a point in aoi that is NOT in any road or building
    terrain_seed = aoi.iloc[0].geometry.representative_point()
    regions.append([terrain_seed.x, terrain_seed.y, 0, 0]) # ID 0 for Terrain

    return coords, regions, segments

# # -- create CityJSON
def doVcBndGeomRd(tris, tri_attr, final_verts_3d, gdf_blds, extent, minz, maxz, jparams, gt_forward, rb, crs):
    # 1. Establish Global Offset for Precision
    x_off, y_off, z_off = float(extent[0]), float(extent[1]), float(minz)

    cm = {
    "type": "CityJSON",
    "version": "1.1",
    "transform": {
        "scale": [0.001, 0.001, 0.001], 
        "translate": [x_off, y_off, z_off]
    },
    "metadata": {
        "title": jparams['cjsn_title'],
        "referenceDate": jparams['cjsn_referenceDate'],
        # Reference system should be the EPSG integer or a standard URN
        "referenceSystem": f"urn:ogc:def:crs:EPSG::{crs}",
        # Format: [minx, miny, minz, maxx, maxy, maxz]
        "geographicalExtent": [
            float(extent[0]), float(extent[1]), float(minz), 
            float(extent[2]), float(extent[3]), float(maxz)
        ],
        "datasetPointOfContact": {
            "contactName": jparams['cjsn_contactName'],
            "emailAddress": jparams['cjsn_emailAddress'],
            "contactType": jparams['cjsn_contactType'],
            "website": jparams['cjsn_website']
        },
        "+metadata-extended": {
            "lineage": [
                {
                    "featureIDs": ["TINRelief"],
                    "source": [{
                        "description": jparams['cjsn_+meta-description'],
                        "sourceSpatialResolution": jparams['cjsn_+meta-sourceSpatialResolution'],
                        "sourceReferenceSystem": jparams['cjsn_+meta-sourceReferenceSystem'],
                        "sourceCitation": jparams['cjsn_+meta-sourceCitation'],
                    }],
                    "processStep": {
                        "description": "Processing of raster DEM using osm_LoD1_3DCityModel workflow",
                        "processor": {
                            "contactName": jparams['cjsn_contactName'],
                            "contactType": jparams['cjsn_contactType'],
                            "website": jparams['cjsn_website']
                        }
                    }
                },
                {
                    "featureIDs": ["Building"],
                    "source": [{
                        "description": "OpenStreetMap contributors",
                        "sourceReferenceSystem": "urn:ogc:def:crs:EPSG::4326",
                        "sourceCitation": "https://www.openstreetmap.org",
                    }],
                    "processStep": {
                        "description": "Processing of building vector contributions using osm_LoD1_3DCityModel workflow",
                        "processor": {
                            "contactName": jparams['cjsn_contactName'],
                            "contactType": jparams['cjsn_contactType'],
                            "website": "https://github.com/AdrianKriger/osm_LoD1_3DCityModel"
                        }
                    }
                }
            ]
        }
    },
    "CityObjects": {},
    "vertices": []
    }

    vertex_lookup = {}
    id_registry = set()

    def get_v_idx(x, y, z):
        pt = (round((float(x) - x_off) * 1000), 
              round((float(y) - y_off) * 1000), 
              round((float(z) - z_off) * 1000))
        if pt not in vertex_lookup:
            vertex_lookup[pt] = len(cm['vertices'])
            cm['vertices'].append(list(pt))
        return vertex_lookup[pt]

    def extrude_roof_ground_v2(oring, irings, height, reverse, allsurfaces):
        o_idx = [get_v_idx(p[0], p[1], height) for p in oring]
        if reverse: o_idx.reverse()
        face = [o_idx]
        for iring in irings:
            i_idx = [get_v_idx(p[0], p[1], height) for p in iring]
            if reverse: i_idx.reverse()
            face.append(i_idx)
        allsurfaces.append(face)

    def extrude_walls_v2(ring, height, ground, allsurfaces):
        for i in range(len(ring)):
            p1, p2 = ring[i], ring[(i + 1) % len(ring)]
            v1, v2 = get_v_idx(p1[0], p1[1], ground), get_v_idx(p2[0], p2[1], ground)
            v3, v4 = get_v_idx(p2[0], p2[1], height), get_v_idx(p1[0], p1[1], height)
            allsurfaces.append([[v1, v2, v3, v4]])

    def get_unique_id(base_name):
        if base_name not in id_registry:
            id_registry.add(base_name); return base_name
        suffix = 1
        while f"{base_name}_{suffix}" in id_registry: suffix += 1
        new_id = f"{base_name}_{suffix}"
        id_registry.add(new_id); return new_id

    tri_attr_int = tri_attr.astype(int).flatten()

    # --- 2. TERRAIN ---
    terrain_tris = tris[tri_attr_int == 0]
    t_boundaries = []
    for t in terrain_tris:
        v1, v2, v3 = final_verts_3d[t[0]], final_verts_3d[t[1]], final_verts_3d[t[2]]
        t_boundaries.append([[get_v_idx(*v1), get_v_idx(*v2), get_v_idx(*v3)]])
    
    cm['CityObjects']['terrain01'] = {
        'type': 'TINRelief',
        'geometry': [{'type': 'CompositeSurface', 'lod': 1, 'boundaries': t_boundaries}]
    }

    # --- 3. BUILDINGS ---
    for i, row in gdf_blds.iterrows():
        attrs = row.drop('geometry').dropna().to_dict()
        building_type = row.get('building')

        # Determine Ground Height
        if building_type in ['bridge', 'roof']:
            # FALLBACK: Use Raster Query for non-grounded structures
            centroid = row.geometry.centroid
            ground_height = float(rasterQuery2(centroid.x, centroid.y, gt_forward, rb))
        else:
            # STANDARD: Use Mesh harvesting
            # Use mesh-harvested height for standard buildings
            tag = 100000 + i
            indices = np.where(tri_attr_int == tag)[0]
            ground_height = float(np.min(final_verts_3d[np.unique(tris[indices]), 2])) if len(indices) > 0 else float(minz)

        # Pull raw inputs from GDF
        bld_h_input = pd.to_numeric(row.get('building_height', 4.0), errors='coerce')
        min_h_input = pd.to_numeric(row.get('min_height'), errors='coerce')
        
        # Logic for Attributes and Geometry Heights
        b_z, t_z = ground_height, ground_height + bld_h_input # Defaults
        bottom_bridge_height = None
        bottom_roof_height = None

        if building_type == 'bridge':
            bottom_bridge_height = round(ground_height + min_h_input, 2)
            b_z = bottom_bridge_height
            t_z = b_z + bld_h_input
            
        elif building_type == 'roof':
            bottom_roof_height = round(ground_height + min_h_input + 1.3, 2)
            b_z = bottom_roof_height
            t_z = b_z + bld_h_input # Thin surface for roofs
        
        # Prepare attributes for the City Model
        attrs.update({
            'ground_height': round(ground_height, 2),
            'building_height': round(bld_h_input, 2),
            'roof_height': round(t_z, 2),
            'bottom_bridge_height': bottom_bridge_height,
            'bottom_roof_height': bottom_roof_height
        })
        
        # ID Generation
        final_id = get_unique_id(f"osm_{attrs.get('osm_id', i)}")

        # Geometry Generation
        footprint = sg.polygon.orient(row.geometry, 1) 
        oring = list(footprint.exterior.coords)[:-1]
        irings = [list(h.coords)[:-1] for h in footprint.interiors]
        
        surfaces = []
        # Walls
        extrude_walls_v2(oring, t_z, b_z, surfaces)
        for iring in irings:
            extrude_walls_v2(iring, t_z, b_z, surfaces)

        # Top and Bottom
        extrude_roof_ground_v2(oring, irings, t_z, False, surfaces) # Top
        extrude_roof_ground_v2(oring, irings, b_z, True, surfaces)  # Bottom

        # Apply the "Strictly > 0" filter immediately
        final_attrs = {k: v for k, v in attrs.items() 
                       if (not isinstance(v, (int, float, np.number)) and v is not None and v != "") 
                       or (isinstance(v, (int, float, np.number)) and not np.isnan(v) and v > 0)}
        
        cm['CityObjects'][final_id] = {
            'type': 'Building', 
            'attributes': final_attrs,
            'geometry': [{'type': 'Solid', 'lod': 1, 'boundaries': [surfaces]}]
        }
  
    return cm

def output_cityjson(extent, minz, maxz, tris, tri_attr, final_verts_3d, dis, jparams,  gt_forward, rb, crs):

    """
    basic function to produce LoD1 City Model
    - buildings and terrain
    """
    cm = doVcBndGeomRd(tris, tri_attr, final_verts_3d, dis, extent, minz, maxz, jparams,  gt_forward, rb, crs)

    # The Robust Encoder for NumPy and Shapely types
    class CityJSONEncoder(json.JSONEncoder):
        def default(self, obj):
            # 1. Handle NumPy types
            if isinstance(obj, (np.integer, np.int64)): return int(obj)
            if isinstance(obj, (np.floating, np.float64)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            
            # 2. Handle Shapely Geometry (The fix for your error)
            # If a Point/Polygon accidentally leaked into attributes, convert to string
            if hasattr(obj, 'wkt'): 
                return str(obj.wkt)
                
            return super(CityJSONEncoder, self).default(obj)

    # Write to the output path defined in jparams
    with open(jparams['cjsn_out'], "w") as fout:
        json.dump(cm, fout, cls=CityJSONEncoder)

    cm_obj = cityjson.load(jparams['cjsn_out'])
    cityjson.save(cm_obj, jparams['cjsn_solid'])

def read_vsimem_geojson(vsimem_path):
    with fiona.open(vsimem_path, driver="GeoJSON") as src:
        #crs = src.crs
        records = []
        for feat in src:
            geom = shape(feat["geometry"])
            props = feat["properties"]
            props["geometry"] = geom
            records.append(props)
    df = GeoDataFrameLite(records)
    df.crs = "EPSG:4326"
    
    return df

def extract_boundaries_by_name(input_pbf, jparams):
    """
    Extract boundaries from an OSM PBF:
        1. Try to extract by name within boundary/place types (neighbourhood, suburb, town, etc.)
        2. If nothing is found, fallback to amenities (e.g., university, research_institute)
    Parameters:
        input_pbf (str): path to OSM PBF
        jparams (dict): must contain 'FocusArea' key for boundary name
    Returns:
        GeoDataFrame
    """
    gdal.UseExceptions()
    gdal.SetConfigOption("OGR_GEOMETRY_ACCEPT_UNCLOSED_RING", "NO")
    
    geojson_vsimem = "/vsimem/boundaries.geojson"
    boundary_name = jparams.get("FocusArea")
    place_types = ["neighbourhood", "suburb", "quarter", "borough", "village", "town", "city"]
    amenity_list = ["university", "research_institute"]

    # --- Try boundary name first, restricted to place types ---
    place_filter = " OR ".join([f"place = '{p}'" for p in place_types])
    where_filter = f"name = '{boundary_name}' AND ({place_filter})"
    gdal.VectorTranslate(
        geojson_vsimem,
        input_pbf,
        format="GeoJSON",
        layers=["multipolygons"],
        options=["-where", where_filter, "-makevalid"]
    )
    gdf = read_vsimem_geojson(geojson_vsimem)
    gdal.Unlink(geojson_vsimem)

    if len(gdf) > 0:
        return gdf

    # --- Fallback to amenities ---
    amenity_filter = " OR ".join([f"amenity = '{a}'" for a in amenity_list])
    where_filter = f"name = '{boundary_name}' AND ({amenity_filter})"
    gdal.VectorTranslate(
        geojson_vsimem,
        input_pbf,
        format="GeoJSON",
        layers=["multipolygons"],
        options=["-where", where_filter, "-makevalid"]
    )
    gdf = read_vsimem_geojson(geojson_vsimem)
    gdal.Unlink(geojson_vsimem)

    return gdf

def _harvestSolar(input_pbf, minx, miny, maxx, maxy, epsg):
    gdal.UseExceptions()
    gdal.SetConfigOption("OGR_GEOMETRY_ACCEPT_UNCLOSED_RING", "NO") 
    gdal.SetConfigOption("OGR_INTERLEAVED_READING", "YES")

    # Define the single, correct SQL filter for solar generators (rooftop renewable)
    sql_where_solar_generator = """
        other_tags LIKE '%"power"=>"generator"%' AND 
        other_tags LIKE '%"generator:source"=>"solar"%'
    """
    all_solar_dfs = [] # Use list to hold DataFrames

    # --- 1. Extract from the 'multipolygons' layer (for complex features) ---
    layer_name_poly = "multipolygons"
    geojson_vsimem_poly = f"/vsimem/solar_temp_{layer_name_poly}.geojson"
    #print(f"1. Checking layer: {layer_name_poly}...")

    gdal.VectorTranslate(
        geojson_vsimem_poly,
        input_pbf,
        format="GeoJSON",
        layers=[layer_name_poly],
        options=[
            "-where", sql_where_solar_generator,
            "-makevalid",
            "-spat", str(minx), str(miny), str(maxx), str(maxy),
        #"-nlt", "POLYGON"
        ]
    )
    gdf_multipolygons = read_vsimem_geojson(geojson_vsimem_poly)
    if not gdf_multipolygons.empty:
        all_solar_dfs.append(gdf_multipolygons)
        #print(f"   -> Found {len(gdf_multipolygons)} features in {layer_name_poly}.")

    # --- 2. Extract from the 'lines' layer (for simple closed ways/polygons) ---
    layer_name_lines = "lines"
    geojson_vsimem_lines = f"/vsimem/solar_temp_{layer_name_lines}.geojson"
    #print(f"2. Checking layer: {layer_name_lines}...")

    gdal.VectorTranslate(
        geojson_vsimem_lines,
        input_pbf,
        format="GeoJSON",
        layers=[layer_name_lines], # This processes closed Ways as Polygons
        options=[
            "-where", sql_where_solar_generator,
            "-makevalid",
            "-spat", str(minx), str(miny), str(maxx), str(maxy),
            "-nlt", "POLYGON"
        ]
    )
    gdf_lines = read_vsimem_geojson(geojson_vsimem_lines)
    if not gdf_lines.empty:
        all_solar_dfs.append(gdf_lines)
        #print(f"   -> Found {len(gdf_lines)} features in {layer_name_lines}.")

    # --- 3. Combine and Finalize ---
    if all_solar_dfs:
        # Concatenate results using pandas.concat
        gdf_solar_combined = pd.concat(all_solar_dfs, ignore_index=True)
        gdf_solar_combined.crs = "EPSG:4326"
    
        # Project the combined GeoDataFrame (assuming .to_crs is a method on the returned structure)
        gdf_solar_combined = gdf_solar_combined.to_crs(epsg)
    
        #print(f"\n✅ Total solar generator polygons harvested: {len(gdf_solar_combined)}")
    #else:
        # Create an empty DataFrame if nothing was found
        #gdf_solar_combined = pd.DataFrame() 
        #print("\033[1m No solar generator features found.\033[0m No rooftop solar are mapped in", jparams['FocusArea'])

    # Cleanup VSI Memory
    gdal.Unlink(geojson_vsimem_poly)
    gdal.Unlink(geojson_vsimem_lines)

    return gdf_solar_combined

def calculate_azimuth_from_geometry(polygon):
    """
    Calculates the azimuth (angle from North, clockwise, 0-180) 
    of the minimum rotated bounding box for a given Shapely Polygon.
    """
    if not polygon or polygon.geom_type not in ['Polygon', 'MultiPolygon']:
        return 0.0

    min_rect = polygon.minimum_rotated_rectangle
    
    if min_rect.geom_type != 'Polygon':
        return 0.0
        
    coords = np.array(min_rect.exterior.coords)
    
    segment1 = coords[1] - coords[0]
    segment2 = coords[2] - coords[1]
    
    len1 = np.linalg.norm(segment1)
    len2 = np.linalg.norm(segment2)

    long_segment = segment1 if len1 >= len2 else segment2
        
    if np.linalg.norm(long_segment) == 0:
        return 0.0

    angle_rad = np.arctan2(long_segment[1], long_segment[0])
    angle_deg = np.degrees(angle_rad)
    
    # Convert angle (from X-axis CCW) to Azimuth (from North CW)
    azimuth = 90.0 - angle_deg
    azimuth = azimuth % 360.0
        
    # Constrain to 0-180 range
    if azimuth > 180.0:
        azimuth -= 180.0

    return azimuth

def _with_solar(gdf_buildings, gdf_solar, epsg):
    """
    Efficient Dual Join: Performs both building-centric and solar-centric joins
    in one loop.

    Returns: (gdf_buildings_modified, gdf_solar_modified)
    """

    n_bld = len(gdf_buildings["geometry"])
    n_sol = len(gdf_solar["geometry"])
    
    # CRITICAL: IDs for both DataFrames
    BLD_ID_COLUMN = "osm_id"  # Assuming 'osm_id' is the unique ID in the building layer
    SOLAR_ID_COLUMN = "osm_id"
    
    # --- BUILDING-CENTRIC OUTPUT (for blds.df) ---
    solar_id_lists = [[] for _ in range(n_bld)]  # List of solar IDs for each building
    solar_m = [[] for _ in range(n_bld)]
    has_solar = [False] * n_bld

    # --- SOLAR-CENTRIC OUTPUT (for gdf_solar.df) ---
    bld_id_lists = [[] for _ in range(n_sol)]  # List of building IDs for each solar panel
    #bld_id_lists =[None] * n_sol
    
    # Brute-force intersection
    for i in range(n_bld):
        b_geom = gdf_buildings["geometry"].iloc[i]
        bld_id = gdf_buildings[BLD_ID_COLUMN].iloc[i] # Get the building ID

        for j in range(n_sol):
            s_geom = gdf_solar["geometry"].iloc[j]
            sol_id = gdf_solar[SOLAR_ID_COLUMN].iloc[j] # Get the solar ID
            s_m = gdf_solar['generator:method'].iloc[j] # Get the building ID

            if b_geom.contains(s_geom):
                # 1. Building-Centric Logic (Attaches solar ID to building)
                has_solar[i] = True
                solar_id_lists[i].append(sol_id)
                solar_m[i].append(s_m)

                # 2. Solar-Centric Logic (Attaches building ID to solar panel)
                bld_id_lists[j].append(bld_id)
                #if bld_id_lists[j] is None:
                #    bld_id_lists[j] = bld_id

    # --- Finalize Lists (Replace [] with None) ---
    for i in range(n_bld):
        if not solar_id_lists[i]:
            solar_id_lists[i] = None

    for j in range(n_sol):
        if not bld_id_lists[j]:
            bld_id_lists[j] = None

    # --- CREATE OUTPUT DataFrames ---

    # 1. Modified Building DataFrame (blds)
    #gdf_blds_out = dict(gdf_buildings)
    gdf_buildings["children"] = solar_id_lists 
    gdf_buildings["has_solar"] = has_solar
    gdf_buildings["method"] = solar_m
    #blds_out = GeoDataFrameLite(gdf_blds_out)
    #blds_out.crs = epsg

    # 2. Modified Solar DataFrame (gdf_solar)
    # We must use the solar DataFrame as the base to keep all geometry/attributes
    #gdf_solar_out = dict(gdf_solar)
    # The new column holds the list of intersecting building IDs
    gdf_solar["parent"] = bld_id_lists 
    #gdf_solar_out = GeoDataFrameLite(gdf_solar_out)
    #gdf_solar_out.crs = epsg
    gdf_solar = gdf_solar.rename(columns={'generator:method': 'method'})
    gdf_solar['area'] = gdf_solar['geometry'].apply(lambda geom: geom.area)
    gdf_solar['azimuth'] = gdf_solar['geometry'].apply(calculate_azimuth_from_geometry)

    return gdf_buildings, gdf_solar

def process_osm_geoms(vsimem_path, input_pbf, layer_name, sql_where, aoi_geom):
    """Generic helper to extract, parse tags, and clip OSM data"""
    
    # 1. Setup Buffering Logic
    # 0.01 degrees is ~1.1km. We use this for both searching and clipping for these types.
    if "highway IS NOT NULL" in sql_where or "leisure IS NOT NULL" in sql_where:
        effective_geom = aoi_geom.buffer(0.01)
    else:
        effective_geom = aoi_geom

    minx, miny, maxx, maxy = effective_geom.bounds
    
    gdal.UseExceptions()
    gdal.SetConfigOption("OGR_GEOMETRY_ACCEPT_UNCLOSED_RING", "NO")
    
    # 2. Extract from PBF
    gdal.VectorTranslate(
        vsimem_path,
        input_pbf,
        format="GeoJSON",
        layers=[layer_name],
        options=[
            "-where", sql_where,
            "-spat", str(minx), str(miny), str(maxx), str(maxy)
        ]
    )
    
    gdf = read_vsimem_geojson(vsimem_path)
    gdal.Unlink(vsimem_path)
    
    if gdf.empty:
        return gdf

    # 3. Safe Tag Parsing
    def safe_convert(tag_string):
        if isinstance(tag_string, str):
            try:
                formatted_string = "{" + tag_string.replace("=>", ":").replace("\n", " ") + "}"
                return json.loads(formatted_string)
            except:
                return {}
        return {}

    gdf["tags"] = gdf["other_tags"].apply(safe_convert)
    
    # Use reset_index to ensure the concat joins on the correct rows
    tags_df = pd.json_normalize(gdf['tags'])
    gdf = pd.concat([gdf.reset_index(drop=True), tags_df.reset_index(drop=True)], axis=1)
    
    if 'other_tags' in gdf.columns: 
        gdf = gdf.drop(columns=['other_tags'])
    if 'tags' in gdf.columns: 
        gdf = gdf.drop(columns=['tags'])

    # 4. Standardize IDs
    if 'osm_id' in gdf.columns:
        if 'osm_way_id' in gdf.columns:
            gdf['osm_id'] = [o if pd.notna(o) else w for o, w in zip(gdf['osm_id'], gdf['osm_way_id'])]
            gdf = gdf.drop(columns=['osm_way_id'])
    elif 'osm_way_id' in gdf.columns:
        gdf = gdf.rename(columns={'osm_way_id': 'osm_id'})

    # 5. Geometry Clip to the EFFECTIVE geometry (either AOI or 1km Buffer)
    gdf = gdf[gdf.geometry.apply(lambda x: x.intersects(effective_geom))]
    
    gdf.crs = "EPSG:4326"
    gdf = gdf.fillna("")
    
    return gdf

def show_interactive_html(html_path, width="100%", height=500):
    """
    Display a local MapLibre HTML file inside a Jupyter notebook.
    """
    if not os.path.exists(html_path):
        raise FileNotFoundError(html_path)

    display(
        IFrame(
            src=html_path,
            width=width,
            height=height,
            extras=['allow-scripts', 'allow-same-origin']
        )
    )

def gdf_to_geojson(obj, geom_col="geometry"):
    """
    Accepts:
    - pandas DataFrame with shapely geometries
    - GeoJSON dict
    - None

    Returns a GeoJSON FeatureCollection dict.
    """

    # Case 1: None
    if obj is None:
        return {"type": "FeatureCollection", "features": []}

    # Case 2: already GeoJSON
    if isinstance(obj, dict):
        # minimal sanity check
        if obj.get("type") == "FeatureCollection":
            return obj
        else:
            raise ValueError("Dict provided is not a GeoJSON FeatureCollection")

    # Case 3: pandas DataFrame
    if hasattr(obj, "iterrows"):
        if len(obj) == 0:
            return {"type": "FeatureCollection", "features": []}

        features = []
        for _, row in obj.iterrows():
            geom = row[geom_col]
            if geom is None or geom.is_empty:
                continue

            props = {
                k: v for k, v in row.items()
                if k != geom_col
            }

            features.append({
                "type": "Feature",
                "geometry": shapely.geometry.mapping(geom),
                "properties": props
            })

        return {
            "type": "FeatureCollection",
            "features": features
        }

    # Anything else is an error
    raise TypeError(f"Unsupported layer type: {type(obj)}")

def create_maplibre_3Dviz(
    result_dir,
    buildings_gdf,
    roads_gdf=None,
    water_gdf=None,
    green_gdf=None,
    brt_gdf=None,
    center=None,
    zoom=16,
    offline=False,
    local_js_path="../data/maplibre-gl.js",
    local_css_path="../data/maplibre-gl.css",
    remote_style="https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json"
):
    # 1. Determine map center safely
    if center is None:
        try:
            # Union all geometries to find the true center
            all_geoms = [g for g in buildings_gdf['geometry'] if g is not None]
            combined = shapely.ops.unary_union(all_geoms)
            center = [combined.centroid.x, combined.centroid.y]
        except:
            center = [0, 0] # Fallback if data is empty

    # 2. ROBUST HARVESTING & CULLING (Handles Pandas + Shapely + Empty Data)
    def harvest_and_cull(df, columns_to_keep):
        # Shield: Check if it's a valid DataFrame-like object with data
        if df is None or not hasattr(df, 'empty') or df.empty:
            return {"type": "FeatureCollection", "features": []}
        
        # Keep ONLY the columns needed to minimize file size
        existing_cols = [c for c in columns_to_keep if c in df.columns]
        # Ensure we have geometry
        if 'geometry' not in df.columns:
            return {"type": "FeatureCollection", "features": []}
            
        temp = df[existing_cols + ['geometry']].copy()
        
        # COORDINATE CULLING: Round to 5 decimal places for efficiency
        def truncate_geom(g):
            if g is None: return None
            return shapely.wkt.loads(shapely.wkt.dumps(g, rounding_precision=5))
            
        temp['geometry'] = temp['geometry'].apply(truncate_geom)
        temp = temp.fillna("")

        # Handle building height specifically
        if 'building_height' in temp.columns:
            temp['building_height'] = pd.to_numeric(temp['building_height'], errors='coerce').fillna(10)
        
        # Manual GeoJSON conversion (Stable for home-baked DFs)
        features = []
        for _, row in temp.iterrows():
            if not row['geometry']: continue
            feat = {
                "type": "Feature",
                "properties": {c: row[c] for c in existing_cols},
                "geometry": shapely.geometry.mapping(row['geometry'])
            }
            features.append(feat)
            
        return {"type": "FeatureCollection", "features": features}

    # Harvest Data 
    building_data = harvest_and_cull(buildings_gdf, ['building_height', 'fill_color', 'osm_id', 'address', 'building', 'plus_code'])
    road_data = harvest_and_cull(roads_gdf, ['highway'])
    water_data = harvest_and_cull(water_gdf, ['natural', 'waterway'])
    green_data = harvest_and_cull(green_gdf, ['leisure', 'name'])
    brt_data = harvest_and_cull(brt_gdf, ['colour'])

    # 3. Asset Logic
    if offline:
        try:
            with open(local_js_path, 'r', encoding='utf-8') as f:
                js_content = f"<script>{f.read()}</script>"
            with open(local_css_path, 'r', encoding='utf-8') as f:
                css_content = f"<style>{f.read()}</style>"
        except FileNotFoundError:
            # Fallback to CDN if local files aren't found to prevent a broken map
            js_content = '<script src="https://unpkg.com/maplibre-gl@3.6.2/dist/maplibre-gl.js"></script>'
            css_content = '<link href="https://unpkg.com/maplibre-gl@3.6.2/dist/maplibre-gl.css" rel="stylesheet" />'
        
        style_js = json.dumps({
            "version": 8, "sources": {}, 
            "layers": [{"id":"bg","type":"background","paint":{"background-color":"#0e0e0e"}}]
        })
    else:
        js_content = '<script src="https://unpkg.com/maplibre-gl@3.6.2/dist/maplibre-gl.js"></script>'
        css_content = '<link href="https://unpkg.com/maplibre-gl@3.6.2/dist/maplibre-gl.css" rel="stylesheet" />'
        style_js = f"'{remote_style}'"

    # 4. HTML Template
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8" />
    <title>geo3D. 3D City Models for Geography and Sustainable Development Education</title>
    {css_content}
    {js_content}
    <style>
        body {{ margin: 0; padding: 0; }}
        #map {{ position: absolute; top: 0; bottom: 0; width: 100%; background: #0e0e0e; }}
        .popup-content {{ font-family: sans-serif; font-size: 12px; line-height: 1.5; color: #333; }}
        hr {{ border: 0; border-top: 1px solid #eee; margin: 5px 0; }}
    </style>
</head>
<body>
<div id="map"></div>
<script>
    const map = new maplibregl.Map({{
        container: 'map',
        style: {style_js},
        center: {json.dumps(center)},
        zoom: {zoom},
        pitch: 60,
        bearing: -10,
        antialias: true
    }});

    map.on('load', () => {{
        map.addControl(new maplibregl.NavigationControl({{ visualizePitch: true }}), 'top-right');

        map.addSource('roads', {{ type: 'geojson', data: {json.dumps(road_data)} }});
        map.addSource('buildings', {{ type: 'geojson', data: {json.dumps(building_data)} }});
        map.addSource('water', {{ type: 'geojson', data: {json.dumps(water_data)} }});
        map.addSource('parks', {{ type: 'geojson', data: {json.dumps(green_data)} }});
        map.addSource('bus', {{ type: 'geojson', data: {json.dumps(brt_data)} }});

        // Layer: Water
        map.addLayer({{
            'id': 'water', 'type': 'fill', 'source': 'water',
            'paint': {{ 'fill-color': '#01579b', 'fill-opacity': 1 }}
        }});

        // Layer: Parks
        map.addLayer({{
        'id': 'parks', 'type': 'fill', 'source': 'parks',
            'paint': {{ 'fill-color': '#66bb6a', 'fill-opacity': 0.4 }}
        }});

        // Layer: Road Hierarchy
        const roadLayers = [
            {{ id: 'road-motorway', filter: ['==', 'highway', 'motorway'], color: '#666666', width: [8, 1, 14, 6] }},
            {{ id: 'road-primary', filter: ['==', 'highway', 'primary'], color: '#8b949e', width: [8, 0.75, 14, 4] }},
            {{ id: 'road-secondary', filter: ['==', 'highway', 'secondary'], color: '#6e7681', width: [9, 0.5, 14, 3] }},
            {{ id: 'road-tertiary', filter: ['==', 'highway', 'tertiary'], color: '#5a5f66', width: [10, 0.4, 14, 2.5] }},
            {{ id: 'road-residential', filter: ['==', 'highway', 'residential'], color: '#444c56', width: [11, 0.3, 16, 2] }},
            {{ id: 'road-service', filter: ['in', 'highway', 'service', 'track', 'minor', 'motorway_link'], color: '#363b42', width: [12, 0.25, 16, 1] }}
        ];

        roadLayers.forEach(layer => {{
            map.addLayer({{
                'id': layer.id, 'type': 'line', 'source': 'roads',
                'filter': layer.filter,
                'paint': {{
                    'line-color': layer.color,
                    'line-width': ['interpolate', ['linear'], ['zoom'], ...layer.width]
                }}
            }});
        }});

        // Bus Routes (Styled by the RGB column you created)
        map.addLayer({{
            id: 'bus-layer',
            type: 'line',
            source: 'bus',
            layout: {{ 'line-join': 'round', 'line-cap': 'round' }},
            paint: {{
                'line-color': [
                    'case',
                    ['has', 'colour'],
                    ['rgb', ['at', 0, ['get', 'colour']], ['at', 1, ['get', 'colour']], ['at', 2, ['get', 'colour']]],
                    '#FF4500'
                ],
                'line-width': 3
            }}
        }});

        // Layer: 3D Buildings
        map.addLayer({{
            'id': '3d-buildings', 'type': 'fill-extrusion', 'source': 'buildings',
            'paint': {{
                'fill-extrusion-color': ['case', ['has', 'fill_color'], ['get', 'fill_color'], '#4a4e5a'],
                'fill-extrusion-height': ['coalesce', ['to-number', ['get', 'building_height']], 10],
                'fill-extrusion-opacity': 0.6
            }}
        }});

        // Popup Logic with Plus Code
        map.on('click', '3d-buildings', (e) => {{
            const p = e.features[0].properties;
            const content = '<div class="popup-content">' +
                '<strong>Building:</strong> ' + (p.building || 'N/A') + '<br><hr>' +
                '<strong>Address:</strong> ' + (p.address || 'N/A') + '<br>' +
                '<strong>Height:</strong> ' + (p.building_height || 0) + 'm<br>' +
                '<strong>Plus Code:</strong> ' + (p.plus_code || 'N/A') +
                '</div>';
            new maplibregl.Popup().setLngLat(e.lngLat).setHTML(content).addTo(map);
        }});
    }});
</script>
</body>
</html>
"""
    with open(result_dir, "w", encoding="utf-8") as f:
        f.write(html_content)
    return result_dir