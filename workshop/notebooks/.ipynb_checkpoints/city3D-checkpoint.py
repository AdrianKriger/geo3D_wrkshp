# -*- coding: utf-8 -*-
# env/geo3D_gthbRepo02
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

import math
import json
import fiona
import copy
import requests
from typing import Optional, Any, Union

import numpy as np
import pandas as pd
#import geopandas as gpd

import shapely.geometry as sg
from shapely.geometry import Point, LineString, Polygon, MultiPolygon, LinearRing, MultiPolygon, MultiLineString, MultiPoint, shape, mapping
from shapely.ops import snap, transform

import pyproj 
from pyproj import CRS, Transformer 
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

from osgeo import gdal

from openlocationcode import openlocationcode as olc

from cjio import cityjson, geom_help

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

#def coords_to_list(coords):
#    """Recursively convert tuples in coordinates to lists for PyDeck."""
#    if isinstance(coords, (float, int)):
#        return coords
#    return [coords_to_list(c) for c in coords]

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

        #with open(output_path, "w", encoding="utf-8") as f:
        #    json.dump(geojson, f, indent=2, ensure_ascii=False)
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

#def process_geometry(geometry):
#    """Return a valid Polygon or None (skip if not area)."""
#    return _ensure_polygon(geometry)

def safe_json_value(val):
    """Convert pandas/NumPy values safely for JSON serialization."""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return None
    if isinstance(val, (np.integer, np.floating)):
        return float(val)
    return val

def process_and_write_geojson(gdf, jparams=None): #, output_file='./data/fp_j.geojson'):

    """
    Consolidates building height calculations and GeoJSON writing.
    Accepts a single GeoDataFrame as input.
    
    Args:
        gdf (gpd.GeoDataFrameLite): The input GeoDataFrameLite with building data.
        output_file (str): The path to the output GeoJSON file.
    """
    #df = gdf.set_crs(WGS84[5:])
    #src_crs = getattr(gdf, "crs", None)
    crs = gdf.crs

    # 1. Filter out rows with missing 'building:levels'
    if 'building:levels' not in gdf.columns:
        print("Warning: 'building:levels' column not found. Skipping height calculations.")
        return
        
    filtered_gdf = gdf[gdf['building:levels'].notna() & (gdf['building:levels'] != '')].copy()
    
    if filtered_gdf.empty:
        print("No buildings with valid 'building:levels' found. No GeoJSON will be created.")
        return

    # Process geometry to ensure all are polygons
    filtered_gdf['geometry'] = filtered_gdf['geometry'].apply(process_geometry)
    filtered_gdf['footprint'] = filtered_gdf['geometry'].apply(lambda g: mapping(g)["coordinates"]) 
    filtered_gdf = filtered_gdf[filtered_gdf['geometry'].notna()]

    # 2. Add new columns using vectorized operations (faster than iterrows)
    #print("Calculating building heights and processing data...")
    height_df = filtered_gdf.apply(calculate_heights, axis=1, result_type='expand')
    height_cols_to_drop = [col for col in height_df.columns if col in filtered_gdf.columns]
    if height_cols_to_drop:
        filtered_gdf = filtered_gdf.drop(columns=height_cols_to_drop)
    #filtered_gdf = pd.concat([filtered_gdf.reset_index(drop=True), height_df.reset_index(drop=True)], axis=1)
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
        'amenity', 'social_facility', 'operator', 'building_height', 'roof_height',
        'ground_height', 'bottom_bridge_height', 'bottom_roof_height',
        'plus_code', 'footprint', 'geometry'
    ]
    
    # Ensure the output columns are unique before reindexing
    final_output_cols = [c for c in output_cols if c in filtered_gdf.columns]
    
    #final_gdf = filtered_gdf.filter(items=output_cols).copy()
    #final_gdf = filtered_gdf.reindex(columns=final_output_cols, fill_value=None).copy()
    #filtered_gdf = filtered_gdf.loc[:, ~filtered_gdf.columns.duplicated()]
    final_gdf = filtered_gdf[final_output_cols].copy()

    # -- Only write GeoJSON if jparams provided
    if jparams is not None:
        features = []
        for idx, row in final_gdf.iterrows():
            geom = row["geometry"]
            if geom is None or geom.is_empty:
                continue
            feature = {
                "id": str(idx),
                "type": "Feature",
                "properties": {
                    k: safe_json_value(v)
                    for k, v in row.items()
                    if k != "geometry"
                },
                "geometry": mapping(geom)
            }
            features.append(feature)
    
        geojson_dict = {
            "type": "FeatureCollection",
            "features": features
        }
    
        # Include CRS if known
        if crs:
            try:
                crs_name = crs.to_string()
            except Exception:
                crs_name = str(crs)
            geojson_dict["crs"] = {
                "type": "name",
                "properties": {"name": crs_name}
            }

        # Write clean, valid JSON
        if jparams is not None:
            fout = jparams["osm_bldings"]
            with open(fout, "w", encoding="utf-8") as f:
                json.dump(geojson_dict, f, indent=2, ensure_ascii=False)
    
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
    #address_parts = [row.get(key) for key in address_keys if row.get(key) not in [None, ""]]
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

def calculate_heights(row, storeyheight=2.8):
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
        roof_height = round(building_height + ground_height, 2)
    
    elif building_type == 'bridge':
        #min_height = row.get('min_height', 0)
        bottom_bridge_height = round(min_height + ground_height, 2)
        # Note: 'building_height' remains the default calculation
        
    elif building_type == 'roof':
        bottom_roof_height = round(levels * storeyheight + ground_height, 2)
        roof_height = round(bottom_roof_height + 1.3, 2)
        # Note: No 'building_height' for roofs
        
    return {
        'ground_height': ground_height,
        'building_height': building_height,
        'roof_height': roof_height,
        'bottom_bridge_height': bottom_bridge_height,
        'bottom_roof_height': bottom_roof_height
    }

#def explode_shapely(df, geometry_col="geometry"):
#    """
#    Explode a pandas.DataFrame with Shapely Multi-geometries into single-part geometries.
#    
#    Parameters
#    ----------
#    df : pandas.DataFrame
#        DataFrame with a column of Shapely geometries.
#    geometry_col : str
#        Name of the column with geometries (default "geometry").
#        
#    Returns
#    -------
#    pandas.DataFrame
#        A new DataFrame where Multi-geometries have been split into single geometries,
#        preserving attributes.
#    """
#    rows = []
#    for idx, row in df.iterrows():
#        geom = row[geometry_col]
#        if geom is None:
#            continue
#        if geom.geom_type.startswith("Multi") or geom.geom_type == "GeometryCollection":
#            for part in geom.geoms:
#                new_row = row.copy()
#                new_row[geometry_col] = part
#                rows.append(new_row)
#        else:
#            rows.append(row)
#    return pd.DataFrame(rows).reset_index(drop=True)

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

def gdf_to_geojson(df):
    """
    Convert a DataFrame with a 'geometry' column to a GeoJSON dict.
    All columns except 'geometry' are included in 'properties'.
    """
    def coords_to_list(coords):
        """Recursively convert tuples to lists for JSON."""
        if isinstance(coords, (float, int)):
            return coords
        return [coords_to_list(c) for c in coords]
    
    features = []
    for idx, row in df.iterrows():
        geom = mapping(row['geometry'])
        geom["coordinates"] = coords_to_list(geom["coordinates"])
        
        # Include all columns except 'geometry' in properties
        props = row.drop(labels='geometry').to_dict()
        
        features.append({
            "type": "Feature",
            "geometry": geom,
            "properties": props
        })
    
    geojson_dict = {
        "type": "FeatureCollection",
        "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
        "features": features
    }
    
    return geojson_dict

#def shapely_to_coords(geom):
#    """
#    Convert a Shapely geometry into a pydeck-ready coordinate list.
#
#    Returns:
#        - Point: [(x, y)]
#        - LineString: [(x1, y1), (x2, y2), ...]
#        - Polygon: [ [exterior coords], [hole1 coords], ... ]
#        - MultiPolygon: [ [outer1 coords + holes], [outer2 coords + holes], ... ]
#        - MultiLineString: [line1 coords, line2 coords, ...]
#    """
#    if isinstance(geom, Point):
#        return [(geom.x, geom.y)]
#    
#    elif isinstance(geom, LineString):
#        return list(geom.coords)
#    
#    elif isinstance(geom, Polygon):
#        coords = [list(geom.exterior.coords)]
#        for hole in geom.interiors:
#            coords.append(list(hole.coords))
#        return coords
#    
#    elif isinstance(geom, MultiPolygon):
#        all_polys = []
#        for poly in geom.geoms:
#            poly_coords = [list(poly.exterior.coords)]
#            for hole in poly.interiors:
#                poly_coords.append(list(hole.coords))
#            all_polys.append(poly_coords)
#        return all_polys
#    
#    elif isinstance(geom, MultiLineString):
#        return [list(line.coords) for line in geom.geoms]
#    
#    else:
#        raise TypeError(f"Unsupported geometry type: {type(geom)}")
#
#def coords(geom):
#    if geom is None:
#        return []
#    if geom.geom_type == "Point":
#        return [(geom.x, geom.y)]
#    elif geom.geom_type == "LineString":
#        return list(geom.coords)
#    elif geom.geom_type == "Polygon":
#        return list(geom.exterior.coords)
#    elif geom.geom_type.startswith("Multi"):
#        # flatten Multi geometries
#        result = []
#        for part in geom.geoms:
#            result.extend(coords(part))
#        return result
#    else:
#        return []

#def getBldVertices(dis, gt_forward, rb):
#    """
#    retrieve vertices from building footprints ~ without duplicates 
#    - these vertices already have a z attribute
#    """  
#    all_coords = []
#    min_zbld = []
#    dps = 3
#    segs = set()
#    
#    for ids, row in dis.iterrows():
#        oring = list(row.geometry.exterior.coords)
#        
#        if row.geometry.exterior.is_ccw == False:
#            #-- to get proper orientation of the normals
#            oring.reverse()
#        
#        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
#        all_coords.extend(coords_rounded)
#        zbld = [z for x, y, z in coords_rounded]
#        min_zbld.append(min(zbld))
#        
#        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
#        
#        for interior in row.geometry.interiors:
#            iring = list(interior.coords)
#            
#            if interior.is_ccw == True:
#                #-- to get proper orientation of the normals
#                iring.reverse() 
#            
#            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
#            all_coords.extend(coords_rounded)
#            
#            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
#    
#    c = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
#    
#    ac = pd.DataFrame(all_coords, 
#                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
#        
#    return ac, c, min_zbld 

def rasterQuery2(mx, my, gt_forward, rb):
    
    px = int((mx - gt_forward[0]) / gt_forward[1])
    py = int((my - gt_forward[3]) / gt_forward[5])

    intval = rb.ReadAsArray(px, py, 1, 1)

    return intval[0][0]

##- 
#def getAOIVertices(aoi, gt_forward, rb): 
#    """
#    retrieve vertices from aoi ~ without duplicates 
#    - these vertices are assigned a z attribute
#    """   
#    aoi_coords = []
#    dps = 3
#    segs = set()
#    
#    for ids, row in aoi.iterrows():
#        oring = list(row.geometry.exterior.coords)
#        
#        if row.geometry.exterior.is_ccw == False:
#            #-- to get proper orientation of the normals
#            oring.reverse()
#        
#        coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in oring]
#        aoi_coords.extend(coords_rounded)
#        
#        segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
#        
#        for interior in row.geometry.interiors:
#            iring = list(interior.coords)
#            
#            if interior.is_ccw == True:
#                #-- to get proper orientation of the normals
#                iring.reverse() 
#            
#            coords_rounded = [(round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2)) for x, y in iring]
#            aoi_coords.extend(coords_rounded)
#            
#            segs.update({(x1, y1, x2, y2) if (x1 < x2) else (x2, y2, x1, y1) for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])})
#    
#    ca = pd.DataFrame.from_dict({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
#    
#    acoi = pd.DataFrame(aoi_coords, 
#                      columns=["x", "y", "z"]).sort_values(by="z", ascending=False).drop_duplicates(subset=["x", "y"]).reset_index(drop=True)
#    
#    return acoi, ca

def getBldVertices(dis, gt_forward, rb):
    """
    Retrieve vertices from building footprints without duplicates.
    Vertices already have a z attribute.
    Works without GeoPandas.
    """
    all_coords = []
    min_zbld = []
    dps = 3
    segs = set()

    for geom in dis['geometry'].values:
        # Exterior
        oring = list(geom.exterior.coords)
        if not geom.exterior.is_ccw:
            oring.reverse()

        coords_rounded = [
            (round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2))
            for x, y in oring
        ]
        all_coords.extend(coords_rounded)
        zbld = [z for x, y, z in coords_rounded]
        min_zbld.append(min(zbld))

        segs.update({
            (x1, y1, x2, y2) if x1 < x2 else (x2, y2, x1, y1)
            for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])
        })

        # Interiors
        for interior in geom.interiors:
            iring = list(interior.coords)
            if interior.is_ccw:
                iring.reverse()
            coords_rounded = [
                (round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2))
                for x, y in iring
            ]
            all_coords.extend(coords_rounded)
            segs.update({
                (x1, y1, x2, y2) if x1 < x2 else (x2, y2, x1, y1)
                for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])
            })

    c = pd.DataFrame({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    ac = pd.DataFrame(all_coords, columns=["x", "y", "z"]).sort_values(by="z", ascending=False)\
        .drop_duplicates(subset=["x", "y"]).reset_index(drop=True)

    return ac, c, min_zbld


def getAOIVertices(aoi, gt_forward, rb):
    """
    Retrieve vertices from AOI without duplicates.
    Vertices are assigned a z attribute.
    Works without GeoPandas.
    """
    aoi_coords = []
    dps = 3
    segs = set()

    for geom in aoi['geometry'].values:
        # Exterior
        oring = list(geom.exterior.coords)
        if not geom.exterior.is_ccw:
            oring.reverse()
        coords_rounded = [
            (round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2))
            for x, y in oring
        ]
        aoi_coords.extend(coords_rounded)
        segs.update({
            (x1, y1, x2, y2) if x1 < x2 else (x2, y2, x1, y1)
            for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])
        })

        # Interiors
        for interior in geom.interiors:
            iring = list(interior.coords)
            if interior.is_ccw:
                iring.reverse()
            coords_rounded = [
                (round(x, dps), round(y, dps), round(float(rasterQuery2(x, y, gt_forward, rb)), 2))
                for x, y in iring
            ]
            aoi_coords.extend(coords_rounded)
            segs.update({
                (x1, y1, x2, y2) if x1 < x2 else (x2, y2, x1, y1)
                for (x1, y1, z1), (x2, y2, z2) in zip(coords_rounded[:-1], coords_rounded[1:])
            })

    ca = pd.DataFrame({"coords": list(segs)}).groupby("coords").size().reset_index(name="count")
    acoi = pd.DataFrame(aoi_coords, columns=["x", "y", "z"]).sort_values(by="z", ascending=False)\
        .drop_duplicates(subset=["x", "y"]).reset_index(drop=True)

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
#def doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result, crs): 
def doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, crs): 
    
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
    #"referenceSystem": jparams['cjsn_referenceSystem'],
    "referenceSystem": f"https://www.opengis.net/def/crs/EPSG/0/{crs}",
    "geographicalExtent": [
        extent[0],
        extent[1],
        minz ,
        extent[2],
        extent[3],
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

        #poly = list(result[lsattributes[i]['osm_id']].values())
        footprint = geom
        footprint = sg.polygon.orient(footprint, 1)

        #-- one building
        oneb = {}
        oneb['type'] = 'Building'
        oneb['attributes'] = {}
        for (k, v) in list(lsattributes[i].items()):
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
            #edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            #extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_bridge_height'], 
            #              allsurfaces, cm, edges)
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_bridge_height'], allsurfaces, cm)
            count = count + 1

        if lsattributes[i]['building'] == 'roof':
            #edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            #extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_roof_height'], 
            #              allsurfaces, cm, edges)
            extrude_walls(oring, lsattributes[i]['roof_height'], lsattributes[i]['bottom_roof_height'], allsurfaces, cm)
            count = count + 1

        if lsattributes[i]['building'] != 'bridge' and lsattributes[i]['building'] != 'roof':
            #new_edges = [[ele for ele in sub if ele <= lsattributes[i]['roof_height']] for sub in poly]
            #new_edges = [[min_zbld[i-count]] + sub_list for sub_list in new_edges]
            #extrude_walls(oring, lsattributes[i]['roof_height'], min_zbld[i-count], 
            #              allsurfaces, cm, new_edges)
            extrude_walls(oring, lsattributes[i]['roof_height'], min_zbld[i-count], allsurfaces, cm)
       
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
            #extrude_int_walls(iring, lsattributes[i]['roof_height'], min_zbld[i-count], allsurfaces, cm)
            extrude_walls(iring, lsattributes[i]['roof_height'], min_zbld[i-count], allsurfaces, cm)

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
        
#- new
def extrude_roof_ground(orng, irngs, height, reverse, allsurfaces, cm):
    oring = copy.deepcopy(orng)
    irings = copy.deepcopy(irngs)
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

def extrude_walls(ring, height, ground, allsurfaces, cm):
    #-- each edge become a wall, ie a rectangle
    for (j, v) in enumerate(ring[:-1]):
        l = []
        cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), ground])
        #values.append(0)
        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), height])
        cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   height])
        t = len(cm['vertices'])
        allsurfaces.append([[t-4, t-3, t-2, t-1]])    
        #allsurfaces.append([t-4, t-3, t-2, t-1])    

    #-- last-first edge
    l = []
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), ground])
    #values.append(0)
    cm['vertices'].append([round(ring[0][0], dps),  round(ring[0][1], dps),  ground])
    cm['vertices'].append([round(ring[0][0], dps),  round(ring[0][1], dps),  height])
    #values.append(0)
    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), height])
    t = len(cm['vertices'])
    allsurfaces.append([[t-4, t-3, t-2, t-1]])
    #allsurfaces.append([t-4, t-3, t-2, t-1])

#def extrude_roof_ground(orng, irngs, height, reverse, allsurfaces, cm):
#    oring = copy.deepcopy(orng)
#    irings = copy.deepcopy(irngs)
#    #irings2 = []
#    if reverse == True:
#        oring.reverse()
#        for each in irings:
#            each.reverse()
#    for (i, pt) in enumerate(oring):
#        cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
#        oring[i] = (len(cm['vertices']) - 1)
#    for (i, iring) in enumerate(irings):
#        for (j, pt) in enumerate(iring):
#            cm['vertices'].append([round(pt[0], dps), round(pt[1], dps), height])
#            irings[i][j] = (len(cm['vertices']) - 1)
#    output = []
#    output.append(oring)
#    for each in irings:
#        output.append(each)
#    allsurfaces.append(output)
#    
#def extrude_walls(ring, height, ground, allsurfaces, cm, edges):  
#    #-- each edge become a wall, ie a rectangle
#    for (j, v) in enumerate(ring[:-1]):
#        #- if iether the left or right vertex has more than 2 heights [grnd and roof] incident:
#        if len(edges[j]) > 2 or len(edges[j+1]) > 2:
#            cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), edges[j][0]])
#            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
#            c = 0
#            #- traverse up [grnd-roof]:
#            for i, o in enumerate(edges[j+1][1:]):
#                cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), o])
#                c = c + 1
#            #- traverse down [roof-grnd]:
#            for i in edges[j][::-1][:-1]:
#                cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), i])
#                c = c + 1
#            t = len(cm['vertices'])
#            c = c + 2
#            b = c
#            l = []
#            for i in range(c):
#                l.append(t-b)
#                b = b - 1 
#            allsurfaces.append([l])
#
#        #- if iether the left and right vertex has only 2 heights [grnd and roof] incident: 
#        if len(edges[j]) == 2 and len(edges[j+1]) == 2:
#            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][0]])
#            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][0]])
#            cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), edges[j+1][1]])
#            cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   edges[j][1]])
#            t = len(cm['vertices'])
#            allsurfaces.append([[t-4, t-3, t-2, t-1]])
#    
#    #- last edge polygon
#    if len(edges[-1]) == 2 and len(edges[0]) == 2:
#        cm['vertices'].append([round(ring[-1][0], dps),  round(ring[-1][1], dps), edges[-1][0]]) 
#        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
#        cm['vertices'].append([round(ring[0][0], dps),  round(ring[0][1], dps),  edges[0][1]])
#        cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), edges[-1][1]])
#        t = len(cm['vertices'])
#        allsurfaces.append([[t-4, t-3, t-2, t-1]])
#        
#    #- last edge polygon   
#    if len(edges[-1]) > 2 or len(edges[0]) > 2:
#        c = 0
#        cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   edges[-1][0]])
#        cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), edges[0][0]])
#        for i, o in enumerate(edges[0][1:]):
#            cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), o])
#            c = c + 1
#        for i in edges[-1][::-1][:-1]:
#            cm['vertices'].append([round(ring[-1][0], dps),   round(ring[-1][1], dps),   i])
#            c = c + 1
#        t = len(cm['vertices'])
#        c = c + 2
#        b = c
#        l = []
#        for i in range(c): 
#            l.append(t-b)
#            b = b - 1 
#        allsurfaces.append([l])
#               
#def extrude_int_walls(ring, height, ground, allsurfaces, cm):
#    #-- each edge become a wall, ie a rectangle
#    for (j, v) in enumerate(ring[:-1]):
#        #l = []
#        cm['vertices'].append([round(ring[j][0], dps),   round(ring[j][1], dps),   ground])
#        #values.append(0)
#        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), ground])
#        #values.append(0)
#        cm['vertices'].append([round(ring[j+1][0], dps), round(ring[j+1][1], dps), height])
#        cm['vertices'].append([round(ring[j][0], dps), round(ring[j][1], dps), height])
#        t = len(cm['vertices'])
#        allsurfaces.append([[t-4, t-3, t-2, t-1]])    
#    #-- last-first edge
#    #l = []
#    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), ground])
#    #values.append(0)
#    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), ground])
#    cm['vertices'].append([round(ring[0][0], dps), round(ring[0][1], dps), height])
#    #values.append(0)
#    cm['vertices'].append([round(ring[-1][0], dps), round(ring[-1][1], dps), height])
#    t = len(cm['vertices'])
#    allsurfaces.append([[t-4, t-3, t-2, t-1]])
    
#def output_cityjson(extent, minz, maxz, TerrainT, pts, jparams, min_zbld, acoi, result, crs):
def output_cityjson(extent, minz, maxz, TerrainT, pts, jparams, min_zbld, acoi, crs):

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
    #cm = doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, result, crs)    
    cm = doVcBndGeomRd(lsgeom, lsattributes, extent, minz, maxz, TerrainT, pts, acoi, jparams, min_zbld, crs)    
    json_str = json.dumps(cm)#, indent=2)
    fout = open(jparams['cjsn_out'], "w")                 
    fout.write(json_str)  
    ##- close fiona object
    c.close() 
    #clean cityjson
    cm = cityjson.load(jparams['cjsn_out'])               
    cityjson.save(cm, jparams['cjsn_solid']) 

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

    # Ensure 'type' column exists and is valid
    #if "type" not in ts.columns:
    #    df["type"] = "MultiPolygon"
    #else:
    #    df["type"] = df["type"].fillna("MultiPolygon")  # fill None / NaN
    #    df["type"] = df["type"].replace({"multipolygon": "MultiPolygon"})  # normalize lowercase    
        #df["type"] = df.apply(
        #lambda row: row["geometry"].geom_type if row["geometry"] is not None else "MultiPolygon", axis=1
        #)
    
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
    #gdf = gpd.read_file(geojson_vsimem)
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
    #gdf = gpd.read_file(geojson_vsimem)
    gdf = read_vsimem_geojson(geojson_vsimem)
    gdal.Unlink(geojson_vsimem)

    return gdf

