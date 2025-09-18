# vector_tileserver02.py
#########################
# helper functions to render .mbtiles as a basemap inside pydeck

# author: arkriger - 2023 - 2025
# github: https://github.com/AdrianKriger/geo3D

#########################

from flask import Flask, Response, abort
import sqlite3
import os

app = Flask(__name__)

# Path to your MBTiles file
AREA_NAME = None
MBTILES_PATH = None

def configure(area_name, port=8000):
    global AREA_NAME, MBTILES_PATH, PORT
    AREA_NAME = area_name
    MBTILES_PATH = os.path.join(os.path.dirname(__file__), "data", f"{AREA_NAME}.mbtiles")
    PORT = port

#def run(port=8000):
    #app.run(port=port)#, host="0.0.0.0")

def get_connection():
    return sqlite3.connect(MBTILES_PATH)

@app.route("/<int:z>/<int:x>/<int:y>.pbf")
def get_tile(z, x, y):
    conn = get_connection()
    cur = conn.cursor()
    tms_y = (1 << z) - 1 - y  # convert XYZ to TMS for MBTiles
    cur.execute(
        "SELECT tile_data FROM tiles WHERE zoom_level=? AND tile_column=? AND tile_row=?",
        (z, x, tms_y)
    )
    row = cur.fetchone()
    conn.close()
    if row is None:
        abort(404)
    return Response(
        row[0],
        mimetype="application/x-protobuf",
        headers={'Access-Control-Allow-Origin': '*'}
    )

@app.route("/tilejson.json")
def tilejson():
    """Serve minimal TileJSON metadata for PyDeck"""
    return {
        "tilejson": "2.2.0",
        "name": AREA_NAME,
        "scheme": "xyz",
        "tiles": [f"http://localhost:{PORT}/{{z}}/{{x}}/{{y}}.pbf"],
        "minzoom": 0,
        "maxzoom": 14
    }
