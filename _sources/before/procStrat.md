# Processing Options, Strategies and Requirements

```{contents}
:local:
```
## 1. Processing Options

<p align="center"><b>There are two processing options</b></p>


| Village | Suburb |
| :-----: | :-----: |
| village is designed for extremely focused analysis at a **neighbourhood** level.  <br /> These are areas with no more than 2 500 buildings | For **larger** areas with more <br /> than 2 500 buildings; such as suburbs, census wards or tracts; please execute suburb. <br /><br />*The suburb processing option caters for regions with no internet access.*|
|  <br /> village harvests [OSM contributions](https://www.openstreetmap.org/about) via [overpass-turbo](https://wiki.openstreetmap.org/wiki/Overpass_turbo) in [GeoJSON](https://geojson.org/) format | an <strong>osm.pbf is necessary</strong>  <br /> With more substantial volumes of data;<br />suburbs extracts the necessary building outlines from the [osm.pbf format](https://wiki.openstreetmap.org/wiki/PBF_Format) (Protocolbuffer Binary Format)|

%
<p align="center"><b>within each processing option; there are two strategies</b></p>

## 2. Processing Strategies

| osm_LoD1_3DCityModel | interactiveOnly |
| :-----: | :-----: |
| <strong>Product:</strong> An <a href="https://www.ogc.org/standard/citygml/">Open Geospatial Consortium (OGC) </a> standard LoD1 3D model adhering to the <a href="https://www.iso.org/standard/66175.html">International Standards Organisation (ISO19107) </a> scheme for 3D primitives. <br><br>These simulation-ready models provide quantitative insights with such use cases as estimating noise propagation, energy demand and wind comfort factor. | <strong>Product:</strong> A colorised pseudo-3D HTML-based visualisation. <br><br><br> The visualisation is for user interaction, navigation, and sharing, promoting community engagement and understanding |
<strong>a raster DEM is necessary</strong> | 

## 3. Data Sources

### osm.pbf

OSM Proto-buff files can be harvested from a number of sources. Either crop an area directly from OpenStreetMap with the [official tool](https://www.openstreetmap.org/export#map=3/0.70/22.15), select a predefined area [from any number of providers](https://wiki.openstreetmap.org/wiki/Planet.osm), such as [Geofabrik](https://download.geofabrik.de), or... download your own.  
e.g.: Provincial extracts for South Africa are available here: [http://download.openstreetmap.fr/extracts/africa/south_africa/](http://download.openstreetmap.fr/extracts/africa/south_africa/)

### raster DEM

For the purposes of *geo3D* a raster Digital Elevation Model (DEM) is the surface of the earth free of man-made and natural features. These are available from a number of on-line sources; most notably [OpenTopography](https://opentopography.org). Almost all government agencies responsible for spatial data; open source their elevation data. 
%
In South Africa this service this via the [Chief Directorate: National Geo-spatial Intelligence](https://ngi.dalrrd.gov.za) CD:NGI [Geospatial portal](http://www.cdngiportal.co.za/cdngiportal/).

***be aware***:  *geo3D*; is a grass-roots learning tool for fine scale analysis at a local level. The very structure and methodology of *geo3D* reinforces this key theme; neighbourhood-level, place- and inquiry-based learning where users have directly contact with, and knowledge of the, environment under investigation. Fine-scale elevation models (with a resolution less than 10-meter) come at a computational cost.