# What does *geo3D* do and why?

The following section aims to explain what *geo3D* produces, the creation process, the necessary ingredients and the usefulness of the product.

```{contents}
:local:
```

## What does *geo3D* do?
The primary product is a topologically correct Level of Detail 1 (LoD1) 3D City Model<sup>*</sup>. Secondary products include an application of spatial data science and an HTML-based visualisation. 

%Our mission is to empower high school learning and community engagement by fostering effective communication and advocacy at the grassroots level.
 
<sup>* ***the goal is a [Open Geospatial Consortium (OGC)](https://www.ogc.org/standard/citygml/) standard model that conforms to the [ISO 19107 spatial schema for 3D primitives](https://www.iso.org/standard/66175.html) [connecting and planar surfaces, correct orientation of the surfaces and watertight volumes]. If the result you achieve is not; you are welcome to raise an [issue](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/issues).***


### `geo3D` products

#### LoD1 3D City Models

##### How? *extrusion!*

An extremely well documented method of producing 3D Models is through extrusion. With extrusion; 2D features are <span style="color: red; opacity: 0.8;"><em>lifted</em></span> from an existing surface creating a volumetric 3D object.  <code><b>geo3D</b></code> infers the height with which to <span style="color: red; opacity: 0.8;"><em>lifted</em></span> 2D features from OSM contributions. 

The OSM tag `building:level` is taken as a [proxy for the height of a building](https://wiki.openstreetmap.org/wiki/Key:building:levels). The calculation is simply `building:level * 2.8 + 1.3`. If a structure does not have a `building:level` tag no LoD1 model is created.
 &nbsp; &nbsp;
% <figure><center>
%  <img src="/before/extrusion_tuDelft.png" width="650" height="350">
%  	%<figcaption>Figure 1</b> - The <code>osm_LoD1_3DCityModel</code> process.  Image adapted from {cite}`ledoux2021`.</figcaption>
%</center></figure>

%![The osm_LoD1_3DCityModel process](./extrusion_tuDelft.png)
%%{width="650" height="350"}
%**Figure 1** - The `osm_LoD1_3DCityModel` process.  *Image adapted from* {cite}`ledoux2021`.

```{figure} extrusion_tuDelft.png
---
#height: 150px
name: extrusion-fig
---
The `osm_LoD1_3DCityModel` process.  *Image adapted from* {cite}`ledoux2021`.
```

% <figure><center>
%  <img src="{{site.baseurl | prepend: site.url}}extrusion_tuDelft.png" alt="alt text" width="650" height="350">
%  <!-- <figcaption>Fig.1 - <code><b>The osm_LoD1_3DCityModel</b></code> process. <span style="color:blue;opacity:0.8;"><em>--image TUDelft</em></span>.</figcaption> -->
%  <figcaption>Fig.1 - The <code><b>osm_LoD1_3DCityModel</b></code> process. <em>--image adapted from</em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
%</center></figure>

%<!--
%<p align="center">
%  <img src="{{site.baseurl | prepend: site.url}}/img/extrusion_tuDelft.png" alt="alt text" width="650" height="350">
% </p> 
%<p align="center"> 
%    Fig 1. The osm_LoD1_3DCityModel process. <span style="color:blue"><em>--image TUDelft</em></span>.
%</p> -->

Fig 1 illustrates the process where the OSM *proxy `building:level` height*  is added to the raster DEM to create a 3D topologically connected surface ~ containing 2D polygons as 3D objects.

%
%<!--### Triangulated MultiSurfaces
%
%MultiSurface outputs are the walls and rooves of buildings, along with the terrain, as a collection of connected triangles. This surface is created in the [Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) format. An accompanying [material file](https://en.wikipedia.org/wiki/Wavefront_.obj_file#Material_template_library) (.mtl) to associate objects with a respective color is [available](https://github.com/AdrianKriger/osm_LoD1_3DCityModel/blob/main/village_campus/result/osm_LoD1_3DCityModel.mtl). 
%
%<figure><center>
%  <img src="{{site.baseurl | prepend: site.url}}/img/objects_horizontal_view_multisurface_tuDelft.png" alt="alt text" width="650" height="350">
%  <figcaption>Fig.2 - illustrates a horizontal view of the 2.75D surface with the exterior of all features together. <em>--image </em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
%</center></figure>   
%-->
##### CityJSON

The LoD1 3D City Model (buildings and terrain) is encoded in the light-weight JSON-based CityJSON format. In the [CityJSON](https://www.cityjson.org/) format these are [Building](https://www.cityjson.org/specs/1.0.1/#building) [City Objects](https://www.cityjson.org/specs/1.0.1/#cityjson-object) separate from the [TINRelief](https://www.cityjson.org/specs/1.0.1/#tinrelief) (ground)
%
%![The osm_LoD1_3DCityModel process](./objects_horizontal_view_solid_tuDelft.png)
%%{width="650" height="350"}
%**Figure 2** -  solid Building CityObjects connected to the terrain. *Image adapted from* {cite}`ledoux2021`.

```{figure} objects_horizontal_view_solid_tuDelft.png
---
#height: 150px
name: profile-fig
---
Solid Building CityObjects connected to the terrain. *Image adapted from* {cite}`ledoux2021`.
```
#### Spatial Data Science

The Jupyter environment allows for extensive customisation and deep analysis through *spatial data science* {cite}`Gran2021`. [*geo3D*](https://github.com/AdrianKriger/geo3D/tree/main) estimates population and calculates Building Volume Per Capita (BVPC).

While estimating population is well-known; BVPC is a fairly recent evaluation. BVPC attempts to integrate with several SDG’s most notably SGD11: developing sustainable cities and communities. BVPC represents the cubic meter of building per person. The value quantifies the amount of living space each person has in their home. BVPC can be used to evaluate overcrowding, and as a proxy for economic and housing inequality {cite}`ghosh2020`.

#### Interactive visualisation 

Further; *geo3D* can produce a pseudo-3D HTML-based visualisation serves to facilitate engagement and understanding at a neighbourhood level. The `iframe` below illustrates building stock differentiated by colour. A school, housing, retail, healthcare and community focused facilities are easily identified while the tooltips highlight the underlying data. Additional features unique to an aoi can also be included. Here farmland, streams, recreational spaces and bus rapid transit routes have been added *- you are thus limited only through data and your imagination*. 

%**To navigate on a laptop without a mouse**:
%
%- `trackpad left-click drag-left` and `-right`;
%- `Ctrl left-click drag-up`, `-down`, `-left` and `-right` to rotate and so-on and
%- `+` next to Backspace zoom-in and `-` next to `+` zoom-out.
%
%<iframe src="_static//interactiveOnly.html" style="width=400 height=250 border: none;"></iframe>
%
%```{iframe} _static/interactiveOnly.html
%:width: 400
%:height: 250
%:frameborder: 0
%```
%
%The visualisation above employs the default [Carto Dark Matter](https://github.com/CartoDB/basemap-styles) basemap. [Pydeck](https://deckgl.readthedocs.io/en/latest/index.html) supports a number of [map_styles](https://deckgl.readthedocs.io/en/latest/deck.html) including the extensive [mapbox gallery](https://www.mapbox.com/gallery/) and [Maptiler](https://www.maptiler.com/) urls (e.g.: `https://api.maptiler.com/maps/{style}/style.json?key={your API key}`).


%<figure><center>
%  <img src="{{site.baseurl | prepend: site.url}}/img/objects_horizontal_view_solid_tuDelft.png" alt="alt text" width="650" height="350">
%  <figcaption>Fig.2 - solid Building CityObjects's connected to the terrain. <em>--image </em><cite><a href="https://github.com/tudelft3d/3dfier"> 3dfier</a></cite> by the<cite><a href="https://3d.bk.tudelft.nl/"> 3D geoinformation research group at TU Delft</a></cite>.</figcaption>
%</center></figure>

%### Interactive Visualisation with Spatial Data Science
%
%A dynamic visualisation and spatial analysis is possible. [Interactive Visualisation](https://adriankriger.github.io/geo3D/docs/interactive/) and [Spatial Data Science](https://adriankriger.github.io/geo3D/docs/spatial/) discusses this further.
%
## Is it useful?

A LoD1 City Model, while basic, offers many advantages over 2D datasets. These may be used for shadow analyses, line of sight predictions, advanced flood simulation, or more advanced quantitative evaluations such as estimating wind comfort factor and simulating noise propagation.

As *geo3D* illustrates population estimation and the calculation of BVPC are also possible.

With the coming revolution in air traffic control, to accommodate newer forms of air services (delivery drones and urban air mobility), an accurate *digital* representation of the built environment will become crucial. A 3D City Model is one component for the effective air space management of the future.

Challenges do exist. Of primary concern are errors in the source data that propagate to the generated 3D model. Care must be taken to ensure the quality of both the vector building outlines and raster DEM.