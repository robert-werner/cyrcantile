# cyrcantile/_base.pxd

from libc.math cimport M_PI, log, tan, atan, exp, sinh, floor, sin
from libc.stdlib cimport malloc, free

cdef struct MinMax:
    int min
    int max

cdef struct Tile:
    int x
    int y
    int z

cdef struct LngLat:
    double lng
    double lat

cdef struct LngLatBbox:
    double west
    double south
    double east
    double north

cdef struct Bbox:
    double left
    double bottom
    double right
    double top

cdef struct XY:
    double x
    double y

cdef double radians(double degrees)
cdef double degrees(double radians)
cdef MinMax minmax(int zoom)
cdef LngLat truncate_lnglat(double lng, double lat)

cdef XY _xy_c(double lng, double lat, bint truncate)
cdef Bbox _xy_bounds_c(Tile tile)
cdef LngLat _lnglat_c(double x, double y, bint truncate)
cdef LngLat _ul_c(Tile tile)
cdef bint _valid_c(Tile tile)
cdef Tile* _neighbours_c(Tile tile, int* count)
cdef Tile _parent_c(Tile tile, int zoom=*)
cdef list _children_c(Tile tile, int target_zoom)
cdef tuple _merge(object merge_set)
cdef object _simplify_c(object tiles)
cdef LngLatBbox _bounds_c(Tile tile)
cdef int _getBboxZoom_c(unsigned int cell0,
                        unsigned int cell1,
                        unsigned int cell2,
                        unsigned int cell3)
cdef unsigned int _rshift_c(unsigned int val, int n)
cdef Tile _tile_c(double lng, double lat, int zoom, bint truncate)
cdef Tile _bounding_tile_c(LngLatBbox bbox, bint truncate)
cdef str _quadkey_c(Tile tile)
cdef Tile _quadkey_to_tile_c(str qk)
cdef dict _feature_c(Tile t, object fid, dict props,
                     bint mercator, bint has_buffer,
                     double buffer, bint has_precision, int precision)
cdef void _collect_coords(object obj, list out)
cdef LngLatBbox _geojson_bounds_c(object obj)
