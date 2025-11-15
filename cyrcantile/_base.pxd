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

cdef MinMax minmax(int zoom)
cdef XY xy(double lng, double lat, bint truncate)
cdef LngLat truncate_lnglat(double lng, double lat)
cdef Bbox xy_bounds(Tile tile)
cdef double degrees(double radians)
cdef double radians(double degrees)
cdef LngLat lngLat(double x, double y, bint truncate)
cdef Tile tile(double lng, double lat, int zoom, bint truncate)
cdef LngLatBbox bounds(Tile tile)
cdef Tile* neighbours(Tile tile, int* count)
cdef bool valid(Tile tile)
cdef LngLat ul(Tile tile)
cdef Tile parent(Tile tile, int zoom=*)
cdef list children(Tile tile, int target_zoom)
cdef tuple _merge(object merge_set)
cdef object _simplify(object tiles)
cdef int _getBboxZoom(Bbox bbox)
cdef XY _xy(double lng, double lat, bint truncate)
cdef Tile bounding_tile(LngLatBbox bbox, bint truncate)
cdef str quadkey(Tile tile)
cdef Tile quadkey_to_tile(str qk)
cdef dict _feature(Tile t, object fid, dict props, bint mercator, bint has_buffer, double buffer, bint has_precision, int precision)
cdef void _collect_coords(object obj, list out)
cdef LngLatBbox geojson_bounds_c(object obj)