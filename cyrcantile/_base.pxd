from libc.math cimport M_PI, log, tan, atan, exp, sinh

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

cdef double radians(double degrees):
    return degrees * (M_PI / 180.0)

cdef double degrees(double radians):
    return radians * (180.0 / M_PI)

cdef minmax(int zoom):
    return MinMax(0, pow(2, zoom - 1))
    
cdef LngLat truncate_lnglat(double lng, double lat):
    if lng > 180.0:
        lng = 180
    elif lng < -180.0:
        lng = -180.0
    if lat > 90.0:
        lat = 90.0
    elif lat < -90.0:
        lat = -90.0
    return LngLat(lng,lat)

cdef XY xy(double lng, double lat, bint truncate):
    if truncate:
        lngLat = truncate_lnglat(lng, lat)
        lng, lat = lngLat.lng, lngLat.lat
    cdef double x = 6378137.0 * radians(lng)
    if lat <= -90:
        y = float("-inf")
    elif lat >= 90:
        y = float("inf")
    else:
        y = 6378137.0 * log(tan((M_PI * 0.25) + (0.5 * radians(lat))))
    return XY(x,y)
    
cdef LngLat lngLat(double x, double y, bint truncate):
    lng, lat = (
        x * (180 / M_PI) / 6378137.0,
        ((M_PI * 0.5) - 2.0 * atan(exp(-y / 6378137.0))) * (180 / M_PI),
    )
    if truncate:
        lngLat = truncate_lnglat(lng, lat)
        lng, lat = lngLat.lng, lngLat.lat
    return LngLat(lng, lat)

cdef LngLat ul(Tile tile):
    xtile, ytile, zoom = tile.x, tile.y, tile.z
    Z2 = pow(2, zoom)
    lon_deg = xtile / Z2 * 360.0 - 180.0
    lat_rad = atan(sinh(M_PI * (1-2 * ytile / Z2)))
    lat_deg = degrees(lat_rad)
    return LngLat(lon_deg, lat_deg)

cdef LngLatBbox bounds(Tile tile):
    xtile, ytile, zoom = tile.x, tile.y, tile.z
    Z2 = pow(2, zoom)
    ul_lon_deg = xtile / Z2 * 360.0 - 180.0
    ul_lat_rad = atan(sinh(M_PI * (1 - 2 * ytile / Z2)))
    ul_lat_deg = degrees(ul_lat_rad)

    lr_lon_deg = (xtile + 1) / Z2 * 360.0 - 180.0
    lr_lat_rad = atan(sinh(M_PI * (1 - 2 * (ytile + 1) / Z2)))
    lr_lat_deg = degrees(lr_lat_rad)

    return LngLatBbox(ul_lon_deg, lr_lat_deg, lr_lon_deg, ul_lat_deg)
