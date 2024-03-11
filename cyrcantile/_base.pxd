from libc.math cimport M_PI, log, tan, atan, exp, sinh

cdef double radians(double degrees):
    return degrees * (M_PI / 180.0)

cdef double degrees(double radians):
    return radians * (180.0 / M_PI)


cdef struct Tile:
    int x
    int y
    int z

cdef struct LngLat:
    float lng
    float lat

cdef struct LngLatBbox:
    float west
    float south
    float east
    float north

cdef struct Bbox:
    float left
    float bottom
    float right
    float top

cdef struct XY:
    float x
    float y

    
cdef LngLat truncate_lnglat(float lng, float lat):
    if lng > 180.0:
        lng = 180
    elif lng < -180.0:
        lng = -180.0
    if lat > 90.0:
        lat = 90.0
    elif lat < -90.0:
        lat = -90.0
    return LngLat(lng,lat)

cdef XY xy(float lng, float lat, bint truncate):
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
    
cdef LngLat lngLat(float x, float y, bint truncate):
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
