from libc.math cimport M_PI, log, tan, atan, exp, sinh, floor, sin

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

cdef int _getBboxZoom(Bbox bbox):
    MAX_ZOOM = 28
    for z in range(0, MAX_ZOOM):
        mask = 1 << (32 - (z + 1))
        if (bbox.left & mask) != (bbox.right & mask) or (bbox.bottom & mask) != (bbox.top & mask):
            return z
    return MAX_ZOOM


cdef rshift(double val, int n):
    return (val % 0x100000000) >> n

cdef XY _xy(double lng, double lat, bint truncate):

    if truncate:
        lng, lat = truncate_lnglat(lng, lat)

    x = lng / 360.0 + 0.5
    sinlat = sin(radians(lat))

    y = 0.5 - 0.25 * log((1.0 + sinlat) / (1.0 - sinlat)) / M_PI
    
    return XY(x,y)

cdef Tile tile(double lng, double lat, int zoom, bint truncate):
    xy_ = _xy(lng, lat, truncate=truncate)
    x, y = xy_.x, xy_.y
    
    Z2 = pow(2, zoom)

    if x <= 0:
        xtile = 0
    elif x >= 1:
        xtile = int(Z2 - 1)
    else:
        # To address loss of precision in round-tripping between tile
        # and lng/lat, points within EPSILON of the right side of a tile
        # are counted in the next tile over.
        xtile = int(floor((x + 1e-14) * Z2))

    if y <= 0:
        ytile = 0
    elif y >= 1:
        ytile = int(Z2 - 1)
    else:
        ytile = int(floor((y + 1e-14) * Z2))

    return Tile(xtile, ytile, zoom)


cdef Tile bounding_tile(LngLatBbox bbox, bint truncate):

    w, s, e, n = bbox.west, bbox.south, bbox.east, bbox.north

    if truncate:
        w, s = truncate_lnglat(w, s)
        e, n = truncate_lnglat(e, n)

    e = e - 1e-11
    s = s + 1e-11

    tmin = tile(w, n, 32, 0)
    tmax = tile(e, s, 32, 0)

    cell = tmin[:2] + tmax[:2]
    z = _getBboxZoom(Bbox(tmin[:2], tmax[:2]))

    if z == 0:
        return Tile(0, 0, 0)

    x = rshift(cell[0], (32 - z))
    y = rshift(cell[1], (32 - z))

    return Tile(x, y, z)
