# cyrcantile/_base.pyx

import operator

from libc.math cimport M_PI, log, tan, atan, exp, sinh, floor, sin
from libc.stdlib cimport malloc, free

from ._base cimport (
    MinMax, Tile, LngLat, LngLatBbox, Bbox, XY,
    radians, degrees, minmax, truncate_lnglat,
    _xy_c, _xy_bounds_c, _lnglat_c, _ul_c, _valid_c,
    _neighbours_c, _parent_c, _children_c, _merge, _simplify_c,
    _bounds_c, _getBboxZoom_c, _rshift_c, _tile_c,
    _bounding_tile_c, _quadkey_c, _quadkey_to_tile_c,
    _feature_c, _collect_coords, _geojson_bounds_c,
)

# Константы
R2D = 180 / M_PI
RE = 6378137.0
CE = 2 * M_PI * RE
EPSILON = 1e-14
LL_EPSILON = 1e-11


# -------------------------
# Исключения (Python)
# -------------------------

class MercantileError(Exception):
    """Base exception"""


class InvalidLatitudeError(MercantileError):
    """Raised when math errors occur beyond ~85 degrees N or S"""


class InvalidZoomError(MercantileError):
    """Raised when a zoom level is invalid"""


class ParentTileError(MercantileError):
    """Raised when a parent tile cannot be determined"""


class QuadKeyError(MercantileError):
    """Raised when errors occur in computing or parsing quad keys"""


class TileArgParsingError(MercantileError):
    """Raised when errors occur in parsing a function's tile arg(s)"""


class TileError(MercantileError):
    """Raised when a tile can't be determined"""


# --------------------------------
# Парсер аргумента *tile (Python)
# --------------------------------

def _parse_tile_arg(*args):
    """
    Порт mercantile._parse_tile_arg: принимает либо (Tile,),
    либо (x, y, z) / x, y, z, возвращает (x, y, z).
    """
    if len(args) == 1:
        args = args[0]
    if len(args) == 3:
        x, y, z = args
        return (int(x), int(y), int(z))
    else:
        raise TileArgParsingError(
            "the tile argument may have 1 or 3 values. "
            "Note that zoom is a keyword-only argument"
        )


# ----------------------------------------------------
# Реализации cdef-ядра (в формате, согласованном с .pxd)
# ----------------------------------------------------

cdef double radians(double degrees):
    return degrees * (M_PI / 180.0)


cdef double degrees(double radians):
    return radians * (180.0 / M_PI)


cdef MinMax minmax(int zoom):
    return MinMax(0, 1 << (zoom - 1))


cdef LngLat truncate_lnglat(double lng, double lat):
    if lng > 180.0:
        lng = 180
    elif lng < -180.0:
        lng = -180.0
    if lat > 90.0:
        lat = 90.0
    elif lat < -90.0:
        lat = -90.0
    return LngLat(lng, lat)


cdef XY _xy_c(double lng, double lat, bint truncate):
    cdef LngLat lngLat
    if truncate:
        lngLat = truncate_lnglat(lng, lat)
        lng, lat = lngLat.lng, lngLat.lat

    cdef double x = RE * radians(lng)
    cdef double y
    if lat <= -90:
        y = float("-inf")
    elif lat >= 90:
        y = float("inf")
    else:
        y = RE * log(tan((M_PI * 0.25) + (0.5 * radians(lat))))

    return XY(x, y)


cdef Bbox _xy_bounds_c(Tile tile):
    cdef double tile_size = CE / (1 << tile.z)

    cdef double left = tile.x * tile_size - CE / 2
    cdef double right = left + tile_size

    cdef double top = CE / 2 - tile.y * tile_size
    cdef double bottom = top - tile_size

    return Bbox(left, bottom, right, top)


cdef LngLat _lnglat_c(double x, double y, bint truncate):
    cdef double lng = x * R2D / RE
    cdef double lat = ((M_PI * 0.5) - 2.0 * atan(exp(-y / RE))) * R2D
    cdef LngLat lngLat
    if truncate:
        lngLat = truncate_lnglat(lng, lat)
        lng, lat = lngLat.lng, lngLat.lat
    return LngLat(lng, lat)


cdef LngLat _ul_c(Tile tile):
    cdef int xtile = tile.x
    cdef int ytile = tile.y
    cdef int zoom = tile.z
    cdef double Z2 = 1 << zoom
    cdef double lon_deg = xtile / Z2 * 360.0 - 180.0
    cdef double lat_rad = atan(sinh(M_PI * (1 - 2 * ytile / Z2)))
    cdef double lat_deg = degrees(lat_rad)
    return LngLat(lon_deg, lat_deg)


cdef bint _valid_c(Tile tile):
    cdef bint validx = 0 <= tile.x <= (1 << tile.z) - 1
    cdef bint validy = 0 <= tile.y <= (1 << tile.z) - 1
    cdef bint validz = 0 <= tile.z
    return validx and validy and validz


cdef Tile* _neighbours_c(Tile tile, int* count):
    cdef int i, j
    cdef int lo, hi
    cdef Tile* arr
    cdef int n = 0

    lo, hi = minmax(tile.z).min, minmax(tile.z).max

    arr = <Tile*> malloc(8 * sizeof(Tile))
    if arr == NULL:
        count[0] = 0
        return NULL

    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            if i == 0 and j == 0:
                continue
            if tile.x + i < 0 or tile.y + j < 0:
                continue
            if tile.x + i > hi or tile.y + j > hi:
                continue

            arr[n].x = tile.x + i
            arr[n].y = tile.y + j
            arr[n].z = tile.z
            n += 1

    count[0] = n
    return arr


cdef Tile _parent_c(Tile tile, int zoom=-1):
    cdef int x = tile.x
    cdef int y = tile.y
    cdef int z = tile.z
    cdef int target_zoom
    cdef Tile return_tile = tile
    cdef int xtile, ytile, ztile

    if z == 0:
        raise ParentTileError("the parent of a zoom-0 tile is undefined")

    if zoom != -1:
        if z <= zoom or zoom != <int>zoom:
            raise InvalidZoomError(
                "zoom must be an integer and less than that of the input tile"
            )
        target_zoom = zoom
    else:
        target_zoom = z - 1

    if x != <int>x or y != <int>y or z != <int>z:
        raise ParentTileError("the parent of a non-integer tile is undefined")

    while return_tile.z > target_zoom:
        xtile = return_tile.x
        ytile = return_tile.y
        ztile = return_tile.z

        if xtile % 2 == 0 and ytile % 2 == 0:
            return_tile = Tile(xtile // 2, ytile // 2, ztile - 1)
        elif xtile % 2 == 0:
            return_tile = Tile(xtile // 2, (ytile - 1) // 2, ztile - 1)
        elif xtile % 2 != 0 and ytile % 2 == 0:
            return_tile = Tile((xtile - 1) // 2, ytile // 2, ztile - 1)
        else:
            return_tile = Tile((xtile - 1) // 2, (ytile - 1) // 2, ztile - 1)

    return return_tile


cdef list _children_c(Tile tile, int target_zoom):
    cdef list tiles = [(tile.x, tile.y, tile.z)]
    cdef int xtile, ytile, ztile

    while tiles[0][2] < target_zoom:
        xtile, ytile, ztile = tiles.pop(0)
        tiles += [
            (xtile * 2,     ytile * 2,     ztile + 1),
            (xtile * 2 + 1, ytile * 2,     ztile + 1),
            (xtile * 2 + 1, ytile * 2 + 1, ztile + 1),
            (xtile * 2,     ytile * 2 + 1, ztile + 1),
        ]

    return [Tile(x, y, z) for x, y, z in tiles]


cdef tuple _merge(object merge_set):
    cdef dict upwards_merge = {}
    cdef object tile, tile_parent
    cdef object supertile, children
    cdef list current_tileset = []
    cdef bint changed = False

    for tile in merge_set:
        tile_parent = parent(tile)
        children = upwards_merge.get(tile_parent)
        if children is None:
            children = set()
            upwards_merge[tile_parent] = children
        (<set>children).add(tile)

    for supertile, children in upwards_merge.items():
        if len(children) == 4:
            current_tileset.append(supertile)
            changed = True
        else:
            current_tileset.extend(children)

    return current_tileset, changed


cdef object _simplify_c(object tiles):
    cdef object root_set = set()
    cdef object tile, supertile
    cdef int x, y, z, i
    cdef bint is_new_tile
    cdef bint is_merging = True
    cdef list sorted_tiles
    cdef tuple res

    sorted_tiles = sorted(tiles, key=operator.itemgetter(2))
    for tile in sorted_tiles:
        x, y, z = tile
        is_new_tile = True
        for i in range(z):
            supertile = parent(tile, zoom=i)
            if supertile in root_set:
                is_new_tile = False
                break
        if is_new_tile:
            (<set>root_set).add(tile)

    while is_merging:
        res = _merge(root_set)
        root_set = res[0]
        is_merging = <bint>res[1]

    return root_set


cdef LngLatBbox _bounds_c(Tile tile):
    cdef int xtile = tile.x
    cdef int ytile = tile.y
    cdef int zoom = tile.z
    cdef double Z2 = 1 << zoom

    cdef double ul_lon_deg = xtile / Z2 * 360.0 - 180.0
    cdef double ul_lat_rad = atan(sinh(M_PI * (1 - 2 * ytile / Z2)))
    cdef double ul_lat_deg = degrees(ul_lat_rad)

    cdef double lr_lon_deg = (xtile + 1) / Z2 * 360.0 - 180.0
    cdef double lr_lat_rad = atan(sinh(M_PI * (1 - 2 * (ytile + 1) / Z2)))
    cdef double lr_lat_deg = degrees(lr_lat_rad)

    return LngLatBbox(ul_lon_deg, lr_lat_deg, lr_lon_deg, ul_lat_deg)


cdef int _getBboxZoom_c(unsigned int cell0,
                         unsigned int cell1,
                         unsigned int cell2,
                         unsigned int cell3):
    cdef int MAX_ZOOM = 28
    cdef int z
    cdef unsigned int mask

    for z in range(0, MAX_ZOOM):
        mask = 1 << (32 - (z + 1))
        if (cell0 & mask) != (cell2 & mask) or (cell1 & mask) != (cell3 & mask):
            return z
    return MAX_ZOOM


cdef unsigned int _rshift_c(unsigned int val, int n):
    return val >> n


cdef Tile _tile_c(double lng, double lat, int zoom, bint truncate):
    cdef double x, y
    if truncate:
        lng, lat = truncate_lnglat(lng, lat)

    x = lng / 360.0 + 0.5
    cdef double sinlat = sin(radians(lat))

    try:
        y = 0.5 - 0.25 * log((1.0 + sinlat) / (1.0 - sinlat)) / M_PI
    except (ValueError, ZeroDivisionError):
        raise InvalidLatitudeError(f"Y can not be computed: lat={lat!r}")

    cdef double Z2 = 1 << zoom
    cdef int xtile
    cdef int ytile

    if x <= 0:
        xtile = 0
    elif x >= 1:
        xtile = <int>(Z2 - 1)
    else:
        xtile = <int>floor((x + EPSILON) * Z2)

    if y <= 0:
        ytile = 0
    elif y >= 1:
        ytile = <int>(Z2 - 1)
    else:
        ytile = <int>floor((y + EPSILON) * Z2)

    return Tile(xtile, ytile, zoom)


cdef Tile _bounding_tile_c(LngLatBbox bbox, bint truncate):
    cdef double w = bbox.west
    cdef double s = bbox.south
    cdef double e = bbox.east
    cdef double n = bbox.north
    cdef LngLat tmp

    if truncate:
        tmp = truncate_lnglat(w, s)
        w, s = tmp.lng, tmp.lat
        tmp = truncate_lnglat(e, n)
        e, n = tmp.lng, tmp.lat

    e = e - LL_EPSILON
    s = s + LL_EPSILON

    cdef Tile tmin
    cdef Tile tmax
    try:
        tmin = _tile_c(w, n, 32, 0)
        tmax = _tile_c(e, s, 32, 0)
    except InvalidLatitudeError:
        return Tile(0, 0, 0)

    cdef unsigned int cell0 = <unsigned int>tmin.x
    cdef unsigned int cell1 = <unsigned int>tmin.y
    cdef unsigned int cell2 = <unsigned int>tmax.x
    cdef unsigned int cell3 = <unsigned int>tmax.y

    cdef int z = _getBboxZoom_c(cell0, cell1, cell2, cell3)

    if z == 0:
        return Tile(0, 0, 0)

    cdef int x = <int>_rshift_c(cell0, 32 - z)
    cdef int y = <int>_rshift_c(cell1, 32 - z)

    return Tile(x, y, z)


cdef str _quadkey_c(Tile tile):
    cdef int xtile = tile.x
    cdef int ytile = tile.y
    cdef int zoom = tile.z
    cdef list qk = []
    cdef int z
    cdef int digit
    cdef int mask

    for z in range(zoom, 0, -1):
        digit = 0
        mask = 1 << (z - 1)
        if xtile & mask:
            digit += 1
        if ytile & mask:
            digit += 2
        qk.append(str(digit))
    return "".join(qk)


cdef Tile _quadkey_to_tile_c(str qk):
    if len(qk) == 0:
        return Tile(0, 0, 0)
    cdef int xtile = 0
    cdef int ytile = 0
    cdef int i
    cdef str digit
    cdef int mask
    for i, digit in enumerate(reversed(qk)):
        mask = 1 << i
        if digit == "1":
            xtile = xtile | mask
        elif digit == "2":
            ytile = ytile | mask
        elif digit == "3":
            xtile = xtile | mask
            ytile = ytile | mask
        elif digit != "0":
            raise QuadKeyError(f"Unexpected quadkey digit: {digit!r}")
    return Tile(xtile, ytile, i + 1)


cdef dict _feature_c(Tile t,
                     object fid,
                     dict props,
                     bint mercator,
                     bint has_buffer,
                     double buffer,
                     bint has_precision,
                     int precision):
    cdef LngLatBbox bb
    cdef double west, south, east, north
    cdef XY p
    cdef dict geom
    cdef dict feat
    cdef list coords
    cdef list ring
    cdef list bbox
    cdef str xyz

    bb = _bounds_c(t)
    west, south, east, north = bb.west, bb.south, bb.east, bb.north

    if mercator:
        p = _xy_c(west, south, 0)
        west, south = p.x, p.y
        p = _xy_c(east, north, 0)
        east, north = p.x, p.y

    if has_buffer:
        west  -= buffer
        south -= buffer
        east  += buffer
        north += buffer

    if has_precision and precision >= 0:
        west  = round(west,  precision)
        south = round(south, precision)
        east  = round(east,  precision)
        north = round(north, precision)

    bbox = [
        west if west < east else east,
        south if south < north else north,
        west if west > east else east,
        south if south > north else north,
    ]

    ring = [
        [west, south],
        [west, north],
        [east, north],
        [east, south],
        [west, south],
    ]
    coords = [ring]
    geom = {
        "type": "Polygon",
        "coordinates": coords,
    }

    xyz = str((t.x, t.y, t.z))

    feat = {
        "type": "Feature",
        "bbox": bbox,
        "id": xyz,
        "geometry": geom,
        "properties": {"title": "XYZ tile %s" % xyz},
    }

    if props is not None:
        feat["properties"].update(props)

    if fid is not None:
        feat["id"] = fid

    return feat


cdef void _collect_coords(object obj, list out):
    cdef object coordinates
    cdef object e

    if isinstance(obj, (tuple, list)):
        coordinates = obj
    elif isinstance(obj, dict) and "features" in obj:
        coordinates = [feat["geometry"]["coordinates"] for feat in obj["features"]]
    elif isinstance(obj, dict) and "geometry" in obj:
        coordinates = obj["geometry"]["coordinates"]
    else:
        if isinstance(obj, dict):
            coordinates = obj.get("coordinates", obj)
        else:
            coordinates = obj

    for e in coordinates:
        if isinstance(e, (float, int)):
            out.append(tuple(coordinates))
            break
        else:
            _collect_coords(e, out)


cdef LngLatBbox _geojson_bounds_c(object obj):
    cdef double w = 180.0
    cdef double s = 90.0
    cdef double e = -180.0
    cdef double n = -90.0
    cdef object coord
    cdef double lng, lat
    cdef LngLatBbox bbox

    for coord in _coords(obj):
        lng = <double>coord[0]
        lat = <double>coord[1]

        if lng < w:
            w = lng
        if lat < s:
            s = lat
        if lng > e:
            e = lng
        if lat > n:
            n = lat

    bbox.west = w
    bbox.south = s
    bbox.east = e
    bbox.north = n
    return bbox


# ---------------------------------
# Python-обёртки (публичный API)
# ---------------------------------

def xy(double lng, double lat, bint truncate=False):
    """
    Совместимо с mercantile.xy: возвращает (x, y).
    """
    cdef XY p = _xy_c(lng, lat, truncate)
    return p.x, p.y


def xy_bounds(*tile):
    t = _parse_tile_arg(*tile)
    cdef Bbox b = _xy_bounds_c(Tile(t[0], t[1], t[2]))
    return b.left, b.bottom, b.right, b.top


def lnglat(double x, double y, bint truncate=False):
    return _lnglat_c(x, y, truncate)


def ul(*tile):
    t = _parse_tile_arg(*tile)
    return _ul_c(Tile(t[0], t[1], t[2]))


def valid(Tile tile):
    return _valid_c(tile)


def neighbors(*tile, **kwargs):
    t = _parse_tile_arg(*tile)
    cdef Tile tt = Tile(t[0], t[1], t[2])

    cdef int count = 0
    cdef Tile* arr = _neighbours_c(tt, &count)
    cdef list tiles = []
    cdef int i
    cdef Tile v

    if arr != NULL:
        for i in range(count):
            v = arr[i]
            if _valid_c(v):
                tiles.append(Tile(v.x, v.y, v.z))
        free(arr)

    return tiles


def simplify(tiles):
    return _simplify_c(tiles)


def bounds(*tile):
    t = _parse_tile_arg(*tile)
    return _bounds_c(Tile(t[0], t[1], t[2]))


def tile(double lng, double lat, int zoom, bint truncate=False):
    return _tile_c(lng, lat, zoom, truncate)


def bounding_tile(*bbox, **kwds):
    if len(bbox) == 2:
        bbox += bbox

    w, s, e, n = bbox
    truncate = bool(kwds.get("truncate"))

    if truncate:
        w, s = truncate_lnglat(w, s)
        e, n = truncate_lnglat(e, n)

    e = e - LL_EPSILON
    s = s + LL_EPSILON

    try:
        tmin = tile(w, n, 32, False)
        tmax = tile(e, s, 32, False)
    except InvalidLatitudeError:
        return Tile(0, 0, 0)

    cell = tmin[:2] + tmax[:2]
    z = _getBboxZoom_c(
        <unsigned int>cell[0],
        <unsigned int>cell[1],
        <unsigned int>cell[2],
        <unsigned int>cell[3],
    )

    if z == 0:
        return Tile(0, 0, 0)

    x = _rshift_c(<unsigned int>cell[0], (32 - z))
    y = _rshift_c(<unsigned int>cell[1], (32 - z))

    return Tile(<int>x, <int>y, z)


def quadkey(*tile):
    t = _parse_tile_arg(*tile)
    return _quadkey_c(Tile(t[0], t[1], t[2]))


def quadkey_to_tile(str qk):
    return _quadkey_to_tile_c(qk)


def parent(*tile, **kwargs):
    tile_t = _parse_tile_arg(*tile)
    x, y, z = tile_t

    if z == 0:
        return None

    zoom = kwargs.get("zoom", None)

    if zoom is not None and (z <= zoom or zoom != int(zoom)):
        raise InvalidZoomError(
            "zoom must be an integer and less than that of the input tile"
        )

    if x != int(x) or y != int(y) or z != int(z):
        raise ParentTileError("the parent of a non-integer tile is undefined")

    target_zoom = z - 1 if zoom is None else zoom

    return_tile = tile_t
    while return_tile[2] > target_zoom:
        xtile, ytile, ztile = return_tile
        if xtile % 2 == 0 and ytile % 2 == 0:
            return_tile = Tile(xtile // 2, ytile // 2, ztile - 1)
        elif xtile % 2 == 0:
            return_tile = Tile(xtile // 2, (ytile - 1) // 2, ztile - 1)
        elif not xtile % 2 == 0 and ytile % 2 == 0:
            return_tile = Tile((xtile - 1) // 2, ytile // 2, ztile - 1)
        else:
            return_tile = Tile((xtile - 1) // 2, (ytile - 1) // 2, ztile - 1)

    return return_tile


def children(*tile, **kwargs):
    tile_t = _parse_tile_arg(*tile)

    zoom = kwargs.get("zoom", None)

    xtile, ytile, ztile = tile_t

    if zoom is not None and (ztile > zoom or zoom != int(zoom)):
        raise InvalidZoomError(
            "zoom must be an integer and greater than that of the input tile"
        )

    target_zoom = zoom if zoom is not None else ztile + 1

    tiles = [Tile(xtile, ytile, ztile)]

    while tiles[0].z < target_zoom:
        t = tiles.pop(0)
        xtile, ytile, ztile = t.x, t.y, t.z
        tiles += [
            Tile(xtile * 2,     ytile * 2,     ztile + 1),
            Tile(xtile * 2 + 1, ytile * 2,     ztile + 1),
            Tile(xtile * 2 + 1, ytile * 2 + 1, ztile + 1),
            Tile(xtile * 2,     ytile * 2 + 1, ztile + 1),
        ]

    return tiles


def rshift(val, n):
    return (val % 0x100000000) >> n


def feature(
    tile,
    fid=None,
    props=None,
    projected="geographic",
    buffer=None,
    precision=None,
):
    _bounds =  bounds(tile)

    if projected == "mercator":
        xys = xy(_bounds['west'], _bounds['south'], truncate=False)
        xyn = xy(_bounds['east'], _bounds['north'], truncate=False)
        west, south = xys
        east, north = xyn

    if buffer:
        _bounds['west'] -= buffer
        _bounds['south'] -= buffer
        _bounds['east'] += buffer
        _bounds['north'] += buffer

    if precision and precision >= 0:
        west, south, east, north = (
            round(v, precision) for v in (_bounds['west'], _bounds['south'], _bounds['east'], _bounds['north'])
        )

    bbox = [min(_bounds['west'], _bounds['east']), min(_bounds['south'], _bounds['north']), max(_bounds['west'], _bounds['east']), max(_bounds['south'], _bounds['north'])]
    geom = {
        "type": "Polygon",
        "coordinates": [
            [[_bounds['west'], _bounds['south']], [_bounds['west'], _bounds['north']], [_bounds['east'], _bounds['north']], [_bounds['east'], _bounds['south']], [_bounds['west'], _bounds['south']]]
        ],
    }

    xyz = str(tile)
    feat = {
        "type": "Feature",
        "bbox": bbox,
        "id": xyz,
        "geometry": geom,
        "properties": {"title": "XYZ tile %s" % xyz},
    }

    if props:
        feat["properties"].update(props)

    if fid is not None:
        feat["id"] = fid

    return feat


def _coords(obj):
    cdef list acc = []
    cdef object xy_

    _collect_coords(obj, acc)
    for xy_ in acc:
        yield xy_[:2]


def geojson_bounds(obj):
    return _geojson_bounds_c(obj)
