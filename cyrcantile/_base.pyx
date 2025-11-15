import operator

from libc.math cimport M_PI, log, tan, atan, exp, sinh, floor, sin

from libc.stdlib cimport malloc

R2D = 180 / M_PI
RE = 6378137.0
CE = 2 * M_PI * RE
EPSILON = 1e-14
LL_EPSILON = 1e-11

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


cdef double radians(double degrees):
    return degrees * (M_PI / 180.0)

cdef double degrees(double radians):
    return radians * (180.0 / M_PI)

cdef MinMax minmax(int zoom):
    return MinMax(0, 1 << (zoom - 1)) # here we bit-shift, so the value is int

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
    return XY(x, y)

cdef Bbox xy_bounds(Tile tile):
    tile_size = CE / 2 ** tile.z

    left = tile.x * tile_size - CE / 2
    right = left + tile_size

    top = CE / 2 - tile.y * tile_size
    bottom = top - tile_size

    return Bbox(left, bottom, right, top)

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
    Z2 = 1 << zoom
    lon_deg = xtile / Z2 * 360.0 - 180.0
    lat_rad = atan(sinh(M_PI * (1 - 2 * ytile / Z2)))
    lat_deg = degrees(lat_rad)
    return LngLat(lon_deg, lat_deg)

cdef bool valid(Tile tile):
    validx = 0 <= tile.x <= 2 ** tile.z - 1
    validy = 0 <= tile.y <= 2 ** tile.z - 1
    validz = 0 <= tile.z
    return validx and validy and validz

cdef Tile* neighbours(Tile tile, int* count):
    cdef int i, j
    cdef int lo, hi
    cdef Tile* arr
    cdef int n = 0

    lo, hi = minmax(tile.z).min, minmax(tile.z).max

    # максимум 8 соседей (вокруг тайла в квадрате 3x3, без центра)
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

            # записываем в C-массив
            arr[n].x = tile.x + i
            arr[n].y = tile.y + j
            arr[n].z = tile.z
            n += 1

    count[0] = n
    return arr

cdef Tile parent(Tile tile, int zoom=-1):
    """
    Cython-версия parent: работает с Tile и целым zoom.
    zoom == -1 означает 'использовать z-1' как в оригинальной функции.
    """
    cdef int x = tile.x
    cdef int y = tile.y
    cdef int z = tile.z
    cdef int target_zoom
    cdef Tile return_tile
    cdef int xtile, ytile, ztile

    if z == 0:
        # поведение можно изменить, если нужно именно None
        raise ParentTileError("the parent of a zoom-0 tile is undefined")

    if zoom != -1:
        # zoom задан явно
        if z <= zoom or zoom != <int>zoom:
            raise InvalidZoomError(
                "zoom must be an integer and less than that of the input tile"
            )
        target_zoom = zoom
    else:
        # zoom не задан: берём непосредственного родителя
        target_zoom = z - 1

    if x != <int>x or y != <int>y or z != <int>z:
        raise ParentTileError("the parent of a non-integer tile is undefined")

    return_tile = tile
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

cdef list children(Tile tile, int target_zoom):
    cdef list tiles = [(tile.x, tile.y, tile.z)]
    cdef int xtile, ytile, ztile

    # строим дерево до нужного зума
    while tiles[0][2] < target_zoom:
        xtile, ytile, ztile = tiles.pop(0)
        tiles += [
            (xtile * 2, ytile * 2, ztile + 1),
            (xtile * 2 + 1, ytile * 2, ztile + 1),
            (xtile * 2 + 1, ytile * 2 + 1, ztile + 1),
            (xtile * 2, ytile * 2 + 1, ztile + 1),
        ]

    # превращаем кортежи в Tile
    return [Tile(x, y, z) for x, y, z in tiles]

cdef tuple _merge(object merge_set):
    """
    Cython-версия merge: на вход любой итерируемый набор тайлов (set/list),
    на выход (current_tileset, changed).
    """
    cdef dict upwards_merge = {}          # parent_tile -> set(children)
    cdef object tile, tile_parent
    cdef object supertile, children
    cdef list current_tileset = []
    cdef bint changed = False

    # Собираем словарь parent -> set(children)
    for tile in merge_set:
        tile_parent = parent(tile)        # Python-обёртка, можно заменить на _parent(...)
        children = upwards_merge.get(tile_parent)
        if children is None:
            children = set()
            upwards_merge[tile_parent] = children
        (<set>children).add(tile)

    # Если у родителя ровно 4 ребёнка — схлопываем в родителя,
    # иначе оставляем детей как есть.
    for supertile, children in upwards_merge.items():
        if len(children) == 4:
            current_tileset.append(supertile)
            changed = True
        else:
            current_tileset.extend(children)

    return current_tileset, changed


cdef object _simplify(object tiles):
    """
    Cython-ядро simplify: принимает любую последовательность тайлов (list/set),
    возвращает упрощённый tileset (список).
    """
    cdef object root_set = set()   # сначала set, потом станет list – как в оригинале
    cdef object tile, supertile
    cdef int x, y, z, i
    cdef bint is_new_tile
    cdef bint is_merging = True
    cdef list sorted_tiles
    cdef tuple res

    # сначала отсеиваем тайлы, полностью покрытые родителями
    sorted_tiles = sorted(tiles, key=operator.itemgetter(2))
    for tile in sorted_tiles:
        x, y, z = tile
        is_new_tile = True

        # перебираем всех родителей по зумам 0..z-1
        for i in range(z):
            supertile = parent(tile, zoom=i)
            if supertile in root_set:
                is_new_tile = False
                break

        if is_new_tile:
            (<set>root_set).add(tile)

    # многократно применяем merge, пока есть, что схлопывать
    while is_merging:
        res = _merge(root_set)
        root_set = res[0]
        is_merging = <bint>res[1]

    return root_set

cdef LngLatBbox bounds(Tile tile):
    xtile, ytile, zoom = tile.x, tile.y, tile.z
    Z2 = 1 << zoom
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

    return XY(x, y)

cdef Tile tile(double lng, double lat, int zoom, bint truncate):
    xy_ = _xy(lng, lat, truncate=truncate)
    x, y = xy_.x, xy_.y

    Z2 = 1 << zoom

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
    cdef double w = bbox.west
    cdef double s = bbox.south
    cdef double e = bbox.east
    cdef double n = bbox.north
    cdef LngLat tmp          # объявляем здесь, на уровне функции

    if truncate:
        tmp = truncate_lnglat(w, s)
        w, s = tmp.lng, tmp.lat
        tmp = truncate_lnglat(e, n)
        e, n = tmp.lng, tmp.lat

    e = e - 1e-11
    s = s + 1e-11

    cdef Tile tmin = tile(w, n, 32, 0)
    cdef Tile tmax = tile(e, s, 32, 0)

    # cell[0] = x_min, cell[1] = y_min, cell[2] = x_max, cell[3] = y_max
    cdef unsigned int cell0 = <unsigned int> tmin.x
    cdef unsigned int cell1 = <unsigned int> tmin.y
    cdef unsigned int cell2 = <unsigned int> tmax.x
    cdef unsigned int cell3 = <unsigned int> tmax.y

    cdef Bbox b
    b.left   = cell0
    b.bottom = cell1
    b.right  = cell2
    b.top    = cell3

    cdef int z = _getBboxZoom(b)

    if z == 0:
        return Tile(0, 0, 0)

    cdef int x = <int> rshift(cell0, 32 - z)
    cdef int y = <int> rshift(cell1, 32 - z)

    return Tile(x, y, z)

cdef str quadkey(Tile tile):
    xtile, ytile, zoom = tile.x, tile.y, tile.z
    qk = []
    for z in range(zoom, 0, -1):
        digit = 0
        mask = 1 << (z - 1)
        if xtile & mask:
            digit += 1
        if ytile & mask:
            digit += 2
        qk.append(str(digit))
    return "".join(qk)

cdef Tile quadkey_to_tile(str qk):
    """Get the tile corresponding to a quadkey

    Parameters
    ----------
    qk : str
        A quadkey string.

    Returns
    -------
    Tile

    """
    if len(qk) == 0:
        return Tile(0, 0, 0)
    xtile, ytile = 0, 0
    for i, digit in enumerate(reversed(qk)):
        mask = 1 << i
        if digit == "1":
            xtile = xtile | mask
        elif digit == "2":
            ytile = ytile | mask
        elif digit == "3":
            xtile = xtile | mask
            ytile = ytile | mask
    return Tile(xtile, ytile, i + 1)


cdef dict _feature(Tile t,
                   object fid,
                   dict props,
                   bint mercator,
                   bint has_buffer,
                   double buffer,
                   bint has_precision,
                   int precision):
    """
    Внутреннее Cython-ядро: работает с Tile и уже разобранными флагами.
    """
    cdef LngLatBbox bb
    cdef double west, south, east, north
    cdef double w2, s2
    cdef XY p
    cdef dict geom
    cdef dict feat
    cdef list coords
    cdef list ring
    cdef list bbox
    cdef str xyz

    # исходный bbox в географических координатах
    bb = bounds(t)
    west, south, east, north = bb.west, bb.south, bb.east, bb.north

    if mercator:
        # перекладываем в Web Mercator через xy(...)
        p = xy(west, south, False)
        west, south = p.x, p.y
        p = xy(east, north, False)
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

    # bbox = [min(west, east), min(south, north), max(west, east), max(south, north)]
    bbox = [
        west if west < east else east,
        south if south < north else north,
        west if west > east else east,
        south if south > north else north,
    ]

    # geometry
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
    """
    Внутреннее Cython-ядро: рекурсивно собирает все (lng, lat) в список out.
    """
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
            # дошли до списка чисел: считаем его координатной парой/кортежом
            out.append(tuple(coordinates))
            break
        else:
            _collect_coords(e, out)


def _coords(obj):
    """
    Генератор с тем же API, что и исходный _coords, но использующий
    cdef-ядро для обхода структуры.
    """
    cdef list acc = []
    cdef object xy

    _collect_coords(obj, acc)

    for xy in acc:
        # xy уже кортеж, но берём только (lng, lat), как в исходнике
        yield xy[:2]


cdef LngLatBbox geojson_bounds_c(object obj):
    cdef double w = 180.0
    cdef double s = 90.0
    cdef double e = -180.0
    cdef double n = -90.0
    cdef object coord
    cdef double lng, lat
    cdef LngLatBbox bbox

    # эквивалент reduce(func, _coords(obj), (180, 90, -180, -90))
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