from mathutils import Vector, Euler, geometry


def extractPerimeterEdges(obj):
    el = []
    perifEdgeList = []
    for p in obj.data.polygons:
        for e in p.edge_keys:
            if e[0] > e[1]:
                el.append((e[1], e[0]))
            else:
                el.append(e)
    while el:
        e = el.pop()
        if e in el:
            el.remove(e)
        else:
            perifEdgeList.append(e)
    return perifEdgeList


def cubic(p, t):
    return p[0] * (1.0 - t) ** 3.0 + 3.0 * p[1] * t * (1.0 - t) ** 2.0 \
           + 3.0 * p[2] * (t ** 2.0) * (1.0 - t) + p[3] * t ** 3.0


def getbezpoints(spl, mt, seg=0):
    points = spl.bezier_points
    p0 = mt * points[seg].co
    p1 = mt * points[seg].handle_right
    p2 = mt * points[seg + 1].handle_left
    p3 = mt * points[seg + 1].co
    return p0, p1, p2, p3


def getnurbspoints(spl, mw):
    pts = []
    ws = []
    for p in spl.points:
        v = Vector(p.co[0:3]) * mw
        pts.append(v)
        ws.append(p.weight)
    return pts, ws


def knots(n, order, type=0):  # 0 uniform 1 endpoints 2 bezier
        kv = []

        t = n + order
        if type == 0:
            for i in range(0, t):
                kv.append(1.0 * i)

        elif type == 1:
            k = 0.0
            for i in range(1, t + 1):
                kv.append(k)
                if i >= order and i <= n:
                    k += 1.0
        elif type == 2:
            if order == 4:
                k = 0.34
                for a in range(0, t):
                    if a >= order and a <= n: k += 0.5
                    kv.append(floor(k))
                    k += 1.0 / 3.0

            elif order == 3:
                k = 0.6
                for a in range(0, t):
                    if n >= a >= order: k += 0.5
                    kv.append(floor(k))

        # normalize the knot vector
        for i in range(0, len(kv)):
            kv[i] = kv[i] / kv[-1]

        return kv


def B(i, k, t, knots):
    ret = 0
    if k > 0:
        n1 = (t - knots[i]) * B(i, k - 1, t, knots)
        d1 = knots[i + k] - knots[i]
        n2 = (knots[i + k + 1] - t) * B(i + 1, k - 1, t, knots)
        d2 = knots[i + k + 1] - knots[i + 1]
        if d1 > 0.0001 or d1 < -0.0001:
            a = n1 / d1
        else:
            a = 0
        if d2 > 0.0001 or d2 < -0.0001:
            b = n2 / d2
        else:
            b = 0
        ret = a + b
        # print "B i = %d, k = %d, ret = %g, a = %g, b = %g\n"%(i,k,ret,a,b)
    else:
        if knots[i] <= t <= knots[i+1]:
            ret = 1
        else:
            ret = 0
    return ret


def C(t, order, points, weights, knots):
    # c = Point([0,0,0])
    c = Vector()
    rational = 0
    i = 0
    while i < len(points):
        b = B(i, order, t, knots)
        p = points[i] * (b * weights[i])
        c = c + p
        rational = rational + b * weights[i]
        i = i + 1

    return c * (1.0 / rational)

    # Return the coordinate of a point in one object
    # at t percentage of the entire object extension


def calct(obj, t):
    if obj.type == 'CURVE':
        spl = None
        mw = obj.matrix_world
        if obj.data.splines.active is None:
            if len(obj.data.splines) > 0:
                spl = obj.data.splines[0]
        else:
            spl = obj.data.splines.active
        if spl is None:
            return False
        if spl.type == "BEZIER":
            points = spl.bezier_points
            nsegs = len(points) - 1
            d = 1.0 / nsegs
            seg = int(t / d)
            t1 = t / d - int(t / d)
            if t == 1:
                seg -= 1
                t1 = 1.0
            p = getbezpoints(spl, mw, seg)
            coord = cubic(p, t1)
            return coord
        elif spl.type == "NURBS":
            data = getnurbspoints(spl, mw)
            pts = data[0]
            ws = data[1]
            order = spl.order_u
            n = len(pts)
            ctype = spl.use_endpoint_u
            kv = knots(n, order, ctype)
            coord = C(t, order - 1, pts, ws, kv)
            return coord
    elif obj.type == 'MESH':
        if t == 1:
            t = 0.999
        mw = obj.matrix_world
        vNum = len(obj.data.vertices)
        vidx = int(vNum * t)
        t1 = vNum * t - int(vNum * t)
        print(
            "CALCT: The object", obj.name, " has num vertices=", vNum, ' index calculated is=', vidx, "with a t=",
            t,
            "and t1=", t1)
        if vidx == 0:
            coord = mw * (
                    (obj.data.vertices[vidx + 1].co - obj.data.vertices[vidx].co) * t1 + obj.data.vertices[vidx].co)
        else:
            coord = mw * (
                    (obj.data.vertices[vidx].co - obj.data.vertices[vidx - 1].co) * t1 + obj.data.vertices[vidx - 1].co)
        print("CALCT: returning the coordinate", coord)
        return coord
    else:
        assert false


def intc(objs, i, stepPercent, spanPercent, tipo=3, tension=0.0, bias=0.0):
    ncurves = len(objs)
    # if 2 curves go to linear interpolation regardless the one you choose
    if ncurves < 3:
        return intl(objs, i, stepPercent, spanPercent)
    else:
        # calculates the points to be interpolated on each curve
        if i == 0:
            p0 = calct(objs[i], stepPercent)
            p1 = p0
            p2 = calct(objs[i + 1], stepPercent)
            p3 = calct(objs[i + 2], stepPercent)
        else:
            if ncurves - 2 == i:
                p0 = calct(objs[i - 1], stepPercent)
                p1 = calct(objs[i], stepPercent)
                p2 = calct(objs[i + 1], stepPercent)
                p3 = p2
            else:
                p0 = calct(objs[i - 1], stepPercent)
                p1 = calct(objs[i], stepPercent)
                p2 = calct(objs[i + 1], stepPercent)
                p3 = calct(objs[i + 2], stepPercent)

    # calculates the interpolation between those points
    # i used methods from this page: http://paulbourke.net/miscellaneous/interpolation/
    if tipo == 0:
        # linear
        return intl(objs, i, stepPercent, spanPercent)
    elif tipo == 1:
        # natural cubic
        t2 = spanPercent * spanPercent
        a0 = p3 - p2 - p0 + p1
        a1 = p0 - p1 - a0
        a2 = p2 - p0
        a3 = p1
        return a0 * spanPercent * t2 + a1 * t2 + a2 * spanPercent + a3
    elif tipo == 2:
        # catmull it seems to be working. ill leave it for now.
        t2 = spanPercent * spanPercent
        a0 = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3
        a1 = p0 - 2.5 * p1 + 2 * p2 - 0.5 * p3
        a2 = -0.5 * p0 + 0.5 * p2
        a3 = p1
        return a0 * spanPercent * spanPercent + a1 * t2 + a2 * spanPercent + a3
    elif tipo == 3:
        # hermite
        tr2 = spanPercent * spanPercent
        tr3 = tr2 * spanPercent
        m0 = (p1 - p0) * (1 + bias) * (1 - tension) / 2
        m0 += (p2 - p1) * (1 - bias) * (1 - tension) / 2
        m1 = (p2 - p1) * (1 + bias) * (1 - tension) / 2
        m1 += (p3 - p2) * (1 - bias) * (1 - tension) / 2
        a0 = 2 * tr3 - 3 * tr2 + 1
        a1 = tr3 - 2 * tr2 + spanPercent
        a2 = tr3 - tr2
        a3 = -2 * tr3 + 3 * tr2

        return a0 * p1 + a1 * m0 + a2 * m1 + a3 * p2


def intl(objs, i, t, tr):
    p1 = calct(objs[i], t)
    p2 = calct(objs[i + 1], t)
    r = p1 + (p2 - p1) * tr
    return r


def loft(objs, steps, spans, interpolation=0, tension=0.0, bias=0.5):
    verts = []

    # for each object
    for i in range(0, len(objs)):
        # For each step
        for j in range(0, steps + 1):
            # t = percentage of steps according to j
            stepPercent = 1.0 * j / steps
            # verts filled in with the coordinate at the
            # point t of the curve
            verts.append(calct(objs[i], stepPercent))
        temp2 = []
        if i < len(objs) - 1:
            for l in range(1, spans):
                spanPercent = 1.0 * l / spans
                for k in range(0, steps + 1):
                    stepPercent = 1.0 * k / steps
                    if interpolation:
                        pos = intc(objs, i, stepPercent, spanPercent, tipo=interpolation, tension=tension,
                                        bias=bias)
                    else:
                        pos = intl(objs, i, stepPercent, spanPercent)

                    temp2.append(pos)
            verts.extend(temp2)
    return verts


def profile(x, m, p):
    # =====================================================================================
    # The profile calculator
    #
    # inputs
    #   x  : point in the range 0 to 1 to calculate
    #   m : camber percentage
    #   p : camber position percentage
    # =====================================================================================
    if (x < 0):
        return 0
    if (x > 1):
        return 0
    if (x < p):
        y = (m / (p * p)) * (2 * p * x - x * x);
    else:
        y = (m / ((1 - p) * (1 - p))) * ((1 - 2 * p) + 2 * p * x - x * x);
    # debug(("Profile return Z="+str(y))
    return y


def curveProfile(x, vxs):
    for i in vxs:
        if i[0] > x:
            return i[1]
    return 0
