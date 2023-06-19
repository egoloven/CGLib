from ..models.point import Point
from ..models.vector import Vector
import numpy as np
import math

# Utils

def make_Uy(sites):
    xs = np.array([site.point.x for site in sites], dtype=float)
    ys = np.array([site.point.y for site in sites], dtype=float)
    return [sites[index] for index in np.lexsort((xs, ys))]

def make_Ux(sites):
    xs = np.array([site.point.x for site in sites], dtype=float)
    ys = np.array([site.point.y for site in sites], dtype=float)
    return [sites[index] for index in np.lexsort((ys, xs))]

def make_U(sites):
    xs = np.array([site.point.x for site in sites], dtype=float)
    ys = np.array([site.point.y for site in sites], dtype=float)
    return {Uy_index: Ux_index for Ux_index, Uy_index in enumerate(np.lexsort((xs, ys)))}

def print_tree(root, indent=0):
    if root is None:
        return

    print(" " * indent + str(root))
    print_tree(root.left, indent + 4)
    print_tree(root.right, indent + 4)

# classes

class Site:
    def __init__(self, point, name=None):
        self.name = name
        self.point = point
        self.edges = []

    def get_name(self):
        if self.name is None:
            return f"({self.point.x}; {self.point.y})"
        return self.name

    def __str__(self):
        return self.get_name()

    __repr__ = __str__

class Edge:
    def __init__(self, left_site=None, right_site=None, start_point=None, end_point=None, infinity=1):
        self.left_site = left_site
        self.right_site = right_site
        self.start_point = start_point
        self.end_point = end_point
        self.infinity = infinity

    def next_site(self, site):
        if self.left_site.point.x == site.point.x and self.left_site.point.y == site.point.y:
            return self.right_site
        return self.left_site

    def __str__(self):
        return f"({self.left_site, self.right_site})"

    __repr__ = __str__

class VoronoiDiagram:
    def __init__(self, Ux):
        self.Ux = Ux
        self.Uy = make_Uy(self.Ux)
        self.U = make_U(self.Ux)
        self.left = None
        self.right = None
        self.edges = []

    def compute_diagram(self):
        if len(self.Ux) == 2:
            self.edges += compute_two_points(self.Ux)
        elif len(self.Ux) == 3:
            self.edges += compute_three_points(self.Ux)
        else:
            median = len(self.Ux) // 2
            self.left = VoronoiDiagram(self.Ux[:median])
            self.right = VoronoiDiagram(self.Ux[median:])
            self.left.compute_diagram()
            self.right.compute_diagram()
            self.merge()
            self.clean()

    def merge(self):
        upper_tangent = find_upper_tangent(self.left, self.right)
        # lower_tangent = find_lower_tangent(self.left, self.right)

        sigma_in = compute_sigma_in(upper_tangent)
        # sigma_out = compute_sigma_out(lower_tangent)
        edge = sigma_in
        left_site = upper_tangent[0]
        right_site = upper_tangent[1]
        print(self)
        while edge is not None:
            print(left_site, right_site)
            self.edges.append(edge)
            edge, left_site, right_site = compute_next_edge(edge, left_site, right_site)

    def clean(self):
        new_left_edges = []
        new_right_edges = []
        for edge in self.left.edges:
            if edge.end_point.x - edge.start_point.x > 0 and edge.infinity == 0.5:
                edge.left_site.edges.remove(edge)
                edge.right_site.edges.remove(edge)

            else:
                new_left_edges.append(edge)

        for edge in self.right.edges:
            if edge.end_point.x - edge.start_point.x < 0 and edge.infinity == 0.5:
                edge.left_site.edges.remove(edge)
                edge.right_site.edges.remove(edge)

            else:
                new_right_edges.append(edge)

        self.left.edges = new_left_edges
        self.right.edges = new_right_edges
        self.edges += self.left.edges
        self.edges += self.right.edges

    def __str__(self):
        name = 'sites:'
        first = True
        for i, site in enumerate(self.Ux):
            if first:
                name += str(site)
                first = False
            else:
                name += ',' + str(site)

        name += f' edges: {len(self.edges)}'
        return name

# functions

def find_intersection_point(breach, edge):
    breach_infinity = breach.infinity
    edge_infinity = edge.infinity

    if breach.start_point == edge.start_point or breach.start_point == edge.end_point:
        return None

    breach_k, breach_b = k_b([Site(breach.start_point), Site(breach.end_point)])
    edge_k, edge_b = k_b([Site(edge.start_point), Site(edge.end_point)])

    # locate x, y coords of intersection

    if breach_k is None:
        if edge_k is None:
            return None
        x = breach_b
        y = edge_k * x + edge_b
    elif breach_k == 0:
        if edge_k == 0:
            return None
        y = breach_b
        x = (y - edge_b) / edge_k
    elif edge_k is None:
        x = edge_b
        y = breach_k * x + breach_b
    elif edge_k == 0:
        y = edge_b
        x = (y - breach_b) / breach_k
    elif math.isclose(edge_k, breach_k):
        return None
    else:
        x = - (breach_b - edge_b) / (breach_k - edge_k)
        y = breach_k * x + breach_b

    # check if intersection can exist

    if Point(x, y) == breach.start_point:
        return None
    if Point(x, y) == edge.start_point:
        return None

    if breach_infinity == 1:
        if edge_infinity == 1:
            return Point(x, y)

        if edge_infinity == 0.5:
            edge_end_point = edge.end_point
            edge_start_point = edge.start_point

            first = Vector([edge_end_point.x - edge_start_point.x, edge_end_point.y - edge_start_point.y])
            second = Vector([x - edge_start_point.x, y - edge_start_point.y])

            first.normalize()
            second.normalize()

            if compare_two_vectors(first, second):
                return Point(x, y)

            return None

        if edge_infinity == 0:
            edge_end_point = edge.end_point
            edge_start_point = edge.start_point

            if min(edge_end_point.x, edge_start_point.x) <= x <= max(edge_end_point.x, edge_start_point.x) and min(edge_end_point.y, edge_start_point.y) <= y <= max(edge_end_point.y, edge_start_point.y):
                return Point(x, y)

            return None

    if  breach_infinity == 0.5:
        if edge_infinity == 1:
            breach_end_point = breach.end_point
            breach_start_point = breach.start_point

            first = Vector([breach_end_point.x - breach_start_point.x, breach_end_point.y - breach_start_point.y])
            second = Vector([x - breach_start_point.x, y - breach_start_point.y])

            first.normalize()
            second.normalize()

            if compare_two_vectors(first, second):
                return Point(x, y)

            return None

        if edge_infinity == 0.5:
            breach_end_point = breach.end_point
            breach_start_point = breach.start_point
            edge_end_point = edge.end_point
            edge_start_point = edge.start_point

            first = Vector([breach_end_point.x - breach_start_point.x, breach_end_point.y - breach_start_point.y])
            second = Vector([x - breach_start_point.x, y - breach_start_point.y])
            third = Vector([edge_end_point.x - edge_start_point.x, edge_end_point.y - edge_start_point.y])
            fourth = Vector([x - edge_start_point.x, y - edge_start_point.y])

            first.normalize()
            second.normalize()
            third.normalize()
            fourth.normalize()

            if compare_two_vectors(first, second) and compare_two_vectors(third, fourth):
                return Point(x, y)

            return None

        if edge_infinity == 0:
            breach_end_point = breach.end_point
            breach_start_point = breach.start_point

            first = Vector([breach_end_point.x - breach_start_point.x, breach_end_point.y - breach_start_point.y])
            second = Vector([x - breach_start_point.x, y - breach_start_point.y])

            first.normalize()
            second.normalize()

            if compare_two_vectors(first, second):
                edge_end_point = edge.end_point
                edge_start_point = edge.start_point

                if min(edge_end_point.x, edge_start_point.x) <= x <= max(edge_end_point.x, edge_start_point.x) and min(edge_end_point.y, edge_start_point.y) <= y <= max(edge_end_point.y, edge_start_point.y):
                    return Point(x, y)

                return None

    return None

def update_edges(intersection_point, intersection_edge, edge, left_site):
    if intersection_point.y < edge.start_point.y:
        edge.end_point = edge.start_point
        edge.start_point = intersection_point
        edge.infinity -= 0.5
    else:
        dx = intersection_point.x - edge.start_point.x
        dy = intersection_point.y - edge.start_point.y
        edge.start_point = intersection_point
        edge.end_point = Point(intersection_point.x + dx, intersection_point.y + dy)
        edge.infinity -= 0.5

    if intersection_edge.infinity == 0.5:
        dx = intersection_edge.end_point.x - intersection_edge.start_point.x
        dy = intersection_edge.end_point.y - intersection_edge.start_point.y
        if intersection_edge in left_site.edges:
            a, b, c = a_b_c([Site(edge.start_point), Site(edge.end_point)])

            if intersection_edge.start_point.y * a + intersection_edge.start_point.x * b + c > 0:
                intersection_edge.end_point = intersection_point
                intersection_edge.infinity = 0

            else:
                intersection_edge.start_point = intersection_point

        else:
            a, b, c = a_b_c([Site(edge.start_point), Site(edge.end_point)])

            if intersection_edge.start_point.y * a + intersection_edge.start_point.x * b + c < 0:
                intersection_edge.end_point = intersection_point
                intersection_edge.infinity = 0

            else:
                intersection_edge.start_point = intersection_point

    elif intersection_edge.infinity == 1:
        intersection_edge.start_point = intersection_point
        dx = intersection_edge.end_point.x - intersection_edge.start_point.x
        dy = intersection_edge.end_point.y - intersection_edge.start_point.y

        if intersection_edge in left_site.edges:
            if dx > 0:
                intersection_edge.end_point = Point(intersection_edge.start_point.x - dx, intersection_edge.start_point.y - dy)

        else:
            if dx < 0:
                intersection_edge.end_point = Point(intersection_edge.start_point.x + dx, intersection_edge.start_point.y + dy)

        intersection_edge.infinity = 0.5

    else:
        dx = intersection_edge.end_point.x - intersection_edge.start_point.x
        dy = intersection_edge.end_point.y - intersection_edge.start_point.y

        if intersection_edge in left_site.edges:
            if dx < 0:
                intersection_edge.start_point = intersection_point

            else:
                intersection_edge.end_point = intersection_point

        else:
            if dx < 0:
                intersection_edge.end_point = intersection_point

            else:
                intersection_edge.start_point = intersection_point


def compute_next_edge(edge, left_site, right_site):
    edges = left_site.edges + right_site.edges
    intersection_points = []
    for e in edges:
        point = find_intersection_point(edge, e)
        if point is not None:
            intersection_points.append([point, e])

    intersection_points.sort(key=lambda x: x[0].y)
    if intersection_points == []:
        edge.infinity = 0.5
        left_site.edges.append(edge)
        right_site.edges.append(edge)
        return None, left_site, right_site

    intersection_point, intersection_edge = intersection_points[-1]
    update_edges(intersection_point, intersection_edge, edge, left_site)
    left_site.edges.append(edge)
    right_site.edges.append(edge)

    if intersection_edge in left_site.edges:
        next_site = intersection_edge.next_site(left_site)
        left_site = next_site
    else:
        next_site = intersection_edge.next_site(right_site)
        right_site = next_site
    
    end_point = Point((left_site.point.x + right_site.point.x) / 2, (left_site.point.y + right_site.point.y) / 2)
    dx = end_point.x - intersection_point.x
    dy = end_point.y - intersection_point.y
    if dy < 0:
        next_edge = Edge(left_site=left_site,
                         right_site=right_site,
                         start_point=intersection_point,
                         end_point=end_point,
                         infinity=0.5)
    else:
        end_point = Point(intersection_point.x - dx, intersection_point.y - dy)
        next_edge = Edge(left_site=left_site,
                         right_site=right_site,
                         start_point=intersection_point,
                         end_point=end_point,
                         infinity=0.5)
    return next_edge, left_site, right_site

def compare_two_vectors(vector1, vector2):
    for coord1, coord2 in zip(vector1.coords, vector2.coords):
        if not math.isclose(coord1, coord2):
            return False
    return True

def compute_sigma_in(upper_tangent):
    left_site, right_site = upper_tangent
    k, b = bisector_k_b(upper_tangent)
    if k == None:
        start_y = 500
        start_x = (left_site.point.x + right_site.point.x) / 2
    elif k == 0:
        start_x = 0
        start_y = b
    else:
        start_y = 500
        start_x = (start_y - b) / k
    
    end_x = (left_site.point.x + right_site.point.x) / 2
    end_y = (left_site.point.y + right_site.point.y) / 2
    edge = Edge(left_site=left_site,
                right_site=right_site,
                start_point=Point(start_x, start_y),
                end_point=Point(end_x, end_y),
                infinity=1)

    return edge

def find_most_left(sites):
    most_left = 0
    for i, site in enumerate(sites[1:]):
        if site.point.y != sites[most_left].point.y:
            break
        most_left = i + 1

    return most_left

def find_upper_tangent(left, right):
    leftUy = left.Uy
    rightUy = right.Uy

    leftUy.reverse()
    rightUy.reverse()
    left_max_y = leftUy[0].point.y
    right_max_y = rightUy[0].point.y

    if left_max_y > right_max_y:
        pot_left = []
        i = 0
        while i < len(leftUy) and leftUy[i].point.y > right_max_y:
            pot_left.append(leftUy[i])
            i += 1

        left_site = pot_left[0]
        a, b, c = a_b_c([left_site, rightUy[0]])
        for site in pot_left[1:]:
            if (site.point.y * a + site.point.x * b + c) > 0:
                left_site = site
                a, b, c = a_b_c([left_site, rightUy[0]])

        right_site = rightUy[0]
        for site in right.Ux[right.U[len(right.Uy) - 1]:]:
            if (site.point.y * a + site.point.x * b + c) > 0:
                right_site = site
                a, b, c = a_b_c([left_site, right_site])
        return [left_site, right_site]

    if left_max_y < right_max_y:
        pot_left = []
        i = 0
        while i < len(rightUy) and rightUy[i].point.y > left_max_y:
            pot_left.append(rightUy[i])
            i += 1

        right_site = pot_left[0]
        lind = find_most_left(leftUy)
        left_site = leftUy[lind]
        a, b, c = a_b_c([left_site, right_site])
        for site in pot_left[1:]:
            if (site.point.y * a + site.point.x * b + c) > 0:
                left_site = site
                a, b, c = a_b_c([left_site, right_site])

        for site in left.Ux[left.U[len(left.Uy) - lind - 1]:]:
            if (site.point.y * a + site.point.x * b + c) > 0:
                left_site = site
                a, b, c = a_b_c([left_site, right_site])
        return [left_site, right_site]

    return [leftUy[find_most_left(leftUy)], rightUy[0]]

def k_b(sites):
    point1, point2 = sites[0].point, sites[1].point
    if point1.x == point2.x:
        return [None, point1.x]
    if point1.y == point2.y:
        return [0, point1.y]

    k = (point2.y - point1.y) / (point2.x - point1.x)
    b = - k * point1.x + point1.y
    return [k, b]

def a_b_c(sites):
    point1 = sites[0].point
    point2 = sites[1].point
    return [point2.x - point1.x,
            - (point2.y - point1.y),
            - (point1.y * (point2.x - point1.x)) + (point1.x * (point2.y - point1.y))]

def bisector_k_b(sites):
    site1, site2 = sites
    point1, point2 = site1.point, site2.point

    if point1.x == point2.x:
        return [0, (point1.x + point2.y) / 2]
    if point1.y == point2.y:
        return [None, (point1.x + point2.y) / 2]

    k = - (point2.x - point1.x) / (point2.y - point1.y)
    x_mid, y_mid = (point1.x + point2.x) / 2, (point1.y + point2.y) / 2
    b = - k * x_mid + y_mid
    return [k, b]

def find_circlecenter(sites):
    ab = [sites[0], sites[1]]
    bc = [sites[1], sites[2]]

    k1, b1 = bisector_k_b(ab)
    k2, b2 = bisector_k_b(bc)

    if k1 is None:
        if k2 is None:
            return None
        out_x = b1
        out_y = k2 * out_x + b2
        return Point(out_x, out_y)

    if k2 is None:
        out_x = b2
        out_y = k1 * out_x + b1
        return Point(out_x, out_y)

    if k1 == 0:
        if k2 == 0:
            return None
        out_y = b1
        out_x = (out_y - b2) / k2
        return Point(out_x, out_y)

    if k2 == 0:
        out_y = b2
        out_x = (out_y - b1) / k1
        return Point(out_x, out_y)

    out_x = (b2 - b1) / (k1 - k2)
    out_y = k1 * out_x + b1
    return Point(out_x, out_y)

def compute_two_points(sites):
    site1, site2 = sites
    mid_point_x = (site1.point.x + site2.point.x) / 2
    mid_point_y = (site1.point.y + site2.point.y) / 2

    k, b = bisector_k_b(sites)
    if k is None:
        second_point_x = mid_point_x
        second_point_y = mid_point_y + 1
    elif k == 0:
        second_point_x = mid_point_x + 1
        second_point_y = mid_point_y
    else:
        second_point_x = mid_point_x + 1
        second_point_y = k * second_point_x + b

    edge = Edge(left_site=site1,
                right_site=site2,
                start_point=Point(mid_point_x, mid_point_y),
                end_point=Point(second_point_x, second_point_y),
                infinity=1)
    site1.edges.append(edge)
    site2.edges.append(edge)
    return [edge]

def orientation(segment, point):
    # segment = list[Site]
    # point = Point
    a, b, c = a_b_c(segment)
    line = a * point.y + b * point.x + c
    if line > 0:
        return 1
    if line < 0:
        return -1
    return 0

def compute_three_points(sites):
    edges = []
    sites_length = len(sites)

    intersection_point = find_circlecenter(sites)

    if intersection_point is None:
        edges += compute_two_points([sites[0], sites[1]])
        edges += compute_two_points([sites[1], sites[2]])
        return edges

    for i in range(sites_length):
        mid_x = (sites[i].point.x + sites[(i + 1) % sites_length].point.x) / 2
        mid_y = (sites[i].point.y + sites[(i + 1) % sites_length].point.y) / 2

        third_x = sites[(i + 2) % sites_length].point.x
        third_y = sites[(i + 2) % sites_length].point.y

        third = Point(third_x, third_y)

        third_orientation = orientation([sites[i], sites[(i + 1) % sites_length]], third)
        intersection_orientation = orientation([sites[i], sites[(i + 1) % sites_length]], intersection_point)

        if third_orientation == intersection_orientation:
            dx = mid_x - intersection_point.x
            dy = mid_y - intersection_point.y
        else:
            dx = intersection_point.x - mid_x
            dy = intersection_point.y - mid_y

        start_x = intersection_point.x
        start_y = intersection_point.y
        end_x = intersection_point.x + dx
        end_y = intersection_point.y + dy

        edge = Edge(left_site=sites[i],
                    right_site=sites[(i + 1) % sites_length],
                    start_point=Point(start_x, start_y),
                    end_point=Point(end_x, end_y),
                    infinity=0.5)
        edges.append(edge)
        sites[i].edges.append(edge)
        sites[(i + 1) % sites_length].edges.append(edge)

    return edges

if __name__ == '__main__':
    file = open('save', 'r')
    points_string = file.read()
    points_string_array = points_string.split('\n')
    sites = []
    for point in points_string_array:
        p = point.split(',')
        sites.append(Site(point=Point(int(p[1]), int(p[2])), name=p[0]))

    Ux = make_Ux(sites)
    VD = VoronoiDiagram(Ux)
    VD.compute_diagram()
    print_tree(VD)