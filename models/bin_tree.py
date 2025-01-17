from .bin_tree_node import Node, NodeWithParent
from .point import Point


class BinTree:
    def __init__(self, root: Node):
        self.root = root

    def __eq__(self, other):
        return self.root == other.root
    
    @property
    def nodes(self):
        """
            Returns the tree represented as left-to-right list of tuples
            with nodes' data, and the data of their left and right children.
        """
        return self._nodes(self.root)
    
    def _nodes(self, node, result=None):
        if result is None:
            result = []
        
        left_data = node.left.data if node.left else None
        right_data = node.right.data if node.right else None

        if node:
            result.append((node.data, left_data, right_data))
        if node.left:
            self._nodes(node.left, result)
        if node.right:
            self._nodes(node.right, result)
        
        return result


class KdTree(BinTree):
    def __init__(self, root: Node, x_range, y_range):
        super().__init__(root)
        self.x_range = x_range
        self.y_range = y_range
        self.partition = []
        self.search_list = []

    def make_tree(self, points, node: Node, vertical=True):
        med = len(points) // 2
        part = (points[med], vertical)
        
        if all(p[0] != part[0] for p in self.partition):
            self.partition.append(part)

        if med == 0:
            return

        if vertical:
            sort_key = lambda p: p.y
        else:
            sort_key = lambda p: p.x
        
        list_l = sorted(points[:med], key=sort_key)
        list_r = sorted(points[-med:], key=sort_key)
        left, right = list_l[med // 2], list_r[med // 2]

        node.left = Node(left)
        if node.data != right:
            node.right = Node(right)

        self.make_tree(list_l, node.left, not vertical)
        self.make_tree(list_r, node.right, not vertical)

    def region_search(self, node: Node, vertical=True):
        if vertical:
            left, right, coord = self.x_range[0], self.x_range[1], node.data.x
        else:
            left, right, coord = self.y_range[0], self.y_range[1], node.data.y

        dots = []
        to_add = self.dot_in_region(node.data)

        if to_add:
            dots.append(node.data)

        intersection = left <= coord <= right
        self.search_list.append((node.data, to_add, intersection))

        if node.left and left < coord:
            dots.extend(self.region_search(node.left, not vertical))
        if node.right and coord < right:
            dots.extend(self.region_search(node.right, not vertical))

        return dots

    def dot_in_region(self, dot):
        return (
            self.x_range[0] <= dot.x
            and dot.x <= self.x_range[1]
            and self.y_range[0] <= dot.y
            and dot.y <= self.y_range[1]
        )


class ChainsBinTree(BinTree):
    def make_tree(self, list, node):
        mid = len(list) // 2
        if mid == 0:
            return

        list_l = list[:mid]
        list_r = list[-mid:]
        left, right = list_l[mid // 2], list_r[mid // 2]

        node.left = NodeWithParent(left, node)
        if node.data != right:
            node.right = NodeWithParent(right, node)

        self.make_tree(list_l, node.left)
        self.make_tree(list_r, node.right)

    @staticmethod
    def _point_in_edge(edge, point):
        return edge.v1[1] <= point.y and edge.v2[1] >= point.y

    @staticmethod
    def _location_against_edge(point, edge):
        return Point.direction(edge.v1.point, edge.v2.point, point)

    def search_point(self, point):
        '''Returns a pair of chains the point is between'''
        current_node = self.root
        left_parent, right_parent = None, None

        while current_node:
            edge = list(
                filter(lambda e: ChainsBinTree._point_in_edge(e, point), current_node.data)
            )[0]
            location = ChainsBinTree._location_against_edge(point, edge)

            if location > 0:
                if current_node.right is not None:
                    current_node = current_node.right
                    left_parent = current_node.parent
                else:
                    return (current_node.data, right_parent.data)

            elif location < 0:
                if current_node.left is not None:
                    current_node = current_node.left
                    right_parent = current_node.parent
                else:
                    return (left_parent.data, current_node.data)
            else:
                return (current_node.data, None)
