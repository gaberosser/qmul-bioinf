from apps.py2cytoscape import cyrest
from apps.py2cytoscape.data.cyrest_client import CyRestClient
from apps.py2cytoscape.data import style

import networkx as nx
from utils import log


class CytoscapeSession(object):
    def __init__(self, clear_existing=True):
        self.logger = log.get_console_logger(self.__class__.__name__)

        # functional API - the python bindings are incomplete here?
        self.cy = CyRestClient()

        if clear_existing:
            # reset the session (in case something is already loaded)
            self.cy.session.delete()

        # command API - the python bindings are much better
        self.cy_cmd = cyrest.cyclient()

        # collections added to the session
        self.name_to_id = {}
        self.collections = {}
        self.auto_net_name = 1

    def add_networkx_graph(self, graph, name=None):
        old_net = None
        if name is None:
            name = self.auto_net_name
            self.auto_net_name += 1
        elif name in self.collections:
            # overwrite if clashing
            self.logger.warn("Net already exists with name %s. Overwriting.", name)
            # get the network before deleting
            old_net = self.collections[name]
            self.cy.network.delete(old_net.obj)

        try:

            this_style = self.cy.style.create(name)
            this_obj = self.cy.network.create_from_networkx(graph, collection=name)

            self.name_to_id[name] = this_obj.get_id()

            self.collections[name] = CytoNet(this_obj, this_style, cy_session=self.cy)
            return self.collections[name]

        except Exception:
            if old_net is not None:
                self.logger.warn("Restoring previous graph")
                self.cy.network.create(name=name, collection=name, data=old_net.obj.to_json())
            raise

    def apply_layout(self, name=None, layout_type='force_directed', **kwargs):
        """
        Apply the specified layout to either a single network OR all networks
        :param name: If specified, this is a string giving a single network name, otherwise we iterate over all of them.
        :param layout_type:
        :param kwargs: Passed to the layout call.
        :return:
        """
        if name is None:
            targets = self.collections.keys()
        else:
            targets = [name]
        layout_func = getattr(self.cy_cmd.layout, layout_type)
        for t in targets:
            layout_func(network=t, **kwargs)


class CytoNet(object):
    allowed_types = {'String', 'Double'}

    def __init__(self, obj, style, cy_session=None):
        self.style = style
        self.obj = obj
        if cy_session is None:
            cy_session = CyRestClient()
        self.cy_session = cy_session

    def apply_style(self):
        if self.cy_session is None:
            raise AttributeError("Unable to apply styles without knowledge of the Cytoscape session")
        self.cy_session.style.apply(self.style, network=self.obj)

    def _create_passthrough_mapping(self, column, prop, col_type='String'):
        if col_type not in self.allowed_types:
            raise TypeError("col_type must be one of %s.", ', '.join(self.allowed_types))
        self.style.create_passthrough_mapping(column=column, vp=prop, col_type=col_type)
        self.apply_style()

    def _create_linear_mapping(
        self,
        column,
        prop,
        x=None,
        y=None,
        col_type='Double',
        col_belongs_to='node'
    ):
        """
        Create a continuous linear mapping.
        :param column:
        :param prop:
        :param x: Length 2, containing (xmin, xmax) for the interpolation. If not supplied, xmin=0, xmax=max(data)
        :param y: Length 2, containing (ymin, ymax) for the interpolation. If not supplied, ymin=10, yxmax=50
        :param col_type:
        :return:
        """
        if col_belongs_to not in {'edge', 'node'}:
            raise TypeError("col_belongs_to must be one of %s.", ', '.join(['edge', 'node']))

        if col_type not in self.allowed_types:
            raise TypeError("col_type must be one of %s.", ', '.join(self.allowed_types))

        if x is None:
            if col_belongs_to == 'edge':
                dat = self.obj.get_edge_column(column)
                x = (dat.min(), dat.max())
            else:
                dat = self.obj.get_node_column(column)
                x = (dat.min(), dat.max())
        if y is None:
            y = (10., 50.)
        slope = style.StyleUtil.create_slope(min=x[0], max=x[1], values=y)
        self.style.create_continuous_mapping(column=column, vp=prop, col_type='Double', points=slope)
        self.apply_style()

    def passthrough_node_fill(self, column, col_type='String'):
        self._create_passthrough_mapping(column, 'NODE_FILL_COLOR', col_type=col_type)

    def passthrough_node_label(self, column, col_type='String'):
        self._create_passthrough_mapping(column, 'NODE_LABEL', col_type=col_type)

    def passthrough_node_size_linear(self, column, xmax=None, xmin=0, ymin=10, ymax=50):
        if xmax is None:
            dat = self.obj.get_node_column(column)
            xmax = dat.max()
        self._create_linear_mapping(column, 'NODE_SIZE', x=(xmin, xmax), y=(ymin, ymax), col_belongs_to='node')

    def passthrough_edge_width_linear(self, column, xmax=None, xmin=None, ymin=0, ymax=5):
        if xmax is None:
            dat = self.obj.get_edge_column(column)
            xmax = dat.max()
        if xmin is None:
            dat = self.obj.get_edge_column(column)
            xmin = dat.min()
        self._create_linear_mapping(column, 'EDGE_WIDTH', x=(xmin, xmax), y=(ymin, ymax), col_belongs_to='edge')

    def update_style_defaults(self, prop_dict):
        self.style.update_defaults(prop_dict)

    def set_node_border_width(self, width):
        self.update_style_defaults({'NODE_BORDER_WIDTH': width})

    def set_node_label_colour(self, colour):
        self.update_style_defaults({'NODE_LABEL_COLOR': colour})

    def set_edge_colour(self, colour):
        self.update_style_defaults({'EDGE_STROKE_UNSELECTED_PAINT': colour})
