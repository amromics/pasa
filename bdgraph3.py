import networkx as nx
from networkx.classes.multigraph import MultiGraph

class BDGraph(MultiGraph):

    """
    Class represent the bidirected graph which is MultiGraph except each edge obtaining
    bi-direction feature.
    For a node: 0=out 1=in (direction w.r.t each tip node of an edge)
    As a consequence, the direction code for an edge can be: 0=00 1=01 2=10 3=11,
    representing 4 types of an bidirected edge:

    Topology   |     Edge      |  BDir  |   Reading
    -----------|---------------|-------|----------------
    u->---->-v | (u_out, v_in) |  01=1 |   u+v+ or v-u-
    u->----<-v | (u_out, v_out)|  00=0 |   u+v- or v+u-
    u-<---->-v | (u_in, v_in)  |  11=3 |   u-v+ or v-u+
    u-<----<-v | (u_in, v_out) |  10=2 |   u-v- or v+u+
    -----------|---------------|-------|----------------

    This class inherit MultiGraph with an additional *bdir* attribute for an edge encoded as
    abovementioned.
    There should be 3-level of granularity for a BDGraph:
        Level 1: edges with same (node pairs + bdir) will be merged.
        Level 2: edges with same (node pairs + bdir + length) will be merged.
        Level 3: edges with same (node pairs + bdir + seq) will be merged.

    """

    def __init__(self, granularity=1, incoming_graph_data=None, **attr):
        # customized dicts here
        MultiGraph.__init__(self, incoming_graph_data, **attr)
        self.granLevel = granularity
    #Overwrite methods
    ############################################################
    def new_edge_key(self, u, v, direction):
        """Returns an unused keypair for edges between nodes `u` and `v` with
        pre-defined bi-direction code `direction`.

        The nodes `u` and `v` do not need to be already in the graph.
        `direction` must be a valid code of (0, 1, 2 or 3)
        Notes
        -----
        Bi-directed edge, not like directed edge, can be traversed 2-ways: u -> v and v -> u
        The corresponding reading for node will be different but still pointing to the same edge.
        Thus for and edge we have a keypair instead of a single key code (also twin keys)

        Parameters
        ----------
        u, v : nodes
        direction: bi-direction code
        Returns
        -------
        keypair : int,int
        """
        if  direction not in (0,1,2,3):
            raise ValueError("Invalid direction identifier: must be 0, 1, 2 or 3!")

        try:
            keydict = self._adj[u][v]
        except KeyError:
            keydict = self.edge_key_dict_factory()

        key_u = len(keydict)//4 + direction
        while key_u in keydict:
            key_u += 4

        key_v = key_u + direction%2 - direction//2 #brilliant!

        return key_u, key_v

    # Add an edge with bi-direction
    def add_edge(self, u, v, direction=None, key=None, **attr):
        """
        Adding a bidirected edge with specified bi-direction to the graph.
        `direction` is the bi-direciton code from `u` to `v`

        This function will add the counterpart from `v` to `u` as well and return
        the pair of twin keys (pointing to the same edge object)
        """

        if  direction not in (None,0,1,2,3):
            raise ValueError("Invalid direction identifier: must be 0, 1, 2, 3 or None!")
        if (direction is None) == (key is None):
            raise ValueError("Must be exactly 1 None argument among direction and key!")

        # add nodes
        if u not in self._adj:
            self._adj[u] = self.adjlist_inner_dict_factory()
            self._node[u] = self.node_attr_dict_factory()
        if v not in self._adj:
            self._adj[v] = self.adjlist_inner_dict_factory()
            self._node[v] = self.node_attr_dict_factory()

        if key is None:
            key_u, key_v = self.new_edge_key(u, v, direction)
        else:
            direction = key%4
            key_u, key_v = key, key + direction%2 - direction//2

        #if (key_u in self._adj[u][v]) != (key_v in self._adj[v][u]):
        #    raise KeyError("Keypair's existence conflict!")

        #u->v
        if v in self._adj[u]:
            keydict = self._adj[u][v]
            datadict = keydict.get(key_u, self.edge_attr_dict_factory())
            datadict.update(attr)
            keydict[key_u] = datadict
        else:
            # selfloops work this way without special treatment
            datadict = self.edge_attr_dict_factory()
            datadict.update(attr)
            keydict = self.edge_key_dict_factory()
            keydict[key_u] = datadict
            self._adj[u][v] = keydict
            #self._adj[v][u] = keydict

        #v->u
        if u in self._adj[v]:
            keydict = self._adj[v][u]
            datadict = keydict.get(key_v, self.edge_attr_dict_factory())
            datadict.update(attr)
            keydict[key_v] = datadict
        else:
            # selfloops work this way without special treatment
            datadict = self.edge_attr_dict_factory()
            datadict.update(attr)
            keydict = self.edge_key_dict_factory()
            keydict[key_v] = datadict
            self._adj[v][u] = keydict

        return key_u, key_v

    #TODO: adapting...
    def add_edges_from(self, ebunch_to_add, **attr):
        """Add all the edges in ebunch_to_add.

        Parameters
        ----------
        ebunch_to_add : container of edges
            Each edge given in the container will be added to the
            graph. The edges can be:

                - 2-tuples (u, v) or
                - 3-tuples (u, v, d) for an edge data dict d, or
                - 3-tuples (u, v, k) for not iterable key k, or
                - 4-tuples (u, v, k, d) for an edge with data and key k

        attr : keyword arguments, optional
            Edge data (or labels or objects) can be assigned using
            keyword arguments.

        Returns
        -------
        A list of edge keys assigned to the edges in `ebunch`.

        See Also
        --------
        add_edge : add a single edge
        add_weighted_edges_from : convenient way to add weighted edges

        Notes
        -----
        Adding the same edge twice has no effect but any edge data
        will be updated when each duplicate edge is added.

        Edge attributes specified in an ebunch take precedence over
        attributes specified via keyword arguments.

        Default keys are generated using the method ``new_edge_key()``.
        This method can be overridden by subclassing the base class and
        providing a custom ``new_edge_key()`` method.

        Examples
        --------
        >>> G = nx.Graph()  # or DiGraph, MultiGraph, MultiDiGraph, etc
        >>> G.add_edges_from([(0, 1), (1, 2)])  # using a list of edge tuples
        >>> e = zip(range(0, 3), range(1, 4))
        >>> G.add_edges_from(e)  # Add the path graph 0-1-2-3

        Associate data to edges

        >>> G.add_edges_from([(1, 2), (2, 3)], weight=3)
        >>> G.add_edges_from([(3, 4), (1, 4)], label="WN2898")
        """
        keylist = []
        for e in ebunch_to_add:
            ne = len(e)
            if ne == 4:
                u, v, key, dd = e
            elif ne == 3:
                u, v, dd = e
                key = None
            elif ne == 2:
                u, v = e
                dd = {}
                key = None
            else:
                msg = f"Edge tuple {e} must be a 2-tuple, 3-tuple or 4-tuple."
                raise NetworkXError(msg)
            ddd = {}
            ddd.update(attr)
            try:
                ddd.update(dd)
            except (TypeError, ValueError):
                if ne != 3:
                    raise
                key = dd  # ne == 3 with 3rd value not dict, must be a key
            key = self.add_edge(u, v, key)
            self[u][v][key].update(ddd)
            keylist.append(key)
        return keylist

