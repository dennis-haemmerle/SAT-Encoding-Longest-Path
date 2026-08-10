import networkx as nx
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


class IncrementalSimplePathEncoder:

    def __init__(self, G: nx.Graph):
        self.G = G
        self.vpool = IDPool()
        self.solver = Solver(name="Cadical195")
        self.max_k = -1

    def to_set(self, val):
        if val is None:
            return set()
        elif isinstance(val, (list, set)):
            return set(val)
        else:
            return {val}

    def extend_to(self, k: int):
        if k <= self.max_k:
            return

        for i in range(self.max_k + 1, k + 1):
            # 1. Position i is occupied by exactly one node.
            lits = [self.vpool.id((v, i)) for v in self.G.nodes()]
            block = CardEnc.equals(lits=lits, bound=1, vpool=self.vpool, encoding=EncType.seqcounter)
            for clause in block.clauses:
                self.solver.add_clause(clause)

            # 2. Each node appears at most once.
            for v in self.G.nodes():
                # New position i is connected with all previous positions.
                for j in range(i):
                    self.solver.add_clause([-self.vpool.id((v, i)), -self.vpool.id((v, j))])

            # 3. Require that positions i-1 and i are connected by an edge.
            if i > 0:
                for u in self.G.nodes():
                    self.solver.add_clause([-self.vpool.id((u, i - 1))] + [self.vpool.id((v, i)) for v in self.G.neighbors(u)])

            self.max_k = i

    def solve(self, k: int, assumptions: list):
        if k < 0 or k > self.G.number_of_edges():
            return False

        self.extend_to(k)

        return bool(self.solver.solve(assumptions=assumptions))

    def longest_simple_path(self, start=None, end=None, symmetry=None):
        self.start = self.to_set(start)
        self.end = self.to_set(end)

        if start:
            self.solver.add_clause([self.vpool.id((v, 0)) for v in self.start])

        # Symmetry breaking
        if symmetry is not None and self.G.is_directed():
            orbit_groups = symmetry.get("orbit_groups", {})
            for orbit in orbit_groups.values():
                if len(orbit) <= 1 or any(v in orbit for v in self.end):
                    continue

                # Only the representative of each orbit is allowed to be the start node.
                starters = [v for v in orbit if v in self.start] if self.start else list(orbit)
                if not starters:
                    continue

                representative = min(starters)

                for v in orbit:
                    if v != representative:
                        self.solver.add_clause([-self.vpool.id((v, 0))])

        longest_path = []

        for k in range(1, min(self.G.number_of_nodes(), self.G.number_of_edges() + 1)):
            assumptions = []

            if self.end:
                if len(self.end) == 1:
                    assumptions.append(self.vpool.id((next(iter(self.end)), k)))
                else:
                    for v in self.G.nodes():
                        if v not in self.end:
                            assumptions.append(-self.vpool.id((v, k)))

            if self.solve(k, assumptions):
                model = set(self.solver.get_model())  # type: ignore
                longest_path = [v for i in range(k + 1) for v in self.G.nodes() if self.vpool.id((v, i)) in model]
            else:
                break
        return longest_path

    def delete(self):
        self.solver.delete()
