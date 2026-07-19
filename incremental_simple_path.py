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

    def solve(self, k: int):
        if k < 0 or k > self.G.number_of_edges():
            return False

        self.extend_to(k)

        return bool(self.solver.solve())

    def longest_simple_path(self):
        longest_path = []
        for k in range(1, min(self.G.number_of_nodes(), self.G.number_of_edges() + 1)):
            if self.solve(k):
                model = set(self.solver.get_model())  # type: ignore
                longest_path = [v for i in range(k + 1) for v in self.G.nodes() if self.vpool.id((v, i)) in model]
            else:
                break
        return longest_path

    def delete(self):
        self.solver.delete()
