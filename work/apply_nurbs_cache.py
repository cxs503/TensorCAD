"""Apply the derivative-net cache patch to tensorcad/nurbs.py (idempotent)."""
from pathlib import Path

path = Path("/srv/tensor/agent/work/m22464/tensorcad/nurbs.py")
text = path.read_text()

old_post = """        object.__setattr__(self, "weights", checked)
        object.__setattr__(self, "knots", (knots_u, knots_v))
"""
new_post = """        object.__setattr__(self, "weights", checked)
        object.__setattr__(self, "knots", (knots_u, knots_v))
        # Derivative control nets are a function of the (immutable) control net and the
        # knot vectors only, so they are built once per surface instead of per evaluation.
        object.__setattr__(self, "_net_cache", {})
"""
assert text.count(old_post) == 1
text = text.replace(old_post, new_post)

anchor = """    # -- evaluation ---------------------------------------------------------
    def derivatives(
"""
helpers = '''    # -- evaluation ---------------------------------------------------------
    def _homogeneous_net(self) -> tuple:
        """The 4-D control net ``(w x, w y, w z, w)``; built once."""
        cache = self._net_cache
        net = cache.get("h")
        if net is None:
            net = tuple(
                tuple(
                    (
                        self.weights[i][j] * self.control_net[i][j][0],
                        self.weights[i][j] * self.control_net[i][j][1],
                        self.weights[i][j] * self.control_net[i][j][2],
                        self.weights[i][j],
                    )
                    for j in range(len(self.control_net[i]))
                )
                for i in range(len(self.control_net))
            )
            cache["h"] = net
        return net

    def _levels_u(self, order: int) -> list:
        """Control nets differentiated ``0..min(order, p)`` times along ``u``."""
        limit = min(order, self.degree[0])
        cache = self._net_cache
        levels = cache.get(("u",))
        if levels is None:
            levels = [(self._homogeneous_net(), self.degree[0], tuple(self.knots[0]))]
            cache[("u",)] = levels
        while len(levels) <= limit:
            reduced, degree, reduced_knots = _differentiate_axis(
                levels[-1][0], 0, levels[-1][1], levels[-1][2]
            )
            levels.append((reduced, degree, reduced_knots))
        return levels

    def _levels_v(self, k: int, order: int) -> list:
        """Control nets differentiated ``0..min(order-k, q)`` times along ``v`` for one ``u``-level."""
        limit = min(order - k, self.degree[1])
        cache = self._net_cache
        key = ("v", k)
        levels = cache.get(key)
        if levels is None:
            base_net = self._levels_u(k)[k][0]
            levels = [(base_net, self.degree[1], tuple(self.knots[1]))]
            cache[key] = levels
        while len(levels) <= limit:
            reduced, degree, reduced_knots = _differentiate_axis(
                levels[-1][0], 1, levels[-1][1], levels[-1][2]
            )
            levels.append((reduced, degree, reduced_knots))
        return levels

    def derivatives(
'''
assert text.count(anchor) == 1
text = text.replace(anchor, helpers)

old_net = '''        net: tuple[tuple[tuple[float, ...], ...], ...] = tuple(
            tuple(
                (
                    self.weights[i][j] * self.control_net[i][j][0],
                    self.weights[i][j] * self.control_net[i][j][1],
                    self.weights[i][j] * self.control_net[i][j][2],
                    self.weights[i][j],
                )
                for j in range(len(self.control_net[i]))
            )
            for i in range(len(self.control_net))
        )
        levels_u: list[tuple[object, int, tuple[float, ...]]] = [(net, p, tuple(self.knots[0]))]
        for _ in range(1, min(order, p) + 1):
            reduced, degree, reduced_knots = _differentiate_axis(
                levels_u[-1][0], 0, levels_u[-1][1], levels_u[-1][2]
            )
            levels_u.append((reduced, degree, reduced_knots))
'''
new_net = '''        levels_u = self._levels_u(order)
'''
assert text.count(old_net) == 1
text = text.replace(old_net, new_net)

old_v = '''            base_net, base_degree, base_knots = levels_u[k]
            levels_v: list[tuple[object, int, tuple[float, ...]]] = [
                (base_net, q, tuple(self.knots[1]))
            ]
            for _ in range(1, min(order - k, q) + 1):
                reduced, degree, reduced_knots = _differentiate_axis(
                    levels_v[-1][0], 1, levels_v[-1][1], levels_v[-1][2]
                )
                levels_v.append((reduced, degree, reduced_knots))
'''
new_v = '''            base_net, base_degree, base_knots = levels_u[k]
            levels_v = self._levels_v(k, order)
'''
assert text.count(old_v) == 1
text = text.replace(old_v, new_v)

path.write_text(text)
print("patched", path)