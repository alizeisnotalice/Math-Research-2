// The stored integers bound support at the grid knots, with the interpolation
// error already added. All interpolation and comparisons here are outward.
struct AngleMesh {
    const L* knots = nullptr;
    const L* values = nullptr;
    int panels = 0;

    L upper(int body, int direction, int normals, L lower, L upper) const {
        const L* x = knots + (body - 1) * (panels + 1);
        if (lower < x[0] || upper > x[panels] || lower > upper) std::abort();
        const L* values0 = values + (body - 1) * (panels + 1) * normals + direction;
        const auto at = [&](int j) { return values0[j * normals]; };
        const auto panel = [&](L angle) {
            return std::min(panels - 1,
                            int(std::upper_bound(x, x + panels + 1, angle) - x) - 1);
        };
        const auto interp = [&](L angle, int j) {
            const L width = x[j + 1] - x[j];
            const W weighted = W(at(j)) * (x[j + 1] - angle) +
                               W(at(j + 1)) * (angle - x[j]);
            return L(weighted >= 0 ? (weighted + width - 1) / width : weighted / width);
        };
        const int left = panel(lower);
        const int right = panel(upper);
        L bound = std::max(interp(lower, left), interp(upper, right));
        for (int j = left + 1; j <= right; ++j) bound = std::max(bound, at(j));
        return bound;
    }
};
