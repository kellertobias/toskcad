//! Declared blind-bore wall and floor reuse the complete stock cylinder proof.
use super::*;
impl Solid {
    /// Drill a plain cap-normal blind bore with a fully protected wall/floor.
    /// The stock must have no prior curved cuts; their separation cannot be
    /// inferred from the base rectilinear boundary's minimum-feature guards.
    /// Returns the actual entry mouth in the current stock placement.
    pub fn rust_connection_blind_bore(
        &self,
        tool: &CircularPocketTool,
        minimum_wall_mm: f64,
        minimum_floor_mm: f64,
        policy: GeometryPolicy,
    ) -> Result<(Self, Point), GeometryError> {
        policy.validate()?;
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return Err(invalid(
                "A blind receiving bore requires the stock's exact geometry policy",
            ));
        }
        policy.length(minimum_wall_mm)?;
        policy.length(minimum_floor_mm)?;
        let protected_radius = (tool.radius_mm + minimum_wall_mm).next_up();
        policy.length(protected_radius)?;
        let start = circular_pocket_body::tool_start(tool, policy)?;
        let finish = std::array::from_fn(|i| start[i] + tool.displacement_mm[i]);
        let stock = self.orthogonal_stock()?;
        let (body, frame) = stock.orthogonal_body().expect("orthogonal stock promotion");
        if body.geometry.solid_count() != 1 || body.has_curved_cuts() {
            return Err(unsupported(
                "A blind receiving bore requires one orthogonal stock solid without prior curved cuts",
            ));
        }
        let [from, to] =
            connection_bore_axis(frame, tool, start, finish, protected_radius, policy)?;
        let pocket = body
            .identify_blind(from, to, tool.radius_mm, policy)
            .map_err(circular_pocket_body::pocket_error)?;
        let axis = (0..3)
            .max_by(|&a, &b| (to[a] - from[a]).abs().total_cmp(&(to[b] - from[b]).abs()))
            .expect("three axes");
        let positive = to[axis] > from[axis];
        // Both sums round outward from the actual retained pocket floor,
        // so neither the radius nor the declared floor can lose material to
        // rounded endpoint arithmetic. The existing intrinsic reserve remains.
        let mut protected_end = to;
        protected_end[axis] = if positive {
            ((pocket.center_mm[axis] + pocket.depth_mm).next_up() + minimum_floor_mm).next_up()
        } else {
            ((pocket.center_mm[axis] - pocket.depth_mm).next_down() - minimum_floor_mm).next_down()
        };
        let witness = body
            .identify_blind(from, protected_end, protected_radius, policy)
            .map_err(circular_pocket_body::pocket_error)?;
        if witness.face != pocket.face || witness.center_mm != pocket.center_mm {
            return Err(unsupported(
                "Blind bore and declared wall/floor must share the actual entry face",
            ));
        }
        let mouth = transform_point(frame, pocket.center_mm);
        let result = stock.with_orthogonal_round_cuts(
            vec![pocket],
            vec![tool.origin.clone()],
            vec![],
            vec![],
        )?;
        Ok((result, mouth))
    }
}
