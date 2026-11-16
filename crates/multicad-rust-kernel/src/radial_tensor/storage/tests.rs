use super::*;
fn material(input: BoundaryInput) -> RadialBoundary {
    BoundaryOwnership::new(input, 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
#[test]
fn actual_source_retained_charge_counts_capacity_without_mutating_geometry() {
    let input = crate::radial_tensor::tests::fixture().input().clone();
    let ordinary = material(input.clone());
    let mut enlarged = input;
    let mut vertices = Vec::with_capacity(enlarged.vertices.len() + 4096);
    vertices.extend_from_slice(&enlarged.vertices);
    let capacity = vertices.capacity();
    enlarged.vertices = vertices;
    let roomy = material(enlarged);
    assert_eq!(roomy.input(), ordinary.input());
    let extra = (capacity - ordinary.input().vertices.capacity()) * std::mem::size_of::<Point>();
    assert_eq!(
        roomy.retained_budget_bytes(),
        ordinary.retained_budget_bytes() + extra
    );
    assert!(ordinary.retained_budget_bytes() > std::mem::size_of::<RadialBoundary>());
}
#[test]
fn pole_owner_identity_counts_shared_allocation_once_but_not_equal_new_data() {
    let body = crate::radial_tensor::tests::material();
    let tensor = body
        .input()
        .surfaces
        .iter()
        .find_map(|s| {
            if let Surface::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .unwrap();
    let clone = tensor.clone();
    let distinct = RetainedTensor::new(tensor.spans(), tensor.poles().to_vec()).unwrap();
    let mut owners = PoleOwners::new();
    let charge = owners.allocation(tensor);
    assert_eq!(owners.allocation(&clone), 0);
    assert_eq!(owners.allocation(&distinct), charge);
    assert!(charge >= std::mem::size_of_val(tensor.poles()));
    assert_eq!(owners.count, 2);
}
