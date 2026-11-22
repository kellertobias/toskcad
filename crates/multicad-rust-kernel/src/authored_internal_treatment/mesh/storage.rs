//! Actual retained allocation capacities, independent of geometric provenance.
use super::*;
pub(super) fn bytes(parts: &[(usize, usize)], base: usize, limit: usize) -> Result<usize, Error> {
    let mut bytes = base;
    for &(capacity, element) in parts {
        bytes = bytes
            .checked_add(capacity.checked_mul(element).ok_or(Error::Limit)?)
            .ok_or(Error::Limit)?;
    }
    if bytes > limit {
        return Err(Error::Limit);
    }
    Ok(bytes)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_capacity_budget_counts_headers_and_refuses_oversized_storage() {
        let mut actual = Vec::<u64>::new();
        actual.try_reserve_exact(1024).unwrap();
        let actual_bytes = actual.capacity() * std::mem::size_of::<u64>();
        assert_eq!(
            bytes(
                &[(actual.capacity(), std::mem::size_of::<u64>())],
                24,
                actual_bytes + 24
            )
            .unwrap(),
            actual_bytes + 24
        );
        assert_eq!(
            bytes(
                &[(actual.capacity(), std::mem::size_of::<u64>())],
                24,
                actual_bytes
            ),
            Err(Error::Limit)
        );
    }
}
