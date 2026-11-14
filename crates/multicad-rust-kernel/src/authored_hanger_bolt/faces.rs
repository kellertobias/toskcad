//! Actual source face roles, independent of sampling or finite topology export.
use super::*;
use crate::authored_internal_treatment::ThreadBandRole;
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceKey {
    MetricTip,
    MetricShoulder,
    WoodShoulder,
    WoodTip,
    PlainWall,
    MetricThread { start: u32, band: ThreadBandRole },
    WoodThread { start: u32, band: ThreadBandRole },
}
const CAPS: [FaceKey; 5] = [
    FaceKey::MetricTip,
    FaceKey::MetricShoulder,
    FaceKey::WoodShoulder,
    FaceKey::WoodTip,
    FaceKey::PlainWall,
];
const BANDS: [ThreadBandRole; 4] = [
    ThreadBandRole::Lower,
    ThreadBandRole::Crest,
    ThreadBandRole::Upper,
    ThreadBandRole::Root,
];
impl Hanger {
    pub fn face_count(&self) -> usize {
        5 + 4 * (self.metric().spec().starts + self.wood().spec().starts) as usize
    }
    pub fn face_key(&self, index: usize) -> Option<FaceKey> {
        if index < 5 {
            return Some(CAPS[index]);
        }
        let k = index - 5;
        let metric = 4 * self.metric().spec().starts as usize;
        if k < metric {
            Some(FaceKey::MetricThread {
                start: (k / 4) as u32,
                band: BANDS[k % 4],
            })
        } else {
            let k = k - metric;
            (k < 4 * self.wood().spec().starts as usize).then(|| FaceKey::WoodThread {
                start: (k / 4) as u32,
                band: BANDS[k % 4],
            })
        }
    }
    pub fn face_index(&self, key: FaceKey) -> Option<u32> {
        match key {
            FaceKey::MetricThread { start, band } if start < self.metric().spec().starts => {
                Some(5 + 4 * start + band_index(band))
            }
            FaceKey::WoodThread { start, band } if start < self.wood().spec().starts => {
                Some(5 + 4 * self.metric().spec().starts + 4 * start + band_index(band))
            }
            _ => CAPS.iter().position(|&x| x == key).map(|x| x as u32),
        }
    }
}
fn band_index(band: ThreadBandRole) -> u32 {
    match band {
        ThreadBandRole::Lower => 0,
        ThreadBandRole::Crest => 1,
        ThreadBandRole::Upper => 2,
        ThreadBandRole::Root => 3,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_role_roundtrip_uses_actual_starts_and_survives_dimension_and_hand_edits() {
        for left in [false, true] {
            let source = Hanger::new(super::super::tests::fixture(left, false), None).unwrap();
            assert_eq!(source.face_count(), 13);
            for i in 0..source.face_count() {
                let key = source.face_key(i).unwrap();
                assert_eq!(source.face_index(key), Some(i as u32));
            }
            assert_eq!(source.face_key(source.face_count()), None);
            assert_eq!(source.face_key(usize::MAX), None);
            assert_eq!(
                source.face_index(FaceKey::MetricThread {
                    start: 1,
                    band: ThreadBandRole::Root
                }),
                None
            );
            assert_eq!(
                source.face_index(FaceKey::WoodThread {
                    start: u32::MAX,
                    band: ThreadBandRole::Crest
                }),
                None
            );
            let mut spec = source.spec();
            spec.plain_length_mm = 3.1;
            let edited = Hanger::new(spec, None).unwrap();
            for i in 0..source.face_count() {
                assert_eq!(source.face_key(i), edited.face_key(i));
            }
        }
    }
}
